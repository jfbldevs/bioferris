use std::any::TypeId;
use std::ops::Range;

use align_core::Alphabet;

use super::{
    AnnotatedRecord, BioIoError, Feature, Format, IoErrorKind, Location, LocationExpression,
    MetadataField, RecordPart, parse_document, writer,
};
use crate::{Dna, DnaSequence, NucleotideSequenceExt, Record, Rna, RnaSequence};

/// Whether explicitly reviewed, removable annotation losses are permitted.
///
/// Discarding never authorizes clipping a feature or interpreting an unsupported
/// location. A feature with an unreviewed qualifier is removed as a whole.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AnnotationLossPolicy {
    /// Reject a transformation that would discard annotation (the default).
    #[default]
    Reject,
    /// Permit reported metadata and whole-feature removal.
    Discard,
}

/// Checked identity, date and explicit loss policy for a derived local record.
///
/// No wall-clock date or database accession is inferred. A date is required for
/// `GenBank`; for EMBL a supplied date adds a local-derivation `DT` field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnotationTransformOptions {
    local_id: String,
    loss_policy: AnnotationLossPolicy,
    date: Option<String>,
}

impl AnnotationTransformOptions {
    /// Creates strict options with a distinct caller-selected local name.
    ///
    /// # Errors
    /// Returns `InvalidHeader` unless the name starts with `LOCAL_`, has a
    /// nonempty ASCII alphanumeric/underscore suffix, and is at most 16 bytes.
    /// Reuse of a source identity is checked when the transformation is called.
    pub fn new(local_id: impl Into<String>) -> Result<Self, BioIoError> {
        let local_id = local_id.into();
        let suffix = local_id.strip_prefix("LOCAL_").unwrap_or_default();
        if suffix.is_empty()
            || local_id.len() > 16
            || !suffix
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
        {
            return Err(transform_error(
                IoErrorKind::InvalidHeader,
                "local identity must be LOCAL_ followed by ASCII letters, digits or underscores, with a total length of 7..=16 bytes",
            ));
        }
        Ok(Self {
            local_id,
            loss_policy: AnnotationLossPolicy::Reject,
            date: None,
        })
    }

    /// Selects whether explicitly reported metadata/whole-feature loss is allowed.
    #[must_use]
    pub const fn with_loss_policy(mut self, policy: AnnotationLossPolicy) -> Self {
        self.loss_policy = policy;
        self
    }

    /// Sets a caller-supplied date for this local derivation, not database history.
    ///
    /// # Errors
    /// Returns `InvalidHeader` unless this is a real Gregorian date in exactly
    /// `DD-MMM-YYYY` form, with an uppercase English month and year 0001..=9999.
    pub fn with_date(mut self, date: impl Into<String>) -> Result<Self, BioIoError> {
        let date = date.into();
        validate_date(&date)?;
        self.date = Some(date);
        Ok(self)
    }

    /// Returns the requested local locus name, without EMBL's local `.1` suffix.
    #[must_use]
    pub fn local_id(&self) -> &str {
        &self.local_id
    }

    /// Returns the explicit annotation-loss policy.
    #[must_use]
    pub const fn loss_policy(&self) -> AnnotationLossPolicy {
        self.loss_policy
    }

    /// Returns the explicitly supplied local derivation date, if any.
    #[must_use]
    pub fn date(&self) -> Option<&str> {
        self.date.as_deref()
    }
}

/// The operation and original-input coordinate system of a derived record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnotationOperation {
    /// Copies a strict zero-based half-open source interval.
    Slice {
        /// Inclusive source start.
        start: usize,
        /// Exclusive source end.
        end: usize,
    },
    /// Reverse-complements a DNA or RNA record without changing its length.
    ReverseComplement,
}

/// One explicit semantic change, using indexes in the unmodified source model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnnotationChange {
    /// Replaces, removes or adds one physical metadata field.
    MetadataChanged {
        /// Original zero-based metadata index; `None` for a new field.
        index: Option<usize>,
        /// Metadata field name.
        key: String,
        /// Original payload, or `None` for an addition.
        before: Option<String>,
        /// Derived payload, or `None` for a removal.
        after: Option<String>,
        /// Why this change is necessary or explicitly permitted.
        reason: String,
    },
    /// Changes a retained feature's canonical location.
    LocationChanged {
        /// Original zero-based feature index.
        feature_index: usize,
        /// Source location text, including any original wrapping.
        before: String,
        /// Canonical derived location text.
        after: String,
        /// Whether this is the special whole-record `source` coverage update.
        source_coverage: bool,
    },
    /// Removes an entire feature under an explicit discard policy.
    FeatureRemoved {
        /// Original zero-based feature index.
        index: usize,
        /// Source feature key.
        key: String,
        /// Source location text.
        location: String,
        /// The loss that required removal.
        reason: String,
    },
}

/// Ordered, typed provenance of an annotation-preserving transformation.
///
/// Changes to existing entries follow source layout order. New metadata fields
/// follow those changes in their derived layout order. This is a Rust readout,
/// not a versioned serialized interchange schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnotationTransformReport {
    /// Original derived record ID, including its version if present.
    pub source_id: String,
    /// Actual output record ID (`LOCAL_name.1` for modern EMBL).
    pub derived_id: String,
    /// The requested operation in source coordinates.
    pub operation: AnnotationOperation,
    /// Every changed location and replaced/removed/added metadata or feature.
    pub changes: Vec<AnnotationChange>,
}

/// A derived, reparsed record and the explicit report explaining its changes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnotatedTransform<A: Alphabet> {
    record: AnnotatedRecord<A>,
    report: AnnotationTransformReport,
}

impl<A: Alphabet> AnnotatedTransform<A> {
    /// Borrows the derived record, whose raw text has also been regenerated.
    #[must_use]
    pub const fn record(&self) -> &AnnotatedRecord<A> {
        &self.record
    }

    /// Borrows the ordered transformation report.
    #[must_use]
    pub const fn report(&self) -> &AnnotationTransformReport {
        &self.report
    }

    /// Consumes the result and returns the derived record.
    #[must_use]
    pub fn into_record(self) -> AnnotatedRecord<A> {
        self.record
    }
}

impl<A: Alphabet> AnnotatedRecord<A> {
    /// Copies a strict interval while mapping complete, reviewed annotations.
    ///
    /// Exact ranges, joins and complements are mapped without reordering joined
    /// segments. Partially retained ordinary features always fail: discard only
    /// permits fully disjoint features and features with unreviewed qualifiers
    /// to be removed as whole features. A simple whole-record `source` feature
    /// becomes the new nonempty coverage; removing it for an empty slice requires
    /// discard. A proper slice of a circular record becomes linear.
    ///
    /// Source case is preserved. The original record is never changed; the
    /// result has a new local identity and regenerated metadata/raw text.
    ///
    /// # Errors
    /// Returns `LengthMismatch` for invalid bounds; `InvalidHeader` for identity
    /// reuse or a missing required `GenBank` date; `Unsupported` for unknown
    /// topology, classification or any unsupported location; `DataLoss` for partial features or
    /// policy-rejected losses. Rendering/reparse errors propagate explicitly.
    pub fn slice(
        &self,
        range: Range<usize>,
        options: AnnotationTransformOptions,
    ) -> Result<AnnotatedTransform<A>, BioIoError> {
        if range.start > range.end || range.end > self.record.sequence().len() {
            return Err(transform_error(
                IoErrorKind::LengthMismatch,
                format!(
                    "slice [{}..{}) is outside the record length {}",
                    range.start,
                    range.end,
                    self.record.sequence().len()
                ),
            ));
        }
        transform(
            self,
            AnnotationOperation::Slice {
                start: range.start,
                end: range.end,
            },
            options,
        )
    }

    /// Reverse-complements DNA/RNA and maps every exact retained location.
    ///
    /// Written join order is retained and double complements are normalized.
    /// Complete retained non-source feature extraction is unchanged. Source
    /// letter case is reversed with its corresponding nucleotide. The source
    /// object is not mutated, and topology is retained.
    ///
    /// # Errors
    /// Returns `Unsupported` for protein/custom alphabets, unknown topology/classification or
    /// unsupported locations; `InvalidHeader` for reused identity or a missing
    /// required `GenBank` date; `DataLoss` for policy-rejected annotation loss.
    /// Rendering/reparse errors propagate explicitly.
    pub fn reverse_complement(
        &self,
        options: AnnotationTransformOptions,
    ) -> Result<AnnotatedTransform<A>, BioIoError> {
        if !is_nucleotide::<A>() {
            return Err(transform_error(
                IoErrorKind::Unsupported,
                "annotated reverse complement supports only the built-in DNA and RNA alphabets",
            ));
        }
        transform(self, AnnotationOperation::ReverseComplement, options)
    }
}

fn transform<A: Alphabet>(
    source: &AnnotatedRecord<A>,
    operation: AnnotationOperation,
    options: AnnotationTransformOptions,
) -> Result<AnnotatedTransform<A>, BioIoError> {
    let derived_id = validate_identity(source, &options)?;
    let topology = topology(source)?;
    let original_length = source.record.sequence().len();
    let derived_topology = match operation {
        AnnotationOperation::Slice { start, end } if start != 0 || end != original_length => {
            "linear"
        }
        _ => topology,
    };
    // Validate every original location before allowing a discard decision. A
    // suspicious qualifier must not conceal unsupported or clipped coordinates.
    for feature in &source.features {
        validate_location(
            feature.location.expression(),
            original_length,
            is_nucleotide::<A>(),
        )?;
        if let AnnotationOperation::Slice { start, end } = operation {
            if !whole_source(feature, original_length)
                && relation(feature.location.expression(), &(start..end)) == Relation::Partial
            {
                return Err(transform_error(
                    IoErrorKind::DataLoss,
                    format!(
                        "slice would partially retain feature {} at {}; feature clipping is unsupported even with discard",
                        feature.key,
                        feature.location.raw()
                    ),
                ));
            }
        }
    }
    validate_header_tail(source, options.loss_policy)?;
    let symbols = transformed_symbols(source, &operation)?;
    let description = match operation {
        AnnotationOperation::Slice { start, end } => {
            format!("BioFerris slice [{start},{end}) of {}.", source.record.id())
        }
        AnnotationOperation::ReverseComplement => {
            format!("BioFerris reverse complement of {}.", source.record.id())
        }
    };
    let record = Record::from_ascii(format!("{derived_id} {description}"), &symbols, None)?;
    let mut model = AnnotatedRecord {
        format: source.format,
        record,
        locus: String::new(),
        molecule: source.molecule.clone(),
        accessions: Vec::new(),
        metadata: Vec::new(),
        features: Vec::new(),
        raw_record: String::new(),
        layout: Vec::new(),
    };
    let mut changes = Vec::new();
    populate_model(
        source,
        &mut model,
        &options,
        &operation,
        derived_topology,
        &description,
        &mut changes,
    )?;
    model.locus = options.local_id;
    let rendered = writer::render(&model, source.format)?;
    let mut records = parse_document(&rendered, source.format)?;
    if records.len() != 1 {
        return Err(transform_error(
            IoErrorKind::InvalidStructure,
            "transformed output did not contain exactly one record",
        ));
    }
    let record = records.pop().ok_or_else(|| {
        transform_error(
            IoErrorKind::InvalidStructure,
            "transformed output is unexpectedly empty",
        )
    })?;
    Ok(AnnotatedTransform {
        record,
        report: AnnotationTransformReport {
            source_id: source.record.id().to_owned(),
            derived_id,
            operation,
            changes,
        },
    })
}

fn validate_identity<A: Alphabet>(
    source: &AnnotatedRecord<A>,
    options: &AnnotationTransformOptions,
) -> Result<String, BioIoError> {
    if source.format == Format::GenBank && options.date.is_none() {
        return Err(transform_error(
            IoErrorKind::InvalidHeader,
            "GenBank transformations require an explicit local derivation date via with_date",
        ));
    }
    let derived_id = match source.format {
        Format::GenBank => options.local_id.clone(),
        Format::Embl => format!("{}.1", options.local_id),
    };
    for identity in [source.locus.as_str(), source.record.id()]
        .into_iter()
        .chain(source.accessions.iter().map(String::as_str))
    {
        if identity == options.local_id || identity == derived_id {
            return Err(transform_error(
                IoErrorKind::InvalidHeader,
                "derived local identity must differ from every source locus, ID and accession",
            ));
        }
    }
    Ok(derived_id)
}

fn topology<A: Alphabet>(record: &AnnotatedRecord<A>) -> Result<&'static str, BioIoError> {
    let header = record.metadata.first().ok_or_else(|| {
        transform_error(
            IoErrorKind::InvalidStructure,
            "record has no header metadata",
        )
    })?;
    let token = match record.format {
        Format::GenBank => header
            .value
            .split_whitespace()
            .nth(if record.molecule == "protein" { 3 } else { 4 }),
        Format::Embl => header.value.split(';').nth(2).map(str::trim),
    };
    match token {
        Some("linear") => Ok("linear"),
        Some("circular") => Ok("circular"),
        _ => Err(transform_error(
            IoErrorKind::Unsupported,
            "annotated transformations require an explicit linear or circular source topology",
        )),
    }
}

fn transformed_symbols<A: Alphabet>(
    source: &AnnotatedRecord<A>,
    operation: &AnnotationOperation,
) -> Result<String, BioIoError> {
    match *operation {
        AnnotationOperation::Slice { start, end } => source
            .record
            .slice(start..end)
            .map(|record| record.source_sequence().to_owned())
            .map_err(|error| transform_error(IoErrorKind::LengthMismatch, error.to_string())),
        AnnotationOperation::ReverseComplement => {
            let bytes = source.record.source_sequence().as_bytes();
            let reversed = if TypeId::of::<A>() == TypeId::of::<Dna>() {
                DnaSequence::from_ascii(bytes)
                    .and_then(|sequence| sequence.reverse_complement())
                    .map(|sequence| sequence.to_string())
            } else if TypeId::of::<A>() == TypeId::of::<Rna>() {
                RnaSequence::from_ascii(bytes)
                    .and_then(|sequence| sequence.reverse_complement())
                    .map(|sequence| sequence.to_string())
            } else {
                return Err(transform_error(
                    IoErrorKind::Unsupported,
                    "annotated reverse complement supports only the built-in DNA and RNA alphabets",
                ));
            }
            .map_err(|error| transform_error(IoErrorKind::InvalidSequence, error.to_string()))?;
            Ok(reversed
                .bytes()
                .zip(source.record.source_sequence().bytes().rev())
                .map(|(byte, original)| {
                    char::from(if original.is_ascii_lowercase() {
                        byte.to_ascii_lowercase()
                    } else {
                        byte
                    })
                })
                .collect())
        }
    }
}

fn is_nucleotide<A: Alphabet>() -> bool {
    TypeId::of::<A>() == TypeId::of::<Dna>() || TypeId::of::<A>() == TypeId::of::<Rna>()
}

fn validate_header_tail<A: Alphabet>(
    source: &AnnotatedRecord<A>,
    policy: AnnotationLossPolicy,
) -> Result<(), BioIoError> {
    if source.format == Format::GenBank {
        let header = source.metadata.first().ok_or_else(layout_error)?;
        let known_fields = if source.molecule == "protein" { 6 } else { 7 };
        if header.value.split_whitespace().count() > known_fields {
            require_discard(policy, "unreviewed trailing LOCUS header content")?;
        }
    }
    Ok(())
}

fn validate_location(
    expression: &LocationExpression,
    length: usize,
    nucleotide: bool,
) -> Result<(), BioIoError> {
    match expression {
        LocationExpression::Range(range) => {
            if range.start >= range.end || range.end > length {
                return Err(transform_error(
                    IoErrorKind::LengthMismatch,
                    "feature location is outside the original record",
                ));
            }
            Ok(())
        }
        LocationExpression::Join(parts) => {
            for part in parts {
                validate_location(part, length, nucleotide)?;
            }
            Ok(())
        }
        LocationExpression::Complement(part) => {
            if !nucleotide {
                return Err(transform_error(
                    IoErrorKind::Unsupported,
                    "complemented features require the built-in DNA or RNA alphabet, including during slicing",
                ));
            }
            validate_location(part, length, nucleotide)
        }
        LocationExpression::Unsupported { reason, .. } => Err(transform_error(
            IoErrorKind::Unsupported,
            format!("cannot transform an unsupported feature location: {reason}"),
        )),
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Relation {
    Contained,
    Disjoint,
    Partial,
}

fn relation(expression: &LocationExpression, slice: &Range<usize>) -> Relation {
    if slice.is_empty() {
        return Relation::Disjoint;
    }
    match expression {
        LocationExpression::Range(range) => {
            if range.start >= slice.start && range.end <= slice.end {
                Relation::Contained
            } else if range.end <= slice.start || range.start >= slice.end {
                Relation::Disjoint
            } else {
                Relation::Partial
            }
        }
        LocationExpression::Join(parts) => {
            let mut combined = None;
            for part in parts {
                let current = relation(part, slice);
                if current == Relation::Partial
                    || combined.is_some_and(|previous| previous != current)
                {
                    return Relation::Partial;
                }
                combined = Some(current);
            }
            combined.unwrap_or(Relation::Partial)
        }
        LocationExpression::Complement(part) => relation(part, slice),
        LocationExpression::Unsupported { .. } => Relation::Partial,
    }
}

fn whole_source(feature: &Feature, length: usize) -> bool {
    feature.key == "source"
        && matches!(feature.location.expression(), LocationExpression::Range(range) if *range == (0..length))
}

fn map_location(
    expression: &LocationExpression,
    operation: &AnnotationOperation,
    length: usize,
) -> Result<LocationExpression, BioIoError> {
    match expression {
        LocationExpression::Range(range) => match *operation {
            AnnotationOperation::Slice { start, .. } => Ok(LocationExpression::Range(
                checked_sub(range.start, start)?..checked_sub(range.end, start)?,
            )),
            AnnotationOperation::ReverseComplement => Ok(LocationExpression::Complement(Box::new(
                LocationExpression::Range(
                    checked_sub(length, range.end)?..checked_sub(length, range.start)?,
                ),
            ))),
        },
        LocationExpression::Join(parts) => parts
            .iter()
            .map(|part| map_location(part, operation, length))
            .collect::<Result<Vec<_>, _>>()
            .map(LocationExpression::Join),
        LocationExpression::Complement(part) => {
            let mapped = map_location(part, operation, length)?;
            Ok(match mapped {
                LocationExpression::Complement(inner) => *inner,
                other => LocationExpression::Complement(Box::new(other)),
            })
        }
        LocationExpression::Unsupported { reason, .. } => Err(transform_error(
            IoErrorKind::Unsupported,
            format!("cannot transform an unsupported feature location: {reason}"),
        )),
    }
}

fn mapped_feature(
    feature: &Feature,
    index: usize,
    operation: &AnnotationOperation,
    original_length: usize,
    output_length: usize,
    policy: AnnotationLossPolicy,
    format: Format,
) -> Result<(Option<Feature>, Option<AnnotationChange>), BioIoError> {
    let coverage = whole_source(feature, original_length);
    let outside = match *operation {
        AnnotationOperation::Slice { start, end } => {
            if coverage {
                output_length == 0
            } else {
                relation(feature.location.expression(), &(start..end)) == Relation::Disjoint
            }
        }
        AnnotationOperation::ReverseComplement => false,
    };
    let unreviewed: Vec<_> = feature
        .qualifiers
        .iter()
        .filter(|qualifier| !reviewed_qualifier(&qualifier.name))
        .map(|qualifier| qualifier.name.as_str())
        .collect();
    let removal_reason = if outside {
        Some(if coverage {
            "whole-record source coverage cannot be retained for an empty slice".to_owned()
        } else {
            "feature is wholly outside the selected interval".to_owned()
        })
    } else if !unreviewed.is_empty() {
        Some(format!(
            "feature has unreviewed qualifiers: {}; the complete feature must be removed together with dependent qualifiers",
            unreviewed.join(", ")
        ))
    } else {
        None
    };
    if let Some(reason) = removal_reason {
        require_discard(policy, &reason)?;
        return Ok((
            None,
            Some(AnnotationChange::FeatureRemoved {
                index,
                key: feature.key.clone(),
                location: feature.location.raw().to_owned(),
                reason,
            }),
        ));
    }
    let expression = if coverage {
        LocationExpression::Range(0..output_length)
    } else {
        map_location(feature.location.expression(), operation, original_length)?
    };
    let text = writer::render_location(&expression, format)?;
    // Reparse to enforce the same nesting/resource/coordinate grammar as input.
    let location = Location::parse(&text)?;
    let change = (text != feature.location.raw()).then(|| AnnotationChange::LocationChanged {
        feature_index: index,
        before: feature.location.raw().to_owned(),
        after: text,
        source_coverage: coverage,
    });
    Ok((
        Some(Feature {
            key: feature.key.clone(),
            location,
            qualifiers: feature.qualifiers.clone(),
            raw_lines: Vec::new(),
        }),
        change,
    ))
}

fn reviewed_qualifier(name: &str) -> bool {
    matches!(
        name,
        "gene"
            | "gene_synonym"
            | "locus_tag"
            | "product"
            | "protein_id"
            | "codon_start"
            | "transl_table"
            | "translation"
            | "pseudo"
            | "pseudogene"
            | "organism"
            | "mol_type"
            | "strain"
            | "isolate"
            | "db_xref"
            | "culture_collection"
            | "plasmid"
            | "organelle"
    )
}

struct MetadataPlan {
    replacement: Option<MetadataField>,
    reason: &'static str,
}

#[allow(clippy::too_many_arguments)]
fn populate_model<A: Alphabet>(
    source: &AnnotatedRecord<A>,
    output: &mut AnnotatedRecord<A>,
    options: &AnnotationTransformOptions,
    operation: &AnnotationOperation,
    topology: &str,
    description: &str,
    changes: &mut Vec<AnnotationChange>,
) -> Result<(), BioIoError> {
    let mut description_seen = false;
    let has_description = source
        .metadata
        .iter()
        .any(|field| is_description(source.format, &field.key));
    let has_feature_header = source.metadata.iter().any(|field| match source.format {
        Format::GenBank => field.key == "FEATURES",
        Format::Embl => field.key == "FH" && recognized_feature_header(field, source.format),
    });
    let mut inserted_feature_header = has_feature_header;
    let mut added = Vec::new();
    for part in &source.layout {
        match *part {
            RecordPart::Metadata(index) => {
                let field = source.metadata.get(index).ok_or_else(layout_error)?;
                let plan = metadata_plan(
                    field,
                    index,
                    output,
                    options,
                    topology,
                    description,
                    &mut description_seen,
                )?;
                if plan.replacement.as_ref().is_none_or(|after| {
                    after.value != field.value || after.continuation != field.continuation
                }) {
                    changes.push(AnnotationChange::MetadataChanged {
                        index: Some(index),
                        key: field.key.clone(),
                        before: Some(field.value.clone()),
                        after: plan.replacement.as_ref().map(|after| after.value.clone()),
                        reason: plan.reason.to_owned(),
                    });
                }
                if is_sequence_header(source.format, &field.key) && !inserted_feature_header {
                    add_feature_header(output, &mut added);
                    inserted_feature_header = true;
                }
                if let Some(replacement) = plan.replacement {
                    push_metadata(output, replacement);
                }
                if index == 0 {
                    if !has_description {
                        let key = description_key(source.format);
                        add_metadata(
                            output,
                            metadata(key, description, source.format),
                            "describe the explicit local derivation",
                            &mut added,
                        );
                    }
                    if source.format == Format::Embl {
                        if let Some(date) = options.date() {
                            add_metadata(
                                output,
                                metadata(
                                    "DT",
                                    &format!("{date} (BioFerris local derivation)"),
                                    source.format,
                                ),
                                "record the caller-supplied local derivation date, not a database update",
                                &mut added,
                            );
                        }
                    }
                }
            }
            RecordPart::Feature(index) => {
                let feature = source.features.get(index).ok_or_else(layout_error)?;
                let (mapped, change) = mapped_feature(
                    feature,
                    index,
                    operation,
                    source.record.sequence().len(),
                    output.record.sequence().len(),
                    options.loss_policy,
                    source.format,
                )?;
                if let Some(change) = change {
                    changes.push(change);
                }
                if let Some(feature) = mapped {
                    if !inserted_feature_header {
                        add_feature_header(output, &mut added);
                        inserted_feature_header = true;
                    }
                    output
                        .layout
                        .push(RecordPart::Feature(output.features.len()));
                    output.features.push(feature);
                }
            }
        }
    }
    changes.extend(added);
    Ok(())
}

fn metadata_plan<A: Alphabet>(
    field: &MetadataField,
    index: usize,
    output: &AnnotatedRecord<A>,
    options: &AnnotationTransformOptions,
    topology: &str,
    description: &str,
    description_seen: &mut bool,
) -> Result<MetadataPlan, BioIoError> {
    let format = output.format;
    let make = |value: String, reason| MetadataPlan {
        replacement: Some(metadata(&field.key, &value, format)),
        reason,
    };
    if index == 0 {
        return Ok(make(
            local_header(output, options, topology, field)?,
            "replace database identity, length and topology with an explicit local derivation header, preserving source classification and removing any explicitly discarded opaque tail",
        ));
    }
    if is_description(format, &field.key) {
        let replacement = if *description_seen {
            None
        } else {
            *description_seen = true;
            Some(metadata(&field.key, description, format))
        };
        return Ok(MetadataPlan {
            replacement,
            reason: "replace the original complete-record description with explicit derivation provenance",
        });
    }
    let identity = match format {
        Format::GenBank => matches!(field.key.as_str(), "ACCESSION" | "VERSION"),
        Format::Embl => matches!(field.key.as_str(), "AC" | "SV"),
    };
    if identity {
        return Ok(MetadataPlan {
            replacement: None,
            reason: "remove the original database accession/version; the derived record has an explicit local identity",
        });
    }
    let preserved = match format {
        Format::GenBank => matches!(field.key.as_str(), "SOURCE" | "ORGANISM"),
        Format::Embl => matches!(field.key.as_str(), "OS" | "OC" | "OG"),
    };
    if preserved || (format == Format::Embl && field.key == "XX" && field.value.is_empty()) {
        return Ok(MetadataPlan {
            replacement: Some(field.clone()),
            reason: "preserve reviewed source metadata",
        });
    }
    if is_sequence_header(format, &field.key) {
        return Ok(make(
            sequence_header(field, output, options.loss_policy)?,
            "regenerate the literal sequence header, length and composition for the derived record, removing any explicitly discarded opaque content",
        ));
    }
    let is_feature_header = match format {
        Format::GenBank => field.key == "FEATURES",
        Format::Embl => field.key == "FH",
    };
    if is_feature_header {
        if recognized_feature_header(field, format) {
            return Ok(make(
                field.value.clone(),
                "regenerate the recognized feature table heading",
            ));
        }
        require_discard(
            options.loss_policy,
            &format!("unrecognized {} feature heading", field.key),
        )?;
        return Ok(MetadataPlan {
            replacement: (format == Format::GenBank)
                .then(|| metadata("FEATURES", "Location/Qualifiers", format)),
            reason: "discard unrecognized feature heading content and regenerate the supported heading",
        });
    }
    if format == Format::Embl && field.key == "DT" && is_local_date(&field.value) {
        return Ok(MetadataPlan {
            replacement: None,
            reason: "remove the previous local derivation date; only the new caller-supplied date applies",
        });
    }
    require_discard(
        options.loss_policy,
        &format!("unreviewed record metadata {} at index {index}", field.key),
    )?;
    Ok(MetadataPlan {
        replacement: None,
        reason: "discard unreviewed source metadata that cannot be certified for the derived record",
    })
}

fn sequence_header<A: Alphabet>(
    field: &MetadataField,
    output: &AnnotatedRecord<A>,
    policy: AnnotationLossPolicy,
) -> Result<String, BioIoError> {
    match output.format {
        Format::GenBank => {
            if !field.value.is_empty() {
                require_discard(policy, "unreviewed ORIGIN header content")?;
            }
            Ok(String::new())
        }
        Format::Embl => {
            match writer::sq_value(&field.value, output.record.source_sequence(), output.format) {
                Ok(value) => Ok(value),
                Err(error) if error.kind() == IoErrorKind::Unsupported => {
                    require_discard(policy, "unreviewed SQ sequence-header content")?;
                    writer::sq_value(
                        "Sequence 0 BP;",
                        output.record.source_sequence(),
                        output.format,
                    )
                }
                Err(error) => Err(error),
            }
        }
    }
}

fn local_header<A: Alphabet>(
    output: &AnnotatedRecord<A>,
    options: &AnnotationTransformOptions,
    topology: &str,
    source_header: &MetadataField,
) -> Result<String, BioIoError> {
    let length = output.record.sequence().len();
    let division = source_division(source_header, output.format, output.molecule == "protein")?;
    match output.format {
        Format::GenBank => {
            let date = options.date().ok_or_else(|| {
                transform_error(IoErrorKind::InvalidHeader, "missing local derivation date")
            })?;
            let unit_and_molecule = if output.molecule == "protein" {
                "aa".to_owned()
            } else {
                format!("bp {}", output.molecule)
            };
            Ok(format!(
                "{} {length} {unit_and_molecule} {topology} {division} {date}",
                options.local_id
            ))
        }
        Format::Embl => {
            let data_class = source_header.value.split(';').nth(4).map(str::trim);
            let Some(
                data_class
                @ ("STD" | "PAT" | "EST" | "GSS" | "HTC" | "HTG" | "WGS" | "TSA" | "STS"),
            ) = data_class
            else {
                return Err(transform_error(
                    IoErrorKind::Unsupported,
                    "annotated transformations require a recognized sequence-bearing EMBL data class; missing, unknown and CON classes are unsupported",
                ));
            };
            Ok(format!(
                "{}; SV 1; {topology}; {}; {data_class}; {division}; {length} BP.",
                options.local_id, output.molecule
            ))
        }
    }
}

fn source_division(
    header: &MetadataField,
    format: Format,
    protein: bool,
) -> Result<&str, BioIoError> {
    // Division codes describe the source, not whether this file was generated
    // locally. In particular, SYN must never be inferred from LOCAL_ identity.
    // The format-specific profiles follow the NCBI field guide and ENA manual.
    let (division, recognized): (Option<&str>, &[&str]) = match format {
        Format::GenBank => (
            header
                .value
                .split_whitespace()
                .nth(if protein { 4 } else { 5 }),
            &[
                "PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN", "UNA", "EST",
                "PAT", "STS", "GSS", "HTG", "HTC", "ENV",
            ],
        ),
        Format::Embl => (
            header.value.split(';').nth(5).map(str::trim),
            &[
                "PHG", "ENV", "FUN", "HUM", "INV", "MAM", "VRT", "MUS", "PLN", "PRO", "ROD", "SYN",
                "TGN", "UNC", "VRL",
            ],
        ),
    };
    division.filter(|code| recognized.contains(code)).ok_or_else(|| {
        transform_error(
            IoErrorKind::Unsupported,
            "annotated transformations require a recognized source division in the selected format; no classification is inferred",
        )
    })
}

fn metadata(key: &str, value: &str, format: Format) -> MetadataField {
    let width = match format {
        Format::GenBank => 12,
        Format::Embl => 5,
    };
    MetadataField {
        key: key.to_owned(),
        value: value.to_owned(),
        continuation: false,
        raw_line: format!("{key:<width$}{value}"),
    }
}

fn push_metadata<A: Alphabet>(record: &mut AnnotatedRecord<A>, field: MetadataField) {
    record
        .layout
        .push(RecordPart::Metadata(record.metadata.len()));
    record.metadata.push(field);
}

fn add_metadata<A: Alphabet>(
    record: &mut AnnotatedRecord<A>,
    field: MetadataField,
    reason: &str,
    changes: &mut Vec<AnnotationChange>,
) {
    changes.push(AnnotationChange::MetadataChanged {
        index: None,
        key: field.key.clone(),
        before: None,
        after: Some(field.value.clone()),
        reason: reason.to_owned(),
    });
    push_metadata(record, field);
}

fn add_feature_header<A: Alphabet>(
    record: &mut AnnotatedRecord<A>,
    changes: &mut Vec<AnnotationChange>,
) {
    let (key, value) = match record.format {
        Format::GenBank => ("FEATURES", "Location/Qualifiers"),
        Format::Embl => ("FH", "Key             Location/Qualifiers"),
    };
    add_metadata(
        record,
        metadata(key, value, record.format),
        "add the supported feature table heading",
        changes,
    );
}

fn recognized_feature_header(field: &MetadataField, format: Format) -> bool {
    if field.continuation {
        return false;
    }
    match format {
        Format::GenBank => field.value == "Location/Qualifiers",
        Format::Embl => {
            field.value.is_empty()
                || field.value.split_whitespace().collect::<Vec<_>>()
                    == ["Key", "Location/Qualifiers"]
        }
    }
}

fn is_description(format: Format, key: &str) -> bool {
    key == description_key(format)
}

fn description_key(format: Format) -> &'static str {
    match format {
        Format::GenBank => "DEFINITION",
        Format::Embl => "DE",
    }
}

fn is_sequence_header(format: Format, key: &str) -> bool {
    match format {
        Format::GenBank => key == "ORIGIN",
        Format::Embl => key == "SQ",
    }
}

fn is_local_date(value: &str) -> bool {
    value
        .strip_suffix(" (BioFerris local derivation)")
        .is_some_and(|date| validate_date(date).is_ok())
}

fn validate_date(date: &str) -> Result<(), BioIoError> {
    let invalid = || {
        transform_error(
            IoErrorKind::InvalidHeader,
            "local derivation date must be a real DD-MMM-YYYY Gregorian date with year 0001..=9999",
        )
    };
    if date.len() != 11
        || !date.is_ascii()
        || date.as_bytes()[2] != b'-'
        || date.as_bytes()[6] != b'-'
        || !date.as_bytes()[..2].iter().all(u8::is_ascii_digit)
        || !date.as_bytes()[7..].iter().all(u8::is_ascii_digit)
    {
        return Err(invalid());
    }
    let day = date[..2].parse::<u16>().map_err(|_| invalid())?;
    let year = date[7..].parse::<u16>().map_err(|_| invalid())?;
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = match &date[3..6] {
        "JAN" | "MAR" | "MAY" | "JUL" | "AUG" | "OCT" | "DEC" => 31,
        "APR" | "JUN" | "SEP" | "NOV" => 30,
        "FEB" if leap => 29,
        "FEB" => 28,
        _ => return Err(invalid()),
    };
    if year == 0 || day == 0 || day > days {
        return Err(invalid());
    }
    Ok(())
}

fn require_discard(policy: AnnotationLossPolicy, reason: &str) -> Result<(), BioIoError> {
    if policy == AnnotationLossPolicy::Reject {
        return Err(transform_error(
            IoErrorKind::DataLoss,
            format!("{reason}; an explicit annotation discard policy is required"),
        ));
    }
    Ok(())
}

fn checked_sub(value: usize, subtract: usize) -> Result<usize, BioIoError> {
    value.checked_sub(subtract).ok_or_else(coordinate_error)
}

fn coordinate_error() -> BioIoError {
    transform_error(
        IoErrorKind::LengthMismatch,
        "feature coordinate arithmetic overflow",
    )
}

fn layout_error() -> BioIoError {
    transform_error(
        IoErrorKind::InvalidStructure,
        "invalid original record layout",
    )
}

fn transform_error(kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, "annotated-transform", message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_identity_and_dates_are_checked_without_a_record() {
        for id in [
            "",
            "LOCAL_",
            "OTHER_name",
            "LOCAL_has-dash",
            "LOCAL_é",
            "LOCAL_12345678901",
        ] {
            assert_eq!(
                AnnotationTransformOptions::new(id).unwrap_err().kind(),
                IoErrorKind::InvalidHeader
            );
        }
        for id in ["LOCAL_a", "LOCAL_1234567890", "LOCAL__"] {
            assert_eq!(AnnotationTransformOptions::new(id).unwrap().local_id(), id);
        }
        for date in ["01-JAN-0001", "29-FEB-2000", "29-FEB-2024", "31-DEC-9999"] {
            assert_eq!(
                AnnotationTransformOptions::new("LOCAL_a")
                    .unwrap()
                    .with_date(date)
                    .unwrap()
                    .date(),
                Some(date)
            );
        }
        for date in [
            "1-JAN-2026",
            "01-jan-2026",
            "00-JAN-2026",
            "01-JAN-0000",
            "29-FEB-1900",
            "29-FEB-2025",
            "31-APR-2026",
            "01-XYZ-2026",
            "01-JAN-10000",
        ] {
            assert_eq!(
                validate_date(date).unwrap_err().kind(),
                IoErrorKind::InvalidHeader
            );
        }
    }

    #[test]
    fn joined_mapping_preserves_written_order_and_normalizes_complements() {
        let source = Location::parse("complement(join(2..4,8..9))").unwrap();
        let mapped = map_location(
            source.expression(),
            &AnnotationOperation::ReverseComplement,
            12,
        )
        .unwrap();
        assert_eq!(
            writer::render_location(&mapped, Format::GenBank).unwrap(),
            "complement(join(complement(9..11),complement(4..5)))"
        );
        let reversed = map_location(&mapped, &AnnotationOperation::ReverseComplement, 12).unwrap();
        assert_eq!(
            writer::render_location(&reversed, Format::GenBank).unwrap(),
            "complement(join(2..4,8..9))"
        );
        assert_eq!(relation(source.expression(), &(1..9)), Relation::Contained);
        assert_eq!(relation(source.expression(), &(4..7)), Relation::Disjoint);
        assert_eq!(relation(source.expression(), &(1..7)), Relation::Partial);
    }

    fn options(policy: AnnotationLossPolicy) -> AnnotationTransformOptions {
        AnnotationTransformOptions::new("LOCAL_checked")
            .unwrap()
            .with_date("16-OCT-2026")
            .unwrap()
            .with_loss_policy(policy)
    }

    #[test]
    fn opaque_headers_are_loss_checked_before_regeneration() {
        for (locus_tail, origin) in [(" opaque", "ORIGIN"), ("", "ORIGIN      opaque")] {
            let text = format!(
                "LOCUS       ORIGINAL 4 bp DNA linear SYN 15-OCT-2026{locus_tail}\n{origin}\n        1 aCGt\n//\n"
            );
            let source = super::super::parse_genbank::<Dna>(&text).unwrap().remove(0);
            assert_eq!(
                source
                    .slice(0..4, options(AnnotationLossPolicy::Reject))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::DataLoss
            );
            let result = source
                .slice(0..4, options(AnnotationLossPolicy::Discard))
                .unwrap();
            assert!(!result.record().raw_record().contains("opaque"));
            assert!(result.report().changes.iter().any(|change| matches!(
                change,
                AnnotationChange::MetadataChanged { index: Some(_), before: Some(before), after: Some(after), .. }
                    if before.contains("opaque") && !after.contains("opaque")
            )));
        }
    }

    #[test]
    fn opaque_sq_requires_explicit_loss_policy() {
        let text = "ID   ORIGINAL; SV 1; linear; genomic DNA; STD; SYN; 4 BP.\nSQ   Sequence 4 BP; MD5 opaque;\n     aCGt 4\n//\n";
        let source = super::super::parse_embl::<Dna>(text).unwrap().remove(0);
        assert_eq!(
            source
                .slice(0..4, options(AnnotationLossPolicy::Reject))
                .unwrap_err()
                .kind(),
            IoErrorKind::DataLoss
        );
        let result = source
            .slice(0..4, options(AnnotationLossPolicy::Discard))
            .unwrap();
        assert!(
            result
                .record()
                .raw_record()
                .contains("SQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;")
        );
        assert!(!result.record().raw_record().contains("MD5 opaque"));
    }

    #[test]
    fn protein_complements_are_unsupported_before_discard_classification() {
        let text = "LOCUS       ORIGINAL 4 aa linear SYN 15-OCT-2026\nFEATURES             Location/Qualifiers\n     misc_feature    complement(2..3)\nORIGIN\n        1 MELK\n//\n";
        let source = super::super::parse_genbank::<crate::Protein>(text)
            .unwrap()
            .remove(0);
        for policy in [AnnotationLossPolicy::Reject, AnnotationLossPolicy::Discard] {
            for range in [0..4, 0..1, 2..2] {
                assert_eq!(
                    source.slice(range, options(policy)).unwrap_err().kind(),
                    IoErrorKind::Unsupported
                );
            }
        }
    }
}
