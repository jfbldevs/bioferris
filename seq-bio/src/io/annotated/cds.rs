//! Annotation-aware complete-CDS translation and optional reference validation.

use std::fmt;

use align_core::Alphabet;

use super::{BioIoError, Feature};
use crate::{
    ConfiguredTranslationExt, GeneticCode, PartialCodonPolicy, ProteinSequence, ReadingFrame,
    Sequence, Translation, TranslationConfig, TranslationError, TranslationMode,
};

/// How an annotated `/translation` is checked against the computed protein.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CdsReferencePolicy {
    /// Check exact equality when present; permit an absent reference (default).
    #[default]
    CheckIfPresent,
    /// Require a reference and exact equality after ASCII whitespace removal.
    RequireMatch,
    /// Do not interpret `/translation`; biological CDS checks still apply.
    Ignore,
}

/// Reference validation performed for a successful CDS translation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CdsReferenceStatus {
    /// No `/translation` was present under `CheckIfPresent`.
    Absent,
    /// The normalized reference matched every computed residue and its length.
    Matched,
    /// The caller explicitly disabled reference interpretation.
    Ignored,
}

/// Annotation-aware translation policies; no annotation is silently overridden.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CdsTranslationOptions {
    /// Table used only when `/transl_table` is absent (default NCBI table 1).
    pub default_genetic_code: GeneticCode,
    /// Whether to compare the source `/translation` (default `CheckIfPresent`).
    pub reference_policy: CdsReferencePolicy,
}

/// An extracted complete-CDS product and its annotation-validation metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CdsTranslation {
    translation: Translation,
    feature_length: usize,
    codon_start: u8,
    ribosomal_slippage: bool,
    reference_status: CdsReferenceStatus,
    reference_translation: Option<String>,
}

impl CdsTranslation {
    /// Returns the checked product; its initial residue is M and stop is omitted.
    #[must_use]
    pub fn peptide(&self) -> &str {
        self.translation.as_str()
    }

    /// Returns the shared translation result and selected code/frame metadata.
    #[must_use]
    pub const fn translation(&self) -> &Translation {
        &self.translation
    }

    /// Returns extracted bases before the codon-start offset, including stop.
    #[must_use]
    pub const fn feature_length(&self) -> usize {
        self.feature_length
    }

    /// Returns the one-based `/codon_start` value, defaulting to 1.
    #[must_use]
    pub const fn codon_start(&self) -> u8 {
        self.codon_start
    }

    /// Returns the annotation-selected genetic code or explicit fallback.
    #[must_use]
    pub const fn genetic_code(&self) -> GeneticCode {
        self.translation.genetic_code()
    }

    /// Returns whether the CDS declares `/ribosomal_slippage`. Its location was
    /// translated exactly as annotated, including any frameshift it encodes.
    #[must_use]
    pub const fn ribosomal_slippage(&self) -> bool {
        self.ribosomal_slippage
    }

    /// Returns whether a source translation was absent, matched or ignored.
    #[must_use]
    pub const fn reference_status(&self) -> CdsReferenceStatus {
        self.reference_status
    }

    /// Returns a normalized, validated matching reference, when checked/present.
    #[must_use]
    pub fn reference_translation(&self) -> Option<&str> {
        self.reference_translation.as_deref()
    }
}

/// A failed interpretation, extraction, translation or exact reference check.
#[derive(Debug)]
#[non_exhaustive]
pub enum CdsTranslationError {
    /// The feature key was not the case-sensitive INSDC key `CDS`.
    NotCds {
        /// Actual feature key.
        key: String,
    },
    /// A qualifier requires biological semantics outside this translation profile.
    UnsupportedQualifier {
        /// Qualifier name without the initial slash.
        name: String,
    },
    /// An interpretation qualifier occurred more than once.
    DuplicateQualifier {
        /// Repeated qualifier name.
        name: String,
    },
    /// A qualifier was a flag, invalid number or invalid reference protein.
    InvalidQualifier {
        /// Invalid qualifier name.
        name: String,
        /// Explanation of the rejected value.
        reason: String,
    },
    /// `RequireMatch` was selected but `/translation` was absent.
    MissingReference,
    /// The normalized source protein disagreed with the computed product.
    ReferenceMismatch {
        /// Normalized source `/translation`.
        reference: String,
        /// Computed product, with no terminal stop.
        computed: String,
        /// Zero-based first differing residue, or shorter length for a prefix.
        first_difference: usize,
    },
    /// The location could not be extracted without inventing sequence semantics.
    Extraction(BioIoError),
    /// The selected table or complete-CDS translation failed validation.
    ///
    /// Any nucleotide error offset refers to the extracted coding-oriented
    /// feature, not the original genomic sequence.
    Translation(TranslationError),
}

impl CdsTranslationError {
    /// Returns a stable category; inspect typed variants for structured details.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NotCds { .. } => "invalid_cds_feature",
            Self::UnsupportedQualifier { .. } => "unsupported_cds_qualifier",
            Self::DuplicateQualifier { .. } | Self::InvalidQualifier { .. } => {
                "invalid_cds_qualifier"
            }
            Self::MissingReference => "missing_cds_translation",
            Self::ReferenceMismatch { .. } => "cds_translation_mismatch",
            Self::Extraction(error) => error.code(),
            Self::Translation(TranslationError::UnsupportedGeneticCode { .. }) => {
                "unsupported_genetic_code"
            }
            Self::Translation(_) => "cds_translation_error",
        }
    }
}

impl fmt::Display for CdsTranslationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotCds { key } => write!(formatter, "expected CDS feature, got {key}"),
            Self::UnsupportedQualifier { name } => {
                write!(formatter, "CDS translation does not support /{name}")
            }
            Self::DuplicateQualifier { name } => {
                write!(formatter, "CDS translation requires at most one /{name}")
            }
            Self::InvalidQualifier { name, reason } => {
                write!(formatter, "invalid CDS /{name}: {reason}")
            }
            Self::MissingReference => {
                formatter.write_str("CDS reference policy requires /translation")
            }
            Self::ReferenceMismatch {
                reference,
                computed,
                first_difference,
            } => write!(
                formatter,
                "CDS /translation differs at residue {first_difference}: reference {reference}, computed {computed}"
            ),
            Self::Extraction(error) => write!(formatter, "CDS extraction failed: {error}"),
            Self::Translation(error) => write!(formatter, "CDS translation failed: {error}"),
        }
    }
}

impl std::error::Error for CdsTranslationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Extraction(error) => Some(error),
            Self::Translation(error) => Some(error),
            _ => None,
        }
    }
}

impl Feature {
    /// Extracts and translates a complete CDS, checking `/translation` if present.
    ///
    /// Exact local joins/complements select biological coding orientation before
    /// `/codon_start` is applied. `/transl_table` takes precedence over the default
    /// NCBI table 1. Reference wrapping whitespace is ignored; residues are not
    /// corrected. The parent feature and sequence are immutable.
    ///
    /// # Errors
    /// Rejects non-CDS features, unsupported exceptions/pseudogenes, invalid or
    /// duplicate interpretation qualifiers, unsupported extraction/tables,
    /// incomplete CDS products and a mismatching source `/translation`.
    /// See [`Self::translate_cds_with`] for explicit reference/fallback policies.
    pub fn translate_cds<A: Alphabet>(
        &self,
        parent: &Sequence<A>,
    ) -> Result<CdsTranslation, CdsTranslationError>
    where
        Sequence<A>: ConfiguredTranslationExt,
    {
        self.translate_cds_with(parent, CdsTranslationOptions::default())
    }

    /// Translates this CDS using annotation semantics and explicit validation.
    ///
    /// `/transl_table` always takes precedence over `default_genetic_code`.
    /// `reference_policy` controls only the source protein comparison: `Ignore`
    /// does not disable complete-CDS, location or biological exception checks.
    /// Translation error offsets use the extracted coding-oriented sequence.
    ///
    /// # Errors
    /// Returns [`CdsTranslationError`] for the failures documented by
    /// [`Self::translate_cds`], or a missing reference under `RequireMatch`.
    pub fn translate_cds_with<A: Alphabet>(
        &self,
        parent: &Sequence<A>,
        options: CdsTranslationOptions,
    ) -> Result<CdsTranslation, CdsTranslationError>
    where
        Sequence<A>: ConfiguredTranslationExt,
    {
        if self.key != "CDS" {
            return Err(CdsTranslationError::NotCds {
                key: self.key.clone(),
            });
        }
        for qualifier in &self.qualifiers {
            if matches!(
                qualifier.name.as_str(),
                "pseudo" | "pseudogene" | "transl_except" | "exception" | "trans_splicing"
            ) {
                return Err(CdsTranslationError::UnsupportedQualifier {
                    name: qualifier.name.clone(),
                });
            }
        }
        // A frameshift lives in the location, usually as a join whose parts repeat
        // or skip bases; the location is translated exactly as annotated, and the
        // complete-CDS and reference checks reject a product read without its shift.
        let ribosomal_slippage = flag_qualifier(self, "ribosomal_slippage")?;
        let codon_start = decimal_qualifier(self, "codon_start")?.unwrap_or(1);
        let reading_frame = match codon_start {
            1 => ReadingFrame::Zero,
            2 => ReadingFrame::One,
            3 => ReadingFrame::Two,
            _ => return Err(invalid("codon_start", "expected 1, 2 or 3")),
        };
        let genetic_code = match decimal_qualifier(self, "transl_table")? {
            Some(id) => GeneticCode::try_from(id).map_err(CdsTranslationError::Translation)?,
            None => options.default_genetic_code,
        };
        let coding = self
            .extract(parent)
            .map_err(CdsTranslationError::Extraction)?;
        let translation = coding
            .translate_with(TranslationConfig {
                genetic_code,
                reading_frame,
                mode: TranslationMode::Cds,
                partial_codons: PartialCodonPolicy::Reject,
            })
            .map_err(CdsTranslationError::Translation)?;
        let (reference_status, reference_translation) =
            compare_reference(self, translation.as_str(), options.reference_policy)?;
        Ok(CdsTranslation {
            translation,
            feature_length: coding.len(),
            codon_start: match reading_frame {
                ReadingFrame::Zero => 1,
                ReadingFrame::One => 2,
                ReadingFrame::Two => 3,
            },
            ribosomal_slippage,
            reference_status,
            reference_translation,
        })
    }
}

fn invalid(name: &str, reason: impl Into<String>) -> CdsTranslationError {
    CdsTranslationError::InvalidQualifier {
        name: name.to_owned(),
        reason: reason.into(),
    }
}

fn valued_qualifier<'a>(
    feature: &'a Feature,
    name: &str,
) -> Result<Option<&'a str>, CdsTranslationError> {
    let mut qualifiers = feature.qualifiers.iter().filter(|q| q.name == name);
    let Some(qualifier) = qualifiers.next() else {
        return Ok(None);
    };
    if qualifiers.next().is_some() {
        return Err(CdsTranslationError::DuplicateQualifier {
            name: name.to_owned(),
        });
    }
    qualifier
        .value
        .as_deref()
        .map(Some)
        .ok_or_else(|| invalid(name, "a value is required"))
}

/// Whether `/name` is present, as a flag without a value and at most once.
fn flag_qualifier(feature: &Feature, name: &str) -> Result<bool, CdsTranslationError> {
    let mut qualifiers = feature.qualifiers.iter().filter(|q| q.name == name);
    let Some(qualifier) = qualifiers.next() else {
        return Ok(false);
    };
    if qualifiers.next().is_some() {
        return Err(CdsTranslationError::DuplicateQualifier {
            name: name.to_owned(),
        });
    }
    if qualifier.value.is_some() {
        return Err(invalid(name, "expected a flag without a value"));
    }
    Ok(true)
}

fn decimal_qualifier(feature: &Feature, name: &str) -> Result<Option<u16>, CdsTranslationError> {
    valued_qualifier(feature, name)?
        .map(|value| {
            let value = value.trim_matches(|c: char| c.is_ascii_whitespace());
            if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_digit()) {
                return Err(invalid(name, "expected one ASCII decimal integer"));
            }
            value
                .parse::<u16>()
                .map_err(|_| invalid(name, "integer exceeds 65535"))
        })
        .transpose()
}

fn compare_reference(
    feature: &Feature,
    computed: &str,
    policy: CdsReferencePolicy,
) -> Result<(CdsReferenceStatus, Option<String>), CdsTranslationError> {
    if policy == CdsReferencePolicy::Ignore {
        return Ok((CdsReferenceStatus::Ignored, None));
    }
    let Some(value) = valued_qualifier(feature, "translation")? else {
        return if policy == CdsReferencePolicy::RequireMatch {
            Err(CdsTranslationError::MissingReference)
        } else {
            Ok((CdsReferenceStatus::Absent, None))
        };
    };
    let reference: String = value.chars().filter(|c| !c.is_ascii_whitespace()).collect();
    if reference.is_empty() || !reference.bytes().all(|byte| byte.is_ascii_uppercase()) {
        return Err(invalid(
            "translation",
            "expected a nonempty uppercase stop-free protein",
        ));
    }
    ProteinSequence::from_ascii(&reference)
        .map_err(|error| invalid("translation", error.to_string()))?;
    if reference != computed {
        let first_difference = reference
            .bytes()
            .zip(computed.bytes())
            .position(|(left, right)| left != right)
            .unwrap_or_else(|| reference.len().min(computed.len()));
        return Err(CdsTranslationError::ReferenceMismatch {
            reference,
            computed: computed.to_owned(),
            first_difference,
        });
    }
    Ok((CdsReferenceStatus::Matched, Some(reference)))
}
