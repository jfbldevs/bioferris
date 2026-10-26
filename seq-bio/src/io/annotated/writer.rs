//! Model-based, preflight-checked writers for the supported annotated profiles.

use std::io::Write;
use std::path::Path;

use align_core::Alphabet;

use super::{
    AnnotatedRecord, BioIoError, Feature, Format, IoErrorKind, LocationExpression, MetadataField,
    Qualifier, ReadLimits, RecordPart, is_component_name, parse_document,
};
use crate::io::io_error;
use crate::publication::{PublishError, publish_new};

/// Writes one same-format `GenBank` record from its supported model.
///
/// Sequence rows and length counters are regenerated; source case, ordered
/// metadata/features and decoded qualifier values are preserved. Raw document
/// caches are not copied. The complete bounded output is reparsed and compared
/// before the first write. This is semantic preservation, not byte identity or
/// certification of opaque location expressions.
/// The caller retains responsibility for flushing its writer.
///
/// # Errors
/// Rejects cross-format conversion, unrepresentable model values and default
/// reader/output limits before writing. An I/O failure can leave partial stream
/// output; use [`write_genbank_new`] for create-only named publication.
pub fn write_genbank<A: Alphabet>(
    mut writer: impl Write,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    let output = render(record, Format::GenBank)?;
    writer
        .write_all(output.as_bytes())
        .map_err(|error| io_error("genbank", error))
}

/// Writes one same-format modern EMBL record with checked semantic preservation.
///
/// Regenerates sequence rows and standard SQ counts, preserving literal source
/// case. U is counted as `other`, never substituted for T. SQ input must contain
/// only the supported length-only or standard A/C/G/T/other count grammar.
/// Metadata, features and decoded qualifier values retain their model order.
/// The caller retains responsibility for flushing its writer.
///
/// # Errors
/// Rejects cross-format conversion, opaque SQ tails, unrepresentable model
/// values and default limits before writing. I/O failures may leave partial
/// stream output; [`write_embl_new`] publishes only a complete new named file.
pub fn write_embl<A: Alphabet>(
    mut writer: impl Write,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    let output = render(record, Format::Embl)?;
    writer
        .write_all(output.as_bytes())
        .map_err(|error| io_error("embl", error))
}

/// Preflights and publishes a new same-format `GenBank` file without overwriting.
///
/// The destination directory must be trusted. Publication uses a temporary
/// file in that directory and a create-only hard link, with no copy fallback.
///
/// # Errors
/// Returns the errors of [`write_genbank`], existing-destination errors and
/// filesystem/publication failures. Failed publication leaves prior output
/// unchanged; temporary output is removed on ordinary error paths.
pub fn write_genbank_new<A: Alphabet>(
    path: impl AsRef<Path>,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    write_new(path.as_ref(), record, Format::GenBank)
}

/// Preflights and publishes a new same-format modern EMBL file without replacing it.
///
/// The destination directory must be trusted; there is no overwrite or copy
/// fallback. A complete document is validated before temporary output is opened.
///
/// # Errors
/// Returns the errors of [`write_embl`], existing-destination errors and
/// filesystem/publication failures. Prior destinations are preserved.
pub fn write_embl_new<A: Alphabet>(
    path: impl AsRef<Path>,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    write_new(path.as_ref(), record, Format::Embl)
}

fn write_new<A: Alphabet>(
    path: &Path,
    record: &AnnotatedRecord<A>,
    format: Format,
) -> Result<(), BioIoError> {
    let output = render(record, format)?;
    publish_new(path, |writer| writer.write_all(output.as_bytes())).map_err(|error| match error {
        PublishError::Io(error) | PublishError::Write(error) => io_error(format.name(), error),
        PublishError::DestinationExists(path) => BioIoError::new(
            IoErrorKind::Io,
            format.name(),
            format!(
                "destination already exists; output is create-only: {}",
                path.display()
            ),
        ),
        PublishError::TemporaryUnavailable => BioIoError::new(
            IoErrorKind::Io,
            format.name(),
            "cannot reserve a temporary output path",
        ),
    })
}

fn error(format: Format, kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, format.name(), message)
}

fn unsupported(format: Format, message: impl Into<String>) -> BioIoError {
    error(format, IoErrorKind::Unsupported, message)
}

struct Output {
    text: String,
    limits: ReadLimits,
    format: Format,
}

impl Output {
    fn line(&mut self, parts: &[&str]) -> Result<(), BioIoError> {
        let mut bytes = 1_usize;
        for part in parts {
            if part.contains(['\r', '\n']) {
                return Err(unsupported(
                    self.format,
                    "embedded line ending in a physical field",
                ));
            }
            bytes = bytes.checked_add(part.len()).ok_or_else(|| {
                error(
                    self.format,
                    IoErrorKind::LimitExceeded,
                    "output byte count overflow",
                )
            })?;
        }
        if bytes > self.limits.max_line_bytes {
            return Err(error(
                self.format,
                IoErrorKind::LimitExceeded,
                "physical output line byte limit exceeded",
            ));
        }
        if bytes > self.limits.max_total_bytes.saturating_sub(self.text.len()) {
            return Err(error(
                self.format,
                IoErrorKind::LimitExceeded,
                "annotated output byte limit exceeded",
            ));
        }
        self.text.try_reserve(bytes).map_err(|failure| {
            error(self.format, IoErrorKind::LimitExceeded, failure.to_string())
        })?;
        for part in parts {
            self.text.push_str(part);
        }
        self.text.push('\n');
        Ok(())
    }
}

pub(super) fn render<A: Alphabet>(
    record: &AnnotatedRecord<A>,
    format: Format,
) -> Result<String, BioIoError> {
    if record.format != format {
        return Err(unsupported(
            format,
            "annotated writing requires the record's source format",
        ));
    }
    let limits = ReadLimits::default();
    if record.record.sequence().len() > limits.max_record_symbols {
        return Err(error(
            format,
            IoErrorKind::LimitExceeded,
            "sequence exceeds the record symbol limit",
        ));
    }
    if record.record.qualities().is_some() {
        return Err(error(
            format,
            IoErrorKind::DataLoss,
            "annotated output cannot preserve quality scores",
        ));
    }
    let mut output = Output {
        text: String::new(),
        limits,
        format,
    };
    let mut metadata_count = 0;
    let mut feature_count = 0;
    for part in &record.layout {
        match *part {
            RecordPart::Metadata(index) if index == metadata_count => {
                let field = record
                    .metadata
                    .get(index)
                    .ok_or_else(|| unsupported(format, "metadata layout index is out of bounds"))?;
                write_metadata(&mut output, field, record)?;
                metadata_count += 1;
            }
            RecordPart::Feature(index) if index == feature_count => {
                let feature = record
                    .features
                    .get(index)
                    .ok_or_else(|| unsupported(format, "feature layout index is out of bounds"))?;
                write_feature(&mut output, feature)?;
                feature_count += 1;
            }
            _ => {
                return Err(unsupported(
                    format,
                    "record layout does not preserve source order",
                ));
            }
        }
    }
    if metadata_count != record.metadata.len() || feature_count != record.features.len() {
        return Err(error(
            format,
            IoErrorKind::DataLoss,
            "record layout omits metadata or features",
        ));
    }
    write_sequence(&mut output, record.record.source_sequence())?;
    output.line(&["//"])?;
    let reparsed = parse_document::<A>(&output.text, format).map_err(|failure| {
        error(
            format,
            IoErrorKind::DataLoss,
            format!("rendered document cannot preserve its model: {failure}"),
        )
    })?;
    let [parsed] = reparsed.as_slice() else {
        return Err(error(
            format,
            IoErrorKind::DataLoss,
            "rendering did not produce exactly one record",
        ));
    };
    verify_semantics(record, parsed, format)?;
    Ok(output.text)
}

fn write_metadata<A: Alphabet>(
    output: &mut Output,
    field: &MetadataField,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    let value = metadata_value(field, record, output.format)?;
    let width = if output.format == Format::GenBank {
        12
    } else {
        5
    };
    if field.key.is_empty() || field.key.len() > width || field.key.contains(['\r', '\n']) {
        return Err(unsupported(
            output.format,
            "unrepresentable metadata field name",
        ));
    }
    let prefix = if field.continuation {
        " ".repeat(width)
    } else {
        let indent = if output.format == Format::GenBank {
            field
                .raw_line
                .get(..width)
                .filter(|prefix| prefix.trim() == field.key)
                .map_or(0, |prefix| {
                    prefix.bytes().take_while(u8::is_ascii_whitespace).count()
                })
        } else {
            0
        };
        if indent + field.key.len() > width {
            return Err(unsupported(
                output.format,
                "metadata indentation exceeds its field width",
            ));
        }
        format!(
            "{}{}{}",
            " ".repeat(indent),
            field.key,
            " ".repeat(width - indent - field.key.len())
        )
    };
    if !field.continuation && value.is_empty() {
        output.line(&[prefix.trim_end()])
    } else if output.format == Format::GenBank
        && !field.continuation
        && field.key == "FEATURES"
        && value == "Location/Qualifiers"
    {
        output.line(&["FEATURES             ", &value])
    } else {
        output.line(&[&prefix, &value])
    }
}

fn metadata_value<A: Alphabet>(
    field: &MetadataField,
    record: &AnnotatedRecord<A>,
    format: Format,
) -> Result<String, BioIoError> {
    if field.value.len() > ReadLimits::default().max_line_bytes {
        return Err(error(
            format,
            IoErrorKind::LimitExceeded,
            "metadata exceeds physical line byte limit",
        ));
    }
    match (format, field.key.as_str(), field.continuation) {
        (Format::GenBank, "LOCUS", false) => {
            replace_token(&field.value, 1, record.record.sequence().len(), format)
        }
        (Format::Embl, "ID", false) => {
            let (prefix, tail) = field
                .value
                .rsplit_once(';')
                .ok_or_else(|| unsupported(format, "invalid ID payload"))?;
            Ok(format!(
                "{prefix};{}",
                replace_token(tail, 0, record.record.sequence().len(), format)?
            ))
        }
        (Format::Embl, "SQ", false) => {
            sq_value(&field.value, record.record.source_sequence(), format)
        }
        _ => Ok(field.value.clone()),
    }
}

fn replace_token(
    value: &str,
    ordinal: usize,
    replacement: usize,
    format: Format,
) -> Result<String, BioIoError> {
    let mut token_start = None;
    let mut token_number = 0;
    for (offset, character) in value
        .char_indices()
        .chain(std::iter::once((value.len(), ' ')))
    {
        if character.is_whitespace() {
            if let Some(start) = token_start.take() {
                if token_number == ordinal {
                    if value[start..offset].parse::<usize>().is_err() {
                        return Err(unsupported(
                            format,
                            "header length is not an unsigned integer",
                        ));
                    }
                    return Ok(format!(
                        "{}{replacement}{}",
                        &value[..start],
                        &value[offset..]
                    ));
                }
                token_number += 1;
            }
        } else if token_start.is_none() {
            token_start = Some(offset);
        }
    }
    Err(unsupported(format, "header lacks its sequence length"))
}

pub(super) fn sq_value(value: &str, source: &str, format: Format) -> Result<String, BioIoError> {
    let fields: Vec<_> = value.split_whitespace().collect();
    let length_only = fields.len() == 3;
    let standard_counts = fields.len() == 13
        && ["A;", "C;", "G;", "T;", "other;"]
            .iter()
            .enumerate()
            .all(|(index, label)| {
                fields[4 + index * 2] == *label && fields[3 + index * 2].parse::<usize>().is_ok()
            });
    if !(length_only || standard_counts)
        || fields.first() != Some(&"Sequence")
        || fields
            .get(1)
            .is_none_or(|length| length.parse::<usize>().is_err())
        || fields.get(2) != Some(&"BP;")
    {
        return Err(unsupported(
            format,
            "SQ must contain only length or standard A/C/G/T/other counts; opaque tails cannot be discarded",
        ));
    }
    let mut counts = [0_usize; 5];
    for symbol in source.bytes() {
        let index = match symbol.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        counts[index] += 1;
    }
    Ok(format!(
        "Sequence {} BP; {} A; {} C; {} G; {} T; {} other;",
        source.len(),
        counts[0],
        counts[1],
        counts[2],
        counts[3],
        counts[4]
    ))
}

fn write_feature(output: &mut Output, feature: &Feature) -> Result<(), BioIoError> {
    if !is_component_name(&feature.key, 15) {
        return Err(unsupported(output.format, "unrepresentable feature key"));
    }
    let location = render_location(&feature.location.expression, output.format)?;
    let start = if output.format == Format::GenBank {
        "     "
    } else {
        "FT   "
    };
    let prefix = format!("{start}{:15} ", feature.key);
    output.line(&[&prefix, &location])?;
    let continuation = format!("{start}                ");
    for qualifier in &feature.qualifiers {
        write_qualifier(output, &continuation, qualifier)?;
    }
    Ok(())
}

fn push_location(output: &mut String, value: &str, format: Format) -> Result<(), BioIoError> {
    let maximum = ReadLimits::default().max_line_bytes - 22;
    if value.len() > maximum.saturating_sub(output.len()) {
        return Err(error(
            format,
            IoErrorKind::LimitExceeded,
            "rendered location exceeds physical line limit",
        ));
    }
    output
        .try_reserve(value.len())
        .map_err(|failure| error(format, IoErrorKind::LimitExceeded, failure.to_string()))?;
    output.push_str(value);
    Ok(())
}

pub(super) fn render_location(
    expression: &LocationExpression,
    format: Format,
) -> Result<String, BioIoError> {
    let mut output = String::new();
    write_location(expression, &mut output, 0, format)?;
    Ok(output)
}

fn write_location(
    expression: &LocationExpression,
    output: &mut String,
    depth: usize,
    format: Format,
) -> Result<(), BioIoError> {
    if depth > super::MAX_LOCATION_DEPTH {
        return Err(error(
            format,
            IoErrorKind::LimitExceeded,
            "location nesting exceeds 64 levels",
        ));
    }
    match expression {
        LocationExpression::Range(range) => {
            if range.start >= range.end {
                return Err(unsupported(
                    format,
                    "empty or reversed feature ranges are not representable",
                ));
            }
            let start = range.start + 1;
            let text = if start == range.end {
                start.to_string()
            } else {
                format!("{start}..{}", range.end)
            };
            push_location(output, &text, format)
        }
        LocationExpression::Join(parts) => {
            if parts.len() < 2 {
                return Err(unsupported(format, "join requires at least two locations"));
            }
            push_location(output, "join(", format)?;
            for (index, part) in parts.iter().enumerate() {
                if index > 0 {
                    push_location(output, ",", format)?;
                }
                write_location(part, output, depth + 1, format)?;
            }
            push_location(output, ")", format)
        }
        LocationExpression::Complement(inner) => {
            push_location(output, "complement(", format)?;
            write_location(inner, output, depth + 1, format)?;
            push_location(output, ")", format)
        }
        LocationExpression::Unsupported { expression, .. } => {
            push_location(output, expression, format)
        }
    }
}

fn escaped(value: &str, format: Format) -> Result<String, BioIoError> {
    let additional = value.bytes().filter(|byte| *byte == b'"').count();
    let bytes = value.len().checked_add(additional).ok_or_else(|| {
        error(
            format,
            IoErrorKind::LimitExceeded,
            "qualifier length overflow",
        )
    })?;
    if bytes > ReadLimits::default().max_line_bytes {
        return Err(error(
            format,
            IoErrorKind::LimitExceeded,
            "escaped qualifier exceeds physical line limit",
        ));
    }
    let mut result = String::new();
    result
        .try_reserve(bytes)
        .map_err(|failure| error(format, IoErrorKind::LimitExceeded, failure.to_string()))?;
    for character in value.chars() {
        result.push(character);
        if character == '"' {
            result.push('"');
        }
    }
    Ok(result)
}

fn write_qualifier(
    output: &mut Output,
    prefix: &str,
    qualifier: &Qualifier,
) -> Result<(), BioIoError> {
    if !is_component_name(&qualifier.name, 20) {
        return Err(unsupported(output.format, "unrepresentable qualifier name"));
    }
    let Some(value) = &qualifier.value else {
        return output.line(&[prefix, "/", &qualifier.name]);
    };
    let quoted = qualifier
        .raw
        .split_once('=')
        .is_some_and(|(_, payload)| payload.starts_with('"'));
    let mut lines = value.split('\n').peekable();
    let mut first = true;
    while let Some(line) = lines.next() {
        let content = if quoted {
            escaped(line, output.format)?
        } else {
            line.to_owned()
        };
        let opening = if first {
            format!("/{}={}", qualifier.name, if quoted { "\"" } else { "" })
        } else {
            String::new()
        };
        let closing = if quoted && lines.peek().is_none() {
            "\""
        } else {
            ""
        };
        output.line(&[prefix, &opening, &content, closing])?;
        first = false;
    }
    Ok(())
}

fn write_sequence(output: &mut Output, source: &str) -> Result<(), BioIoError> {
    if !source.is_ascii() {
        return Err(error(
            output.format,
            IoErrorKind::InvalidSequence,
            "sequence must contain ASCII symbols",
        ));
    }
    for offset in (0..source.len()).step_by(60) {
        let end = (offset + 60).min(source.len());
        let mut groups = String::new();
        for start in (offset..end).step_by(10) {
            if start != offset {
                groups.push(' ');
            }
            groups.push_str(&source[start..(start + 10).min(end)]);
        }
        match output.format {
            Format::GenBank => output.line(&[&format!("{:9} ", offset + 1), &groups])?,
            Format::Embl => output.line(&[&format!("     {groups:65} {end:9}")])?,
        }
    }
    Ok(())
}

fn verify_semantics<A: Alphabet>(
    original: &AnnotatedRecord<A>,
    parsed: &AnnotatedRecord<A>,
    format: Format,
) -> Result<(), BioIoError> {
    let same_record = original.record == parsed.record
        && original.format == parsed.format
        && original.locus == parsed.locus
        && original.molecule == parsed.molecule
        && original.accessions == parsed.accessions
        && original.layout == parsed.layout;
    let same_features = original.features.len() == parsed.features.len()
        && original
            .features
            .iter()
            .zip(&parsed.features)
            .all(|(left, right)| {
                left.key == right.key
                    && left.location.expression == right.location.expression
                    && left.qualifiers.len() == right.qualifiers.len()
                    && left
                        .qualifiers
                        .iter()
                        .zip(&right.qualifiers)
                        .all(|(left, right)| left.name == right.name && left.value == right.value)
            });
    if !same_record || !same_features || original.metadata.len() != parsed.metadata.len() {
        return Err(error(
            format,
            IoErrorKind::DataLoss,
            "annotated rendering changed record, feature or layout semantics",
        ));
    }
    for (left, right) in original.metadata.iter().zip(&parsed.metadata) {
        if left.key != right.key
            || left.continuation != right.continuation
            || metadata_value(left, original, format)? != right.value
        {
            return Err(error(
                format,
                IoErrorKind::DataLoss,
                "annotated rendering changed metadata semantics",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Dna, Rna};

    const GENBANK: &str = "LOCUS       TEST 4 bp DNA linear SYN 16-OCT-2026\nDEFINITION  Writer example.\nACCESSION   TEST\nVERSION     TEST.1\nCOMMENT     original comment\nFEATURES             Location/Qualifiers\n     misc_feature    1..4\n                     /note=\"first \"\"quoted\"\"\n                     /second line\"\n                     /pseudo\n                     /number=2\nORIGIN\n        1 aCgT\n//\n";

    #[test]
    fn model_values_replace_raw_caches_without_reflowing_qualifiers() {
        let mut record = super::super::parse_genbank::<Dna>(GENBANK)
            .unwrap()
            .remove(0);
        record.raw_record = "not an annotated document".into();
        record.metadata[4].value = "changed comment".into();
        record.metadata[4].raw_line = "COMMENT     stale comment".into();
        record.features[0].raw_lines = vec!["stale feature cache".into()];
        record.features[0].location.raw = "stale location cache".into();
        record.features[0].qualifiers[2].value = Some("3".into());
        let rendered = render(&record, Format::GenBank).unwrap();
        assert!(rendered.contains("COMMENT     changed comment\n"));
        assert!(
            rendered
                .contains("/note=\"first \"\"quoted\"\"\n                     /second line\"\n")
        );
        assert!(rendered.contains("/pseudo\n                     /number=3\n"));
        assert!(!rendered.contains("stale"));
        assert!(rendered.contains("FEATURES             Location/Qualifiers\n"));
        assert!(rendered.contains("ORIGIN\n"));
        assert!(rendered.contains("        1 aCgT\n"));
    }

    #[test]
    fn cross_format_and_semantic_mismatch_fail_before_external_writes() {
        let mut record = super::super::parse_genbank::<Dna>(GENBANK)
            .unwrap()
            .remove(0);
        let mut bytes = Vec::new();
        assert_eq!(
            write_embl(&mut bytes, &record).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
        assert!(bytes.is_empty());
        record.metadata[1].value = "Different derived description.".into();
        assert_eq!(
            write_genbank(&mut bytes, &record).unwrap_err().kind(),
            IoErrorKind::DataLoss
        );
        assert!(bytes.is_empty());
    }

    #[test]
    fn embl_counts_are_literal_and_opaque_sq_tails_are_not_dropped() {
        let text =
            "ID   RNA; SV 1; linear; RNA; STD; SYN; 5 BP.\nSQ   Sequence 5 BP;\n     aUcgn 5\n//\n";
        let record = super::super::parse_embl::<Rna>(text).unwrap().remove(0);
        let rendered = render(&record, Format::Embl).unwrap();
        assert!(rendered.contains("SQ   Sequence 5 BP; 1 A; 1 C; 1 G; 0 T; 2 other;\n"));
        let opaque = text.replace("Sequence 5 BP;", "Sequence 5 BP; checksum opaque;");
        let record = super::super::parse_embl::<Rna>(&opaque).unwrap().remove(0);
        let mut bytes = Vec::new();
        assert_eq!(
            write_embl(&mut bytes, &record).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
        assert!(bytes.is_empty());
    }

    #[test]
    fn empty_metadata_fields_use_unpadded_canonical_lines() {
        let text =
            "ID   EMPTY; SV 1; linear; DNA; STD; SYN; 0 BP.\nXX\nFH\nSQ   Sequence 0 BP;\n//\n";
        let record = super::super::parse_embl::<Dna>(text).unwrap().remove(0);
        let rendered = render(&record, Format::Embl).unwrap();
        assert!(rendered.contains("\nXX\nFH\nSQ   Sequence 0 BP;"));
        assert!(rendered.lines().all(|line| !line.ends_with(' ')));
    }

    #[test]
    fn physical_line_and_total_limits_are_inclusive_and_preflighted() {
        let mut output = Output {
            text: String::new(),
            format: Format::GenBank,
            limits: ReadLimits {
                max_line_bytes: 4,
                max_total_bytes: 5,
                max_record_symbols: 10,
            },
        };
        output.line(&["abc"]).unwrap();
        output.line(&[""]).unwrap();
        assert_eq!(output.text, "abc\n\n");
        assert_eq!(
            output.line(&[""]).unwrap_err().kind(),
            IoErrorKind::LimitExceeded
        );
        assert_eq!(
            output.line(&["abcd"]).unwrap_err().kind(),
            IoErrorKind::LimitExceeded
        );
        assert_eq!(output.text, "abc\n\n");
    }

    #[test]
    fn stream_failures_keep_io_category_and_may_leave_a_valid_prefix() {
        struct Failing {
            written: Vec<u8>,
        }
        impl Write for Failing {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                if !self.written.is_empty() {
                    return Err(std::io::Error::other("deliberate stream failure"));
                }
                self.written.extend_from_slice(&bytes[..4]);
                Ok(4)
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let record = super::super::parse_genbank::<Dna>(GENBANK)
            .unwrap()
            .remove(0);
        let mut writer = Failing {
            written: Vec::new(),
        };
        let failure = write_genbank(&mut writer, &record).unwrap_err();
        assert_eq!(failure.kind(), IoErrorKind::Io);
        assert_eq!(writer.written, b"LOCU");
        assert!(std::error::Error::source(&failure).is_some());
    }
}
