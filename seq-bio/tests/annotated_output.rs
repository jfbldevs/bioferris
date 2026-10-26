//! Independent annotated-output preservation, coordinate and publication checks.

use std::fs;
use std::io::{self, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::io::annotated::{
    AnnotatedFormat, AnnotatedRecord, AnnotationChange, AnnotationLossPolicy, AnnotationOperation,
    AnnotationTransformOptions, LocationExpression, parse_embl, parse_genbank, write_embl,
    write_embl_new, write_genbank, write_genbank_new,
};
use bioferris::io::{BioIoError, IoErrorKind};
use bioferris::prelude::*;
use proptest::prelude::*;

const FORMATS: [AnnotatedFormat; 2] = [AnnotatedFormat::GenBank, AnnotatedFormat::Embl];

const LITERAL_GENBANK: &str = r#"LOCUS       ORIGINAL 18 bp DNA linear SYN 15-OCT-2026
DEFINITION  A wrapped
            description.
ACCESSION   TEST0001 SECONDARY
VERSION     TEST0001.2
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            unclassified synthetic sequences.
DBLINK      BioProject: PRJNA000000
COMMENT     Uninterpreted metadata remains ordered.
FEATURES             Location/Qualifiers
     source          1..18
                     /organism="synthetic construct"
                     /mol_type="other DNA"
     CDS             complement(join(1..3,
                     10..12))
                     /gene="quoted_gene"
                     /note="first ""quoted""
                     second /literal"
                     /db_xref="DB:1"
                     /db_xref="DB:2"
                     /codon_start=2
                     /transl_table=11
                     /pseudo
     misc_feature    <1..>18
     misc_feature    REMOTE.1:1..3
ORIGIN
        1 aTggcttaac ccgggttt
//
"#;

const LITERAL_EMBL: &str = r#"ID   TEST0001; SV 2; linear; genomic DNA; STD; SYN; 18 BP.
XX
AC   TEST0001; SECONDARY;
DE   A wrapped
DE   description.
OS   synthetic construct
OC   unclassified synthetic sequences.
DR   MD5; this-summary-is-not-authoritative.
FH   Key             Location/Qualifiers
FH
FT   source          1..18
FT                   /organism="synthetic construct"
FT   CDS             complement(join(1..3,
FT                   10..12))
FT                   /note="first ""quoted""
FT                   second /literal"
FT                   /db_xref="DB:1"
FT                   /db_xref="DB:2"
FT                   /pseudo
XX
SQ   Sequence 18 BP; 99 A; 98 C; 97 G; 96 T; 95 other;
     aTggcttaac ccgggttt 18
//
"#;

struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-annotated-output-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn count(&self) -> usize {
        fs::read_dir(&self.0).unwrap().count()
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn parse<A: Alphabet>(format: AnnotatedFormat, text: &str) -> AnnotatedRecord<A> {
    let mut records = match format {
        AnnotatedFormat::GenBank => parse_genbank(text),
        AnnotatedFormat::Embl => parse_embl(text),
    }
    .unwrap();
    assert_eq!(records.len(), 1);
    records.remove(0)
}

fn write<A: Alphabet>(
    format: AnnotatedFormat,
    writer: impl Write,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    match format {
        AnnotatedFormat::GenBank => write_genbank(writer, record),
        AnnotatedFormat::Embl => write_embl(writer, record),
    }
}

fn write_new<A: Alphabet>(
    format: AnnotatedFormat,
    path: impl AsRef<Path>,
    record: &AnnotatedRecord<A>,
) -> Result<(), BioIoError> {
    match format {
        AnnotatedFormat::GenBank => write_genbank_new(path, record),
        AnnotatedFormat::Embl => write_embl_new(path, record),
    }
}

fn render<A: Alphabet>(record: &AnnotatedRecord<A>) -> String {
    let mut bytes = Vec::new();
    write(record.format(), &mut bytes, record).unwrap();
    String::from_utf8(bytes).unwrap()
}

fn semantic_metadata<A: Alphabet>(record: &AnnotatedRecord<A>) -> Vec<(&str, &str, bool)> {
    record
        .metadata()
        .iter()
        .filter(|field| !matches!(field.key.as_str(), "LOCUS" | "ID" | "SQ"))
        .map(|field| (field.key.as_str(), field.value.as_str(), field.continuation))
        .collect()
}

fn assert_semantics<A: Alphabet>(before: &AnnotatedRecord<A>, after: &AnnotatedRecord<A>) {
    assert_eq!(before.format(), after.format());
    assert_eq!(before.record(), after.record());
    assert_eq!(before.locus(), after.locus());
    assert_eq!(before.declared_molecule(), after.declared_molecule());
    assert_eq!(before.accessions(), after.accessions());
    assert_eq!(semantic_metadata(before), semantic_metadata(after));
    assert_eq!(before.features().len(), after.features().len());
    for (old, new) in before.features().iter().zip(after.features()) {
        assert_eq!(old.key, new.key);
        assert_eq!(old.location.expression(), new.location.expression());
        assert_eq!(
            old.qualifiers
                .iter()
                .map(|qualifier| (&qualifier.name, &qualifier.value))
                .collect::<Vec<_>>(),
            new.qualifiers
                .iter()
                .map(|qualifier| (&qualifier.name, &qualifier.value))
                .collect::<Vec<_>>()
        );
    }
}

fn assert_canonical_round_trip<A: Alphabet>(record: &AnnotatedRecord<A>) -> String {
    let original = record.clone();
    let output = render(record);
    let parsed = parse(record.format(), &output);
    assert_semantics(record, &parsed);
    assert_eq!(parsed.raw_record(), output);
    assert_eq!(render(&parsed), output);
    assert_eq!(record, &original);
    assert!(output.ends_with("\n//\n"));
    output
}

fn feature(format: AnnotatedFormat, key: &str, location: &str, qualifiers: &[&str]) -> String {
    let prefix = match format {
        AnnotatedFormat::GenBank => "     ",
        AnnotatedFormat::Embl => "FT   ",
    };
    let mut text = format!("{prefix}{key:<15} {location}\n");
    for qualifier in qualifiers {
        text.push_str(&format!("{prefix}                {qualifier}\n"));
    }
    text
}

fn document(format: AnnotatedFormat, molecule: &str, sequence: &str, features: &str) -> String {
    let length = sequence.len();
    match format {
        AnnotatedFormat::GenBank => {
            let unit = if molecule == "protein" { "aa" } else { "bp" };
            let token = if molecule == "protein" { "" } else { molecule };
            let row = if sequence.is_empty() {
                String::new()
            } else {
                format!("        1 {sequence}\n")
            };
            format!(
                "LOCUS       ORIGINAL {length} {unit} {token} linear SYN 15-OCT-2026\nDEFINITION  Original annotation.\nFEATURES             Location/Qualifiers\n{features}ORIGIN\n{row}//\n"
            )
        }
        AnnotatedFormat::Embl => {
            let row = if sequence.is_empty() {
                String::new()
            } else {
                format!("     {sequence} {length}\n")
            };
            format!(
                "ID   ORIGINAL; SV 1; linear; {molecule}; STD; SYN; {length} BP.\nDE   Original annotation.\nFH   Key             Location/Qualifiers\nFH\n{features}SQ   Sequence {length} BP;\n{row}//\n"
            )
        }
    }
}

fn options(id: &str) -> AnnotationTransformOptions {
    AnnotationTransformOptions::new(id)
        .unwrap()
        .with_date("16-OCT-2026")
        .unwrap()
}

fn discard(id: &str) -> AnnotationTransformOptions {
    options(id).with_loss_policy(AnnotationLossPolicy::Discard)
}

fn local_id(format: AnnotatedFormat, base: &str) -> String {
    match format {
        AnnotatedFormat::GenBank => base.to_owned(),
        AnnotatedFormat::Embl => format!("{base}.1"),
    }
}

#[test]
fn same_format_writers_preserve_literal_semantics_not_original_layout() {
    for (format, text) in [
        (AnnotatedFormat::GenBank, LITERAL_GENBANK),
        (AnnotatedFormat::Embl, LITERAL_EMBL),
    ] {
        let record = parse::<Dna>(format, text);
        assert_eq!(record.record().id(), "TEST0001.2");
        assert_eq!(
            record.record().description(),
            Some("A wrapped description.")
        );
        assert_eq!(record.record().source_sequence(), "aTggcttaacccgggttt");
        assert_eq!(record.accessions(), ["TEST0001", "SECONDARY"]);
        let output = assert_canonical_round_trip(&record);
        assert_ne!(output, text);
        assert!(output.contains("/db_xref=\"DB:1\""));
        assert!(output.contains("/db_xref=\"DB:2\""));
        assert!(output.contains("/pseudo\n"));
        let parsed = parse::<Dna>(format, &output);
        let note = parsed.features()[1]
            .qualifiers
            .iter()
            .find(|qualifier| qualifier.name == "note")
            .unwrap();
        assert_eq!(
            note.value.as_deref(),
            Some("first \"quoted\"\nsecond /literal")
        );
    }
}

#[test]
fn embl_summary_is_regenerated_from_symbols_not_cached_counts() {
    let record = parse::<Dna>(AnnotatedFormat::Embl, LITERAL_EMBL);
    let output = assert_canonical_round_trip(&record);
    assert!(output.contains("SQ   Sequence 18 BP; 3 A; 4 C; 5 G; 6 T; 0 other;\n"));
    assert!(!output.contains("99 A"));
    for (molecule, sequence, summary) in [
        ("DNA", "aCgTnR", "1 A; 1 C; 1 G; 1 T; 2 other;"),
        ("RNA", "aCgUnR", "1 A; 1 C; 1 G; 0 T; 3 other;"),
    ] {
        let text = document(AnnotatedFormat::Embl, molecule, sequence, "");
        let output = if molecule == "DNA" {
            assert_canonical_round_trip(&parse::<Dna>(AnnotatedFormat::Embl, &text))
        } else {
            assert_canonical_round_trip(&parse::<Rna>(AnnotatedFormat::Embl, &text))
        };
        assert!(output.contains(summary), "{output}");
    }
}

#[test]
fn real_public_fixtures_round_trip_without_relaxing_unsupported_locations() {
    for (format, input, expected_id, expected_length) in [
        (
            AnnotatedFormat::GenBank,
            include_str!("fixtures/annotated-dq067560.1.gb"),
            "DQ067560.1",
            534,
        ),
        (
            AnnotatedFormat::Embl,
            include_str!("fixtures/annotated-dq067560.1.embl"),
            "DQ067560.1",
            534,
        ),
        (
            AnnotatedFormat::Embl,
            include_str!("fixtures/annotated-m77789.2.embl"),
            "M77789.2",
            2686,
        ),
    ] {
        let record = parse::<Dna>(format, input);
        assert_eq!(record.record().id(), expected_id);
        assert_eq!(record.record().sequence().len(), expected_length);
        let output = assert_canonical_round_trip(&record);
        let reopened = parse::<Dna>(format, &output);
        for (before, after) in record.features().iter().zip(reopened.features()) {
            match before.extract(record.record().sequence()) {
                Ok(expected) => assert_eq!(
                    after.extract(reopened.record().sequence()).unwrap(),
                    expected
                ),
                Err(error) => assert_eq!(
                    after
                        .extract(reopened.record().sequence())
                        .unwrap_err()
                        .kind(),
                    error.kind()
                ),
            }
        }
    }
}

#[test]
fn writers_support_empty_sequences_rna_and_genbank_protein_without_conversion() {
    for format in FORMATS {
        assert_canonical_round_trip(&parse::<Dna>(format, &document(format, "DNA", "", "")));
        let nucleic = parse::<Rna>(format, &document(format, "mRNA", "aUgCryN", ""));
        assert_eq!(nucleic.record().source_sequence(), "aUgCryN");
        assert_canonical_round_trip(&nucleic);
    }
    let format = AnnotatedFormat::GenBank;
    let peptide = parse::<Protein>(format, &document(format, "protein", "mKBZJXUO", ""));
    assert_eq!(peptide.record().source_sequence(), "mKBZJXUO");
    assert_canonical_round_trip(&peptide);
}

#[test]
fn unsupported_location_text_is_retained_by_writers_but_not_certified() {
    for format in FORMATS {
        for location in [
            "<1..>5",
            "REMOTE.1:1..3",
            "2^3",
            "order(1,5)",
            "join(1,OTHER:2)",
        ] {
            let text = document(
                format,
                "DNA",
                "aCgTN",
                &feature(format, "misc_feature", location, &[]),
            );
            let record = parse::<Dna>(format, &text);
            let output = assert_canonical_round_trip(&record);
            let parsed = parse::<Dna>(format, &output);
            assert_eq!(
                parsed.features()[0]
                    .extract(parsed.record().sequence())
                    .unwrap_err()
                    .kind(),
                IoErrorKind::Unsupported
            );
        }
    }
}

#[test]
fn cross_format_output_is_rejected_before_writing_or_creating_anything() {
    let directory = Directory::new();
    for format in FORMATS {
        let record = parse::<Dna>(format, &document(format, "DNA", "AC", ""));
        let other = match format {
            AnnotatedFormat::GenBank => AnnotatedFormat::Embl,
            AnnotatedFormat::Embl => AnnotatedFormat::GenBank,
        };
        let mut bytes = Vec::new();
        assert_eq!(
            write(other, &mut bytes, &record).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
        assert!(bytes.is_empty());
        assert_eq!(
            write_new(other, directory.path("must-not-exist"), &record)
                .unwrap_err()
                .kind(),
            IoErrorKind::Unsupported
        );
        assert_eq!(directory.count(), 0);
    }
}

struct FailingWriter {
    bytes: Vec<u8>,
    remaining: usize,
}

impl Write for FailingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.remaining == 0 {
            return Err(io::Error::other("independent injected writer failure"));
        }
        let count = bytes.len().min(self.remaining);
        self.bytes.extend_from_slice(&bytes[..count]);
        self.remaining -= count;
        Ok(count)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn late_stream_write_errors_are_observable_without_changing_the_source() {
    for format in FORMATS {
        let record = parse::<Dna>(format, &document(format, "DNA", "aCgTrYsW", ""));
        let original = record.clone();
        let expected = render(&record);
        let remaining = expected.len() - 3;
        let mut writer = FailingWriter {
            bytes: Vec::new(),
            remaining,
        };
        let error = write(format, &mut writer, &record).unwrap_err();
        assert_eq!(error.kind(), IoErrorKind::Io);
        assert_eq!(error.code(), "io_error");
        assert_eq!(writer.bytes, expected.as_bytes()[..remaining]);
        assert_eq!(record, original);
    }
}

#[test]
fn named_outputs_are_create_only_and_leave_no_partial_file_on_failure() {
    for format in FORMATS {
        let directory = Directory::new();
        let record = parse::<Dna>(format, &document(format, "DNA", "aCgT", ""));
        let path = directory.path("résultat.data");
        write_new(format, &path, &record).unwrap();
        assert_eq!(fs::read(&path).unwrap(), render(&record).as_bytes());
        let previous = fs::read(&path).unwrap();
        assert_eq!(
            write_new(format, &path, &record).unwrap_err().kind(),
            IoErrorKind::Io
        );
        assert_eq!(fs::read(&path).unwrap(), previous);
        assert!(write_new(format, &directory.0, &record).is_err());
        assert!(write_new(format, directory.path("missing-parent/output"), &record).is_err());
        assert_eq!(directory.count(), 1);
    }
}

#[cfg(unix)]
#[test]
fn dangling_named_output_links_are_preserved_without_following_them() {
    use std::os::unix::fs::symlink;

    let directory = Directory::new();
    let destination = directory.path("dangling");
    let target = directory.path("absent");
    symlink(&target, &destination).unwrap();
    for format in FORMATS {
        let record = parse::<Dna>(format, &document(format, "DNA", "AC", ""));
        assert_eq!(
            write_new(format, &destination, &record).unwrap_err().kind(),
            IoErrorKind::Io
        );
        assert_eq!(fs::read_link(&destination).unwrap(), target);
        assert!(!target.exists());
        assert_eq!(directory.count(), 1);
    }
}

fn complement(byte: u8, rna: bool) -> u8 {
    let upper = match byte.to_ascii_uppercase() {
        b'A' => {
            if rna {
                b'U'
            } else {
                b'T'
            }
        }
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        other => panic!("invalid oracle nucleotide {other}"),
    };
    if byte.is_ascii_lowercase() {
        upper.to_ascii_lowercase()
    } else {
        upper
    }
}

fn reverse_complement(text: &str, rna: bool) -> String {
    text.bytes()
        .rev()
        .map(|byte| char::from(complement(byte, rna)))
        .collect()
}

#[derive(Clone)]
enum OracleLocation {
    Span(Range<usize>),
    Join(Vec<Self>),
    Complement(Box<Self>),
}

impl OracleLocation {
    fn text(&self) -> String {
        match self {
            Self::Span(range) => format!("{}..{}", range.start + 1, range.end),
            Self::Join(parts) => format!(
                "join({})",
                parts.iter().map(Self::text).collect::<Vec<_>>().join(",")
            ),
            Self::Complement(inner) => format!("complement({})", inner.text()),
        }
    }

    fn extract(&self, sequence: &str, rna: bool) -> String {
        match self {
            Self::Span(range) => sequence[range.clone()].to_ascii_uppercase(),
            Self::Join(parts) => parts
                .iter()
                .map(|part| part.extract(sequence, rna))
                .collect(),
            Self::Complement(inner) => reverse_complement(&inner.extract(sequence, rna), rna),
        }
    }
}

fn assert_raw_cache<A: Alphabet>(record: &AnnotatedRecord<A>) {
    let parsed = parse::<A>(record.format(), record.raw_record());
    assert_eq!(&parsed, record);
    assert_eq!(render(record), record.raw_record());
    for feature in record.features() {
        for line in &feature.raw_lines {
            assert!(record.raw_record().lines().any(|raw| raw == line));
        }
    }
}

#[test]
fn local_identity_and_calendar_options_are_validated_without_guessed_dates() {
    for id in ["LOCAL_a", "LOCAL_1234567890", "LOCAL_A_B"] {
        let configured = AnnotationTransformOptions::new(id).unwrap();
        assert_eq!(configured.local_id(), id);
        assert_eq!(configured.loss_policy(), AnnotationLossPolicy::Reject);
        assert_eq!(configured.date(), None);
    }
    for id in [
        "",
        "LOCAL_",
        "local_a",
        "LOCAL_a-b",
        "LOCAL_a.b",
        "LOCAL_é",
        "LOCAL_a b",
        "LOCAL_12345678901",
        "ACCESSION.1",
    ] {
        assert_eq!(
            AnnotationTransformOptions::new(id).unwrap_err().kind(),
            IoErrorKind::InvalidHeader,
            "{id}"
        );
    }
    for date in ["01-JAN-0001", "29-FEB-2000", "29-FEB-2024", "31-DEC-9999"] {
        let configured = options("LOCAL_date").with_date(date).unwrap();
        assert_eq!(configured.date(), Some(date));
    }
    for date in [
        "",
        "2026-10-16",
        "1-OCT-2026",
        "00-JAN-2026",
        "32-JAN-2026",
        "31-APR-2026",
        "29-FEB-1900",
        "29-FEB-2026",
        "01-XYZ-2026",
        "01-JAN-0000",
        "01-JAN-10000",
        "16-OCT-2026\n",
    ] {
        assert_eq!(
            options("LOCAL_date").with_date(date).unwrap_err().kind(),
            IoErrorKind::InvalidHeader,
            "{date}"
        );
    }
    for format in FORMATS {
        let entry = parse::<Dna>(format, &document(format, "DNA", "AC", ""));
        let result = entry.slice(0..2, AnnotationTransformOptions::new("LOCAL_date").unwrap());
        match format {
            AnnotatedFormat::GenBank => {
                assert_eq!(result.unwrap_err().kind(), IoErrorKind::InvalidHeader);
            }
            AnnotatedFormat::Embl => {
                let transformed = result.unwrap();
                assert!(
                    !transformed
                        .record()
                        .metadata()
                        .iter()
                        .any(|field| field.key == "DT")
                );
                assert_raw_cache(transformed.record());
            }
        }
    }
}

#[test]
fn new_identity_description_and_report_are_explicit_even_for_full_slice() {
    for format in FORMATS {
        let input = document(format, "DNA", "aCgTN", "");
        let entry = parse::<Dna>(format, &input);
        let original = entry.clone();
        let source_id = entry.record().id().to_owned();
        let transformed = entry.slice(0..5, options("LOCAL_copy")).unwrap();
        let derived = transformed.record();
        assert_eq!(derived.record().id(), local_id(format, "LOCAL_copy"));
        assert_eq!(derived.locus(), "LOCAL_copy");
        assert!(derived.accessions().is_empty());
        assert_eq!(derived.record().source_sequence(), "aCgTN");
        assert_eq!(
            derived.record().description(),
            Some(format!("BioFerris slice [0,5) of {source_id}.").as_str())
        );
        assert_eq!(transformed.report().source_id, source_id);
        assert_eq!(transformed.report().derived_id, derived.record().id());
        assert_eq!(
            transformed.report().operation,
            AnnotationOperation::Slice { start: 0, end: 5 }
        );
        assert!(!transformed.report().changes.is_empty());
        assert_raw_cache(derived);
        assert_eq!(transformed.clone().into_record(), *derived);
        assert_eq!(entry, original);
        let reused = derived.slice(0..5, options("LOCAL_copy")).unwrap_err();
        assert_eq!(reused.kind(), IoErrorKind::InvalidHeader);
    }
}

#[test]
fn every_empty_and_nonempty_sequence_slice_is_strict_and_case_preserving() {
    let sequence = "aCgTrY";
    for format in FORMATS {
        let entry = parse::<Dna>(format, &document(format, "DNA", sequence, ""));
        let original = entry.clone();
        for start in 0..=sequence.len() {
            for end in start..=sequence.len() {
                let result = entry.slice(start..end, options("LOCAL_slice")).unwrap();
                assert_eq!(
                    result.record().record().source_sequence(),
                    &sequence[start..end]
                );
                assert_eq!(
                    result.record().record().sequence().to_string(),
                    sequence[start..end].to_ascii_uppercase()
                );
                assert!(result.record().features().is_empty());
                assert_raw_cache(result.record());
            }
        }
        for (start, end) in [(4, 3), (0, 7), (7, 7), (usize::MAX, usize::MAX)] {
            assert_eq!(
                entry
                    .slice(start..end, options("LOCAL_bad"))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::LengthMismatch
            );
        }
        assert_eq!(entry, original);
    }
}

#[test]
fn exhaustive_single_feature_slices_distinguish_contained_outside_and_partial() {
    let sequence = "aCgTr";
    for format in FORMATS {
        for first in 0..sequence.len() {
            for last in first + 1..=sequence.len() {
                let location = format!("{}..{last}", first + 1);
                let text = document(
                    format,
                    "DNA",
                    sequence,
                    &feature(format, "gene", &location, &["/gene=\"example\""]),
                );
                let entry = parse::<Dna>(format, &text);
                let original = entry.clone();
                for start in 0..=sequence.len() {
                    for end in start..=sequence.len() {
                        let contained = first >= start && last <= end;
                        let intersects = start < end && first < end && start < last;
                        for policy in [AnnotationLossPolicy::Reject, AnnotationLossPolicy::Discard]
                        {
                            let result = entry.slice(
                                start..end,
                                options("LOCAL_feature").with_loss_policy(policy),
                            );
                            if contained {
                                let transformed = result.unwrap();
                                let retained = &transformed.record().features()[0];
                                assert_eq!(
                                    retained.location.expression(),
                                    &LocationExpression::Range(first - start..last - start)
                                );
                                assert_eq!(
                                    retained
                                        .extract(transformed.record().record().sequence())
                                        .unwrap()
                                        .to_string(),
                                    sequence[first..last].to_ascii_uppercase()
                                );
                                assert_eq!(
                                    retained.qualifiers[0].value.as_deref(),
                                    Some("example")
                                );
                            } else if !intersects && policy == AnnotationLossPolicy::Discard {
                                let transformed = result.unwrap();
                                assert!(transformed.record().features().is_empty());
                                assert!(transformed.report().changes.iter().any(|change| matches!(change,
                                    AnnotationChange::FeatureRemoved { index: 0, key, .. } if key == "gene")));
                            } else {
                                assert_eq!(
                                    result.unwrap_err().kind(),
                                    IoErrorKind::DataLoss,
                                    "feature {first}..{last}, slice {start}..{end}, policy {policy:?}"
                                );
                            }
                        }
                    }
                }
                assert_eq!(entry, original);
            }
        }
    }
}

#[test]
fn whole_record_source_coverage_resizes_but_empty_removal_requires_permission() {
    for format in FORMATS {
        let sequence = "aCgTr";
        let text = document(
            format,
            "DNA",
            sequence,
            &feature(
                format,
                "source",
                "1..5",
                &[
                    "/organism=\"synthetic construct\"",
                    "/mol_type=\"other DNA\"",
                ],
            ),
        );
        let entry = parse::<Dna>(format, &text);
        for start in 0..=sequence.len() {
            for end in start..=sequence.len() {
                if start == end {
                    assert_eq!(
                        entry
                            .slice(start..end, options("LOCAL_empty"))
                            .unwrap_err()
                            .kind(),
                        IoErrorKind::DataLoss
                    );
                    let result = entry.slice(start..end, discard("LOCAL_empty")).unwrap();
                    assert!(result.record().features().is_empty());
                    assert!(result.record().record().sequence().is_empty());
                } else {
                    let result = entry.slice(start..end, options("LOCAL_region")).unwrap();
                    let source = &result.record().features()[0];
                    assert_eq!(
                        source.location.expression(),
                        &LocationExpression::Range(0..end - start)
                    );
                    assert_eq!(source.qualifiers, entry.features()[0].qualifiers);
                    if start != 0 || end != sequence.len() {
                        assert!(result.report().changes.iter().any(|change| matches!(
                            change,
                            AnnotationChange::LocationChanged {
                                feature_index: 0,
                                source_coverage: true,
                                ..
                            }
                        )));
                    }
                }
            }
        }
        let reversed = entry.reverse_complement(options("LOCAL_reverse")).unwrap();
        assert_eq!(
            reversed.record().features()[0].location.expression(),
            &LocationExpression::Range(0..5)
        );
        for location in ["complement(1..5)", "join(1..2,3..5)", "2..5"] {
            let nonspecial = parse::<Dna>(
                format,
                &document(
                    format,
                    "DNA",
                    sequence,
                    &feature(format, "source", location, &[]),
                ),
            );
            assert_eq!(
                nonspecial
                    .slice(2..4, discard("LOCAL_partial"))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::DataLoss
            );
        }
    }
}

#[test]
fn joined_features_cannot_be_partially_clipped_or_partially_discarded() {
    for format in FORMATS {
        for location in [
            "join(2..3,7..8)",
            "complement(join(2..3,7..8))",
            "join(complement(7..8),2..3)",
        ] {
            let entry = parse::<Dna>(
                format,
                &document(
                    format,
                    "DNA",
                    "aCgTrYsWkM",
                    &feature(format, "CDS", location, &["/translation=\"XX\""]),
                ),
            );
            for range in [0..5, 2..8, 3..8, 2..7, 6..10] {
                assert_eq!(
                    entry
                        .slice(range, discard("LOCAL_clip"))
                        .unwrap_err()
                        .kind(),
                    IoErrorKind::DataLoss
                );
            }
            let contained = entry.slice(1..8, options("LOCAL_join")).unwrap();
            assert_eq!(
                contained.record().features()[0]
                    .extract(contained.record().record().sequence())
                    .unwrap(),
                entry.features()[0]
                    .extract(entry.record().sequence())
                    .unwrap()
            );
            assert_eq!(
                contained.record().features()[0].qualifiers[0]
                    .value
                    .as_deref(),
                Some("XX")
            );
            let outside = entry.slice(3..6, discard("LOCAL_gap")).unwrap();
            assert!(outside.record().features().is_empty());
        }
    }
}

fn check_reverse_locations<A: Alphabet>(molecule: &str, sequence: &str, rna: bool) {
    for format in FORMATS {
        for first_start in 0..sequence.len() {
            for first_end in first_start + 1..=sequence.len() {
                for second_start in 0..sequence.len() {
                    for second_end in second_start + 1..=sequence.len() {
                        let first = OracleLocation::Span(first_start..first_end);
                        let second = OracleLocation::Span(second_start..second_end);
                        for location in [
                            first.clone(),
                            OracleLocation::Join(vec![first.clone(), second.clone()]),
                            OracleLocation::Complement(Box::new(OracleLocation::Join(vec![
                                first.clone(),
                                second.clone(),
                            ]))),
                            OracleLocation::Join(vec![
                                OracleLocation::Complement(Box::new(first.clone())),
                                second.clone(),
                            ]),
                            OracleLocation::Join(vec![
                                first.clone(),
                                OracleLocation::Complement(Box::new(second.clone())),
                            ]),
                            OracleLocation::Complement(Box::new(OracleLocation::Complement(
                                Box::new(first),
                            ))),
                        ] {
                            let input = document(
                                format,
                                molecule,
                                sequence,
                                &feature(
                                    format,
                                    "gene",
                                    &location.text(),
                                    &["/gene=\"preserved\""],
                                ),
                            );
                            let entry = parse::<A>(format, &input);
                            let expected = location.extract(sequence, rna);
                            let result = entry.reverse_complement(options("LOCAL_rc")).unwrap();
                            let derived = result.record();
                            assert_eq!(
                                derived.record().source_sequence(),
                                reverse_complement(sequence, rna)
                            );
                            assert_eq!(derived.features().len(), 1);
                            assert_eq!(
                                derived.features()[0]
                                    .extract(derived.record().sequence())
                                    .unwrap()
                                    .to_string(),
                                expected,
                                "{}",
                                location.text()
                            );
                            assert_eq!(
                                result.report().operation,
                                AnnotationOperation::ReverseComplement
                            );
                            assert_eq!(entry.raw_record(), input);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn exhaustive_join_and_complement_reverse_mapping_preserves_biological_extraction() {
    check_reverse_locations::<Dna>("DNA", "aCgTr", false);
    check_reverse_locations::<Rna>("RNA", "aCgUr", true);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn generated_case_and_ambiguity_match_independent_reverse_oracle(raw in prop::collection::vec(any::<u8>(), 0..128)) {
        let symbols = b"ACGTRYSWKMBDHVN";
        let sequence: String = raw.iter().map(|byte| {
            let value = symbols[usize::from(*byte) % symbols.len()];
            char::from(if byte & 1 == 0 { value } else { value.to_ascii_lowercase() })
        }).collect();
        for format in FORMATS {
            let entry = parse::<Dna>(format, &document(format, "DNA", &sequence, ""));
            let result = entry.reverse_complement(options("LOCAL_rc")).unwrap();
            prop_assert_eq!(result.record().record().source_sequence(), reverse_complement(&sequence, false));
            assert_raw_cache(result.record());
            let restored = result.record().reverse_complement(options("LOCAL_again")).unwrap();
            prop_assert_eq!(restored.record().record().source_sequence(), sequence.as_str());
            prop_assert_eq!(entry.record().source_sequence(), sequence.as_str());
        }
    }
}

#[test]
fn unsupported_locations_cannot_be_guessed_even_with_discard_permission() {
    for format in FORMATS {
        for location in [
            "<1..>5",
            "REMOTE.1:1..3",
            "2^3",
            "order(1,5)",
            "join(1,OTHER:2)",
        ] {
            let entry = parse::<Dna>(
                format,
                &document(
                    format,
                    "DNA",
                    "aCgTN",
                    &feature(format, "misc_feature", location, &[]),
                ),
            );
            let original = entry.clone();
            for policy in [AnnotationLossPolicy::Reject, AnnotationLossPolicy::Discard] {
                for range in [0..5, 1..4, 0..0] {
                    assert_eq!(
                        entry
                            .slice(range, options("LOCAL_unknown").with_loss_policy(policy))
                            .unwrap_err()
                            .kind(),
                        IoErrorKind::Unsupported
                    );
                }
                assert_eq!(
                    entry
                        .reverse_complement(options("LOCAL_unknown").with_loss_policy(policy))
                        .unwrap_err()
                        .kind(),
                    IoErrorKind::Unsupported
                );
            }
            assert_eq!(entry, original);
        }
    }
}

#[test]
fn protein_slice_is_supported_but_reverse_complement_is_never_fabricated() {
    let format = AnnotatedFormat::GenBank;
    let text = document(
        format,
        "protein",
        "mKBZJXUO",
        &feature(format, "Region", "2..4", &["/product=\"region\""]),
    );
    let entry = parse::<Protein>(format, &text);
    let result = entry.slice(1..5, options("LOCAL_peptide")).unwrap();
    assert_eq!(result.record().record().source_sequence(), "KBZJ");
    assert_eq!(
        result.record().features()[0].location.expression(),
        &LocationExpression::Range(0..3)
    );
    assert_canonical_round_trip(result.record());
    for configured in [options("LOCAL_rc"), discard("LOCAL_rc")] {
        assert_eq!(
            entry.reverse_complement(configured).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
    }
}

fn with_metadata(document: &str, format: AnnotatedFormat, metadata: &str) -> String {
    let marker = match format {
        AnnotatedFormat::GenBank => "FEATURES             Location/Qualifiers\n",
        AnnotatedFormat::Embl => "FH   Key             Location/Qualifiers\n",
    };
    document.replacen(marker, &format!("{metadata}{marker}"), 1)
}

#[test]
fn reviewed_qualifiers_are_retained_in_order_with_values_and_flags() {
    let qualifiers = [
        "/gene=\"demo\"",
        "/gene_synonym=\"alias\"",
        "/locus_tag=\"tag_1\"",
        "/product=\"product\"",
        "/protein_id=\"PROTEIN.1\"",
        "/codon_start=2",
        "/transl_table=11",
        "/translation=\"XX\"",
        "/pseudo",
        "/pseudogene=\"unknown\"",
        "/organism=\"synthetic construct\"",
        "/mol_type=\"other DNA\"",
        "/strain=\"strain\"",
        "/isolate=\"isolate\"",
        "/db_xref=\"DB:1\"",
        "/db_xref=\"DB:2\"",
        "/culture_collection=\"CC:1\"",
        "/plasmid=\"plasmid\"",
        "/organelle=\"plastid\"",
    ];
    for format in FORMATS {
        let entry = parse::<Dna>(
            format,
            &document(
                format,
                "DNA",
                "aCgTrY",
                &feature(format, "gene", "2..5", &qualifiers),
            ),
        );
        let expected = entry.features()[0]
            .qualifiers
            .iter()
            .map(|q| (&q.name, &q.value))
            .collect::<Vec<_>>();
        for transformed in [
            entry.slice(1..5, options("LOCAL_keep")).unwrap(),
            entry.reverse_complement(options("LOCAL_keep")).unwrap(),
        ] {
            let actual = transformed.record().features()[0]
                .qualifiers
                .iter()
                .map(|q| (&q.name, &q.value))
                .collect::<Vec<_>>();
            assert_eq!(actual, expected);
            assert_raw_cache(transformed.record());
            assert!(
                !transformed
                    .report()
                    .changes
                    .iter()
                    .any(|change| matches!(change, AnnotationChange::FeatureRemoved { .. }))
            );
        }
    }
}

#[test]
fn unreviewed_qualifiers_drop_the_entire_feature_only_when_authorized() {
    for format in FORMATS {
        for qualifier in [
            "/note=\"coordinates 1..3\"",
            "/transl_except=\"(pos:1..3,aa:Sec)\"",
            "/anticodon=\"(pos:1..3,aa:Leu)\"",
            "/unknown=\"opaque\"",
            "/Gene=\"case-sensitive\"",
        ] {
            let annotations = format!(
                "{}{}",
                feature(format, "gene", "1..3", &["/gene=\"retained\""]),
                feature(format, "CDS", "1..6", &["/translation=\"XX\"", qualifier])
            );
            let entry = parse::<Dna>(format, &document(format, "DNA", "aCgTrY", &annotations));
            let original = entry.clone();
            assert_eq!(
                entry
                    .slice(0..6, options("LOCAL_reject"))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::DataLoss
            );
            assert_eq!(
                entry
                    .reverse_complement(options("LOCAL_reject"))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::DataLoss
            );
            for result in [
                entry.slice(0..6, discard("LOCAL_discard")).unwrap(),
                entry.reverse_complement(discard("LOCAL_discard")).unwrap(),
            ] {
                assert_eq!(result.record().features().len(), 1);
                assert_eq!(result.record().features()[0].key, "gene");
                assert_eq!(
                    result.record().features()[0].qualifiers[0].value.as_deref(),
                    Some("retained")
                );
                let removals: Vec<_> = result
                    .report()
                    .changes
                    .iter()
                    .filter_map(|change| match change {
                        AnnotationChange::FeatureRemoved {
                            index,
                            key,
                            location,
                            reason,
                        } => Some((*index, key.as_str(), location.as_str(), reason)),
                        _ => None,
                    })
                    .collect();
                assert_eq!(removals.len(), 1);
                assert_eq!(
                    (removals[0].0, removals[0].1, removals[0].2),
                    (1, "CDS", "1..6")
                );
                assert!(!removals[0].3.is_empty());
                assert!(!result.record().raw_record().contains(qualifier));
                assert_raw_cache(result.record());
            }
            // Discard cannot hide a partial-feature coordinate problem merely
            // because the same feature also contains an unreviewed qualifier.
            assert_eq!(
                entry
                    .slice(0..5, discard("LOCAL_partial"))
                    .unwrap_err()
                    .kind(),
                IoErrorKind::DataLoss
            );
            assert_eq!(entry, original);
        }
    }
}

#[test]
fn unknown_metadata_and_database_history_require_explicit_reported_loss() {
    for format in FORMATS {
        let (safe, unsafe_metadata, unsafe_keys) = match format {
            AnnotatedFormat::GenBank => (
                "SOURCE      synthetic construct\n  ORGANISM  synthetic construct\n            unclassified synthetic sequences.\n",
                "COMMENT     opaque note\nREFERENCE   1  (bases 1 to 6)\n  AUTHORS   Example,A.\n  TITLE     Source publication\n  JOURNAL   Unpublished\n   PUBMED   12345\nDBLINK      BioProject: PRJNA000000\n",
                vec![
                    "COMMENT",
                    "REFERENCE",
                    "AUTHORS",
                    "TITLE",
                    "JOURNAL",
                    "PUBMED",
                    "DBLINK",
                ],
            ),
            AnnotatedFormat::Embl => (
                "OS   synthetic construct\nOC   unclassified synthetic sequences.\nOG   Plastid\n",
                "CC   opaque note\nRN   [1]\nRA   Example A.;\nRT   \"Source publication\";\nRL   Unpublished.\nDR   MD5; abcdef0123456789.\nDT   01-JAN-2000 (Rel. 1, Created)\n",
                vec!["CC", "RN", "RA", "RT", "RL", "DR", "DT"],
            ),
        };
        let base = document(format, "DNA", "aCgTrY", "");
        let input = with_metadata(&base, format, &format!("{safe}{unsafe_metadata}"));
        let entry = parse::<Dna>(format, &input);
        let original = entry.clone();
        assert_eq!(
            entry
                .slice(0..6, options("LOCAL_reject"))
                .unwrap_err()
                .kind(),
            IoErrorKind::DataLoss
        );
        assert_eq!(
            entry
                .reverse_complement(options("LOCAL_reject"))
                .unwrap_err()
                .kind(),
            IoErrorKind::DataLoss
        );
        for result in [
            entry.slice(1..5, discard("LOCAL_clean")).unwrap(),
            entry.reverse_complement(discard("LOCAL_clean")).unwrap(),
        ] {
            for (index, field) in entry
                .metadata()
                .iter()
                .enumerate()
                .filter(|(_, field)| unsafe_keys.contains(&field.key.as_str()))
            {
                assert!(result.report().changes.iter().any(|change| matches!(change,
                    AnnotationChange::MetadataChanged { index: Some(original_index), key, before: Some(value), after: None, reason }
                    if *original_index == index && key == &field.key && value == &field.value && !reason.is_empty())), "missing removal report for {field:?}");
            }
            for field in entry.metadata().iter().filter(|field| {
                matches!(
                    field.key.as_str(),
                    "SOURCE" | "ORGANISM" | "OS" | "OC" | "OG"
                )
            }) {
                assert!(
                    result
                        .record()
                        .metadata()
                        .iter()
                        .any(|retained| retained.key == field.key
                            && retained.value == field.value
                            && retained.continuation == field.continuation)
                );
            }
            assert!(!result.record().raw_record().contains("abcdef0123456789"));
            assert!(!result.record().raw_record().contains("01-JAN-2000"));
            assert!(!result.record().raw_record().contains("Source publication"));
            if format == AnnotatedFormat::Embl {
                let dates = result
                    .record()
                    .metadata()
                    .iter()
                    .filter(|field| field.key == "DT")
                    .collect::<Vec<_>>();
                assert_eq!(dates.len(), 1);
                assert_eq!(dates[0].value, "16-OCT-2026 (BioFerris local derivation)");
            }
            assert_raw_cache(result.record());
        }
        assert_eq!(entry, original);
    }
}

#[test]
fn accession_and_version_metadata_are_not_reused_as_derived_identity() {
    for format in FORMATS {
        let metadata = match format {
            AnnotatedFormat::GenBank => "ACCESSION   ORIGINAL SECONDARY\nVERSION     ORIGINAL.1\n",
            AnnotatedFormat::Embl => "AC   ORIGINAL; SECONDARY;\nSV   ORIGINAL.1\n",
        };
        let text = with_metadata(&document(format, "DNA", "aCgT", ""), format, metadata);
        let entry = parse::<Dna>(format, &text);
        assert_eq!(entry.record().id(), "ORIGINAL.1");
        let result = entry.slice(0..4, options("LOCAL_identity")).unwrap();
        assert_eq!(
            result.record().record().id(),
            local_id(format, "LOCAL_identity")
        );
        assert!(result.record().accessions().is_empty());
        for field in result.record().metadata() {
            assert!(!matches!(
                field.key.as_str(),
                "ACCESSION" | "VERSION" | "AC" | "SV"
            ));
        }
        assert_eq!(result.report().source_id, "ORIGINAL.1");
        assert_eq!(
            result.record().record().description(),
            Some("BioFerris slice [0,4) of ORIGINAL.1.")
        );
        assert_raw_cache(result.record());
    }
}

#[test]
fn circular_topology_is_retained_only_for_whole_sequence_operations() {
    for format in FORMATS {
        let circular = document(format, "DNA", "aCgTrY", "").replacen("linear", "circular", 1);
        let entry = parse::<Dna>(format, &circular);
        for result in [
            entry.slice(0..6, options("LOCAL_whole")).unwrap(),
            entry.reverse_complement(options("LOCAL_whole")).unwrap(),
        ] {
            assert!(
                result
                    .record()
                    .raw_record()
                    .lines()
                    .next()
                    .unwrap()
                    .contains("circular")
            );
            assert_raw_cache(result.record());
        }
        for range in [0..5, 1..6, 1..5, 0..0, 6..6] {
            let result = entry.slice(range, options("LOCAL_part")).unwrap();
            let header = result.record().raw_record().lines().next().unwrap();
            assert!(header.contains("linear"));
            assert!(!header.contains("circular"));
        }
        let unknown = document(format, "DNA", "AC", "").replacen("linear", "unknown", 1);
        let entry = parse::<Dna>(format, &unknown);
        assert_canonical_round_trip(&entry);
        assert_eq!(
            entry
                .slice(0..2, discard("LOCAL_unknown"))
                .unwrap_err()
                .kind(),
            IoErrorKind::Unsupported
        );
        assert_eq!(
            entry
                .reverse_complement(discard("LOCAL_unknown"))
                .unwrap_err()
                .kind(),
            IoErrorKind::Unsupported
        );
    }
}

#[test]
fn source_classification_is_preserved_without_marking_natural_records_synthetic() {
    for (format, divisions) in [
        (
            AnnotatedFormat::GenBank,
            "PRI ROD MAM VRT INV PLN BCT VRL PHG SYN UNA EST PAT STS GSS HTG HTC ENV",
        ),
        (
            AnnotatedFormat::Embl,
            "PHG ENV FUN HUM INV MAM VRT MUS PLN PRO ROD SYN TGN UNC VRL",
        ),
    ] {
        for division in divisions.split_whitespace() {
            let text = document(format, "DNA", "aCgT", "").replacen("SYN", division, 1);
            let entry = parse::<Dna>(format, &text);
            for result in [
                entry.slice(1..3, options("LOCAL_part")).unwrap(),
                entry.reverse_complement(discard("LOCAL_rc")).unwrap(),
            ] {
                let header = &result.record().metadata()[0].value;
                let retained = match format {
                    AnnotatedFormat::GenBank => header.split_whitespace().nth(5).unwrap(),
                    AnnotatedFormat::Embl => header.split(';').nth(5).unwrap().trim(),
                };
                assert_eq!(retained, division);
                assert_raw_cache(result.record());
            }
            assert_eq!(entry.raw_record(), text);
        }
    }
    for data_class in [
        "STD", "PAT", "EST", "GSS", "HTC", "HTG", "WGS", "TSA", "STS",
    ] {
        let text = document(AnnotatedFormat::Embl, "DNA", "ACGT", "").replacen(
            "STD; SYN",
            &format!("{data_class}; PRO"),
            1,
        );
        let entry = parse::<Dna>(AnnotatedFormat::Embl, &text);
        for result in [
            entry.slice(1..3, options("LOCAL_part")).unwrap(),
            entry.reverse_complement(options("LOCAL_rc")).unwrap(),
        ] {
            assert!(
                result.record().metadata()[0]
                    .value
                    .contains(&format!("; {data_class}; PRO;"))
            );
            assert_raw_cache(result.record());
        }
    }
    let text = document(AnnotatedFormat::GenBank, "protein", "MELK", "").replacen("SYN", "PRI", 1);
    let protein = parse::<Protein>(AnnotatedFormat::GenBank, &text);
    let result = protein.slice(1..3, options("LOCAL_peptide")).unwrap();
    assert_eq!(
        result.record().metadata()[0].value,
        "LOCAL_peptide 2 aa linear PRI 16-OCT-2026"
    );
}

#[test]
fn missing_unknown_and_cross_format_classifications_are_not_inferred_even_with_discard() {
    for format in FORMATS {
        let original = document(format, "DNA", "ACGT", "");
        let invalid = match format {
            AnnotatedFormat::GenBank => vec![
                original.replacen(" SYN", "", 1),
                original.replacen(" SYN 15-OCT-2026", "", 1),
                original.replacen("SYN", "ZZZ", 1),
                original.replacen("SYN", "PRO", 1),
            ],
            AnnotatedFormat::Embl => vec![
                original.replacen("SYN", "", 1),
                original.replacen("SYN", "ZZZ", 1),
                original.replacen("SYN", "BCT", 1),
                original.replacen("STD", "", 1),
                original.replacen("STD", "ZZZ", 1),
                original.replacen("STD", "CON", 1),
            ],
        };
        for text in invalid {
            let entry = parse::<Dna>(format, &text);
            assert_canonical_round_trip(&entry);
            for policy in [AnnotationLossPolicy::Reject, AnnotationLossPolicy::Discard] {
                let selected = options("LOCAL_checked").with_loss_policy(policy);
                assert_eq!(
                    entry.slice(0..4, selected.clone()).unwrap_err().kind(),
                    IoErrorKind::Unsupported
                );
                assert_eq!(
                    entry.reverse_complement(selected).unwrap_err().kind(),
                    IoErrorKind::Unsupported
                );
            }
            assert_eq!(entry.raw_record(), text);
        }
    }
}

#[test]
fn opaque_sq_information_is_rejected_before_any_stream_or_path_output() {
    let format = AnnotatedFormat::Embl;
    let ordinary = document(format, "DNA", "aCgTrY", "");
    let directory = Directory::new();
    let previous = directory.path("previous.embl");
    fs::write(&previous, b"preserved previous content").unwrap();
    for payload in [
        "Sequence 6 BP; opaque trailing metadata",
        "Sequence 6 BP; 1 A;",
        "Sequence 6 BP; 1 A; 1 C; 1 G; 1 T; 2 other; extra",
    ] {
        let input = ordinary.replace("Sequence 6 BP;", payload);
        let entry = parse::<Dna>(format, &input);
        let mut output = Vec::new();
        assert_eq!(
            write(format, &mut output, &entry).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
        assert!(output.is_empty());
        assert_eq!(
            write_new(format, directory.path("absent.embl"), &entry)
                .unwrap_err()
                .kind(),
            IoErrorKind::Unsupported
        );
        assert_eq!(
            write_new(format, &previous, &entry).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
        assert_eq!(fs::read(&previous).unwrap(), b"preserved previous content");
        assert_eq!(directory.count(), 1);
    }
}

#[test]
fn quoted_empty_unicode_and_multiline_qualifiers_remain_distinct() {
    for format in FORMATS {
        let prefix = if format == AnnotatedFormat::GenBank {
            "                     "
        } else {
            "FT                   "
        };
        let mut annotations = feature(
            format,
            "misc_feature",
            "1..4",
            &["/note=\"\"", "/note=\"β and \"\"quotes\"\"\"", "/pseudo"],
        );
        annotations.push_str(&format!(
            "{prefix}/experiment=unquoted\n{prefix}continuation\n"
        ));
        let entry = parse::<Dna>(format, &document(format, "DNA", "aCgT", &annotations));
        let output = assert_canonical_round_trip(&entry);
        let parsed = parse::<Dna>(format, &output);
        let qualifiers = &parsed.features()[0].qualifiers;
        assert_eq!(qualifiers.len(), 4);
        assert_eq!(qualifiers[0].value.as_deref(), Some(""));
        assert_eq!(qualifiers[1].value.as_deref(), Some("β and \"quotes\""));
        assert_eq!(qualifiers[2].value, None);
        assert_eq!(
            qualifiers[3].value.as_deref(),
            Some("unquoted\ncontinuation")
        );
    }
}

#[test]
fn opaque_structural_payloads_require_discard_and_report_their_replacement() {
    let ordinary_record = document(AnnotatedFormat::GenBank, "DNA", "aCgTrY", "");
    let standard_entry = document(AnnotatedFormat::Embl, "DNA", "aCgTrY", "");
    for (format, input, key, sentinel) in [
        (
            AnnotatedFormat::GenBank,
            ordinary_record.replace("ORIGIN\n", "ORIGIN      unreviewed-origin-payload\n"),
            "ORIGIN",
            "unreviewed-origin-payload",
        ),
        (
            AnnotatedFormat::GenBank,
            ordinary_record.replace("15-OCT-2026\n", "15-OCT-2026 unreviewed-header-tail\n"),
            "LOCUS",
            "unreviewed-header-tail",
        ),
        (
            AnnotatedFormat::Embl,
            standard_entry.replace(
                "SQ   Sequence 6 BP;",
                "SQ   Sequence 6 BP; unreviewed-summary-tail",
            ),
            "SQ",
            "unreviewed-summary-tail",
        ),
    ] {
        let entry = parse::<Dna>(format, &input);
        let original = entry.clone();
        let (index, field) = entry
            .metadata()
            .iter()
            .enumerate()
            .find(|(_, field)| field.key == key)
            .unwrap();
        assert!(field.value.contains(sentinel));
        assert_eq!(
            entry
                .slice(0..6, options("LOCAL_reject"))
                .unwrap_err()
                .kind(),
            IoErrorKind::DataLoss
        );
        assert_eq!(
            entry
                .reverse_complement(options("LOCAL_reject"))
                .unwrap_err()
                .kind(),
            IoErrorKind::DataLoss
        );
        for result in [
            entry.slice(0..6, discard("LOCAL_clean")).unwrap(),
            entry.reverse_complement(discard("LOCAL_clean")).unwrap(),
        ] {
            assert!(!result.record().raw_record().contains(sentinel));
            let replacement = result
                .record()
                .metadata()
                .iter()
                .find(|field| field.key == key)
                .unwrap();
            assert!(result.report().changes.iter().any(|change| matches!(change,
                AnnotationChange::MetadataChanged { index: Some(position), key: changed_key, before: Some(before), after: Some(after), reason }
                if *position == index && changed_key == key && before == &field.value && after == &replacement.value && !reason.is_empty())));
            assert_raw_cache(result.record());
        }
        assert_eq!(entry, original);
    }
}

#[test]
fn protein_complement_locations_cannot_survive_or_be_discarded_by_slice() {
    let format = AnnotatedFormat::GenBank;
    for location in [
        "complement(1..2)",
        "join(1,complement(3))",
        "complement(complement(1..2))",
    ] {
        let text = document(
            format,
            "protein",
            "mKBZJX",
            &feature(format, "Region", location, &[]),
        );
        let entry = parse::<Protein>(format, &text);
        // Preserving the opaque-to-protein extraction expression in a same-
        // format document is distinct from claiming a meaningful transform.
        assert_canonical_round_trip(&entry);
        for policy in [AnnotationLossPolicy::Reject, AnnotationLossPolicy::Discard] {
            for range in [0..6, 4..6, 0..0] {
                assert_eq!(
                    entry
                        .slice(range, options("LOCAL_protein").with_loss_policy(policy))
                        .unwrap_err()
                        .kind(),
                    IoErrorKind::Unsupported
                );
            }
        }
    }
}

#[test]
fn generated_description_and_section_headers_have_addition_reports() {
    for format in FORMATS {
        let input = document(format, "DNA", "aCgT", "")
            .replace("DEFINITION  Original annotation.\n", "")
            .replace("DE   Original annotation.\n", "")
            .replace("FEATURES             Location/Qualifiers\n", "")
            .replace("FH   Key             Location/Qualifiers\nFH\n", "");
        let entry = parse::<Dna>(format, &input);
        assert!(entry.record().description().is_none());
        let result = entry.slice(1..3, options("LOCAL_added")).unwrap();
        let (description_key, section_key) = match format {
            AnnotatedFormat::GenBank => ("DEFINITION", "FEATURES"),
            AnnotatedFormat::Embl => ("DE", "FH"),
        };
        for key in [description_key, section_key] {
            assert!(result.report().changes.iter().any(|change| matches!(change,
                AnnotationChange::MetadataChanged { index: None, key: changed, before: None, after: Some(_), reason }
                if changed == key && !reason.is_empty())));
        }
        assert_raw_cache(result.record());
    }
}

#[test]
fn location_nesting_boundary_is_checked_and_valid_depth_can_be_transformed() {
    for format in FORMATS {
        let valid = format!("{}1..2{}", "complement(".repeat(64), ")".repeat(64));
        let input = document(format, "DNA", "aCgT", &feature(format, "gene", &valid, &[]));
        let entry = parse::<Dna>(format, &input);
        assert_canonical_round_trip(&entry);
        let result = entry.reverse_complement(options("LOCAL_depth")).unwrap();
        assert_eq!(
            result.record().features()[0]
                .extract(result.record().record().sequence())
                .unwrap()
                .to_string(),
            "AC"
        );
        let invalid = format!("complement({valid})");
        let input = document(
            format,
            "DNA",
            "aCgT",
            &feature(format, "gene", &invalid, &[]),
        );
        let error = match format {
            AnnotatedFormat::GenBank => parse_genbank::<Dna>(&input).unwrap_err(),
            AnnotatedFormat::Embl => parse_embl::<Dna>(&input).unwrap_err(),
        };
        assert_eq!(error.kind(), IoErrorKind::LimitExceeded);
    }
}

#[test]
fn tiny_canonical_serializations_match_independent_literal_outputs() {
    let expected_genbank = concat!(
        "LOCUS       ORIGINAL 4 bp DNA linear SYN 15-OCT-2026\n",
        "DEFINITION  Original annotation.\n",
        "FEATURES             Location/Qualifiers\n",
        "     gene            1\n",
        "                     /gene=\"example\"\n",
        "ORIGIN\n",
        "        1 aCgT\n",
        "//\n",
    );
    let expected_embl = concat!(
        "ID   ORIGINAL; SV 1; linear; DNA; STD; SYN; 4 BP.\n",
        "DE   Original annotation.\n",
        "FH   Key             Location/Qualifiers\n",
        "FH\n",
        "FT   gene            1\n",
        "FT                   /gene=\"example\"\n",
        "SQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n",
        "     aCgT",
        "          ",
        "          ",
        "          ",
        "          ",
        "          ",
        "          ",
        "          ",
        "4\n",
        "//\n",
    );
    for (format, expected) in [
        (AnnotatedFormat::GenBank, expected_genbank),
        (AnnotatedFormat::Embl, expected_embl),
    ] {
        let input = document(
            format,
            "DNA",
            "aCgT",
            &feature(format, "gene", "1..1", &["/gene=\"example\""]),
        );
        let record = parse::<Dna>(format, &input);
        assert_eq!(render(&record), expected);
    }
}
