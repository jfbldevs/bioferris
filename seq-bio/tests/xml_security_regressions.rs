//! Focused regressions for the quick-xml 0.41.0 security update and the shared
//! line reader's retry of transient `Interrupted` reads.
//!
//! Inputs are committed fixtures or small finite expansions of them; no hostile
//! unbounded XML is generated.

use std::fmt::Write as _;
use std::io::{self, BufReader, Read};

use bioferris::io::blast::{
    BlastTabularOptions, BlastXmlFormat, BlastXmlOptions, parse_blast_tabular, parse_blast_xml,
};
use bioferris::io::{FastaReader, FastqReader, IoErrorKind};
use bioferris::{Dna, Record};

const FASTA: &[u8] = b">sample mixed case\r\nAc\r\ngTN\r\n>second\nRY\n";
const FASTQ: &[u8] = b"@sample\nACGT\n+\n!I~+\n@second\nRY\n+\n!~\n";
const LEGACY_XML: &str = include_str!("fixtures/blast-example.xml");
const XML2: &str = include_str!("fixtures/blast-example-xml2.xml");
const TABULAR: &[u8] = include_bytes!("fixtures/blast-example.tab6");

/// Yields one byte per call and fails once with `Interrupted` at `interrupt_at`.
struct InterruptOnce {
    bytes: &'static [u8],
    position: usize,
    interrupt_at: usize,
    interrupted: bool,
}

impl Read for InterruptOnce {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if !self.interrupted && self.position == self.interrupt_at {
            self.interrupted = true;
            return Err(io::ErrorKind::Interrupted.into());
        }
        let Some(&byte) = self.bytes.get(self.position) else {
            return Ok(0);
        };
        let Some(slot) = destination.first_mut() else {
            return Ok(0);
        };
        *slot = byte;
        self.position += 1;
        Ok(1)
    }
}

fn interrupted(bytes: &'static [u8], interrupt_at: usize) -> BufReader<InterruptOnce> {
    BufReader::with_capacity(
        4,
        InterruptOnce {
            bytes,
            position: 0,
            interrupt_at,
            interrupted: false,
        },
    )
}

fn fasta(reader: impl io::BufRead) -> Result<Vec<Record<Dna>>, IoErrorKind> {
    FastaReader::<_, Dna>::new(reader)
        .collect::<Result<_, _>>()
        .map_err(|error| error.kind())
}

fn fastq(reader: impl io::BufRead) -> Result<Vec<Record<Dna>>, IoErrorKind> {
    FastqReader::<_, Dna>::new(reader)
        .collect::<Result<_, _>>()
        .map_err(|error| error.kind())
}

#[test]
fn interrupted_line_reads_retry_without_losing_a_partial_line() {
    let sequences = fasta(FASTA).unwrap();
    let reads = fastq(FASTQ).unwrap();
    assert_eq!(sequences.len(), 2);
    assert_eq!(reads.len(), 2);
    for position in 0..=FASTA.len() {
        assert_eq!(
            fasta(interrupted(FASTA, position)),
            Ok(sequences.clone()),
            "FASTA interruption at byte {position}"
        );
    }
    for position in 0..=FASTQ.len() {
        assert_eq!(
            fastq(interrupted(FASTQ, position)),
            Ok(reads.clone()),
            "FASTQ interruption at byte {position}"
        );
    }
    let expected = parse_blast_tabular(TABULAR, BlastTabularOptions::default()).unwrap();
    for position in 0..=TABULAR.len() {
        let observed = parse_blast_tabular(
            interrupted(TABULAR, position),
            BlastTabularOptions::default(),
        )
        .unwrap_or_else(|error| panic!("tabular interruption at {position}: {error}"));
        assert_eq!(observed.queries(), expected.queries());
    }
}

fn xml_options(format: BlastXmlFormat) -> BlastXmlOptions {
    BlastXmlOptions {
        format,
        ..BlastXmlOptions::default()
    }
}

/// Inserts `count` distinct prefixed namespace declarations into the root tag.
fn with_root_declarations(input: &str, root: &str, count: usize) -> String {
    let mut declarations = String::new();
    for index in 0..count {
        write!(declarations, " xmlns:n{index}=\"urn:n{index}\"").unwrap();
    }
    let start = input.find(root).unwrap() + root.len();
    let mut expanded = input.to_owned();
    expanded.insert_str(start, &declarations);
    expanded
}

#[test]
fn blast_xml_namespace_declarations_are_capped_at_256_per_element() {
    // XML2 already declares its default namespace on the root element.
    for (format, input, root, existing) in [
        (BlastXmlFormat::Legacy, LEGACY_XML, "<BlastOutput", 0),
        (BlastXmlFormat::Xml2, XML2, "<BlastXML2", 1),
    ] {
        let expected = parse_blast_xml(input.as_bytes(), xml_options(format)).unwrap();
        let accepted = with_root_declarations(input, root, 256 - existing);
        assert_eq!(
            parse_blast_xml(accepted.as_bytes(), xml_options(format))
                .unwrap()
                .queries(),
            expected.queries(),
            "{format:?}: 256 declarations"
        );
        let rejected = with_root_declarations(input, root, 257 - existing);
        let error = parse_blast_xml(rejected.as_bytes(), xml_options(format)).unwrap_err();
        assert_eq!(
            error.kind(),
            IoErrorKind::InvalidStructure,
            "{format:?}: {error}"
        );
    }
}

#[test]
fn blast_xml_duplicate_root_attributes_are_rejected() {
    for (format, input, root) in [
        (BlastXmlFormat::Legacy, LEGACY_XML, "<BlastOutput"),
        (BlastXmlFormat::Xml2, XML2, "<BlastXML2"),
    ] {
        let start = input.find(root).unwrap() + root.len();
        let mut duplicate = input.to_owned();
        duplicate.insert_str(start, " xmlns:n=\"urn:a\" xmlns:n=\"urn:b\"");
        let error = parse_blast_xml(duplicate.as_bytes(), xml_options(format)).unwrap_err();
        assert_eq!(
            error.kind(),
            IoErrorKind::InvalidStructure,
            "{format:?}: {error}"
        );
    }
}

#[test]
fn documented_namespace_limit_error_matches_the_library() {
    let guide = include_str!("../../docs/testing/security.md");
    let marker = "<!-- security-namespace-limit-output -->\n```text\n";
    assert_eq!(guide.matches(marker).count(), 1);
    let expected = guide
        .split_once(marker)
        .unwrap()
        .1
        .split_once("```")
        .unwrap()
        .0;
    let input = with_root_declarations(LEGACY_XML, "<BlastOutput", 257);
    let error = parse_blast_xml(input.as_bytes(), xml_options(BlastXmlFormat::Legacy)).unwrap_err();
    assert_eq!(
        (error.code(), error.record(), error.line()),
        ("invalid_structure", Some(1), Some(2))
    );
    assert_eq!(format!("kind: {:?}\n{error}\n", error.kind()), expected);
}
