//! Streaming GenBank/EMBL records: equality with the document readers, error
//! context and policy, limits, paths and the documented example.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::io::annotated::{
    AnnotatedFormat, AnnotatedReader, AnnotatedRecord, RecordErrorPolicy, embl, genbank,
    parse_embl, parse_genbank,
};
use bioferris::io::{BioIoError, IoErrorKind, ReadLimits};
use bioferris::prelude::Alphabet;
use bioferris::{Dna, Rna};
use flate2::{Compression, write::GzEncoder};

const PUBLIC_GENBANK: &str = include_str!("fixtures/annotated-dq067560.1.gb");
const PUBLIC_EMBL: &str = include_str!("fixtures/annotated-dq067560.1.embl");
const STREAM: &str = include_str!("fixtures/annotated-stream.gb");

fn genbank_record(name: &str, sequence: &str, molecule: &str) -> String {
    format!(
        "LOCUS       {name} {} bp {molecule} linear SYN 22-OCT-2026\nORIGIN\n        1 {sequence}\n//\n",
        sequence.len()
    )
}

fn embl_record(name: &str, sequence: &str, molecule: &str) -> String {
    let length = sequence.len();
    format!(
        "ID   {name}; SV 1; linear; {molecule}; STD; SYN; {length} BP.\nSQ   Sequence {length} BP;\n     {sequence} {length}\n//\n"
    )
}

fn stream<A: Alphabet>(
    input: &str,
    format: AnnotatedFormat,
    policy: RecordErrorPolicy,
) -> Vec<Result<AnnotatedRecord<A>, BioIoError>> {
    AnnotatedReader::new(input.as_bytes(), format)
        .with_record_errors(policy)
        .collect()
}

fn context(error: &BioIoError) -> (IoErrorKind, &str, Option<usize>, Option<usize>) {
    (error.kind(), error.format(), error.record(), error.line())
}

#[test]
fn streamed_records_equal_the_document_readers() {
    let second = genbank_record("SECOND", "acgt", "DNA").replace('\n', "\r\n");
    let third = genbank_record("THIRD", "ggcc", "DNA");
    let input = format!("\n{PUBLIC_GENBANK}\n\n{second}{}", third.trim_end());
    let streamed = stream::<Dna>(&input, AnnotatedFormat::GenBank, RecordErrorPolicy::Stop)
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(streamed, parse_genbank::<Dna>(&input).unwrap());
    assert_eq!(streamed.len(), 3);
    assert_eq!(streamed[1].raw_record(), second);
    assert_eq!(streamed[2].raw_record(), third.trim_end());

    let input = format!(
        "{PUBLIC_EMBL}{}\n{}",
        embl_record("SECOND", "acgt", "genomic DNA"),
        embl_record("THIRD", "ggcc", "genomic DNA")
    );
    let streamed = stream::<Dna>(&input, AnnotatedFormat::Embl, RecordErrorPolicy::Stop)
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(streamed, parse_embl::<Dna>(&input).unwrap());
    assert_eq!(streamed.len(), 3);

    for empty in ["", " \n\r\n"] {
        assert!(stream::<Dna>(empty, AnnotatedFormat::GenBank, RecordErrorPolicy::Stop).is_empty());
    }
}

#[test]
fn framing_errors_match_the_document_readers_and_always_end_iteration() {
    let complete = genbank_record("TEST", "atg", "DNA");
    let open = complete.replace("//\n", "");
    for input in [
        format!("{complete}{open}"),
        format!("{open}{complete}"),
        format!("{complete}>TEST\nATG\n"),
        format!("{complete}//\n"),
    ] {
        let expected = parse_genbank::<Dna>(&input).unwrap_err();
        for policy in [RecordErrorPolicy::Stop, RecordErrorPolicy::Continue] {
            let mut reader =
                AnnotatedReader::<_, Dna>::new(input.as_bytes(), AnnotatedFormat::GenBank)
                    .with_record_errors(policy);
            let mut error = None;
            for item in reader.by_ref() {
                if let Err(failure) = item {
                    error = Some(failure);
                    break;
                }
            }
            let error = error.expect("a framing error");
            assert_eq!(context(&error), context(&expected), "{input:?}");
            assert!(reader.next().is_none(), "{input:?}");
        }
    }
}

#[test]
fn record_errors_stop_by_default_and_continue_on_request() {
    let stopped = stream::<Dna>(STREAM, AnnotatedFormat::GenBank, RecordErrorPolicy::Stop);
    assert_eq!(stopped.len(), 2);
    assert_eq!(stopped[0].as_ref().unwrap().record().id(), "DEMO1.1");
    let error = stopped[1].as_ref().unwrap_err();
    assert_eq!(
        context(error),
        (IoErrorKind::InvalidSequence, "genbank", Some(2), Some(11))
    );

    let continued = stream::<Dna>(
        STREAM,
        AnnotatedFormat::GenBank,
        RecordErrorPolicy::Continue,
    );
    let outcomes = continued
        .iter()
        .map(|item| match item {
            Ok(entry) => entry.record().id().to_owned(),
            Err(error) => format!("error {}", error.record().unwrap()),
        })
        .collect::<Vec<_>>();
    assert_eq!(outcomes, ["DEMO1.1", "error 2", "DEMO3.1"]);

    // A second pass with the other molecule reads the record the first skipped.
    let rna = stream::<Rna>(
        STREAM,
        AnnotatedFormat::GenBank,
        RecordErrorPolicy::Continue,
    );
    let outcomes = rna
        .iter()
        .map(|item| match item {
            Ok(entry) => entry.record().id().to_owned(),
            Err(error) => format!("error {}", error.record().unwrap()),
        })
        .collect::<Vec<_>>();
    assert_eq!(outcomes, ["error 1", "DEMO2.1", "error 3"]);
}

#[test]
fn limits_bound_lines_totals_and_record_symbols() {
    let input = format!(
        "{}{}{}",
        genbank_record("LONG", "acgtacgta", "DNA"),
        genbank_record("SHORT", "acgt", "DNA"),
        genbank_record("LONGER", "acgtacgtacgt", "DNA")
    );
    let symbols = ReadLimits {
        max_record_symbols: 5,
        ..ReadLimits::default()
    };
    let outcomes =
        AnnotatedReader::<_, Dna>::with_limits(input.as_bytes(), AnnotatedFormat::GenBank, symbols)
            .with_record_errors(RecordErrorPolicy::Continue)
            .map(|item| match item {
                Ok(entry) => entry.record().id().to_owned(),
                Err(error) => format!("{} {}", error.code(), error.record().unwrap()),
            })
            .collect::<Vec<_>>();
    assert_eq!(outcomes, ["limit_exceeded 1", "SHORT", "limit_exceeded 3"]);

    let lines = ReadLimits {
        max_line_bytes: 40,
        ..ReadLimits::default()
    };
    let items =
        AnnotatedReader::<_, Dna>::with_limits(input.as_bytes(), AnnotatedFormat::GenBank, lines)
            .with_record_errors(RecordErrorPolicy::Continue)
            .collect::<Vec<_>>();
    assert_eq!(items.len(), 1);
    let error = items[0].as_ref().unwrap_err();
    assert_eq!(
        context(error),
        (IoErrorKind::LimitExceeded, "genbank", Some(1), Some(1))
    );

    let total = ReadLimits {
        max_total_bytes: input.find("LOCUS       LONGER").unwrap() + 10,
        ..ReadLimits::default()
    };
    let items =
        AnnotatedReader::<_, Dna>::with_limits(input.as_bytes(), AnnotatedFormat::GenBank, total)
            .with_record_errors(RecordErrorPolicy::Continue)
            .collect::<Vec<_>>();
    assert_eq!(items.len(), 3);
    assert!(items[0].is_ok() && items[1].is_ok());
    let error = items[2].as_ref().unwrap_err();
    assert_eq!(
        context(error),
        (IoErrorKind::LimitExceeded, "genbank", Some(3), Some(9))
    );
}

struct TempInput(PathBuf);

impl TempInput {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let serial = NEXT.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "seq-bio-annotated-stream-{}-{serial}.data",
            std::process::id()
        ));
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        file.write_all(bytes).unwrap();
        Self(path)
    }
}

impl Drop for TempInput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[test]
fn path_readers_stream_plain_and_gzip_and_report_late_corruption() {
    let plain = TempInput::new(PUBLIC_GENBANK.as_bytes());
    let records = genbank::<Dna>(&plain.0)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records, parse_genbank::<Dna>(PUBLIC_GENBANK).unwrap());

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(PUBLIC_EMBL.as_bytes()).unwrap();
    let mut compressed = encoder.finish().unwrap();
    let valid = TempInput::new(&compressed);
    let records = embl::<Dna>(&valid.0)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records, parse_embl::<Dna>(PUBLIC_EMBL).unwrap());

    // The CRC is checked only once the stream reaches its end, after the
    // complete record has been returned.
    let crc_index = compressed.len() - 8;
    compressed[crc_index] ^= 1;
    let corrupt = TempInput::new(&compressed);
    let items = embl::<Dna>(&corrupt.0).unwrap().collect::<Vec<_>>();
    assert_eq!(items.len(), 2);
    assert!(items[0].is_ok());
    assert_eq!(items[1].as_ref().unwrap_err().kind(), IoErrorKind::Io);

    let invalid_utf8 = TempInput::new(&[0xff]);
    let items = genbank::<Dna>(&invalid_utf8.0).unwrap().collect::<Vec<_>>();
    assert_eq!(items.len(), 1);
    assert_eq!(
        context(items[0].as_ref().unwrap_err()),
        (IoErrorKind::InvalidStructure, "genbank", Some(1), Some(1))
    );
    assert!(genbank::<Dna>(plain.0.with_extension("missing")).is_err());
}

fn marked<'a>(body: &'a str, name: &str) -> &'a str {
    let start = format!("<!-- {name} -->");
    assert_eq!(body.matches(&start).count(), 1, "unique marker: {name}");
    body.split_once(&start).unwrap().1
}

fn fence<'a>(body: &'a str, language: &str) -> &'a str {
    let marker = format!("```{language}\n");
    let (_, rest) = body.split_once(&marker).expect("missing code fence");
    let (content, _) = rest.split_once("```").expect("unclosed code fence");
    content
}

#[test]
fn documented_stream_example_uses_the_actual_input_program_and_output() {
    let guide = include_str!("../../docs/api/annotated-input-v1.md");
    assert_eq!(
        fence(marked(guide, "annotated-stream-input"), "text"),
        STREAM
    );
    let program = include_str!("../examples/stream_annotations.rs");
    let (_, program) = program
        .split_once("\n\n")
        .expect("example module documentation");
    assert_eq!(
        fence(marked(guide, "annotated-stream-program"), "rust").trim(),
        program.trim()
    );
    let command = fence(marked(guide, "annotated-stream-command"), "bash");
    assert_eq!(command, "cargo bio stream_annotations\n");

    // Run the real checked-in alias from the toolkit root, as documented.
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root");
    let output = Command::new(env!("CARGO"))
        .args(["bio", "stream_annotations", "--offline", "--locked"])
        .current_dir(root)
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    assert_eq!(
        output.stdout,
        fence(marked(guide, "annotated-stream-output"), "text").as_bytes()
    );
}
