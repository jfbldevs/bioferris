//! Independent record, byte-span, strict-schema and stale-source index checks.

use std::fs::{self, File, FileTimes, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::io::index::{
    IndexAlphabet, IndexError, IndexOptions, IndexedFormat, SequenceIndex, SidecarFormat,
    SourceCheck,
};
use bioferris::io::{FastaReader, FastqReader, ReadLimits};
use bioferris::prelude::*;
use flate2::{Compression, write::GzEncoder};
use proptest::prelude::*;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-index-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: impl AsRef<Path>) -> PathBuf {
        self.0.join(name)
    }

    fn write(&self, name: impl AsRef<Path>, bytes: &[u8]) -> PathBuf {
        let path = self.path(name);
        fs::write(&path, bytes).unwrap();
        path
    }

    fn manifest(&self, name: &str, value: &Value) -> PathBuf {
        self.write(name, &serde_json::to_vec(value).unwrap())
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

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn sequential<A: IndexAlphabet>(input: &[u8], format: IndexedFormat) -> Vec<Record<A>> {
    match format {
        IndexedFormat::Fasta => FastaReader::new(input).collect::<Result<_, _>>().unwrap(),
        IndexedFormat::FastqSanger => FastqReader::new(input).collect::<Result<_, _>>().unwrap(),
    }
}

fn build<A: IndexAlphabet>(source: &Path, format: IndexedFormat) -> SequenceIndex<A> {
    SequenceIndex::build(source, format, IndexOptions::default()).unwrap()
}

fn saved<A: IndexAlphabet>(index: &SequenceIndex<A>, path: &Path) -> Value {
    index.save_new(path).unwrap();
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

fn check_records<A: IndexAlphabet>(index: &mut SequenceIndex<A>, expected: &[Record<A>]) {
    assert_eq!(index.len(), expected.len());
    assert_eq!(index.is_empty(), expected.is_empty());
    assert_eq!(
        index.ids().collect::<Vec<_>>(),
        expected.iter().map(Record::id).collect::<Vec<_>>()
    );
    assert_eq!(index.ids().len(), expected.len());
    // Seek backwards, forwards, and to the same key repeatedly rather than
    // accidentally passing through sequential-only cursor behavior.
    for record in expected.iter().rev().chain(expected) {
        assert_eq!(index.get(record.id()).unwrap().as_ref(), Some(record));
    }
    for absent in ["", "missing-key", " ", "--json"] {
        assert!(index.get(absent).unwrap().is_none());
    }
}

fn check_round_trip<A: IndexAlphabet>(input: &[u8], format: IndexedFormat, expected: &[Record<A>]) {
    let directory = Directory::new();
    let source = directory.write("source.data", input);
    let sidecar = directory.path("index.json");
    assert_eq!(sequential::<A>(input, format), expected);
    let mut index = build::<A>(&source, format);
    assert_eq!(index.format(), format);
    assert_eq!(index.source_len(), u64::try_from(input.len()).unwrap());
    assert_eq!(index.source_sha256(), digest(input));
    check_records(&mut index, expected);
    let manifest = saved(&index, &sidecar);
    assert_eq!(manifest["format_id"], "bioferris-sequence-index");
    assert_eq!(manifest["schema_version"], 1);
    assert_eq!(manifest["molecule"], A::MOLECULE);
    assert_eq!(manifest["source_format"], format.name());
    assert_eq!(manifest["key_policy"], "record-id-first-token-v1");
    assert_eq!(manifest["source_len"], input.len());
    assert_eq!(manifest["source_sha256"], digest(input));
    let mut reopened =
        SequenceIndex::<A>::open(&source, &sidecar, IndexOptions::default()).unwrap();
    check_records(&mut reopened, expected);
    // The binary sidecar reopens to the same records under either source check.
    let binary = directory.path("index.bfidx");
    index.save_new_as(&binary, SidecarFormat::Binary).unwrap();
    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut reopened =
            SequenceIndex::<A>::open_with(&source, &binary, IndexOptions::default(), check)
                .unwrap();
        assert_eq!(reopened.format(), format);
        assert_eq!(reopened.source_sha256(), digest(input));
        check_records(&mut reopened, expected);
    }
    assert_eq!(fs::read(&source).unwrap(), input);
}

#[test]
fn defaults_and_supported_alphabet_identifiers_are_explicit() {
    let options = IndexOptions::default();
    assert_eq!(options.max_records, 100_000);
    assert_eq!(options.max_id_bytes, 16 * 1024 * 1024);
    assert_eq!(options.max_index_bytes, 64 * 1024 * 1024);
    assert_eq!(options.read_limits.max_line_bytes, 8 * 1024 * 1024);
    assert_eq!(options.read_limits.max_record_symbols, 100_000_000);
    assert_eq!(options.read_limits.max_total_bytes, 1024 * 1024 * 1024);
    assert_eq!(Dna::MOLECULE, "dna");
    assert_eq!(Rna::MOLECULE, "rna");
    assert_eq!(Protein::MOLECULE, "protein");
    assert_eq!(IndexedFormat::Fasta.name(), "fasta");
    assert_eq!(IndexedFormat::FastqSanger.name(), "fastq-sanger");
}

#[test]
fn empty_files_and_fasta_blank_only_files_index_without_inventing_records() {
    for format in [IndexedFormat::Fasta, IndexedFormat::FastqSanger] {
        check_round_trip::<Dna>(b"", format, &[]);
        check_round_trip::<Rna>(b"", format, &[]);
        check_round_trip::<Protein>(b"", format, &[]);
    }
    check_round_trip::<Dna>(b"\n \t\r\n\r\n", IndexedFormat::Fasta, &[]);
}

#[test]
fn headers_case_descriptions_and_qualities_equal_sequential_records() {
    let fasta = "\r\n>  α\tfirst description  \r\na C\tG\r\ntN\r\n\r\n>β empty\r\n>third\r\nrysw";
    let dna_expected = [
        Record::<Dna>::from_ascii("  α\tfirst description  ", "aCGtN", None).unwrap(),
        Record::from_ascii("β empty", "", None).unwrap(),
        Record::from_ascii("third", "rysw", None).unwrap(),
    ];
    check_round_trip(fasta.as_bytes(), IndexedFormat::Fasta, &dna_expected);
    let rna_input = fasta.replace('t', "u");
    let rna_expected = sequential::<Rna>(rna_input.as_bytes(), IndexedFormat::Fasta);
    check_round_trip(rna_input.as_bytes(), IndexedFormat::Fasta, &rna_expected);
    let protein_input = ">protein\tpreserved metadata\nmK\nBZJXUO\n>empty\n";
    let protein_expected = [
        Record::<Protein>::from_ascii("protein\tpreserved metadata", "mKBZJXUO", None).unwrap(),
        Record::from_ascii("empty", "", None).unwrap(),
    ];
    check_round_trip(
        protein_input.as_bytes(),
        IndexedFormat::Fasta,
        &protein_expected,
    );

    let reads_input = "@  read\tdescription  \r\naC\r\ngT\r\n+  read\tdescription  \r\n@\r\n+!~\r\n@second\r\nN\r\n+\r\nI";
    let expected = [
        Record::<Dna>::from_ascii("  read\tdescription  ", "aCgT", Some(vec![31, 10, 0, 93]))
            .unwrap(),
        Record::from_ascii("second", "N", Some(vec![40])).unwrap(),
    ];
    check_round_trip(
        reads_input.as_bytes(),
        IndexedFormat::FastqSanger,
        &expected,
    );
    let rna_input = reads_input.replace('T', "U");
    let rna_expected = sequential::<Rna>(rna_input.as_bytes(), IndexedFormat::FastqSanger);
    check_round_trip(
        rna_input.as_bytes(),
        IndexedFormat::FastqSanger,
        &rna_expected,
    );
    let protein_input = b"@peptide\nmKBZ\n+peptide\n@+!~\n";
    let protein_expected =
        [Record::<Protein>::from_ascii("peptide", "mKBZ", Some(vec![31, 10, 0, 93])).unwrap()];
    check_round_trip(protein_input, IndexedFormat::FastqSanger, &protein_expected);
}

fn append_line(output: &mut Vec<u8>, line: &[u8], crlf: bool) {
    output.extend_from_slice(line);
    output.extend_from_slice(if crlf { b"\r\n" } else { b"\n" });
}

fn generated<A: IndexAlphabet>(
    raw: &[Vec<u8>],
    format: IndexedFormat,
    wrapping: usize,
    crlf: bool,
    final_newline: bool,
) -> (Vec<u8>, Vec<Record<A>>) {
    let mut bytes = Vec::new();
    let mut records = Vec::new();
    if format == IndexedFormat::Fasta {
        append_line(&mut bytes, b" \t", crlf);
    }
    for (ordinal, raw) in raw.iter().enumerate() {
        let header = format!("  gene_{ordinal}\tannotation α {ordinal}  ");
        let sequence: Vec<_> = raw
            .iter()
            .map(|byte| {
                let symbol = A::SYMBOLS[usize::from(*byte) % A::SYMBOLS.len()];
                if byte & 1 == 0 {
                    symbol
                } else {
                    symbol.to_ascii_lowercase()
                }
            })
            .collect();
        let quality: Vec<_> = raw.iter().map(|byte| byte % 94).collect();
        let prefix = if format == IndexedFormat::Fasta {
            '>'
        } else {
            '@'
        };
        append_line(&mut bytes, format!("{prefix}{header}").as_bytes(), crlf);
        let width = (wrapping + ordinal) % 13 + 1;
        for line in sequence.chunks(width) {
            if format == IndexedFormat::Fasta {
                let spaced: Vec<_> = line.iter().flat_map(|byte| [*byte, b' ', b'\t']).collect();
                append_line(&mut bytes, &spaced, crlf);
            } else {
                append_line(&mut bytes, line, crlf);
            }
        }
        if format == IndexedFormat::Fasta {
            append_line(&mut bytes, b"", crlf);
        } else {
            let repeated = if ordinal % 2 == 0 {
                format!("+{header}")
            } else {
                "+".to_owned()
            };
            append_line(&mut bytes, repeated.as_bytes(), crlf);
            for line in quality.chunks((width + 4) % 11 + 1) {
                let encoded: Vec<_> = line.iter().map(|score| score + 33).collect();
                append_line(&mut bytes, &encoded, crlf);
            }
        }
        records.push(
            Record::from_ascii(
                header,
                sequence,
                (format == IndexedFormat::FastqSanger).then_some(quality),
            )
            .unwrap(),
        );
    }
    if !final_newline && bytes.ends_with(b"\n") {
        bytes.truncate(bytes.len() - if crlf { 2 } else { 1 });
    }
    (bytes, records)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn generated_variable_wrapping_preserves_all_alphabets(
        raw in prop::collection::vec(prop::collection::vec(any::<u8>(), 1..65), 0..6),
        wrapping in 0_usize..40,
        crlf in any::<bool>(),
        final_newline in any::<bool>(),
    ) {
        for format in [IndexedFormat::Fasta, IndexedFormat::FastqSanger] {
            let (input, expected) = generated::<Dna>(&raw, format, wrapping, crlf, final_newline);
            check_round_trip(&input, format, &expected);
            let (input, expected) = generated::<Rna>(&raw, format, wrapping, crlf, final_newline);
            check_round_trip(&input, format, &expected);
            let (input, expected) = generated::<Protein>(&raw, format, wrapping, crlf, final_newline);
            check_round_trip(&input, format, &expected);
        }
    }
}

#[test]
fn persisted_fasta_spans_follow_literal_raw_boundaries_not_buffer_read_ahead() {
    let directory = Directory::new();
    let prefix = b"\r\n \t\n";
    let first = b">alpha full\r\na C\r\ngT\n\n";
    let second = ">β\tsecond\nNN\n".as_bytes();
    let third = b">empty";
    let input = [prefix.as_slice(), first, second, third].concat();
    let source = directory.write("source.fasta", &input);
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    let value = saved(&index, &directory.path("index.json"));
    let mut offset = prefix.len();
    for (entry, (raw, id, length, line, record)) in
        value["entries"].as_array().unwrap().iter().zip([
            (first.as_slice(), "alpha", 4, 3, 1),
            (second, "β", 2, 7, 2),
            (third.as_slice(), "empty", 0, 9, 3),
        ])
    {
        assert_eq!(
            *entry,
            json!({"id":id,"offset":offset,"byte_length":raw.len(),
            "sequence_length":length,"record":record,"header_line":line,"sha256":digest(raw)})
        );
        offset += raw.len();
    }
    assert_eq!(offset, input.len());
    assert_eq!(value["entries"].as_array().unwrap().len(), 3);
}

#[test]
fn persisted_fastq_spans_count_wrapped_quality_bytes_and_crlf_exactly() {
    let directory = Directory::new();
    let first = b"@a first\r\naC\r\ngT\r\n+a first\r\n@+\r\n!~\r\n";
    let second = b"@b\nN\n+\n+";
    let input = [first.as_slice(), second].concat();
    let source = directory.write("source.fastq", &input);
    let index = build::<Dna>(&source, IndexedFormat::FastqSanger);
    let value = saved(&index, &directory.path("index.json"));
    assert_eq!(
        value["entries"],
        json!([
            {"id":"a","offset":0,"byte_length":first.len(),"sequence_length":4,
                "record":1,"header_line":1,"sha256":digest(first)},
            {"id":"b","offset":first.len(),"byte_length":second.len(),"sequence_length":1,
                "record":2,"header_line":7,"sha256":digest(second)},
        ])
    );
}

#[test]
fn keys_are_first_tokens_case_sensitive_and_source_ordered() {
    let input = b">z last alphabetically\nA\n>A uppercase\nC\n>a lowercase\nG\n";
    let directory = Directory::new();
    let source = directory.write("source.fasta", input);
    let mut index = build::<Dna>(&source, IndexedFormat::Fasta);
    assert_eq!(index.ids().collect::<Vec<_>>(), ["z", "A", "a"]);
    assert_eq!(index.get("A").unwrap().unwrap().source_sequence(), "C");
    assert_eq!(index.get("a").unwrap().unwrap().source_sequence(), "G");
    assert!(index.get("z last alphabetically").unwrap().is_none());
    assert!(index.get(" z").unwrap().is_none());
}

#[test]
fn duplicate_first_tokens_are_rejected_instead_of_replaced() {
    for (format, input) in [
        (
            IndexedFormat::Fasta,
            ">same first\nA\n> other\nC\n>  same\tsecond  \nG\n",
        ),
        (
            IndexedFormat::FastqSanger,
            "@same first\nA\n+\n!\n@other\nC\n+\n+\n@ same second\nG\n+\n~\n",
        ),
    ] {
        let directory = Directory::new();
        let source = directory.write("source", input.as_bytes());
        let error = SequenceIndex::<Dna>::build(&source, format, IndexOptions::default())
            .err()
            .unwrap();
        assert_eq!(error.code(), "duplicate_record_id");
        assert!(
            matches!(error, IndexError::DuplicateId {id, first_record:1, duplicate_record:3} if id == "same")
        );
        assert_eq!(directory.count(), 1);
    }
}

#[test]
fn real_phix174_fixture_matches_sequential_records_and_table_11_translation() {
    let input = include_bytes!("fixtures/phix174-e.fasta");
    let expected = sequential::<Dna>(input, IndexedFormat::Fasta);
    assert_eq!(expected.len(), 1);
    assert_eq!(expected[0].sequence().len(), 276);
    check_round_trip(input, IndexedFormat::Fasta, &expected);
    let directory = Directory::new();
    let source = directory.write("phix174.fasta", input);
    let mut index = build::<Dna>(&source, IndexedFormat::Fasta);
    let record = index.get(expected[0].id()).unwrap().unwrap();
    let protein = record
        .sequence()
        .translate_with(TranslationConfig {
            genetic_code: GeneticCode::BacterialArchaealPlastid,
            mode: TranslationMode::Cds,
            ..TranslationConfig::default()
        })
        .unwrap();
    assert_eq!(protein.len(), 91);
    assert_eq!(record, expected[0]);
}

fn assert_limit(error: &IndexError, resource: &str, limit: usize, observed: usize) {
    assert_eq!(error.code(), "limit_exceeded");
    let IndexError::LimitExceeded {
        resource: actual,
        limit: bound,
        observed: count,
    } = error
    else {
        panic!("expected an index limit error: {error}");
    };
    assert_eq!(*actual, resource);
    assert_eq!(*bound, u64::try_from(limit).unwrap());
    assert_eq!(*count, u64::try_from(observed).unwrap());
}

#[test]
fn entry_and_utf8_id_budgets_are_aggregate_inclusive_and_checked_on_open() {
    let directory = Directory::new();
    let source = directory.write(
        "source",
        ">a description\nA\n>β long metadata ignored by ID budget\nC\n".as_bytes(),
    );
    let sidecar = directory.path("index.json");
    let index = SequenceIndex::<Dna>::build(
        &source,
        IndexedFormat::Fasta,
        IndexOptions {
            max_records: 2,
            max_id_bytes: 3,
            ..IndexOptions::default()
        },
    )
    .unwrap();
    index.save_new(&sidecar).unwrap();
    SequenceIndex::<Dna>::open(
        &source,
        &sidecar,
        IndexOptions {
            max_records: 2,
            max_id_bytes: 3,
            ..IndexOptions::default()
        },
    )
    .unwrap();
    for (resource, records, id_bytes, limit, observed) in [
        ("index_records", 1, 3, 1, 2),
        ("index_id_bytes", 2, 2, 2, 3),
        ("index_records", 0, 3, 0, 1),
        ("index_id_bytes", 2, 0, 0, 1),
    ] {
        let options = || IndexOptions {
            max_records: records,
            max_id_bytes: id_bytes,
            ..IndexOptions::default()
        };
        let error = SequenceIndex::<Dna>::build(&source, IndexedFormat::Fasta, options())
            .err()
            .unwrap();
        assert_limit(&error, resource, limit, observed);
        let error = SequenceIndex::<Dna>::open(&source, &sidecar, options())
            .err()
            .unwrap();
        assert_eq!(error.code(), "limit_exceeded");
    }
}

#[test]
fn zero_entry_and_id_budgets_allow_only_empty_collections() {
    let directory = Directory::new();
    let source = directory.write("empty", b"");
    let sidecar = directory.path("empty.json");
    let options = || IndexOptions {
        max_records: 0,
        max_id_bytes: 0,
        read_limits: ReadLimits {
            max_line_bytes: 0,
            max_record_symbols: 0,
            max_total_bytes: 0,
        },
        ..IndexOptions::default()
    };
    let index = SequenceIndex::<Dna>::build(&source, IndexedFormat::Fasta, options()).unwrap();
    assert!(index.is_empty());
    assert_eq!(
        index.source_sha256(),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    index.save_new(&sidecar).unwrap();
    assert!(
        SequenceIndex::<Dna>::open(&source, &sidecar, options())
            .unwrap()
            .is_empty()
    );
    let empty_record = directory.write("empty-record", b">one\n");
    assert_eq!(
        SequenceIndex::<Dna>::build(
            &empty_record,
            IndexedFormat::Fasta,
            IndexOptions {
                max_records: 0,
                ..IndexOptions::default()
            }
        )
        .err()
        .unwrap()
        .code(),
        "limit_exceeded"
    );
}

#[test]
fn serialized_byte_budgets_include_the_entire_compact_sidecar() {
    let directory = Directory::new();
    let source = directory.write("source", b">one\naCgT\n");
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    let sidecar = directory.path("index.json");
    let value = saved(&index, &sidecar);
    let bytes = fs::read(&sidecar).unwrap();
    let exact = bytes.len();
    assert!(!bytes.ends_with(b"\n"));
    let exact_index = SequenceIndex::<Dna>::build(
        &source,
        IndexedFormat::Fasta,
        IndexOptions {
            max_index_bytes: exact,
            ..IndexOptions::default()
        },
    )
    .unwrap();
    exact_index.save_new(directory.path("exact.json")).unwrap();
    assert_eq!(fs::read(directory.path("exact.json")).unwrap(), bytes);
    SequenceIndex::<Dna>::open(
        &source,
        &sidecar,
        IndexOptions {
            max_index_bytes: exact,
            ..IndexOptions::default()
        },
    )
    .unwrap();
    let mut empty_manifest = value.clone();
    empty_manifest["entries"] = json!([]);
    let empty_manifest_bytes = serde_json::to_vec(&empty_manifest).unwrap().len();
    // A zero limit fails before records are parsed; observed bytes describe the
    // attempted empty manifest, not a future index whose records are not read.
    for (limit, observed) in [(0, empty_manifest_bytes), (exact - 1, exact)] {
        let error = SequenceIndex::<Dna>::build(
            &source,
            IndexedFormat::Fasta,
            IndexOptions {
                max_index_bytes: limit,
                ..IndexOptions::default()
            },
        )
        .err()
        .unwrap();
        assert_limit(&error, "index_bytes", limit, observed);
        assert_eq!(
            SequenceIndex::<Dna>::open(
                &source,
                &sidecar,
                IndexOptions {
                    max_index_bytes: limit,
                    ..IndexOptions::default()
                }
            )
            .err()
            .unwrap()
            .code(),
            "limit_exceeded"
        );
    }
    let pretty = serde_json::to_vec_pretty(&value).unwrap();
    assert!(pretty.len() > exact);
    let pretty_path = directory.write("pretty.json", &pretty);
    SequenceIndex::<Dna>::open(
        &source,
        &pretty_path,
        IndexOptions {
            max_index_bytes: pretty.len(),
            ..IndexOptions::default()
        },
    )
    .unwrap();
    assert_eq!(
        SequenceIndex::<Dna>::open(
            &source,
            &pretty_path,
            IndexOptions {
                max_index_bytes: pretty.len() - 1,
                ..IndexOptions::default()
            }
        )
        .err()
        .unwrap()
        .code(),
        "limit_exceeded"
    );
}

#[test]
fn existing_line_record_and_total_byte_limits_remain_effective() {
    let directory = Directory::new();
    let input = b">a\r\naC\r\n";
    let source = directory.write("source", input);
    let exact = ReadLimits {
        max_line_bytes: 4,
        max_record_symbols: 2,
        max_total_bytes: input.len(),
    };
    let index = SequenceIndex::<Dna>::build(
        &source,
        IndexedFormat::Fasta,
        IndexOptions {
            read_limits: exact,
            ..IndexOptions::default()
        },
    )
    .unwrap();
    let sidecar = directory.path("index.json");
    index.save_new(&sidecar).unwrap();
    for limits in [
        ReadLimits {
            max_line_bytes: 3,
            ..exact
        },
        ReadLimits {
            max_record_symbols: 1,
            ..exact
        },
        ReadLimits {
            max_total_bytes: input.len() - 1,
            ..exact
        },
    ] {
        assert_eq!(
            SequenceIndex::<Dna>::build(
                &source,
                IndexedFormat::Fasta,
                IndexOptions {
                    read_limits: limits,
                    ..IndexOptions::default()
                }
            )
            .err()
            .unwrap()
            .code(),
            "limit_exceeded"
        );
        assert_eq!(
            SequenceIndex::<Dna>::open(
                &source,
                &sidecar,
                IndexOptions {
                    read_limits: limits,
                    ..IndexOptions::default()
                }
            )
            .err()
            .unwrap()
            .code(),
            "limit_exceeded"
        );
    }
}

#[test]
fn late_invalid_records_propagate_existing_reader_errors() {
    for (format, input, code) in [
        (
            IndexedFormat::Fasta,
            ">good\nAC\n>bad\nAZ\n",
            "invalid_sequence",
        ),
        (
            IndexedFormat::Fasta,
            ">good\nAC\n>  \nA\n",
            "invalid_header",
        ),
        (
            IndexedFormat::FastqSanger,
            "@good\nAC\n+\n!+\n@bad\nT\n+\n ",
            "invalid_quality",
        ),
        (
            IndexedFormat::FastqSanger,
            "@good\nAC\n+\n!+\n@bad\nAC\n+\n!",
            "truncated_input",
        ),
        (
            IndexedFormat::FastqSanger,
            "@good\nAC\n+\n!+\n@bad\nA\n+wrong\n!",
            "invalid_header",
        ),
        (
            IndexedFormat::FastqSanger,
            "@empty\n+\n",
            "unsupported_format_feature",
        ),
    ] {
        let directory = Directory::new();
        let source = directory.write("source", input.as_bytes());
        let error = SequenceIndex::<Dna>::build(&source, format, IndexOptions::default())
            .err()
            .unwrap();
        assert_eq!(error.code(), code);
        assert!(matches!(error, IndexError::Input(_)));
        assert_eq!(directory.count(), 1);
    }
}

#[test]
fn every_sanger_quality_value_survives_indexed_retrieval() {
    let mut input = b"@qualities exact range\n".to_vec();
    input.extend(std::iter::repeat_n(b'A', 94));
    input.extend_from_slice(b"\n+\n");
    input.extend(33..=126);
    let expected = [Record::<Dna>::from_ascii(
        "qualities exact range",
        vec![b'A'; 94],
        Some((0..=93).collect()),
    )
    .unwrap()];
    check_round_trip(&input, IndexedFormat::FastqSanger, &expected);
}

#[test]
fn malformed_sidecar_schema_is_rejected_without_trusting_source_digest() {
    let directory = Directory::new();
    let source = directory.write("source", b">a\nAC\n>b\nGT\n");
    let original = saved(
        &build::<Dna>(&source, IndexedFormat::Fasta),
        &directory.path("valid.json"),
    );
    let mutations = [
        ("schema_version", json!(2)),
        ("schema_version", json!("1")),
        ("format_id", json!("other-index")),
        ("key_policy", json!("full-header")),
        ("molecule", json!("rna")),
        ("source_format", json!("gzip-fasta")),
        ("source_len", json!(-1)),
        ("source_len", json!(1.5)),
        ("source_sha256", json!("not-a-sha256")),
        ("source_sha256", json!("A".repeat(64))),
        ("source_path", json!("/untrusted/manifest/path")),
        ("entries", json!(null)),
    ];
    for (ordinal, (field, value)) in mutations.into_iter().enumerate() {
        let mut changed = original.clone();
        changed[field] = value;
        let path = directory.manifest(&format!("bad-{ordinal}.json"), &changed);
        assert_eq!(
            SequenceIndex::<Dna>::open(&source, path, IndexOptions::default())
                .err()
                .unwrap()
                .code(),
            "invalid_index",
            "field {field}"
        );
    }
    for field in [
        "format_id",
        "schema_version",
        "molecule",
        "source_format",
        "key_policy",
        "source_len",
        "source_sha256",
        "entries",
    ] {
        let mut changed = original.clone();
        changed.as_object_mut().unwrap().remove(field);
        let path = directory.manifest(&format!("missing-{field}.json"), &changed);
        assert_eq!(
            SequenceIndex::<Dna>::open(&source, path, IndexOptions::default())
                .err()
                .unwrap()
                .code(),
            "invalid_index"
        );
    }
    let compact = original.to_string();
    let duplicate_top = format!("{{\"schema_version\":1,{}", &compact[1..]);
    let duplicate_entry = compact.replacen("\"id\":\"a\"", "\"id\":\"a\",\"id\":\"a\"", 1);
    // A derived Serde struct can accept a positional sequence even though the
    // persisted contract requires named JSON object fields at both levels.
    let positional_manifest = json!([
        original["format_id"],
        original["schema_version"],
        original["molecule"],
        original["source_format"],
        original["key_policy"],
        original["source_len"],
        original["source_sha256"],
        original["entries"]
    ]);
    let mut positional_entry = original.clone();
    let entry = &original["entries"][0];
    positional_entry["entries"][0] = json!([
        entry["id"],
        entry["offset"],
        entry["byte_length"],
        entry["sequence_length"],
        entry["record"],
        entry["header_line"],
        entry["sha256"]
    ]);
    for (ordinal, raw) in [
        duplicate_top,
        duplicate_entry,
        positional_manifest.to_string(),
        positional_entry.to_string(),
        format!("{compact} trailing"),
        "{".to_owned(),
    ]
    .iter()
    .enumerate()
    {
        let path = directory.write(format!("raw-{ordinal}.json"), raw.as_bytes());
        assert_eq!(
            SequenceIndex::<Dna>::open(&source, path, IndexOptions::default())
                .err()
                .unwrap()
                .code(),
            "invalid_index"
        );
    }
}

#[test]
fn forged_offsets_lengths_ids_and_metadata_fail_even_with_matching_source_hash() {
    let directory = Directory::new();
    let input = b">a\nAC\nGT\n>b\nTT\n";
    let source = directory.write("source", input);
    let original = saved(
        &build::<Dna>(&source, IndexedFormat::Fasta),
        &directory.path("valid.json"),
    );
    let mut variants = Vec::new();
    for (field, value) in [
        ("id", json!("forged")),
        ("id", json!("")),
        ("id", json!("two tokens")),
        ("offset", json!(u64::MAX)),
        ("offset", json!(-1)),
        ("byte_length", json!(u64::MAX)),
        ("byte_length", json!(0)),
        ("sequence_length", json!(3)),
        ("record", json!(0)),
        ("record", json!(2)),
        ("header_line", json!(0)),
        ("header_line", json!(2)),
        ("sha256", json!("0".repeat(64))),
        ("unknown", json!(true)),
    ] {
        let mut changed = original.clone();
        changed["entries"][0][field] = value;
        variants.push(changed);
    }
    let mut truncated = original.clone();
    truncated["entries"][0]["byte_length"] = json!(6);
    truncated["entries"][0]["sequence_length"] = json!(2);
    truncated["entries"][0]["sha256"] = json!(digest(&input[..6]));
    variants.push(truncated);
    let mut shifted = original.clone();
    shifted["entries"][0]["offset"] = json!(1);
    shifted["entries"][0]["byte_length"] = json!(8);
    shifted["entries"][0]["sha256"] = json!(digest(&input[1..9]));
    variants.push(shifted);
    let mut omitted = original.clone();
    omitted["entries"].as_array_mut().unwrap().pop();
    variants.push(omitted);
    let mut empty = original.clone();
    empty["entries"] = json!([]);
    variants.push(empty);
    let mut reordered = original.clone();
    reordered["entries"].as_array_mut().unwrap().reverse();
    variants.push(reordered);
    let mut duplicate = original.clone();
    duplicate["entries"][1]["id"] = json!("a");
    variants.push(duplicate);
    for (ordinal, changed) in variants.iter().enumerate() {
        assert_eq!(changed["source_sha256"], original["source_sha256"]);
        let path = directory.manifest(&format!("forged-{ordinal}.json"), changed);
        let error = SequenceIndex::<Dna>::open(&source, path, IndexOptions::default())
            .err()
            .unwrap();
        assert_eq!(error.code(), "invalid_index", "forgery {ordinal}: {error}");
    }
}

#[test]
fn opening_a_different_molecule_or_source_never_silently_reuses_an_index() {
    let directory = Directory::new();
    let source = directory.write("source", b">a\nACGT\n");
    let sidecar = directory.path("index.json");
    let original = saved(&build::<Dna>(&source, IndexedFormat::Fasta), &sidecar);
    assert_eq!(
        SequenceIndex::<Rna>::open(&source, &sidecar, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "invalid_index"
    );
    assert_eq!(
        SequenceIndex::<Protein>::open(&source, &sidecar, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "invalid_index"
    );
    let another = directory.write("another", b">a\nTGCA\n");
    assert_eq!(
        SequenceIndex::<Dna>::open(&another, &sidecar, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "stale_source"
    );
    let mut wrong_digest = original;
    wrong_digest["source_sha256"] = json!("0".repeat(64));
    let corrupted = directory.manifest("wrong-digest.json", &wrong_digest);
    assert_eq!(
        SequenceIndex::<Dna>::open(&source, corrupted, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "stale_source"
    );
}

#[test]
fn append_truncate_and_same_size_source_changes_invalidate_saved_indexes() {
    let original = b">a\nACGT\n>b\nTT\n";
    let mut appended = original.to_vec();
    appended.extend_from_slice(b">c\nN\n");
    let mut changed_case = original.to_vec();
    changed_case[3] = b'a';
    let mut changed_header = original.to_vec();
    changed_header[1] = b'z';
    for changed in [
        appended,
        original[..original.len() - 1].to_vec(),
        changed_case,
        changed_header,
    ] {
        let directory = Directory::new();
        let source = directory.write("source", original);
        let original_time = fs::metadata(&source).unwrap().modified().unwrap();
        let sidecar = directory.path("index.json");
        build::<Dna>(&source, IndexedFormat::Fasta)
            .save_new(&sidecar)
            .unwrap();
        fs::write(&source, &changed).unwrap();
        // Restoring the timestamp rules out size/mtime-only verification, even
        // for a case change that leaves the normalized biological sequence equal.
        File::options()
            .write(true)
            .open(&source)
            .unwrap()
            .set_times(FileTimes::new().set_modified(original_time))
            .unwrap();
        assert_eq!(
            fs::metadata(&source).unwrap().modified().unwrap(),
            original_time
        );
        assert_ne!(digest(original), digest(&changed));
        let error = SequenceIndex::<Dna>::open(&source, &sidecar, IndexOptions::default())
            .err()
            .unwrap();
        assert_eq!(error.code(), "stale_source");
    }
}

#[test]
fn quality_only_changes_are_detected_without_losing_valid_fastq_structure() {
    let directory = Directory::new();
    let source = directory.write("source.fastq", b"@a\naCgT\n+\n!@+~\n");
    let sidecar = directory.path("index.json");
    let index = build::<Dna>(&source, IndexedFormat::FastqSanger);
    index.save_new(&sidecar).unwrap();
    fs::write(&source, b"@a\naCgT\n+\n!A+~\n").unwrap();
    assert_eq!(
        sequential::<Dna>(&fs::read(&source).unwrap(), IndexedFormat::FastqSanger)[0]
            .source_sequence(),
        "aCgT"
    );
    assert_eq!(
        SequenceIndex::<Dna>::open(&source, &sidecar, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "stale_source"
    );
}

#[test]
fn changes_after_open_are_errors_even_for_missing_or_unchanged_selected_records() {
    let directory = Directory::new();
    let source = directory.write("source", b">a\nAC\n>b\nGT\n");
    let sidecar = directory.path("index.json");
    build::<Dna>(&source, IndexedFormat::Fasta)
        .save_new(&sidecar)
        .unwrap();
    let mut index = SequenceIndex::<Dna>::open(&source, &sidecar, IndexOptions::default()).unwrap();
    assert_eq!(index.get("a").unwrap().unwrap().source_sequence(), "AC");
    OpenOptions::new()
        .append(true)
        .open(&source)
        .unwrap()
        .write_all(b">c\nN\n")
        .unwrap();
    assert_eq!(index.get("a").unwrap_err().code(), "stale_source");
    assert_eq!(index.get("not-present").unwrap_err().code(), "stale_source");
}

#[test]
fn same_size_selected_span_mutation_after_open_cannot_pass_as_a_cached_record() {
    let directory = Directory::new();
    let source = directory.write("source", b">a\nACGT\n");
    let original_time = fs::metadata(&source).unwrap().modified().unwrap();
    let mut index = build::<Dna>(&source, IndexedFormat::Fasta);
    let mut writer = OpenOptions::new().write(true).open(&source).unwrap();
    writer.seek(SeekFrom::Start(3)).unwrap();
    writer.write_all(b"T").unwrap();
    writer.flush().unwrap();
    writer
        .set_times(FileTimes::new().set_modified(original_time))
        .unwrap();
    assert_eq!(index.get("a").unwrap_err().code(), "stale_source");
}

#[cfg(unix)]
#[test]
fn replacement_and_unlink_of_the_original_path_invalidate_the_retained_handle() {
    let directory = Directory::new();
    let bytes = b">a\nACGT\n";
    let source = directory.write("source", bytes);
    let mut index = build::<Dna>(&source, IndexedFormat::Fasta);
    fs::rename(&source, directory.path("old-source")).unwrap();
    fs::write(&source, bytes).unwrap();
    assert_eq!(index.get("a").unwrap_err().code(), "stale_source");
    let mut second = build::<Dna>(&source, IndexedFormat::Fasta);
    fs::remove_file(&source).unwrap();
    assert_eq!(second.get("a").unwrap_err().code(), "stale_source");
    assert_eq!(second.get("absent").unwrap_err().code(), "stale_source");
}

#[test]
fn persistence_is_deterministic_and_does_not_embed_source_paths() {
    let directory = Directory::new();
    let input = b">read description\naCgT\n";
    let first = directory.write("first-source", input);
    let second = directory.write("second-source", input);
    let first_index = build::<Dna>(&first, IndexedFormat::Fasta);
    let second_index = build::<Dna>(&second, IndexedFormat::Fasta);
    first_index.save_new(directory.path("first.json")).unwrap();
    second_index
        .save_new(directory.path("second.json"))
        .unwrap();
    assert_eq!(
        fs::read(directory.path("first.json")).unwrap(),
        fs::read(directory.path("second.json")).unwrap()
    );
    let mut relocated = SequenceIndex::<Dna>::open(
        &second,
        directory.path("first.json"),
        IndexOptions::default(),
    )
    .unwrap();
    assert_eq!(
        relocated.get("read").unwrap().unwrap().source_sequence(),
        "aCgT"
    );
}

#[test]
fn unicode_source_and_sidecar_paths_are_supported() {
    let directory = Directory::new();
    let source = directory.write("séquences-🧬.fasta", ">α description\naCgT\n".as_bytes());
    let sidecar = directory.path("índice-🧬.json");
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    index.save_new(&sidecar).unwrap();
    let mut reopened =
        SequenceIndex::<Dna>::open(&source, &sidecar, IndexOptions::default()).unwrap();
    assert_eq!(
        reopened.get("α").unwrap().unwrap().header(),
        "α description"
    );
}

#[cfg(unix)]
#[test]
fn non_utf8_paths_are_not_lossily_serialized_or_rejected() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;
    let directory = Directory::new();
    let source = directory.write(
        OsString::from_vec(b"source-\xff.fasta".to_vec()),
        b">a\nAC\n",
    );
    let sidecar = directory.path(OsString::from_vec(b"index-\xfe.json".to_vec()));
    assert!(source.to_str().is_none());
    assert!(sidecar.to_str().is_none());
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    index.save_new(&sidecar).unwrap();
    let mut reopened =
        SequenceIndex::<Dna>::open(&source, &sidecar, IndexOptions::default()).unwrap();
    assert_eq!(reopened.get("a").unwrap().unwrap().source_sequence(), "AC");
}

#[test]
fn compressed_sources_are_rejected_by_signature_while_plain_gzip_names_work() {
    let directory = Directory::new();
    let input = b">a\nAC\n";
    let source = directory.write("plain.fasta.gz", input);
    let sidecar = directory.path("plain.json");
    build::<Dna>(&source, IndexedFormat::Fasta)
        .save_new(&sidecar)
        .unwrap();
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(input).unwrap();
    let compressed = directory.write("compressed.fasta", &gzip.finish().unwrap());
    let bgzf = directory.write(
        "bgzf.fasta",
        &[
            31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0, 27, 0,
        ],
    );
    for path in [compressed, bgzf] {
        assert_eq!(
            SequenceIndex::<Dna>::build(&path, IndexedFormat::Fasta, IndexOptions::default())
                .err()
                .unwrap()
                .code(),
            "unsupported_index_source"
        );
        assert_eq!(
            SequenceIndex::<Dna>::open(&path, &sidecar, IndexOptions::default())
                .err()
                .unwrap()
                .code(),
            "unsupported_index_source"
        );
    }
}

#[test]
fn directories_and_missing_sources_report_explicit_errors() {
    let directory = Directory::new();
    assert_eq!(
        SequenceIndex::<Dna>::build(&directory.0, IndexedFormat::Fasta, IndexOptions::default())
            .err()
            .unwrap()
            .code(),
        "unsupported_index_source"
    );
    let missing = SequenceIndex::<Dna>::build(
        directory.path("missing"),
        IndexedFormat::Fasta,
        IndexOptions::default(),
    )
    .err()
    .unwrap();
    assert_eq!(missing.code(), "io_error");
    assert!(matches!(missing, IndexError::Io(_)));
    #[cfg(unix)]
    assert_eq!(
        SequenceIndex::<Dna>::build(
            Path::new("/dev/null"),
            IndexedFormat::Fasta,
            IndexOptions::default()
        )
        .err()
        .unwrap()
        .code(),
        "unsupported_index_source"
    );
}

#[test]
fn save_is_create_only_and_failure_never_changes_prior_bytes_or_leaves_temporary_files() {
    let directory = Directory::new();
    let source = directory.write("source", b">a\nAC\n");
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    let destination = directory.write("existing.json", b"user-owned destination\n");
    assert_eq!(
        index.save_new(&destination).unwrap_err().code(),
        "output_error"
    );
    assert_eq!(fs::read(&destination).unwrap(), b"user-owned destination\n");
    assert_eq!(index.save_new(&source).unwrap_err().code(), "output_error");
    assert_eq!(fs::read(&source).unwrap(), b">a\nAC\n");
    assert_eq!(
        index.save_new(&directory.0).unwrap_err().code(),
        "output_error"
    );
    assert!(
        index
            .save_new(directory.path("missing-parent/index.json"))
            .is_err()
    );
    assert!(
        index
            .save_new(directory.path("existing.json/index.json"))
            .is_err()
    );
    assert_eq!(directory.count(), 2);
    let new_path = directory.path("new.json");
    index.save_new(&new_path).unwrap();
    let bytes = fs::read(&new_path).unwrap();
    assert_eq!(
        index.save_new(&new_path).unwrap_err().code(),
        "output_error"
    );
    assert_eq!(fs::read(&new_path).unwrap(), bytes);
    assert_eq!(directory.count(), 3);
}

#[cfg(unix)]
#[test]
fn dangling_output_symlinks_are_never_followed_or_replaced() {
    use std::os::unix::fs::symlink;
    let directory = Directory::new();
    let source = directory.write("source", b">a\nAC\n");
    let index = build::<Dna>(&source, IndexedFormat::Fasta);
    let output = directory.path("index.json");
    symlink("absent-target", &output).unwrap();
    assert_eq!(index.save_new(&output).unwrap_err().code(), "output_error");
    assert!(
        fs::symlink_metadata(&output)
            .unwrap()
            .file_type()
            .is_symlink()
    );
    assert_eq!(fs::read_link(&output).unwrap(), Path::new("absent-target"));
    assert!(!directory.path("absent-target").exists());
    assert_eq!(directory.count(), 2);
}

#[cfg(unix)]
#[test]
fn index_fifo_probe() {
    let Some(source) = std::env::var_os("BIOFERRIS_INDEX_FIFO_SOURCE") else {
        return;
    };
    let sidecar = std::env::var_os("BIOFERRIS_INDEX_FIFO_SIDECAR").unwrap();
    assert_eq!(
        SequenceIndex::<Dna>::build(
            Path::new(&source),
            IndexedFormat::Fasta,
            IndexOptions::default()
        )
        .err()
        .unwrap()
        .code(),
        "unsupported_index_source"
    );
    assert_eq!(
        SequenceIndex::<Dna>::open(
            Path::new(&source),
            Path::new(&sidecar),
            IndexOptions::default()
        )
        .err()
        .unwrap()
        .code(),
        "unsupported_index_source"
    );
}

#[cfg(unix)]
#[test]
fn fifo_sources_are_rejected_before_an_open_can_block() {
    use std::process::{Command, Stdio};
    use std::time::{Duration, Instant};
    let directory = Directory::new();
    let source = directory.write("regular", b">a\nAC\n");
    let sidecar = directory.path("index.json");
    build::<Dna>(&source, IndexedFormat::Fasta)
        .save_new(&sidecar)
        .unwrap();
    let fifo = directory.path("source.fifo");
    // mkfifo is fixture setup only; all indexing behavior runs in this Rust
    // harness. A timed child makes a preflight regression fail instead of hang.
    assert!(
        Command::new("mkfifo")
            .arg(&fifo)
            .status()
            .unwrap()
            .success()
    );
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "index_fifo_probe", "--nocapture"])
        .env("BIOFERRIS_INDEX_FIFO_SOURCE", &fifo)
        .env("BIOFERRIS_INDEX_FIFO_SIDECAR", &sidecar)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(5);
    while child.try_wait().unwrap().is_none() {
        if Instant::now() >= deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("indexing a FIFO blocked instead of rejecting the source type");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stdout: {}; stderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
