//! Binary sidecars and quick source checks: layout, equivalence with JSON
//! sidecars, structural and checksum validation, budgets, and exactly what a
//! quick check does and does not detect.

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use bioferris::io::index::{
    IndexError, IndexOptions, IndexedFormat, SequenceIndex, SidecarFormat, SourceCheck,
};
use bioferris::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};

const EXAMPLE: &str = "tests/fixtures/index-example.fastq";

struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-binary-index-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.path(name);
        fs::write(&path, bytes).unwrap();
        path
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

fn example(directory: &Directory) -> PathBuf {
    directory.write("reads.fastq", &fs::read(EXAMPLE).unwrap())
}

fn build(source: &Path) -> SequenceIndex<Dna> {
    SequenceIndex::build(source, IndexedFormat::FastqSanger, IndexOptions::default()).unwrap()
}

fn save_binary(index: &SequenceIndex<Dna>, path: &Path) -> Vec<u8> {
    index.save_new_as(path, SidecarFormat::Binary).unwrap();
    fs::read(path).unwrap()
}

fn open(
    source: &Path,
    sidecar: &Path,
    check: SourceCheck,
) -> Result<SequenceIndex<Dna>, IndexError> {
    SequenceIndex::open_with(source, sidecar, IndexOptions::default(), check)
}

fn le64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

/// Recomputes the trailing checksum, so a test reaches the rule it targets.
fn reseal(mut bytes: Vec<u8>) -> Vec<u8> {
    let body = bytes.len() - 32;
    let checksum = Sha256::digest(&bytes[..body]);
    bytes[body..].copy_from_slice(&checksum);
    bytes
}

fn invalid_message(error: &IndexError) -> &str {
    match error {
        IndexError::InvalidIndex(message) => message,
        other => panic!("expected invalid_index, got {other:?}"),
    }
}

fn stale_message(error: &IndexError) -> &str {
    match error {
        IndexError::StaleSource(message) => message,
        other => panic!("expected stale_source, got {other:?}"),
    }
}

fn modified(path: &Path) -> SystemTime {
    fs::metadata(path).unwrap().modified().unwrap()
}

/// Rewrites a file with the same length and puts its modification time back.
fn rewrite_keeping_time(path: &Path, bytes: &[u8]) {
    let original = modified(path);
    assert_eq!(
        fs::metadata(path).unwrap().len(),
        u64::try_from(bytes.len()).unwrap()
    );
    fs::write(path, bytes).unwrap();
    File::options()
        .write(true)
        .open(path)
        .unwrap()
        .set_modified(original)
        .unwrap();
    assert_eq!(modified(path), original);
}

fn records(index: &mut SequenceIndex<Dna>) -> Vec<Record<Dna>> {
    let ids = index.ids().map(str::to_owned).collect::<Vec<_>>();
    ids.iter()
        .map(|id| index.get(id).unwrap().unwrap())
        .collect()
}

#[test]
fn binary_layout_matches_the_documented_bytes() {
    let directory = Directory::new();
    let source = example(&directory);
    let index = build(&source);
    index.save_new(directory.path("reads.json")).unwrap();
    let json: Value =
        serde_json::from_slice(&fs::read(directory.path("reads.json")).unwrap()).unwrap();
    let bytes = save_binary(&index, &directory.path("reads.bfidx"));

    // 88-byte header, 2 x 64-byte entries, "read1\nread2\n", 2 x 8-byte
    // lookup positions and a 32-byte checksum.
    assert_eq!(bytes.len(), 88 + 2 * 64 + 12 + 2 * 8 + 32);
    assert_eq!(&bytes[..8], b"\x89BFIDX\r\n");
    assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 2);
    assert_eq!(
        &bytes[12..16],
        [1, 2, 1, 0],
        "dna, fastq-sanger, key policy, reserved"
    );
    assert_eq!(le64(&bytes, 16), 2, "records");
    assert_eq!(le64(&bytes, 24), 10, "ID bytes without line feeds");
    assert_eq!(le64(&bytes, 32), 67, "source length");
    let time = modified(&source).duration_since(UNIX_EPOCH).unwrap();
    assert_eq!(le64(&bytes, 40), time.as_secs());
    assert_eq!(
        u32::from_le_bytes(bytes[48..52].try_into().unwrap()),
        time.subsec_nanos()
    );
    assert_eq!(&bytes[52..56], [0; 4]);
    assert_eq!(hex(&bytes[56..88]), json["source_sha256"].as_str().unwrap());
    assert_eq!(
        hex(&bytes[56..88]),
        hex(&Sha256::digest(fs::read(&source).unwrap()))
    );
    for (position, entry) in json["entries"].as_array().unwrap().iter().enumerate() {
        let at = 88 + 64 * position;
        assert_eq!(le64(&bytes, at), entry["offset"].as_u64().unwrap());
        assert_eq!(le64(&bytes, at + 8), entry["byte_length"].as_u64().unwrap());
        assert_eq!(
            le64(&bytes, at + 16),
            entry["sequence_length"].as_u64().unwrap()
        );
        assert_eq!(
            le64(&bytes, at + 24),
            entry["header_line"].as_u64().unwrap()
        );
        assert_eq!(
            hex(&bytes[at + 32..at + 64]),
            entry["sha256"].as_str().unwrap()
        );
    }
    // Literal spans: read1 is bytes 0..42 from line 1, read2 bytes 42..67 from line 5.
    assert_eq!(
        [
            le64(&bytes, 88),
            le64(&bytes, 96),
            le64(&bytes, 152),
            le64(&bytes, 160)
        ],
        [0, 42, 42, 25]
    );
    assert_eq!([le64(&bytes, 112), le64(&bytes, 176)], [1, 5]);
    assert_eq!(&bytes[216..228], b"read1\nread2\n");
    assert_eq!([le64(&bytes, 228), le64(&bytes, 236)], [0, 1]);
    assert_eq!(&bytes[244..], Sha256::digest(&bytes[..244]).as_slice());
}

fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;

    bytes.iter().fold(String::new(), |mut text, byte| {
        write!(text, "{byte:02x}").unwrap();
        text
    })
}

#[test]
fn binary_sidecars_reopen_to_the_records_of_json_sidecars() {
    let directory = Directory::new();
    let source = example(&directory);
    let index = build(&source);
    let json = directory.path("reads.json");
    let binary = directory.path("reads.bfidx");
    index.save_new(&json).unwrap();
    index.save_new_as(&binary, SidecarFormat::Binary).unwrap();
    drop(index);

    let mut from_json =
        SequenceIndex::<Dna>::open(&source, &json, IndexOptions::default()).unwrap();
    let expected = records(&mut from_json);
    assert_eq!(expected.len(), 2);
    assert_eq!(expected[1].source_sequence(), "Nn");
    assert_eq!(expected[1].qualities(), Some([31, 10].as_slice()));
    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut from_binary = open(&source, &binary, check).unwrap();
        assert_eq!(from_binary.ids().collect::<Vec<_>>(), ["read1", "read2"]);
        assert_eq!(
            from_binary.ids().rev().collect::<Vec<_>>(),
            ["read2", "read1"]
        );
        assert_eq!(from_binary.len(), 2);
        assert_eq!(from_binary.source_len(), 67);
        assert_eq!(from_binary.format(), IndexedFormat::FastqSanger);
        assert_eq!(from_binary.source_sha256(), from_json.source_sha256());
        assert_eq!(records(&mut from_binary), expected);
        assert!(from_binary.get("absent").unwrap().is_none());
        assert!(from_binary.get("read1 coding").unwrap().is_none());
    }
}

#[test]
fn either_sidecar_saves_the_same_bytes_in_either_format() {
    let directory = Directory::new();
    let source = example(&directory);
    let index = build(&source);
    index.save_new(directory.path("built.json")).unwrap();
    let built_binary = save_binary(&index, &directory.path("built.bfidx"));
    let built_json = fs::read(directory.path("built.json")).unwrap();

    let from_binary = open(&source, &directory.path("built.bfidx"), SourceCheck::Full).unwrap();
    from_binary.save_new(directory.path("again.json")).unwrap();
    assert_eq!(fs::read(directory.path("again.json")).unwrap(), built_json);
    assert_eq!(
        save_binary(&from_binary, &directory.path("again.bfidx")),
        built_binary
    );

    let from_json = open(&source, &directory.path("built.json"), SourceCheck::Full).unwrap();
    assert_eq!(
        save_binary(&from_json, &directory.path("converted.bfidx")),
        built_binary
    );
    // Saving stays create-only for the new format.
    let error = from_json
        .save_new_as(directory.path("built.bfidx"), SidecarFormat::Binary)
        .unwrap_err();
    assert_eq!(error.code(), "output_error");
    assert_eq!(
        fs::read(directory.path("built.bfidx")).unwrap(),
        built_binary
    );
}

/// IDs that need JSON escaping, and the sidecar `096a948b` wrote for them,
/// before JSON was written from the binary-ready representation.
const ESCAPES: &[u8] = b"@q\"uote x\nACGT\n+\nIIII\n@back\\slash\nAC\n+\nII\n@\xc2\xb5\nA\n+\nI\n";
const ESCAPES_JSON: &str = concat!(
    r#"{"format_id":"bioferris-sequence-index","schema_version":1,"molecule":"dna","#,
    r#""source_format":"fastq-sanger","key_policy":"record-id-first-token-v1","source_len":52,"#,
    r#""source_sha256":"89f236b6b731c169a18a1d7ac527aa5ba6c65b5cf68196db33b031a924441747","#,
    r#""entries":[{"id":"q\"uote","offset":0,"byte_length":22,"sequence_length":4,"record":1,"#,
    r#""header_line":1,"sha256":"086ead9243dfa6d729628559f7dd812c951a4db9a30464b0bfb6deb15e02835b"},"#,
    r#"{"id":"back\\slash","offset":22,"byte_length":20,"sequence_length":2,"record":2,"#,
    r#""header_line":5,"sha256":"b9f09565b3fb8c4297678250be6402ca3b5089bd40ada469d474429e12137c02"},"#,
    r#"{"id":"µ","offset":42,"byte_length":10,"sequence_length":1,"record":3,"header_line":9,"#,
    r#""sha256":"d512e4c83e805859331dab92ce699e1093d3d9e5d60319c7f107287b688fb0b3"}]}"#,
);

#[test]
fn json_sidecars_keep_the_bytes_of_earlier_versions() {
    let directory = Directory::new();
    let source = directory.write("escapes.fastq", ESCAPES);
    let index = build(&source);
    index.save_new(directory.path("built.json")).unwrap();
    assert_eq!(
        fs::read_to_string(directory.path("built.json")).unwrap(),
        ESCAPES_JSON
    );
    save_binary(&index, &directory.path("escapes.bfidx"));
    for (name, sidecar) in [
        ("from-json.json", "built.json"),
        ("from-binary.json", "escapes.bfidx"),
    ] {
        let reopened = open(&source, &directory.path(sidecar), SourceCheck::Full).unwrap();
        assert_eq!(
            reopened.ids().collect::<Vec<_>>(),
            ["q\"uote", "back\\slash", "µ"]
        );
        reopened.save_new(directory.path(name)).unwrap();
        assert_eq!(
            fs::read_to_string(directory.path(name)).unwrap(),
            ESCAPES_JSON,
            "{sidecar}"
        );
    }
}

#[test]
fn empty_sources_have_header_only_binary_sidecars() {
    let directory = Directory::new();
    let source = directory.write("empty.fastq", b"");
    let bytes = save_binary(&build(&source), &directory.path("empty.bfidx"));
    assert_eq!(bytes.len(), 120);
    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut index = open(&source, &directory.path("empty.bfidx"), check).unwrap();
        assert!(index.is_empty());
        assert_eq!(index.ids().len(), 0);
        assert!(index.get("read1").unwrap().is_none());
    }
}

#[test]
fn quick_checks_need_a_binary_sidecar() {
    let directory = Directory::new();
    let source = example(&directory);
    build(&source)
        .save_new(directory.path("reads.json"))
        .unwrap();
    let error = open(&source, &directory.path("reads.json"), SourceCheck::Quick)
        .err()
        .unwrap();
    assert_eq!(
        invalid_message(&error),
        "quick source checks need a binary sidecar, which records the source modification time"
    );
    open(&source, &directory.path("reads.json"), SourceCheck::Full).unwrap();
}

/// One hundred four-base reads `r000` to `r099`; a quick check samples
/// positions 0, 3, 6, 9, ... 99, so `r001` is not sampled and `r003` is.
fn hundred_reads(sequence_of: impl Fn(usize) -> &'static str) -> Vec<u8> {
    (0..100)
        .flat_map(|read| format!("@r{read:03}\n{}\n+\nIIII\n", sequence_of(read)).into_bytes())
        .collect()
}

#[test]
fn quick_checks_detect_length_time_and_sampled_record_changes() {
    let directory = Directory::new();
    let original = hundred_reads(|_| "ACGT");
    let source = directory.write("reads.fastq", &original);
    let sidecar = directory.path("reads.bfidx");
    save_binary(&build(&source), &sidecar);
    open(&source, &sidecar, SourceCheck::Quick).unwrap();

    // A different length.
    let mut longer = original.clone();
    longer.extend_from_slice(b"\n");
    fs::write(&source, &longer).unwrap();
    let error = open(&source, &sidecar, SourceCheck::Quick).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "source byte length differs from the index"
    );

    // The same bytes with a new modification time: quick refuses, full accepts.
    fs::write(&source, &original).unwrap();
    File::options()
        .write(true)
        .open(&source)
        .unwrap()
        .set_modified(UNIX_EPOCH + std::time::Duration::from_secs(1_000_000_000))
        .unwrap();
    let error = open(&source, &sidecar, SourceCheck::Quick).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "source modification time differs from the index"
    );
    open(&source, &sidecar, SourceCheck::Full).unwrap();

    // A same-length change to a sampled record, with the time restored.
    let source = directory.write("sampled.fastq", &original);
    let sidecar = directory.path("sampled.bfidx");
    save_binary(&build(&source), &sidecar);
    rewrite_keeping_time(
        &source,
        &hundred_reads(|read| if read == 3 { "ACGA" } else { "ACGT" }),
    );
    let error = open(&source, &sidecar, SourceCheck::Quick).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "record span bytes differ from the verified SHA-256"
    );
}

#[test]
fn quick_checks_miss_unsampled_changes_that_keep_length_and_time() {
    let directory = Directory::new();
    let source = directory.write("reads.fastq", &hundred_reads(|_| "ACGT"));
    let sidecar = directory.path("reads.bfidx");
    save_binary(&build(&source), &sidecar);
    rewrite_keeping_time(
        &source,
        &hundred_reads(|read| if read == 1 { "ACGA" } else { "ACGT" }),
    );

    // The documented limit: the quick check passes, and the index describes the
    // indexed bytes rather than the current file.
    let mut index = open(&source, &sidecar, SourceCheck::Quick).unwrap();
    assert_eq!(index.len(), 100);
    // A lookup still never returns the changed record...
    let error = index.get("r001").unwrap_err();
    assert_eq!(
        stale_message(&error),
        "record span bytes differ from the verified SHA-256"
    );
    // ...and unchanged records remain retrievable.
    assert_eq!(
        index.get("r002").unwrap().unwrap().source_sequence(),
        "ACGT"
    );
    let error = open(&source, &sidecar, SourceCheck::Full).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "source SHA-256 differs from the index"
    );
}

#[test]
fn damaged_binary_sidecars_fail_their_checksum_or_length() {
    let directory = Directory::new();
    let source = example(&directory);
    let valid = save_binary(&build(&source), &directory.path("valid.bfidx"));
    let mut flipped = valid.clone();
    flipped[120] ^= 1;
    let mut truncated = valid.clone();
    truncated.pop();
    let mut extended = valid.clone();
    extended.push(0);
    let cases = [
        (
            flipped,
            "binary sidecar checksum does not match its contents",
        ),
        (
            truncated,
            "binary sidecar length does not match its record and ID counts",
        ),
        (
            extended,
            "binary sidecar length does not match its record and ID counts",
        ),
        (
            valid[..100].to_vec(),
            "binary sidecar is shorter than its header and checksum",
        ),
    ];
    for (ordinal, (bytes, message)) in cases.into_iter().enumerate() {
        let path = directory.write(&format!("damaged-{ordinal}.bfidx"), &bytes);
        for check in [SourceCheck::Full, SourceCheck::Quick] {
            let error = open(&source, &path, check).err().unwrap();
            assert_eq!(invalid_message(&error), message, "case {ordinal}");
        }
    }
}

#[test]
fn structural_rules_apply_even_with_a_valid_checksum() {
    let directory = Directory::new();
    let source = example(&directory);
    let valid = save_binary(&build(&source), &directory.path("valid.bfidx"));
    let forged = |offset: usize, replacement: &[u8]| {
        let mut bytes = valid.clone();
        bytes[offset..offset + replacement.len()].copy_from_slice(replacement);
        reseal(bytes)
    };
    let cases = [
        (
            forged(8, &3u32.to_le_bytes()),
            "unsupported format ID, schema version or key policy",
        ),
        (
            forged(14, &[2]),
            "unsupported format ID, schema version or key policy",
        ),
        (
            forged(12, &[2]),
            "index molecule \"rna\" does not match dna",
        ),
        (forged(12, &[9]), "unknown molecule code in binary sidecar"),
        (
            forged(13, &[9]),
            "unknown source format code in binary sidecar",
        ),
        (
            forged(15, &[1]),
            "reserved header bytes or the modification time are invalid",
        ),
        (
            forged(48, &1_000_000_000u32.to_le_bytes()),
            "reserved header bytes or the modification time are invalid",
        ),
        (
            forged(160, &0u64.to_le_bytes()),
            "entry coordinates, order or record metadata are inconsistent",
        ),
        (
            forged(152, &0u64.to_le_bytes()),
            "entry coordinates, order or record metadata are inconsistent",
        ),
        (
            forged(88, &u64::MAX.to_le_bytes()),
            "record span overflows its coordinate type",
        ),
        (
            forged(216, b"rea 1"),
            "entry ID must be a nonempty first-token identifier",
        ),
        (forged(216, &[0xff]), "record IDs must be UTF-8"),
        (
            forged(216, b"rea\nd1\nead2\n"),
            "the ID section must hold one line-feed-terminated ID per record",
        ),
        (
            forged(222, b"read1"),
            "sidecar contains duplicate record IDs",
        ),
        (
            forged(228, &[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            "lookup positions are not in ascending ID order",
        ),
        (
            forged(236, &2u64.to_le_bytes()),
            "lookup position is outside the entry table",
        ),
    ];
    for (ordinal, (bytes, message)) in cases.into_iter().enumerate() {
        let path = directory.write(&format!("forged-{ordinal}.bfidx"), &bytes);
        for check in [SourceCheck::Full, SourceCheck::Quick] {
            let error = open(&source, &path, check).err().unwrap();
            assert_eq!(invalid_message(&error), message, "case {ordinal}");
        }
    }
}

#[test]
fn full_checks_reject_consistent_forgeries_that_quick_checks_cannot_all_see() {
    let directory = Directory::new();
    let source = example(&directory);
    let valid = save_binary(&build(&source), &directory.path("valid.bfidx"));
    let forge = |name: &str, offset: usize, replacement: &[u8]| {
        let mut bytes = valid.clone();
        bytes[offset..offset + replacement.len()].copy_from_slice(replacement);
        directory.write(name, &reseal(bytes))
    };
    let differ = "saved entries differ from the fully rederived source records";

    // A saved span digest: the source still hashes as saved, the span does not.
    let digest = forge("digest.bfidx", 120, &[0; 32]);
    assert_eq!(
        invalid_message(&open(&source, &digest, SourceCheck::Full).err().unwrap()),
        differ
    );
    let error = open(&source, &digest, SourceCheck::Quick).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "record span bytes differ from the verified SHA-256"
    );

    // An identifier that keeps the lookup order.
    let id = forge("id.bfidx", 222, b"read3");
    assert_eq!(
        invalid_message(&open(&source, &id, SourceCheck::Full).err().unwrap()),
        differ
    );
    let error = open(&source, &id, SourceCheck::Quick).err().unwrap();
    assert_eq!(
        invalid_message(&error),
        "retrieved record ID or sequence length differs from the index"
    );

    // Header line numbers and the whole-source digest are only verified in full.
    let line = forge("line.bfidx", 176, &6u64.to_le_bytes());
    assert_eq!(
        invalid_message(&open(&source, &line, SourceCheck::Full).err().unwrap()),
        differ
    );
    open(&source, &line, SourceCheck::Quick).unwrap();
    let whole = forge("whole.bfidx", 56, &[7; 32]);
    let error = open(&source, &whole, SourceCheck::Full).err().unwrap();
    assert_eq!(
        stale_message(&error),
        "source SHA-256 differs from the index"
    );
    let quick = open(&source, &whole, SourceCheck::Quick).unwrap();
    assert_eq!(quick.source_sha256(), "07".repeat(32));
}

#[test]
fn a_different_source_is_refused_by_both_checks() {
    let directory = Directory::new();
    let source = example(&directory);
    let sidecar = directory.path("reads.bfidx");
    save_binary(&build(&source), &sidecar);
    let mut other = fs::read(&source).unwrap();
    other[15] = b'c'; // atgGCT -> acgGCT in read1, the same length.
    let other_path = directory.write("other.fastq", &other);
    File::options()
        .write(true)
        .open(&other_path)
        .unwrap()
        .set_modified(modified(&source))
        .unwrap();
    let error = open(&other_path, &sidecar, SourceCheck::Full)
        .err()
        .unwrap();
    assert_eq!(
        stale_message(&error),
        "source SHA-256 differs from the index"
    );
    let error = open(&other_path, &sidecar, SourceCheck::Quick)
        .err()
        .unwrap();
    assert_eq!(
        stale_message(&error),
        "record span bytes differ from the verified SHA-256"
    );
    let error = SequenceIndex::<Rna>::open(&source, &sidecar, IndexOptions::default())
        .err()
        .unwrap();
    assert_eq!(
        invalid_message(&error),
        "index molecule \"dna\" does not match rna"
    );
}

#[test]
fn binary_budgets_bound_records_ids_and_sidecar_bytes() {
    let directory = Directory::new();
    let source = example(&directory);
    let sidecar = directory.path("reads.bfidx");
    let bytes = save_binary(&build(&source), &sidecar);
    let options = |records: usize, id_bytes: usize, index_bytes: usize| IndexOptions {
        max_records: records,
        max_id_bytes: id_bytes,
        max_index_bytes: index_bytes,
        ..IndexOptions::default()
    };
    let size = u64::try_from(bytes.len()).unwrap();
    let exact = options(2, 10, bytes.len());
    for (limits, resource, limit, observed) in [
        (options(1, 10, bytes.len()), "index_records", 1, 2),
        (options(2, 9, bytes.len()), "index_id_bytes", 9, 10),
        (
            options(2, 10, bytes.len() - 1),
            "index_bytes",
            size - 1,
            size,
        ),
    ] {
        for check in [SourceCheck::Full, SourceCheck::Quick] {
            let error = SequenceIndex::<Dna>::open_with(&source, &sidecar, limits, check)
                .err()
                .unwrap();
            assert!(
                matches!(error, IndexError::LimitExceeded { resource: found, limit: l, observed: o }
                    if found == resource && l == limit && o == observed),
                "{error:?}"
            );
        }
    }
    let index =
        SequenceIndex::<Dna>::open_with(&source, &sidecar, exact, SourceCheck::Full).unwrap();
    // The JSON form is larger than this index's sidecar budget allows.
    let before = directory.count();
    let error = index.save_new(directory.path("reads.json")).unwrap_err();
    assert!(
        matches!(error, IndexError::LimitExceeded { resource: "index_bytes", limit, .. }
            if limit == size),
        "{error:?}"
    );
    assert_eq!(
        directory.count(),
        before,
        "no sidecar or temporary file is left"
    );
    assert_eq!(save_binary(&index, &directory.path("copy.bfidx")), bytes);
}

#[test]
fn sidecar_changes_after_opening_are_detected() {
    let directory = Directory::new();
    let source = example(&directory);
    let sidecar = directory.path("reads.bfidx");
    save_binary(&build(&source), &sidecar);
    let mut index = open(&source, &sidecar, SourceCheck::Full).unwrap();
    assert!(index.get("read1").unwrap().is_some());
    File::options()
        .append(true)
        .open(&sidecar)
        .unwrap()
        .write_all(b"x")
        .unwrap();
    for id in ["read1", "absent"] {
        let error = index.get(id).unwrap_err();
        assert_eq!(
            invalid_message(&error),
            "sidecar identity, size or timestamps differ from the verified file"
        );
    }
}
