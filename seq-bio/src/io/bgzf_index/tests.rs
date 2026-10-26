use super::*;
use crate::io::bgzf::{BGZF_EOF_MARKER, BgzfWriteOptions, BgzfWriter, write_bgzf_new};
use crate::io::gzip::{GzipWriteOptions, GzipWriter, write_gzip_new};
use crate::io::{LossPolicy, write_record_fasta, write_record_fastq};
use crate::{Dna, Protein, Rna};
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-bgzf-index-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

fn encode(input: &[u8], size: usize) -> Vec<u8> {
    let mut writer = BgzfWriter::with_options(
        Vec::new(),
        BgzfWriteOptions {
            block_size: size,
            ..Default::default()
        },
    )
    .unwrap();
    writer.write_all(input).unwrap();
    writer.finish().unwrap()
}

fn build(
    directory: &Directory,
    input: &[u8],
    format: IndexedFormat,
    size: usize,
) -> BgzfSequenceIndex<Dna> {
    let source = directory.path("source.bgz");
    std::fs::write(&source, encode(input, size)).unwrap();
    BgzfSequenceIndex::build(source, format, BgzfIndexOptions::default()).unwrap()
}

#[test]
fn fasta_headers_case_empty_records_crlf_wrap_and_source_order_survive() {
    let directory = Directory::new();
    let input = b"\r\n>z description\r\naCg\r\nTN\r\n>a empty\r\n>m\r\nTT\r\n";
    let mut index = build(&directory, input, IndexedFormat::Fasta, 7);
    assert_eq!(index.ids().collect::<Vec<_>>(), ["z", "a", "m"]);
    assert_eq!(
        index
            .get("z")
            .unwrap()
            .unwrap()
            .source_sequence()
            .as_bytes(),
        b"aCgTN"
    );
    assert_eq!(index.get("z").unwrap().unwrap().header(), "z description");
    assert!(index.get("a").unwrap().unwrap().sequence().is_empty());
    assert!(index.get("Z").unwrap().is_none());
    assert_eq!(index.uncompressed_len(), input.len() as u64);
    assert_eq!(index.format(), IndexedFormat::Fasta);
    assert_eq!(index.len(), 3);
    assert!(!index.is_empty());
    let sidecar = directory.path("index.json");
    index.save_new(&sidecar).unwrap();
    let saved = std::fs::read(&sidecar).unwrap();
    assert!(
        !String::from_utf8(saved.clone())
            .unwrap()
            .contains("source.bgz")
    );
    let mut opened = BgzfSequenceIndex::<Dna>::open(
        directory.path("source.bgz"),
        &sidecar,
        BgzfIndexOptions::default(),
    )
    .unwrap();
    let second = directory.path("second.json");
    opened.save_new(&second).unwrap();
    assert_eq!(saved, std::fs::read(second).unwrap());
    for key in ["m", "a", "z", "z", "absent"] {
        assert_eq!(
            opened
                .get(key)
                .unwrap()
                .map(|record| record.header().to_owned()),
            index
                .get(key)
                .unwrap()
                .map(|record| record.header().to_owned())
        );
    }
}

#[test]
fn fastq_wrapping_header_and_exact_sanger_qualities_survive_block_crossings() {
    let directory = Directory::new();
    let input =
        b"@r mixed case\r\naCg\r\nTN\r\n+r mixed case\r\n!I~\r\n@+\r\n@q\r\nAA\r\n+\r\n~~\r\n";
    let mut index = build(&directory, input, IndexedFormat::FastqSanger, 3);
    let record = index.get("r").unwrap().unwrap();
    assert_eq!(record.source_sequence().as_bytes(), b"aCgTN");
    assert_eq!(record.header(), "r mixed case");
    assert_eq!(record.qualities().unwrap(), &[0, 40, 93, 31, 10]);
    assert_eq!(
        index.get("q").unwrap().unwrap().qualities().unwrap(),
        &[93, 93]
    );
    let sidecar = directory.path("index.json");
    index.save_new(&sidecar).unwrap();
    let mut opened = BgzfSequenceIndex::<Dna>::open(
        directory.path("source.bgz"),
        sidecar,
        BgzfIndexOptions::default(),
    )
    .unwrap();
    assert_eq!(
        opened.get("r").unwrap().unwrap().qualities(),
        record.qualities()
    );
}

// A directly constructed RFC 1951 final stored block, independent of the writer.
fn stored(input: &[u8]) -> Vec<u8> {
    let len = u16::try_from(input.len()).unwrap();
    let size = 18 + 5 + input.len() + 8;
    let mut bytes = vec![31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0];
    bytes.extend_from_slice(&u16::try_from(size - 1).unwrap().to_le_bytes());
    bytes.push(1);
    bytes.extend_from_slice(&len.to_le_bytes());
    bytes.extend_from_slice(&(!len).to_le_bytes());
    bytes.extend_from_slice(input);
    bytes.extend_from_slice(&crc32fast::hash(input).to_le_bytes());
    bytes.extend_from_slice(&u32::from(len).to_le_bytes());
    bytes
}

#[test]
fn embedded_empty_members_map_canonically_and_all_equivalent_ends_retrieve() {
    let directory = Directory::new();
    let mut raw = BGZF_EOF_MARKER.to_vec();
    raw.extend_from_slice(&stored(b">a\nAC"));
    let boundary = raw.len() as u64;
    raw.extend_from_slice(&BGZF_EOF_MARKER);
    raw.extend_from_slice(&stored(b"\n>b\nTT\n"));
    let last_empty = raw.len() as u64;
    raw.extend_from_slice(&stored(b""));
    raw.extend_from_slice(&BGZF_EOF_MARKER);
    std::fs::write(directory.path("source.bgz"), &raw).unwrap();
    let mut index = BgzfSequenceIndex::<Dna>::build(
        directory.path("source.bgz"),
        IndexedFormat::Fasta,
        BgzfIndexOptions::default(),
    )
    .unwrap();
    assert_eq!(index.block_count(), 6);
    assert_eq!(index.virtual_offset(0).unwrap().raw(), 0);
    assert_eq!(
        index.virtual_offset(5).unwrap(),
        BgzfVirtualOffset::new(boundary, 0).unwrap()
    );
    assert_eq!(
        index.virtual_offset(index.uncompressed_len()).unwrap(),
        BgzfVirtualOffset::new(last_empty, 0).unwrap()
    );
    for block in &index.manifest.blocks {
        assert_eq!(
            index
                .uncompressed_offset(BgzfVirtualOffset::new(block.compressed_offset, 0).unwrap())
                .unwrap(),
            block.uncompressed_offset
        );
        assert_eq!(
            index
                .uncompressed_offset(
                    BgzfVirtualOffset::new(block.compressed_offset, block.uncompressed_size)
                        .unwrap()
                )
                .unwrap(),
            block.uncompressed_offset + u64::from(block.uncompressed_size)
        );
    }
    assert_eq!(
        index
            .uncompressed_offset(BgzfVirtualOffset::new(raw.len() as u64, 0).unwrap())
            .unwrap(),
        index.uncompressed_len()
    );
    assert_eq!(
        index
            .get("a")
            .unwrap()
            .unwrap()
            .source_sequence()
            .as_bytes(),
        b"AC"
    );
    assert_eq!(
        index
            .get("b")
            .unwrap()
            .unwrap()
            .source_sequence()
            .as_bytes(),
        b"TT"
    );
    index.save_new(directory.path("index.json")).unwrap();
    BgzfSequenceIndex::<Dna>::open(
        directory.path("source.bgz"),
        directory.path("index.json"),
        BgzfIndexOptions::default(),
    )
    .unwrap();
}

#[test]
fn full_65536_byte_block_end_uses_next_block_not_wrapped_data_offset() {
    let directory = Directory::new();
    let mut input = b">a\n".to_vec();
    input.resize(65_535, b'A');
    input.extend_from_slice(b"\n>b\nTT\n");
    let mut index = build(&directory, &input, IndexedFormat::Fasta, 65_536);
    assert_eq!(index.manifest.blocks[0].uncompressed_size, 65_536);
    let boundary = index.virtual_offset(65_536).unwrap();
    assert_eq!(boundary.uncompressed_offset(), 0);
    assert_eq!(
        boundary.compressed_offset(),
        u64::from(index.manifest.blocks[0].compressed_size)
    );
    assert_eq!(index.manifest.entries[0].virtual_end, boundary.raw());
    assert_eq!(index.manifest.entries[1].virtual_start, boundary.raw());
    assert_eq!(index.get("a").unwrap().unwrap().sequence().len(), 65_532);
    assert_eq!(
        index
            .get("b")
            .unwrap()
            .unwrap()
            .source_sequence()
            .as_bytes(),
        b"TT"
    );
}

#[test]
fn every_logical_position_round_trips_under_many_block_layouts() {
    for size in [1, 2, 3, 4, 7, 17, 64, 256, 65_536] {
        let directory = Directory::new();
        let input = b"\n>r description\r\naCgT\nNN\n>x\nTTT\n>empty\n";
        let mut index = build(&directory, input, IndexedFormat::Fasta, size);
        for position in 0..=input.len() as u64 {
            let virtual_offset = index.virtual_offset(position).unwrap();
            assert_eq!(
                index.uncompressed_offset(virtual_offset).unwrap(),
                position,
                "block size {size}, position {position}"
            );
        }
        assert_eq!(
            index
                .get("r")
                .unwrap()
                .unwrap()
                .source_sequence()
                .as_bytes(),
            b"aCgTNN"
        );
        assert!(index.virtual_offset(input.len() as u64 + 1).is_err());
        assert!(
            index
                .uncompressed_offset(BgzfVirtualOffset::new(1, 0).unwrap())
                .is_err()
        );
        assert!(
            index
                .uncompressed_offset(BgzfVirtualOffset::new(0, u16::MAX.into()).unwrap())
                .is_err()
        );
    }
}

#[test]
fn empty_document_has_one_required_member_and_all_zero_record_budgets_work() {
    let directory = Directory::new();
    std::fs::write(directory.path("source.bgz"), BGZF_EOF_MARKER).unwrap();
    let mut options = BgzfIndexOptions::default();
    options.index.max_records = 0;
    options.index.max_id_bytes = 0;
    options.index.read_limits.max_total_bytes = 0;
    options.max_blocks = 1;
    options.max_compressed_bytes = 28;
    let mut index = BgzfSequenceIndex::<Dna>::build(
        directory.path("source.bgz"),
        IndexedFormat::Fasta,
        options,
    )
    .unwrap();
    assert!(index.is_empty());
    assert!(index.get("absent").unwrap().is_none());
    assert_eq!(index.block_count(), 1);
    assert_eq!(index.virtual_offset(0).unwrap().raw(), 0);
}

#[test]
fn exact_inclusive_and_one_below_limits_are_enforced() {
    let directory = Directory::new();
    let input = b">alpha\nAC\n>beta\nTT\n";
    let index = build(&directory, input, IndexedFormat::Fasta, 8);
    let mut exact = BgzfIndexOptions {
        max_blocks: index.block_count(),
        max_compressed_bytes: index.source_len(),
        ..Default::default()
    };
    exact.index.max_records = 2;
    exact.index.max_id_bytes = 9;
    exact.index.max_index_bytes =
        usize::try_from(serialized_size(&index.manifest).unwrap()).unwrap();
    exact.index.read_limits.max_total_bytes = input.len();
    BgzfSequenceIndex::<Dna>::build(directory.path("source.bgz"), IndexedFormat::Fasta, exact)
        .unwrap();
    for field in 0..6 {
        let mut options = exact;
        match field {
            0 => options.max_blocks -= 1,
            1 => options.max_compressed_bytes -= 1,
            2 => options.index.max_records -= 1,
            3 => options.index.max_id_bytes -= 1,
            4 => options.index.max_index_bytes -= 1,
            _ => options.index.read_limits.max_total_bytes -= 1,
        }
        let error = BgzfSequenceIndex::<Dna>::build(
            directory.path("source.bgz"),
            IndexedFormat::Fasta,
            options,
        )
        .unwrap_err();
        assert_eq!(error.code(), "limit_exceeded", "field {field}: {error}");
    }
    for field in 0..6 {
        let mut options = exact;
        match field {
            0 => options.max_blocks = 0,
            1 => options.max_compressed_bytes = 0,
            2 => options.index.max_records = 0,
            3 => options.index.max_id_bytes = 0,
            4 => options.index.max_index_bytes = 0,
            _ => options.index.read_limits.max_total_bytes = 0,
        }
        assert_eq!(
            BgzfSequenceIndex::<Dna>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                options
            )
            .unwrap_err()
            .code(),
            "limit_exceeded"
        );
    }
}

#[test]
fn duplicate_ids_and_ordinary_gzip_and_missing_marker_fail() {
    let directory = Directory::new();
    for input in [
        b">r one\nA\n>r two\nT\n".as_slice(),
        b"@r one\nA\n+\n!\n@r two\nT\n+\nI\n",
    ] {
        std::fs::write(directory.path("source.bgz"), encode(input, 4)).unwrap();
        let format = if input[0] == b'>' {
            IndexedFormat::Fasta
        } else {
            IndexedFormat::FastqSanger
        };
        let error = BgzfSequenceIndex::<Dna>::build(
            directory.path("source.bgz"),
            format,
            BgzfIndexOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            BgzfIndexError::Index(IndexError::DuplicateId {
                first_record: 1,
                duplicate_record: 2,
                ..
            })
        ));
    }
    let mut gzip = GzipWriter::new(Vec::new());
    gzip.write_all(b">r\nAC\n").unwrap();
    std::fs::write(directory.path("source.bgz"), gzip.finish().unwrap()).unwrap();
    assert_eq!(
        BgzfSequenceIndex::<Dna>::build(
            directory.path("source.bgz"),
            IndexedFormat::Fasta,
            BgzfIndexOptions::default()
        )
        .unwrap_err()
        .code(),
        "invalid_bgzf_header"
    );
    std::fs::write(directory.path("source.bgz"), stored(b">r\nAC\n")).unwrap();
    assert_eq!(
        BgzfSequenceIndex::<Dna>::build(
            directory.path("source.bgz"),
            IndexedFormat::Fasta,
            BgzfIndexOptions::default()
        )
        .unwrap_err()
        .code(),
        "missing_bgzf_eof"
    );
}

#[test]
fn malformed_biological_and_compressed_inputs_never_create_an_index() {
    let directory = Directory::new();
    for input in [b"AC\n".as_slice(), b">\nAC\n", b">r\nAX\n"] {
        std::fs::write(directory.path("source.bgz"), encode(input, 4)).unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
    for input in [
        b"@r\nAC\n+\n!\n".as_slice(),
        b"@r\nAC\n+q\n!!\n",
        b"@r\nAC\n+\n\x7f!\n",
    ] {
        std::fs::write(directory.path("source.bgz"), encode(input, 4)).unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::build(
                directory.path("source.bgz"),
                IndexedFormat::FastqSanger,
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
    let valid = encode(b">r\nAC\n", 4);
    for length in 0..valid.len() {
        std::fs::write(directory.path("source.bgz"), &valid[..length]).unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
}

#[test]
fn forged_sidecars_are_rejected_even_with_correct_source_sha256() {
    let directory = Directory::new();
    let index = build(&directory, b">a\nAC\n>b\nTT\n", IndexedFormat::Fasta, 5);
    let original = serde_json::to_value(&index.manifest).unwrap();
    for field in [
        "offset",
        "byte_length",
        "sequence_length",
        "virtual_start",
        "virtual_end",
        "record",
        "header_line",
    ] {
        let mut value = original.clone();
        let number = value["entries"][0][field].as_u64().unwrap();
        value["entries"][0][field] = (number + 1).into();
        std::fs::write(
            directory.path("index.json"),
            serde_json::to_vec(&value).unwrap(),
        )
        .unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::open(
                directory.path("source.bgz"),
                directory.path("index.json"),
                BgzfIndexOptions::default()
            )
            .is_err(),
            "field {field}"
        );
    }
    for field in [
        "compressed_offset",
        "compressed_size",
        "uncompressed_offset",
        "uncompressed_size",
    ] {
        let mut value = original.clone();
        value["blocks"][0][field] = (value["blocks"][0][field].as_u64().unwrap() + 1).into();
        std::fs::write(
            directory.path("index.json"),
            serde_json::to_vec(&value).unwrap(),
        )
        .unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::open(
                directory.path("source.bgz"),
                directory.path("index.json"),
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
    for mutation in 0..5 {
        let mut value = original.clone();
        match mutation {
            0 => value["entries"][0]["sha256"] = "a".repeat(64).into(),
            1 => value["entries"][0]["id"] = "different".into(),
            2 => value["blocks"][0]["canonical_eof_marker"] = true.into(),
            3 => value["uncompressed_len"] = (index.uncompressed_len() + 1).into(),
            _ => value["molecule"] = "rna".into(),
        }
        std::fs::write(
            directory.path("index.json"),
            serde_json::to_vec(&value).unwrap(),
        )
        .unwrap();
        assert!(
            BgzfSequenceIndex::<Dna>::open(
                directory.path("source.bgz"),
                directory.path("index.json"),
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
}

#[test]
fn strict_object_schema_refuses_duplicates_unknown_fields_and_positional_arrays() {
    let directory = Directory::new();
    let index = build(&directory, b">a\nAC\n", IndexedFormat::Fasta, 5);
    let original = serde_json::to_string(&index.manifest).unwrap();
    let positional = serde_json::to_value(&index.manifest)
        .unwrap()
        .as_object()
        .unwrap()
        .values()
        .cloned()
        .collect::<Vec<_>>();
    for content in [
        original.replacen('{', "{\"schema_version\":1,", 1),
        original.replacen('{', "{\"unexpected\":1,", 1),
        original.replace("\"id\":\"a\"", "\"id\":\"a\",\"id\":\"a\""),
        original.replace(
            "\"compressed_offset\":0",
            "\"compressed_offset\":0,\"unexpected\":1",
        ),
        serde_json::to_string(&positional).unwrap(),
    ] {
        std::fs::write(directory.path("index.json"), content).unwrap();
        assert_eq!(
            BgzfSequenceIndex::<Dna>::open(
                directory.path("source.bgz"),
                directory.path("index.json"),
                BgzfIndexOptions::default()
            )
            .unwrap_err()
            .code(),
            "invalid_index"
        );
    }
}

#[test]
fn stale_source_replacement_and_missing_lookup_are_checked() {
    let directory = Directory::new();
    let mut index = build(&directory, b">r\nAC\n", IndexedFormat::Fasta, 5);
    index.save_new(directory.path("index.json")).unwrap();
    std::fs::rename(directory.path("source.bgz"), directory.path("previous.bgz")).unwrap();
    std::fs::write(directory.path("source.bgz"), encode(b">r\nAC\n", 5)).unwrap();
    for key in ["r", "missing"] {
        assert_eq!(index.get(key).unwrap_err().code(), "stale_source");
    }
    assert_eq!(
        index
            .save_new(directory.path("new.json"))
            .unwrap_err()
            .code(),
        "stale_source"
    );
    assert!(!directory.path("new.json").exists());
    std::fs::write(directory.path("source.bgz"), encode(b">r\nTT\n", 5)).unwrap();
    assert_eq!(
        BgzfSequenceIndex::<Dna>::open(
            directory.path("source.bgz"),
            directory.path("index.json"),
            BgzfIndexOptions::default()
        )
        .unwrap_err()
        .code(),
        "stale_source"
    );
}

#[test]
fn create_only_sidecars_and_compressed_biological_writers_preserve_metadata() {
    let directory = Directory::new();
    let record =
        Record::<Dna>::from_ascii("r description", b"aCgT", Some(vec![0, 20, 40, 93])).unwrap();
    write_bgzf_new(
        directory.path("record.bgz"),
        BgzfWriteOptions::default(),
        |output| write_record_fastq(output, &record),
    )
    .unwrap();
    let mut index = BgzfSequenceIndex::<Dna>::build(
        directory.path("record.bgz"),
        IndexedFormat::FastqSanger,
        BgzfIndexOptions::default(),
    )
    .unwrap();
    let retrieved = index.get("r").unwrap().unwrap();
    assert_eq!(retrieved.header(), record.header());
    assert_eq!(
        retrieved.source_sequence().as_bytes(),
        record.source_sequence().as_bytes()
    );
    assert_eq!(retrieved.qualities(), record.qualities());
    write_gzip_new(
        directory.path("record.gz"),
        GzipWriteOptions::default(),
        |output| write_record_fastq(output, &record),
    )
    .unwrap();
    let records =
        FastqReader::<_, Dna>::new(super::super::open_input(directory.path("record.gz")).unwrap())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
    assert_eq!(
        records[0].source_sequence().as_bytes(),
        record.source_sequence().as_bytes()
    );
    assert_eq!(records[0].qualities(), record.qualities());
    let failed = write_bgzf_new(
        directory.path("loss.bgz"),
        BgzfWriteOptions::default(),
        |output| write_record_fasta(output, &record, LossPolicy::Reject),
    );
    assert!(failed.is_err());
    assert!(!directory.path("loss.bgz").exists());
    index.save_new(directory.path("index.json")).unwrap();
    let before = std::fs::read(directory.path("index.json")).unwrap();
    assert!(index.save_new(directory.path("index.json")).is_err());
    assert_eq!(std::fs::read(directory.path("index.json")).unwrap(), before);
    assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 3);
}

#[test]
fn explicit_rna_and_protein_alphabets_do_not_guess_from_file_extensions() {
    let directory = Directory::new();
    for (content, molecule) in [
        (b">rna\naCuG\n".as_slice(), "rna"),
        (b">protein\nMaFW\n".as_slice(), "protein"),
    ] {
        std::fs::write(directory.path("source.bgz"), encode(content, 4)).unwrap();
        if molecule == "rna" {
            let mut index = BgzfSequenceIndex::<Rna>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                BgzfIndexOptions::default(),
            )
            .unwrap();
            assert_eq!(
                index
                    .get("rna")
                    .unwrap()
                    .unwrap()
                    .source_sequence()
                    .as_bytes(),
                b"aCuG"
            );
        } else {
            let mut index = BgzfSequenceIndex::<Protein>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                BgzfIndexOptions::default(),
            )
            .unwrap();
            assert_eq!(
                index
                    .get("protein")
                    .unwrap()
                    .unwrap()
                    .source_sequence()
                    .as_bytes(),
                b"MaFW"
            );
        }
        assert!(
            BgzfSequenceIndex::<Dna>::build(
                directory.path("source.bgz"),
                IndexedFormat::Fasta,
                BgzfIndexOptions::default()
            )
            .is_err()
        );
    }
}

#[test]
fn retained_native_seqio_indexes_match_headers_case_qualities_raw_spans_and_blocks() {
    let observations = include_bytes!("../../../tests/fixtures/bgzf-index-biopython-1.88.json");
    assert_eq!(
        format!("{:x}", Sha256::digest(observations)),
        "376c24835a9a9697eb750f41470f031210c97471af6883dfbc7872ed2e4c343c"
    );
    let reference: serde_json::Value = serde_json::from_slice(observations).unwrap();
    assert_eq!(reference["reference"]["biopython"], "1.88");
    let fixtures: [&[u8]; 3] = [
        include_bytes!("../../../tests/fixtures/bgzf-native-index.fasta.bgz"),
        include_bytes!("../../../tests/fixtures/bgzf-native-index.fastq.bgz"),
        include_bytes!("../../../tests/fixtures/bgzf-native-boundary.fasta.bgz"),
    ];
    let mut compared = 0;
    for (case, bytes) in reference["cases"].as_array().unwrap().iter().zip(fixtures) {
        let directory = Directory::new();
        std::fs::write(directory.path("source.bgz"), bytes).unwrap();
        let format = if case["format"] == "fasta" {
            IndexedFormat::Fasta
        } else {
            IndexedFormat::FastqSanger
        };
        let mut index = BgzfSequenceIndex::<Dna>::build(
            directory.path("source.bgz"),
            format,
            BgzfIndexOptions::default(),
        )
        .unwrap();
        assert_eq!(
            index.source_sha256(),
            case["source_sha256"].as_str().unwrap()
        );
        assert_eq!(
            index.ids().collect::<Vec<_>>(),
            case["ids"]
                .as_array()
                .unwrap()
                .iter()
                .map(|id| id.as_str().unwrap())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            index.block_count(),
            case["blocks"].as_array().unwrap().len()
        );
        for (block, expected) in index
            .manifest
            .blocks
            .iter()
            .zip(case["blocks"].as_array().unwrap())
        {
            assert_eq!(
                serde_json::json!([
                    block.compressed_offset,
                    block.compressed_size,
                    block.uncompressed_offset,
                    block.uncompressed_size
                ]),
                *expected
            );
        }
        for expected in case["records"].as_array().unwrap() {
            let id = expected["id"].as_str().unwrap();
            let record = index.get(id).unwrap().unwrap();
            assert_eq!(record.header(), expected["header"].as_str().unwrap());
            assert_eq!(
                record.sequence().len() as u64,
                expected["sequence_length"].as_u64().unwrap()
            );
            assert_eq!(
                format!("{:x}", Sha256::digest(record.source_sequence().as_bytes())),
                expected["sequence_sha256"].as_str().unwrap()
            );
            assert_eq!(
                serde_json::to_value(record.qualities()).unwrap(),
                expected["qualities"]
            );
            let entry = index
                .manifest
                .entries
                .iter()
                .find(|entry| entry.id == id)
                .unwrap();
            assert_eq!(entry.byte_length, expected["raw_length"].as_u64().unwrap());
            assert_eq!(entry.sha256, expected["raw_sha256"].as_str().unwrap());
            compared += 1;
        }
        index.save_new(directory.path("index.json")).unwrap();
        let reopened = BgzfSequenceIndex::<Dna>::open(
            directory.path("source.bgz"),
            directory.path("index.json"),
            BgzfIndexOptions::default(),
        )
        .unwrap();
        assert_eq!(reopened.source_sha256(), index.source_sha256());
    }
    assert_eq!(compared, 7);
}

#[test]
fn thousands_of_short_records_share_members_without_changing_raw_span_hashes() {
    let directory = Directory::new();
    let mut input = String::from("\n");
    for ordinal in 0..5_000 {
        use std::fmt::Write;
        writeln!(input, ">r{ordinal} demo\nacgT").unwrap();
    }
    let mut index = build(&directory, input.as_bytes(), IndexedFormat::Fasta, 65_536);
    assert_eq!(index.len(), 5_000);
    assert!(index.block_count() < 10);
    for entry in &index.manifest.entries {
        let start = usize::try_from(entry.offset).unwrap();
        let end = usize::try_from(entry.offset + entry.byte_length).unwrap();
        assert_eq!(
            entry.sha256,
            format!("{:x}", Sha256::digest(&input.as_bytes()[start..end]))
        );
    }
    for id in ["r4999", "r0", "r2500", "r8"] {
        let record = index.get(id).unwrap().unwrap();
        assert_eq!(record.id(), id);
        assert_eq!(record.source_sequence(), "acgT");
    }
    index.save_new(directory.path("index.json")).unwrap();
    let reopened = BgzfSequenceIndex::<Dna>::open(
        directory.path("source.bgz"),
        directory.path("index.json"),
        BgzfIndexOptions::default(),
    )
    .unwrap();
    assert_eq!(reopened.len(), 5_000);
}
