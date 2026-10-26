//! Record-preserving streaming, quality and compression contracts.

use std::io::{self, BufReader, Cursor, Read, Write};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use bioferris::io::{
    FastaReader, FastqReader, IoErrorKind, LossPolicy, ReadLimits, decode_input, one_record,
    write_record_fasta, write_record_fastq,
};
use bioferris::prelude::*;
use flate2::{Compression, write::GzEncoder};

#[test]
fn record_validation_and_case_preservation() {
    let record =
        Record::<Dna>::from_ascii("read-1  mixed case", "aCgT", Some(vec![0, 40, 62, 93])).unwrap();
    assert_eq!(record.id(), "read-1");
    assert_eq!(record.description(), Some("mixed case"));
    assert_eq!(record.header(), "read-1  mixed case");
    assert_eq!(record.source_sequence(), "aCgT");
    assert_eq!(record.sequence().to_string(), "ACGT");
    assert_eq!(record.qualities(), Some([0, 40, 62, 93].as_slice()));
    for (header, sequence, qualities, kind) in [
        (" ", "A", None, IoErrorKind::InvalidHeader),
        ("x\n", "A", None, IoErrorKind::InvalidHeader),
        ("x", "U", None, IoErrorKind::InvalidSequence),
        ("x", "A", Some(vec![94]), IoErrorKind::InvalidQuality),
        ("x", "A", Some(vec![]), IoErrorKind::LengthMismatch),
    ] {
        assert_eq!(
            Record::<Dna>::from_ascii(header, sequence, qualities)
                .unwrap_err()
                .kind(),
            kind
        );
    }
}

#[test]
fn streaming_real_phix174_cds_retains_identity_and_matches_refseq_protein() {
    let record = one_record(FastaReader::<_, Dna>::new(Cursor::new(include_bytes!(
        "fixtures/phix174-e.fasta"
    ))))
    .unwrap();
    assert_eq!(record.id(), "NC_001422.1:568-843");
    assert_eq!(
        record.description(),
        Some("Escherichia phage phiX174, complete genome")
    );
    assert_eq!(record.sequence().len(), 276);
    assert!(record.qualities().is_none());
    // Independent RefSeq NP_040709.1 annotation; fixture provenance is in README.
    let protein = concat!(
        "MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNL",
        "RKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE"
    );
    assert_eq!(record.sequence().translate_cds().unwrap().as_str(), protein);
}

#[test]
fn fasta_matches_legacy_sequence_semantics_without_changing_header_identity() {
    let input = ">a description\r\nac gt\r\n\r\n>b\r\nNN\r\n>empty\r\n";
    let old = parse_fasta::<Dna>(input.as_bytes()).unwrap();
    let records = FastaReader::<_, Dna>::new(Cursor::new(input))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), old.len());
    for (new, old) in records.iter().zip(&old) {
        assert_eq!(new.sequence(), old.sequence());
    }
    assert_eq!(old[0].identifier(), "a description");
    assert_eq!(records[0].id(), "a");
    assert_eq!(records[0].header(), "a description");
    assert_eq!(records[0].source_sequence(), "acgt");
    assert!(FastaReader::<_, Dna>::new(Cursor::new("")).next().is_none());
    assert!(parse_fasta::<Dna>(b"".as_slice()).is_err());
}

struct Counted {
    input: Cursor<Vec<u8>>,
    bytes: Arc<AtomicUsize>,
}
impl Read for Counted {
    fn read(&mut self, target: &mut [u8]) -> io::Result<usize> {
        let length = self.input.read(target)?;
        self.bytes.fetch_add(length, Ordering::SeqCst);
        Ok(length)
    }
}

#[test]
fn fasta_is_lazy_and_errors_are_fused_with_source_context() {
    let bytes = Arc::new(AtomicUsize::new(0));
    let input = Counted {
        input: Cursor::new(b">first\nAC\n>bad\nU\n>last\nGT\n".to_vec()),
        bytes: bytes.clone(),
    };
    let mut records = FastaReader::<_, Dna>::new(BufReader::with_capacity(1, input));
    assert_eq!(bytes.load(Ordering::SeqCst), 0);
    assert_eq!(records.next().unwrap().unwrap().id(), "first");
    assert_eq!(bytes.load(Ordering::SeqCst), 15);
    let error = records.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::InvalidSequence);
    assert_eq!((error.record(), error.line()), (Some(2), Some(4)));
    assert!(records.next().is_none());
    assert!(records.next().is_none());
    let mut bad_header = FastaReader::<_, Dna>::new(Cursor::new(">ok\nA\n>\nG\n"));
    assert!(bad_header.next().unwrap().is_ok());
    assert_eq!(
        bad_header.next().unwrap().unwrap_err().kind(),
        IoErrorKind::InvalidHeader
    );
}

#[test]
fn fastq_wrapping_is_length_delimited_and_round_trips_scores() {
    let input =
        "@read description\r\naC\r\ngT\r\n+read description\r\n@+\r\n!~\r\n@second\nN\n+\nI\n";
    let records = FastqReader::<_, Dna>::new(Cursor::new(input))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].source_sequence(), "aCgT");
    assert_eq!(records[0].qualities(), Some([31, 10, 0, 93].as_slice()));
    let mut written = Vec::new();
    for record in &records {
        write_record_fastq(&mut written, record).unwrap();
    }
    let reloaded = FastqReader::<_, Dna>::new(Cursor::new(&written))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records, reloaded);
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "@read description\naCgT\n+\n@+!~\n@second\nN\n+\nI\n"
    );
}

#[test]
fn fastq_malformed_inputs_fail_once() {
    for (input, kind) in [
        ("read\nA\n+\n!\n", IoErrorKind::InvalidHeader),
        ("@r\nA\n+other\n!\n", IoErrorKind::InvalidHeader),
        ("@r\nA\n", IoErrorKind::Truncated),
        ("@r\nAA\n+\n!\n", IoErrorKind::Truncated),
        ("@r\nA\n+\n!!\n", IoErrorKind::LengthMismatch),
        ("@r\nA\n+\n \n", IoErrorKind::InvalidQuality),
        ("@r\nA C\n+\n!!!\n", IoErrorKind::InvalidSequence),
        ("@r\n+\n\n", IoErrorKind::Unsupported),
        ("@r\nA\n+\n\n", IoErrorKind::InvalidQuality),
    ] {
        let mut records = FastqReader::<_, Dna>::new(Cursor::new(input));
        let error = records.next().unwrap().unwrap_err();
        assert_eq!(error.kind(), kind, "{input:?}: {error}");
        assert_eq!(error.record(), Some(1));
        assert!(error.line().is_some());
        assert!(records.next().is_none());
    }
}

#[test]
fn invalid_utf8_lines_fail_before_their_symbols_are_checked() {
    // 0xff and a lone 0xc3 are invalid UTF-8; `J` would be an invalid DNA symbol.
    for (fasta, input, line, message) in [
        (
            true,
            &b">r\xff\nA\n"[..],
            1,
            "invalid utf-8 sequence of 1 bytes from index 2",
        ),
        (
            true,
            &b">r\nAC\xffJ\n"[..],
            2,
            "invalid utf-8 sequence of 1 bytes from index 2",
        ),
        (
            true,
            &b">r\nA\n>s\nAJ\xc3"[..],
            4,
            "incomplete utf-8 byte sequence from index 2",
        ),
        (
            false,
            &b"@r\nA\xff\n+\n!!\n"[..],
            2,
            "invalid utf-8 sequence of 1 bytes from index 1",
        ),
        (
            false,
            &b"@r\nA\n+\xff\n!\n"[..],
            3,
            "invalid utf-8 sequence of 1 bytes from index 1",
        ),
        (
            false,
            &b"@r\nA\n+\n\xff\n"[..],
            4,
            "invalid utf-8 sequence of 1 bytes from index 0",
        ),
    ] {
        let error = if fasta {
            FastaReader::<_, Dna>::new(Cursor::new(input))
                .find_map(Result::err)
                .unwrap()
        } else {
            FastqReader::<_, Dna>::new(Cursor::new(input))
                .next()
                .unwrap()
                .unwrap_err()
        };
        assert_eq!(
            error.kind(),
            IoErrorKind::InvalidStructure,
            "{input:?}: {error}"
        );
        assert_eq!(error.line(), Some(line), "{input:?}");
        assert!(error.to_string().contains(message), "{input:?}: {error}");
    }
}

#[test]
fn explicit_limits_apply_before_large_record_and_line_accumulation() {
    for limits in [
        ReadLimits {
            max_line_bytes: 3,
            ..ReadLimits::default()
        },
        ReadLimits {
            max_total_bytes: 5,
            ..ReadLimits::default()
        },
        ReadLimits {
            max_record_symbols: 2,
            ..ReadLimits::default()
        },
    ] {
        let error = FastaReader::<_, Dna>::with_limits(Cursor::new(">r\nACG\n"), limits)
            .next()
            .unwrap()
            .unwrap_err();
        assert_eq!(error.kind(), IoErrorKind::LimitExceeded);
    }
    let limits = ReadLimits {
        max_record_symbols: 2,
        ..ReadLimits::default()
    };
    assert_eq!(
        FastqReader::<_, Dna>::with_limits(Cursor::new("@r\nACG\n+\n!!!\n"), limits)
            .next()
            .unwrap()
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

struct ShortReads {
    bytes: Cursor<Vec<u8>>,
    interrupted: bool,
}
impl Read for ShortReads {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if !self.interrupted {
            self.interrupted = true;
            return Err(io::Error::from(io::ErrorKind::Interrupted));
        }
        let length = output.len().min(1);
        self.bytes.read(&mut output[..length])
    }
}

#[test]
fn compression_detection_handles_short_reads_interruptions_and_short_files() {
    let reader = ShortReads {
        bytes: Cursor::new(gzip(b">short\nAC\n")),
        interrupted: false,
    };
    let record = one_record(FastaReader::<_, Dna>::new(decode_input(reader).unwrap())).unwrap();
    assert_eq!(record.sequence().to_string(), "AC");
    for input in [b"".as_slice(), b">", b"\n"] {
        let reader = ShortReads {
            bytes: Cursor::new(input.to_vec()),
            interrupted: false,
        };
        let mut output = String::new();
        decode_input(reader)
            .unwrap()
            .read_to_string(&mut output)
            .unwrap();
        assert_eq!(output.as_bytes(), input);
    }
}

#[test]
fn gzip_magic_members_crc_and_truncation_are_checked() {
    let mut compressed = gzip(b">a\nAC\n");
    compressed.extend(gzip(b">b\nGT\n"));
    let reader = decode_input(Cursor::new(compressed)).unwrap();
    let records = FastaReader::<_, Dna>::new(reader)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 2);
    for compressed in [
        {
            let mut data = gzip(b">a\nAC\n");
            let index = data.len() - 8;
            data[index] ^= 1;
            data
        },
        {
            let mut data = gzip(b">a\nAC\n");
            data.truncate(data.len() - 4);
            data
        },
    ] {
        let reader = decode_input(Cursor::new(compressed)).unwrap();
        let error = FastaReader::<_, Dna>::new(reader)
            .collect::<Result<Vec<_>, _>>()
            .unwrap_err();
        assert_eq!(error.kind(), IoErrorKind::Io);
    }
    let reader = decode_input(Cursor::new(gzip(b">large\nAAAAAAAAAAAA\n"))).unwrap();
    let limits = ReadLimits {
        max_total_bytes: 10,
        ..ReadLimits::default()
    };
    assert_eq!(
        FastaReader::<_, Dna>::with_limits(reader, limits)
            .next()
            .unwrap()
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
}

#[test]
fn singleton_and_explicit_loss_contracts() {
    let read = |input| one_record(FastaReader::<_, Dna>::new(Cursor::new(input)));
    assert!(read(">a\nAC\n").is_ok());
    assert_eq!(read("").unwrap_err().kind(), IoErrorKind::LengthMismatch);
    assert_eq!(
        read(">a\nAC\n>b\nGT\n").unwrap_err().kind(),
        IoErrorKind::LengthMismatch
    );
    assert_eq!(
        read(">a\nAC\n>b\nU\n").unwrap_err().kind(),
        IoErrorKind::InvalidSequence
    );
    let record = Record::<Dna>::from_ascii("a desc", "ac", Some(vec![0, 40])).unwrap();
    let mut output = Vec::new();
    assert_eq!(
        write_record_fasta(&mut output, &record, LossPolicy::Reject)
            .unwrap_err()
            .kind(),
        IoErrorKind::DataLoss
    );
    assert!(output.is_empty());
    write_record_fasta(&mut output, &record, LossPolicy::Allow).unwrap();
    assert_eq!(output, b">a desc\nac\n");
    let plain = read(">a\nAC\n").unwrap();
    let mut output = Vec::new();
    assert_eq!(
        write_record_fastq(&mut output, &plain).unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    assert!(output.is_empty());
}

proptest::proptest! {
    #[test]
    fn arbitrary_dna_fastq_roundtrip(sequence in "[ACGTacgt]{1,200}", score in 0u8..=93) {
        let record = Record::<Dna>::from_ascii("read metadata", &sequence, Some(vec![score; sequence.len()])).unwrap();
        let mut output = Vec::new();
        write_record_fastq(&mut output, &record).unwrap();
        let reloaded = one_record(FastqReader::<_, Dna>::new(Cursor::new(output))).unwrap();
        proptest::prop_assert_eq!(record, reloaded);
    }
}
