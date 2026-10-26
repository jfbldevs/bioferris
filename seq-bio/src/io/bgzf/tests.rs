//! Independent format fixtures and public reference interoperability checks.

use super::*;
use std::cell::RefCell;
use std::io::Cursor;
use std::rc::Rc;

fn independent_crc(bytes: &[u8]) -> u32 {
    let mut crc = !0_u32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ if crc & 1 == 1 { 0xedb8_8320 } else { 0 };
        }
    }
    !crc
}

fn member(data: &[u8], deflate: &[u8]) -> Vec<u8> {
    let size = 18 + deflate.len() + 8;
    let mut bytes = vec![31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0];
    bytes.extend_from_slice(&u16::try_from(size - 1).unwrap().to_le_bytes());
    bytes.extend_from_slice(deflate);
    bytes.extend_from_slice(&independent_crc(data).to_le_bytes());
    bytes.extend_from_slice(&u32::try_from(data.len()).unwrap().to_le_bytes());
    bytes
}

fn stored(data: &[u8]) -> Vec<u8> {
    let length = u16::try_from(data.len()).unwrap();
    let mut deflate = vec![1]; // Final RFC 1951 stored block with zero padding.
    deflate.extend_from_slice(&length.to_le_bytes());
    deflate.extend_from_slice(&(!length).to_le_bytes());
    deflate.extend_from_slice(data);
    member(data, &deflate)
}

fn with_marker(mut data: Vec<u8>) -> Vec<u8> {
    data.extend_from_slice(&BGZF_EOF_MARKER);
    data
}

fn decoded(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    let mut reader = BgzfReader::new(bytes);
    let mut output = Vec::new();
    let mut buffer = [0; 8192];
    loop {
        let amount = reader.read_bytes(&mut buffer)?;
        if amount == 0 {
            return Ok(output);
        }
        output.extend_from_slice(&buffer[..amount]);
    }
}

fn compressed(data: &[u8], options: BgzfWriteOptions) -> Vec<u8> {
    let mut writer = BgzfWriter::with_options(Vec::new(), options).unwrap();
    writer.write_all(data).unwrap();
    writer.finish().unwrap()
}

fn random_bytes(length: usize) -> Vec<u8> {
    let mut state = 0xa341_316c_u32;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state.to_le_bytes()[0]
        })
        .collect()
}

#[test]
fn independent_stored_deflate_crc_and_empty_marker_match_the_specification() {
    assert_eq!(independent_crc(b"123456789"), 0xcbf4_3926);
    let bytes = with_marker(stored(b"123456789\0\r\n"));
    assert_eq!(decoded(&bytes).unwrap(), b"123456789\0\r\n");
    assert_eq!(decoded(&BGZF_EOF_MARKER).unwrap(), b"");
    let mut reader = BgzfReader::new(bytes.as_slice());
    assert_eq!(reader.eof_status(), BgzfEofStatus::Unchecked);
    assert_eq!(reader.finish().unwrap(), BgzfEofStatus::Present);
    assert!(reader.current_block().unwrap().canonical_eof_marker);
}

#[test]
fn virtual_offsets_match_public_biopython_observations_and_reject_overflow() {
    assert_eq!(
        BgzfVirtualOffset::new(55_074, 126).unwrap().raw(),
        3_609_329_790
    );
    let maximum = BgzfVirtualOffset::new(BGZF_MAX_COMPRESSED_OFFSET, 65_535).unwrap();
    assert_eq!(maximum.raw(), u64::MAX);
    assert_eq!(maximum.uncompressed_offset(), 65_535);
    assert_eq!(maximum.compressed_offset(), BGZF_MAX_COMPRESSED_OFFSET);
    assert_eq!(u64::from(maximum), u64::MAX);
    assert_eq!(BgzfVirtualOffset::from_raw(u64::MAX), maximum);
    assert!(matches!(
        BgzfVirtualOffset::new(1 << 48, 0),
        Err(BgzfError::InvalidVirtualOffset { .. })
    ));
    assert!(matches!(
        BgzfVirtualOffset::new(0, 65_536),
        Err(BgzfError::InvalidVirtualOffset { .. })
    ));
    assert!(BgzfVirtualOffset::new(1, 0).unwrap() > BgzfVirtualOffset::new(0, 65_535).unwrap());
}

#[test]
fn strict_and_legacy_eof_policies_are_explicit() {
    assert!(matches!(
        decoded(b""),
        Err(BgzfError::MissingEof { offset: 0 })
    ));
    let bytes = stored(b"ABC");
    let mut strict = BgzfReader::new(bytes.as_slice());
    let mut output = [0; 10];
    assert_eq!(strict.read_bytes(&mut output).unwrap(), 3);
    assert_eq!(&output[..3], b"ABC");
    assert!(matches!(
        strict.read_bytes(&mut output),
        Err(BgzfError::MissingEof { offset: 34 })
    ));
    assert_eq!(strict.eof_status(), BgzfEofStatus::Missing);
    assert!(matches!(
        strict.virtual_position(),
        Err(BgzfError::FailedState)
    ));
    for input in [b"".as_slice(), bytes.as_slice()] {
        let mut reader = BgzfReader::with_options(
            input,
            BgzfReadOptions {
                eof_policy: BgzfEofPolicy::AllowMissing,
            },
        );
        assert_eq!(reader.finish().unwrap(), BgzfEofStatus::Missing);
    }
}

#[test]
fn embedded_empty_members_never_hide_later_data_or_invalid_suffixes() {
    let mut bytes = with_marker(stored(b"first"));
    bytes.extend_from_slice(&stored(b"second"));
    bytes.extend_from_slice(&BGZF_EOF_MARKER);
    assert_eq!(decoded(&bytes).unwrap(), b"firstsecond");
    bytes.extend_from_slice(&[0; 20]);
    assert!(matches!(
        decoded(&bytes),
        Err(BgzfError::InvalidHeader { .. })
    ));
    let empty = with_marker(with_marker(stored(b"")));
    assert_eq!(decoded(&empty).unwrap(), b"");
}

#[test]
fn every_partial_member_is_truncated_without_exposing_its_bytes() {
    let bytes = stored(b"ACGT");
    for length in 1..bytes.len() {
        let mut reader = BgzfReader::new(&bytes[..length]);
        let mut output = [0xa5; 100];
        assert!(
            matches!(
                reader.read_bytes(&mut output),
                Err(BgzfError::Truncated { offset: 0, .. })
            ),
            "prefix length={length}"
        );
        assert_eq!(output, [0xa5; 100]);
        assert!(matches!(
            reader.read_bytes(&mut output),
            Err(BgzfError::FailedState)
        ));
    }
    let mut reader = BgzfReader::new(&bytes[..3]);
    assert!(matches!(
        reader.finish(),
        Err(BgzfError::Truncated {
            part: "header",
            expected: 11,
            actual: 2,
            ..
        })
    ));
}

#[test]
fn invalid_magic_method_flags_and_zero_files_are_refused() {
    let good = stored(b"abc");
    for (index, value) in [(0, 0), (1, 0), (2, 0), (3, 0), (3, 12), (3, 5)] {
        let mut bytes = good.clone();
        bytes[index] = value;
        assert!(matches!(
            decoded(&bytes),
            Err(BgzfError::InvalidHeader { offset: 0, .. })
        ));
    }
    assert!(matches!(
        decoded(&[0; 64]),
        Err(BgzfError::InvalidHeader {
            observed: [0, 0, 0, 0],
            ..
        })
    ));
}

fn replace_extra(original: &[u8], extra: &[u8]) -> Vec<u8> {
    let size = original.len() - 6 + extra.len();
    let mut header = original[..12].to_vec();
    header[10..12].copy_from_slice(&u16::try_from(extra.len()).unwrap().to_le_bytes());
    header.extend_from_slice(extra);
    if let Some(index) = extra
        .windows(4)
        .position(|window| window == b"BC\x02\x00")
        .filter(|index| index + 6 <= extra.len())
    {
        header[12 + index + 4..12 + index + 6]
            .copy_from_slice(&u16::try_from(size - 1).unwrap().to_le_bytes());
    }
    header.extend_from_slice(&original[18..]);
    header
}

#[test]
fn extra_fields_support_other_subfields_and_refuse_bad_bc_structure() {
    let good = stored(b"abc");
    let extra = b"XY\x03\x00oneBC\x02\x00\x00\x00ZZ\x00\x00";
    assert_eq!(
        decoded(&with_marker(replace_extra(&good, extra))).unwrap(),
        b"abc"
    );
    for extra in [
        b"".as_slice(),
        b"XY\x00\x00",
        b"B",
        b"BC\x03\x00abc",
        b"BC\x02\x00\x00",
        b"BC\x02\x00\x00\x00BC\x02\x00\x00\x00",
    ] {
        assert!(
            matches!(
                decoded(&replace_extra(&good, extra)),
                Err(BgzfError::InvalidExtra { offset: 0, .. })
            ),
            "extra={extra:?}"
        );
    }
    let mut huge = good.clone();
    huge[10..12].copy_from_slice(&u16::MAX.to_le_bytes());
    assert!(matches!(
        decoded(&huge),
        Err(BgzfError::InvalidBlockSize { .. })
    ));
}

#[test]
fn compressed_size_and_payload_consumption_are_checked_exactly() {
    let mut too_small = stored(b"abc");
    too_small[16..18].copy_from_slice(&0_u16.to_le_bytes());
    assert!(matches!(
        decoded(&too_small),
        Err(BgzfError::InvalidBlockSize {
            compressed_size: 1,
            ..
        })
    ));
    let mut too_large = stored(b"abc");
    too_large[16..18].copy_from_slice(&u16::MAX.to_le_bytes());
    assert!(matches!(
        decoded(&too_large),
        Err(BgzfError::Truncated {
            part: "payload/trailer",
            ..
        })
    ));
    let valid = stored(b"abc");
    let mut payload = valid[18..valid.len() - 8].to_vec();
    payload.push(0);
    assert!(matches!(
        decoded(&member(b"abc", &payload)),
        Err(BgzfError::Deflate { .. })
    ));
    payload.truncate(payload.len() - 2);
    assert!(matches!(
        decoded(&member(b"abc", &payload)),
        Err(BgzfError::Deflate { .. })
    ));
    assert!(matches!(
        decoded(&member(b"", &[7])),
        Err(BgzfError::Deflate { .. })
    ));
}

#[test]
fn crc_and_isize_errors_retain_actual_and_declared_values() {
    let mut bytes = stored(b"abc");
    let trailer = bytes.len() - 8;
    bytes[trailer] ^= 1;
    assert!(
        matches!(decoded(&bytes), Err(BgzfError::CrcMismatch { offset: 0, expected, actual }) if actual == 0x3524_41c2 && expected == actual ^ 1)
    );
    bytes[trailer] ^= 1;
    bytes[trailer + 4..].copy_from_slice(&4_u32.to_le_bytes());
    assert!(matches!(
        decoded(&bytes),
        Err(BgzfError::SizeMismatch {
            offset: 0,
            expected: 4,
            actual: 3
        })
    ));
    bytes[trailer + 4..].copy_from_slice(&65_537_u32.to_le_bytes());
    assert!(matches!(
        decoded(&bytes),
        Err(BgzfError::UncompressedLimit {
            offset: 0,
            observed: 65_537
        })
    ));
}

#[test]
fn decompression_expansion_is_bounded_independently_of_forged_isize() {
    let data = vec![0; 200_000];
    let mut deflater = flate2::write::DeflateEncoder::new(Vec::new(), Compression::default());
    deflater.write_all(&data).unwrap();
    let payload = deflater.finish().unwrap();
    let bytes = member(b"", &payload); // Falsely small ISIZE cannot defeat expansion bounds.
    assert!(matches!(
        decoded(&bytes),
        Err(BgzfError::UncompressedLimit {
            observed: 65_537,
            ..
        })
    ));
}

#[test]
fn maximum_compressed_and_uncompressed_blocks_are_supported() {
    let exact_compressed = stored(&vec![b'G'; 65_505]);
    assert_eq!(exact_compressed.len(), 65_536);
    assert_eq!(
        decoded(&with_marker(exact_compressed)).unwrap().len(),
        65_505
    );
    let data = vec![b'A'; 65_536];
    let bytes = compressed(
        &data,
        BgzfWriteOptions {
            block_size: 65_536,
            ..Default::default()
        },
    );
    let mut reader = BgzfReader::new(bytes.as_slice());
    let mut first = vec![0; 65_535];
    reader.read_exact(&mut first).unwrap();
    assert_eq!(
        reader.virtual_position().unwrap().uncompressed_offset(),
        65_535
    );
    assert_eq!(reader.current_block().unwrap().uncompressed_size, 65_536);
    reader.read_exact(&mut [0]).unwrap();
    assert_eq!(reader.virtual_position().unwrap().uncompressed_offset(), 0);
    assert_eq!(
        reader.virtual_position().unwrap().compressed_offset(),
        reader.current_block().unwrap().compressed_size as u64
    );
    assert_eq!(reader.finish().unwrap(), BgzfEofStatus::Present);
}

#[test]
fn seek_validates_targets_normalizes_boundaries_and_recovers_failed_state() {
    let first = stored(b"abc");
    let second = stored(b"DEF");
    let mut bytes = first.clone();
    bytes.extend_from_slice(&second);
    let bytes = with_marker(bytes);
    let mut reader = BgzfReader::new(Cursor::new(bytes.clone()));
    let expected = BgzfVirtualOffset::new(first.len() as u64, 0).unwrap();
    assert_eq!(
        reader
            .seek_virtual(BgzfVirtualOffset::new(0, 3).unwrap())
            .unwrap(),
        expected
    );
    let mut out = [0; 2];
    reader.read_exact(&mut out).unwrap();
    assert_eq!(&out, b"DE");
    assert!(matches!(
        reader.seek_virtual(BgzfVirtualOffset::new(0, 4).unwrap()),
        Err(BgzfError::InvalidSeek {
            block_length: Some(3),
            ..
        })
    ));
    assert!(matches!(
        reader.read_bytes(&mut out),
        Err(BgzfError::FailedState)
    ));
    assert_eq!(
        reader
            .seek_virtual(BgzfVirtualOffset::new(0, 1).unwrap())
            .unwrap()
            .raw(),
        1
    );
    reader.read_exact(&mut out).unwrap();
    assert_eq!(&out, b"bc");
    assert!(matches!(
        reader.seek_virtual(BgzfVirtualOffset::new(bytes.len() as u64 + 1, 0).unwrap()),
        Err(BgzfError::InvalidSeek { .. })
    ));
    assert!(
        reader
            .seek_virtual(BgzfVirtualOffset::new(1, 0).unwrap())
            .is_err()
    );
    let end = BgzfVirtualOffset::new(bytes.len() as u64, 0).unwrap();
    assert_eq!(reader.seek_virtual(end).unwrap(), end);
    assert_eq!(reader.read_bytes(&mut out).unwrap(), 0);
    assert_eq!(reader.eof_status(), BgzfEofStatus::Present);
}

#[test]
fn any_valid_final_empty_member_is_supported_including_eof_seeks() {
    let mut empty = stored(b"");
    empty[4] = 42; // Valid noncanonical MTIME and stored rather than fixed DEFLATE.
    let empty = replace_extra(&empty, b"XY\x01\x00xBC\x02\x00\x00\x00");
    let mut bytes = stored(b"abc");
    bytes.extend_from_slice(&empty);
    assert_eq!(decoded(&bytes).unwrap(), b"abc");
    let mut reader = BgzfReader::new(Cursor::new(bytes.clone()));
    reader
        .seek_virtual(BgzfVirtualOffset::new(bytes.len() as u64, 0).unwrap())
        .unwrap();
    assert_eq!(reader.finish().unwrap(), BgzfEofStatus::Present);
}

#[test]
fn finish_after_seek_validates_the_suffix_without_claiming_unread_prefix() {
    let mut invalid_prefix = stored(b"bad");
    let crc = invalid_prefix.len() - 8;
    invalid_prefix[crc] ^= 1;
    let valid_start = invalid_prefix.len();
    invalid_prefix.extend_from_slice(&stored(b"good"));
    let bytes = with_marker(invalid_prefix);
    assert!(matches!(
        decoded(&bytes),
        Err(BgzfError::CrcMismatch { .. })
    ));
    let mut reader = BgzfReader::new(Cursor::new(bytes));
    reader
        .seek_virtual(BgzfVirtualOffset::new(valid_start as u64, 0).unwrap())
        .unwrap();
    assert_eq!(reader.finish().unwrap(), BgzfEofStatus::Present);
}

struct ShortReader {
    bytes: Cursor<Vec<u8>>,
    calls: usize,
}
impl Read for ShortReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        self.calls += 1;
        if self.calls % 7 == 0 {
            return Err(io::ErrorKind::Interrupted.into());
        }
        let amount = output.len().min(1);
        self.bytes.read(&mut output[..amount])
    }
}

#[test]
fn read_only_short_interrupted_streams_and_bufread_use_the_same_checked_core() {
    let bytes = with_marker(stored(b"line1\nline2\r\n\0last"));
    let input = ShortReader {
        bytes: Cursor::new(bytes),
        calls: 0,
    };
    let mut reader = BgzfReader::new(input);
    let mut line = String::new();
    assert_eq!(reader.read_line(&mut line).unwrap(), 6);
    assert_eq!(line, "line1\n");
    assert_eq!(reader.fill_buf().unwrap(), b"line2\r\n\0last");
    reader.consume(7);
    let mut rest = Vec::new();
    reader.read_to_end(&mut rest).unwrap();
    assert_eq!(rest, b"\0last");
    assert_eq!(reader.eof_status(), BgzfEofStatus::Present);
    assert_eq!(reader.read_bytes(&mut []).unwrap(), 0);
}

#[test]
fn standard_io_errors_retain_typed_integrity_causes_and_do_not_leak_bad_bytes() {
    let mut bytes = stored(b"abc");
    bytes[20] ^= 1;
    let mut reader = BgzfReader::new(bytes.as_slice());
    let mut output = Vec::new();
    let error = reader.read_to_end(&mut output).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert!(
        error
            .get_ref()
            .unwrap()
            .downcast_ref::<BgzfError>()
            .is_some()
    );
    assert!(output.is_empty());
    let mut reader = BgzfReader::new([31_u8].as_slice());
    let error = reader.fill_buf().unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
    assert!(matches!(
        error.get_ref().unwrap().downcast_ref::<BgzfError>(),
        Some(BgzfError::Truncated { offset: 0, .. })
    ));
}

#[derive(Clone, Default)]
struct Sink {
    bytes: Rc<RefCell<Vec<u8>>>,
    fail_after: Option<usize>,
    flush_error: bool,
    calls: usize,
}
impl Write for Sink {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.calls += 1;
        if self.calls % 11 == 0 {
            return Err(io::ErrorKind::Interrupted.into());
        }
        let length = self.bytes.borrow().len();
        if self.fail_after.is_some_and(|limit| length >= limit) {
            return Err(io::ErrorKind::PermissionDenied.into());
        }
        let count = bytes
            .len()
            .min(3)
            .min(self.fail_after.map_or(usize::MAX, |limit| limit - length));
        self.bytes.borrow_mut().extend_from_slice(&bytes[..count]);
        Ok(count)
    }
    fn flush(&mut self) -> io::Result<()> {
        if self.flush_error {
            Err(io::ErrorKind::BrokenPipe.into())
        } else {
            Ok(())
        }
    }
}

#[test]
fn writers_finalize_empty_outputs_and_preserve_all_bytes_with_short_writes() {
    assert_eq!(
        BgzfWriter::new(Vec::new()).finish().unwrap(),
        BGZF_EOF_MARKER
    );
    let input = random_bytes(200_003);
    let sink = Sink::default();
    let captured = sink.bytes.clone();
    let mut writer = BgzfWriter::new(sink);
    assert_eq!(writer.write_bytes(&[]).unwrap(), 0);
    writer.write_all(&input).unwrap();
    writer.finish().unwrap();
    let bytes = captured.borrow();
    assert_eq!(decoded(&bytes).unwrap(), input);
    let mut gunzip = flate2::read::MultiGzDecoder::new(bytes.as_slice());
    let mut independently_decoded = Vec::new();
    gunzip.read_to_end(&mut independently_decoded).unwrap();
    assert_eq!(independently_decoded, input);
}

#[test]
fn every_compression_level_and_extreme_candidate_size_is_supported() {
    for level in 0..=9 {
        let input = random_bytes(65_537);
        let bytes = compressed(
            &input,
            BgzfWriteOptions {
                block_size: 65_536,
                compression_level: level,
            },
        );
        assert_eq!(decoded(&bytes).unwrap(), input, "level={level}");
        let mut position = 0;
        let mut slice = bytes.as_slice();
        let mut total = 0;
        while let Some(block) = read_block(&mut slice, &mut position).unwrap() {
            assert!(block.info.compressed_size <= 65_536);
            assert!(block.info.uncompressed_size <= 65_536);
            total += block.data.len();
        }
        assert_eq!(total, input.len());
    }
    let input: Vec<_> = (0..=255).collect();
    let bytes = compressed(
        &input,
        BgzfWriteOptions {
            block_size: 1,
            ..Default::default()
        },
    );
    assert_eq!(decoded(&bytes).unwrap(), input);
}

#[test]
fn flushing_positions_are_stable_and_drop_never_performs_hidden_finalization() {
    let sink = Sink::default();
    let captured = sink.bytes.clone();
    let mut writer = BgzfWriter::new(sink);
    writer.write_all(b"first").unwrap();
    assert!(captured.borrow().is_empty());
    let position = writer.virtual_position().unwrap();
    assert_eq!(position.uncompressed_offset(), 0);
    assert_eq!(position.compressed_offset(), captured.borrow().len() as u64);
    assert!(matches!(
        decoded(&captured.borrow()),
        Err(BgzfError::MissingEof { .. })
    ));
    writer.write_all(b"second").unwrap();
    writer.finish().unwrap();
    let mut reader = BgzfReader::new(Cursor::new(captured.borrow().clone()));
    reader.seek_virtual(position).unwrap();
    let mut remaining = Vec::new();
    reader.read_to_end(&mut remaining).unwrap();
    assert_eq!(remaining, b"second");
    let sink = Sink::default();
    let captured = sink.bytes.clone();
    {
        let mut writer = BgzfWriter::new(sink);
        writer.write_all(b"unfinalized").unwrap();
    }
    assert!(captured.borrow().is_empty());
}

#[test]
fn writer_options_write_failures_and_flush_failures_are_observable() {
    for options in [
        BgzfWriteOptions {
            compression_level: 10,
            ..Default::default()
        },
        BgzfWriteOptions {
            block_size: 0,
            ..Default::default()
        },
        BgzfWriteOptions {
            block_size: 65_537,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            BgzfWriter::with_options(Vec::new(), options),
            Err(BgzfError::InvalidOptions { .. })
        ));
    }
    let sink = Sink {
        fail_after: Some(5),
        ..Default::default()
    };
    let mut writer = BgzfWriter::new(sink);
    writer.write_all(b"data").unwrap();
    assert!(matches!(
        writer.flush_blocks(),
        Err(BgzfError::Io {
            operation: "write block",
            offset: 0,
            ..
        })
    ));
    assert!(matches!(
        writer.write_bytes(b"more"),
        Err(BgzfError::FailedState)
    ));
    let sink = Sink {
        flush_error: true,
        ..Default::default()
    };
    assert!(matches!(
        BgzfWriter::new(sink).finish(),
        Err(BgzfError::Io {
            operation: "finalize output",
            ..
        })
    ));
    let sink = Sink {
        flush_error: true,
        ..Default::default()
    };
    let mut writer = BgzfWriter::new(sink);
    let error = writer.flush().unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
    assert!(matches!(
        error.get_ref().unwrap().downcast_ref::<BgzfError>(),
        Some(BgzfError::Io { .. })
    ));
    assert!(matches!(
        writer.virtual_position(),
        Err(BgzfError::FailedState)
    ));
}

#[test]
fn explicit_origins_and_checked_address_ends_do_not_wrap() {
    let mut reader = BgzfReader::with_origin(
        BGZF_EOF_MARKER.as_slice(),
        BGZF_MAX_COMPRESSED_OFFSET - 27,
        BgzfReadOptions::default(),
    )
    .unwrap();
    assert!(matches!(reader.finish(), Err(BgzfError::AddressOverflow)));
    let writer = BgzfWriter::with_origin(
        Vec::new(),
        BGZF_MAX_COMPRESSED_OFFSET - 27,
        BgzfWriteOptions::default(),
    )
    .unwrap();
    assert!(matches!(writer.finish(), Err(BgzfError::AddressOverflow)));
    let writer = BgzfWriter::with_origin(
        Vec::new(),
        BGZF_MAX_COMPRESSED_OFFSET - 28,
        BgzfWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(writer.finish().unwrap(), BGZF_EOF_MARKER);
    assert!(BgzfReader::with_origin(b"".as_slice(), 1 << 48, BgzfReadOptions::default()).is_err());
    assert!(BgzfWriter::with_origin(Vec::new(), 1 << 48, BgzfWriteOptions::default()).is_err());
    let bytes = with_marker(stored(b"abc"));
    let mut prefixed = vec![0; 10];
    prefixed.extend_from_slice(&bytes);
    let mut source = Cursor::new(prefixed);
    source.set_position(10);
    let mut reader = BgzfReader::with_origin(source, 10, BgzfReadOptions::default()).unwrap();
    assert_eq!(reader.virtual_position().unwrap().compressed_offset(), 10);
    reader
        .seek_virtual(BgzfVirtualOffset::new(10, 1).unwrap())
        .unwrap();
    let mut remainder = Vec::new();
    reader.read_to_end(&mut remainder).unwrap();
    assert_eq!(remainder, b"bc");
}

#[test]
fn public_biopython_bam_fixture_matches_observed_blocks_bytes_and_positions() {
    let bytes = include_bytes!("../../../tests/fixtures/bgzf-biopython-ex1.bam");
    let data = decoded(bytes).unwrap();
    assert_eq!(data.len(), 456_614);
    assert_eq!(independent_crc(&data), 1_043_922_902);
    assert_eq!(&data[..4], b"BAM\x01");
    let expected = [
        (0, 18_239, 65_536),
        (18_239, 18_223, 65_536),
        (36_462, 18_017, 65_536),
        (54_479, 17_342, 65_536),
        (71_821, 17_715, 65_536),
        (89_536, 17_728, 65_536),
        (107_264, 17_292, 63_398),
        (124_556, 28, 0),
    ];
    let mut slice = bytes.as_slice();
    let mut position = 0;
    for (offset, compressed, uncompressed) in expected {
        let block = read_block(&mut slice, &mut position).unwrap().unwrap();
        assert_eq!(
            (
                block.info.compressed_offset,
                block.info.compressed_size,
                block.data.len()
            ),
            (offset, compressed, uncompressed)
        );
    }
    assert!(read_block(&mut slice, &mut position).unwrap().is_none());
    let mut reader = BgzfReader::new(Cursor::new(bytes));
    for (amount, expected_position) in [
        (1, 1),
        (65_534, 65_535),
        (1, 1_195_311_104),
        (1, 1_195_311_105),
        (65_535, 2_389_573_632),
        (1, 2_389_573_633),
    ] {
        reader.read_exact(&mut vec![0; amount]).unwrap();
        assert_eq!(reader.virtual_position().unwrap().raw(), expected_position);
    }
    reader.finish().unwrap();
    assert_eq!(reader.virtual_position().unwrap().raw(), 8_164_737_024);
}

#[test]
fn native_biopython_writer_fixture_is_decoded_byte_for_byte() {
    let bytes = include_bytes!("../../../tests/fixtures/bgzf-biopython-generated.bgz");
    let mut expected: Vec<u8> = (0..140_000)
        .map(|index| u8::try_from((index * 17) % 256).unwrap())
        .collect();
    expected.extend_from_slice(b"\0\0>\r\n");
    assert_eq!(decoded(bytes).unwrap(), expected);
    assert_eq!(independent_crc(&expected), 3_774_126_258);
}

#[test]
fn exhaustive_single_byte_mutations_never_publish_unverified_data_or_panic() {
    let original = with_marker(stored(b"ACGT"));
    let mut mutations = 0;
    for index in 0..original.len() {
        for value in 0..=255_u8 {
            if value == original[index] {
                continue;
            }
            let mut changed = original.clone();
            changed[index] = value;
            if let Ok(output) = decoded(&changed) {
                assert_eq!(output, b"ACGT", "byte {index} changed to {value}");
            }
            mutations += 1;
        }
    }
    assert_eq!(mutations, 16_065);
}

struct FaultReader {
    bytes: Cursor<Vec<u8>>,
    fail_at: usize,
}
impl Read for FaultReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let position = usize::try_from(self.bytes.position()).unwrap();
        if position == self.fail_at {
            return Err(io::ErrorKind::PermissionDenied.into());
        }
        let amount = output.len().min(self.fail_at - position);
        self.bytes.read(&mut output[..amount])
    }
}

#[test]
fn underlying_read_failures_preserve_typed_cause_and_precise_compressed_offset() {
    let bytes = with_marker(stored(b"abc"));
    for fail_at in [0, 1, 11, 12, 17, 18, 20, 30] {
        let input = FaultReader {
            bytes: Cursor::new(bytes.clone()),
            fail_at,
        };
        let mut reader = BgzfReader::new(input);
        let error = reader.read_to_end(&mut Vec::new()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::PermissionDenied);
        let typed = error
            .get_ref()
            .unwrap()
            .downcast_ref::<BgzfError>()
            .unwrap();
        assert_eq!(typed.code(), "bgzf_io_error");
        assert!(
            matches!(typed, BgzfError::Io { offset, source, .. } if *offset == fail_at as u64 && source.kind() == io::ErrorKind::PermissionDenied)
        );
        assert!(std::error::Error::source(typed).is_some());
    }
}

#[test]
fn all_zero_uncompressed_bytes_are_preserved_as_data() {
    for length in [0, 1, 65_535, 65_536, 65_537, 131_073] {
        let input = vec![0; length];
        assert_eq!(
            decoded(&compressed(&input, BgzfWriteOptions::default())).unwrap(),
            input
        );
    }
}

#[derive(Debug)]
struct ZeroWriter;
impl Write for ZeroWriter {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Ok(0)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn zero_length_sink_writes_are_reported_instead_of_silently_finalizing() {
    let error = BgzfWriter::new(ZeroWriter).finish().unwrap_err();
    assert!(
        matches!(error, BgzfError::Io { source, offset: 0, .. } if source.kind() == io::ErrorKind::WriteZero)
    );
}

#[test]
fn legacy_eof_seeks_are_explicit_and_within_eof_offsets_are_refused() {
    let bytes = stored(b"abc");
    let target = BgzfVirtualOffset::new(bytes.len() as u64, 0).unwrap();
    let mut strict = BgzfReader::new(Cursor::new(bytes.clone()));
    assert!(matches!(
        strict.seek_virtual(target),
        Err(BgzfError::MissingEof { .. })
    ));
    let mut legacy = BgzfReader::with_options(
        Cursor::new(bytes.clone()),
        BgzfReadOptions {
            eof_policy: BgzfEofPolicy::AllowMissing,
        },
    );
    assert_eq!(legacy.seek_virtual(target).unwrap(), target);
    assert_eq!(legacy.eof_status(), BgzfEofStatus::Missing);
    assert!(matches!(
        legacy.seek_virtual(BgzfVirtualOffset::new(bytes.len() as u64, 1).unwrap()),
        Err(BgzfError::InvalidSeek { .. })
    ));
}
