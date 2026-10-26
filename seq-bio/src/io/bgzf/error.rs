use std::fmt;
use std::io;

use super::BgzfVirtualOffset;

/// A BGZF operation failed with original compressed positions retained.
#[derive(Debug)]
#[non_exhaustive]
pub enum BgzfError {
    /// An underlying stream operation failed.
    Io {
        /// Operation that failed.
        operation: &'static str,
        /// Compressed address at the operation.
        offset: u64,
        /// Original I/O failure.
        source: io::Error,
    },
    /// A member stopped before its declared component was complete.
    Truncated {
        /// Compressed member start.
        offset: u64,
        /// `header`, `extra fields`, `payload/trailer` or `terminal window`.
        part: &'static str,
        /// Required component bytes.
        expected: usize,
        /// Bytes obtained before physical EOF.
        actual: usize,
    },
    /// Magic, method or flags differ from the declared BGZF profile.
    InvalidHeader {
        /// Compressed member start.
        offset: u64,
        /// Observed ID1, ID2, CM and FLG bytes.
        observed: [u8; 4],
    },
    /// Required extra subfields are missing, duplicated or malformed.
    InvalidExtra {
        /// Compressed member start.
        offset: u64,
        /// Stable diagnostic explanation.
        message: &'static str,
    },
    /// Declared compressed size cannot contain the header and trailer.
    InvalidBlockSize {
        /// Compressed member start.
        offset: u64,
        /// Observed/required compressed size.
        compressed_size: usize,
        /// Complete header size including extra fields.
        header_size: usize,
    },
    /// Raw DEFLATE was invalid, incomplete or had unused payload bytes.
    Deflate {
        /// Compressed member start.
        offset: u64,
        /// Original backend error or checked end/consumption failure.
        message: String,
    },
    /// Declared or expanded uncompressed size exceeds 65,536 bytes.
    UncompressedLimit {
        /// Compressed member start.
        offset: u64,
        /// Declared or observed uncompressed byte count.
        observed: u64,
    },
    /// Trailer ISIZE differs from the actual decoded byte count.
    SizeMismatch {
        /// Compressed member start.
        offset: u64,
        /// Trailer ISIZE.
        expected: u32,
        /// Actual raw-DEFLATE output bytes.
        actual: u64,
    },
    /// Trailer CRC32 differs from the decoded bytes.
    CrcMismatch {
        /// Compressed member start.
        offset: u64,
        /// Trailer CRC32.
        expected: u32,
        /// Computed CRC32.
        actual: u32,
    },
    /// Strict policy requires a validated final empty member.
    MissingEof {
        /// Physical EOF address.
        offset: u64,
    },
    /// A virtual-position component exceeds its allocated width.
    InvalidVirtualOffset {
        /// Rejected compressed block address.
        compressed: u64,
        /// Rejected within-block data offset.
        uncompressed: u32,
    },
    /// Compressed address arithmetic exceeds the 48-bit position range.
    AddressOverflow,
    /// A seek target is beyond the source or its decoded data block.
    InvalidSeek {
        /// Requested virtual position.
        target: BgzfVirtualOffset,
        /// Complete compressed source length.
        source_length: u64,
        /// Decoded target block length, when available.
        block_length: Option<usize>,
    },
    /// Compression level or candidate block size is invalid.
    InvalidOptions {
        /// Rejected level.
        compression_level: u32,
        /// Rejected candidate data block size.
        block_size: usize,
    },
    /// A previous operation failed; reader recovery requires a successful seek.
    FailedState,
}

impl BgzfError {
    /// Returns a stable machine-readable category independent of display wording.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Io { .. } => "bgzf_io_error",
            Self::Truncated { .. } => "truncated_bgzf",
            Self::InvalidHeader { .. } => "invalid_bgzf_header",
            Self::InvalidExtra { .. } => "invalid_bgzf_extra",
            Self::InvalidBlockSize { .. } => "invalid_bgzf_block_size",
            Self::Deflate { .. } => "invalid_bgzf_deflate",
            Self::UncompressedLimit { .. } => "bgzf_uncompressed_limit",
            Self::SizeMismatch { .. } => "bgzf_size_mismatch",
            Self::CrcMismatch { .. } => "bgzf_crc_mismatch",
            Self::MissingEof { .. } => "missing_bgzf_eof",
            Self::InvalidVirtualOffset { .. } => "invalid_bgzf_virtual_offset",
            Self::AddressOverflow => "bgzf_address_overflow",
            Self::InvalidSeek { .. } => "invalid_bgzf_seek",
            Self::InvalidOptions { .. } => "invalid_bgzf_options",
            Self::FailedState => "bgzf_failed_state",
        }
    }
}

impl fmt::Display for BgzfError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io {
                operation,
                offset,
                source,
            } => write!(
                formatter,
                "BGZF {operation} failed at compressed byte {offset}: {source}"
            ),
            Self::Truncated {
                offset,
                part,
                expected,
                actual,
            } => write!(
                formatter,
                "BGZF block at {offset} has truncated {part}: expected {expected} bytes, obtained {actual}"
            ),
            Self::InvalidHeader { offset, observed } => write!(
                formatter,
                "BGZF block at {offset} has unsupported magic/method/flags {observed:?}"
            ),
            Self::InvalidExtra { offset, message } => write!(
                formatter,
                "BGZF block at {offset} has invalid extra fields: {message}"
            ),
            Self::InvalidBlockSize {
                offset,
                compressed_size,
                header_size,
            } => write!(
                formatter,
                "BGZF block at {offset} has compressed size {compressed_size} incompatible with header size {header_size}"
            ),
            Self::Deflate { offset, message } => write!(
                formatter,
                "BGZF raw DEFLATE failed at block {offset}: {message}"
            ),
            Self::UncompressedLimit { offset, observed } => write!(
                formatter,
                "BGZF block at {offset} expands/declares {observed} bytes; maximum is 65536"
            ),
            Self::SizeMismatch {
                offset,
                expected,
                actual,
            } => write!(
                formatter,
                "BGZF block at {offset} declares ISIZE {expected}, decoded {actual}"
            ),
            Self::CrcMismatch {
                offset,
                expected,
                actual,
            } => write!(
                formatter,
                "BGZF block at {offset} declares CRC32 {expected:08x}, computed {actual:08x}"
            ),
            Self::MissingEof { offset } => write!(
                formatter,
                "BGZF physical EOF at {offset} lacks a validated final empty member"
            ),
            Self::InvalidVirtualOffset {
                compressed,
                uncompressed,
            } => write!(
                formatter,
                "BGZF virtual offset components exceed 48+16 bits: compressed={compressed}, uncompressed={uncompressed}"
            ),
            Self::AddressOverflow => {
                formatter.write_str("BGZF compressed address exceeds its 48-bit position range")
            }
            Self::InvalidSeek {
                target,
                source_length,
                block_length,
            } => write!(
                formatter,
                "BGZF virtual target {} is invalid for compressed length {source_length}, decoded block length {block_length:?}",
                target.raw()
            ),
            Self::InvalidOptions {
                compression_level,
                block_size,
            } => write!(
                formatter,
                "BGZF requires compression level 0..=9 and block size 1..=65536; received {compression_level}, {block_size}"
            ),
            Self::FailedState => formatter.write_str("BGZF stream is in a failed state"),
        }
    }
}

impl std::error::Error for BgzfError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}
