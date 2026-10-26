//! Checked BGZF byte streams, explicit finalization and virtual positions.
//!
//! BGZF is a block-framed gzip profile, not a biological record format. Readers
//! validate each whole block before exposing its bytes. Call `finish` to verify
//! the remaining suffix and EOF policy; dropping a handle does not verify it.

use std::fmt;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use flate2::{Compress, Compression, Decompress, FlushCompress, FlushDecompress, Status};

use crate::publication::{PublishError, publish_new};

#[path = "bgzf/error.rs"]
mod error;
pub use error::BgzfError;

/// Maximum compressed member size, including the gzip header and trailer.
pub const BGZF_MAX_BLOCK_SIZE: usize = 65_536;
/// Maximum decompressed data size of one BGZF member.
pub const BGZF_MAX_UNCOMPRESSED_SIZE: usize = 65_536;
/// Largest compressed block address representable by a virtual offset.
pub const BGZF_MAX_COMPRESSED_OFFSET: u64 = (1_u64 << 48) - 1;
/// Canonical empty EOF member recommended by the SAM/BAM specification.
pub const BGZF_EOF_MARKER: [u8; 28] = [
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

const HEADER_SIZE: usize = 18;
const TRAILER_SIZE: usize = 8;

/// A comparable BGZF position, encoded as a 48-bit address and 16-bit data offset.
///
/// No addition/subtraction is defined. A value does not itself prove that its
/// compressed component is a block boundary in a particular source.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct BgzfVirtualOffset(u64);

impl BgzfVirtualOffset {
    /// Checks both components before constructing the standard packed position.
    ///
    /// # Errors
    /// Returns `InvalidVirtualOffset` if either component exceeds its bit width.
    pub fn new(compressed: u64, uncompressed: u32) -> Result<Self, BgzfError> {
        if compressed > BGZF_MAX_COMPRESSED_OFFSET || uncompressed > u32::from(u16::MAX) {
            return Err(BgzfError::InvalidVirtualOffset {
                compressed,
                uncompressed,
            });
        }
        Ok(Self((compressed << 16) | u64::from(uncompressed)))
    }

    /// Interprets a packed `u64`; physical boundary validity remains unchecked.
    #[must_use]
    pub const fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    /// Returns the standard packed representation.
    #[must_use]
    pub const fn raw(self) -> u64 {
        self.0
    }

    /// Returns the compressed block address.
    #[must_use]
    pub const fn compressed_offset(self) -> u64 {
        self.0 >> 16
    }

    /// Returns the within-block uncompressed byte offset.
    #[must_use]
    pub const fn uncompressed_offset(self) -> u16 {
        let bytes = self.0.to_le_bytes();
        u16::from_le_bytes([bytes[0], bytes[1]])
    }
}

impl From<BgzfVirtualOffset> for u64 {
    fn from(value: BgzfVirtualOffset) -> Self {
        value.raw()
    }
}

/// Whether physical EOF must immediately follow a validated empty member.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BgzfEofPolicy {
    /// Refuse a missing terminal empty member, including an empty physical file.
    #[default]
    Require,
    /// Accept a missing marker for legacy streams and report `Missing` status.
    AllowMissing,
}

/// Observed terminal-member status, not proof of unread-prefix validity.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BgzfEofStatus {
    /// Physical EOF has not yet been observed on this traversal.
    #[default]
    Unchecked,
    /// Physical EOF immediately follows a validated empty member.
    Present,
    /// Physical EOF lacks a terminal empty member under explicit legacy permission.
    Missing,
}

/// Reader policies; individual block buffers always obey the format bounds.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BgzfReadOptions {
    /// Strict terminal-empty-member checking by default.
    pub eof_policy: BgzfEofPolicy,
}

/// Structural information about one completely validated block.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BgzfBlockInfo {
    /// Compressed byte address of the member start.
    pub compressed_offset: u64,
    /// Complete compressed member size, including framing.
    pub compressed_size: usize,
    /// Completely validated uncompressed byte count.
    pub uncompressed_size: usize,
    /// Whether the complete encoded member equals the canonical 28-byte marker.
    pub canonical_eof_marker: bool,
}

pub(crate) struct Block {
    pub(crate) info: BgzfBlockInfo,
    pub(crate) data: Vec<u8>,
}

/// A bounded streaming BGZF decoder that retains one validated data block.
///
/// New readers begin with compressed address zero relative to the supplied
/// stream start. For seeking, the underlying seek coordinate system must agree;
/// use `with_origin` for an explicitly prepositioned seekable source.
#[derive(Debug)]
pub struct BgzfReader<R> {
    inner: R,
    options: BgzfReadOptions,
    data: Vec<u8>,
    cursor: usize,
    current: Option<BgzfBlockInfo>,
    next_block: u64,
    last_empty: bool,
    eof_status: BgzfEofStatus,
    physical_eof: bool,
    failed: bool,
}

impl<R: Read> BgzfReader<R> {
    /// Creates a lazy reader with strict EOF policy; no bytes are read yet.
    #[must_use]
    pub fn new(inner: R) -> Self {
        Self::with_options(inner, BgzfReadOptions::default())
    }

    /// Creates a lazy reader with an explicit EOF policy and origin zero.
    #[must_use]
    pub fn with_options(inner: R, options: BgzfReadOptions) -> Self {
        Self {
            inner,
            options,
            data: Vec::new(),
            cursor: 0,
            current: None,
            next_block: 0,
            last_empty: false,
            eof_status: BgzfEofStatus::Unchecked,
            physical_eof: false,
            failed: false,
        }
    }

    /// Creates a reader at a caller-established compressed source address.
    ///
    /// The source must already be positioned at `origin`. This method declares
    /// that position; it neither seeks nor validates an unread prefix.
    ///
    /// # Errors
    /// Returns `InvalidVirtualOffset` for an origin outside 48 bits.
    pub fn with_origin(inner: R, origin: u64, options: BgzfReadOptions) -> Result<Self, BgzfError> {
        BgzfVirtualOffset::new(origin, 0)?;
        let mut result = Self::with_options(inner, options);
        result.next_block = origin;
        Ok(result)
    }

    /// Returns the current canonical virtual position without reading more data.
    ///
    /// Consuming a whole block always reports the next block at offset zero,
    /// including after a 65,536-byte block and after a terminal empty block.
    ///
    /// # Errors
    /// Returns `FailedState` after a failed decode/seek.
    pub fn virtual_position(&self) -> Result<BgzfVirtualOffset, BgzfError> {
        self.check_state()?;
        if self.cursor == self.data.len() {
            BgzfVirtualOffset::new(self.next_block, 0)
        } else {
            let start = self
                .current
                .map_or(self.next_block, |info| info.compressed_offset);
            let within = u32::try_from(self.cursor).map_err(|_| BgzfError::AddressOverflow)?;
            BgzfVirtualOffset::new(start, within)
        }
    }

    /// Returns information about the most recently validated member.
    #[must_use]
    pub const fn current_block(&self) -> Option<&BgzfBlockInfo> {
        self.current.as_ref()
    }

    /// Returns the EOF status established on the current traversal.
    #[must_use]
    pub const fn eof_status(&self) -> BgzfEofStatus {
        self.eof_status
    }

    /// Borrows the underlying stream without allowing untracked mutation.
    #[must_use]
    pub const fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns the underlying stream without draining or validating unread data.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Copies bytes from one validated data block, returning zero only at EOF or
    /// for an empty output buffer. Embedded empty members are skipped.
    ///
    /// # Errors
    /// Returns typed framing, decompression, integrity, EOF, address or I/O errors.
    /// No bytes from a failed block are exposed. Earlier reads may have returned
    /// a valid prefix. A subsequent read refuses a failed state.
    pub fn read_bytes(&mut self, output: &mut [u8]) -> Result<usize, BgzfError> {
        self.check_state()?;
        if output.is_empty() {
            return Ok(0);
        }
        self.ensure_data()?;
        let amount = output.len().min(self.data.len() - self.cursor);
        output[..amount].copy_from_slice(&self.data[self.cursor..self.cursor + amount]);
        self.cursor += amount;
        Ok(amount)
    }

    /// Drains and validates the remaining suffix, returning its EOF status.
    ///
    /// After a seek, unread earlier blocks are not checked by this call.
    ///
    /// # Errors
    /// Returns the same typed errors as [`Self::read_bytes`].
    pub fn finish(&mut self) -> Result<BgzfEofStatus, BgzfError> {
        loop {
            self.ensure_data()?;
            if self.physical_eof {
                return Ok(self.eof_status);
            }
            self.cursor = self.data.len();
        }
    }

    fn check_state(&self) -> Result<(), BgzfError> {
        if self.failed {
            Err(BgzfError::FailedState)
        } else {
            Ok(())
        }
    }

    fn ensure_data(&mut self) -> Result<(), BgzfError> {
        self.check_state()?;
        while self.cursor == self.data.len() && !self.physical_eof {
            match read_block(&mut self.inner, &mut self.next_block) {
                Ok(Some(block)) => {
                    self.last_empty = block.data.is_empty();
                    self.data = block.data;
                    self.cursor = 0;
                    self.current = Some(block.info);
                }
                Ok(None) => {
                    self.physical_eof = true;
                    self.eof_status = if self.last_empty {
                        BgzfEofStatus::Present
                    } else {
                        BgzfEofStatus::Missing
                    };
                    if !self.last_empty && self.options.eof_policy == BgzfEofPolicy::Require {
                        self.failed = true;
                        return Err(BgzfError::MissingEof {
                            offset: self.next_block,
                        });
                    }
                }
                Err(error) => {
                    self.failed = true;
                    return Err(error);
                }
            }
        }
        Ok(())
    }
}

impl<R: Read + Seek> BgzfReader<R> {
    /// Seeks to and validates a BGZF block and its within-block offset.
    ///
    /// A decodable target is not authentication of source-wide block boundaries.
    /// Physical EOF is accepted at offset zero if its terminal member satisfies
    /// the selected EOF policy. Exact short-block ends are normalized. Successful
    /// seeks recover a failed reader; failed seeks leave it failed.
    ///
    /// # Errors
    /// Returns typed seek, framing, integrity, EOF, address or underlying I/O errors.
    pub fn seek_virtual(
        &mut self,
        target: BgzfVirtualOffset,
    ) -> Result<BgzfVirtualOffset, BgzfError> {
        self.failed = true;
        let length = self
            .inner
            .seek(SeekFrom::End(0))
            .map_err(|source| BgzfError::Io {
                operation: "seek source end",
                offset: target.compressed_offset(),
                source,
            })?;
        let compressed = target.compressed_offset();
        let within = usize::from(target.uncompressed_offset());
        if compressed > length || (compressed == length && within != 0) {
            return Err(BgzfError::InvalidSeek {
                target,
                source_length: length,
                block_length: None,
            });
        }
        if compressed == length {
            let present = terminal_empty(&mut self.inner, length)?;
            if !present && self.options.eof_policy == BgzfEofPolicy::Require {
                return Err(BgzfError::MissingEof { offset: length });
            }
            self.inner
                .seek(SeekFrom::Start(length))
                .map_err(|source| BgzfError::Io {
                    operation: "seek physical EOF",
                    offset: length,
                    source,
                })?;
            self.data.clear();
            self.cursor = 0;
            self.current = None;
            self.next_block = length;
            self.last_empty = present;
            self.physical_eof = true;
            self.eof_status = if present {
                BgzfEofStatus::Present
            } else {
                BgzfEofStatus::Missing
            };
        } else {
            self.inner
                .seek(SeekFrom::Start(compressed))
                .map_err(|source| BgzfError::Io {
                    operation: "seek target block",
                    offset: compressed,
                    source,
                })?;
            let mut next = compressed;
            let block = read_block(&mut self.inner, &mut next)?.ok_or(BgzfError::InvalidSeek {
                target,
                source_length: length,
                block_length: None,
            })?;
            if within > block.data.len() {
                return Err(BgzfError::InvalidSeek {
                    target,
                    source_length: length,
                    block_length: Some(block.data.len()),
                });
            }
            self.last_empty = block.data.is_empty();
            self.data = block.data;
            self.cursor = within;
            self.current = Some(block.info);
            self.next_block = next;
            self.physical_eof = false;
            self.eof_status = BgzfEofStatus::Unchecked;
        }
        self.failed = false;
        self.virtual_position()
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        self.read_bytes(output).map_err(io_error)
    }
}

impl<R: Read> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.ensure_data().map_err(io_error)?;
        Ok(&self.data[self.cursor..])
    }

    fn consume(&mut self, amount: usize) {
        self.cursor += amount.min(self.data.len() - self.cursor);
    }
}

fn io_error(error: BgzfError) -> io::Error {
    let kind = match &error {
        BgzfError::Io { source, .. } => source.kind(),
        BgzfError::Truncated { .. } | BgzfError::MissingEof { .. } => io::ErrorKind::UnexpectedEof,
        BgzfError::FailedState => io::ErrorKind::Other,
        BgzfError::InvalidOptions { .. } | BgzfError::InvalidVirtualOffset { .. } => {
            io::ErrorKind::InvalidInput
        }
        _ => io::ErrorKind::InvalidData,
    };
    io::Error::new(kind, error)
}

fn read_exact<R: Read>(
    reader: &mut R,
    output: &mut [u8],
    position: &mut u64,
    block: u64,
    part: &'static str,
) -> Result<(), BgzfError> {
    let mut actual = 0;
    while actual < output.len() {
        match reader.read(&mut output[actual..]) {
            Ok(0) => {
                return Err(BgzfError::Truncated {
                    offset: block,
                    part,
                    expected: output.len(),
                    actual,
                });
            }
            Ok(amount) => {
                actual += amount;
                *position = next_address(*position, amount)?;
            }
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => {
                return Err(BgzfError::Io {
                    operation: "read block",
                    offset: *position,
                    source,
                });
            }
        }
    }
    Ok(())
}

pub(crate) fn read_block<R: Read>(
    reader: &mut R,
    position: &mut u64,
) -> Result<Option<Block>, BgzfError> {
    let start = *position;
    let mut fixed = [0_u8; 12];
    loop {
        match reader.read(&mut fixed[..1]) {
            Ok(0) => return Ok(None),
            Ok(_) => {
                *position = next_address(*position, 1)?;
                break;
            }
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => {
                return Err(BgzfError::Io {
                    operation: "read block start",
                    offset: start,
                    source,
                });
            }
        }
    }
    read_exact(reader, &mut fixed[1..], position, start, "header")?;
    if fixed[..4] != [31, 139, 8, 4] {
        return Err(BgzfError::InvalidHeader {
            offset: start,
            observed: [fixed[0], fixed[1], fixed[2], fixed[3]],
        });
    }
    let xlen = usize::from(u16::from_le_bytes([fixed[10], fixed[11]]));
    let header_length = 12 + xlen;
    if header_length > BGZF_MAX_BLOCK_SIZE - TRAILER_SIZE {
        return Err(BgzfError::InvalidBlockSize {
            offset: start,
            compressed_size: header_length + TRAILER_SIZE,
            header_size: header_length,
        });
    }
    let mut raw = Vec::with_capacity(header_length);
    raw.extend_from_slice(&fixed);
    raw.resize(header_length, 0);
    read_exact(reader, &mut raw[12..], position, start, "extra fields")?;
    let size = block_size(&raw[12..], start)?;
    if size < header_length + TRAILER_SIZE {
        return Err(BgzfError::InvalidBlockSize {
            offset: start,
            compressed_size: size,
            header_size: header_length,
        });
    }
    next_address(start, size)?;
    raw.resize(size, 0);
    read_exact(
        reader,
        &mut raw[header_length..],
        position,
        start,
        "payload/trailer",
    )?;
    decode_block(&raw, header_length, start).map(Some)
}

fn decode_block(raw: &[u8], header_length: usize, start: u64) -> Result<Block, BgzfError> {
    let size = raw.len();
    let trailer = size - TRAILER_SIZE;
    let expected_crc = le_u32(&raw[trailer..trailer + 4]);
    let expected_size = le_u32(&raw[trailer + 4..]);
    if u64::from(expected_size) > BGZF_MAX_UNCOMPRESSED_SIZE as u64 {
        return Err(BgzfError::UncompressedLimit {
            offset: start,
            observed: u64::from(expected_size),
        });
    }
    let payload = &raw[header_length..trailer];
    let mut data = vec![0_u8; BGZF_MAX_UNCOMPRESSED_SIZE + 1];
    let mut decoder = Decompress::new(false);
    let status = decoder
        .decompress(payload, &mut data, FlushDecompress::Finish)
        .map_err(|error| BgzfError::Deflate {
            offset: start,
            message: error.to_string(),
        })?;
    let actual_size = decoder.total_out();
    if actual_size > BGZF_MAX_UNCOMPRESSED_SIZE as u64 {
        return Err(BgzfError::UncompressedLimit {
            offset: start,
            observed: actual_size,
        });
    }
    if status != Status::StreamEnd || decoder.total_in() != payload.len() as u64 {
        return Err(BgzfError::Deflate {
            offset: start,
            message: "raw DEFLATE must finish and consume the complete payload".to_owned(),
        });
    }
    if actual_size != u64::from(expected_size) {
        return Err(BgzfError::SizeMismatch {
            offset: start,
            expected: expected_size,
            actual: actual_size,
        });
    }
    data.truncate(usize::try_from(actual_size).map_err(|_| BgzfError::AddressOverflow)?);
    let actual_crc = crc32fast::hash(&data);
    if actual_crc != expected_crc {
        return Err(BgzfError::CrcMismatch {
            offset: start,
            expected: expected_crc,
            actual: actual_crc,
        });
    }
    Ok(Block {
        info: BgzfBlockInfo {
            compressed_offset: start,
            compressed_size: size,
            uncompressed_size: data.len(),
            canonical_eof_marker: raw == BGZF_EOF_MARKER,
        },
        data,
    })
}

fn block_size(extra: &[u8], offset: u64) -> Result<usize, BgzfError> {
    let mut position = 0;
    let mut found = None;
    while position < extra.len() {
        if extra.len() - position < 4 {
            return Err(BgzfError::InvalidExtra {
                offset,
                message: "truncated extra subfield header",
            });
        }
        let length = usize::from(u16::from_le_bytes([
            extra[position + 2],
            extra[position + 3],
        ]));
        let end = position + 4 + length;
        if end > extra.len() {
            return Err(BgzfError::InvalidExtra {
                offset,
                message: "extra subfield exceeds XLEN",
            });
        }
        if &extra[position..position + 2] == b"BC" {
            if found.is_some() || length != 2 {
                return Err(BgzfError::InvalidExtra {
                    offset,
                    message: "BC must occur exactly once with SLEN=2",
                });
            }
            found = Some(
                usize::from(u16::from_le_bytes([
                    extra[position + 4],
                    extra[position + 5],
                ])) + 1,
            );
        }
        position = end;
    }
    found.ok_or(BgzfError::InvalidExtra {
        offset,
        message: "missing BC subfield",
    })
}

fn le_u32(bytes: &[u8]) -> u32 {
    u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

fn next_address(position: u64, amount: usize) -> Result<u64, BgzfError> {
    position
        .checked_add(u64::try_from(amount).map_err(|_| BgzfError::AddressOverflow)?)
        .filter(|value| *value <= BGZF_MAX_COMPRESSED_OFFSET)
        .ok_or(BgzfError::AddressOverflow)
}

fn terminal_empty<R: Read + Seek>(reader: &mut R, length: u64) -> Result<bool, BgzfError> {
    let count = usize::try_from(length.min(BGZF_MAX_BLOCK_SIZE as u64))
        .map_err(|_| BgzfError::AddressOverflow)?;
    let start = length - count as u64;
    reader
        .seek(SeekFrom::Start(start))
        .map_err(|source| BgzfError::Io {
            operation: "seek terminal member",
            offset: start,
            source,
        })?;
    let mut tail = vec![0; count];
    let mut position = start;
    read_exact(reader, &mut tail, &mut position, start, "terminal window")?;
    for index in 0..tail.len().saturating_sub(11) {
        if tail.get(index..index + 4) != Some(&[31, 139, 8, 4]) {
            continue;
        }
        let xlen = usize::from(u16::from_le_bytes([tail[index + 10], tail[index + 11]]));
        let Some(extra) = tail.get(index + 12..index + 12 + xlen) else {
            continue;
        };
        let Ok(size) = block_size(extra, start + index as u64) else {
            continue;
        };
        if size != tail.len() - index {
            continue;
        }
        let mut slice = &tail[index..];
        let mut address = start + index as u64;
        if let Ok(Some(block)) = read_block(&mut slice, &mut address) {
            if block.data.is_empty() && address == length {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Compression options checked before any writer output is produced.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BgzfWriteOptions {
    /// Compression level in 0..=9; default 6 uses the existing Rust backend.
    pub compression_level: u32,
    /// Maximum candidate input block size in 1..=65,536; default 65,280.
    pub block_size: usize,
}

impl Default for BgzfWriteOptions {
    fn default() -> Self {
        Self {
            compression_level: 6,
            block_size: 65_280,
        }
    }
}

impl BgzfWriteOptions {
    /// Checks the compression level and candidate block size.
    ///
    /// # Errors
    /// Returns `InvalidOptions` for unsupported values.
    pub fn validate(self) -> Result<(), BgzfError> {
        if self.compression_level > 9
            || !(1..=BGZF_MAX_UNCOMPRESSED_SIZE).contains(&self.block_size)
        {
            return Err(BgzfError::InvalidOptions {
                compression_level: self.compression_level,
                block_size: self.block_size,
            });
        }
        Ok(())
    }
}

/// A buffered byte writer that finalizes only through explicit `finish`.
///
/// Dropping it performs no hidden I/O and may leave an incomplete output. Named
/// file publication must separately use the existing create-only transaction.
#[derive(Debug)]
pub struct BgzfWriter<W> {
    inner: W,
    options: BgzfWriteOptions,
    data: Vec<u8>,
    position: u64,
    failed: bool,
}

/// A create-only BGZF publication failure retaining the application's own error.
#[derive(Debug)]
#[non_exhaustive]
pub enum BgzfOutputError<E> {
    /// Compression, framing or explicit finalization failed.
    Compression(BgzfError),
    /// The biological/application writer callback failed.
    Write(E),
    /// Temporary output, synchronization or publication failed.
    Io(io::Error),
    /// Any existing file, directory or symlink occupies the destination.
    DestinationExists(PathBuf),
    /// No same-directory temporary filename could be reserved.
    TemporaryUnavailable,
}

impl<E: fmt::Display> fmt::Display for BgzfOutputError<E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Compression(error) => error.fmt(formatter),
            Self::Write(error) => write!(formatter, "BGZF data writer failed: {error}"),
            Self::Io(error) => error.fmt(formatter),
            Self::DestinationExists(path) => write!(
                formatter,
                "destination already exists; output is create-only: {}",
                path.display()
            ),
            Self::TemporaryUnavailable => {
                formatter.write_str("could not reserve temporary BGZF output")
            }
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for BgzfOutputError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Compression(error) => Some(error),
            Self::Write(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::DestinationExists(_) | Self::TemporaryUnavailable => None,
        }
    }
}

/// Explicitly finalizes a new BGZF file before create-only publication.
///
/// Existing FASTA/FASTQ writers can write directly to the supplied `Write` sink.
/// The destination directory must be trusted. No existing destination is replaced.
/// Failed callbacks/finalization leave no named partial output or temporary file.
///
/// # Errors
/// Returns distinct callback, compression, existing-path and publication errors.
pub fn write_bgzf_new<E>(
    path: impl AsRef<Path>,
    options: BgzfWriteOptions,
    write: impl FnOnce(&mut dyn Write) -> Result<(), E>,
) -> Result<(), BgzfOutputError<E>> {
    options.validate().map_err(BgzfOutputError::Compression)?;
    publish_new(path.as_ref(), |output| {
        let mut writer =
            BgzfWriter::with_options(output, options).map_err(BgzfOutputError::Compression)?;
        write(&mut writer).map_err(BgzfOutputError::Write)?;
        writer.finish().map_err(BgzfOutputError::Compression)?;
        Ok(())
    })
    .map_err(|error| match error {
        PublishError::Write(error) => error,
        PublishError::Io(error) => BgzfOutputError::Io(error),
        PublishError::DestinationExists(path) => BgzfOutputError::DestinationExists(path),
        PublishError::TemporaryUnavailable => BgzfOutputError::TemporaryUnavailable,
    })
}

impl<W: Write> BgzfWriter<W> {
    /// Creates an origin-zero writer with default options and no initial I/O.
    #[must_use]
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            options: BgzfWriteOptions::default(),
            data: Vec::new(),
            position: 0,
            failed: false,
        }
    }

    /// Validates explicit options before buffering or writing any data.
    ///
    /// # Errors
    /// Returns `InvalidOptions` for unsupported values.
    pub fn with_options(inner: W, options: BgzfWriteOptions) -> Result<Self, BgzfError> {
        Self::with_origin(inner, 0, options)
    }

    /// Declares a caller-established output address and validates options.
    ///
    /// No seek or append validation occurs; the sink must already be positioned.
    ///
    /// # Errors
    /// Returns `InvalidOptions` or `InvalidVirtualOffset`.
    pub fn with_origin(
        inner: W,
        origin: u64,
        options: BgzfWriteOptions,
    ) -> Result<Self, BgzfError> {
        options.validate()?;
        BgzfVirtualOffset::new(origin, 0)?;
        Ok(Self {
            inner,
            options,
            data: Vec::new(),
            position: origin,
            failed: false,
        })
    }

    /// Borrows the underlying output without permitting untracked writes.
    #[must_use]
    pub const fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Buffers up to one candidate block, returning the number accepted.
    ///
    /// Full prior buffers are emitted before accepting new input. Accepted bytes
    /// can remain buffered until a later write, flush or finish. No implicit final
    /// EOF marker is emitted by this method.
    ///
    /// # Errors
    /// Returns failed-state, compression, address or underlying write errors.
    pub fn write_bytes(&mut self, input: &[u8]) -> Result<usize, BgzfError> {
        self.check_state()?;
        if input.is_empty() {
            return Ok(0);
        }
        if self.data.len() == self.options.block_size {
            self.emit_pending()?;
        }
        let amount = input.len().min(self.options.block_size - self.data.len());
        self.data.extend_from_slice(&input[..amount]);
        Ok(amount)
    }

    /// Emits pending blocks and flushes the sink without writing a terminal marker.
    ///
    /// # Errors
    /// Returns compression, address, failed-state or underlying output errors.
    pub fn flush_blocks(&mut self) -> Result<(), BgzfError> {
        self.check_state()?;
        self.emit_pending()?;
        if let Err(source) = self.inner.flush() {
            self.failed = true;
            return Err(BgzfError::Io {
                operation: "flush output",
                offset: self.position,
                source,
            });
        }
        Ok(())
    }

    /// Flushes pending data and returns a stable next-block position.
    ///
    /// Calling this can affect block boundaries; it avoids publishing an offset
    /// inside a pending incompressible candidate that later needs splitting.
    ///
    /// # Errors
    /// Returns the same errors as [`Self::flush_blocks`].
    pub fn virtual_position(&mut self) -> Result<BgzfVirtualOffset, BgzfError> {
        self.flush_blocks()?;
        BgzfVirtualOffset::new(self.position, 0)
    }

    /// Emits all data, the canonical empty marker, flushes and returns the sink.
    ///
    /// # Errors
    /// Returns compression, address, failed-state or underlying output errors.
    /// Failure can leave a partial output; no successful finalized sink is returned.
    pub fn finish(mut self) -> Result<W, BgzfError> {
        self.check_state()?;
        self.emit_pending()?;
        self.publish(&BGZF_EOF_MARKER)?;
        self.inner.flush().map_err(|source| BgzfError::Io {
            operation: "finalize output",
            offset: self.position,
            source,
        })?;
        Ok(self.inner)
    }

    fn check_state(&self) -> Result<(), BgzfError> {
        if self.failed {
            Err(BgzfError::FailedState)
        } else {
            Ok(())
        }
    }

    fn emit_pending(&mut self) -> Result<(), BgzfError> {
        let result = self.emit_pending_inner();
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn emit_pending_inner(&mut self) -> Result<(), BgzfError> {
        while !self.data.is_empty() {
            let mut amount = self.data.len();
            let encoded = loop {
                if let Some(encoded) = encode_block(
                    &self.data[..amount],
                    self.options.compression_level,
                    self.position,
                )? {
                    break encoded;
                }
                if amount == 1 {
                    return Err(BgzfError::Deflate {
                        offset: self.position,
                        message: "a one-byte input did not fit a BGZF block".to_owned(),
                    });
                }
                amount = (amount / 2).max(1);
            };
            self.publish(&encoded)?;
            self.data.drain(..amount);
        }
        Ok(())
    }

    fn publish(&mut self, encoded: &[u8]) -> Result<(), BgzfError> {
        let next = next_address(self.position, encoded.len())?;
        if let Err(source) = self.inner.write_all(encoded) {
            self.failed = true;
            return Err(BgzfError::Io {
                operation: "write block",
                offset: self.position,
                source,
            });
        }
        self.position = next;
        Ok(())
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        self.write_bytes(input).map_err(io_error)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_blocks().map_err(io_error)
    }
}

fn encode_block(data: &[u8], level: u32, offset: u64) -> Result<Option<Vec<u8>>, BgzfError> {
    let mut compressed = vec![0; BGZF_MAX_BLOCK_SIZE - HEADER_SIZE - TRAILER_SIZE];
    let mut compressor = Compress::new(Compression::new(level), false);
    let status = compressor
        .compress(data, &mut compressed, FlushCompress::Finish)
        .map_err(|error| BgzfError::Deflate {
            offset,
            message: error.to_string(),
        })?;
    if status != Status::StreamEnd {
        return Ok(None);
    }
    if compressor.total_in() != data.len() as u64 {
        return Err(BgzfError::Deflate {
            offset,
            message: "compressor did not consume the complete input".to_owned(),
        });
    }
    compressed
        .truncate(usize::try_from(compressor.total_out()).map_err(|_| BgzfError::AddressOverflow)?);
    let size = HEADER_SIZE + compressed.len() + TRAILER_SIZE;
    let bsize = u16::try_from(size - 1).map_err(|_| BgzfError::AddressOverflow)?;
    let mut encoded = Vec::with_capacity(size);
    encoded.extend_from_slice(&[31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0]);
    encoded.extend_from_slice(&bsize.to_le_bytes());
    encoded.extend_from_slice(&compressed);
    encoded.extend_from_slice(&crc32fast::hash(data).to_le_bytes());
    encoded.extend_from_slice(
        &u32::try_from(data.len())
            .map_err(|_| BgzfError::AddressOverflow)?
            .to_le_bytes(),
    );
    Ok(Some(encoded))
}

#[cfg(test)]
#[path = "bgzf/tests.rs"]
mod tests;
