//! Bounded biological file readers and explicit, checked conversions.
//!
//! The legacy `read_fasta`/`parse_fasta` exports retain their original behavior.
//! New streaming readers are [`fasta`](crate::io::fasta),
//! [`fastq`](crate::io::fastq), [`FastaReader`](crate::io::FastaReader)
//! and [`FastqReader`](crate::io::FastqReader).

use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;

pub use align::advanced::formats::*;

pub mod alignment;
pub mod annotated;
pub mod bgzf;
pub mod bgzf_index;
pub mod blast;
pub mod graph;
pub mod gzip;
pub mod index;
pub mod result_json;
mod sequences;
pub use sequences::{
    FastaReader, FastqReader, LossPolicy, fasta, fastq, one_record, write_record_fasta,
    write_record_fastq,
};

/// Stable categories shared by the new file import APIs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IoErrorKind {
    /// A filesystem, stream or decompression operation failed.
    Io,
    /// A required identifier or header is invalid.
    InvalidHeader,
    /// A residue is invalid for the selected molecule.
    InvalidSequence,
    /// A quality character or value is invalid.
    InvalidQuality,
    /// Declared and observed dimensions disagree.
    LengthMismatch,
    /// Required input was missing at the end of a stream.
    Truncated,
    /// Valid syntax or semantics is outside the supported profile.
    Unsupported,
    /// An explicit resource limit was exceeded.
    LimitExceeded,
    /// A structural invariant or grammar rule was violated.
    InvalidStructure,
    /// The requested conversion would discard information.
    DataLoss,
    /// A JSON schema/version or required field is invalid.
    InvalidSchema,
}

impl IoErrorKind {
    /// Returns the stable machine-readable category code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Io => "io_error",
            Self::InvalidHeader => "invalid_header",
            Self::InvalidSequence => "invalid_sequence",
            Self::InvalidQuality => "invalid_quality",
            Self::LengthMismatch => "length_mismatch",
            Self::Truncated => "truncated_input",
            Self::Unsupported => "unsupported_format_feature",
            Self::LimitExceeded => "limit_exceeded",
            Self::InvalidStructure => "invalid_structure",
            Self::DataLoss => "data_loss",
            Self::InvalidSchema => "invalid_schema",
        }
    }
}

/// A structured input error with optional one-based record/line context.
#[derive(Debug)]
pub struct BioIoError {
    kind: IoErrorKind,
    format: &'static str,
    message: String,
    record: Option<usize>,
    line: Option<usize>,
    source: Option<std::io::Error>,
}

impl BioIoError {
    /// Constructs a contextualizable file error.
    pub fn new(kind: IoErrorKind, format: &'static str, message: impl Into<String>) -> Self {
        Self {
            kind,
            format,
            message: message.into(),
            record: None,
            line: None,
            source: None,
        }
    }
    /// Attaches one-based source record and line numbers.
    #[must_use]
    pub const fn at(mut self, record: usize, line: usize) -> Self {
        self.record = Some(record);
        self.line = Some(line);
        self
    }
    /// Returns the typed category.
    #[must_use]
    pub const fn kind(&self) -> IoErrorKind {
        self.kind
    }
    /// Returns the stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.kind.code()
    }
    /// Returns the source format or processing layer.
    #[must_use]
    pub const fn format(&self) -> &'static str {
        self.format
    }
    /// Returns the one-based record ordinal, when known.
    #[must_use]
    pub const fn record(&self) -> Option<usize> {
        self.record
    }
    /// Returns the one-based line number, when known.
    #[must_use]
    pub const fn line(&self) -> Option<usize> {
        self.line
    }
    /// Returns the explanatory message; control flow should use `kind`/`code`.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
    pub(crate) fn in_format(mut self, format: &'static str) -> Self {
        self.format = format;
        self
    }
}

impl fmt::Display for BioIoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.format, self.message)?;
        if let (Some(record), Some(line)) = (self.record, self.line) {
            write!(f, " (record {record}, line {line})")?;
        }
        Ok(())
    }
}
impl std::error::Error for BioIoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_ref()
            .map(|error| error as &(dyn std::error::Error + 'static))
    }
}

/// Maximum input sizes, enforced before buffering a line or record.
#[derive(Clone, Copy, Debug)]
pub struct ReadLimits {
    /// Bytes per physical line including its line ending (default 8 MiB).
    pub max_line_bytes: usize,
    /// Symbols per biological record (default 100 million).
    pub max_record_symbols: usize,
    /// Total decompressed source bytes (default 1 GiB).
    pub max_total_bytes: usize,
}
impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_line_bytes: 8 * 1024 * 1024,
            max_record_symbols: 100_000_000,
            max_total_bytes: 1024 * 1024 * 1024,
        }
    }
}

/// Opens plain or gzip data by content signature, not filename extension.
///
/// Concatenated gzip members are decoded; integrity is checked as the entire
/// input is consumed. Stopping a reader early cannot validate an unread trailer.
///
/// # Errors
/// Returns an I/O error if the path cannot be opened or inspected. Later
/// decompression errors are returned by the consuming parser.
pub fn open_input(path: impl AsRef<Path>) -> Result<Box<dyn BufRead>, BioIoError> {
    let file = File::open(path).map_err(|e| io_error("input", e))?;
    decode_input(file)
}

/// Wraps a stream with plain/gzip detection, including concatenated members.
///
/// # Errors
/// Returns initial stream inspection errors. Later decompression errors surface
/// while the returned stream is read. Consume to EOF to verify gzip integrity.
pub fn decode_input(mut reader: impl Read + 'static) -> Result<Box<dyn BufRead>, BioIoError> {
    let mut signature = [0u8; 2];
    let mut count = 0;
    while count < signature.len() {
        match reader.read(&mut signature[count..]) {
            Ok(0) => break,
            Ok(length) => count += length,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(io_error("input", error)),
        }
    }
    let gzip = count == 2 && signature == [0x1f, 0x8b];
    let reader = BufReader::new(Cursor::new(signature[..count].to_vec()).chain(reader));
    if gzip {
        Ok(Box::new(BufReader::new(
            flate2::bufread::MultiGzDecoder::new(reader),
        )))
    } else {
        Ok(Box::new(reader))
    }
}

pub(crate) fn io_error(format: &'static str, error: std::io::Error) -> BioIoError {
    let mut failure = BioIoError::new(IoErrorKind::Io, format, error.to_string());
    failure.source = Some(error);
    failure
}

/// Raw-byte boundaries retained internally for validated sequence indexing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RecordSpan {
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) header_line: usize,
}

pub(crate) struct Lines<R> {
    reader: R,
    pub(crate) limits: ReadLimits,
    pub(crate) line: usize,
    total: usize,
    last_start: usize,
    format: &'static str,
}

impl<R: BufRead> Lines<R> {
    pub(crate) fn new(reader: R, limits: ReadLimits, format: &'static str) -> Self {
        Self {
            reader,
            limits,
            line: 0,
            total: 0,
            last_start: 0,
            format,
        }
    }

    pub(crate) fn byte_position(&self) -> u64 {
        self.total as u64
    }

    pub(crate) fn line_start(&self) -> u64 {
        self.last_start as u64
    }

    pub(crate) fn next_line(&mut self) -> Result<Option<String>, BioIoError> {
        let Some(mut bytes) = self.next_physical_bytes()? else {
            return Ok(None);
        };
        if bytes.last() == Some(&b'\n') {
            bytes.pop();
            if bytes.last() == Some(&b'\r') {
                bytes.pop();
            }
        }
        String::from_utf8(bytes)
            .map(Some)
            .map_err(|e| BioIoError::new(IoErrorKind::InvalidStructure, self.format, e.to_string()))
    }

    /// Returns the next line with its original `\n` or `\r\n` ending, if any.
    pub(crate) fn next_physical_line(&mut self) -> Result<Option<String>, BioIoError> {
        let Some(bytes) = self.next_physical_bytes()? else {
            return Ok(None);
        };
        String::from_utf8(bytes)
            .map(Some)
            .map_err(|e| BioIoError::new(IoErrorKind::InvalidStructure, self.format, e.to_string()))
    }

    /// Reads the next line into `line`, without its `\n` or `\r\n` ending, and
    /// rejects invalid UTF-8 exactly as [`Self::next_line`] does. The buffer is
    /// reused, so reading allocates only when a line outgrows it. Returns
    /// `false` at the end of the input.
    pub(crate) fn next_line_into(&mut self, line: &mut Vec<u8>) -> Result<bool, BioIoError> {
        line.clear();
        if !self.read_physical(line)? {
            return Ok(false);
        }
        if line.last() == Some(&b'\n') {
            line.pop();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
        }
        if !line.is_ascii() {
            std::str::from_utf8(line).map_err(|e| {
                BioIoError::new(IoErrorKind::InvalidStructure, self.format, e.to_string())
            })?;
        }
        Ok(true)
    }

    fn next_physical_bytes(&mut self) -> Result<Option<Vec<u8>>, BioIoError> {
        let mut bytes = Vec::new();
        Ok(self.read_physical(&mut bytes)?.then_some(bytes))
    }

    /// Appends the next physical line, with its ending, to the empty `bytes`.
    fn read_physical(&mut self, bytes: &mut Vec<u8>) -> Result<bool, BioIoError> {
        self.last_start = self.total;
        loop {
            let available = match self.reader.fill_buf() {
                Ok(available) => available,
                // Interrupted is transient; retain any already-consumed prefix.
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(io_error(self.format, error)),
            };
            if available.is_empty() {
                break;
            }
            let newline = find_newline(available);
            let count = newline.map_or(available.len(), |index| index + 1);
            if count > self.limits.max_line_bytes.saturating_sub(bytes.len())
                || count > self.limits.max_total_bytes.saturating_sub(self.total)
            {
                return Err(BioIoError::new(
                    IoErrorKind::LimitExceeded,
                    self.format,
                    "line or total decompressed byte limit exceeded",
                ));
            }
            bytes.try_reserve(count).map_err(|e| {
                BioIoError::new(IoErrorKind::LimitExceeded, self.format, e.to_string())
            })?;
            bytes.extend_from_slice(&available[..count]);
            self.reader.consume(count);
            self.total += count;
            if newline.is_some() {
                break;
            }
        }
        if bytes.is_empty() {
            return Ok(false);
        }
        self.line += 1;
        Ok(true)
    }
}

/// The position of the first `\n`, testing eight bytes at a time.
fn find_newline(bytes: &[u8]) -> Option<usize> {
    const ONES: u64 = 0x0101_0101_0101_0101;
    const HIGHS: u64 = 0x8080_8080_8080_8080;
    const NEWLINES: u64 = 0x0a0a_0a0a_0a0a_0a0a;
    let mut chunks = bytes.chunks_exact(8);
    let mut offset = 0;
    for chunk in &mut chunks {
        let mut word = [0; 8];
        word.copy_from_slice(chunk);
        let word = u64::from_le_bytes(word) ^ NEWLINES;
        // Nonzero exactly when some byte of `word` is zero, that is, a newline.
        if word.wrapping_sub(ONES) & !word & HIGHS != 0 {
            return chunk
                .iter()
                .position(|byte| *byte == b'\n')
                .map(|index| offset + index);
        }
        offset += 8;
    }
    chunks
        .remainder()
        .iter()
        .position(|byte| *byte == b'\n')
        .map(|index| offset + index)
}

// Document readers intentionally buffer a bounded document; sequence readers
// below do not call this helper and only retain the current record.
pub(crate) fn read_text(
    mut reader: impl BufRead,
    format: &'static str,
) -> Result<String, BioIoError> {
    let limit = ReadLimits::default().max_total_bytes;
    let mut bytes = Vec::new();
    reader
        .by_ref()
        .take(limit as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| io_error(format, e))?;
    if bytes.len() > limit {
        return Err(BioIoError::new(
            IoErrorKind::LimitExceeded,
            format,
            "document byte limit exceeded",
        ));
    }
    String::from_utf8(bytes)
        .map_err(|e| BioIoError::new(IoErrorKind::InvalidStructure, format, e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::find_newline;

    #[test]
    fn word_wise_newline_search_equals_a_byte_scan() {
        // Every newline position and none, in buffers of 0 to 40 bytes whose
        // other bytes include values next to `\n` and high bytes.
        for length in 0..=40 {
            for newline in (0..length).map(Some).chain([None]) {
                for filler in [0x00, 0x09, 0x0b, 0x8a, 0xff, b'A'] {
                    let mut bytes = vec![filler; length];
                    if let Some(position) = newline {
                        bytes[position] = b'\n';
                        if position + 3 < length {
                            bytes[position + 3] = b'\n';
                        }
                    }
                    assert_eq!(
                        find_newline(&bytes),
                        bytes.iter().position(|byte| *byte == b'\n'),
                        "{bytes:?}"
                    );
                }
            }
        }
    }
}
