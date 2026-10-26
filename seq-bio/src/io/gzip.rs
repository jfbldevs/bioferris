//! Deterministic gzip output with explicit, observable finalization.
//!
//! The raw compressor owns no I/O handle. Dropping this writer frees it and the
//! sink without emitting a hidden trailer or suppressing a finalization error.

use std::fmt;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use flate2::{Compress, Compression, FlushCompress, Status};

use crate::publication::{PublishError, publish_new};

const OUTPUT_BYTES: usize = 16 * 1024;

/// Ordinary gzip compression policy; the header never includes a clock or path.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GzipWriteOptions {
    /// Compression level in 0..=9, default 6.
    pub compression_level: u32,
}

impl Default for GzipWriteOptions {
    fn default() -> Self {
        Self {
            compression_level: 6,
        }
    }
}

impl GzipWriteOptions {
    /// Validates the compression level before any output is written.
    ///
    /// # Errors
    /// Returns `InvalidCompressionLevel` for levels above 9.
    pub fn validate(self) -> Result<(), GzipError> {
        if self.compression_level > 9 {
            return Err(GzipError::InvalidCompressionLevel {
                level: self.compression_level,
            });
        }
        Ok(())
    }
}

/// A single-member gzip writer using a bounded compressed-output scratch buffer.
///
/// Existing FASTA/FASTQ writers can write into it through `Write`. Call `finish`
/// and check its result. Dropping without finishing emits no hidden output.
#[derive(Debug)]
#[must_use = "call finish and check its result to finalize the gzip member"]
pub struct GzipWriter<W> {
    inner: W,
    compressor: Compress,
    crc: crc32fast::Hasher,
    input_bytes: u64,
    options: GzipWriteOptions,
    header_written: bool,
    failed: bool,
}

impl<W: Write> GzipWriter<W> {
    /// Creates an unwritten origin-independent gzip stream with default options.
    pub fn new(inner: W) -> Self {
        Self::configured(inner, GzipWriteOptions::default())
    }

    /// Creates a writer after checking explicit compression options.
    ///
    /// # Errors
    /// Returns `InvalidCompressionLevel` without writing any bytes.
    pub fn with_options(inner: W, options: GzipWriteOptions) -> Result<Self, GzipError> {
        options.validate()?;
        Ok(Self::configured(inner, options))
    }

    fn configured(inner: W, options: GzipWriteOptions) -> Self {
        Self {
            inner,
            compressor: Compress::new(Compression::new(options.compression_level), false),
            crc: crc32fast::Hasher::new(),
            input_bytes: 0,
            options,
            header_written: false,
            failed: false,
        }
    }

    /// Borrows the sink without allowing untracked writes.
    #[must_use]
    pub const fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns accepted uncompressed bytes, including data buffered by DEFLATE.
    #[must_use]
    pub const fn input_bytes(&self) -> u64 {
        self.input_bytes
    }

    /// Compresses a nonempty prefix and returns the number of bytes accepted.
    ///
    /// # Errors
    /// Returns typed compression, sink, overflow or failed-state errors. A sink
    /// may retain a prefix before a failure; this writer then refuses further use.
    pub fn write_bytes(&mut self, input: &[u8]) -> Result<usize, GzipError> {
        let result = self.write_inner(input);
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn write_inner(&mut self, input: &[u8]) -> Result<usize, GzipError> {
        self.check_state()?;
        if input.is_empty() {
            return Ok(0);
        }
        self.write_header()?;
        let mut output = [0; OUTPUT_BYTES];
        loop {
            let before_in = self.compressor.total_in();
            let before_out = self.compressor.total_out();
            self.compressor
                .compress(input, &mut output, FlushCompress::None)
                .map_err(|error| GzipError::Deflate(error.to_string()))?;
            let consumed = usize::try_from(
                self.compressor
                    .total_in()
                    .checked_sub(before_in)
                    .ok_or(GzipError::CountOverflow)?,
            )
            .map_err(|_| GzipError::CountOverflow)?;
            let produced = usize::try_from(
                self.compressor
                    .total_out()
                    .checked_sub(before_out)
                    .ok_or(GzipError::CountOverflow)?,
            )
            .map_err(|_| GzipError::CountOverflow)?;
            self.publish(&output[..produced], "write compressed data")?;
            if consumed > 0 {
                self.input_bytes = self
                    .input_bytes
                    .checked_add(u64::try_from(consumed).map_err(|_| GzipError::CountOverflow)?)
                    .ok_or(GzipError::CountOverflow)?;
                self.crc.update(&input[..consumed]);
                return Ok(consumed);
            }
            if produced == 0 {
                return Err(GzipError::Deflate("compressor made no progress".into()));
            }
        }
    }

    /// Emits buffered compressed bytes and flushes the sink without ending gzip.
    ///
    /// # Errors
    /// Returns compression, sink, overflow or failed-state errors.
    pub fn flush_compressed(&mut self) -> Result<(), GzipError> {
        let result = self.flush_inner();
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn flush_inner(&mut self) -> Result<(), GzipError> {
        self.check_state()?;
        self.write_header()?;
        self.pump(false)?;
        self.inner.flush().map_err(|source| GzipError::Io {
            operation: "flush output",
            source,
        })
    }

    /// Completes DEFLATE, writes CRC32/ISIZE, flushes and returns the sink.
    ///
    /// ISIZE is the uncompressed length modulo 2^32 as required by gzip. Header
    /// MTIME is zero, OS is 255 and optional text/extra fields are absent.
    ///
    /// # Errors
    /// Returns compression, sink, overflow or failed-state errors. No finalized
    /// sink is returned on failure, and Drop performs no later output attempts.
    pub fn finish(mut self) -> Result<W, GzipError> {
        self.check_state()?;
        self.write_header()?;
        self.pump(true)?;
        let mut trailer = [0; 8];
        trailer[..4].copy_from_slice(&self.crc.clone().finalize().to_le_bytes());
        trailer[4..].copy_from_slice(&self.input_bytes.to_le_bytes()[..4]);
        self.publish(&trailer, "write gzip trailer")?;
        self.inner.flush().map_err(|source| GzipError::Io {
            operation: "finalize output",
            source,
        })?;
        Ok(self.inner)
    }

    fn check_state(&self) -> Result<(), GzipError> {
        if self.failed {
            Err(GzipError::FailedState)
        } else {
            Ok(())
        }
    }

    fn write_header(&mut self) -> Result<(), GzipError> {
        if !self.header_written {
            let xfl = match self.options.compression_level {
                9 => 2,
                1 => 4,
                _ => 0,
            };
            self.publish(&[31, 139, 8, 0, 0, 0, 0, 0, xfl, 255], "write gzip header")?;
            self.header_written = true;
        }
        Ok(())
    }

    fn publish(&mut self, bytes: &[u8], operation: &'static str) -> Result<(), GzipError> {
        self.inner
            .write_all(bytes)
            .map_err(|source| GzipError::Io { operation, source })
    }

    fn pump(&mut self, finish: bool) -> Result<(), GzipError> {
        let mut output = [0; OUTPUT_BYTES];
        loop {
            let before = self.compressor.total_out();
            let flush = if finish {
                FlushCompress::Finish
            } else {
                FlushCompress::Sync
            };
            let status = self
                .compressor
                .compress(&[], &mut output, flush)
                .map_err(|error| GzipError::Deflate(error.to_string()))?;
            let produced = usize::try_from(
                self.compressor
                    .total_out()
                    .checked_sub(before)
                    .ok_or(GzipError::CountOverflow)?,
            )
            .map_err(|_| GzipError::CountOverflow)?;
            self.publish(&output[..produced], "flush compressed data")?;
            if finish && status == Status::StreamEnd {
                return Ok(());
            }
            if !finish && produced < output.len() {
                return Ok(());
            }
            if produced == 0 {
                return Err(GzipError::Deflate("compressor did not finish".into()));
            }
        }
    }
}

impl<W: Write> Write for GzipWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.write_bytes(bytes).map_err(as_io)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.flush_compressed().map_err(as_io)
    }
}

fn as_io(error: GzipError) -> io::Error {
    let kind = match &error {
        GzipError::Io { source, .. } => source.kind(),
        GzipError::InvalidCompressionLevel { .. } => io::ErrorKind::InvalidInput,
        _ => io::ErrorKind::Other,
    };
    io::Error::new(kind, error)
}

/// An explicit gzip computation or finalization failure.
#[derive(Debug)]
#[non_exhaustive]
pub enum GzipError {
    /// Compression level is not supported.
    InvalidCompressionLevel {
        /// Rejected compression level.
        level: u32,
    },
    /// An underlying sink operation failed.
    Io {
        /// Operation that failed.
        operation: &'static str,
        /// Original I/O error.
        source: io::Error,
    },
    /// The raw DEFLATE backend failed or made no progress.
    Deflate(String),
    /// An input/output count cannot be represented.
    CountOverflow,
    /// A previous write or flush failed.
    FailedState,
}

impl GzipError {
    /// Returns the stable machine-readable error category.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidCompressionLevel { .. } => "invalid_gzip_options",
            Self::Io { .. } => "gzip_io_error",
            Self::Deflate(_) => "gzip_deflate_error",
            Self::CountOverflow => "gzip_count_overflow",
            Self::FailedState => "gzip_failed_state",
        }
    }
}

impl fmt::Display for GzipError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCompressionLevel { level } => write!(
                formatter,
                "gzip compression level must be in 0..=9, received {level}"
            ),
            Self::Io { operation, source } => {
                write!(formatter, "gzip {operation} failed: {source}")
            }
            Self::Deflate(message) => write!(formatter, "gzip DEFLATE failed: {message}"),
            Self::CountOverflow => formatter.write_str("gzip input/output count overflow"),
            Self::FailedState => formatter.write_str("gzip writer is in a failed state"),
        }
    }
}

impl std::error::Error for GzipError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// A create-only gzip publication failure preserving the callback's own error.
#[derive(Debug)]
#[non_exhaustive]
pub enum GzipOutputError<E> {
    /// The compression engine or explicit finalization failed.
    Compression(GzipError),
    /// The biological/application writer callback failed.
    Write(E),
    /// Temporary file, synchronization or hard-link publication failed.
    Io(io::Error),
    /// A file, directory or symlink already occupies the requested destination.
    DestinationExists(PathBuf),
    /// No temporary publication filename could be reserved.
    TemporaryUnavailable,
}

impl<E: fmt::Display> fmt::Display for GzipOutputError<E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Compression(error) => error.fmt(formatter),
            Self::Write(error) => write!(formatter, "gzip data writer failed: {error}"),
            Self::Io(error) => error.fmt(formatter),
            Self::DestinationExists(path) => write!(
                formatter,
                "destination already exists; output is create-only: {}",
                path.display()
            ),
            Self::TemporaryUnavailable => {
                formatter.write_str("could not reserve temporary gzip output")
            }
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for GzipOutputError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Compression(error) => Some(error),
            Self::Write(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::DestinationExists(_) | Self::TemporaryUnavailable => None,
        }
    }
}

/// Writes and explicitly finalizes a new gzip file before create-only publication.
///
/// The callback receives a standard `Write` sink, suitable for existing biological
/// writers. The destination directory must be trusted. No existing destination
/// is overwritten, and failed temporary output is removed by the shared publisher.
///
/// # Errors
/// Returns separate callback, compression, existing-path and publication failures.
pub fn write_gzip_new<E>(
    path: impl AsRef<Path>,
    options: GzipWriteOptions,
    write: impl FnOnce(&mut dyn Write) -> Result<(), E>,
) -> Result<(), GzipOutputError<E>> {
    options.validate().map_err(GzipOutputError::Compression)?;
    publish_new(path.as_ref(), |output| {
        let mut writer =
            GzipWriter::with_options(output, options).map_err(GzipOutputError::Compression)?;
        write(&mut writer).map_err(GzipOutputError::Write)?;
        writer.finish().map_err(GzipOutputError::Compression)?;
        Ok(())
    })
    .map_err(|error| match error {
        PublishError::Write(error) => error,
        PublishError::Io(error) => GzipOutputError::Io(error),
        PublishError::DestinationExists(path) => GzipOutputError::DestinationExists(path),
        PublishError::TemporaryUnavailable => GzipOutputError::TemporaryUnavailable,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::io::Read;
    use std::rc::Rc;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn decompress(bytes: &[u8]) -> io::Result<Vec<u8>> {
        let mut result = Vec::new();
        flate2::read::MultiGzDecoder::new(bytes).read_to_end(&mut result)?;
        Ok(result)
    }

    fn encoded(input: &[u8], level: u32) -> Vec<u8> {
        let mut writer = GzipWriter::with_options(
            Vec::new(),
            GzipWriteOptions {
                compression_level: level,
            },
        )
        .unwrap();
        writer.write_all(input).unwrap();
        writer.finish().unwrap()
    }

    #[test]
    fn header_trailer_and_deterministic_bytes_follow_rfc1952() {
        for level in 0..=9 {
            let bytes = encoded(b"123456789", level);
            let xfl = match level {
                1 => 4,
                9 => 2,
                _ => 0,
            };
            assert_eq!(&bytes[..10], &[31, 139, 8, 0, 0, 0, 0, 0, xfl, 255]);
            assert_eq!(
                &bytes[bytes.len() - 8..],
                &[0x26, 0x39, 0xf4, 0xcb, 9, 0, 0, 0]
            );
            assert_eq!(decompress(&bytes).unwrap(), b"123456789");
            assert_eq!(encoded(b"123456789", level), bytes);
        }
    }

    #[test]
    fn empty_all_byte_and_large_inputs_round_trip_at_every_level() {
        let all_bytes: Vec<u8> = (0..=255).collect();
        let large: Vec<_> = (0..300_000)
            .map(|index| u8::try_from(index % 251).unwrap())
            .collect();
        for level in 0..=9 {
            for input in [b"".as_slice(), all_bytes.as_slice(), large.as_slice()] {
                assert_eq!(
                    decompress(&encoded(input, level)).unwrap(),
                    input,
                    "level={level}, length={}",
                    input.len()
                );
            }
        }
    }

    #[test]
    fn repeated_flushes_leave_the_stream_open_and_preserve_exact_data() {
        let mut writer = GzipWriter::new(Vec::new());
        writer.flush_compressed().unwrap();
        writer.flush_compressed().unwrap();
        writer.write_all(b"first\0").unwrap();
        writer.flush_compressed().unwrap();
        writer.write_all(b"second\r\n").unwrap();
        writer.flush_compressed().unwrap();
        writer.flush_compressed().unwrap();
        assert_eq!(writer.input_bytes(), 14);
        let output = writer.finish().unwrap();
        assert_eq!(decompress(&output).unwrap(), b"first\0second\r\n");
    }

    #[derive(Clone, Debug, Default)]
    struct Sink {
        bytes: Rc<RefCell<Vec<u8>>>,
        writes: Rc<RefCell<usize>>,
        fail_after: Option<usize>,
        fail_flush: bool,
        interrupted: bool,
    }

    impl Write for Sink {
        fn write(&mut self, input: &[u8]) -> io::Result<usize> {
            *self.writes.borrow_mut() += 1;
            if self.interrupted && *self.writes.borrow() % 7 == 0 {
                return Err(io::ErrorKind::Interrupted.into());
            }
            let present = self.bytes.borrow().len();
            if self.fail_after.is_some_and(|limit| present >= limit) {
                return Err(io::ErrorKind::PermissionDenied.into());
            }
            let count = input
                .len()
                .min(3)
                .min(self.fail_after.map_or(usize::MAX, |limit| limit - present));
            self.bytes.borrow_mut().extend_from_slice(&input[..count]);
            Ok(count)
        }
        fn flush(&mut self) -> io::Result<()> {
            if self.fail_flush {
                Err(io::ErrorKind::BrokenPipe.into())
            } else {
                Ok(())
            }
        }
    }

    #[test]
    fn short_and_interrupted_writes_finish_without_losing_bytes() {
        let sink = Sink {
            interrupted: true,
            ..Default::default()
        };
        let captured = sink.bytes.clone();
        let mut writer = GzipWriter::new(sink);
        assert_eq!(writer.write_bytes(&[]).unwrap(), 0);
        assert!(captured.borrow().is_empty());
        let data: Vec<u8> = (0..100_000)
            .map(|index| u8::try_from(index % 251).unwrap())
            .collect();
        writer.write_all(&data).unwrap();
        writer.finish().unwrap();
        assert_eq!(decompress(&captured.borrow()).unwrap(), data);
    }

    #[test]
    fn drop_performs_no_hidden_write_or_finalization() {
        let sink = Sink::default();
        let bytes = sink.bytes.clone();
        let shared_calls = sink.writes.clone();
        let mut writer = GzipWriter::new(sink);
        writer.write_all(b"unfinalized").unwrap();
        writer.flush_compressed().unwrap();
        let before = bytes.borrow().clone();
        let calls = *shared_calls.borrow();
        drop(writer);
        assert_eq!(*shared_calls.borrow(), calls);
        assert_eq!(*bytes.borrow(), before);
        assert!(decompress(&before).is_err());
    }

    #[test]
    fn write_flush_and_finalization_failures_preserve_causes_and_failed_state() {
        let mut writer = GzipWriter::new(Sink {
            fail_after: Some(0),
            ..Default::default()
        });
        assert!(matches!(
            writer.write_bytes(b"abc"),
            Err(GzipError::Io {
                operation: "write gzip header",
                ..
            })
        ));
        assert!(matches!(writer.finish(), Err(GzipError::FailedState)));
        let mut writer = GzipWriter::new(Sink {
            fail_flush: true,
            ..Default::default()
        });
        let error = writer.flush().unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
        assert!(matches!(
            error.get_ref().unwrap().downcast_ref::<GzipError>(),
            Some(GzipError::Io { .. })
        ));
        assert!(matches!(
            writer.write_bytes(b"abc"),
            Err(GzipError::FailedState)
        ));
        let writer = GzipWriter::new(Sink {
            fail_flush: true,
            ..Default::default()
        });
        assert!(matches!(
            writer.finish(),
            Err(GzipError::Io {
                operation: "finalize output",
                ..
            })
        ));
        let writer = GzipWriter::new(Sink {
            fail_after: Some(12),
            ..Default::default()
        });
        assert!(matches!(writer.finish(), Err(GzipError::Io { .. })));
    }

    #[test]
    fn invalid_options_and_count_overflow_do_not_wrap() {
        assert!(matches!(
            GzipWriter::with_options(
                Vec::new(),
                GzipWriteOptions {
                    compression_level: 10
                }
            ),
            Err(GzipError::InvalidCompressionLevel { level: 10 })
        ));
        let mut writer = GzipWriter::new(Vec::new());
        writer.input_bytes = u64::MAX;
        assert!(matches!(
            writer.write_bytes(b"a"),
            Err(GzipError::CountOverflow)
        ));
        assert!(matches!(
            writer.write_bytes(b"b"),
            Err(GzipError::FailedState)
        ));
        let mut writer = GzipWriter::new(Vec::new());
        writer.input_bytes = 1 << 32;
        let output = writer.finish().unwrap();
        assert_eq!(&output[output.len() - 4..], &[0; 4]);
    }

    struct Directory(PathBuf);
    impl Directory {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "bioferris-gzip-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Directory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn create_only_publication_refuses_existing_paths_and_cleans_callback_errors() {
        let directory = Directory::new();
        let path = directory.0.join("records.fa.gz");
        write_gzip_new(&path, GzipWriteOptions::default(), |output| {
            output.write_all(b">a\nACGT\n")
        })
        .unwrap();
        let original = std::fs::read(&path).unwrap();
        assert_eq!(decompress(&original).unwrap(), b">a\nACGT\n");
        let error = write_gzip_new(&path, GzipWriteOptions::default(), |output| {
            output.write_all(b"replacement")
        })
        .unwrap_err();
        assert!(matches!(error, GzipOutputError::DestinationExists(_)));
        assert_eq!(std::fs::read(&path).unwrap(), original);
        let failed = directory.0.join("failed.gz");
        let error = write_gzip_new(&failed, GzipWriteOptions::default(), |output| {
            output.write_all(b"prefix")?;
            Err(io::Error::other("callback failure"))
        })
        .unwrap_err();
        assert!(matches!(error, GzipOutputError::Write(_)));
        assert!(!failed.exists());
        assert_eq!(std::fs::read_dir(&directory.0).unwrap().count(), 1);
    }

    #[cfg(unix)]
    #[test]
    fn create_only_publication_refuses_dangling_symlinks() {
        let directory = Directory::new();
        let path = directory.0.join("dangling.gz");
        std::os::unix::fs::symlink(directory.0.join("missing"), &path).unwrap();
        assert!(matches!(
            write_gzip_new(&path, GzipWriteOptions::default(), |output| output
                .write_all(b"data")),
            Err(GzipOutputError::DestinationExists(_))
        ));
        assert!(std::fs::symlink_metadata(&path).unwrap().is_symlink());
    }

    #[test]
    fn existing_record_writers_compose_without_losing_case_headers_or_qualities() {
        let record =
            crate::Record::<crate::Dna>::from_ascii("read1 coding", b"atgGCT", Some(vec![40; 6]))
                .unwrap();
        let mut writer = GzipWriter::new(Vec::new());
        crate::io::write_record_fastq(&mut writer, &record).unwrap();
        let output = writer.finish().unwrap();
        assert_eq!(
            decompress(&output).unwrap(),
            b"@read1 coding\natgGCT\n+\nIIIIII\n"
        );
        let directory = Directory::new();
        let path = directory.0.join("forbidden.fa.gz");
        let error = write_gzip_new(&path, GzipWriteOptions::default(), |output| {
            crate::io::write_record_fasta(output, &record, crate::io::LossPolicy::Reject)
        })
        .unwrap_err();
        assert!(matches!(error, GzipOutputError::Write(source) if source.code() == "data_loss"));
        assert!(!path.exists());
    }
}
