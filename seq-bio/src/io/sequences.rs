use std::io::{BufRead, Write};
use std::iter::FusedIterator;
use std::marker::PhantomData;
use std::path::Path;

use crate::Record;
use crate::prelude::Alphabet;
use crate::record::split_header;

use super::{BioIoError, IoErrorKind, Lines, ReadLimits, RecordSpan, io_error, open_input};

/// A lazy FASTA iterator retaining only the current record and one header.
///
/// A record keeps its source text and its normalized sequence, so memory peaks
/// at about twice the largest record.
///
/// Blank lines and ASCII whitespace within sequence lines are ignored. Empty
/// sequences are allowed; an empty document yields no records. After any error
/// the iterator is fused. The legacy eager `read_fasta` is unchanged.
pub struct FastaReader<R, A: Alphabet> {
    lines: Lines<R>,
    buffer: Vec<u8>,
    pending: Option<(String, usize, u64)>,
    last_span: Option<RecordSpan>,
    ordinal: usize,
    finished: bool,
    alphabet: PhantomData<A>,
}

impl<R: BufRead, A: Alphabet> FastaReader<R, A> {
    /// Constructs a reader with default resource limits; performs no reading.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self::with_limits(reader, ReadLimits::default())
    }

    /// Constructs a reader with explicit limits; performs no reading.
    #[must_use]
    pub fn with_limits(reader: R, limits: ReadLimits) -> Self {
        Self {
            lines: Lines::new(reader, limits, "fasta"),
            buffer: Vec::new(),
            pending: None,
            last_span: None,
            ordinal: 0,
            finished: false,
            alphabet: PhantomData,
        }
    }

    fn read_record(&mut self) -> Result<Option<Record<A>>, BioIoError> {
        self.ordinal += 1;
        let (header, header_line, start) = if let Some(pending) = self.pending.take() {
            pending
        } else {
            loop {
                let Some(line) = self.line()? else {
                    return Ok(None);
                };
                if line.trim().is_empty() {
                    continue;
                }
                let Some(header) = line.strip_prefix('>') else {
                    return Err(self.error(IoErrorKind::InvalidHeader, "expected a '>' header"));
                };
                break (header.to_owned(), self.lines.line, self.lines.line_start());
            }
        };
        // Validate identity immediately without consuming the following record.
        let parts = split_header(&header)
            .map_err(|e| e.in_format("fasta").at(self.ordinal, header_line))?;
        let mut sequence = Vec::new();
        while self.line_into()? {
            if self.buffer.first() == Some(&b'>') {
                let header = self.buffer_text(1)?;
                self.pending = Some((header, self.lines.line, self.lines.line_start()));
                break;
            }
            // A line of valid symbols that fits the limit is copied whole; any
            // other line takes the byte loop, which reports the same errors.
            if self.buffer.len()
                <= self
                    .lines
                    .limits
                    .max_record_symbols
                    .saturating_sub(sequence.len())
                && self
                    .buffer
                    .iter()
                    .all(|byte| !byte.is_ascii_whitespace() && A::index_of(*byte).is_some())
            {
                sequence.extend_from_slice(&self.buffer);
                continue;
            }
            for byte in self
                .buffer
                .iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
            {
                if A::index_of(byte).is_none() {
                    return Err(self.error(
                        IoErrorKind::InvalidSequence,
                        format!("invalid {} byte {byte}", A::NAME),
                    ));
                }
                if sequence.len() >= self.lines.limits.max_record_symbols {
                    return Err(
                        self.error(IoErrorKind::LimitExceeded, "record symbol limit exceeded")
                    );
                }
                sequence.push(byte);
            }
        }
        let record = Record::from_checked_header(header, &parts, sequence, None)
            .map_err(|e| e.in_format("fasta").at(self.ordinal, header_line))?;
        self.last_span = Some(RecordSpan {
            start,
            end: self
                .pending
                .as_ref()
                .map_or_else(|| self.lines.byte_position(), |pending| pending.2),
            header_line,
        });
        Ok(Some(record))
    }

    pub(crate) fn record_span(&self) -> Option<RecordSpan> {
        self.last_span
    }

    fn line(&mut self) -> Result<Option<String>, BioIoError> {
        let line = self.lines.line + 1;
        self.lines.next_line().map_err(|e| e.at(self.ordinal, line))
    }

    /// Reads the next line into the reused buffer, as [`Self::line`] would.
    fn line_into(&mut self) -> Result<bool, BioIoError> {
        let line = self.lines.line + 1;
        self.lines
            .next_line_into(&mut self.buffer)
            .map_err(|e| e.at(self.ordinal, line))
    }

    /// The buffered line from byte `from` on, which reading checked is UTF-8.
    fn buffer_text(&self, from: usize) -> Result<String, BioIoError> {
        std::str::from_utf8(&self.buffer[from..])
            .map(str::to_owned)
            .map_err(|e| self.error(IoErrorKind::InvalidStructure, e.to_string()))
    }

    fn error(&self, kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
        BioIoError::new(kind, "fasta", message).at(self.ordinal, self.lines.line)
    }
}

impl<R: BufRead, A: Alphabet> Iterator for FastaReader<R, A> {
    type Item = Result<Record<A>, BioIoError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.last_span = None;
        if self.finished {
            return None;
        }
        match self.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}
impl<R: BufRead, A: Alphabet> FusedIterator for FastaReader<R, A> {}

/// A lazy Sanger FASTQ (Phred+33) reader with length-delimited qualities.
///
/// A record keeps its source text, normalized sequence and scores; the reader
/// keeps only the current record.
///
/// Wrapped sequence/quality lines and CRLF are supported. Empty reads, blank
/// lines and whitespace inside sequence/quality lines are rejected. A nonempty
/// `+` header must repeat the complete `@` header exactly. Phred+64/Solexa are not
/// detected or converted: the caller must select the correct encoding.
pub struct FastqReader<R, A: Alphabet> {
    lines: Lines<R>,
    buffer: Vec<u8>,
    last_span: Option<RecordSpan>,
    ordinal: usize,
    finished: bool,
    alphabet: PhantomData<A>,
}

impl<R: BufRead, A: Alphabet> FastqReader<R, A> {
    /// Constructs a reader with default resource limits; performs no reading.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self::with_limits(reader, ReadLimits::default())
    }

    /// Constructs a reader with explicit limits; performs no reading.
    #[must_use]
    pub fn with_limits(reader: R, limits: ReadLimits) -> Self {
        Self {
            lines: Lines::new(reader, limits, "fastq-sanger"),
            buffer: Vec::new(),
            last_span: None,
            ordinal: 0,
            finished: false,
            alphabet: PhantomData,
        }
    }

    fn read_record(&mut self) -> Result<Option<Record<A>>, BioIoError> {
        self.ordinal += 1;
        let start = self.lines.byte_position();
        if !self.line_into()? {
            return Ok(None);
        }
        if self.buffer.first() != Some(&b'@') {
            return Err(self.error(IoErrorKind::InvalidHeader, "expected an '@' header"));
        }
        let header = self.buffer_text(1)?;
        let header_line = self.lines.line;
        let parts = split_header(&header)
            .map_err(|e| e.in_format("fastq-sanger").at(self.ordinal, header_line))?;
        let mut sequence = Vec::new();
        loop {
            self.required_line_into("missing '+' separator")?;
            if self.buffer.first() == Some(&b'+') {
                let repeated = &self.buffer[1..];
                if !repeated.is_empty() && repeated != header.as_bytes() {
                    return Err(self.error(
                        IoErrorKind::InvalidHeader,
                        "'+' header must be empty or exactly repeat the '@' header",
                    ));
                }
                break;
            }
            if self.buffer.is_empty() {
                return Err(self.error(IoErrorKind::InvalidSequence, "empty sequence line"));
            }
            if self.buffer.len()
                > self
                    .lines
                    .limits
                    .max_record_symbols
                    .saturating_sub(sequence.len())
            {
                return Err(self.error(IoErrorKind::LimitExceeded, "record symbol limit exceeded"));
            }
            if let Some(byte) = self
                .buffer
                .iter()
                .copied()
                .find(|byte| A::index_of(*byte).is_none())
            {
                return Err(self.error(
                    IoErrorKind::InvalidSequence,
                    format!("invalid {} byte {byte}", A::NAME),
                ));
            }
            sequence.extend_from_slice(&self.buffer);
        }
        if sequence.is_empty() {
            return Err(self.error(
                IoErrorKind::Unsupported,
                "empty FASTQ reads are not supported",
            ));
        }
        let mut qualities = Vec::with_capacity(sequence.len());
        while qualities.len() < sequence.len() {
            self.required_line_into("quality data ended before the sequence length")?;
            if self.buffer.is_empty() || self.buffer.iter().any(|byte| !(33..=126).contains(byte)) {
                return Err(self.error(
                    IoErrorKind::InvalidQuality,
                    "quality characters must be ASCII 33..=126",
                ));
            }
            if self.buffer.len() > sequence.len() - qualities.len() {
                return Err(self.error(
                    IoErrorKind::LengthMismatch,
                    "quality count exceeds sequence length",
                ));
            }
            qualities.extend(self.buffer.iter().map(|byte| byte - 33));
        }
        let record = Record::from_checked_header(header, &parts, sequence, Some(qualities))
            .map_err(|e| e.in_format("fastq-sanger").at(self.ordinal, header_line))?;
        self.last_span = Some(RecordSpan {
            start,
            end: self.lines.byte_position(),
            header_line,
        });
        Ok(Some(record))
    }

    pub(crate) fn record_span(&self) -> Option<RecordSpan> {
        self.last_span
    }

    /// Reads the next line into the reused buffer.
    fn line_into(&mut self) -> Result<bool, BioIoError> {
        let line = self.lines.line + 1;
        self.lines
            .next_line_into(&mut self.buffer)
            .map_err(|e| e.at(self.ordinal, line))
    }

    fn required_line_into(&mut self, message: &str) -> Result<(), BioIoError> {
        if self.line_into()? {
            return Ok(());
        }
        Err(
            BioIoError::new(IoErrorKind::Truncated, "fastq-sanger", message)
                .at(self.ordinal, self.lines.line + 1),
        )
    }

    /// The buffered line from byte `from` on, which reading checked is UTF-8.
    fn buffer_text(&self, from: usize) -> Result<String, BioIoError> {
        std::str::from_utf8(&self.buffer[from..])
            .map(str::to_owned)
            .map_err(|e| self.error(IoErrorKind::InvalidStructure, e.to_string()))
    }

    fn error(&self, kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
        BioIoError::new(kind, "fastq-sanger", message).at(self.ordinal, self.lines.line)
    }
}

impl<R: BufRead, A: Alphabet> Iterator for FastqReader<R, A> {
    type Item = Result<Record<A>, BioIoError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.last_span = None;
        if self.finished {
            return None;
        }
        match self.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(error) => {
                self.finished = true;
                Some(Err(error))
            }
        }
    }
}
impl<R: BufRead, A: Alphabet> FusedIterator for FastqReader<R, A> {}

/// Opens a plain/gzip FASTA file as a lazy record iterator.
///
/// # Errors
/// Returns opening errors immediately; parse/decompression errors are items.
pub fn fasta<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<FastaReader<Box<dyn BufRead>, A>, BioIoError> {
    Ok(FastaReader::new(open_input(path)?))
}

/// Opens a plain/gzip Sanger FASTQ file as a lazy record iterator.
///
/// # Errors
/// Returns opening errors immediately; parse/decompression errors are items.
pub fn fastq<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<FastqReader<Box<dyn BufRead>, A>, BioIoError> {
    Ok(FastqReader::new(open_input(path)?))
}

/// Consumes a reader that must contain exactly one record, including its EOF.
///
/// # Errors
/// Rejects empty/multiple-record input and propagates any encountered input error.
pub fn one_record<A: Alphabet>(
    mut records: impl Iterator<Item = Result<Record<A>, BioIoError>>,
) -> Result<Record<A>, BioIoError> {
    let first = records.next().ok_or_else(|| {
        BioIoError::new(
            IoErrorKind::LengthMismatch,
            "records",
            "expected one record, found none",
        )
    })??;
    if let Some(next) = records.next() {
        next?;
        return Err(BioIoError::new(
            IoErrorKind::LengthMismatch,
            "records",
            "expected one record, found more than one",
        ));
    }
    Ok(first)
}

/// Explicit permission to drop quality scores when writing FASTA.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LossPolicy {
    /// Reject a conversion that would discard quality scores.
    Reject,
    /// Permit discarding quality scores; never fabricate missing information.
    Allow,
}

/// Writes one FASTA record with 80-symbol lines, preserving header and case.
///
/// # Errors
/// Rejects quality loss unless explicitly allowed, before writing any bytes.
/// A writer I/O failure may leave partial output; path-level atomicity is the
/// caller's responsibility. Line wrapping and source line endings are not kept.
pub fn write_record_fasta<A: Alphabet>(
    mut writer: impl Write,
    record: &Record<A>,
    loss: LossPolicy,
) -> Result<(), BioIoError> {
    if record.qualities().is_some() && loss == LossPolicy::Reject {
        return Err(BioIoError::new(
            IoErrorKind::DataLoss,
            "fasta",
            "FASTA cannot retain quality scores; opt in to discarding them",
        ));
    }
    writeln!(writer, ">{}", record.header()).map_err(|e| io_error("fasta", e))?;
    if record.source_sequence().is_empty() {
        writer.write_all(b"\n").map_err(|e| io_error("fasta", e))?;
    }
    for chunk in record.source_sequence().as_bytes().chunks(80) {
        writer
            .write_all(chunk)
            .and_then(|()| writer.write_all(b"\n"))
            .map_err(|e| io_error("fasta", e))?;
    }
    Ok(())
}

/// Writes one Sanger FASTQ record with 80-symbol sequence and quality lines.
///
/// # Errors
/// Rejects absent quality scores or an empty read before writing any bytes.
/// Does not fabricate scores. Writer I/O failures may leave partial output.
pub fn write_record_fastq<A: Alphabet>(
    mut writer: impl Write,
    record: &Record<A>,
) -> Result<(), BioIoError> {
    let scores = record.qualities().ok_or_else(|| {
        BioIoError::new(
            IoErrorKind::DataLoss,
            "fastq-sanger",
            "FASTQ requires quality scores; none are available",
        )
    })?;
    if scores.is_empty() {
        return Err(BioIoError::new(
            IoErrorKind::Unsupported,
            "fastq-sanger",
            "empty FASTQ reads are not supported",
        ));
    }
    writeln!(writer, "@{}", record.header()).map_err(|e| io_error("fastq-sanger", e))?;
    for chunk in record.source_sequence().as_bytes().chunks(80) {
        writer
            .write_all(chunk)
            .and_then(|()| writer.write_all(b"\n"))
            .map_err(|e| io_error("fastq-sanger", e))?;
    }
    writeln!(writer, "+").map_err(|e| io_error("fastq-sanger", e))?;
    for chunk in scores.chunks(80) {
        let encoded: Vec<u8> = chunk.iter().map(|score| score + 33).collect();
        writer
            .write_all(&encoded)
            .and_then(|()| writer.write_all(b"\n"))
            .map_err(|e| io_error("fastq-sanger", e))?;
    }
    Ok(())
}

#[cfg(test)]
mod span_tests {
    use std::io::{BufReader, Cursor};

    use crate::Dna;

    use super::*;

    #[test]
    fn fasta_spans_preserve_raw_offsets_across_lookahead_and_empty_records() {
        let leading = "\r\n \t\n";
        let first = ">récord description\r\nac G\r\n\r\n";
        let second = ">empty\r\n";
        let third = ">last\nT";
        let text = format!("{leading}{first}{second}{third}");
        let mut reader =
            FastaReader::<_, Dna>::new(BufReader::with_capacity(1, Cursor::new(text.as_bytes())));
        assert_eq!(reader.record_span(), None);
        let mut start = leading.len() as u64;
        for (raw, line, id, sequence) in [
            (first, 3, "récord", "acG"),
            (second, 6, "empty", ""),
            (third, 7, "last", "T"),
        ] {
            let record = reader.next().unwrap().unwrap();
            assert_eq!(record.id(), id);
            assert_eq!(record.source_sequence(), sequence);
            let end = start + raw.len() as u64;
            assert_eq!(
                reader.record_span(),
                Some(RecordSpan {
                    start,
                    end,
                    header_line: line,
                })
            );
            assert_eq!(
                &text.as_bytes()[usize::try_from(start).unwrap()..usize::try_from(end).unwrap()],
                raw.as_bytes()
            );
            start = end;
        }
        assert_eq!(start, text.len() as u64);
        assert!(reader.next().is_none());
        assert_eq!(reader.record_span(), None);
        assert!(reader.next().is_none());
        assert_eq!(reader.record_span(), None);
    }

    #[test]
    fn fastq_spans_are_quality_length_delimited_not_header_delimited() {
        let first = "@one desc\r\nac\r\nGT\r\n+one desc\r\n@+\r\n!!\r\n";
        let second = "@two\nN\n+\nI";
        let text = format!("{first}{second}");
        let mut reader =
            FastqReader::<_, Dna>::new(BufReader::with_capacity(1, Cursor::new(text.as_bytes())));
        assert_eq!(reader.record_span(), None);
        let one = reader.next().unwrap().unwrap();
        assert_eq!(one.source_sequence(), "acGT");
        assert_eq!(one.qualities(), Some([31, 10, 0, 0].as_slice()));
        assert_eq!(
            reader.record_span(),
            Some(RecordSpan {
                start: 0,
                end: first.len() as u64,
                header_line: 1,
            })
        );
        assert_eq!(reader.next().unwrap().unwrap().id(), "two");
        assert_eq!(
            reader.record_span(),
            Some(RecordSpan {
                start: first.len() as u64,
                end: text.len() as u64,
                header_line: 7,
            })
        );
        assert!(reader.next().is_none());
        assert_eq!(reader.record_span(), None);
    }

    #[test]
    fn failed_or_empty_reads_do_not_expose_a_previous_record_span() {
        let mut sequences = FastaReader::<_, Dna>::new(Cursor::new(">ok\nA\n>bad\n!\n"));
        assert!(sequences.next().unwrap().is_ok());
        assert!(sequences.record_span().is_some());
        assert!(sequences.next().unwrap().is_err());
        assert_eq!(sequences.record_span(), None);
        assert!(sequences.next().is_none());

        let mut reads = FastqReader::<_, Dna>::new(Cursor::new("@ok\nA\n+\n!\n@bad\nA\n+\n"));
        assert!(reads.next().unwrap().is_ok());
        assert!(reads.record_span().is_some());
        assert!(reads.next().unwrap().is_err());
        assert_eq!(reads.record_span(), None);
        assert!(reads.next().is_none());

        let mut empty = FastaReader::<_, Dna>::new(Cursor::new("\r\n \t\n"));
        assert!(empty.next().is_none());
        assert_eq!(empty.record_span(), None);
    }
}
