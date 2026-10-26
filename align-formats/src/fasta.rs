use core::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use align_core::{Alphabet, Sequence, SequenceError};
use align_msa::{MsaError, MultipleAlignment, SequenceRecord};

const OUTPUT_LINE_WIDTH: usize = 80;

/// A malformed FASTA document or failed FASTA input/output operation.
#[derive(Debug)]
#[non_exhaustive]
pub enum FastaError {
    /// Opening, reading, or writing the underlying stream failed.
    Io(io::Error),
    /// Sequence content appeared before the first header.
    SequenceBeforeHeader {
        /// One-based input line.
        line: usize,
    },
    /// A header contained no identifier.
    EmptyHeader {
        /// One-based input line.
        line: usize,
    },
    /// A record contains a symbol outside the selected alphabet.
    InvalidSequence {
        /// Header of the affected record.
        identifier: String,
        /// Alphabet validation failure.
        source: SequenceError,
    },
    /// A parsed header violates the shared record contract.
    InvalidRecord {
        /// One-based header line.
        line: usize,
        /// Record validation failure.
        source: MsaError,
    },
    /// The stream contained no FASTA records.
    EmptyInput,
}

impl fmt::Display for FastaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "I/O error: {error}"),
            Self::SequenceBeforeHeader { line } => {
                write!(
                    formatter,
                    "sequence data occurs before a FASTA header at line {line}"
                )
            }
            Self::EmptyHeader { line } => {
                write!(formatter, "FASTA header at line {line} has no identifier")
            }
            Self::InvalidSequence { identifier, source } => {
                write!(
                    formatter,
                    "record {identifier:?} has an invalid sequence: {source}"
                )
            }
            Self::InvalidRecord { line, source } => {
                write!(
                    formatter,
                    "FASTA header at line {line} is invalid: {source}"
                )
            }
            Self::EmptyInput => formatter.write_str("FASTA input contains no records"),
        }
    }
}

impl std::error::Error for FastaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::InvalidSequence { source, .. } => Some(source),
            Self::InvalidRecord { source, .. } => Some(source),
            Self::SequenceBeforeHeader { .. } | Self::EmptyHeader { .. } | Self::EmptyInput => None,
        }
    }
}

impl From<io::Error> for FastaError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// Reads and validates every record in a FASTA file.
///
/// The complete trimmed header is retained as the stable record identifier.
/// ASCII whitespace in sequence lines is ignored and symbols are normalized by
/// [`Sequence::from_ascii`].
///
/// # Errors
///
/// Returns [`FastaError`] for failed I/O, malformed structure, an empty document,
/// an invalid identifier, or a symbol outside alphabet `A`.
pub fn read_fasta<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<Vec<SequenceRecord<A>>, FastaError> {
    let input = BufReader::new(File::open(path)?);
    parse_fasta(input)
}

/// Parses and validates every record from a buffered FASTA stream.
///
/// # Errors
///
/// Returns [`FastaError`] for failed I/O, malformed structure, an empty document,
/// an invalid identifier, or a symbol outside alphabet `A`.
pub fn parse_fasta<A: Alphabet>(
    reader: impl BufRead,
) -> Result<Vec<SequenceRecord<A>>, FastaError> {
    let mut records = Vec::new();
    let mut current: Option<(String, usize)> = None;
    let mut symbols = Vec::new();

    for (zero_based_line, line) in reader.lines().enumerate() {
        let line_number = zero_based_line + 1;
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some((identifier, header_line)) = current.take() {
                records.push(build_record(identifier, header_line, &symbols)?);
                symbols.clear();
            }
            let identifier = header.trim();
            if identifier.is_empty() {
                return Err(FastaError::EmptyHeader { line: line_number });
            }
            current = Some((identifier.to_owned(), line_number));
        } else if line.trim().is_empty() {
            continue;
        } else if current.is_none() {
            return Err(FastaError::SequenceBeforeHeader { line: line_number });
        } else {
            symbols.extend(line.bytes().filter(|byte| !byte.is_ascii_whitespace()));
        }
    }

    if let Some((identifier, header_line)) = current {
        records.push(build_record(identifier, header_line, &symbols)?);
    }
    if records.is_empty() {
        return Err(FastaError::EmptyInput);
    }
    Ok(records)
}

fn build_record<A: Alphabet>(
    identifier: String,
    header_line: usize,
    symbols: &[u8],
) -> Result<SequenceRecord<A>, FastaError> {
    let sequence = Sequence::from_ascii(symbols).map_err(|source| FastaError::InvalidSequence {
        identifier: identifier.clone(),
        source,
    })?;
    SequenceRecord::new(identifier, sequence).map_err(|source| FastaError::InvalidRecord {
        line: header_line,
        source,
    })
}

/// Writes ungapped records to a FASTA file with 80 symbols per line.
///
/// # Errors
///
/// Returns [`FastaError::Io`] if the file cannot be created or written.
pub fn write_fasta<A: Alphabet>(
    path: impl AsRef<Path>,
    records: &[SequenceRecord<A>],
) -> Result<(), FastaError> {
    let mut output = BufWriter::new(File::create(path)?);
    write_fasta_to(&mut output, records)?;
    output.flush()?;
    Ok(())
}

/// Writes ungapped records to a FASTA stream with 80 symbols per line.
///
/// # Errors
///
/// Returns [`FastaError::Io`] if the stream cannot be written.
pub fn write_fasta_to<A: Alphabet>(
    mut writer: impl Write,
    records: &[SequenceRecord<A>],
) -> Result<(), FastaError> {
    for record in records {
        writeln!(writer, ">{}", record.identifier())?;
        write_wrapped(&mut writer, record.sequence().to_ascii_string().as_bytes())?;
    }
    Ok(())
}

/// Writes a rectangular multiple alignment to a FASTA file.
///
/// # Errors
///
/// Returns [`FastaError::Io`] if the file cannot be created or written.
pub fn write_alignment_fasta<A: Alphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FastaError> {
    let mut output = BufWriter::new(File::create(path)?);
    write_alignment_fasta_to(&mut output, alignment)?;
    output.flush()?;
    Ok(())
}

/// Writes a rectangular multiple alignment to a FASTA stream.
///
/// # Errors
///
/// Returns [`FastaError::Io`] if the stream cannot be written.
pub fn write_alignment_fasta_to<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FastaError> {
    for row in alignment.rows() {
        writeln!(writer, ">{}", row.identifier())?;
        write_wrapped(&mut writer, row.to_ascii_with_gaps().as_bytes())?;
    }
    Ok(())
}

fn write_wrapped(mut writer: impl Write, symbols: &[u8]) -> Result<(), io::Error> {
    if symbols.is_empty() {
        writer.write_all(b"\n")?;
        return Ok(());
    }
    for chunk in symbols.chunks(OUTPUT_LINE_WIDTH) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Convenient `.write_fasta(path)` output for records and alignments.
pub trait FastaWrite {
    /// Writes `self` to `path` in FASTA representation.
    ///
    /// # Errors
    ///
    /// Returns [`FastaError::Io`] if the file cannot be created or written.
    fn write_fasta(&self, path: impl AsRef<Path>) -> Result<(), FastaError>;
}

impl<A: Alphabet> FastaWrite for [SequenceRecord<A>] {
    fn write_fasta(&self, path: impl AsRef<Path>) -> Result<(), FastaError> {
        write_fasta(path, self)
    }
}

impl<A: Alphabet> FastaWrite for Vec<SequenceRecord<A>> {
    fn write_fasta(&self, path: impl AsRef<Path>) -> Result<(), FastaError> {
        write_fasta(path, self)
    }
}

impl<A: Alphabet> FastaWrite for MultipleAlignment<A> {
    fn write_fasta(&self, path: impl AsRef<Path>) -> Result<(), FastaError> {
        write_alignment_fasta(path, self)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use align_core::{Dna, DnaSequence};
    use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

    use super::{FastaError, parse_fasta, write_alignment_fasta_to, write_fasta_to};

    #[test]
    fn parses_multiline_lowercase_fasta() {
        let input = Cursor::new(b">first description\nac gt\nN\n\n>second\nTTAA\n");
        let records = parse_fasta::<Dna>(input).expect("valid FASTA");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].identifier(), "first description");
        assert_eq!(records[0].sequence().to_ascii_string(), "ACGTN");
        assert_eq!(records[1].sequence().to_ascii_string(), "TTAA");
    }

    #[test]
    fn reports_sequence_before_header() {
        let error = parse_fasta::<Dna>(Cursor::new(b"ACGT\n")).expect_err("header is required");
        assert!(matches!(
            error,
            FastaError::SequenceBeforeHeader { line: 1 }
        ));
    }

    #[test]
    fn reports_empty_input_header_and_invalid_symbols() {
        assert!(matches!(
            parse_fasta::<Dna>(Cursor::new(b"\n")),
            Err(FastaError::EmptyInput)
        ));
        assert!(matches!(
            parse_fasta::<Dna>(Cursor::new(b">   \n")),
            Err(FastaError::EmptyHeader { line: 1 })
        ));
        assert!(matches!(
            parse_fasta::<Dna>(Cursor::new(b">bad\nACU\n")),
            Err(FastaError::InvalidSequence { identifier, .. }) if identifier == "bad"
        ));
    }

    #[test]
    fn writes_records_and_alignment() {
        let records = vec![
            SequenceRecord::new("a", DnaSequence::from_ascii("ACGT").unwrap()).unwrap(),
            SequenceRecord::new("b", DnaSequence::from_ascii("AGT").unwrap()).unwrap(),
        ];
        let mut fasta = Vec::new();
        write_fasta_to(&mut fasta, &records).unwrap();
        assert_eq!(String::from_utf8(fasta).unwrap(), ">a\nACGT\n>b\nAGT\n");

        let alignment = MultipleAlignment::try_new(vec![
            AlignedRow::<Dna>::new(
                "a",
                "ACGT".bytes().map(align_core::Symbol::from_ascii).collect(),
            )
            .unwrap(),
            AlignedRow::<Dna>::new(
                "b",
                vec![
                    align_core::Symbol::from_ascii(b'A'),
                    None,
                    align_core::Symbol::from_ascii(b'G'),
                    align_core::Symbol::from_ascii(b'T'),
                ],
            )
            .unwrap(),
        ])
        .unwrap();
        let mut aligned = Vec::new();
        write_alignment_fasta_to(&mut aligned, &alignment).unwrap();
        assert_eq!(String::from_utf8(aligned).unwrap(), ">a\nACGT\n>b\nA-GT\n");
    }
}
