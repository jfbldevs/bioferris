//! Conservative `GenBank` and EMBL documents, writing and checked transforms.
//!
//! [`parse_genbank`], [`read_genbank`] and their EMBL counterparts buffer a
//! bounded document; [`AnnotatedReader`], opened by [`genbank`] or [`embl`],
//! streams one record at a time.
//!
//! Exact locations use zero-based, half-open ranges internally. Fuzzy, remote,
//! between-base and unrecognized locations are retained, not guessed. See the
//! repository's `docs/api/annotated-input-v1.md` for the supported profile.
//!
//! ```
//! use bioferris::{Dna, io::annotated::parse_genbank};
//!
//! let input = r#"LOCUS       DEMO 9 bp DNA linear SYN 15-OCT-2026
//! DEFINITION  Example coding sequence.
//! ACCESSION   DEMO
//! VERSION     DEMO.1
//! FEATURES             Location/Qualifiers
//!      CDS             1..9
//!                      /gene="demo"
//! ORIGIN
//!         1 atggcttaa
//! //
//! "#;
//! let records = parse_genbank::<Dna>(input)?;
//! let entry = &records[0];
//! assert_eq!(entry.record().id(), "DEMO.1");
//! let coding = entry.features()[0].extract(entry.record().sequence())?;
//! assert_eq!(coding.to_string(), "ATGGCTTAA");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::io::{BufRead, Cursor};
use std::iter::FusedIterator;
use std::marker::PhantomData;
use std::ops::Range;
use std::path::Path;

use align_core::Alphabet;

use super::{BioIoError, IoErrorKind, Lines, ReadLimits, open_input, read_text};
use crate::{DnaSequence, NucleotideSequenceExt, Record, RnaSequence, Sequence};

mod cds;
mod transform;
mod writer;

pub use cds::{
    CdsReferencePolicy, CdsReferenceStatus, CdsTranslation, CdsTranslationError,
    CdsTranslationOptions,
};
pub use transform::{
    AnnotatedTransform, AnnotationChange, AnnotationLossPolicy, AnnotationOperation,
    AnnotationTransformOptions, AnnotationTransformReport,
};
pub use writer::{write_embl, write_embl_new, write_genbank, write_genbank_new};

const MAX_LOCATION_DEPTH: usize = 64;

/// A physical metadata line, retained in input order, including unknown fields.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetadataField {
    /// Header field name. Continuations inherit the preceding field's name.
    pub key: String,
    /// Field payload without its fixed-column prefix or outer whitespace.
    pub value: String,
    /// Whether the `GenBank` field name was blank on this physical line.
    pub continuation: bool,
    /// Original physical line, excluding its CRLF or LF ending.
    pub raw_line: String,
}

/// One feature qualifier; repeated names remain separate and ordered.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Qualifier {
    /// Case-preserved qualifier name without its initial slash.
    pub name: String,
    /// Decoded value, or `None` for a flag such as `/pseudo`.
    ///
    /// Quotes and doubled-quote escaping are decoded. Wrapped physical payloads
    /// retain a newline; translation and other whitespace semantics are not guessed.
    pub value: Option<String>,
    /// Original qualifier payloads joined with newlines, including `/` and quotes.
    pub raw: String,
}

/// A location expression with explicitly supported extraction semantics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LocationExpression {
    /// Exact zero-based half-open range; a point is represented by length one.
    Range(Range<usize>),
    /// Concatenate child expressions in their written order.
    Join(Vec<Self>),
    /// Reverse-complement the entire child result, including joined segments.
    Complement(Box<Self>),
    /// Preserved expression whose sequence-extraction semantics are unsupported.
    Unsupported {
        /// Whitespace-normalized location text.
        expression: String,
        /// Explanation of the unsupported semantics.
        reason: String,
    },
}

/// Original location text plus a checked, partial interpretation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Location {
    raw: String,
    expression: LocationExpression,
}

impl Location {
    /// Parses exact points/ranges, `join` and `complement`, preserving others.
    ///
    /// # Errors
    /// Rejects empty or unbalanced expressions, malformed supported operators,
    /// zero/reversed/overflowing exact coordinates and nesting beyond 64 levels.
    /// Unsupported expressions are not certified as valid INSDC syntax.
    pub fn parse(text: &str) -> Result<Self, BioIoError> {
        let compact: String = text.chars().filter(|c| !c.is_ascii_whitespace()).collect();
        let expression = parse_location(&compact, 0)?;
        Ok(Self {
            raw: text.to_owned(),
            expression,
        })
    }

    /// Returns the input text, with newlines preserved across wrapped locations.
    #[must_use]
    pub fn raw(&self) -> &str {
        &self.raw
    }

    /// Returns the interpreted expression without discarding unsupported parts.
    #[must_use]
    pub const fn expression(&self) -> &LocationExpression {
        &self.expression
    }

    /// Extracts a feature from a typed sequence without changing its alphabet.
    ///
    /// # Errors
    /// Returns `Unsupported` for preserved non-exact semantics or protein
    /// complements, `LengthMismatch` for out-of-bounds ranges and `LimitExceeded`
    /// if repeated joins would exceed the default record-size bound.
    pub fn extract<A: Alphabet>(&self, sequence: &Sequence<A>) -> Result<Sequence<A>, BioIoError> {
        let mut output = Vec::new();
        extract_location(&self.expression, sequence, &mut output)?;
        Sequence::from_ascii(output)
            .map_err(|error| location_error(IoErrorKind::InvalidSequence, error.to_string()))
    }
}

/// An ordered feature with its location, repeated qualifiers and original lines.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Feature {
    /// Original feature key, including unknown keys.
    pub key: String,
    /// Location text and supported interpretation.
    pub location: Location,
    /// Qualifiers in source order; repeated names are not overwritten.
    pub qualifiers: Vec<Qualifier>,
    /// Original physical lines, excluding line endings.
    pub raw_lines: Vec<String>,
}

impl Feature {
    /// Extracts this feature using its location; returns nucleotides for a CDS.
    /// Use [`Self::translate_cds`] for annotation-aware protein translation.
    ///
    /// # Errors
    /// Propagates the explicit boundaries documented by [`Location::extract`].
    pub fn extract<A: Alphabet>(&self, sequence: &Sequence<A>) -> Result<Sequence<A>, BioIoError> {
        self.location.extract(sequence)
    }
}

/// A typed record with preserved INSDC metadata and features.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnotatedRecord<A: Alphabet> {
    format: Format,
    record: Record<A>,
    locus: String,
    molecule: String,
    accessions: Vec<String>,
    metadata: Vec<MetadataField>,
    features: Vec<Feature>,
    raw_record: String,
    layout: Vec<RecordPart>,
}

impl<A: Alphabet> AnnotatedRecord<A> {
    /// Returns the source/output profile. Cross-format writing is not implicit.
    #[must_use]
    pub const fn format(&self) -> AnnotatedFormat {
        self.format
    }

    /// Returns the typed sequence and derived identifier/description.
    #[must_use]
    pub const fn record(&self) -> &Record<A> {
        &self.record
    }
    /// Returns the original LOCUS or ID identifier, before accession preference.
    #[must_use]
    pub fn locus(&self) -> &str {
        &self.locus
    }
    /// Returns the declared molecule text; this does not select a genetic code.
    #[must_use]
    pub fn declared_molecule(&self) -> &str {
        &self.molecule
    }
    /// Returns all accessions in source order, including repetitions.
    #[must_use]
    pub fn accessions(&self) -> &[String] {
        &self.accessions
    }
    /// Returns every physical metadata line in input order.
    #[must_use]
    pub fn metadata(&self) -> &[MetadataField] {
        &self.metadata
    }
    /// Returns features in input order.
    #[must_use]
    pub fn features(&self) -> &[Feature] {
        &self.features
    }
    /// Returns source text, or regenerated text for a transformed record.
    ///
    /// Inter-record blank lines are not part of a record. This is retention of
    /// an original document. Writers render the model independently of this cache;
    /// checked transforms regenerate it along with feature and qualifier raw text.
    #[must_use]
    pub fn raw_record(&self) -> &str {
        &self.raw_record
    }
}

/// Parses a bounded, possibly multi-record `GenBank` document.
///
/// # Errors
/// Rejects malformed/truncated records, missing literal sequence, mismatched
/// lengths or molecule/alphabet, invalid symbols and default resource limits.
pub fn parse_genbank<A: Alphabet>(input: &str) -> Result<Vec<AnnotatedRecord<A>>, BioIoError> {
    parse_document(&read_text(Cursor::new(input), "genbank")?, Format::GenBank)
}

/// Reads a plain or gzip `GenBank` file into a bounded document.
///
/// # Errors
/// Returns I/O/decompression errors or the errors documented by [`parse_genbank`].
pub fn read_genbank<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<Vec<AnnotatedRecord<A>>, BioIoError> {
    parse_document(&read_text(open_input(path)?, "genbank")?, Format::GenBank)
}

/// Parses a bounded, possibly multi-record modern EMBL document.
///
/// # Errors
/// Rejects malformed/truncated records, missing literal sequence, mismatched
/// lengths or molecule/alphabet, invalid symbols and default resource limits.
pub fn parse_embl<A: Alphabet>(input: &str) -> Result<Vec<AnnotatedRecord<A>>, BioIoError> {
    parse_document(&read_text(Cursor::new(input), "embl")?, Format::Embl)
}

/// Reads a plain or gzip EMBL file into a bounded document.
///
/// # Errors
/// Returns I/O/decompression errors or the errors documented by [`parse_embl`].
pub fn read_embl<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<Vec<AnnotatedRecord<A>>, BioIoError> {
    parse_document(&read_text(open_input(path)?, "embl")?, Format::Embl)
}

/// Explicit annotated input/output profile; not a conversion instruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnnotatedFormat {
    /// Literal-sequence `GenBank`, including the existing protein profile.
    GenBank,
    /// Modern seven-field, nucleotide EMBL.
    Embl,
}

impl AnnotatedFormat {
    /// Returns the stable profile name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::GenBank => "genbank",
            Self::Embl => "embl",
        }
    }
    fn is_start(self, line: &str) -> bool {
        match self {
            Self::GenBank => line.starts_with("LOCUS       "),
            Self::Embl => line.starts_with("ID   "),
        }
    }
}

type Format = AnnotatedFormat;

/// What an [`AnnotatedReader`] does after a record fails to parse.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RecordErrorPolicy {
    /// Yield the error and end iteration, like the other streaming readers.
    #[default]
    Stop,
    /// Yield the error, then continue with the next record.
    ///
    /// This applies only to an error inside a complete record, from its
    /// `LOCUS`/`ID` line through its `//` terminator. Errors that leave the next
    /// record boundary unknown still end iteration: input outside a record, a
    /// missing terminator, an exceeded line or total byte limit, invalid UTF-8,
    /// and I/O or decompression failures.
    Continue,
}

/// A lazy `GenBank` or EMBL reader that keeps one record at a time.
///
/// Records are delimited exactly as [`parse_genbank`] and [`parse_embl`]
/// delimit a document, and each record is parsed by the same code, so every
/// record equals the corresponding element of their result, raw source text
/// included. Blank lines between records are skipped, and an empty input yields
/// no records. Errors carry one-based record and line numbers; a record's number
/// is its position in the input. After an error the iterator ends, unless
/// [`RecordErrorPolicy::Continue`] lets it move past a record that failed to
/// parse.
///
/// Memory holds the current record's source text and parsed values, not the
/// whole input.
pub struct AnnotatedReader<R, A: Alphabet> {
    lines: Lines<R>,
    format: Format,
    policy: RecordErrorPolicy,
    buffer: String,
    records: usize,
    finished: bool,
    alphabet: PhantomData<A>,
}

impl<R: BufRead, A: Alphabet> AnnotatedReader<R, A> {
    /// Constructs a reader with default resource limits; performs no reading.
    #[must_use]
    pub fn new(reader: R, format: AnnotatedFormat) -> Self {
        Self::with_limits(reader, format, ReadLimits::default())
    }

    /// Constructs a reader with explicit limits; performs no reading.
    ///
    /// `max_line_bytes` and `max_total_bytes` bound the input as it is read;
    /// `max_record_symbols` bounds each record's declared and observed sequence.
    #[must_use]
    pub fn with_limits(reader: R, format: AnnotatedFormat, limits: ReadLimits) -> Self {
        Self {
            lines: Lines::new(reader, limits, format.name()),
            format,
            policy: RecordErrorPolicy::Stop,
            buffer: String::new(),
            records: 0,
            finished: false,
            alphabet: PhantomData,
        }
    }

    /// Selects what happens after a record fails to parse; the default stops.
    #[must_use]
    pub const fn with_record_errors(mut self, policy: RecordErrorPolicy) -> Self {
        self.policy = policy;
        self
    }

    /// Reads the next complete record into the buffer and returns the number of
    /// its first line, or `None` at the end of the input.
    fn next_record_text(&mut self) -> Result<Option<usize>, BioIoError> {
        self.buffer.clear();
        let mut first_line = None;
        loop {
            // A failure belongs to the current record or, between records, the next.
            let record = self.records + usize::from(first_line.is_none());
            let line_number = self.lines.line + 1;
            let Some(physical) = self
                .lines
                .next_physical_line()
                .map_err(|error| error.at(record, line_number))?
            else {
                return match first_line {
                    None => Ok(None),
                    Some(line) => Err(BioIoError::new(
                        IoErrorKind::Truncated,
                        self.format.name(),
                        "record is missing its // terminator",
                    )
                    .at(record, line)),
                };
            };
            let line = physical.trim_end_matches(['\r', '\n']);
            if first_line.is_none() {
                if line.trim().is_empty() {
                    continue;
                }
                if !self.format.is_start(line) {
                    return Err(BioIoError::new(
                        IoErrorKind::InvalidHeader,
                        self.format.name(),
                        "expected a LOCUS/ID record header",
                    )
                    .at(record, line_number));
                }
                self.records += 1;
                first_line = Some(line_number);
            } else if self.format.is_start(line) {
                return Err(BioIoError::new(
                    IoErrorKind::Truncated,
                    self.format.name(),
                    "record is missing its // terminator",
                )
                .at(record, line_number));
            }
            let terminated = line == "//";
            self.buffer.push_str(&physical);
            if terminated {
                return Ok(first_line);
            }
        }
    }
}

impl<R: BufRead, A: Alphabet> Iterator for AnnotatedReader<R, A> {
    type Item = Result<AnnotatedRecord<A>, BioIoError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let first_line = match self.next_record_text() {
            Ok(Some(first_line)) => first_line,
            Ok(None) => {
                self.finished = true;
                return None;
            }
            Err(error) => {
                self.finished = true;
                return Some(Err(error));
            }
        };
        let parsed = parse_record(
            &self.buffer,
            self.format,
            self.records,
            first_line,
            self.lines.limits,
        );
        if parsed.is_err() && self.policy == RecordErrorPolicy::Stop {
            self.finished = true;
        }
        Some(parsed)
    }
}

impl<R: BufRead, A: Alphabet> FusedIterator for AnnotatedReader<R, A> {}

/// Opens a plain or gzip `GenBank` file as a lazy record iterator.
///
/// # Errors
/// Returns opening errors immediately; parse and decompression errors are items.
pub fn genbank<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<AnnotatedReader<Box<dyn BufRead>, A>, BioIoError> {
    Ok(AnnotatedReader::new(
        open_input(path)?,
        AnnotatedFormat::GenBank,
    ))
}

/// Opens a plain or gzip EMBL file as a lazy record iterator.
///
/// # Errors
/// Returns opening errors immediately; parse and decompression errors are items.
pub fn embl<A: Alphabet>(
    path: impl AsRef<Path>,
) -> Result<AnnotatedReader<Box<dyn BufRead>, A>, BioIoError> {
    Ok(AnnotatedReader::new(
        open_input(path)?,
        AnnotatedFormat::Embl,
    ))
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum RecordPart {
    Metadata(usize),
    Feature(usize),
}

fn location_error(kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, "insdc-location", message)
}

fn parse_location(text: &str, depth: usize) -> Result<LocationExpression, BioIoError> {
    if depth > MAX_LOCATION_DEPTH {
        return Err(location_error(
            IoErrorKind::LimitExceeded,
            "location nesting exceeds 64 levels",
        ));
    }
    if text.is_empty() || !text.is_ascii() {
        return Err(location_error(
            IoErrorKind::InvalidStructure,
            "location must be nonempty ASCII",
        ));
    }
    let mut balance = 0_usize;
    for byte in text.bytes() {
        if byte == b'(' {
            balance += 1;
        }
        if byte == b')' {
            balance = balance.checked_sub(1).ok_or_else(|| {
                location_error(
                    IoErrorKind::InvalidStructure,
                    "unbalanced location parentheses",
                )
            })?;
        }
    }
    if balance != 0 {
        return Err(location_error(
            IoErrorKind::InvalidStructure,
            "unbalanced location parentheses",
        ));
    }
    for operator in ["complement", "join"] {
        if let Some(body) = text.strip_prefix(&format!("{operator}(")) {
            let body = body.strip_suffix(')').ok_or_else(|| {
                location_error(
                    IoErrorKind::InvalidStructure,
                    "location operator lacks closing parenthesis",
                )
            })?;
            let parts = split_location_parts(body)?;
            if operator == "complement" {
                if parts.len() != 1 {
                    return Err(location_error(
                        IoErrorKind::InvalidStructure,
                        "complement requires one location",
                    ));
                }
                return Ok(LocationExpression::Complement(Box::new(parse_location(
                    parts[0],
                    depth + 1,
                )?)));
            }
            if parts.len() < 2 {
                return Err(location_error(
                    IoErrorKind::InvalidStructure,
                    "join requires at least two locations",
                ));
            }
            return Ok(LocationExpression::Join(
                parts
                    .into_iter()
                    .map(|part| parse_location(part, depth + 1))
                    .collect::<Result<_, _>>()?,
            ));
        }
    }
    if text.bytes().all(|b| b.is_ascii_digit() || b == b'.') {
        let (start, end) = if let Some((start, end)) = text.split_once("..") {
            (coordinate(start)?, coordinate(end)?)
        } else if text.contains('.') {
            return Ok(unsupported_location(text, "uncertain single-base location"));
        } else {
            let point = coordinate(text)?;
            (point, point)
        };
        if start > end {
            return Err(location_error(
                IoErrorKind::InvalidStructure,
                "location start exceeds its end",
            ));
        }
        return Ok(LocationExpression::Range(start - 1..end));
    }
    let reason = if text.contains(':') {
        "remote references require an explicit resolver"
    } else if text.contains(['<', '>']) {
        "fuzzy boundaries cannot be extracted exactly"
    } else if text.contains('^') {
        "between-base locations are not sequence spans"
    } else if text.starts_with("order(") {
        "order does not imply sequence concatenation"
    } else {
        "unrecognized location expression"
    };
    Ok(unsupported_location(text, reason))
}

fn unsupported_location(text: &str, reason: &str) -> LocationExpression {
    LocationExpression::Unsupported {
        expression: text.to_owned(),
        reason: reason.to_owned(),
    }
}

fn coordinate(text: &str) -> Result<usize, BioIoError> {
    let value = text.parse::<usize>().map_err(|_| {
        location_error(
            IoErrorKind::InvalidStructure,
            "invalid or overflowing exact coordinate",
        )
    })?;
    if value == 0 {
        return Err(location_error(
            IoErrorKind::InvalidStructure,
            "INSDC coordinates begin at one",
        ));
    }
    Ok(value)
}

fn split_location_parts(text: &str) -> Result<Vec<&str>, BioIoError> {
    let mut depth = 0_usize;
    let mut start = 0;
    let mut parts = Vec::new();
    for (index, byte) in text.bytes().enumerate() {
        match byte {
            b'(' => depth += 1,
            b')' => {
                depth = depth.checked_sub(1).ok_or_else(|| {
                    location_error(
                        IoErrorKind::InvalidStructure,
                        "unexpected closing parenthesis",
                    )
                })?;
            }
            b',' if depth == 0 => {
                parts.push(&text[start..index]);
                start = index + 1;
            }
            _ => {}
        }
    }
    parts.push(&text[start..]);
    if depth != 0 || parts.iter().any(|part| part.is_empty()) {
        return Err(location_error(
            IoErrorKind::InvalidStructure,
            "empty or unbalanced location operand",
        ));
    }
    Ok(parts)
}

fn check_range(range: &Range<usize>, length: usize) -> Result<(), BioIoError> {
    if range.start > range.end || range.end > length {
        return Err(location_error(
            IoErrorKind::LengthMismatch,
            format!(
                "feature range {}..{} exceeds sequence length {length}",
                range.start, range.end
            ),
        ));
    }
    Ok(())
}

fn validate_locations(expression: &LocationExpression, length: usize) -> Result<(), BioIoError> {
    match expression {
        LocationExpression::Range(range) => check_range(range, length),
        LocationExpression::Join(parts) => parts
            .iter()
            .try_for_each(|part| validate_locations(part, length)),
        LocationExpression::Complement(inner) => validate_locations(inner, length),
        LocationExpression::Unsupported { .. } => Ok(()),
    }
}

fn extract_location<A: Alphabet>(
    expression: &LocationExpression,
    sequence: &Sequence<A>,
    output: &mut Vec<u8>,
) -> Result<(), BioIoError> {
    match expression {
        LocationExpression::Range(range) => {
            check_range(range, sequence.len())?;
            let length = range.end - range.start;
            if length
                > ReadLimits::default()
                    .max_record_symbols
                    .saturating_sub(output.len())
            {
                return Err(location_error(
                    IoErrorKind::LimitExceeded,
                    "extracted feature exceeds the record symbol limit",
                ));
            }
            output
                .try_reserve(length)
                .map_err(|error| location_error(IoErrorKind::LimitExceeded, error.to_string()))?;
            output.extend(
                sequence.as_slice()[range.clone()]
                    .iter()
                    .map(|symbol| symbol.as_ascii()),
            );
            Ok(())
        }
        LocationExpression::Join(parts) => parts
            .iter()
            .try_for_each(|part| extract_location(part, sequence, output)),
        LocationExpression::Complement(inner) => {
            let start = output.len();
            extract_location(inner, sequence, output)?;
            let reversed = match A::NAME {
                "DNA" => DnaSequence::from_ascii(&output[start..])
                    .and_then(|value| value.reverse_complement())
                    .map(|value| value.to_ascii_string()),
                "RNA" => RnaSequence::from_ascii(&output[start..])
                    .and_then(|value| value.reverse_complement())
                    .map(|value| value.to_ascii_string()),
                _ => {
                    return Err(location_error(
                        IoErrorKind::Unsupported,
                        "complement extraction requires DNA or RNA",
                    ));
                }
            }
            .map_err(|error| location_error(IoErrorKind::InvalidSequence, error.to_string()))?;
            output[start..].copy_from_slice(reversed.as_bytes());
            Ok(())
        }
        LocationExpression::Unsupported { reason, .. } => {
            Err(location_error(IoErrorKind::Unsupported, reason))
        }
    }
}

fn parse_document<A: Alphabet>(
    text: &str,
    format: Format,
) -> Result<Vec<AnnotatedRecord<A>>, BioIoError> {
    let mut records = Vec::new();
    let mut start = None;
    let mut offset = 0;
    for (index, physical) in text.split_inclusive('\n').enumerate() {
        let line = physical.trim_end_matches(['\r', '\n']);
        let record_number = records.len() + 1;
        if physical.len() > ReadLimits::default().max_line_bytes {
            return Err(BioIoError::new(
                IoErrorKind::LimitExceeded,
                format.name(),
                "physical line byte limit exceeded",
            )
            .at(record_number, index + 1));
        }
        if start.is_none() && !line.trim().is_empty() {
            if !format.is_start(line) {
                return Err(BioIoError::new(
                    IoErrorKind::InvalidHeader,
                    format.name(),
                    "expected a LOCUS/ID record header",
                )
                .at(record_number, index + 1));
            }
            start = Some((offset, index + 1));
        } else if start.is_some() && format.is_start(line) {
            return Err(BioIoError::new(
                IoErrorKind::Truncated,
                format.name(),
                "record is missing its // terminator",
            )
            .at(record_number, index + 1));
        }
        offset += physical.len();
        if line == "//" {
            let Some((start_offset, start_line)) = start.take() else {
                return Err(BioIoError::new(
                    IoErrorKind::InvalidStructure,
                    format.name(),
                    "terminator outside a record",
                )
                .at(record_number, index + 1));
            };
            records.push(parse_record(
                &text[start_offset..offset],
                format,
                record_number,
                start_line,
                ReadLimits::default(),
            )?);
        }
    }
    if let Some((_, line)) = start {
        return Err(BioIoError::new(
            IoErrorKind::Truncated,
            format.name(),
            "record is missing its // terminator",
        )
        .at(records.len() + 1, line));
    }
    Ok(records)
}

struct Header {
    locus: String,
    length: usize,
    molecule: String,
    version: Option<String>,
}

fn parse_header<A: Alphabet>(
    line: &str,
    format: Format,
    limits: ReadLimits,
) -> Result<Header, BioIoError> {
    let invalid = || {
        BioIoError::new(
            IoErrorKind::InvalidHeader,
            format.name(),
            "invalid LOCUS/ID name, length, unit or molecule",
        )
    };
    let header = match format {
        Format::GenBank => {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() < 4 {
                return Err(invalid());
            }
            let length = fields[2].parse().map_err(|_| invalid())?;
            let molecule = match fields.get(3) {
                Some(&"aa") => "protein",
                Some(&"bp") => {
                    let molecule = *fields.get(4).ok_or_else(invalid)?;
                    if molecule.eq_ignore_ascii_case("protein") {
                        return Err(invalid());
                    }
                    molecule
                }
                _ => return Err(invalid()),
            };
            Header {
                locus: fields[1].to_owned(),
                length,
                molecule: molecule.to_owned(),
                version: None,
            }
        }
        Format::Embl => {
            let fields: Vec<_> = line[5..].split(';').map(str::trim).collect();
            if fields.len() != 7 {
                return Err(BioIoError::new(
                    IoErrorKind::Unsupported,
                    format.name(),
                    "only modern seven-field EMBL ID lines are supported",
                ));
            }
            let length_unit: Vec<_> = fields[6].split_whitespace().collect();
            if fields[0].is_empty()
                || fields[0].contains(char::is_whitespace)
                || length_unit.len() != 2
                || length_unit[1] != "BP."
            {
                return Err(invalid());
            }
            let version = fields[1].strip_prefix("SV ").ok_or_else(invalid)?;
            if version.is_empty() || !version.bytes().all(|b| b.is_ascii_digit()) {
                return Err(invalid());
            }
            Header {
                locus: fields[0].to_owned(),
                length: length_unit[0].parse().map_err(|_| invalid())?,
                molecule: fields[3].to_owned(),
                version: Some(format!("{}.{version}", fields[0])),
            }
        }
    };
    if header.length > limits.max_record_symbols {
        return Err(BioIoError::new(
            IoErrorKind::LimitExceeded,
            format.name(),
            "declared sequence length exceeds the record symbol limit",
        ));
    }
    let expected = molecule_alphabet(&header.molecule, format)?;
    if expected != A::NAME {
        return Err(BioIoError::new(
            IoErrorKind::InvalidSequence,
            format.name(),
            format!(
                "declared molecule {} does not match requested {} alphabet",
                header.molecule,
                A::NAME
            ),
        ));
    }
    Ok(header)
}

fn molecule_alphabet(molecule: &str, format: Format) -> Result<&'static str, BioIoError> {
    let molecule_upper = molecule.to_ascii_uppercase();
    let molecule_base = ["SS-", "DS-", "MS-"]
        .iter()
        .find_map(|prefix| molecule_upper.strip_prefix(prefix))
        .unwrap_or(&molecule_upper);
    let expected = match molecule_base {
        "DNA" => "DNA",
        "GENOMIC DNA" if matches!(format, Format::Embl) => "DNA",
        "OTHER DNA" if matches!(format, Format::Embl) && molecule_upper == "OTHER DNA" => "DNA",
        "GENOMIC RNA" if matches!(format, Format::Embl) => "RNA",
        "RNA" | "MRNA" | "TRNA" | "RRNA" | "URNA" | "SNRNA" | "SNORNA" | "SCRNA" | "PRE-RNA"
        | "NCRNA" | "TMRNA" | "CRNA" => "RNA",
        "PROTEIN" if molecule == "protein" && matches!(format, Format::GenBank) => "protein",
        _ => {
            return Err(BioIoError::new(
                IoErrorKind::Unsupported,
                format.name(),
                "unsupported declared molecule",
            ));
        }
    };
    Ok(expected)
}

struct RecordBuilder {
    header: Header,
    metadata: Vec<MetadataField>,
    features: Vec<Feature>,
    feature: Option<FeatureBuilder>,
    sequence: String,
    in_features: bool,
    in_sequence: bool,
    last_key: String,
    layout: Vec<RecordPart>,
}

fn parse_record<A: Alphabet>(
    text: &str,
    format: Format,
    number: usize,
    first_line: usize,
    limits: ReadLimits,
) -> Result<AnnotatedRecord<A>, BioIoError> {
    let mut lines = text.lines();
    let first = lines.next().ok_or_else(|| {
        BioIoError::new(IoErrorKind::InvalidHeader, format.name(), "empty record")
    })?;
    let header = parse_header::<A>(first, format, limits).map_err(|e| e.at(number, first_line))?;
    let mut builder = RecordBuilder {
        header,
        metadata: Vec::new(),
        features: Vec::new(),
        feature: None,
        sequence: String::new(),
        in_features: false,
        in_sequence: false,
        last_key: String::new(),
        layout: Vec::new(),
    };
    builder.push_metadata(first, format)?;
    for (offset, line) in lines.enumerate() {
        let context = |error: BioIoError| {
            error
                .in_format(format.name())
                .at(number, first_line + offset + 1)
        };
        if line == "//" {
            break;
        }
        if builder.in_sequence {
            parse_sequence_line::<A>(line, format, limits, &mut builder.sequence)
                .map_err(context)?;
        } else {
            parse_annotation_line(line, format, &mut builder).map_err(context)?;
        }
    }
    finish_record(builder, text, format)
        .map_err(|error| error.in_format(format.name()).at(number, first_line))
}

fn metadata_line(
    line: &str,
    format: Format,
    previous: &mut String,
) -> Result<MetadataField, BioIoError> {
    let width = match format {
        Format::GenBank => 12,
        Format::Embl => 5,
    };
    let (key, payload) = if line.len() >= width {
        (
            line.get(..width).ok_or_else(|| {
                BioIoError::new(
                    IoErrorKind::InvalidHeader,
                    format.name(),
                    "non-ASCII field prefix",
                )
            })?,
            &line[width..],
        )
    } else {
        (line, "")
    };
    let key = key.trim();
    let continuation = key.is_empty();
    if continuation && previous.is_empty() {
        return Err(BioIoError::new(
            IoErrorKind::InvalidHeader,
            format.name(),
            "metadata continuation without a field",
        ));
    }
    if !continuation {
        key.clone_into(previous);
    }
    Ok(MetadataField {
        key: previous.clone(),
        value: payload.trim().to_owned(),
        continuation,
        raw_line: line.to_owned(),
    })
}

fn parse_annotation_line(
    line: &str,
    format: Format,
    builder: &mut RecordBuilder,
) -> Result<(), BioIoError> {
    let sequence_start = match format {
        Format::GenBank => line == "ORIGIN" || line.starts_with("ORIGIN      "),
        Format::Embl => line.starts_with("SQ   "),
    };
    if sequence_start {
        builder.finish_feature()?;
        if let Format::Embl = format {
            let words: Vec<_> = line[5..].split_whitespace().collect();
            if words.len() < 3 || words[0] != "Sequence" || words[2] != "BP;" {
                return Err(BioIoError::new(
                    IoErrorKind::InvalidHeader,
                    format.name(),
                    "invalid SQ sequence header",
                ));
            }
            if words[1].parse::<usize>().ok() != Some(builder.header.length) {
                return Err(BioIoError::new(
                    IoErrorKind::LengthMismatch,
                    format.name(),
                    "SQ length differs from ID length",
                ));
            }
        }
        builder.push_metadata(line, format)?;
        builder.in_sequence = true;
        return Ok(());
    }
    let feature_start = match format {
        Format::GenBank => line.starts_with("FEATURES"),
        Format::Embl => line.starts_with("FH"),
    };
    if feature_start {
        builder.in_features = true;
        builder.push_metadata(line, format)?;
        return Ok(());
    }
    let is_feature = match format {
        Format::GenBank => builder.in_features && line.starts_with("     "),
        Format::Embl => line.starts_with("FT   "),
    };
    if is_feature {
        return builder.feature_line(line, format);
    }
    builder.finish_feature()?;
    if !line.trim().is_empty() {
        builder.push_metadata(line, format)?;
    }
    Ok(())
}

impl RecordBuilder {
    fn push_metadata(&mut self, line: &str, format: Format) -> Result<(), BioIoError> {
        let field = metadata_line(line, format, &mut self.last_key)?;
        self.layout.push(RecordPart::Metadata(self.metadata.len()));
        self.metadata.push(field);
        Ok(())
    }

    fn finish_feature(&mut self) -> Result<(), BioIoError> {
        if let Some(feature) = self.feature.take() {
            self.features.push(feature.finish()?);
        }
        Ok(())
    }

    fn feature_line(&mut self, line: &str, format: Format) -> Result<(), BioIoError> {
        let invalid = || {
            BioIoError::new(
                IoErrorKind::InvalidStructure,
                format.name(),
                "invalid fixed-column feature line",
            )
        };
        if line.len() < 21 || !line.as_bytes()[..21].is_ascii() || line.as_bytes()[20] != b' ' {
            return Err(invalid());
        }
        let key = line[5..20].trim();
        let payload = line[21..].trim();
        if !key.is_empty() {
            if !is_component_name(key, 15) || payload.is_empty() {
                return Err(invalid());
            }
            self.finish_feature()?;
            self.layout.push(RecordPart::Feature(self.features.len()));
            self.feature = Some(FeatureBuilder {
                key: key.to_owned(),
                location: payload.to_owned(),
                qualifiers: Vec::new(),
                qualifier: None,
                raw_lines: vec![line.to_owned()],
            });
            return Ok(());
        }
        let feature = self.feature.as_mut().ok_or_else(invalid)?;
        feature.raw_lines.push(line.to_owned());
        feature.continue_payload(payload)
    }
}

struct FeatureBuilder {
    key: String,
    location: String,
    qualifiers: Vec<Qualifier>,
    qualifier: Option<String>,
    raw_lines: Vec<String>,
}

impl FeatureBuilder {
    fn continue_payload(&mut self, payload: &str) -> Result<(), BioIoError> {
        if let Some(raw) = &mut self.qualifier {
            if quoted_value_is_open(raw)? {
                raw.push('\n');
                raw.push_str(payload);
                return Ok(());
            }
        }
        if payload.starts_with('/') {
            self.finish_qualifier()?;
            self.qualifier = Some(payload.to_owned());
        } else if let Some(raw) = &mut self.qualifier {
            if raw.contains('"') {
                return Err(location_error(
                    IoErrorKind::InvalidStructure,
                    "text after a closed quoted qualifier",
                ));
            }
            raw.push('\n');
            raw.push_str(payload);
        } else {
            self.location.push('\n');
            self.location.push_str(payload);
        }
        Ok(())
    }

    fn finish_qualifier(&mut self) -> Result<(), BioIoError> {
        if let Some(raw) = self.qualifier.take() {
            self.qualifiers.push(parse_qualifier(raw)?);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Feature, BioIoError> {
        self.finish_qualifier()?;
        Ok(Feature {
            key: self.key,
            location: Location::parse(&self.location)?,
            qualifiers: self.qualifiers,
            raw_lines: self.raw_lines,
        })
    }
}

fn quoted_value_is_open(raw: &str) -> Result<bool, BioIoError> {
    let Some((_, value)) = raw.split_once('=') else {
        return Ok(false);
    };
    let Some(value) = value.strip_prefix('"') else {
        return Ok(false);
    };
    let mut bytes = value.bytes().peekable();
    while let Some(byte) = bytes.next() {
        if byte == b'"' {
            if bytes.peek() == Some(&b'"') {
                bytes.next();
            } else if bytes.any(|b| !b.is_ascii_whitespace()) {
                return Err(location_error(
                    IoErrorKind::InvalidStructure,
                    "unexpected text after qualifier closing quote",
                ));
            } else {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn parse_qualifier(raw: String) -> Result<Qualifier, BioIoError> {
    let payload = raw
        .strip_prefix('/')
        .ok_or_else(|| location_error(IoErrorKind::InvalidStructure, "qualifier lacks slash"))?;
    let (name, value) = payload
        .split_once('=')
        .map_or((payload, None), |(key, value)| (key, Some(value)));
    if !is_component_name(name, 20) {
        return Err(location_error(
            IoErrorKind::InvalidStructure,
            "invalid qualifier name",
        ));
    }
    let value = match value {
        Some(value) if value.starts_with('"') => {
            if quoted_value_is_open(&raw)? {
                return Err(location_error(
                    IoErrorKind::Truncated,
                    "unterminated quoted qualifier",
                ));
            }
            let inside = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .ok_or_else(|| {
                    location_error(IoErrorKind::InvalidStructure, "invalid quoted qualifier")
                })?;
            Some(inside.replace("\"\"", "\""))
        }
        Some(value) if value.is_empty() || value.contains('"') => {
            return Err(location_error(
                IoErrorKind::InvalidStructure,
                "invalid unquoted qualifier value",
            ));
        }
        Some(value) => Some(value.to_owned()),
        None => None,
    };
    Ok(Qualifier {
        name: name.to_owned(),
        value,
        raw,
    })
}

fn is_component_name(name: &str, maximum: usize) -> bool {
    name.len() <= maximum
        && name.bytes().any(|byte| byte.is_ascii_alphabetic())
        && name
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'\'' | b'*'))
}

fn parse_sequence_line<A: Alphabet>(
    line: &str,
    format: Format,
    limits: ReadLimits,
    sequence: &mut String,
) -> Result<(), BioIoError> {
    if line.trim().is_empty() {
        return Ok(());
    }
    let mut fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() < 2 {
        return Err(BioIoError::new(
            IoErrorKind::InvalidSequence,
            format.name(),
            "sequence line requires symbols and a position counter",
        ));
    }
    let position = match format {
        Format::GenBank => fields.remove(0),
        Format::Embl => fields.pop().unwrap_or_default(),
    };
    let count = position.parse::<usize>().map_err(|_| {
        BioIoError::new(
            IoErrorKind::InvalidSequence,
            format.name(),
            "invalid sequence position counter",
        )
    })?;
    let additional: usize = fields.iter().map(|field| field.len()).sum();
    if additional > limits.max_record_symbols.saturating_sub(sequence.len()) {
        return Err(BioIoError::new(
            IoErrorKind::LimitExceeded,
            format.name(),
            "sequence exceeds the record symbol limit",
        ));
    }
    let expected = match format {
        Format::GenBank => sequence.len() + 1,
        Format::Embl => sequence.len() + additional,
    };
    if count != expected {
        return Err(BioIoError::new(
            IoErrorKind::LengthMismatch,
            format.name(),
            "sequence position counter is inconsistent",
        ));
    }
    for field in fields {
        if !field
            .bytes()
            .all(|byte| byte.is_ascii_alphabetic() && A::index_of(byte).is_some())
        {
            return Err(BioIoError::new(
                IoErrorKind::InvalidSequence,
                format.name(),
                format!(
                    "sequence row contains a residue outside the {} alphabet",
                    A::NAME
                ),
            ));
        }
        sequence.push_str(field);
    }
    Ok(())
}

fn finish_record<A: Alphabet>(
    mut builder: RecordBuilder,
    raw: &str,
    format: Format,
) -> Result<AnnotatedRecord<A>, BioIoError> {
    builder.finish_feature()?;
    if !builder.in_sequence {
        return Err(BioIoError::new(
            IoErrorKind::Unsupported,
            format.name(),
            "literal ORIGIN/SQ sequence is required; CONTIG-only records are unsupported",
        ));
    }
    if builder.sequence.len() != builder.header.length {
        return Err(BioIoError::new(
            IoErrorKind::LengthMismatch,
            format.name(),
            format!(
                "declared {} symbols, observed {}",
                builder.header.length,
                builder.sequence.len()
            ),
        ));
    }
    let mut description = Vec::new();
    let mut accessions = Vec::new();
    let mut version = builder.header.version;
    for field in &builder.metadata {
        match field.key.as_str() {
            "DEFINITION" | "DE" => description.push(field.value.as_str()),
            "ACCESSION" => accessions.extend(field.value.split_whitespace().map(str::to_owned)),
            "AC" => accessions.extend(
                field
                    .value
                    .split(';')
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(str::to_owned),
            ),
            "VERSION" | "SV" => {
                let value = field.value.split_whitespace().next().ok_or_else(|| {
                    BioIoError::new(
                        IoErrorKind::InvalidHeader,
                        format.name(),
                        "empty version field",
                    )
                })?;
                if version.as_ref().is_some_and(|previous| previous != value) {
                    return Err(BioIoError::new(
                        IoErrorKind::InvalidHeader,
                        format.name(),
                        "conflicting version identifiers",
                    ));
                }
                version = Some(value.to_owned());
            }
            _ => {}
        }
    }
    let id = version
        .or_else(|| accessions.first().cloned())
        .unwrap_or_else(|| builder.header.locus.clone());
    let description = description.join(" ");
    let header = if description.is_empty() {
        id
    } else {
        format!("{id} {description}")
    };
    let record = Record::from_ascii(header, builder.sequence, None)
        .map_err(|error| error.in_format(format.name()))?;
    for feature in &builder.features {
        validate_locations(feature.location.expression(), record.sequence().len())?;
    }
    Ok(AnnotatedRecord {
        format,
        record,
        locus: builder.header.locus,
        molecule: builder.header.molecule,
        accessions,
        metadata: builder.metadata,
        features: builder.features,
        raw_record: raw.to_owned(),
        layout: builder.layout,
    })
}
