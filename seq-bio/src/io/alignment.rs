//! Conservative alignment-file import without silently losing biological data.
//!
//! Imported documents preserve missing data, all-gap columns and annotations.
//! Conversion to the existing computational alignment type is checked separately.
//!
//! ```
//! use bioferris::{Dna, io::alignment::{AlignmentFormat, parse_alignment}};
//!
//! let input = b">reference\nAC-G\n>sample\nACTG\n";
//! let document = parse_alignment::<Dna>(&input[..], AlignmentFormat::Fasta)?;
//! println!("{} rows, {} columns", document.rows().len(), document.width());
//! let alignment = document.to_multiple_alignment()?;
//! assert_eq!(alignment.rows()[0].to_ascii_with_gaps(), "AC-G");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::marker::PhantomData;
use std::path::Path;

use align::prelude::{AlignedRow, MultipleAlignment, PhylipLayout, PhylipNameStyle};
use align_core::{Alphabet, Symbol};

use super::{BioIoError, IoErrorKind, open_input, read_text};

/// An explicitly selected alignment format; no filename/content guessing occurs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentFormat {
    /// Rectangular FASTA, preserving the complete trimmed header as the row name.
    Fasta,
    /// CLUSTAL blocks, with optional consensus lines and checked residue counts.
    Clustal,
    /// PHYLIP with both name convention and layout selected explicitly.
    Phylip {
        /// Fixed ten-byte ASCII names or whitespace-delimited relaxed names.
        names: PhylipNameStyle,
        /// Sequential records or equal-width interleaved blocks.
        layout: PhylipLayout,
    },
    /// One Stockholm 1.0 document, including GF, GS, GC and GR annotations.
    Stockholm,
    /// One NEXUS DATA block, with a noninterleaved DNA/RNA/protein matrix.
    Nexus,
}

impl AlignmentFormat {
    /// Returns the format label used in structured diagnostics.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Fasta => "aligned FASTA",
            Self::Clustal => "CLUSTAL",
            Self::Phylip { .. } => "PHYLIP",
            Self::Stockholm => "Stockholm 1.0",
            Self::Nexus => "NEXUS",
        }
    }
}

/// Scope of a preserved annotation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AnnotationScope {
    /// Free-text file annotation (Stockholm GF).
    File,
    /// Free-text sequence annotation (Stockholm GS).
    Sequence,
    /// One character per alignment column (GC or CLUSTAL consensus).
    Column,
    /// One character per column of a named sequence (Stockholm GR).
    Residue,
    /// Source comment, retained as text rather than interpreted as executable data.
    Comment,
}

/// One annotation; repeated free-text tags remain separate and ordered.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentAnnotation {
    scope: AnnotationScope,
    sequence: Option<String>,
    tag: String,
    value: String,
    line: usize,
}

impl AlignmentAnnotation {
    /// Returns whether this annotation describes a file, row, column or residue.
    #[must_use]
    pub const fn scope(&self) -> AnnotationScope {
        self.scope
    }

    /// Returns the target sequence name for GS/GR, otherwise `None`.
    #[must_use]
    pub fn sequence(&self) -> Option<&str> {
        self.sequence.as_deref()
    }

    /// Returns the original annotation tag (or `comment` / `consensus`).
    #[must_use]
    pub fn tag(&self) -> &str {
        &self.tag
    }

    /// Returns the text, or concatenated column annotation for wrapped input.
    #[must_use]
    pub fn value(&self) -> &str {
        &self.value
    }

    /// Returns the one-based source line of the first fragment.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

/// One named, rectangular row in an imported document.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentRow {
    identifier: String,
    aligned: String,
    line: usize,
}

impl AlignmentRow {
    /// Returns the name, in input row order.
    #[must_use]
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    /// Returns uppercase residues with the original gap/missing-data glyphs.
    #[must_use]
    pub fn aligned(&self) -> &str {
        &self.aligned
    }

    /// Returns the one-based source line where the row first appeared.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

/// A validated document, deliberately broader than [`MultipleAlignment`].
///
/// All-gap columns and explicit missing data are valid imported data. No parser
/// deletes columns or replaces missing data with gaps or alphabet ambiguity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentDocument<A: Alphabet> {
    format: AlignmentFormat,
    rows: Vec<AlignmentRow>,
    annotations: Vec<AlignmentAnnotation>,
    width: usize,
    gaps: Vec<u8>,
    missing: Option<u8>,
    marker: PhantomData<A>,
}

impl<A: Alphabet> AlignmentDocument<A> {
    /// Returns rows in their stable input order.
    #[must_use]
    pub fn rows(&self) -> &[AlignmentRow] {
        &self.rows
    }

    /// Returns the number of columns, including gaps and missing cells.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Returns preserved annotations in first-occurrence order.
    #[must_use]
    pub fn annotations(&self) -> &[AlignmentAnnotation] {
        &self.annotations
    }

    /// Returns the explicitly selected format/profile.
    #[must_use]
    pub const fn format(&self) -> AlignmentFormat {
        self.format
    }

    /// Returns the ASCII glyphs representing gaps in this document.
    #[must_use]
    pub fn gap_symbols(&self) -> &[u8] {
        &self.gaps
    }

    /// Returns the ASCII missing-data glyph, if this profile has one.
    #[must_use]
    pub const fn missing_symbol(&self) -> Option<u8> {
        self.missing
    }

    /// Converts without discarding annotations or changing biological cells.
    ///
    /// # Errors
    ///
    /// Returns `DataLoss` for annotations or missing data, and `InvalidStructure`
    /// for kernel-incompatible data such as an all-gap column. Dot gaps, where
    /// supported, become the computational type's gap cells, not residues.
    pub fn to_multiple_alignment(&self) -> Result<MultipleAlignment<A>, BioIoError> {
        if !self.annotations.is_empty() {
            return Err(failure(
                self.format,
                IoErrorKind::DataLoss,
                "the computational alignment cannot retain annotations; explicitly use to_multiple_alignment_discarding_annotations if intended",
            ));
        }
        self.to_multiple_alignment_discarding_annotations()
    }

    /// Extracts an alignment with explicit permission to omit annotations.
    ///
    /// The source document is borrowed and remains available. This permission
    /// does not authorize changing missing data, deleting columns or repairing
    /// otherwise invalid computational alignments.
    ///
    /// # Errors
    ///
    /// Returns `DataLoss` for missing data or `InvalidStructure` for unsupported
    /// computational structure (including all-gap columns).
    pub fn to_multiple_alignment_discarding_annotations(
        &self,
    ) -> Result<MultipleAlignment<A>, BioIoError> {
        let mut rows = Vec::new();
        for (index, row) in self.rows.iter().enumerate() {
            let mut cells = Vec::new();
            for (column, byte) in row.aligned.bytes().enumerate() {
                if Some(byte) == self.missing {
                    return Err(failure(
                        self.format,
                        IoErrorKind::DataLoss,
                        format!("missing data at zero-based column {column} cannot be represented by the computational alignment"),
                    ).at(index + 1, row.line));
                }
                if self.gaps.contains(&byte) {
                    cells.push(None);
                } else {
                    let symbol = Symbol::<A>::from_ascii(byte).ok_or_else(|| {
                        failure(
                            self.format,
                            IoErrorKind::InvalidSequence,
                            "invalid stored symbol",
                        )
                    })?;
                    cells.push(Some(symbol));
                }
            }
            rows.push(
                AlignedRow::new(row.identifier.clone(), cells).map_err(|error| {
                    failure(
                        self.format,
                        IoErrorKind::InvalidStructure,
                        error.to_string(),
                    )
                    .at(index + 1, row.line)
                })?,
            );
        }
        MultipleAlignment::try_new(rows).map_err(|error| {
            failure(
                self.format,
                IoErrorKind::InvalidStructure,
                error.to_string(),
            )
        })
    }
}

/// Parses one complete alignment document using the common bounded text reader.
///
/// This is an eager alignment reader, not a constant-memory record iterator.
/// The alphabet and PHYLIP dialect must be chosen explicitly. Trailing documents
/// and unsupported constructs are rejected rather than silently skipped.
///
/// # Errors
///
/// Returns structured I/O, syntax, alphabet, dimensions, annotation or limit
/// errors. Source locations use one-based rows/lines where available.
pub fn parse_alignment<A: Alphabet>(
    reader: impl BufRead,
    format: AlignmentFormat,
) -> Result<AlignmentDocument<A>, BioIoError> {
    let text = read_text(reader, format.name())?;
    let raw = match format {
        AlignmentFormat::Fasta => parse_fasta(&text)?,
        AlignmentFormat::Clustal => parse_clustal(&text)?,
        AlignmentFormat::Phylip { names, layout } => parse_phylip(&text, names, layout)?,
        AlignmentFormat::Stockholm => parse_stockholm(&text)?,
        AlignmentFormat::Nexus => parse_nexus::<A>(&text)?,
    };
    raw.validate::<A>(format)
}

/// Reads one alignment from a path through the shared plain/gzip input layer.
///
/// # Errors
///
/// Returns the same validation errors as [`parse_alignment`], or path/input errors.
pub fn read_alignment<A: Alphabet>(
    path: impl AsRef<Path>,
    format: AlignmentFormat,
) -> Result<AlignmentDocument<A>, BioIoError> {
    parse_alignment(open_input(path)?, format)
}

#[derive(Default)]
struct RawDocument {
    rows: Vec<AlignmentRow>,
    annotations: Vec<AlignmentAnnotation>,
    gaps: Vec<u8>,
    missing: Option<u8>,
}

impl RawDocument {
    fn new(gaps: &[u8], missing: Option<u8>) -> Self {
        Self {
            gaps: gaps.to_vec(),
            missing,
            ..Self::default()
        }
    }

    fn validate<A: Alphabet>(
        mut self,
        format: AlignmentFormat,
    ) -> Result<AlignmentDocument<A>, BioIoError> {
        let width = self
            .rows
            .first()
            .ok_or_else(|| {
                failure(
                    format,
                    IoErrorKind::InvalidStructure,
                    "alignment has no rows",
                )
            })?
            .aligned
            .len();
        let mut identifiers = HashSet::new();
        for (index, row) in self.rows.iter_mut().enumerate() {
            if row.identifier.trim().is_empty() || row.identifier.chars().any(char::is_control) {
                return Err(failure(
                    format,
                    IoErrorKind::InvalidHeader,
                    "row identifier is empty or contains a control character",
                )
                .at(index + 1, row.line));
            }
            if !identifiers.insert(row.identifier.clone()) {
                return Err(failure(
                    format,
                    IoErrorKind::InvalidStructure,
                    format!("duplicate row identifier {:?}", row.identifier),
                )
                .at(index + 1, row.line));
            }
            if row.aligned.len() != width {
                return Err(failure(
                    format,
                    IoErrorKind::LengthMismatch,
                    format!("row has {} columns; expected {width}", row.aligned.len()),
                )
                .at(index + 1, row.line));
            }
            let mut normalized = String::with_capacity(width);
            for (column, byte) in row.aligned.bytes().enumerate() {
                if self.gaps.contains(&byte) || Some(byte) == self.missing {
                    normalized.push(char::from(byte));
                } else if let Some(symbol) = Symbol::<A>::from_ascii(byte) {
                    normalized.push(char::from(symbol.as_ascii()));
                } else {
                    return Err(failure(
                        format,
                        IoErrorKind::InvalidSequence,
                        format!(
                            "invalid {} byte {byte} at zero-based alignment column {column}",
                            A::NAME
                        ),
                    )
                    .at(index + 1, row.line));
                }
            }
            row.aligned = normalized;
        }
        for annotation in &self.annotations {
            if let Some(name) = &annotation.sequence {
                if !identifiers.contains(name) {
                    return Err(failure(
                        format,
                        IoErrorKind::InvalidStructure,
                        format!("annotation refers to unknown sequence {name:?}"),
                    )
                    .at(1, annotation.line));
                }
            }
            if matches!(
                annotation.scope,
                AnnotationScope::Column | AnnotationScope::Residue
            ) && (!annotation.value.is_ascii() || annotation.value.len() != width)
            {
                return Err(failure(
                    format,
                    IoErrorKind::LengthMismatch,
                    format!(
                        "annotation {:?} has {} bytes; expected {width} ASCII column characters",
                        annotation.tag,
                        annotation.value.len()
                    ),
                )
                .at(1, annotation.line));
            }
        }
        Ok(AlignmentDocument {
            format,
            rows: self.rows,
            annotations: self.annotations,
            width,
            gaps: self.gaps,
            missing: self.missing,
            marker: PhantomData,
        })
    }
}

fn failure(format: AlignmentFormat, kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, format.name(), message)
}

fn row(identifier: impl Into<String>, aligned: String, line: usize) -> AlignmentRow {
    AlignmentRow {
        identifier: identifier.into(),
        aligned,
        line,
    }
}

fn compact(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_ascii_whitespace())
        .collect()
}

fn parse_fasta(text: &str) -> Result<RawDocument, BioIoError> {
    let format = AlignmentFormat::Fasta;
    let mut result = RawDocument::new(b"-", Some(b'?'));
    for (index, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            result
                .rows
                .push(row(header.trim(), String::new(), index + 1));
        } else if let Some(current) = result.rows.last_mut() {
            if line.starts_with(';') {
                return Err(failure(
                    format,
                    IoErrorKind::Unsupported,
                    "FASTA comment lines are not part of this profile",
                )
                .at(result.rows.len(), index + 1));
            }
            current.aligned.push_str(&compact(line));
        } else {
            return Err(failure(
                format,
                IoErrorKind::InvalidHeader,
                "expected a FASTA header before sequence data",
            )
            .at(1, index + 1));
        }
    }
    Ok(result)
}

struct ClustalFragment {
    identifier: String,
    sequence: String,
    count: Option<usize>,
    line: usize,
}

#[derive(Default)]
struct ClustalBlock {
    rows: Vec<ClustalFragment>,
    consensus: Option<(String, usize)>,
    sequence_column: usize,
}

fn parse_clustal(text: &str) -> Result<RawDocument, BioIoError> {
    let format = AlignmentFormat::Clustal;
    let mut lines = text.lines().enumerate().peekable();
    let (_, header) = lines
        .find(|(_, line)| !line.trim().is_empty())
        .ok_or_else(|| failure(format, IoErrorKind::InvalidHeader, "missing CLUSTAL header"))?;
    if header.split_whitespace().next() != Some("CLUSTAL") {
        return Err(failure(
            format,
            IoErrorKind::InvalidHeader,
            "expected a CLUSTAL header",
        ));
    }
    let mut result = RawDocument::new(b"-", None);
    let mut block = ClustalBlock::default();
    let mut consensus = String::new();
    let mut consensus_line = None;
    while let Some((index, line)) = lines.next() {
        if line.trim().is_empty() {
            flush_clustal(&mut result, &mut block, &mut consensus, &mut consensus_line)?;
            continue;
        }
        if line
            .bytes()
            .all(|byte| byte.is_ascii_whitespace() || b"*:.".contains(&byte))
        {
            clustal_consensus(&mut block, line, index + 1)?;
            continue;
        }
        if block.consensus.is_some() {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "sequence rows must precede consensus; separate blocks with a blank line",
            )
            .at(1, index + 1));
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if !(2..=3).contains(&fields.len()) {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "expected name, sequence and optional cumulative residue count",
            )
            .at(1, index + 1));
        }
        if block.rows.is_empty() {
            let prefix = line.find(fields[0]).unwrap_or(0) + fields[0].len();
            block.sequence_column =
                prefix + line[prefix..].len() - line[prefix..].trim_start().len();
        }
        let count = fields
            .get(2)
            .map(|value| parse_number(value, format, index + 1))
            .transpose()?;
        block.rows.push(ClustalFragment {
            identifier: fields[0].to_owned(),
            sequence: fields[1].to_owned(),
            count,
            line: index + 1,
        });
        if lines.peek().is_none() {
            flush_clustal(&mut result, &mut block, &mut consensus, &mut consensus_line)?;
        }
    }
    flush_clustal(&mut result, &mut block, &mut consensus, &mut consensus_line)?;
    if let Some(line) = consensus_line {
        result.annotations.push(AlignmentAnnotation {
            scope: AnnotationScope::Column,
            sequence: None,
            tag: "consensus".to_owned(),
            value: consensus,
            line,
        });
    }
    Ok(result)
}

fn clustal_consensus(
    block: &mut ClustalBlock,
    line: &str,
    number: usize,
) -> Result<(), BioIoError> {
    let format = AlignmentFormat::Clustal;
    let Some(first) = block.rows.first() else {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            "consensus appears before a sequence block",
        )
        .at(1, number));
    };
    if block.consensus.is_some() || !line.is_ascii() || line.len() < block.sequence_column {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            "duplicate or misaligned consensus line",
        )
        .at(1, number));
    }
    if !line[..block.sequence_column]
        .bytes()
        .all(|byte| byte == b' ' || byte == b'\t')
    {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            "consensus does not align with sequence columns",
        )
        .at(1, number));
    }
    let value = &line[block.sequence_column..];
    if value.len() > first.sequence.len() && !value[first.sequence.len()..].trim().is_empty() {
        return Err(failure(
            format,
            IoErrorKind::LengthMismatch,
            "consensus exceeds block width",
        )
        .at(1, number));
    }
    let mut value = value[..value.len().min(first.sequence.len())].to_owned();
    value.extend(std::iter::repeat_n(' ', first.sequence.len() - value.len()));
    block.consensus = Some((value, number));
    Ok(())
}

fn flush_clustal(
    result: &mut RawDocument,
    block: &mut ClustalBlock,
    consensus: &mut String,
    consensus_line: &mut Option<usize>,
) -> Result<(), BioIoError> {
    let Some(first) = block.rows.first() else {
        return Ok(());
    };
    let format = AlignmentFormat::Clustal;
    let width = first.sequence.len();
    let initial = result.rows.is_empty();
    if !initial && block.rows.len() != result.rows.len() {
        return Err(failure(
            format,
            IoErrorKind::LengthMismatch,
            "CLUSTAL block has a different row count",
        )
        .at(1, first.line));
    }
    for (index, fragment) in block.rows.iter().enumerate() {
        if fragment.sequence.len() != width {
            return Err(failure(
                format,
                IoErrorKind::LengthMismatch,
                "CLUSTAL block fragments must have equal widths",
            )
            .at(index + 1, fragment.line));
        }
        if initial {
            result
                .rows
                .push(row(&fragment.identifier, String::new(), fragment.line));
        }
        let target = &mut result.rows[index];
        if target.identifier != fragment.identifier {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "CLUSTAL blocks must retain sequence names and order",
            )
            .at(index + 1, fragment.line));
        }
        target.aligned.push_str(&fragment.sequence);
        if let Some(count) = fragment.count {
            if target.aligned.bytes().filter(|byte| *byte != b'-').count() != count {
                return Err(failure(
                    format,
                    IoErrorKind::LengthMismatch,
                    "incorrect cumulative CLUSTAL residue count",
                )
                .at(index + 1, fragment.line));
            }
        }
    }
    if let Some((value, line)) = block.consensus.take() {
        consensus_line.get_or_insert(line);
        consensus.push_str(&value);
    } else {
        consensus.extend(std::iter::repeat_n(' ', width));
    }
    block.rows.clear();
    Ok(())
}

fn parse_number(value: &str, format: AlignmentFormat, line: usize) -> Result<usize, BioIoError> {
    if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            format!("invalid unsigned dimension/count {value:?}"),
        )
        .at(1, line));
    }
    value.parse().map_err(|_| {
        failure(
            format,
            IoErrorKind::LimitExceeded,
            "dimension/count exceeds this platform's integer range",
        )
        .at(1, line)
    })
}

fn parse_phylip(
    text: &str,
    names: PhylipNameStyle,
    layout: PhylipLayout,
) -> Result<RawDocument, BioIoError> {
    let format = AlignmentFormat::Phylip { names, layout };
    let lines: Vec<_> = text
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .collect();
    let Some((header_line, header)) = lines.first() else {
        return Err(failure(
            format,
            IoErrorKind::InvalidHeader,
            "missing PHYLIP dimensions",
        ));
    };
    let fields: Vec<_> = header.split_whitespace().collect();
    if fields.len() != 2 {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "PHYLIP header must contain only taxon and column counts",
        )
        .at(1, header_line + 1));
    }
    let count = parse_number(fields[0], format, header_line + 1)?;
    let width = parse_number(fields[1], format, header_line + 1)?;
    if count == 0 {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            "PHYLIP taxon count must be positive",
        ));
    }
    let mut result = RawDocument::new(b"-", Some(b'?'));
    let mut cursor = 1;
    for index in 0..count {
        let (line_number, line) = lines
            .get(cursor)
            .ok_or_else(|| failure(format, IoErrorKind::Truncated, "missing PHYLIP taxon row"))?;
        cursor += 1;
        let (identifier, sequence) = phylip_first_line(line, names, format, line_number + 1)?;
        result.rows.push(row(identifier, sequence, line_number + 1));
        if layout == PhylipLayout::Sequential {
            while result.rows[index].aligned.len() < width {
                let (_, continuation) = lines.get(cursor).ok_or_else(|| {
                    failure(
                        format,
                        IoErrorKind::Truncated,
                        "truncated sequential PHYLIP row",
                    )
                    .at(index + 1, line_number + 1)
                })?;
                cursor += 1;
                result.rows[index].aligned.push_str(&compact(continuation));
            }
        }
    }
    if layout == PhylipLayout::Interleaved {
        equal_fragments(&result.rows, format)?;
        while result.rows[0].aligned.len() < width {
            let mut block_width = None;
            for (index, target) in result.rows.iter_mut().enumerate() {
                let (line_number, line) = lines.get(cursor).ok_or_else(|| {
                    failure(
                        format,
                        IoErrorKind::Truncated,
                        "truncated interleaved PHYLIP block",
                    )
                })?;
                cursor += 1;
                let fragment = compact(line);
                if block_width.is_some_and(|expected| expected != fragment.len()) {
                    return Err(failure(
                        format,
                        IoErrorKind::LengthMismatch,
                        "interleaved PHYLIP block fragments must have equal widths",
                    )
                    .at(index + 1, line_number + 1));
                }
                block_width = Some(fragment.len());
                target.aligned.push_str(&fragment);
            }
        }
    }
    if cursor != lines.len() {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "trailing PHYLIP rows, options, trees or multiple datasets are not accepted",
        )
        .at(1, lines[cursor].0 + 1));
    }
    validate_phylip_rows(&result.rows, format, width)?;
    Ok(result)
}

fn validate_phylip_rows(
    rows: &[AlignmentRow],
    format: AlignmentFormat,
    width: usize,
) -> Result<(), BioIoError> {
    for (index, target) in rows.iter().enumerate() {
        if target.aligned.len() != width {
            return Err(failure(
                format,
                IoErrorKind::LengthMismatch,
                format!(
                    "PHYLIP row has {} columns; header declares {width}",
                    target.aligned.len()
                ),
            )
            .at(index + 1, target.line));
        }
        if target
            .aligned
            .bytes()
            .any(|byte| byte.is_ascii_digit() || byte == b'.')
        {
            return Err(failure(format, IoErrorKind::Unsupported, "numbered sequence data and dot/match notation are not supported by this PHYLIP profile").at(index + 1, target.line));
        }
    }
    Ok(())
}

fn phylip_first_line(
    line: &str,
    names: PhylipNameStyle,
    format: AlignmentFormat,
    line_number: usize,
) -> Result<(String, String), BioIoError> {
    let (name, sequence) = match names {
        PhylipNameStyle::Strict => {
            if line.len() < 10 || !line.as_bytes()[..10].is_ascii() {
                return Err(failure(
                    format,
                    IoErrorKind::InvalidHeader,
                    "strict PHYLIP requires a ten-byte ASCII name field",
                )
                .at(1, line_number));
            }
            (&line[..10], &line[10..])
        }
        PhylipNameStyle::Relaxed => {
            let line = line.trim_start();
            let boundary = line.find(char::is_whitespace).unwrap_or(line.len());
            (&line[..boundary], &line[boundary..])
        }
    };
    Ok((name.trim().to_owned(), compact(sequence)))
}

fn equal_fragments(rows: &[AlignmentRow], format: AlignmentFormat) -> Result<(), BioIoError> {
    let expected = rows.first().map_or(0, |row| row.aligned.len());
    for (index, row) in rows.iter().enumerate() {
        if row.aligned.len() != expected {
            return Err(failure(
                format,
                IoErrorKind::LengthMismatch,
                "interleaved block fragments must have equal widths",
            )
            .at(index + 1, row.line));
        }
    }
    Ok(())
}

fn parse_stockholm(text: &str) -> Result<RawDocument, BioIoError> {
    let format = AlignmentFormat::Stockholm;
    let mut lines = text.lines().enumerate();
    if lines.next().map(|(_, line)| line.trim()) != Some("# STOCKHOLM 1.0") {
        return Err(failure(
            format,
            IoErrorKind::InvalidHeader,
            "expected # STOCKHOLM 1.0 on the first line",
        ));
    }
    let mut result = RawDocument::new(b"-.", Some(b'?'));
    let mut indices = HashMap::<String, usize>::new();
    let mut terminated = false;
    for (index, line) in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if terminated {
            return Err(failure(
                format,
                IoErrorKind::Unsupported,
                "only one Stockholm alignment is accepted; trailing data is not ignored",
            )
            .at(1, index + 1));
        }
        if line == "//" {
            terminated = true;
            continue;
        }
        if line.starts_with("#=") {
            stockholm_annotation(&mut result, line, index + 1)?;
        } else if let Some(comment) = line.strip_prefix('#') {
            result.annotations.push(AlignmentAnnotation {
                scope: AnnotationScope::Comment,
                sequence: None,
                tag: "comment".to_owned(),
                value: comment.to_owned(),
                line: index + 1,
            });
        } else {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() != 2 {
                return Err(failure(
                    format,
                    IoErrorKind::InvalidStructure,
                    "expected one sequence name and one aligned fragment",
                )
                .at(1, index + 1));
            }
            if let Some(existing) = indices.get(fields[0]) {
                result.rows[*existing].aligned.push_str(fields[1]);
            } else {
                indices.insert(fields[0].to_owned(), result.rows.len());
                result
                    .rows
                    .push(row(fields[0], fields[1].to_owned(), index + 1));
            }
        }
    }
    if !terminated {
        return Err(failure(
            format,
            IoErrorKind::Truncated,
            "missing Stockholm // terminator",
        ));
    }
    Ok(result)
}

fn take_field(text: &str) -> Option<(&str, &str)> {
    let text = text.trim_start();
    if text.is_empty() {
        return None;
    }
    let end = text.find(char::is_whitespace).unwrap_or(text.len());
    Some((&text[..end], text[end..].trim_start()))
}

fn stockholm_annotation(
    result: &mut RawDocument,
    line: &str,
    number: usize,
) -> Result<(), BioIoError> {
    let format = AlignmentFormat::Stockholm;
    let malformed = || {
        failure(
            format,
            IoErrorKind::InvalidStructure,
            "incomplete Stockholm annotation",
        )
        .at(1, number)
    };
    let (prefix, rest) = take_field(line).ok_or_else(malformed)?;
    let scope = match prefix {
        "#=GF" => AnnotationScope::File,
        "#=GS" => AnnotationScope::Sequence,
        "#=GC" => AnnotationScope::Column,
        "#=GR" => AnnotationScope::Residue,
        _ => {
            return Err(failure(
                format,
                IoErrorKind::Unsupported,
                format!("unsupported Stockholm markup {prefix:?}"),
            )
            .at(1, number));
        }
    };
    let (sequence, rest) = if matches!(scope, AnnotationScope::Sequence | AnnotationScope::Residue)
    {
        let (sequence, rest) = take_field(rest).ok_or_else(malformed)?;
        (Some(sequence.to_owned()), rest)
    } else {
        (None, rest)
    };
    let (tag, value) = take_field(rest).ok_or_else(malformed)?;
    if value.is_empty() {
        return Err(malformed());
    }
    if matches!(scope, AnnotationScope::Column | AnnotationScope::Residue) {
        if value.bytes().any(|byte| byte.is_ascii_whitespace()) {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "column annotations cannot contain whitespace",
            )
            .at(1, number));
        }
        if let Some(annotation) = result.annotations.iter_mut().find(|annotation| {
            annotation.scope == scope && annotation.sequence == sequence && annotation.tag == tag
        }) {
            annotation.value.push_str(value);
            return Ok(());
        }
    }
    result.annotations.push(AlignmentAnnotation {
        scope,
        sequence,
        tag: tag.to_owned(),
        value: value.to_owned(),
        line: number,
    });
    Ok(())
}

struct NexusCommand {
    text: String,
    line: usize,
}

/// Splits commands without treating quoted semicolons or comments as syntax.
fn nexus_commands(
    text: &str,
    first_line: usize,
) -> Result<(Vec<NexusCommand>, Vec<AlignmentAnnotation>), BioIoError> {
    let format = AlignmentFormat::Nexus;
    let mut characters = text.chars().peekable();
    let mut commands = Vec::new();
    let mut annotations = Vec::new();
    let mut command = String::new();
    let mut line = first_line;
    let mut start_line = first_line;
    let mut quoted = false;
    let mut has_token = false;
    while let Some(character) = characters.next() {
        if character == '\n' {
            line += 1;
        }
        if !has_token && !character.is_whitespace() && !matches!(character, '[' | ';') {
            start_line = line;
            has_token = true;
        }
        if character == '\'' {
            command.push(character);
            if quoted && characters.peek() == Some(&'\'') {
                command.push(characters.next().unwrap_or('\''));
            } else {
                quoted = !quoted;
            }
        } else if character == '[' && !quoted {
            let comment_line = line;
            let mut depth = 1usize;
            let mut value = String::new();
            for next in characters.by_ref() {
                if next == '\n' {
                    line += 1;
                    command.push('\n');
                }
                if next == '[' {
                    depth += 1;
                }
                if next == ']' {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                value.push(next);
            }
            if depth != 0 {
                return Err(
                    failure(format, IoErrorKind::Truncated, "unterminated NEXUS comment")
                        .at(1, comment_line),
                );
            }
            annotations.push(AlignmentAnnotation {
                scope: AnnotationScope::Comment,
                sequence: None,
                tag: "comment".to_owned(),
                value,
                line: comment_line,
            });
        } else if character == ']' && !quoted {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "unmatched NEXUS comment terminator",
            )
            .at(1, line));
        } else if character == ';' && !quoted {
            if has_token {
                commands.push(NexusCommand {
                    text: command.trim_start().to_owned(),
                    line: start_line,
                });
            }
            command.clear();
            has_token = false;
            start_line = line;
        } else {
            command.push(character);
        }
    }
    if quoted || !command.trim().is_empty() {
        return Err(failure(
            format,
            IoErrorKind::Truncated,
            "unterminated NEXUS quote or missing command semicolon",
        )
        .at(1, start_line));
    }
    Ok((commands, annotations))
}

fn nexus_assignments(text: &str, line: usize) -> Result<HashMap<String, String>, BioIoError> {
    let format = AlignmentFormat::Nexus;
    let spaced = text.replace('=', " = ");
    let fields: Vec<_> = spaced.split_whitespace().collect();
    if fields.len() % 3 != 0 {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "expected only key=value NEXUS options in this profile",
        )
        .at(1, line));
    }
    let mut result = HashMap::new();
    for triple in fields.chunks_exact(3) {
        if triple[1] != "="
            || result
                .insert(triple[0].to_ascii_uppercase(), triple[2].to_owned())
                .is_some()
        {
            return Err(failure(
                format,
                IoErrorKind::InvalidStructure,
                "invalid or repeated NEXUS option",
            )
            .at(1, line));
        }
    }
    Ok(result)
}

fn parse_nexus<A: Alphabet>(text: &str) -> Result<RawDocument, BioIoError> {
    let format = AlignmentFormat::Nexus;
    let prefix = text.len() - text.trim_start().len();
    let first_line = 1 + text[..prefix].bytes().filter(|byte| *byte == b'\n').count();
    let text = text.trim_start();
    if text
        .get(..6)
        .is_none_or(|header| !header.eq_ignore_ascii_case("#NEXUS"))
        || text
            .as_bytes()
            .get(6)
            .is_some_and(|byte| !byte.is_ascii_whitespace())
    {
        return Err(failure(
            format,
            IoErrorKind::InvalidHeader,
            "expected #NEXUS header",
        ));
    }
    let (commands, annotations) = nexus_commands(&text[6..], first_line)?;
    let begins_data = commands.first().is_some_and(|command| {
        command
            .text
            .split_whitespace()
            .map(str::to_ascii_uppercase)
            .eq(["BEGIN", "DATA"].map(str::to_owned))
    });
    if commands.len() != 5 || !begins_data {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "expected exactly BEGIN DATA, DIMENSIONS, FORMAT, MATRIX and END commands; extra blocks/commands are not ignored",
        ));
    }
    let dimension_text = command_body(&commands[1], "DIMENSIONS")?;
    let mut dimensions = nexus_assignments(dimension_text, commands[1].line)?;
    let count = dimensions
        .remove("NTAX")
        .ok_or_else(|| failure(format, IoErrorKind::InvalidStructure, "missing NTAX"))?;
    let width = dimensions
        .remove("NCHAR")
        .ok_or_else(|| failure(format, IoErrorKind::InvalidStructure, "missing NCHAR"))?;
    if !dimensions.is_empty() {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "unsupported DIMENSIONS option",
        ));
    }
    let count = parse_number(&count, format, commands[1].line)?;
    let width = parse_number(&width, format, commands[1].line)?;
    let (gap, missing) = nexus_format::<A>(&commands[2])?;
    let matrix = command_body(&commands[3], "MATRIX")?;
    let matrix_prefix = commands[3].text.len() - matrix.len();
    let matrix_line = commands[3].line
        + commands[3].text[..matrix_prefix]
            .bytes()
            .filter(|byte| *byte == b'\n')
            .count();
    let mut result = RawDocument::new(&[gap], Some(missing));
    result.annotations = annotations;
    for (index, line) in matrix.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let number = matrix_line + index;
        let (name, data) = nexus_row(line, number)?;
        let sequence = compact(data);
        if sequence.len() != width {
            return Err(failure(format, IoErrorKind::LengthMismatch, format!("NEXUS row has {} columns; expected {width}; wrapped/interleaved rows are unsupported", sequence.len())).at(result.rows.len() + 1, number));
        }
        result.rows.push(row(name, sequence, number));
    }
    if result.rows.len() != count {
        return Err(failure(
            format,
            IoErrorKind::LengthMismatch,
            format!(
                "NEXUS has {} rows; expected NTAX={count}",
                result.rows.len()
            ),
        ));
    }
    if !commands[4].text.trim().eq_ignore_ascii_case("END")
        && !commands[4].text.trim().eq_ignore_ascii_case("ENDBLOCK")
    {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "expected END after MATRIX",
        ));
    }
    Ok(result)
}

fn nexus_format<A: Alphabet>(command: &NexusCommand) -> Result<(u8, u8), BioIoError> {
    let format = AlignmentFormat::Nexus;
    let mut options = nexus_assignments(command_body(command, "FORMAT")?, command.line)?;
    let datatype = options.remove("DATATYPE").ok_or_else(|| {
        failure(
            format,
            IoErrorKind::Unsupported,
            "FORMAT DATATYPE must be explicit",
        )
    })?;
    if !datatype.eq_ignore_ascii_case(A::NAME) {
        return Err(failure(
            format,
            IoErrorKind::InvalidSequence,
            format!(
                "DATATYPE={datatype} does not match selected alphabet {}",
                A::NAME
            ),
        ));
    }
    if !["DNA", "RNA", "PROTEIN"]
        .iter()
        .any(|allowed| datatype.eq_ignore_ascii_case(allowed))
    {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "only DNA, RNA and PROTEIN datatypes are supported",
        ));
    }
    let gap = nexus_glyph::<A>(options.remove("GAP"), b'-', command.line)?;
    let missing = nexus_glyph::<A>(options.remove("MISSING"), b'?', command.line)?;
    if gap == missing {
        return Err(failure(
            format,
            IoErrorKind::InvalidStructure,
            "gap and missing-data symbols must differ",
        ));
    }
    if let Some(interleave) = options.remove("INTERLEAVE") {
        if !interleave.eq_ignore_ascii_case("NO") {
            return Err(failure(
                format,
                IoErrorKind::Unsupported,
                "interleaved NEXUS is not supported",
            ));
        }
    }
    if !options.is_empty() {
        return Err(failure(
            format,
            IoErrorKind::Unsupported,
            "unsupported FORMAT option (including MATCHCHAR, EQUATE, TOKENS or TRANSPOSE)",
        ));
    }
    Ok((gap, missing))
}

fn command_body<'a>(command: &'a NexusCommand, expected: &str) -> Result<&'a str, BioIoError> {
    let (keyword, rest) = take_field(&command.text).ok_or_else(|| {
        failure(
            AlignmentFormat::Nexus,
            IoErrorKind::InvalidStructure,
            "empty command",
        )
    })?;
    if !keyword.eq_ignore_ascii_case(expected) {
        return Err(failure(
            AlignmentFormat::Nexus,
            IoErrorKind::Unsupported,
            format!("expected {expected}, found {keyword}"),
        )
        .at(1, command.line));
    }
    Ok(rest)
}

fn nexus_glyph<A: Alphabet>(
    value: Option<String>,
    default: u8,
    line: usize,
) -> Result<u8, BioIoError> {
    let Some(value) = value else {
        return Ok(default);
    };
    if value.len() != 1 || !value.is_ascii() {
        return Err(failure(
            AlignmentFormat::Nexus,
            IoErrorKind::Unsupported,
            "GAP and MISSING must be single unquoted ASCII punctuation characters",
        )
        .at(1, line));
    }
    let byte = value.as_bytes()[0];
    if !byte.is_ascii_punctuation() || b"[];'\"=".contains(&byte) || A::index_of(byte).is_some() {
        return Err(failure(
            AlignmentFormat::Nexus,
            IoErrorKind::Unsupported,
            "unsupported GAP or MISSING glyph",
        )
        .at(1, line));
    }
    Ok(byte)
}

fn nexus_row(line: &str, number: usize) -> Result<(String, &str), BioIoError> {
    let format = AlignmentFormat::Nexus;
    if let Some(remainder) = line.strip_prefix('\'') {
        let mut characters = remainder.char_indices().peekable();
        let mut name = String::new();
        while let Some((index, character)) = characters.next() {
            if character == '\'' {
                if characters.peek().is_some_and(|(_, next)| *next == '\'') {
                    characters.next();
                    name.push('\'');
                } else {
                    let data = &remainder[index + 1..];
                    if !data.is_empty() && !data.starts_with(char::is_whitespace) {
                        return Err(failure(
                            format,
                            IoErrorKind::InvalidHeader,
                            "quoted NEXUS name must be separated from sequence data",
                        )
                        .at(1, number));
                    }
                    return Ok((name, data));
                }
            } else {
                name.push(character);
            }
        }
        Err(failure(
            format,
            IoErrorKind::Truncated,
            "unterminated quoted NEXUS name",
        )
        .at(1, number))
    } else {
        let (name, data) = take_field(line).ok_or_else(|| {
            failure(format, IoErrorKind::InvalidHeader, "missing NEXUS name").at(1, number)
        })?;
        if name
            .chars()
            .any(|character| "(){} /\\,:=\"'".contains(character))
        {
            return Err(failure(
                format,
                IoErrorKind::Unsupported,
                "punctuation in NEXUS taxon names must be single-quoted",
            )
            .at(1, number));
        }
        Ok((name.replace('_', " "), data))
    }
}
