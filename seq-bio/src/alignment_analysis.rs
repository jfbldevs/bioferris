//! Immutable alignment analysis without weakening computational alignment rules.
//!
//! The analysis rectangle can retain imported missing cells, gap glyphs, opaque
//! annotations and all-gap columns. Coordinates are zero-based; slices are strict
//! half-open ranges and retain the original source positions.

use std::fmt;
use std::ops::Range;

use align::prelude::{MultipleAlignment, MultipleAlignmentView, PairwiseAlignmentView};
use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};

use crate::Strand;
use crate::io::alignment::{AlignmentAnnotation, AlignmentDocument};

mod sealed {
    pub trait Sealed {
        const UNKNOWN: u8;
        fn ambiguous(byte: u8) -> bool;
    }
}

/// A supported biological alphabet with explicitly defined ambiguity semantics.
///
/// This sealed trait supports DNA, RNA and extended stop-free protein only.
/// Protein residues outside the standard 20, including U and O, are ambiguous
/// for these analysis policies; the source cells themselves are preserved.
pub trait AlignmentAnalysisAlphabet: Alphabet + sealed::Sealed {}

impl AlignmentAnalysisAlphabet for Dna {}
impl AlignmentAnalysisAlphabet for Rna {}
impl AlignmentAnalysisAlphabet for Protein {}

impl sealed::Sealed for Dna {
    const UNKNOWN: u8 = b'N';
    fn ambiguous(byte: u8) -> bool {
        !b"ACGT".contains(&byte)
    }
}

impl sealed::Sealed for Rna {
    const UNKNOWN: u8 = b'N';
    fn ambiguous(byte: u8) -> bool {
        !b"ACGU".contains(&byte)
    }
}

impl sealed::Sealed for Protein {
    const UNKNOWN: u8 = b'X';
    fn ambiguous(byte: u8) -> bool {
        !b"ACDEFGHIKLMNPQRSTVWY".contains(&byte)
    }
}

/// Limits checked before allocating a new alignment rectangle.
///
/// They bound dimensions and cell count, not exact process memory or annotation
/// bytes already retained by an imported document. Zero limits are valid.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentAnalysisLimits {
    /// Maximum row count, including rows with no residues.
    pub max_rows: usize,
    /// Maximum column count, including when there are no rows.
    pub max_columns: usize,
    /// Maximum product of row count and column count.
    pub max_cells: usize,
}

impl Default for AlignmentAnalysisLimits {
    fn default() -> Self {
        Self {
            max_rows: 100_000,
            max_columns: 1_000_000,
            max_cells: 10_000_000,
        }
    }
}

impl AlignmentAnalysisLimits {
    fn check(self, rows: usize, columns: usize) -> Result<(), AlignmentAnalysisError> {
        for (resource, observed, limit) in [
            ("rows", rows, self.max_rows),
            ("columns", columns, self.max_columns),
        ] {
            if observed > limit {
                return Err(AlignmentAnalysisError::LimitExceeded {
                    resource,
                    limit,
                    observed,
                    rows,
                    columns,
                });
            }
        }
        let cells = rows
            .checked_mul(columns)
            .ok_or(AlignmentAnalysisError::CountOverflow)?;
        if cells > self.max_cells {
            return Err(AlignmentAnalysisError::LimitExceeded {
                resource: "cells",
                limit: self.max_cells,
                observed: cells,
                rows,
                columns,
            });
        }
        Ok(())
    }
}

/// One preserved biological, gap or missing-data cell.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentCell<A: AlignmentAnalysisAlphabet> {
    /// A validated uppercase symbol, including literal alphabet ambiguity.
    Residue(Symbol<A>),
    /// An explicit gap, retaining its source ASCII glyph.
    Gap(u8),
    /// Missing data, retaining its source ASCII glyph without guessing a residue.
    Missing(u8),
}

impl<A: AlignmentAnalysisAlphabet> AlignmentCell<A> {
    /// Returns the original normalized residue or source gap/missing glyph.
    #[must_use]
    pub fn as_ascii(self) -> u8 {
        match self {
            Self::Residue(symbol) => symbol.as_ascii(),
            Self::Gap(byte) | Self::Missing(byte) => byte,
        }
    }
}

/// A checked source interval and its orientation in the represented row.
///
/// The source length is the full original sequence length, not the aligned
/// interval length. Declaring reverse orientation does not complement symbols.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RowCoordinates {
    source_length: usize,
    range: Range<usize>,
    strand: Strand,
}

impl RowCoordinates {
    /// Validates an in-bounds source interval; attaching it additionally checks
    /// that its length equals the row's represented residue count.
    ///
    /// # Errors
    /// Returns `InvalidCoordinates` for reversed or out-of-source ranges.
    pub fn new(
        source_length: usize,
        range: Range<usize>,
        strand: Strand,
    ) -> Result<Self, AlignmentAnalysisError> {
        if range.start > range.end || range.end > source_length {
            return Err(AlignmentAnalysisError::InvalidCoordinates {
                row: None,
                source_length,
                range,
                represented_residues: None,
            });
        }
        Ok(Self {
            source_length,
            range,
            strand,
        })
    }

    /// Returns the complete original source length.
    #[must_use]
    pub const fn source_length(&self) -> usize {
        self.source_length
    }

    /// Returns the ascending half-open represented source interval.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// Returns the orientation of the represented residues.
    #[must_use]
    pub const fn strand(&self) -> Strand {
        self.strand
    }

    fn position(&self, ordinal: usize) -> usize {
        match self.strand {
            Strand::Forward => self.range.start + ordinal,
            Strand::Reverse => self.range.end - 1 - ordinal,
        }
    }

    fn slice(&self, start: usize, end: usize) -> Self {
        let range = match self.strand {
            Strand::Forward => self.range.start + start..self.range.start + end,
            Strand::Reverse => self.range.end - end..self.range.end - start,
        };
        Self {
            source_length: self.source_length,
            range,
            strand: self.strand,
        }
    }
}

/// One named immutable row and its optional known source coordinates.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentAnalysisRow<A: AlignmentAnalysisAlphabet> {
    identifier: String,
    cells: Vec<AlignmentCell<A>>,
    coordinates: Option<RowCoordinates>,
}

impl<A: AlignmentAnalysisAlphabet> AlignmentAnalysisRow<A> {
    /// Returns the original row identifier.
    #[must_use]
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    /// Returns the immutable cells in alignment-column order.
    #[must_use]
    pub fn cells(&self) -> &[AlignmentCell<A>] {
        &self.cells
    }

    /// Returns known coordinates; a missing-data row has no inferred origin.
    #[must_use]
    pub const fn coordinates(&self) -> Option<&RowCoordinates> {
        self.coordinates.as_ref()
    }

    /// Renders the row while retaining source gap and missing-data glyphs.
    #[must_use]
    pub fn aligned(&self) -> String {
        self.cells
            .iter()
            .map(|cell| char::from(cell.as_ascii()))
            .collect()
    }
}

/// Explicit handling of opaque annotations on a nontrivial derivative.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum AlignmentAnnotationPolicy {
    /// Refuse the transform if it would retain potentially stale metadata.
    #[default]
    Reject,
    /// Remove all annotations and return their ordered original metadata.
    Discard,
}

/// A successful slice plus every explicitly removed original annotation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentSlice<A: AlignmentAnalysisAlphabet> {
    alignment: AlignmentAnalysis<A>,
    discarded_annotations: Vec<AlignmentAnnotation>,
}

impl<A: AlignmentAnalysisAlphabet> AlignmentSlice<A> {
    /// Borrows the derived rectangle.
    #[must_use]
    pub const fn alignment(&self) -> &AlignmentAnalysis<A> {
        &self.alignment
    }

    /// Returns removed annotations in source order, including scopes, targets,
    /// tags, source lines and complete values.
    #[must_use]
    pub fn discarded_annotations(&self) -> &[AlignmentAnnotation] {
        &self.discarded_annotations
    }

    /// Separates the derivative from its ordered metadata-loss report.
    #[must_use]
    pub fn into_parts(self) -> (AlignmentAnalysis<A>, Vec<AlignmentAnnotation>) {
        (self.alignment, self.discarded_annotations)
    }
}

/// A checked analysis rectangle independent of computational MSA restrictions.
///
/// Empty dimensions, all-gap columns and missing cells are retained. Constructors
/// preserve identifiers and annotations; they never rerun an alignment kernel.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentAnalysis<A: AlignmentAnalysisAlphabet> {
    rows: Vec<AlignmentAnalysisRow<A>>,
    width: usize,
    annotations: Vec<AlignmentAnnotation>,
    limits: AlignmentAnalysisLimits,
}

impl<A: AlignmentAnalysisAlphabet> AlignmentAnalysis<A> {
    /// Copies a canonical MSA under default limits.
    ///
    /// # Errors
    /// Returns `LimitExceeded` or `CountOverflow` before model allocation.
    pub fn from_multiple(alignment: &MultipleAlignment<A>) -> Result<Self, AlignmentAnalysisError> {
        Self::from_multiple_with_limits(alignment, AlignmentAnalysisLimits::default())
    }

    /// Copies a canonical MSA with explicit dimension budgets.
    ///
    /// # Errors
    /// Returns `LimitExceeded` or `CountOverflow` before model allocation.
    pub fn from_multiple_with_limits(
        alignment: &MultipleAlignment<A>,
        limits: AlignmentAnalysisLimits,
    ) -> Result<Self, AlignmentAnalysisError> {
        limits.check(alignment.sequence_count(), alignment.width())?;
        let rows = alignment
            .rows()
            .iter()
            .map(|row| {
                let cells: Vec<_> = row
                    .symbols()
                    .iter()
                    .map(|cell| cell.map_or(AlignmentCell::Gap(b'-'), AlignmentCell::Residue))
                    .collect();
                analysis_row(row.identifier().to_owned(), cells)
            })
            .collect();
        Ok(Self {
            rows,
            width: alignment.width(),
            annotations: Vec::new(),
            limits,
        })
    }

    /// Copies the canonical MSA borrowed from an accepted algorithm result.
    ///
    /// # Errors
    /// Returns the same errors as [`Self::from_multiple`].
    pub fn from_multiple_result<R>(result: &R) -> Result<Self, AlignmentAnalysisError>
    where
        R: MultipleAlignmentView<Alphabet = A> + ?Sized,
    {
        Self::from_multiple(result.alignment())
    }

    /// Copies a result's canonical MSA under explicit dimension budgets.
    ///
    /// # Errors
    /// Returns the same errors as [`Self::from_multiple_with_limits`].
    pub fn from_multiple_result_with_limits<R>(
        result: &R,
        limits: AlignmentAnalysisLimits,
    ) -> Result<Self, AlignmentAnalysisError>
    where
        R: MultipleAlignmentView<Alphabet = A> + ?Sized,
    {
        Self::from_multiple_with_limits(result.alignment(), limits)
    }

    /// Renders an accepted pairwise traceback against both complete originals.
    ///
    /// Counts later use actual symbols, not the traceback's equality labels.
    /// Row identifiers are `first` and `second`; source ranges and lengths are
    /// retained, including for local alignment.
    ///
    /// # Errors
    /// Returns dimension-limit errors or `InvalidPairwise` for a traceback that
    /// cannot render against the supplied original sequences.
    pub fn from_pairwise<R: PairwiseAlignmentView + ?Sized>(
        first: &Sequence<A>,
        second: &Sequence<A>,
        result: &R,
    ) -> Result<Self, AlignmentAnalysisError> {
        Self::from_pairwise_with_limits(first, second, result, AlignmentAnalysisLimits::default())
    }

    /// Renders a pairwise traceback with explicit dimension budgets.
    ///
    /// # Errors
    /// Returns the same errors as [`Self::from_pairwise`].
    pub fn from_pairwise_with_limits<R: PairwiseAlignmentView + ?Sized>(
        first: &Sequence<A>,
        second: &Sequence<A>,
        result: &R,
        limits: AlignmentAnalysisLimits,
    ) -> Result<Self, AlignmentAnalysisError> {
        let alignment = result.alignment();
        limits.check(2, alignment.columns())?;
        let rendered = alignment.render(first, second).map_err(|error| {
            AlignmentAnalysisError::InvalidPairwise {
                message: error.to_string(),
            }
        })?;
        let mut rows = Vec::with_capacity(2);
        for (identifier, text, length, range) in [
            (
                "first",
                rendered.first(),
                first.len(),
                alignment.first_range(),
            ),
            (
                "second",
                rendered.second(),
                second.len(),
                alignment.second_range(),
            ),
        ] {
            let mut row = analysis_row(identifier.to_owned(), parse_cells(text, b"-", None)?);
            row.coordinates = Some(RowCoordinates::new(length, range, Strand::Forward)?);
            rows.push(row);
        }
        Ok(Self {
            rows,
            width: alignment.columns(),
            annotations: Vec::new(),
            limits,
        })
    }

    /// Copies an imported document without discarding its cells or annotations.
    ///
    /// A row containing missing data retains no guessed source coordinates.
    ///
    /// # Errors
    /// Returns dimension-limit errors or `InvalidCell` for inconsistent stored
    /// data; parser validation normally makes `InvalidCell` unreachable.
    pub fn from_document(document: &AlignmentDocument<A>) -> Result<Self, AlignmentAnalysisError> {
        Self::from_document_with_limits(document, AlignmentAnalysisLimits::default())
    }

    /// Copies an imported document with explicit dimension budgets.
    ///
    /// # Errors
    /// Returns the same errors as [`Self::from_document`].
    pub fn from_document_with_limits(
        document: &AlignmentDocument<A>,
        limits: AlignmentAnalysisLimits,
    ) -> Result<Self, AlignmentAnalysisError> {
        limits.check(document.rows().len(), document.width())?;
        let rows = document
            .rows()
            .iter()
            .map(|row| {
                Ok(analysis_row(
                    row.identifier().to_owned(),
                    parse_cells(
                        row.aligned(),
                        document.gap_symbols(),
                        document.missing_symbol(),
                    )?,
                ))
            })
            .collect::<Result<_, AlignmentAnalysisError>>()?;
        Ok(Self {
            rows,
            width: document.width(),
            annotations: document.annotations().to_vec(),
            limits,
        })
    }

    /// Returns rows in their stable alignment order.
    #[must_use]
    pub fn rows(&self) -> &[AlignmentAnalysisRow<A>] {
        &self.rows
    }

    /// Returns the explicit number of rows, including empty rows.
    #[must_use]
    pub fn sequence_count(&self) -> usize {
        self.rows.len()
    }

    /// Returns the explicit number of columns, including when there are no rows.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Returns original opaque annotations in source order.
    #[must_use]
    pub fn annotations(&self) -> &[AlignmentAnnotation] {
        &self.annotations
    }

    /// Returns dimension limits inherited by derivatives.
    #[must_use]
    pub const fn limits(&self) -> AlignmentAnalysisLimits {
        self.limits
    }

    /// Attaches an explicit source origin to one row, returning a new rectangle.
    ///
    /// # Errors
    /// Returns `InvalidRow`, `MissingCoordinates`, or `InvalidCoordinates` if the
    /// row contains missing data or interval length differs from residue count.
    pub fn with_row_coordinates(
        &self,
        row: usize,
        coordinates: RowCoordinates,
    ) -> Result<Self, AlignmentAnalysisError> {
        let selected = self.checked_row(row)?;
        if selected
            .cells
            .iter()
            .any(|cell| matches!(cell, AlignmentCell::Missing(_)))
        {
            return Err(AlignmentAnalysisError::MissingCoordinates { row });
        }
        let residues = residue_count(&selected.cells);
        if coordinates.range.len() != residues {
            return Err(AlignmentAnalysisError::InvalidCoordinates {
                row: Some(row),
                source_length: coordinates.source_length,
                range: coordinates.range,
                represented_residues: Some(residues),
            });
        }
        let mut result = self.clone();
        result.rows[row].coordinates = Some(coordinates);
        Ok(result)
    }

    /// Copies a strict row range while retaining original columns and origins.
    ///
    /// # Errors
    /// Returns `InvalidRange` or `AnnotationLoss` for a nontrivial annotated slice.
    pub fn slice_rows(&self, rows: Range<usize>) -> Result<Self, AlignmentAnalysisError> {
        Ok(self
            .slice_with(rows, 0..self.width, AlignmentAnnotationPolicy::Reject)?
            .alignment)
    }

    /// Copies a strict column range while retaining original source positions.
    ///
    /// # Errors
    /// Returns `InvalidRange` or `AnnotationLoss` for a nontrivial annotated slice.
    pub fn slice_columns(&self, columns: Range<usize>) -> Result<Self, AlignmentAnalysisError> {
        Ok(self
            .slice_with(
                0..self.rows.len(),
                columns,
                AlignmentAnnotationPolicy::Reject,
            )?
            .alignment)
    }

    /// Copies strict row/column ranges with explicit annotation handling.
    ///
    /// Identity slices retain annotations even with discard permission. Every
    /// other annotated transform either refuses or reports all removed metadata.
    /// No generic annotation substring projection is attempted.
    ///
    /// # Errors
    /// Returns `InvalidRange`, `AnnotationLoss`, or dimension arithmetic errors.
    pub fn slice_with(
        &self,
        rows: Range<usize>,
        columns: Range<usize>,
        annotations: AlignmentAnnotationPolicy,
    ) -> Result<AlignmentSlice<A>, AlignmentAnalysisError> {
        checked_range("rows", &rows, self.rows.len())?;
        checked_range("columns", &columns, self.width)?;
        let identity = rows == (0..self.rows.len()) && columns == (0..self.width);
        if !identity
            && !self.annotations.is_empty()
            && annotations == AlignmentAnnotationPolicy::Reject
        {
            return Err(AlignmentAnalysisError::AnnotationLoss {
                annotations: self.annotations.clone(),
            });
        }
        self.limits.check(rows.len(), columns.len())?;
        let selected = self.rows[rows]
            .iter()
            .map(|row| {
                let before = residue_count(&row.cells[..columns.start]);
                let retained = residue_count(&row.cells[columns.clone()]);
                AlignmentAnalysisRow {
                    identifier: row.identifier.clone(),
                    cells: row.cells[columns.clone()].to_vec(),
                    coordinates: row
                        .coordinates
                        .as_ref()
                        .map(|origin| origin.slice(before, before + retained)),
                }
            })
            .collect();
        Ok(AlignmentSlice {
            alignment: Self {
                rows: selected,
                width: columns.len(),
                annotations: if identity {
                    self.annotations.clone()
                } else {
                    Vec::new()
                },
                limits: self.limits,
            },
            discarded_annotations: if identity {
                Vec::new()
            } else {
                self.annotations.clone()
            },
        })
    }

    /// Maps a column to its original source position; gaps return `None`.
    ///
    /// # Errors
    /// Returns invalid row/column errors or `MissingCoordinates` when a source
    /// mapping is unknown. Missing data elsewhere in a row also prevents guessing.
    pub fn column_to_sequence(
        &self,
        row: usize,
        column: usize,
    ) -> Result<Option<usize>, AlignmentAnalysisError> {
        let selected = self.checked_row(row)?;
        let cell =
            selected
                .cells
                .get(column)
                .copied()
                .ok_or(AlignmentAnalysisError::InvalidColumn {
                    column,
                    columns: self.width,
                })?;
        let coordinates = ensure_known_row(selected, row)?;
        match cell {
            AlignmentCell::Residue(_) => Ok(Some(
                coordinates.position(residue_count(&selected.cells[..column])),
            )),
            AlignmentCell::Gap(_) => Ok(None),
            AlignmentCell::Missing(_) => Err(AlignmentAnalysisError::MissingCoordinates { row }),
        }
    }

    /// Maps a valid original source position to a column, or `None` if unaligned.
    ///
    /// # Errors
    /// Returns `InvalidRow`, `MissingCoordinates`, or `InvalidSourcePosition`.
    pub fn sequence_to_column(
        &self,
        row: usize,
        position: usize,
    ) -> Result<Option<usize>, AlignmentAnalysisError> {
        let selected = self.checked_row(row)?;
        let coordinates = ensure_known_row(selected, row)?;
        if position >= coordinates.source_length {
            return Err(AlignmentAnalysisError::InvalidSourcePosition {
                row,
                position,
                source_length: coordinates.source_length,
            });
        }
        if !coordinates.range.contains(&position) {
            return Ok(None);
        }
        let ordinal = match coordinates.strand {
            Strand::Forward => position - coordinates.range.start,
            Strand::Reverse => coordinates.range.end - 1 - position,
        };
        Ok(selected
            .cells
            .iter()
            .enumerate()
            .filter(|(_, cell)| matches!(cell, AlignmentCell::Residue(_)))
            .nth(ordinal)
            .map(|(column, _)| column))
    }

    fn checked_row(&self, row: usize) -> Result<&AlignmentAnalysisRow<A>, AlignmentAnalysisError> {
        self.rows
            .get(row)
            .ok_or(AlignmentAnalysisError::InvalidRow {
                row,
                rows: self.rows.len(),
            })
    }
}

fn analysis_row<A: AlignmentAnalysisAlphabet>(
    identifier: String,
    cells: Vec<AlignmentCell<A>>,
) -> AlignmentAnalysisRow<A> {
    let residues = residue_count(&cells);
    let coordinates = (!cells
        .iter()
        .any(|cell| matches!(cell, AlignmentCell::Missing(_))))
    .then_some(RowCoordinates {
        source_length: residues,
        range: 0..residues,
        strand: Strand::Forward,
    });
    AlignmentAnalysisRow {
        identifier,
        cells,
        coordinates,
    }
}

fn parse_cells<A: AlignmentAnalysisAlphabet>(
    text: &str,
    gaps: &[u8],
    missing: Option<u8>,
) -> Result<Vec<AlignmentCell<A>>, AlignmentAnalysisError> {
    text.bytes()
        .enumerate()
        .map(|(column, byte)| {
            if Some(byte) == missing {
                Ok(AlignmentCell::Missing(byte))
            } else if gaps.contains(&byte) {
                Ok(AlignmentCell::Gap(byte))
            } else {
                Symbol::from_ascii(byte)
                    .map(AlignmentCell::Residue)
                    .ok_or(AlignmentAnalysisError::InvalidCell { column, byte })
            }
        })
        .collect()
}

fn residue_count<A: AlignmentAnalysisAlphabet>(cells: &[AlignmentCell<A>]) -> usize {
    cells
        .iter()
        .filter(|cell| matches!(cell, AlignmentCell::Residue(_)))
        .count()
}

fn ensure_known_row<A: AlignmentAnalysisAlphabet>(
    selected: &AlignmentAnalysisRow<A>,
    row: usize,
) -> Result<&RowCoordinates, AlignmentAnalysisError> {
    selected
        .coordinates
        .as_ref()
        .ok_or(AlignmentAnalysisError::MissingCoordinates { row })
}

fn checked_range(
    axis: &'static str,
    range: &Range<usize>,
    length: usize,
) -> Result<(), AlignmentAnalysisError> {
    if range.start > range.end || range.end > length {
        return Err(AlignmentAnalysisError::InvalidRange {
            axis,
            start: range.start,
            end: range.end,
            length,
        });
    }
    Ok(())
}

/// How ambiguous residues affect statistics or consensus evidence.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum AlignmentAmbiguityPolicy {
    /// Treat each accepted alphabet symbol as a literal distinct residue.
    #[default]
    Literal,
    /// Exclude ambiguity-containing evidence rather than expanding it.
    Ignore,
    /// Refuse the first ambiguous residue in the documented traversal order.
    Reject,
}

/// The exact pair count used as the identity denominator.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum IdentityDenominator {
    /// Compare retained residue/residue pairs only.
    #[default]
    ResiduePairs,
    /// Include single-gap pairs but not double gaps or excluded evidence.
    ResidueAndSingleGap,
    /// Include every unordered row pair in every column, including missing data.
    AllPairs,
}

/// Explicit analysis policies; defaults use literal residue-pair identity.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct AlignmentStatisticsOptions {
    /// Handling of residues outside the canonical alphabet subset.
    pub ambiguity: AlignmentAmbiguityPolicy,
    /// Exact denominator used for the reported identity fraction.
    pub denominator: IdentityDenominator,
}

/// Exact cell counts; `ambiguous` is a subset of `residues`, not another state.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct AlignmentCellCounts {
    /// All represented residue cells, including ambiguous ones.
    pub residues: u64,
    /// Explicit gap cells.
    pub gaps: u64,
    /// Explicit missing-data cells.
    pub missing: u64,
    /// Residues outside the standard concrete nucleotide/20-amino-acid subset.
    pub ambiguous: u64,
}

impl AlignmentCellCounts {
    fn add<A: AlignmentAnalysisAlphabet>(
        &mut self,
        cell: AlignmentCell<A>,
    ) -> Result<(), AlignmentAnalysisError> {
        match cell {
            AlignmentCell::Residue(symbol) => {
                self.residues = add(self.residues, 1)?;
                if A::ambiguous(symbol.as_ascii()) {
                    self.ambiguous = add(self.ambiguous, 1)?;
                }
            }
            AlignmentCell::Gap(_) => self.gaps = add(self.gaps, 1)?,
            AlignmentCell::Missing(_) => self.missing = add(self.missing, 1)?,
        }
        Ok(())
    }
}

/// Disjoint run counts; a wholly gapped nonempty row is leading only.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct AlignmentGapRuns {
    /// Gap runs that start at column zero.
    pub leading: u64,
    /// Runs touching neither alignment boundary; missing cells separate runs.
    pub internal: u64,
    /// Runs ending at the final column, excluding wholly gapped rows.
    pub trailing: u64,
}

/// Exact disjoint unordered-pair categories accumulated across columns.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct AlignmentPairCounts {
    /// All unordered row pairs across all represented columns.
    pub total: u64,
    /// Pairs containing at least one missing cell, classified first.
    pub missing: u64,
    /// Nonmissing pairs containing ambiguity under the ignore policy.
    pub excluded_ambiguity: u64,
    /// Retained pairs with two gaps; never residue identities.
    pub double_gap: u64,
    /// Retained pairs with one residue and one gap.
    pub single_gap: u64,
    /// Retained pairs with two residues, including identities and mismatches.
    pub residue_pairs: u64,
    /// Literal equal-symbol residue pairs.
    pub identities: u64,
    /// Literal unequal-symbol residue pairs.
    pub mismatches: u64,
}

impl AlignmentPairCounts {
    fn accumulate(&mut self, other: Self) -> Result<(), AlignmentAnalysisError> {
        self.total = add(self.total, other.total)?;
        self.missing = add(self.missing, other.missing)?;
        self.excluded_ambiguity = add(self.excluded_ambiguity, other.excluded_ambiguity)?;
        self.double_gap = add(self.double_gap, other.double_gap)?;
        self.single_gap = add(self.single_gap, other.single_gap)?;
        self.residue_pairs = add(self.residue_pairs, other.residue_pairs)?;
        self.identities = add(self.identities, other.identities)?;
        self.mismatches = add(self.mismatches, other.mismatches)?;
        Ok(())
    }
}

/// Exact row metrics, retaining the original row identifier.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentRowStatistics {
    /// Original row identifier.
    pub identifier: String,
    /// Exact residue/gap/missing counts independent of comparison policy.
    pub cells: AlignmentCellCounts,
    /// Exact leading/internal/trailing gap-run partition.
    pub gap_runs: AlignmentGapRuns,
}

/// Exact cell and policy-qualified comparison counts for one column.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentColumnStatistics {
    /// Zero-based column index in the analyzed rectangle.
    pub column: usize,
    /// Exact residue/gap/missing cell counts.
    pub cells: AlignmentCellCounts,
    /// Unordered-pair categories under the selected ambiguity policy.
    pub pairs: AlignmentPairCounts,
}

/// Aggregate analysis plus exact row/column evidence for its identity fraction.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentStatistics {
    /// Explicit row count.
    pub rows: usize,
    /// Explicit column count.
    pub columns: usize,
    /// Exact counts across the complete rectangle.
    pub cells: AlignmentCellCounts,
    /// Sum of the disjoint row gap-run categories.
    pub gap_runs: AlignmentGapRuns,
    /// Row metrics in stable source order.
    pub row_statistics: Vec<AlignmentRowStatistics>,
    /// Column metrics in increasing column order.
    pub column_statistics: Vec<AlignmentColumnStatistics>,
    /// Exact accumulated unordered-pair comparison counts.
    pub pairs: AlignmentPairCounts,
    /// Exact denominator selected by the requested policy.
    pub denominator: u64,
    /// Literal identities divided by denominator, absent when denominator is zero.
    pub identity: Option<f64>,
    /// Policies used for this readout.
    pub options: AlignmentStatisticsOptions,
}

impl<A: AlignmentAnalysisAlphabet> AlignmentAnalysis<A> {
    /// Computes exact cells, gap runs and unordered-pair comparisons.
    ///
    /// Per-column symbol frequencies avoid a quadratic pair matrix. Time is
    /// O(rows * columns + alphabet size * columns), with row/column-sized output.
    /// Ambiguity rejection visits rows then columns, before publishing a result.
    /// Counts converted into the final `f64` fraction may lose integer precision
    /// above 2^53; the exact numerator and denominator remain available.
    ///
    /// # Errors
    /// Returns `RejectedCell` for reject-policy ambiguity or `CountOverflow`.
    pub fn statistics(
        &self,
        options: AlignmentStatisticsOptions,
    ) -> Result<AlignmentStatistics, AlignmentAnalysisError> {
        let mut result = AlignmentStatistics {
            rows: self.rows.len(),
            columns: self.width,
            cells: AlignmentCellCounts::default(),
            gap_runs: AlignmentGapRuns::default(),
            row_statistics: Vec::with_capacity(self.rows.len()),
            column_statistics: Vec::with_capacity(self.width),
            pairs: AlignmentPairCounts::default(),
            denominator: 0,
            identity: None,
            options,
        };
        for (row_index, row) in self.rows.iter().enumerate() {
            let mut counts = AlignmentCellCounts::default();
            for (column, &cell) in row.cells.iter().enumerate() {
                if let AlignmentCell::Residue(symbol) = cell {
                    if options.ambiguity == AlignmentAmbiguityPolicy::Reject
                        && A::ambiguous(symbol.as_ascii())
                    {
                        return Err(AlignmentAnalysisError::RejectedCell {
                            row: row_index,
                            column,
                            byte: symbol.as_ascii(),
                            kind: "ambiguity",
                        });
                    }
                }
                counts.add(cell)?;
                result.cells.add(cell)?;
            }
            let runs = gap_runs(&row.cells)?;
            result.gap_runs.leading = add(result.gap_runs.leading, runs.leading)?;
            result.gap_runs.internal = add(result.gap_runs.internal, runs.internal)?;
            result.gap_runs.trailing = add(result.gap_runs.trailing, runs.trailing)?;
            result.row_statistics.push(AlignmentRowStatistics {
                identifier: row.identifier.clone(),
                cells: counts,
                gap_runs: runs,
            });
        }
        for column in 0..self.width {
            let mut counts = AlignmentCellCounts::default();
            let mut frequencies = vec![0_u64; A::SYMBOLS.len()];
            for row in &self.rows {
                let cell = row.cells[column];
                counts.add(cell)?;
                if let AlignmentCell::Residue(symbol) = cell {
                    if options.ambiguity != AlignmentAmbiguityPolicy::Ignore
                        || !A::ambiguous(symbol.as_ascii())
                    {
                        let frequency = &mut frequencies[symbol.alphabet_index()];
                        *frequency = add(*frequency, 1)?;
                    }
                }
            }
            let total_rows = as_count(self.rows.len())?;
            let nonmissing = total_rows - counts.missing;
            let ignored = if options.ambiguity == AlignmentAmbiguityPolicy::Ignore {
                counts.ambiguous
            } else {
                0
            };
            let retained = nonmissing - ignored;
            let residues = counts.residues - ignored;
            let identities = frequencies
                .into_iter()
                .try_fold(0, |sum, count| add(sum, choose_two(count)?))?;
            let residue_pairs = choose_two(residues)?;
            let pairs = AlignmentPairCounts {
                total: choose_two(total_rows)?,
                missing: choose_two(total_rows)? - choose_two(nonmissing)?,
                excluded_ambiguity: choose_two(nonmissing)? - choose_two(retained)?,
                double_gap: choose_two(counts.gaps)?,
                single_gap: multiply(counts.gaps, residues)?,
                residue_pairs,
                identities,
                mismatches: residue_pairs - identities,
            };
            result.pairs.accumulate(pairs)?;
            result.column_statistics.push(AlignmentColumnStatistics {
                column,
                cells: counts,
                pairs,
            });
        }
        result.denominator = match options.denominator {
            IdentityDenominator::ResiduePairs => result.pairs.residue_pairs,
            IdentityDenominator::ResidueAndSingleGap => {
                add(result.pairs.residue_pairs, result.pairs.single_gap)?
            }
            IdentityDenominator::AllPairs => result.pairs.total,
        };
        result.identity = fraction(result.pairs.identities, result.denominator);
        Ok(result)
    }
}

fn gap_runs<A: AlignmentAnalysisAlphabet>(
    cells: &[AlignmentCell<A>],
) -> Result<AlignmentGapRuns, AlignmentAnalysisError> {
    let mut result = AlignmentGapRuns::default();
    let mut column = 0;
    while column < cells.len() {
        if !matches!(cells[column], AlignmentCell::Gap(_)) {
            column += 1;
            continue;
        }
        let start = column;
        while column < cells.len() && matches!(cells[column], AlignmentCell::Gap(_)) {
            column += 1;
        }
        if start == 0 {
            result.leading = add(result.leading, 1)?;
        } else if column == cells.len() {
            result.trailing = add(result.trailing, 1)?;
        } else {
            result.internal = add(result.internal, 1)?;
        }
    }
    Ok(result)
}

fn as_count(value: usize) -> Result<u64, AlignmentAnalysisError> {
    u64::try_from(value).map_err(|_| AlignmentAnalysisError::CountOverflow)
}

fn add(first: u64, second: u64) -> Result<u64, AlignmentAnalysisError> {
    first
        .checked_add(second)
        .ok_or(AlignmentAnalysisError::CountOverflow)
}

fn multiply(first: u64, second: u64) -> Result<u64, AlignmentAnalysisError> {
    first
        .checked_mul(second)
        .ok_or(AlignmentAnalysisError::CountOverflow)
}

fn choose_two(count: u64) -> Result<u64, AlignmentAnalysisError> {
    if count < 2 {
        return Ok(0);
    }
    if count % 2 == 0 {
        multiply(count / 2, count - 1)
    } else {
        multiply(count, (count - 1) / 2)
    }
}

#[allow(clippy::cast_precision_loss)]
fn fraction(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator != 0).then(|| numerator as f64 / denominator as f64)
}

/// Whether gap cells contribute votes to consensus.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ConsensusGapPolicy {
    /// Exclude gaps from the evidence denominator.
    #[default]
    Ignore,
    /// Count all gap glyphs as votes for the canonical output gap `-`.
    Vote,
    /// Refuse any gap cell.
    Reject,
}

/// Whether missing-data cells are ignored or refused for consensus.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ConsensusMissingPolicy {
    /// Exclude missing cells from the evidence denominator.
    #[default]
    Ignore,
    /// Refuse any missing cell instead of guessing its state.
    Reject,
}

/// How equal maximum consensus vote counts are resolved.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ConsensusTiePolicy {
    /// Emit N for DNA/RNA or X for protein and report the tie.
    #[default]
    Unknown,
    /// Select the lowest uppercase ASCII glyph, including `-` when gaps vote.
    AsciiOrder,
    /// Refuse the tied column and retain its tied glyphs in the error.
    Error,
}

/// Independent consensus policies and an inclusive minimum-support threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConsensusOptions {
    /// Whether to omit, count or reject gaps.
    pub gaps: ConsensusGapPolicy,
    /// Whether to omit or reject missing data.
    pub missing: ConsensusMissingPolicy,
    /// Whether to count literal, omit or reject ambiguity symbols.
    pub ambiguity: AlignmentAmbiguityPolicy,
    /// How equal maximum votes are resolved.
    pub ties: ConsensusTiePolicy,
    /// Finite minimum winning-count/evidence ratio, inclusive, in 0..=1.
    pub min_support: f64,
}

impl Default for ConsensusOptions {
    fn default() -> Self {
        Self {
            gaps: ConsensusGapPolicy::Ignore,
            missing: ConsensusMissingPolicy::Ignore,
            ambiguity: AlignmentAmbiguityPolicy::Literal,
            ties: ConsensusTiePolicy::Unknown,
            min_support: 0.5,
        }
    }
}

impl ConsensusOptions {
    /// Validates the threshold before any consensus computation.
    ///
    /// # Errors
    /// Returns `InvalidMinimumSupport` for NaN, infinity or values outside 0..=1.
    pub fn validate(self) -> Result<(), AlignmentAnalysisError> {
        if !self.min_support.is_finite() || !(0.0..=1.0).contains(&self.min_support) {
            return Err(AlignmentAnalysisError::InvalidMinimumSupport {
                value: self.min_support,
            });
        }
        Ok(())
    }
}

/// Why a consensus column emitted its selected glyph.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConsensusStatus {
    /// A unique residue winner meets the minimum support.
    Accepted,
    /// A selected gap winner meets the minimum support.
    Gap,
    /// Equal maximum votes were represented by the unknown symbol.
    Tie,
    /// Equal maximum votes were resolved by explicit ASCII order.
    TieResolved,
    /// Maximum support did not meet the inclusive minimum threshold.
    BelowThreshold,
    /// Every cell was ignored, or the rectangle has no rows.
    NoEvidence,
}

/// Complete evidence for one consensus column.
#[derive(Clone, Debug, PartialEq)]
pub struct ConsensusColumn {
    /// Zero-based column in the analyzed rectangle.
    pub column: usize,
    /// Emitted uppercase residue, unknown symbol or canonical gap `-`.
    pub symbol: u8,
    /// Reason for the emitted symbol.
    pub status: ConsensusStatus,
    /// Maximum vote count, zero when no evidence exists.
    pub winning_count: u64,
    /// Number of nonignored voting cells.
    pub denominator: u64,
    /// Maximum count divided by voting cells, absent with no evidence.
    pub support: Option<f64>,
    /// All maximum-vote glyphs in ASCII order; empty with no evidence.
    pub winners: Vec<u8>,
}

/// Gapped consensus text and exact per-column evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentConsensus {
    /// One output glyph per input column, including all-gap/no-evidence columns.
    pub sequence: String,
    /// Evidence and decision in column order.
    pub columns: Vec<ConsensusColumn>,
    /// Explicit policies used to obtain this result.
    pub options: ConsensusOptions,
}

impl<A: AlignmentAnalysisAlphabet> AlignmentAnalysis<A> {
    /// Computes consensus with explicit gap, missing, ambiguity and tie rules.
    ///
    /// Missing/gap/ambiguity refusal scans columns then rows. A tied column is
    /// handled before minimum-support testing: error/unknown tie policies remain
    /// observable even if the maximum support is below the threshold. ASCII tie
    /// resolution still must meet the threshold. No-evidence columns emit the
    /// canonical unknown even when the minimum support is zero.
    ///
    /// # Errors
    /// Returns `InvalidMinimumSupport`, `RejectedCell`, `ConsensusTie` or
    /// `CountOverflow`. No partial consensus is returned on failure.
    pub fn consensus(
        &self,
        options: ConsensusOptions,
    ) -> Result<AlignmentConsensus, AlignmentAnalysisError> {
        options.validate()?;
        let mut result = AlignmentConsensus {
            sequence: String::with_capacity(self.width),
            columns: Vec::with_capacity(self.width),
            options,
        };
        for column in 0..self.width {
            let mut votes = [0_u64; 128];
            let mut denominator = 0_u64;
            for (row_index, row) in self.rows.iter().enumerate() {
                let cell = row.cells[column];
                let vote = match cell {
                    AlignmentCell::Gap(byte) => match options.gaps {
                        ConsensusGapPolicy::Ignore => None,
                        ConsensusGapPolicy::Vote => Some(b'-'),
                        ConsensusGapPolicy::Reject => {
                            return Err(AlignmentAnalysisError::RejectedCell {
                                row: row_index,
                                column,
                                byte,
                                kind: "gap",
                            });
                        }
                    },
                    AlignmentCell::Missing(byte) => match options.missing {
                        ConsensusMissingPolicy::Ignore => None,
                        ConsensusMissingPolicy::Reject => {
                            return Err(AlignmentAnalysisError::RejectedCell {
                                row: row_index,
                                column,
                                byte,
                                kind: "missing",
                            });
                        }
                    },
                    AlignmentCell::Residue(symbol) => {
                        let byte = symbol.as_ascii();
                        if A::ambiguous(byte) {
                            match options.ambiguity {
                                AlignmentAmbiguityPolicy::Literal => Some(byte),
                                AlignmentAmbiguityPolicy::Ignore => None,
                                AlignmentAmbiguityPolicy::Reject => {
                                    return Err(AlignmentAnalysisError::RejectedCell {
                                        row: row_index,
                                        column,
                                        byte,
                                        kind: "ambiguity",
                                    });
                                }
                            }
                        } else {
                            Some(byte)
                        }
                    }
                };
                if let Some(byte) = vote {
                    votes[usize::from(byte)] = add(votes[usize::from(byte)], 1)?;
                    denominator = add(denominator, 1)?;
                }
            }
            let winning_count = votes.iter().copied().max().unwrap_or(0);
            let winners: Vec<_> = (0_u8..128)
                .filter(|byte| winning_count > 0 && votes[usize::from(*byte)] == winning_count)
                .collect();
            let support = fraction(winning_count, denominator);
            let tied = winners.len() > 1;
            if tied && options.ties == ConsensusTiePolicy::Error {
                return Err(AlignmentAnalysisError::ConsensusTie { column, winners });
            }
            let (symbol, status) = if denominator == 0 {
                (A::UNKNOWN, ConsensusStatus::NoEvidence)
            } else if tied && options.ties == ConsensusTiePolicy::Unknown {
                (A::UNKNOWN, ConsensusStatus::Tie)
            } else if support.is_some_and(|value| value < options.min_support) {
                (A::UNKNOWN, ConsensusStatus::BelowThreshold)
            } else {
                let symbol = winners[0];
                let status = if tied {
                    ConsensusStatus::TieResolved
                } else if symbol == b'-' {
                    ConsensusStatus::Gap
                } else {
                    ConsensusStatus::Accepted
                };
                (symbol, status)
            };
            result.sequence.push(char::from(symbol));
            result.columns.push(ConsensusColumn {
                column,
                symbol,
                status,
                winning_count,
                denominator,
                support,
                winners,
            });
        }
        Ok(result)
    }
}

/// A checked alignment-analysis operation failed without a partial result.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum AlignmentAnalysisError {
    /// A dimension exceeds its explicit constructor limit.
    LimitExceeded {
        /// `rows`, `columns` or `cells`.
        resource: &'static str,
        /// Inclusive configured resource limit.
        limit: usize,
        /// Observed resource count.
        observed: usize,
        /// Source/derived row count.
        rows: usize,
        /// Source/derived column count.
        columns: usize,
    },
    /// A dimension, pair count or aggregate cannot fit its integer type.
    CountOverflow,
    /// A requested range is reversed or out of bounds.
    InvalidRange {
        /// `rows` or `columns`.
        axis: &'static str,
        /// Inclusive requested start.
        start: usize,
        /// Exclusive requested end.
        end: usize,
        /// Available length on that axis.
        length: usize,
    },
    /// A requested row does not exist.
    InvalidRow {
        /// Requested zero-based row.
        row: usize,
        /// Number of available rows.
        rows: usize,
    },
    /// A requested alignment column does not exist.
    InvalidColumn {
        /// Requested zero-based column.
        column: usize,
        /// Number of available columns.
        columns: usize,
    },
    /// A requested position is beyond its complete original source.
    InvalidSourcePosition {
        /// Zero-based row.
        row: usize,
        /// Requested zero-based original source position.
        position: usize,
        /// Full source length.
        source_length: usize,
    },
    /// An origin is invalid or consumes a different number of residues.
    InvalidCoordinates {
        /// Row when attaching an origin, absent when constructing one.
        row: Option<usize>,
        /// Full declared source length.
        source_length: usize,
        /// Requested ascending source range.
        range: Range<usize>,
        /// Actual represented residues, when attached to a known row.
        represented_residues: Option<usize>,
    },
    /// Missing-data cells prevent inferring the row's original source positions.
    MissingCoordinates {
        /// Row whose mapping is unknown.
        row: usize,
    },
    /// A nontrivial transform would retain potentially stale opaque metadata.
    AnnotationLoss {
        /// Every affected annotation, in original order with complete metadata.
        annotations: Vec<AlignmentAnnotation>,
    },
    /// A pairwise result could not render against the supplied originals.
    InvalidPairwise {
        /// Original checked-render error message.
        message: String,
    },
    /// An adapter encountered an invalid stored biological byte.
    InvalidCell {
        /// Zero-based column.
        column: usize,
        /// Invalid ASCII byte.
        byte: u8,
    },
    /// A requested policy rejects a retained cell.
    RejectedCell {
        /// Zero-based row.
        row: usize,
        /// Zero-based column.
        column: usize,
        /// Original normalized/source glyph.
        byte: u8,
        /// `gap`, `missing` or `ambiguity`.
        kind: &'static str,
    },
    /// Consensus minimum support is not finite or lies outside 0..=1.
    InvalidMinimumSupport {
        /// Rejected threshold, which may be NaN or infinite.
        value: f64,
    },
    /// Error-policy consensus found equal maximum vote counts.
    ConsensusTie {
        /// Zero-based tied column.
        column: usize,
        /// Maximum-vote glyphs in stable ASCII order.
        winners: Vec<u8>,
    },
}

impl AlignmentAnalysisError {
    /// Returns a stable machine-readable category independent of display text.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::LimitExceeded { .. } => "limit_exceeded",
            Self::CountOverflow => "count_overflow",
            Self::InvalidRange { .. } => "invalid_range",
            Self::InvalidRow { .. } => "invalid_row",
            Self::InvalidColumn { .. } => "invalid_column",
            Self::InvalidSourcePosition { .. } => "invalid_source_position",
            Self::InvalidCoordinates { .. } => "invalid_coordinates",
            Self::MissingCoordinates { .. } => "missing_coordinates",
            Self::AnnotationLoss { .. } => "annotation_loss",
            Self::InvalidPairwise { .. } => "invalid_pairwise",
            Self::InvalidCell { .. } => "invalid_cell",
            Self::RejectedCell { .. } => "rejected_cell",
            Self::InvalidMinimumSupport { .. } => "invalid_minimum_support",
            Self::ConsensusTie { .. } => "consensus_tie",
        }
    }
}

impl fmt::Display for AlignmentAnalysisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                limit,
                observed,
                rows,
                columns,
            } => write!(
                formatter,
                "alignment {resource} count {observed} exceeds limit {limit} ({rows} rows, {columns} columns)"
            ),
            Self::CountOverflow => {
                formatter.write_str("alignment analysis count exceeds integer limits")
            }
            Self::InvalidRange {
                axis,
                start,
                end,
                length,
            } => write!(
                formatter,
                "{axis} range {start}..{end} is invalid for length {length}"
            ),
            Self::InvalidRow { row, rows } => {
                write!(formatter, "row {row} is outside alignment row count {rows}")
            }
            Self::InvalidColumn { column, columns } => {
                write!(
                    formatter,
                    "column {column} is outside alignment width {columns}"
                )
            }
            Self::InvalidSourcePosition {
                row,
                position,
                source_length,
            } => write!(
                formatter,
                "source position {position} is invalid for row {row} source length {source_length}"
            ),
            Self::InvalidCoordinates {
                row,
                source_length,
                range,
                represented_residues,
            } => write!(
                formatter,
                "invalid source coordinates {}..{} for source length {source_length}, row {row:?}, represented residues {represented_residues:?}",
                range.start, range.end
            ),
            Self::MissingCoordinates { row } => write!(
                formatter,
                "row {row} has no known residue-to-source mapping; missing data is not guessed"
            ),
            Self::AnnotationLoss { annotations } => write!(
                formatter,
                "slice could invalidate {} annotations; explicit discard permission is required",
                annotations.len()
            ),
            Self::InvalidPairwise { message } => {
                write!(formatter, "invalid pairwise input: {message}")
            }
            Self::InvalidCell { column, byte } => {
                write!(
                    formatter,
                    "invalid stored byte {byte} at alignment column {column}"
                )
            }
            Self::RejectedCell {
                row,
                column,
                byte,
                kind,
            } => write!(
                formatter,
                "policy rejects {kind} byte {byte} at row {row}, column {column}"
            ),
            Self::InvalidMinimumSupport { value } => write!(
                formatter,
                "consensus minimum support must be finite and in 0..=1; received {value}"
            ),
            Self::ConsensusTie { column, winners } => write!(
                formatter,
                "consensus tie at column {column} between ASCII glyphs {winners:?}"
            ),
        }
    }
}

impl std::error::Error for AlignmentAnalysisError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::alignment::{AlignmentFormat, AnnotationScope, parse_alignment};
    use align::prelude::{
        AlignedRow, AlignmentOperation, MafftStrategy, Multiple, Pairwise, PairwiseAlignment,
        SequenceRecord,
    };
    use std::fmt::Write as _;

    fn document<A: AlignmentAnalysisAlphabet>(rows: &[&str]) -> AlignmentDocument<A> {
        let mut input = String::from("# STOCKHOLM 1.0\n");
        if rows.iter().all(|row| row.is_empty()) {
            input.clear();
            for (index, row) in rows.iter().enumerate() {
                writeln!(input, ">r{index}\n{row}").unwrap();
            }
            return parse_alignment(input.as_bytes(), AlignmentFormat::Fasta).unwrap();
        }
        for (index, row) in rows.iter().enumerate() {
            writeln!(input, "r{index} {row}").unwrap();
        }
        input.push_str("//\n");
        parse_alignment(input.as_bytes(), AlignmentFormat::Stockholm).unwrap()
    }

    fn dna_analysis(rows: &[&str]) -> AlignmentAnalysis<Dna> {
        AlignmentAnalysis::from_document(&document(rows)).unwrap()
    }

    fn rendered<A: AlignmentAnalysisAlphabet>(analysis: &AlignmentAnalysis<A>) -> Vec<String> {
        analysis
            .rows()
            .iter()
            .map(AlignmentAnalysisRow::aligned)
            .collect()
    }

    #[test]
    fn imported_rectangle_retains_distinct_cells_and_original_glyphs() {
        let input = b"# STOCKHOLM 1.0\n#=GF ID retained\nr0 A.-?N\nr1 C--?G\n//\n";
        let source = parse_alignment::<Dna>(input.as_slice(), AlignmentFormat::Stockholm).unwrap();
        let analysis = AlignmentAnalysis::from_document(&source).unwrap();
        assert_eq!(analysis.width(), 5);
        assert_eq!(analysis.sequence_count(), 2);
        assert_eq!(rendered(&analysis), ["A.-?N", "C--?G"]);
        assert_eq!(analysis.rows()[0].identifier(), "r0");
        assert_eq!(analysis.rows()[0].cells()[1], AlignmentCell::Gap(b'.'));
        assert_eq!(analysis.rows()[0].cells()[2], AlignmentCell::Gap(b'-'));
        assert_eq!(analysis.rows()[0].cells()[3], AlignmentCell::Missing(b'?'));
        assert!(matches!(
            analysis.rows()[0].cells()[4],
            AlignmentCell::Residue(_)
        ));
        assert_eq!(analysis.annotations(), source.annotations());
        assert!(analysis.rows()[0].coordinates().is_none());
        assert_eq!(source.rows()[0].aligned(), "A.-?N");
    }

    #[test]
    fn nexus_custom_gap_and_missing_symbols_are_preserved() {
        let input = b"#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA GAP=~ MISSING=!;\nMATRIX\na A~!C\nb CG!C\n;\nEND;\n";
        let source = parse_alignment::<Dna>(input.as_slice(), AlignmentFormat::Nexus).unwrap();
        let analysis = AlignmentAnalysis::from_document(&source).unwrap();
        assert_eq!(rendered(&analysis), ["A~!C", "CG!C"]);
        assert_eq!(analysis.rows()[0].cells()[1], AlignmentCell::Gap(b'~'));
        assert_eq!(analysis.rows()[0].cells()[2], AlignmentCell::Missing(b'!'));
        let sliced = analysis.slice_columns(1..3).unwrap();
        assert_eq!(rendered(&sliced), ["~!", "G!"]);
    }

    #[test]
    fn default_limits_are_explicit_and_enforced_before_copying() {
        let source = document::<Dna>(&["AC-G", "ACTG"]);
        assert_eq!(
            AlignmentAnalysisLimits::default(),
            AlignmentAnalysisLimits {
                max_rows: 100_000,
                max_columns: 1_000_000,
                max_cells: 10_000_000,
            }
        );
        for (limits, resource, observed) in [
            (
                AlignmentAnalysisLimits {
                    max_rows: 1,
                    ..Default::default()
                },
                "rows",
                2,
            ),
            (
                AlignmentAnalysisLimits {
                    max_columns: 3,
                    ..Default::default()
                },
                "columns",
                4,
            ),
            (
                AlignmentAnalysisLimits {
                    max_cells: 7,
                    ..Default::default()
                },
                "cells",
                8,
            ),
        ] {
            assert!(matches!(
                AlignmentAnalysis::from_document_with_limits(&source, limits),
                Err(AlignmentAnalysisError::LimitExceeded {
                    resource: actual, observed: count, rows: 2, columns: 4, ..
                }) if actual == resource && count == observed
            ));
        }
        let limits = AlignmentAnalysisLimits {
            max_rows: 2,
            max_columns: 4,
            max_cells: 8,
        };
        let analysis = AlignmentAnalysis::from_document_with_limits(&source, limits).unwrap();
        assert_eq!(analysis.slice_columns(1..3).unwrap().limits(), limits);
        assert_eq!(source.width(), 4);
    }

    #[test]
    fn checked_dimension_and_pair_arithmetic_refuses_overflow() {
        let unlimited = AlignmentAnalysisLimits {
            max_rows: usize::MAX,
            max_columns: usize::MAX,
            max_cells: usize::MAX,
        };
        assert_eq!(
            unlimited.check(usize::MAX, 2),
            Err(AlignmentAnalysisError::CountOverflow)
        );
        assert_eq!(choose_two(0).unwrap(), 0);
        assert_eq!(choose_two(1).unwrap(), 0);
        assert_eq!(choose_two(6).unwrap(), 15);
        assert_eq!(
            choose_two(6_000_000_000).unwrap(),
            17_999_999_997_000_000_000
        );
        assert_eq!(
            choose_two(u64::MAX),
            Err(AlignmentAnalysisError::CountOverflow)
        );
        assert_eq!(add(u64::MAX, 1), Err(AlignmentAnalysisError::CountOverflow));
        assert_eq!(
            multiply(u64::MAX, 2),
            Err(AlignmentAnalysisError::CountOverflow)
        );
    }

    #[test]
    fn strict_slices_preserve_rows_all_gap_columns_and_source() {
        let analysis = dna_analysis(&["AC-G", "ACTG"]);
        let selected = analysis.slice_rows(0..1).unwrap();
        assert_eq!(rendered(&selected), ["AC-G"]);
        assert_eq!(selected.width(), 4);
        assert_eq!(selected.rows()[0].cells()[2], AlignmentCell::Gap(b'-'));
        assert_eq!(
            rendered(&analysis.slice_columns(1..3).unwrap()),
            ["C-", "CT"]
        );
        assert_eq!(rendered(&analysis), ["AC-G", "ACTG"]);
        assert_eq!(analysis.slice_columns(0..4).unwrap(), analysis);
        assert!(matches!(
            analysis.slice_columns(4..5),
            Err(AlignmentAnalysisError::InvalidRange {
                axis: "columns",
                ..
            })
        ));
        assert!(matches!(
            analysis.slice_rows(0..3),
            Err(AlignmentAnalysisError::InvalidRange { axis: "rows", .. })
        ));
        let reversed = Range { start: 3, end: 2 };
        assert!(matches!(
            analysis.slice_columns(reversed),
            Err(AlignmentAnalysisError::InvalidRange { .. })
        ));
    }

    #[test]
    fn empty_rectangles_retain_explicit_dimensions() {
        let source = dna_analysis(&["AC-G", "ACTG"]);
        let no_rows = source.slice_rows(1..1).unwrap();
        assert_eq!((no_rows.sequence_count(), no_rows.width()), (0, 4));
        assert_eq!(no_rows.slice_columns(1..3).unwrap().width(), 2);
        let counts = no_rows
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(counts.cells, AlignmentCellCounts::default());
        assert_eq!(counts.pairs, AlignmentPairCounts::default());
        assert_eq!(counts.identity, None);
        assert_eq!(counts.column_statistics.len(), 4);
        let consensus = no_rows.consensus(ConsensusOptions::default()).unwrap();
        assert_eq!(consensus.sequence, "NNNN");
        assert!(
            consensus
                .columns
                .iter()
                .all(|column| column.status == ConsensusStatus::NoEvidence)
        );
        let no_columns = source.slice_columns(3..3).unwrap();
        assert_eq!((no_columns.sequence_count(), no_columns.width()), (2, 0));
        assert_eq!(rendered(&no_columns), ["", ""]);
        assert_eq!(
            no_columns
                .consensus(ConsensusOptions::default())
                .unwrap()
                .sequence,
            ""
        );
        assert_eq!(no_columns.sequence_to_column(0, 0).unwrap(), None);
        let parsed_empty = dna_analysis(&["", ""]);
        assert_eq!(
            (parsed_empty.sequence_count(), parsed_empty.width()),
            (2, 0)
        );
        assert!(matches!(
            parsed_empty.sequence_to_column(0, 0),
            Err(AlignmentAnalysisError::InvalidSourcePosition {
                source_length: 0,
                ..
            })
        ));
    }

    #[test]
    fn original_ungapped_coordinate_mapping_and_bounds_are_checked() {
        let analysis = dna_analysis(&["AC-G", "ACTG"]);
        assert_eq!(analysis.rows()[0].coordinates().unwrap().range(), 0..3);
        assert_eq!(analysis.column_to_sequence(0, 0).unwrap(), Some(0));
        assert_eq!(analysis.column_to_sequence(0, 2).unwrap(), None);
        assert_eq!(analysis.column_to_sequence(0, 3).unwrap(), Some(2));
        assert_eq!(analysis.sequence_to_column(0, 2).unwrap(), Some(3));
        assert!(matches!(
            analysis.column_to_sequence(2, 0),
            Err(AlignmentAnalysisError::InvalidRow { row: 2, rows: 2 })
        ));
        assert!(matches!(
            analysis.column_to_sequence(0, 4),
            Err(AlignmentAnalysisError::InvalidColumn {
                column: 4,
                columns: 4
            })
        ));
        assert!(matches!(
            analysis.sequence_to_column(0, 3),
            Err(AlignmentAnalysisError::InvalidSourcePosition {
                source_length: 3,
                ..
            })
        ));
        let slice = analysis.slice_columns(1..4).unwrap();
        assert_eq!(slice.column_to_sequence(0, 0).unwrap(), Some(1));
        assert_eq!(slice.sequence_to_column(0, 0).unwrap(), None);
        assert_eq!(slice.rows()[0].coordinates().unwrap().range(), 1..3);
    }

    #[test]
    fn reverse_source_origins_survive_slicing_and_empty_ranges() {
        let analysis = dna_analysis(&["AC-GT"])
            .with_row_coordinates(0, RowCoordinates::new(20, 10..14, Strand::Reverse).unwrap())
            .unwrap();
        let actual: Vec<_> = (0..5)
            .map(|column| analysis.column_to_sequence(0, column).unwrap())
            .collect();
        assert_eq!(actual, [Some(13), Some(12), None, Some(11), Some(10)]);
        let selected = analysis.slice_columns(1..4).unwrap();
        assert_eq!(selected.rows()[0].coordinates().unwrap().range(), 11..13);
        assert_eq!(
            selected.rows()[0].coordinates().unwrap().strand(),
            Strand::Reverse
        );
        assert_eq!(
            selected.rows()[0].coordinates().unwrap().source_length(),
            20
        );
        assert_eq!(selected.column_to_sequence(0, 0).unwrap(), Some(12));
        assert_eq!(selected.column_to_sequence(0, 1).unwrap(), None);
        assert_eq!(selected.column_to_sequence(0, 2).unwrap(), Some(11));
        assert_eq!(selected.sequence_to_column(0, 13).unwrap(), None);
        let gap_only = analysis.slice_columns(2..3).unwrap();
        assert_eq!(gap_only.rows()[0].coordinates().unwrap().range(), 12..12);
        assert_eq!(gap_only.column_to_sequence(0, 0).unwrap(), None);
        for position in 0..20 {
            assert_eq!(gap_only.sequence_to_column(0, position).unwrap(), None);
        }
    }

    #[test]
    fn coordinate_origins_refuse_bad_ranges_lengths_and_missing_rows() {
        assert!(matches!(
            RowCoordinates::new(4, 3..5, Strand::Forward),
            Err(AlignmentAnalysisError::InvalidCoordinates { row: None, .. })
        ));
        assert!(RowCoordinates::new(4, Range { start: 3, end: 2 }, Strand::Reverse).is_err());
        let analysis = dna_analysis(&["AC-G"]);
        let origin = RowCoordinates::new(10, 4..6, Strand::Forward).unwrap();
        assert!(matches!(
            analysis.with_row_coordinates(0, origin),
            Err(AlignmentAnalysisError::InvalidCoordinates {
                row: Some(0),
                represented_residues: Some(3),
                ..
            })
        ));
        let missing = dna_analysis(&["A?CG"]);
        assert!(matches!(
            missing.column_to_sequence(0, 0),
            Err(AlignmentAnalysisError::MissingCoordinates { row: 0 })
        ));
        assert!(matches!(
            missing.sequence_to_column(0, 0),
            Err(AlignmentAnalysisError::MissingCoordinates { row: 0 })
        ));
        assert!(matches!(
            missing.with_row_coordinates(0, RowCoordinates::new(3, 0..3, Strand::Forward).unwrap()),
            Err(AlignmentAnalysisError::MissingCoordinates { row: 0 })
        ));
        let resolved = missing.slice_columns(2..4).unwrap();
        assert!(resolved.rows()[0].coordinates().is_none());
        let explicit = resolved
            .with_row_coordinates(0, RowCoordinates::new(10, 6..8, Strand::Forward).unwrap())
            .unwrap();
        assert_eq!(explicit.column_to_sequence(0, 0).unwrap(), Some(6));
        assert!(missing.rows()[0].coordinates().is_none());
    }

    #[test]
    fn pairwise_adapter_retains_local_origins_and_uses_actual_symbols() {
        let first = Sequence::<Dna>::from_ascii("TTACGCC").unwrap();
        let second = Sequence::<Dna>::from_ascii("GGACTGA").unwrap();
        // Structurally valid labels deliberately contradict the rendered symbols.
        let traceback = PairwiseAlignment::try_new(
            999,
            2..5,
            2..6,
            vec![
                AlignmentOperation::Substitution,
                AlignmentOperation::Match,
                AlignmentOperation::Insertion,
                AlignmentOperation::Substitution,
            ],
        )
        .unwrap();
        let analysis = AlignmentAnalysis::from_pairwise(&first, &second, &traceback).unwrap();
        assert_eq!(rendered(&analysis), ["AC-G", "ACTG"]);
        assert_eq!(analysis.rows()[0].identifier(), "first");
        assert_eq!(analysis.rows()[1].identifier(), "second");
        assert_eq!(analysis.rows()[0].coordinates().unwrap().range(), 2..5);
        assert_eq!(analysis.rows()[0].coordinates().unwrap().source_length(), 7);
        assert_eq!(analysis.sequence_to_column(0, 1).unwrap(), None);
        assert_eq!(analysis.sequence_to_column(0, 4).unwrap(), Some(3));
        let statistics = analysis
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(statistics.pairs.identities, 3);
        assert_eq!(statistics.pairs.mismatches, 0);
        assert_eq!(statistics.pairs.single_gap, 1);
        assert_eq!(statistics.identity, Some(1.0));
        assert_eq!(traceback.score(), 999);
    }

    #[test]
    fn pairwise_adapter_refuses_source_bounds_and_honors_limits() {
        let sequence = Sequence::<Dna>::from_ascii("AC").unwrap();
        let traceback =
            PairwiseAlignment::try_new(0, 3..4, 0..1, vec![AlignmentOperation::Match]).unwrap();
        assert!(matches!(
            AlignmentAnalysis::from_pairwise(&sequence, &sequence, &traceback),
            Err(AlignmentAnalysisError::InvalidPairwise { .. })
        ));
        let limits = AlignmentAnalysisLimits {
            max_rows: 1,
            ..Default::default()
        };
        assert!(matches!(
            AlignmentAnalysis::from_pairwise_with_limits(&sequence, &sequence, &traceback, limits),
            Err(AlignmentAnalysisError::LimitExceeded {
                resource: "rows",
                ..
            })
        ));
    }

    #[test]
    fn canonical_msa_and_result_adapters_share_symbols_and_origins() {
        let rows = vec![
            AlignedRow::new(
                "reference",
                vec![
                    Symbol::<Dna>::from_ascii(b'A'),
                    None,
                    Symbol::from_ascii(b'C'),
                ],
            )
            .unwrap(),
            AlignedRow::new(
                "sample",
                vec![
                    Symbol::<Dna>::from_ascii(b'A'),
                    Symbol::from_ascii(b'G'),
                    Symbol::from_ascii(b'C'),
                ],
            )
            .unwrap(),
        ];
        let source = MultipleAlignment::try_new(rows).unwrap();
        let direct = AlignmentAnalysis::from_multiple(&source).unwrap();
        let by_view = AlignmentAnalysis::from_multiple_result(&&source).unwrap();
        assert_eq!(direct, by_view);
        assert_eq!(rendered(&direct), ["A-C", "AGC"]);
        assert_eq!(direct.rows()[0].coordinates().unwrap().source_length(), 2);
        assert_eq!(direct.slice_rows(0..1).unwrap().width(), 3);
        let limits = AlignmentAnalysisLimits {
            max_cells: 5,
            ..Default::default()
        };
        assert!(matches!(
            AlignmentAnalysis::from_multiple_with_limits(&source, limits),
            Err(AlignmentAnalysisError::LimitExceeded {
                resource: "cells",
                ..
            })
        ));
        assert!(AlignmentAnalysis::from_multiple_result_with_limits(&&source, limits).is_err());
    }

    #[test]
    fn accepted_pairwise_and_tiny_msa_results_adapt_without_realigning() {
        let first = Sequence::<Dna>::from_ascii("TTACGCC").unwrap();
        let second = Sequence::<Dna>::from_ascii("GGACGA").unwrap();
        let run = Pairwise::smith_waterman().run(&first, &second).unwrap();
        let analysis = AlignmentAnalysis::from_pairwise(run.first(), run.second(), &run).unwrap();
        assert_eq!(
            analysis.rows()[0].coordinates().unwrap().range(),
            run.result().first_range()
        );
        assert_eq!(
            analysis.rows()[1].coordinates().unwrap().range(),
            run.result().second_range()
        );
        let rendered_pair = run.result().render(&first, &second).unwrap();
        assert_eq!(analysis.rows()[0].aligned(), rendered_pair.first());
        let records = vec![
            SequenceRecord::new("a", Sequence::<Dna>::from_ascii("ACGT").unwrap()).unwrap(),
            SequenceRecord::new("b", Sequence::<Dna>::from_ascii("AGT").unwrap()).unwrap(),
        ];
        let result = Multiple::mafft(MafftStrategy::FftNs1)
            .align(&records)
            .unwrap();
        let analysis = AlignmentAnalysis::from_multiple_result(&result).unwrap();
        for (adapted, original) in analysis.rows().iter().zip(result.alignment().rows()) {
            assert_eq!(adapted.identifier(), original.identifier());
            assert_eq!(adapted.aligned(), original.to_ascii_with_gaps());
        }
    }

    #[test]
    fn opaque_annotation_identity_preservation_and_discard_reports_are_exact() {
        let input = b"# STOCKHOLM 1.0\n#=GF ID source\n# comment with coordinates 1-4\n#=GS a DE original row\na A-CG\nb ATCG\n#=GC SS_cons <..>\n#=GR a SS <..>\n//\n";
        let source = parse_alignment::<Dna>(input.as_slice(), AlignmentFormat::Stockholm).unwrap();
        let analysis = AlignmentAnalysis::from_document(&source).unwrap();
        assert_eq!(analysis.annotations().len(), 5);
        assert_eq!(analysis.annotations()[0].scope(), AnnotationScope::File);
        assert_eq!(
            analysis.slice_columns(0..4).unwrap().annotations(),
            source.annotations()
        );
        let identity = analysis
            .slice_with(0..2, 0..4, AlignmentAnnotationPolicy::Discard)
            .unwrap();
        assert!(identity.discarded_annotations().is_empty());
        assert_eq!(identity.alignment().annotations(), source.annotations());
        for error in [
            analysis.slice_rows(0..1).unwrap_err(),
            analysis.slice_columns(0..2).unwrap_err(),
        ] {
            match error {
                AlignmentAnalysisError::AnnotationLoss { annotations } => {
                    assert_eq!(annotations, source.annotations());
                }
                other => panic!("unexpected error {other}"),
            }
        }
        let slice = analysis
            .slice_with(0..1, 1..3, AlignmentAnnotationPolicy::Discard)
            .unwrap();
        assert_eq!(rendered(slice.alignment()), ["-C"]);
        assert!(slice.alignment().annotations().is_empty());
        assert_eq!(slice.discarded_annotations(), source.annotations());
        assert_eq!(slice.discarded_annotations()[4].sequence(), Some("a"));
        assert_eq!(slice.discarded_annotations()[4].tag(), "SS");
        assert_eq!(slice.discarded_annotations()[4].value(), "<..>");
        assert_eq!(slice.discarded_annotations()[4].line(), 8);
        let (derivative, report) = slice.into_parts();
        assert_eq!(derivative.width(), 2);
        assert_eq!(report, source.annotations());
        assert_eq!(analysis.annotations(), source.annotations());
    }

    #[test]
    fn exact_statistics_cover_missing_ambiguity_gaps_and_all_denominators() {
        let analysis = dna_analysis(&["A-N?", "A--A", "C--A", "N.-A"]);
        let literal = analysis
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(
            literal.cells,
            AlignmentCellCounts {
                residues: 8,
                gaps: 7,
                missing: 1,
                ambiguous: 2
            }
        );
        assert_eq!(literal.pairs.total, 24);
        assert_eq!(literal.pairs.missing, 3);
        assert_eq!(literal.pairs.excluded_ambiguity, 0);
        assert_eq!(literal.pairs.double_gap, 9);
        assert_eq!(literal.pairs.single_gap, 3);
        assert_eq!(literal.pairs.residue_pairs, 9);
        assert_eq!(literal.pairs.identities, 4);
        assert_eq!(literal.pairs.mismatches, 5);
        assert_eq!(literal.denominator, 9);
        assert_eq!(literal.identity, Some(4.0 / 9.0));
        let ignored = analysis
            .statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(ignored.pairs.excluded_ambiguity, 6);
        assert_eq!(ignored.pairs.residue_pairs, 6);
        assert_eq!(ignored.pairs.single_gap, 0);
        assert_eq!(ignored.pairs.identities, 4);
        assert_eq!(ignored.pairs.mismatches, 2);
        assert_eq!(ignored.denominator, 6);
        let with_gaps = analysis
            .statistics(AlignmentStatisticsOptions {
                denominator: IdentityDenominator::ResidueAndSingleGap,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(with_gaps.denominator, 12);
        let all = analysis
            .statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                denominator: IdentityDenominator::AllPairs,
            })
            .unwrap();
        assert_eq!(all.denominator, 24);
        assert_eq!(all.identity, Some(1.0 / 6.0));
        assert_eq!(all.cells, literal.cells);
        assert_eq!(literal.column_statistics[1].cells.gaps, 4);
        assert_eq!(literal.column_statistics[1].pairs.double_gap, 6);
    }

    #[test]
    fn missing_pair_priority_precedes_excluded_ambiguity() {
        let analysis = dna_analysis(&["?", "N", "-"]);
        let counts = analysis
            .statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            counts.pairs,
            AlignmentPairCounts {
                total: 3,
                missing: 2,
                excluded_ambiguity: 1,
                ..Default::default()
            }
        );
        assert_eq!(counts.identity, None);
    }

    #[test]
    fn gap_runs_partition_all_gap_missing_and_boundary_cases() {
        let analysis = dna_analysis(&["--A--?--C--", "-----------", "ACGTACGTACG"]);
        let result = analysis
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(
            result.row_statistics[0].gap_runs,
            AlignmentGapRuns {
                leading: 1,
                internal: 2,
                trailing: 1
            }
        );
        assert_eq!(
            result.row_statistics[1].gap_runs,
            AlignmentGapRuns {
                leading: 1,
                internal: 0,
                trailing: 0
            }
        );
        assert_eq!(
            result.row_statistics[2].gap_runs,
            AlignmentGapRuns::default()
        );
        assert_eq!(
            result.gap_runs,
            AlignmentGapRuns {
                leading: 2,
                internal: 2,
                trailing: 1
            }
        );
        assert_eq!(result.row_statistics[0].cells.missing, 1);
    }

    #[test]
    fn literal_ambiguity_is_not_resolved_as_compatible_identity() {
        let analysis = dna_analysis(&["RN", "AA"]);
        let result = analysis
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(result.pairs.identities, 0);
        assert_eq!(result.pairs.mismatches, 2);
        let ignored = analysis
            .statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(ignored.pairs.excluded_ambiguity, 2);
        assert_eq!(ignored.identity, None);
        assert!(matches!(
            analysis.statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Reject,
                ..Default::default()
            }),
            Err(AlignmentAnalysisError::RejectedCell {
                row: 0,
                column: 0,
                byte: b'R',
                kind: "ambiguity"
            })
        ));
    }

    #[test]
    fn single_row_and_all_gap_identity_never_invent_evidence() {
        for rows in [&["ACGT"][..], &["----", "...."][..], &["????", "????"][..]] {
            let analysis = dna_analysis(rows);
            let result = analysis
                .statistics(AlignmentStatisticsOptions::default())
                .unwrap();
            assert_eq!(result.pairs.identities, 0);
            assert_eq!(result.denominator, 0);
            assert_eq!(result.identity, None);
        }
    }

    #[test]
    fn default_consensus_has_literal_votes_and_explicit_column_evidence() {
        let analysis = dna_analysis(&["AA-N?", "AC-A?", "AC-A?"]);
        let result = analysis.consensus(ConsensusOptions::default()).unwrap();
        assert_eq!(result.sequence, "ACNAN");
        assert_eq!(result.columns[0].winning_count, 3);
        assert_eq!(result.columns[0].denominator, 3);
        assert_eq!(result.columns[0].status, ConsensusStatus::Accepted);
        assert_eq!(result.columns[1].support, Some(2.0 / 3.0));
        assert_eq!(result.columns[2].status, ConsensusStatus::NoEvidence);
        assert_eq!(result.columns[2].support, None);
        assert!(result.columns[2].winners.is_empty());
        assert_eq!(result.columns[3].denominator, 3);
        assert_eq!(result.columns[3].winners, [b'A']);
        assert_eq!(result.columns[4].status, ConsensusStatus::NoEvidence);
    }

    #[test]
    fn consensus_tie_policies_and_threshold_are_distinguishable() {
        let analysis = dna_analysis(&["AC", "CA"]);
        let unknown = analysis.consensus(ConsensusOptions::default()).unwrap();
        assert_eq!(unknown.sequence, "NN");
        assert_eq!(unknown.columns[0].status, ConsensusStatus::Tie);
        assert_eq!(unknown.columns[0].winners, b"AC");
        let ordered = analysis
            .consensus(ConsensusOptions {
                ties: ConsensusTiePolicy::AsciiOrder,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(ordered.sequence, "AA");
        assert_eq!(ordered.columns[0].status, ConsensusStatus::TieResolved);
        let below = analysis
            .consensus(ConsensusOptions {
                ties: ConsensusTiePolicy::AsciiOrder,
                min_support: 0.500_001,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(below.sequence, "NN");
        assert_eq!(below.columns[0].status, ConsensusStatus::BelowThreshold);
        assert!(
            matches!(analysis.consensus(ConsensusOptions { ties: ConsensusTiePolicy::Error, min_support: 1.0, ..Default::default() }), Err(AlignmentAnalysisError::ConsensusTie { column: 0, winners }) if winners == b"AC")
        );
        let unique = dna_analysis(&["A", "A", "C"]);
        assert_eq!(
            unique
                .consensus(ConsensusOptions {
                    min_support: 2.0 / 3.0,
                    ..Default::default()
                })
                .unwrap()
                .sequence,
            "A"
        );
        assert_eq!(
            unique
                .consensus(ConsensusOptions {
                    min_support: 0.7,
                    ..Default::default()
                })
                .unwrap()
                .columns[0]
                .status,
            ConsensusStatus::BelowThreshold
        );
    }

    #[test]
    fn gap_votes_normalize_source_glyphs_and_expose_gap_ties() {
        let analysis = dna_analysis(&["-.", ".-", "AA"]);
        let voted = analysis
            .consensus(ConsensusOptions {
                gaps: ConsensusGapPolicy::Vote,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(voted.sequence, "--");
        assert_eq!(voted.columns[0].status, ConsensusStatus::Gap);
        assert_eq!(voted.columns[0].winning_count, 2);
        assert_eq!(voted.columns[0].denominator, 3);
        assert_eq!(
            analysis
                .consensus(ConsensusOptions::default())
                .unwrap()
                .sequence,
            "AA"
        );
        let tie = dna_analysis(&["A", "."]);
        let result = tie
            .consensus(ConsensusOptions {
                gaps: ConsensusGapPolicy::Vote,
                ties: ConsensusTiePolicy::AsciiOrder,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.sequence, "-");
        assert_eq!(result.columns[0].status, ConsensusStatus::TieResolved);
        assert_eq!(result.columns[0].winners, b"-A");
        let gaps_only = dna_analysis(&["--", ".."]);
        assert_eq!(
            gaps_only
                .consensus(ConsensusOptions {
                    gaps: ConsensusGapPolicy::Vote,
                    ..Default::default()
                })
                .unwrap()
                .sequence,
            "--"
        );
    }

    #[test]
    fn consensus_ignored_ambiguity_and_missing_reduce_denominators() {
        let analysis = dna_analysis(&["AN?", "NN?", "?A?"]);
        let result = analysis
            .consensus(ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.sequence, "AAN");
        assert_eq!(result.columns[0].denominator, 1);
        assert_eq!(result.columns[1].denominator, 1);
        assert_eq!(result.columns[2].denominator, 0);
        assert_eq!(result.columns[2].status, ConsensusStatus::NoEvidence);
        let literal = analysis.consensus(ConsensusOptions::default()).unwrap();
        assert_eq!(literal.sequence, "NNN");
        assert_eq!(literal.columns[0].status, ConsensusStatus::Tie);
        assert_eq!(literal.columns[1].status, ConsensusStatus::Accepted);
        assert_eq!(literal.columns[1].winning_count, 2);
    }

    #[test]
    fn each_consensus_rejection_reports_original_cell_and_traversal() {
        let analysis = dna_analysis(&["A?N", ".AA"]);
        assert!(matches!(
            analysis.consensus(ConsensusOptions {
                gaps: ConsensusGapPolicy::Reject,
                ..Default::default()
            }),
            Err(AlignmentAnalysisError::RejectedCell {
                row: 1,
                column: 0,
                byte: b'.',
                kind: "gap"
            })
        ));
        assert!(matches!(
            analysis.consensus(ConsensusOptions {
                missing: ConsensusMissingPolicy::Reject,
                ..Default::default()
            }),
            Err(AlignmentAnalysisError::RejectedCell {
                row: 0,
                column: 1,
                byte: b'?',
                kind: "missing"
            })
        ));
        assert!(matches!(
            analysis.consensus(ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Reject,
                ..Default::default()
            }),
            Err(AlignmentAnalysisError::RejectedCell {
                row: 0,
                column: 2,
                byte: b'N',
                kind: "ambiguity"
            })
        ));
    }

    #[test]
    fn consensus_invalid_threshold_precedes_cell_rejection_even_for_empty_input() {
        let analysis = dna_analysis(&["N"]);
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
            let options = ConsensusOptions {
                min_support: value,
                ambiguity: AlignmentAmbiguityPolicy::Reject,
                ..Default::default()
            };
            assert!(matches!(
                analysis.consensus(options),
                Err(AlignmentAnalysisError::InvalidMinimumSupport { .. })
            ));
            assert!(
                analysis
                    .slice_rows(0..0)
                    .unwrap()
                    .consensus(options)
                    .is_err()
            );
        }
        for value in [0.0, -0.0, 1.0] {
            assert!(
                ConsensusOptions {
                    min_support: value,
                    ..Default::default()
                }
                .validate()
                .is_ok()
            );
        }
        let no_evidence = dna_analysis(&["-"])
            .consensus(ConsensusOptions {
                min_support: 0.0,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(no_evidence.sequence, "N");
        assert_eq!(no_evidence.columns[0].status, ConsensusStatus::NoEvidence);
    }

    #[test]
    fn rna_and_extended_protein_analysis_have_explicit_canonical_subsets() {
        let rna = AlignmentAnalysis::from_document(&document::<Rna>(&["AUN", "AUA"])).unwrap();
        let rna_stats = rna
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(rna_stats.cells.ambiguous, 1);
        assert_eq!(
            rna.consensus(ConsensusOptions::default()).unwrap().sequence,
            "AUN"
        );
        let protein =
            AlignmentAnalysis::from_document(&document::<Protein>(&["ABZJXUO", "ACDEFGH"]))
                .unwrap();
        let protein_stats = protein
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap();
        assert_eq!(protein_stats.cells.ambiguous, 6);
        assert_eq!(protein_stats.pairs.identities, 1);
        assert_eq!(
            protein
                .consensus(ConsensusOptions::default())
                .unwrap()
                .sequence,
            "AXXXXXX"
        );
        let excluded = protein
            .statistics(AlignmentStatisticsOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(excluded.pairs.excluded_ambiguity, 6);
        assert_eq!(excluded.denominator, 1);
        let literal =
            AlignmentAnalysis::from_document(&document::<Protein>(&["UU", "OO"])).unwrap();
        assert_eq!(
            literal
                .consensus(ConsensusOptions::default())
                .unwrap()
                .sequence,
            "XX"
        );
    }

    fn explicit_pair_oracle(cells: &[u8], policy: AlignmentAmbiguityPolicy) -> AlignmentPairCounts {
        let mut result = AlignmentPairCounts::default();
        for left in 0..cells.len() {
            for right in left + 1..cells.len() {
                let (a, b) = (cells[left], cells[right]);
                result.total += 1;
                if a == b'?' || b == b'?' {
                    result.missing += 1;
                } else if policy == AlignmentAmbiguityPolicy::Ignore && (a == b'N' || b == b'N') {
                    result.excluded_ambiguity += 1;
                } else if b"-.".contains(&a) && b"-.".contains(&b) {
                    result.double_gap += 1;
                } else if b"-.".contains(&a) || b"-.".contains(&b) {
                    result.single_gap += 1;
                } else {
                    result.residue_pairs += 1;
                    if a == b {
                        result.identities += 1;
                    } else {
                        result.mismatches += 1;
                    }
                }
            }
        }
        result
    }

    #[test]
    fn exhaustive_frequency_counts_match_an_independent_explicit_pair_oracle() {
        let alphabet = b"ACN-.?";
        let mut comparisons = 0;
        for length in 1..=5_u32 {
            for encoded in 0..alphabet.len().pow(length) {
                let mut value = encoded;
                let bytes: Vec<_> = (0..length)
                    .map(|_| {
                        let byte = alphabet[value % alphabet.len()];
                        value /= alphabet.len();
                        byte
                    })
                    .collect();
                let rows: Vec<_> = bytes
                    .iter()
                    .map(|byte| char::from(*byte).to_string())
                    .collect();
                let refs: Vec<_> = rows.iter().map(String::as_str).collect();
                let analysis = dna_analysis(&refs);
                for ambiguity in [
                    AlignmentAmbiguityPolicy::Literal,
                    AlignmentAmbiguityPolicy::Ignore,
                ] {
                    let expected = explicit_pair_oracle(&bytes, ambiguity);
                    for denominator in [
                        IdentityDenominator::ResiduePairs,
                        IdentityDenominator::ResidueAndSingleGap,
                        IdentityDenominator::AllPairs,
                    ] {
                        let actual = analysis
                            .statistics(AlignmentStatisticsOptions {
                                ambiguity,
                                denominator,
                            })
                            .unwrap();
                        assert_eq!(
                            actual.pairs, expected,
                            "cells={bytes:?}, policy={ambiguity:?}"
                        );
                        let expected_denominator = match denominator {
                            IdentityDenominator::ResiduePairs => expected.residue_pairs,
                            IdentityDenominator::ResidueAndSingleGap => {
                                expected.residue_pairs + expected.single_gap
                            }
                            IdentityDenominator::AllPairs => expected.total,
                        };
                        assert_eq!(actual.denominator, expected_denominator);
                        assert_eq!(
                            actual.identity,
                            fraction(expected.identities, expected_denominator)
                        );
                        assert_eq!(
                            expected.total,
                            expected.missing
                                + expected.excluded_ambiguity
                                + expected.double_gap
                                + expected.single_gap
                                + expected.residue_pairs
                        );
                        comparisons += 1;
                    }
                }
            }
        }
        assert_eq!(comparisons, 55_980);
    }

    #[test]
    fn exhaustive_gapped_coordinate_inverses_and_slice_composition() {
        let alphabet = b"AC-";
        let mut rectangles = 0;
        for length in 0..=6_u32 {
            for encoded in 0..alphabet.len().pow(length) {
                let mut value = encoded;
                let bytes: Vec<_> = (0..length)
                    .map(|_| {
                        let byte = alphabet[value % alphabet.len()];
                        value /= alphabet.len();
                        byte
                    })
                    .collect();
                let text = String::from_utf8(bytes.clone()).unwrap();
                let residues = bytes.iter().filter(|&&byte| byte != b'-').count();
                for strand in [Strand::Forward, Strand::Reverse] {
                    let analysis = dna_analysis(&[&text])
                        .with_row_coordinates(
                            0,
                            RowCoordinates::new(20, 7..7 + residues, strand).unwrap(),
                        )
                        .unwrap();
                    let mut ordinal = 0;
                    for (column, &byte) in bytes.iter().enumerate() {
                        if byte == b'-' {
                            assert_eq!(analysis.column_to_sequence(0, column).unwrap(), None);
                        } else {
                            let expected = if strand == Strand::Forward {
                                7 + ordinal
                            } else {
                                6 + residues - ordinal
                            };
                            assert_eq!(
                                analysis.column_to_sequence(0, column).unwrap(),
                                Some(expected)
                            );
                            assert_eq!(
                                analysis.sequence_to_column(0, expected).unwrap(),
                                Some(column)
                            );
                            ordinal += 1;
                        }
                    }
                    for start in 0..=bytes.len() {
                        for end in start..=bytes.len() {
                            let slice = analysis.slice_columns(start..end).unwrap();
                            for offset in 0..end - start {
                                assert_eq!(
                                    slice.column_to_sequence(0, offset).unwrap(),
                                    analysis.column_to_sequence(0, start + offset).unwrap()
                                );
                            }
                            for inner in 0..=end - start {
                                assert_eq!(
                                    slice.slice_columns(inner..end - start).unwrap(),
                                    analysis.slice_columns(start + inner..end).unwrap()
                                );
                            }
                        }
                    }
                    rectangles += 1;
                }
            }
        }
        assert_eq!(rectangles, 2_186);
    }

    #[test]
    fn row_permutation_preserves_aggregate_statistics_and_consensus_votes() {
        let first = dna_analysis(&["AN-?", "AC-A", "CCNA"]);
        let second = dna_analysis(&["CCNA", "AN-?", "AC-A"]);
        for ambiguity in [
            AlignmentAmbiguityPolicy::Literal,
            AlignmentAmbiguityPolicy::Ignore,
        ] {
            let options = AlignmentStatisticsOptions {
                ambiguity,
                ..Default::default()
            };
            let a = first.statistics(options).unwrap();
            let b = second.statistics(options).unwrap();
            assert_eq!(a.pairs, b.pairs);
            assert_eq!(a.column_statistics, b.column_statistics);
            assert_eq!(a.cells, b.cells);
            assert_eq!(a.identity, b.identity);
            let options = ConsensusOptions {
                ambiguity,
                gaps: ConsensusGapPolicy::Vote,
                ..Default::default()
            };
            assert_eq!(
                first.consensus(options).unwrap(),
                second.consensus(options).unwrap()
            );
        }
    }

    fn sorted_ballots(
        cells: &[u8],
        options: ConsensusOptions,
    ) -> Result<Vec<u8>, AlignmentAnalysisError> {
        let mut ballots = Vec::new();
        for (row, &byte) in cells.iter().enumerate() {
            let (reject, retain, normalized, kind) = match byte {
                b'-' | b'.' => (
                    options.gaps == ConsensusGapPolicy::Reject,
                    options.gaps == ConsensusGapPolicy::Vote,
                    b'-',
                    "gap",
                ),
                b'?' => (
                    options.missing == ConsensusMissingPolicy::Reject,
                    false,
                    byte,
                    "missing",
                ),
                b'N' => (
                    options.ambiguity == AlignmentAmbiguityPolicy::Reject,
                    options.ambiguity == AlignmentAmbiguityPolicy::Literal,
                    byte,
                    "ambiguity",
                ),
                _ => (false, true, byte, "residue"),
            };
            if reject {
                return Err(AlignmentAnalysisError::RejectedCell {
                    row,
                    column: 0,
                    byte,
                    kind,
                });
            }
            if retain {
                ballots.push(normalized);
            }
        }
        ballots.sort_unstable();
        Ok(ballots)
    }

    fn sorted_consensus_oracle(
        cells: &[u8],
        options: ConsensusOptions,
    ) -> Result<ConsensusColumn, AlignmentAnalysisError> {
        let ballots = sorted_ballots(cells, options)?;
        let groups: Vec<_> = ballots.chunk_by(|first, second| first == second).collect();
        let winning_count = groups.iter().map(|group| group.len()).max().unwrap_or(0) as u64;
        let winners: Vec<_> = groups
            .iter()
            .filter(|group| group.len() as u64 == winning_count)
            .map(|group| group[0])
            .collect();
        let denominator = ballots.len() as u64;
        let support = fraction(winning_count, denominator);
        let (symbol, status) = if ballots.is_empty() {
            (b'N', ConsensusStatus::NoEvidence)
        } else if winners.len() > 1 && options.ties != ConsensusTiePolicy::AsciiOrder {
            if options.ties == ConsensusTiePolicy::Error {
                return Err(AlignmentAnalysisError::ConsensusTie { column: 0, winners });
            }
            (b'N', ConsensusStatus::Tie)
        } else if support.unwrap() < options.min_support {
            (b'N', ConsensusStatus::BelowThreshold)
        } else if winners.len() > 1 {
            (winners[0], ConsensusStatus::TieResolved)
        } else if winners[0] == b'-' {
            (b'-', ConsensusStatus::Gap)
        } else {
            (winners[0], ConsensusStatus::Accepted)
        };
        Ok(ConsensusColumn {
            column: 0,
            symbol,
            status,
            winning_count,
            denominator,
            support,
            winners,
        })
    }

    #[test]
    fn exhaustive_consensus_policy_cross_product_matches_sorted_ballot_oracle() {
        let alphabet = b"ACN-.?";
        let mut cases = 0;
        for length in 0..=4_u32 {
            for encoded in 0..alphabet.len().pow(length) {
                let mut value = encoded;
                let bytes: Vec<_> = (0..length)
                    .map(|_| {
                        let byte = alphabet[value % alphabet.len()];
                        value /= alphabet.len();
                        byte
                    })
                    .collect();
                let rows: Vec<_> = bytes
                    .iter()
                    .map(|byte| char::from(*byte).to_string())
                    .collect();
                let refs: Vec<_> = rows.iter().map(String::as_str).collect();
                let analysis = if length == 0 {
                    dna_analysis(&["A"]).slice_rows(0..0).unwrap()
                } else {
                    dna_analysis(&refs)
                };
                for gaps in [
                    ConsensusGapPolicy::Ignore,
                    ConsensusGapPolicy::Vote,
                    ConsensusGapPolicy::Reject,
                ] {
                    for missing in [
                        ConsensusMissingPolicy::Ignore,
                        ConsensusMissingPolicy::Reject,
                    ] {
                        for ambiguity in [
                            AlignmentAmbiguityPolicy::Literal,
                            AlignmentAmbiguityPolicy::Ignore,
                            AlignmentAmbiguityPolicy::Reject,
                        ] {
                            for ties in [
                                ConsensusTiePolicy::Unknown,
                                ConsensusTiePolicy::AsciiOrder,
                                ConsensusTiePolicy::Error,
                            ] {
                                for min_support in [0.0, 0.5, 1.0] {
                                    let options = ConsensusOptions {
                                        gaps,
                                        missing,
                                        ambiguity,
                                        ties,
                                        min_support,
                                    };
                                    let actual = analysis.consensus(options).map(|mut result| {
                                        let column = result.columns.remove(0);
                                        assert_eq!(result.sequence.as_bytes(), [column.symbol]);
                                        column
                                    });
                                    assert_eq!(
                                        actual,
                                        sorted_consensus_oracle(&bytes, options),
                                        "cells={bytes:?}, options={options:?}"
                                    );
                                    cases += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(cases, 251_910);
    }

    #[test]
    fn every_typed_error_has_a_stable_code_and_nonempty_display() {
        let cases = [
            (
                AlignmentAnalysisError::LimitExceeded {
                    resource: "cells",
                    limit: 0,
                    observed: 1,
                    rows: 1,
                    columns: 1,
                },
                "limit_exceeded",
            ),
            (AlignmentAnalysisError::CountOverflow, "count_overflow"),
            (
                AlignmentAnalysisError::InvalidRange {
                    axis: "rows",
                    start: 1,
                    end: 0,
                    length: 1,
                },
                "invalid_range",
            ),
            (
                AlignmentAnalysisError::InvalidRow { row: 1, rows: 0 },
                "invalid_row",
            ),
            (
                AlignmentAnalysisError::InvalidColumn {
                    column: 1,
                    columns: 0,
                },
                "invalid_column",
            ),
            (
                AlignmentAnalysisError::InvalidSourcePosition {
                    row: 0,
                    position: 1,
                    source_length: 0,
                },
                "invalid_source_position",
            ),
            (
                AlignmentAnalysisError::InvalidCoordinates {
                    row: None,
                    source_length: 0,
                    range: 0..1,
                    represented_residues: None,
                },
                "invalid_coordinates",
            ),
            (
                AlignmentAnalysisError::MissingCoordinates { row: 0 },
                "missing_coordinates",
            ),
            (
                AlignmentAnalysisError::AnnotationLoss {
                    annotations: Vec::new(),
                },
                "annotation_loss",
            ),
            (
                AlignmentAnalysisError::InvalidPairwise {
                    message: "invalid source".to_owned(),
                },
                "invalid_pairwise",
            ),
            (
                AlignmentAnalysisError::InvalidCell {
                    column: 0,
                    byte: b'*',
                },
                "invalid_cell",
            ),
            (
                AlignmentAnalysisError::RejectedCell {
                    row: 0,
                    column: 0,
                    byte: b'N',
                    kind: "ambiguity",
                },
                "rejected_cell",
            ),
            (
                AlignmentAnalysisError::InvalidMinimumSupport { value: 2.0 },
                "invalid_minimum_support",
            ),
            (
                AlignmentAnalysisError::ConsensusTie {
                    column: 0,
                    winners: b"AC".to_vec(),
                },
                "consensus_tie",
            ),
        ];
        for (error, expected) in cases {
            assert_eq!(error.code(), expected);
            assert!(!error.to_string().is_empty());
            assert!(std::error::Error::source(&error).is_none());
        }
    }

    #[test]
    fn adapter_cell_validation_has_exact_byte_and_column() {
        assert!(matches!(
            parse_cells::<Dna>("A*", b"-", None),
            Err(AlignmentAnalysisError::InvalidCell {
                column: 1,
                byte: b'*'
            })
        ));
    }
}
