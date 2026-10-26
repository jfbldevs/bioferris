use core::mem::size_of;

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer, Symbol,
};

use crate::AlignmentError;
use crate::affine_simd::{AffineBackend, AffineKernelPolicy, SimdBoundaryMode, try_fill_simd};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub(crate) enum BestTrace {
    Start = 0,
    Diagonal = 1,
    Deletion = 2,
    Insertion = 3,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum GapTrace {
    Open,
    Extend,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(transparent)]
pub(crate) struct TraceCell(u8);

impl TraceCell {
    const BEST_MASK: u8 = 0b0000_0011;
    const DELETION_EXTENDS: u8 = 0b0000_0100;
    const INSERTION_EXTENDS: u8 = 0b0000_1000;

    pub(crate) const fn from_bits(bits: u8) -> Self {
        Self(bits)
    }

    pub(crate) const fn bits(self) -> u8 {
        self.0
    }

    pub(crate) const fn with_best(mut self, trace: BestTrace) -> Self {
        self.0 = (self.0 & !Self::BEST_MASK) | trace as u8;
        self
    }

    pub(crate) const fn with_deletion(mut self, trace: GapTrace) -> Self {
        if matches!(trace, GapTrace::Extend) {
            self.0 |= Self::DELETION_EXTENDS;
        } else {
            self.0 &= !Self::DELETION_EXTENDS;
        }
        self
    }

    pub(crate) const fn with_insertion(mut self, trace: GapTrace) -> Self {
        if matches!(trace, GapTrace::Extend) {
            self.0 |= Self::INSERTION_EXTENDS;
        } else {
            self.0 &= !Self::INSERTION_EXTENDS;
        }
        self
    }

    pub(crate) const fn from_traces(
        best: BestTrace,
        deletion: GapTrace,
        insertion: GapTrace,
    ) -> Self {
        let mut bits = best as u8;
        if matches!(deletion, GapTrace::Extend) {
            bits |= Self::DELETION_EXTENDS;
        }
        if matches!(insertion, GapTrace::Extend) {
            bits |= Self::INSERTION_EXTENDS;
        }
        Self(bits)
    }

    pub(crate) const fn best(self) -> Option<BestTrace> {
        match self.0 & Self::BEST_MASK {
            0 => Some(BestTrace::Start),
            1 => Some(BestTrace::Diagonal),
            2 => Some(BestTrace::Deletion),
            3 => Some(BestTrace::Insertion),
            _ => None,
        }
    }

    pub(crate) const fn deletion(self) -> GapTrace {
        if self.0 & Self::DELETION_EXTENDS == 0 {
            GapTrace::Open
        } else {
            GapTrace::Extend
        }
    }

    pub(crate) const fn insertion(self) -> GapTrace {
        if self.0 & Self::INSERTION_EXTENDS == 0 {
            GapTrace::Open
        } else {
            GapTrace::Extend
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TracebackState {
    Best,
    Deletion,
    Insertion,
}

pub(crate) fn every_path_score_fits<A, S>(
    substitution: &S,
    gap: AffineGap,
    maximum_operations: usize,
) -> bool
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    every_path_score_fits_for::<A, S>(
        substitution,
        gap,
        maximum_operations,
        i128::from(Score::MIN),
        i128::from(Score::MAX),
    )
}

pub(crate) fn every_path_score_fits_for<A, S>(
    substitution: &S,
    gap: AffineGap,
    maximum_operations: usize,
    minimum_score: i128,
    maximum_score: i128,
) -> bool
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let Some((minimum_substitution, maximum_substitution)) = substitution.score_bounds() else {
        return false;
    };
    if minimum_substitution > maximum_substitution {
        return false;
    }
    let Ok(maximum_operations) = i128::try_from(maximum_operations) else {
        return false;
    };
    let minimum_increment = minimum_substitution
        .min(gap.open_score())
        .min(gap.extension_score());
    let maximum_increment = maximum_substitution
        .max(gap.open_score())
        .max(gap.extension_score());
    let minimum_path_score = i128::from(minimum_increment) * maximum_operations;
    let maximum_path_score = i128::from(maximum_increment) * maximum_operations;
    minimum_path_score >= minimum_score && maximum_path_score <= maximum_score
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceLayout {
    RowMajor,
    ColumnMajor,
}

#[derive(Clone, Copy)]
pub(crate) enum TraceStorage<'a> {
    Packed(&'a [TraceCell]),
    SimdU8(&'a [u8]),
}

impl TraceStorage<'_> {
    pub(crate) fn get(self, index: usize) -> Option<TraceCell> {
        match self {
            Self::Packed(traces) => traces.get(index).copied(),
            Self::SimdU8(traces) => traces.get(index).copied().map(TraceCell::from_bits),
        }
    }
}

impl TraceLayout {
    pub(crate) fn index(
        self,
        row: usize,
        column: usize,
        rows: usize,
        columns: usize,
    ) -> Option<usize> {
        match self {
            Self::RowMajor => row.checked_mul(columns)?.checked_add(column),
            Self::ColumnMajor => column.checked_mul(rows)?.checked_add(row),
        }
    }
}

#[inline]
pub(crate) fn select_best(
    diagonal: Score,
    deletion: Score,
    insertion: Score,
) -> (Score, BestTrace) {
    let (best_gap, gap_trace) = if deletion >= insertion {
        (deletion, BestTrace::Deletion)
    } else {
        (insertion, BestTrace::Insertion)
    };
    if diagonal >= best_gap {
        (diagonal, BestTrace::Diagonal)
    } else {
        (best_gap, gap_trace)
    }
}

#[inline]
pub(crate) fn select_extended_gap(open: Score, extension: Score) -> (Score, GapTrace) {
    if extension >= open {
        (extension, GapTrace::Extend)
    } else {
        (open, GapTrace::Open)
    }
}

pub(crate) fn checked_candidate(
    base: Score,
    increment: Score,
    row: usize,
    column: usize,
) -> Result<Option<Score>, AlignmentError> {
    if let Some(score) = base.checked_add(increment) {
        return Ok(Some(score));
    }
    if increment > 0 {
        Err(AlignmentError::ScoreOverflow { row, column })
    } else {
        Ok(None)
    }
}

#[inline]
pub(crate) fn select_checked_gap(
    open: Option<Score>,
    extension: Option<Score>,
) -> Option<(Score, GapTrace)> {
    match (open, extension) {
        (Some(open), Some(extension)) if extension >= open => Some((extension, GapTrace::Extend)),
        (Some(open), _) => Some((open, GapTrace::Open)),
        (None, Some(extension)) => Some((extension, GapTrace::Extend)),
        (None, None) => None,
    }
}

pub(crate) fn select_checked_best(
    diagonal: Option<Score>,
    deletion: Option<Score>,
    insertion: Option<Score>,
    row: usize,
    column: usize,
) -> Result<(Score, BestTrace), AlignmentError> {
    let mut selected = diagonal.map(|score| (score, BestTrace::Diagonal));
    for (candidate, trace) in [
        (deletion, BestTrace::Deletion),
        (insertion, BestTrace::Insertion),
    ] {
        if let Some(candidate) = candidate {
            match selected {
                Some((best, _)) if candidate <= best => {}
                _ => selected = Some((candidate, trace)),
            }
        }
    }
    selected.ok_or(AlignmentError::ScoreOverflow { row, column })
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn fill_bounded_cached_row<const DELETION_CAN_EXTEND: bool, A>(
    second: &[Symbol<A>],
    substitution_row: &[Score; 256],
    open_score: Score,
    extension_score: Score,
    mut diagonal_base: Score,
    mut left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
) where
    A: Alphabet,
{
    debug_assert_eq!(best_scores.len(), second.len());
    debug_assert_eq!(deletion_scores.len(), second.len());
    debug_assert_eq!(trace_cells.len(), second.len());

    let Some((second_symbol, remaining_second)) = second.split_first() else {
        return;
    };
    let (best_cell, remaining_best) = best_scores
        .split_first_mut()
        .expect("score row matches non-empty sequence");
    let (deletion_cell, remaining_deletions) = deletion_scores
        .split_first_mut()
        .expect("deletion row matches non-empty sequence");
    let (trace_cell, remaining_traces) = trace_cells
        .split_first_mut()
        .expect("trace row matches non-empty sequence");

    let old_best = *best_cell;
    let deletion_open = old_best + open_score;
    let (deletion_score, deletion_trace) = if DELETION_CAN_EXTEND {
        select_extended_gap(deletion_open, *deletion_cell + extension_score)
    } else {
        (deletion_open, GapTrace::Open)
    };
    *deletion_cell = deletion_score;

    let mut insertion_score = left_best + open_score;
    let diagonal = diagonal_base + substitution_row[second_symbol.alphabet_index()];
    let (best, best_trace) = select_best(diagonal, deletion_score, insertion_score);
    *best_cell = best;
    *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, GapTrace::Open);
    diagonal_base = old_best;
    left_best = best;

    for (((second_symbol, best_cell), deletion_cell), trace_cell) in remaining_second
        .iter()
        .copied()
        .zip(remaining_best)
        .zip(remaining_deletions)
        .zip(remaining_traces)
    {
        let old_best = *best_cell;
        let deletion_open = old_best + open_score;
        let (deletion_score, deletion_trace) = if DELETION_CAN_EXTEND {
            select_extended_gap(deletion_open, *deletion_cell + extension_score)
        } else {
            (deletion_open, GapTrace::Open)
        };
        *deletion_cell = deletion_score;

        let insertion_open = left_best + open_score;
        let (new_insertion_score, insertion_trace) =
            select_extended_gap(insertion_open, insertion_score + extension_score);
        insertion_score = new_insertion_score;

        let diagonal = diagonal_base + substitution_row[second_symbol.alphabet_index()];
        let (best, best_trace) = select_best(diagonal, deletion_score, insertion_score);
        *best_cell = best;
        *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, insertion_trace);
        diagonal_base = old_best;
        left_best = best;
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn fill_bounded_row<const DELETION_CAN_EXTEND: bool, A, S>(
    first_symbol: Symbol<A>,
    second: &[Symbol<A>],
    substitution: &S,
    open_score: Score,
    extension_score: Score,
    mut diagonal_base: Score,
    mut left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    debug_assert_eq!(best_scores.len(), second.len());
    debug_assert_eq!(deletion_scores.len(), second.len());
    debug_assert_eq!(trace_cells.len(), second.len());

    let Some((second_symbol, remaining_second)) = second.split_first() else {
        return;
    };
    let (best_cell, remaining_best) = best_scores
        .split_first_mut()
        .expect("score row matches non-empty sequence");
    let (deletion_cell, remaining_deletions) = deletion_scores
        .split_first_mut()
        .expect("deletion row matches non-empty sequence");
    let (trace_cell, remaining_traces) = trace_cells
        .split_first_mut()
        .expect("trace row matches non-empty sequence");

    let old_best = *best_cell;
    let deletion_open = old_best + open_score;
    let (deletion_score, deletion_trace) = if DELETION_CAN_EXTEND {
        select_extended_gap(deletion_open, *deletion_cell + extension_score)
    } else {
        (deletion_open, GapTrace::Open)
    };
    *deletion_cell = deletion_score;

    let mut insertion_score = left_best + open_score;
    let diagonal = diagonal_base + substitution.score(first_symbol, *second_symbol);
    let (best, best_trace) = select_best(diagonal, deletion_score, insertion_score);
    *best_cell = best;
    *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, GapTrace::Open);
    diagonal_base = old_best;
    left_best = best;

    for (((second_symbol, best_cell), deletion_cell), trace_cell) in remaining_second
        .iter()
        .copied()
        .zip(remaining_best)
        .zip(remaining_deletions)
        .zip(remaining_traces)
    {
        let old_best = *best_cell;
        let deletion_open = old_best + open_score;
        let (deletion_score, deletion_trace) = if DELETION_CAN_EXTEND {
            select_extended_gap(deletion_open, *deletion_cell + extension_score)
        } else {
            (deletion_open, GapTrace::Open)
        };
        *deletion_cell = deletion_score;

        let insertion_open = left_best + open_score;
        let (new_insertion_score, insertion_trace) =
            select_extended_gap(insertion_open, insertion_score + extension_score);
        insertion_score = new_insertion_score;

        let diagonal = diagonal_base + substitution.score(first_symbol, second_symbol);
        let (best, best_trace) = select_best(diagonal, deletion_score, insertion_score);
        *best_cell = best;
        *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, insertion_trace);
        diagonal_base = old_best;
        left_best = best;
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn fill_bounded_substitution_row<const DELETION_CAN_EXTEND: bool, A, S>(
    first_symbol: Symbol<A>,
    second: &[Symbol<A>],
    substitution: &S,
    substitution_row: &mut [Score; 256],
    open_score: Score,
    extension_score: Score,
    diagonal_base: Score,
    left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    if width <= substitution_row.len()
        && substitution.fill_score_row(first_symbol, &mut substitution_row[..width])
    {
        fill_bounded_cached_row::<DELETION_CAN_EXTEND, A>(
            second,
            substitution_row,
            open_score,
            extension_score,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            trace_cells,
        );
    } else {
        fill_bounded_row::<DELETION_CAN_EXTEND, A, S>(
            first_symbol,
            second,
            substitution,
            open_score,
            extension_score,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            trace_cells,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_bounded_cached<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    mut substitution_row: [Score; 256],
    gap: AffineGap,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    if columns > 1 {
        best_scores[1] = open_score;
        traceback[1] = TraceCell::default()
            .with_best(BestTrace::Insertion)
            .with_insertion(GapTrace::Open);
        let mut boundary_score = open_score;
        for column in 2..columns {
            boundary_score += extension_score;
            best_scores[column] = boundary_score;
            traceback[column] = TraceCell::default()
                .with_best(BestTrace::Insertion)
                .with_insertion(GapTrace::Extend);
        }
    }

    let Some((_, remaining_first)) = first.split_first() else {
        return;
    };

    let mut diagonal_base = best_scores[0];
    best_scores[0] += open_score;
    deletion_scores[0] = best_scores[0];
    traceback[columns] = TraceCell::default()
        .with_best(BestTrace::Deletion)
        .with_deletion(GapTrace::Open);
    fill_bounded_cached_row::<false, A>(
        second,
        &substitution_row,
        open_score,
        extension_score,
        diagonal_base,
        best_scores[0],
        &mut best_scores[1..],
        &mut deletion_scores[1..],
        &mut traceback[columns + 1..columns * 2],
    );

    for (first_index, first_symbol) in remaining_first.iter().copied().enumerate() {
        let row = first_index + 2;
        let row_start = row * columns;
        diagonal_base = best_scores[0];
        best_scores[0] += extension_score;
        deletion_scores[0] = best_scores[0];
        traceback[row_start] = TraceCell::default()
            .with_best(BestTrace::Deletion)
            .with_deletion(GapTrace::Extend);
        fill_bounded_substitution_row::<true, A, S>(
            first_symbol,
            second,
            substitution,
            &mut substitution_row,
            open_score,
            extension_score,
            diagonal_base,
            best_scores[0],
            &mut best_scores[1..],
            &mut deletion_scores[1..],
            &mut traceback[row_start + 1..row_start + columns],
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_bounded_scalar<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    if columns > 1 {
        best_scores[1] = open_score;
        traceback[1] = TraceCell::default()
            .with_best(BestTrace::Insertion)
            .with_insertion(GapTrace::Open);
        let mut boundary_score = open_score;
        for column in 2..columns {
            boundary_score += extension_score;
            best_scores[column] = boundary_score;
            traceback[column] = TraceCell::default()
                .with_best(BestTrace::Insertion)
                .with_insertion(GapTrace::Extend);
        }
    }

    let Some((first_symbol, remaining_first)) = first.split_first() else {
        return;
    };

    let mut diagonal_base = best_scores[0];
    best_scores[0] += open_score;
    deletion_scores[0] = best_scores[0];
    traceback[columns] = TraceCell::default()
        .with_best(BestTrace::Deletion)
        .with_deletion(GapTrace::Open);
    fill_bounded_row::<false, A, S>(
        *first_symbol,
        second,
        substitution,
        open_score,
        extension_score,
        diagonal_base,
        best_scores[0],
        &mut best_scores[1..],
        &mut deletion_scores[1..],
        &mut traceback[columns + 1..columns * 2],
    );

    for (first_index, first_symbol) in remaining_first.iter().copied().enumerate() {
        let row = first_index + 2;
        let row_start = row * columns;
        diagonal_base = best_scores[0];
        best_scores[0] += extension_score;
        deletion_scores[0] = best_scores[0];
        traceback[row_start] = TraceCell::default()
            .with_best(BestTrace::Deletion)
            .with_deletion(GapTrace::Extend);
        fill_bounded_row::<true, A, S>(
            first_symbol,
            second,
            substitution,
            open_score,
            extension_score,
            diagonal_base,
            best_scores[0],
            &mut best_scores[1..],
            &mut deletion_scores[1..],
            &mut traceback[row_start + 1..row_start + columns],
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_bounded<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    let mut substitution_row = [0; 256];
    let has_cached_row = first.first().is_some_and(|first_symbol| {
        width <= substitution_row.len()
            && substitution.fill_score_row(*first_symbol, &mut substitution_row[..width])
    });
    if has_cached_row {
        fill_bounded_cached(
            first,
            second,
            substitution,
            substitution_row,
            gap,
            best_scores,
            deletion_scores,
            traceback,
            columns,
        );
    } else {
        fill_bounded_scalar(
            first,
            second,
            substitution,
            gap,
            best_scores,
            deletion_scores,
            traceback,
            columns,
        );
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_checked<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    deletion_representable: &mut [u8],
    traceback: &mut [TraceCell],
    columns: usize,
) -> Result<(), AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    for column in 1..columns {
        let increment = if column == 1 {
            open_score
        } else {
            extension_score
        };
        best_scores[column] = best_scores[column - 1]
            .checked_add(increment)
            .ok_or(AlignmentError::ScoreOverflow { row: 0, column })?;
        traceback[column] = TraceCell::default()
            .with_best(BestTrace::Insertion)
            .with_insertion(if column == 1 {
                GapTrace::Open
            } else {
                GapTrace::Extend
            });
    }

    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = best_scores[0];
        let boundary_increment = if row == 1 {
            open_score
        } else {
            extension_score
        };
        best_scores[0] = best_scores[0]
            .checked_add(boundary_increment)
            .ok_or(AlignmentError::ScoreOverflow { row, column: 0 })?;
        deletion_scores[0] = best_scores[0];
        deletion_representable[0] = 1;
        traceback[row_start] = TraceCell::default()
            .with_best(BestTrace::Deletion)
            .with_deletion(if row == 1 {
                GapTrace::Open
            } else {
                GapTrace::Extend
            });

        let mut insertion_score = 0;
        let mut insertion_representable = false;
        for (second_index, second_symbol) in second.iter().copied().enumerate() {
            let column = second_index + 1;
            let old_best = best_scores[column];

            let deletion_open = checked_candidate(old_best, open_score, row, column)?;
            let deletion_extension = if row > 1 && deletion_representable[column] != 0 {
                checked_candidate(deletion_scores[column], extension_score, row, column)?
            } else {
                None
            };
            let selected_deletion = select_checked_gap(deletion_open, deletion_extension);
            let deletion_candidate = selected_deletion.map(|(score, _)| score);
            let deletion_trace =
                selected_deletion.map_or(GapTrace::Open, |(_, predecessor)| predecessor);
            if let Some((score, _)) = selected_deletion {
                deletion_scores[column] = score;
                deletion_representable[column] = 1;
            } else {
                deletion_scores[column] = 0;
                deletion_representable[column] = 0;
            }

            let insertion_open =
                checked_candidate(best_scores[column - 1], open_score, row, column)?;
            let insertion_extension = if column > 1 && insertion_representable {
                checked_candidate(insertion_score, extension_score, row, column)?
            } else {
                None
            };
            let selected_insertion = select_checked_gap(insertion_open, insertion_extension);
            let insertion_candidate = selected_insertion.map(|(score, _)| score);
            let insertion_trace =
                selected_insertion.map_or(GapTrace::Open, |(_, predecessor)| predecessor);
            if let Some((score, _)) = selected_insertion {
                insertion_score = score;
                insertion_representable = true;
            } else {
                insertion_score = 0;
                insertion_representable = false;
            }

            let substitution_score = substitution.score(first_symbol, second_symbol);
            let diagonal = checked_candidate(diagonal_base, substitution_score, row, column)?;
            let (best, best_trace) = select_checked_best(
                diagonal,
                deletion_candidate,
                insertion_candidate,
                row,
                column,
            )?;
            best_scores[column] = best;
            traceback[row_start + column] = TraceCell::default()
                .with_best(best_trace)
                .with_deletion(deletion_trace)
                .with_insertion(insertion_trace);
            diagonal_base = old_best;
        }
    }
    Ok(())
}

/// Exact memory model for the algorithm-owned element buffers.
///
/// Allocation metadata, allocator rounding, and input sequences are not included.
/// `peak_buffer_bytes` covers two rolling score rows, checked-path validity flags,
/// the packed traceback matrix, and the maximum operation-vector capacity retained
/// by the returned result.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GotohMemoryEstimate {
    cells: usize,
    score_buffer_bytes: usize,
    state_validity_buffer_bytes: usize,
    traceback_matrix_bytes: usize,
    maximum_operation_bytes: usize,
    peak_buffer_bytes: usize,
}

impl GotohMemoryEstimate {
    /// Returns the number of dynamic-programming cells, including boundaries.
    #[must_use]
    pub const fn cells(self) -> usize {
        self.cells
    }

    /// Returns bytes occupied by both rolling score rows.
    #[must_use]
    pub const fn score_buffer_bytes(self) -> usize {
        self.score_buffer_bytes
    }

    /// Returns bytes occupied by checked-path state-validity flags.
    #[must_use]
    pub const fn state_validity_buffer_bytes(self) -> usize {
        self.state_validity_buffer_bytes
    }

    /// Returns bytes occupied by packed traceback cells.
    #[must_use]
    pub const fn traceback_matrix_bytes(self) -> usize {
        self.traceback_matrix_bytes
    }

    /// Returns bytes reserved for the maximum possible operation count.
    #[must_use]
    pub const fn maximum_operation_bytes(self) -> usize {
        self.maximum_operation_bytes
    }

    /// Returns the sum of all simultaneously live algorithm-owned buffers.
    #[must_use]
    pub const fn peak_buffer_bytes(self) -> usize {
        self.peak_buffer_bytes
    }
}

pub(crate) fn reserve_workspace_buffer<T>(
    buffer: &mut Vec<T>,
    required_elements: usize,
    name: &'static str,
) -> Result<(), AlignmentError> {
    let additional = required_elements.saturating_sub(buffer.len());
    buffer
        .try_reserve_exact(additional)
        .map_err(|_| AlignmentError::AllocationFailed {
            buffer: name,
            elements: required_elements,
        })
}

fn resize_workspace_buffer<T: Clone>(
    buffer: &mut Vec<T>,
    required_elements: usize,
    value: T,
    name: &'static str,
) -> Result<(), AlignmentError> {
    reserve_workspace_buffer(buffer, required_elements, name)?;
    buffer.resize(required_elements, value);
    Ok(())
}

/// Reusable scratch storage for exact affine pairwise alignment.
///
/// One workspace can be shared sequentially by [`Gotoh`],
/// [`crate::AffineEndGapAligner`], and [`crate::AffineLocalAligner`] calls, including
/// calls over different supported alphabets. Reuse retains scalar score/validity
/// buffers and, when selected, padded SIMD score columns, substitution profiles,
/// symbol indices, and column-major traceback storage. Every returned
/// [`PairwiseAlignment`] remains independently owned.
///
/// The workspace is deliberately passed as `&mut`, which prevents concurrent use of
/// the same scratch buffers. Use one workspace per concurrent alignment.
#[derive(Debug, Default)]
pub struct AffineWorkspace {
    pub(crate) best_scores: Vec<Score>,
    pub(crate) deletion_scores: Vec<Score>,
    pub(crate) deletion_representable: Vec<u8>,
    pub(crate) traceback: Vec<TraceCell>,
    pub(crate) simd_best_scores: Vec<i32>,
    pub(crate) simd_insertion_scores: Vec<i32>,
    pub(crate) simd_profile: Vec<i32>,
    pub(crate) simd_second_indices: Vec<usize>,
    pub(crate) simd_traceback: Vec<u8>,
    kernel_policy: AffineKernelPolicy,
    last_backend: Option<AffineBackend>,
}

impl AffineWorkspace {
    /// Creates an empty workspace that allocates on its first alignment.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            best_scores: Vec::new(),
            deletion_scores: Vec::new(),
            deletion_representable: Vec::new(),
            traceback: Vec::new(),
            simd_best_scores: Vec::new(),
            simd_insertion_scores: Vec::new(),
            simd_profile: Vec::new(),
            simd_second_indices: Vec::new(),
            simd_traceback: Vec::new(),
            kernel_policy: AffineKernelPolicy::Auto,
            last_backend: None,
        }
    }

    /// Selects the affine recurrence backend for subsequent calls using this workspace.
    pub fn set_kernel_policy(&mut self, policy: AffineKernelPolicy) {
        self.kernel_policy = policy;
    }

    /// Returns the configured affine recurrence backend policy.
    #[must_use]
    pub const fn kernel_policy(&self) -> AffineKernelPolicy {
        self.kernel_policy
    }

    /// Returns the backend used by the most recent completed recurrence fill.
    ///
    /// `None` means that no fill has completed since construction or that the most
    /// recent backend selection failed before execution.
    #[must_use]
    pub const fn last_backend(&self) -> Option<AffineBackend> {
        self.last_backend
    }

    pub(crate) fn clear_last_backend(&mut self) {
        self.last_backend = None;
    }

    pub(crate) fn set_last_backend(&mut self, backend: AffineBackend) {
        self.last_backend = Some(backend);
    }

    /// Fallibly reserves scratch capacity for the supplied maximum dimensions.
    ///
    /// This reserves the checked-arithmetic validity row as well as the buffers used
    /// by the common bounded path. Existing allocations and contents are retained.
    /// A later alignment may reuse the capacity for any smaller compatible shape.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if the dimensions cannot be
    /// represented, or [`AlignmentError::AllocationFailed`] if a reservation fails.
    pub fn reserve(&mut self, first_len: usize, second_len: usize) -> Result<(), AlignmentError> {
        let estimate = Gotoh::<()>::memory_requirement(first_len, second_len)?;
        let columns = second_len + 1;
        reserve_workspace_buffer(
            &mut self.best_scores,
            columns,
            "affine workspace best-score row",
        )?;
        reserve_workspace_buffer(
            &mut self.deletion_scores,
            columns,
            "affine workspace deletion-score row",
        )?;
        reserve_workspace_buffer(
            &mut self.deletion_representable,
            columns,
            "affine workspace state-validity row",
        )?;
        reserve_workspace_buffer(
            &mut self.traceback,
            estimate.cells,
            "affine workspace traceback matrix",
        )
    }

    /// Returns the number of columns retained by both rolling score-row buffers.
    #[must_use]
    pub fn score_row_capacity(&self) -> usize {
        self.best_scores
            .capacity()
            .min(self.deletion_scores.capacity())
    }

    /// Returns the number of checked-path validity entries currently retained.
    #[must_use]
    pub fn state_validity_capacity(&self) -> usize {
        self.deletion_representable.capacity()
    }

    /// Returns the number of packed traceback cells currently retained.
    #[must_use]
    pub fn traceback_capacity(&self) -> usize {
        self.traceback.capacity()
    }

    /// Returns the number of signed-64-bit terminal scores retained for SIMD.
    #[must_use]
    pub fn simd_final_score_row_capacity(&self) -> usize {
        self.best_scores.capacity()
    }

    /// Returns the number of signed-32-bit backing entries retained by both SIMD rows.
    #[must_use]
    pub fn simd_score_row_capacity(&self) -> usize {
        self.simd_best_scores
            .capacity()
            .min(self.simd_insertion_scores.capacity())
    }

    /// Returns the number of signed-32-bit backing profile entries retained.
    #[must_use]
    pub fn simd_profile_capacity(&self) -> usize {
        self.simd_profile.capacity()
    }

    /// Returns the number of second-sequence alphabet indices retained for SIMD.
    #[must_use]
    pub fn simd_second_index_capacity(&self) -> usize {
        self.simd_second_indices.capacity()
    }

    /// Returns the number of packed SIMD traceback bytes retained.
    #[must_use]
    pub fn simd_traceback_capacity(&self) -> usize {
        self.simd_traceback.capacity()
    }

    pub(crate) fn prepare(
        &mut self,
        columns: usize,
        cells: usize,
        checked_arithmetic: bool,
    ) -> Result<(), AlignmentError> {
        resize_workspace_buffer(
            &mut self.best_scores,
            columns,
            0,
            "affine workspace best-score row",
        )?;
        resize_workspace_buffer(
            &mut self.deletion_scores,
            columns,
            0,
            "affine workspace deletion-score row",
        )?;
        if checked_arithmetic {
            resize_workspace_buffer(
                &mut self.deletion_representable,
                columns,
                0,
                "affine workspace state-validity row",
            )?;
            self.deletion_representable.fill(0);
        }
        resize_workspace_buffer(
            &mut self.traceback,
            cells,
            TraceCell::default(),
            "affine workspace traceback matrix",
        )?;

        // The fill overwrites every score/deletion dependency before it is read,
        // except the origin and (for checked arithmetic) the first-row validity
        // flags. Every traceback cell that can be visited is also overwritten except
        // a local zero-score origin. Reset only those live initial conditions.
        self.best_scores[0] = 0;
        self.traceback[0] = TraceCell::default();
        Ok(())
    }

    pub(crate) fn prepare_simd(
        &mut self,
        first_len: usize,
        backing_first_len: usize,
        second_len: usize,
        alphabet_width: usize,
        cells: usize,
    ) -> Result<(), AlignmentError> {
        debug_assert!(backing_first_len >= first_len);
        let rows = backing_first_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len: 0,
            })?;
        let columns = second_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let profile_elements = backing_first_len.checked_mul(alphabet_width).ok_or(
            AlignmentError::DimensionOverflow {
                first_len,
                second_len: alphabet_width,
            },
        )?;
        resize_workspace_buffer(
            &mut self.best_scores,
            columns,
            0,
            "affine workspace SIMD final score row",
        )?;
        resize_workspace_buffer(
            &mut self.simd_best_scores,
            rows,
            0,
            "affine workspace SIMD best-score column",
        )?;
        resize_workspace_buffer(
            &mut self.simd_insertion_scores,
            rows,
            0,
            "affine workspace SIMD insertion-score column",
        )?;
        resize_workspace_buffer(
            &mut self.simd_profile,
            profile_elements,
            0,
            "affine workspace SIMD substitution profile",
        )?;
        resize_workspace_buffer(
            &mut self.simd_second_indices,
            second_len,
            0,
            "affine workspace SIMD second-sequence indices",
        )?;
        resize_workspace_buffer(
            &mut self.simd_traceback,
            cells,
            0,
            "affine workspace SIMD traceback matrix",
        )
    }
}

/// Scalar global Gotoh alignment with one symmetric affine gap score.
#[derive(Clone, Debug)]
pub struct Gotoh<S> {
    substitution: S,
    gap: AffineGap,
}

impl<S> Gotoh<S> {
    /// Creates a global Gotoh aligner with explicit substitution and affine scores.
    #[must_use]
    pub const fn new(substitution: S, gap: AffineGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured affine gap score.
    #[must_use]
    pub const fn gap(&self) -> AffineGap {
        self.gap
    }

    /// Computes exact element-buffer requirements for two sequence lengths.
    ///
    /// With [`Gotoh::align_with_workspace`], this estimate describes the active
    /// elements for the call. A reused [`AffineWorkspace`] may retain greater
    /// capacities from an earlier alignment; allocator metadata and retained excess
    /// capacity are not included.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions or any byte count
    /// cannot be represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<GotohMemoryEstimate, AlignmentError> {
        let rows = first_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let columns = second_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let cells = rows
            .checked_mul(columns)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let maximum_operations =
            first_len
                .checked_add(second_len)
                .ok_or(AlignmentError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let score_buffer_bytes = columns
            .checked_mul(size_of::<Score>())
            .and_then(|bytes| bytes.checked_mul(2))
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let state_validity_buffer_bytes =
            columns
                .checked_mul(size_of::<u8>())
                .ok_or(AlignmentError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let traceback_matrix_bytes =
            cells
                .checked_mul(size_of::<TraceCell>())
                .ok_or(AlignmentError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let maximum_operation_bytes = maximum_operations
            .checked_mul(size_of::<AlignmentOperation>())
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let peak_buffer_bytes = score_buffer_bytes
            .checked_add(state_validity_buffer_bytes)
            .and_then(|bytes| bytes.checked_add(traceback_matrix_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_operation_bytes))
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;

        Ok(GotohMemoryEstimate {
            cells,
            score_buffer_bytes,
            state_validity_buffer_bytes,
            traceback_matrix_bytes,
            maximum_operation_bytes,
            peak_buffer_bytes,
        })
    }

    /// Aligns two complete sequences and returns one deterministic optimal traceback.
    ///
    /// Best-state ties prefer diagonal, deletion, then insertion. Within a gap,
    /// extension wins an opening tie, implementing the paper's straight-ahead
    /// traceback guidance.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation failure,
    /// or an internal traceback/result invariant violation.
    #[allow(clippy::too_many_lines)]
    pub fn align<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<PairwiseAlignment, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let mut workspace = AffineWorkspace::new();
        self.align_with_workspace(first, second, &mut workspace)
    }

    /// Aligns two complete sequences while reusing caller-owned scratch storage.
    ///
    /// The score, ranges, operations, and deterministic tie-breaking are identical
    /// to [`Self::align`]. Only scratch-allocation lifetime differs. The returned
    /// alignment does not borrow from `workspace` and remains valid after its next
    /// use or after it is dropped.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation failure,
    /// or an internal traceback/result invariant violation.
    #[allow(clippy::too_many_lines)]
    pub fn align_with_workspace<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AffineWorkspace,
    ) -> Result<PairwiseAlignment, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let first_len = first.len();
        let second_len = second.len();
        let estimate = Self::memory_requirement(first_len, second_len)?;
        let rows = first_len + 1;
        let columns = second_len + 1;
        let maximum_operations = first_len + second_len;
        let bounded =
            every_path_score_fits::<A, S>(&self.substitution, self.gap, maximum_operations);
        let trace_layout = if try_fill_simd(
            first.as_slice(),
            second.as_slice(),
            &self.substitution,
            self.gap,
            SimdBoundaryMode::Global,
            maximum_operations,
            workspace,
        )?
        .is_some()
        {
            TraceLayout::ColumnMajor
        } else {
            workspace.prepare(columns, estimate.cells, !bounded)?;
            let best_scores = &mut workspace.best_scores;
            let deletion_scores = &mut workspace.deletion_scores;
            let traceback = &mut workspace.traceback;
            if bounded {
                fill_bounded(
                    first.as_slice(),
                    second.as_slice(),
                    &self.substitution,
                    self.gap,
                    best_scores,
                    deletion_scores,
                    traceback,
                    columns,
                );
            } else {
                fill_checked(
                    first.as_slice(),
                    second.as_slice(),
                    &self.substitution,
                    self.gap,
                    best_scores,
                    deletion_scores,
                    &mut workspace.deletion_representable,
                    traceback,
                    columns,
                )?;
            }
            workspace.set_last_backend(AffineBackend::Scalar);
            TraceLayout::RowMajor
        };

        let mut operations = Vec::new();
        operations
            .try_reserve_exact(maximum_operations)
            .map_err(|_| AlignmentError::AllocationFailed {
                buffer: "Gotoh traceback operations",
                elements: maximum_operations,
            })?;
        let mut row = first_len;
        let mut column = second_len;
        let mut state = TracebackState::Best;
        let trace_storage = match trace_layout {
            TraceLayout::RowMajor => TraceStorage::Packed(&workspace.traceback),
            TraceLayout::ColumnMajor => TraceStorage::SimdU8(&workspace.simd_traceback),
        };
        while row != 0 || column != 0 {
            let index = trace_layout
                .index(row, column, rows, columns)
                .ok_or(AlignmentError::TracebackInvariant { row, column })?;
            let trace = trace_storage
                .get(index)
                .ok_or(AlignmentError::TracebackInvariant { row, column })?;
            match state {
                TracebackState::Best => match trace.best() {
                    Some(BestTrace::Diagonal) if row != 0 && column != 0 => {
                        let operation =
                            if first.as_slice()[row - 1] == second.as_slice()[column - 1] {
                                AlignmentOperation::Match
                            } else {
                                AlignmentOperation::Substitution
                            };
                        operations.push(operation);
                        row -= 1;
                        column -= 1;
                    }
                    Some(BestTrace::Deletion) if row != 0 => {
                        state = TracebackState::Deletion;
                    }
                    Some(BestTrace::Insertion) if column != 0 => {
                        state = TracebackState::Insertion;
                    }
                    _ => return Err(AlignmentError::TracebackInvariant { row, column }),
                },
                TracebackState::Deletion if row != 0 => {
                    operations.push(AlignmentOperation::Deletion);
                    row -= 1;
                    state = if matches!(trace.deletion(), GapTrace::Extend) {
                        TracebackState::Deletion
                    } else {
                        TracebackState::Best
                    };
                }
                TracebackState::Insertion if column != 0 => {
                    operations.push(AlignmentOperation::Insertion);
                    column -= 1;
                    state = if matches!(trace.insertion(), GapTrace::Extend) {
                        TracebackState::Insertion
                    } else {
                        TracebackState::Best
                    };
                }
                TracebackState::Deletion | TracebackState::Insertion => {
                    return Err(AlignmentError::TracebackInvariant { row, column });
                }
            }
        }
        if !matches!(state, TracebackState::Best) {
            return Err(AlignmentError::TracebackInvariant { row, column });
        }
        operations.reverse();

        PairwiseAlignment::try_new(
            workspace.best_scores[second_len],
            0..first_len,
            0..second_len,
            operations,
        )
        .map_err(AlignmentError::InvalidResult)
    }
}

#[cfg(test)]
mod tests {
    use align_core::{AffineGap, AlignmentOperation, Dna, MatchMismatch};

    use super::{BestTrace, GapTrace, Gotoh, TraceCell, every_path_score_fits};

    #[test]
    fn packed_trace_cell_uses_one_byte() {
        assert_eq!(core::mem::size_of::<TraceCell>(), 1);
    }

    #[test]
    fn direct_trace_cell_constructor_round_trips_every_state() {
        for best in [
            BestTrace::Start,
            BestTrace::Diagonal,
            BestTrace::Deletion,
            BestTrace::Insertion,
        ] {
            for deletion in [GapTrace::Open, GapTrace::Extend] {
                for insertion in [GapTrace::Open, GapTrace::Extend] {
                    let trace = TraceCell::from_traces(best, deletion, insertion);
                    assert_eq!(trace.best(), Some(best));
                    assert_eq!(trace.deletion(), deletion);
                    assert_eq!(trace.insertion(), insertion);
                }
            }
        }
    }

    #[test]
    fn memory_model_matches_element_layout() {
        let estimate = Gotoh::<MatchMismatch>::memory_requirement(3, 4).expect("valid dimensions");
        assert_eq!(estimate.cells(), 20);
        assert_eq!(estimate.score_buffer_bytes(), 80);
        assert_eq!(estimate.state_validity_buffer_bytes(), 5);
        assert_eq!(estimate.traceback_matrix_bytes(), 20);
        assert_eq!(estimate.maximum_operation_bytes(), 7);
        assert_eq!(estimate.peak_buffer_bytes(), 112);
    }

    #[test]
    fn memory_model_rejects_dimension_overflow() {
        assert!(Gotoh::<MatchMismatch>::memory_requirement(usize::MAX, 0).is_err());
    }

    #[test]
    fn configuration_accessors_are_exact() {
        let gap = AffineGap::new(-5, -1).expect("valid gap");
        let aligner = Gotoh::new(MatchMismatch::new(2, -3), gap);
        assert_eq!(aligner.substitution(), &MatchMismatch::new(2, -3));
        assert_eq!(aligner.gap(), gap);
    }

    #[test]
    fn arithmetic_preflight_handles_biological_and_extreme_ranges() {
        let gap = AffineGap::new(-5, -1).expect("valid gap");
        assert!(every_path_score_fits::<Dna, _>(
            &MatchMismatch::new(2, -3),
            gap,
            1_000
        ));
        assert!(!every_path_score_fits::<Dna, _>(
            &MatchMismatch::new(Score::MAX, -1),
            gap,
            2
        ));
    }

    #[test]
    fn linear_specialization_has_expected_canonical_trace() {
        let first = align_core::DnaSequence::from_ascii(b"AC").expect("valid DNA");
        let second = align_core::DnaSequence::from_ascii(b"A").expect("valid DNA");
        let aligner = Gotoh::new(
            MatchMismatch::new(2, -1),
            AffineGap::new(-2, -2).expect("valid gap"),
        );
        let result = aligner.align(&first, &second).expect("alignment succeeds");
        assert_eq!(result.score(), 0);
        assert_eq!(
            result.operations(),
            &[AlignmentOperation::Match, AlignmentOperation::Deletion]
        );
    }

    use align_core::Score;
}
