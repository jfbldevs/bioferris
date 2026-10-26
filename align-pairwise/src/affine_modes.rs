use align_core::{
    AffineGap, AlignmentOperation, Alphabet, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer, Symbol,
};

use crate::affine_simd::{AffineBackend, SimdBoundaryMode, SimdEndpoint, try_fill_simd};
use crate::gotoh::{
    BestTrace, GapTrace, TraceCell, TraceLayout, TraceStorage, TracebackState, checked_candidate,
    every_path_score_fits, select_best, select_checked_best, select_checked_gap,
    select_extended_gap,
};
use crate::{AffineWorkspace, AlignmentError, FreeEndGaps, Gotoh, GotohMemoryEstimate};

/// Versioned affine four-end profile implemented by [`AffineEndGapAligner`].
pub const AFFINE_END_GAP_V1: &str = "AFFINE-ENDGAP-v1";

/// Versioned affine local profile implemented by [`AffineLocalAligner`].
pub const SMITH_WATERMAN_GOTOH_LOCAL_AFFINE_V1: &str = "SMITH-WATERMAN-GOTOH-local-affine-v1";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BoundaryMode {
    EndGap(FreeEndGaps),
    Local,
}

impl BoundaryMode {
    const fn first_start_is_free(self) -> bool {
        match self {
            Self::EndGap(free_ends) => free_ends.first_start(),
            Self::Local => true,
        }
    }

    const fn second_start_is_free(self) -> bool {
        match self {
            Self::EndGap(free_ends) => free_ends.second_start(),
            Self::Local => true,
        }
    }

    const fn is_local(self) -> bool {
        matches!(self, Self::Local)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Endpoint {
    score: Score,
    row: usize,
    column: usize,
}

impl Endpoint {
    const ORIGIN: Self = Self {
        score: 0,
        row: 0,
        column: 0,
    };

    const fn new(score: Score, row: usize, column: usize) -> Self {
        Self { score, row, column }
    }

    fn observe_local(&mut self, score: Score, row: usize, column: usize) {
        if score > self.score {
            *self = Self::new(score, row, column);
        }
    }

    fn observe_end_gap(&mut self, score: Score, row: usize, column: usize) {
        let candidate_consumed = row + column;
        let selected_consumed = self.row + self.column;
        if score > self.score
            || (score == self.score && candidate_consumed > selected_consumed)
            || (score == self.score && candidate_consumed == selected_consumed && row > self.row)
        {
            *self = Self::new(score, row, column);
        }
    }
}

impl From<SimdEndpoint> for Endpoint {
    fn from(endpoint: SimdEndpoint) -> Self {
        Self::new(endpoint.score, endpoint.row, endpoint.column)
    }
}

#[derive(Clone, Copy)]
struct FillOutcome {
    local_endpoint: Endpoint,
    last_column_endpoint: Endpoint,
}

fn local_best(diagonal: Score, deletion: Score, insertion: Score) -> (Score, BestTrace) {
    let (score, trace) = select_best(diagonal, deletion, insertion);
    if score > 0 {
        (score, trace)
    } else {
        (0, BestTrace::Start)
    }
}

fn checked_local_best(
    diagonal: Option<Score>,
    deletion: Option<Score>,
    insertion: Option<Score>,
    row: usize,
    column: usize,
) -> Result<(Score, BestTrace), AlignmentError> {
    if diagonal.is_none() && deletion.is_none() && insertion.is_none() {
        return Ok((0, BestTrace::Start));
    }
    let (score, trace) = select_checked_best(diagonal, deletion, insertion, row, column)?;
    if score > 0 {
        Ok((score, trace))
    } else {
        Ok((0, BestTrace::Start))
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[inline]
fn fill_bounded_cached_row<const DELETION_CAN_EXTEND: bool, const LOCAL: bool, A>(
    second: &[Symbol<A>],
    substitution_row: &[Score; 256],
    open_score: Score,
    extension_score: Score,
    row: usize,
    mut diagonal_base: Score,
    mut left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
    local_endpoint: &mut Endpoint,
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
    let (best, best_trace) = if LOCAL {
        local_best(diagonal, deletion_score, insertion_score)
    } else {
        select_best(diagonal, deletion_score, insertion_score)
    };
    *best_cell = best;
    *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, GapTrace::Open);
    if LOCAL {
        local_endpoint.observe_local(best, row, 1);
    }
    diagonal_base = old_best;
    left_best = best;

    for (offset, (((second_symbol, best_cell), deletion_cell), trace_cell)) in remaining_second
        .iter()
        .copied()
        .zip(remaining_best)
        .zip(remaining_deletions)
        .zip(remaining_traces)
        .enumerate()
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
        let (best, best_trace) = if LOCAL {
            local_best(diagonal, deletion_score, insertion_score)
        } else {
            select_best(diagonal, deletion_score, insertion_score)
        };
        *best_cell = best;
        *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, insertion_trace);
        if LOCAL {
            local_endpoint.observe_local(best, row, offset + 2);
        }
        diagonal_base = old_best;
        left_best = best;
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[inline]
fn fill_bounded_row<const DELETION_CAN_EXTEND: bool, const LOCAL: bool, A, S>(
    first_symbol: Symbol<A>,
    second: &[Symbol<A>],
    substitution: &S,
    open_score: Score,
    extension_score: Score,
    row: usize,
    mut diagonal_base: Score,
    mut left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
    local_endpoint: &mut Endpoint,
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
    let (best, best_trace) = if LOCAL {
        local_best(diagonal, deletion_score, insertion_score)
    } else {
        select_best(diagonal, deletion_score, insertion_score)
    };
    *best_cell = best;
    *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, GapTrace::Open);
    if LOCAL {
        local_endpoint.observe_local(best, row, 1);
    }
    diagonal_base = old_best;
    left_best = best;

    for (offset, (((second_symbol, best_cell), deletion_cell), trace_cell)) in remaining_second
        .iter()
        .copied()
        .zip(remaining_best)
        .zip(remaining_deletions)
        .zip(remaining_traces)
        .enumerate()
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
        let (best, best_trace) = if LOCAL {
            local_best(diagonal, deletion_score, insertion_score)
        } else {
            select_best(diagonal, deletion_score, insertion_score)
        };
        *best_cell = best;
        *trace_cell = TraceCell::from_traces(best_trace, deletion_trace, insertion_trace);
        if LOCAL {
            local_endpoint.observe_local(best, row, offset + 2);
        }
        diagonal_base = old_best;
        left_best = best;
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn fill_bounded_substitution_row<const DELETION_CAN_EXTEND: bool, const LOCAL: bool, A, S>(
    first_symbol: Symbol<A>,
    second: &[Symbol<A>],
    substitution: &S,
    substitution_row: &mut [Score; 256],
    open_score: Score,
    extension_score: Score,
    row: usize,
    diagonal_base: Score,
    left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    trace_cells: &mut [TraceCell],
    local_endpoint: &mut Endpoint,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    if width <= substitution_row.len()
        && substitution.fill_score_row(first_symbol, &mut substitution_row[..width])
    {
        fill_bounded_cached_row::<DELETION_CAN_EXTEND, LOCAL, A>(
            second,
            substitution_row,
            open_score,
            extension_score,
            row,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            trace_cells,
            local_endpoint,
        );
    } else {
        fill_bounded_row::<DELETION_CAN_EXTEND, LOCAL, A, S>(
            first_symbol,
            second,
            substitution,
            open_score,
            extension_score,
            row,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            trace_cells,
            local_endpoint,
        );
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_bounded_cached_for_mode<const LOCAL: bool, A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    mut substitution_row: [Score; 256],
    gap: AffineGap,
    first_start_is_free: bool,
    second_start_is_free: bool,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) -> FillOutcome
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    if !second_start_is_free && columns > 1 {
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

    let second_len = second.len();
    let mut local_endpoint = Endpoint::ORIGIN;
    let mut last_column_endpoint = Endpoint::new(best_scores[second_len], 0, second_len);

    let Some((_, remaining_first)) = first.split_first() else {
        return FillOutcome {
            local_endpoint,
            last_column_endpoint,
        };
    };

    let mut diagonal_base = best_scores[0];
    if first_start_is_free {
        best_scores[0] = 0;
        deletion_scores[0] = 0;
        traceback[columns] = TraceCell::default();
    } else {
        best_scores[0] += open_score;
        deletion_scores[0] = best_scores[0];
        traceback[columns] = TraceCell::default()
            .with_best(BestTrace::Deletion)
            .with_deletion(GapTrace::Open);
    }
    fill_bounded_cached_row::<false, LOCAL, A>(
        second,
        &substitution_row,
        open_score,
        extension_score,
        1,
        diagonal_base,
        best_scores[0],
        &mut best_scores[1..],
        &mut deletion_scores[1..],
        &mut traceback[columns + 1..columns * 2],
        &mut local_endpoint,
    );
    last_column_endpoint.observe_end_gap(best_scores[second_len], 1, second_len);

    for (first_index, first_symbol) in remaining_first.iter().copied().enumerate() {
        let row = first_index + 2;
        let row_start = row * columns;
        diagonal_base = best_scores[0];
        if first_start_is_free {
            best_scores[0] = 0;
            deletion_scores[0] = 0;
            traceback[row_start] = TraceCell::default();
        } else {
            best_scores[0] += extension_score;
            deletion_scores[0] = best_scores[0];
            traceback[row_start] = TraceCell::default()
                .with_best(BestTrace::Deletion)
                .with_deletion(GapTrace::Extend);
        }
        fill_bounded_substitution_row::<true, LOCAL, A, S>(
            first_symbol,
            second,
            substitution,
            &mut substitution_row,
            open_score,
            extension_score,
            row,
            diagonal_base,
            best_scores[0],
            &mut best_scores[1..],
            &mut deletion_scores[1..],
            &mut traceback[row_start + 1..row_start + columns],
            &mut local_endpoint,
        );
        last_column_endpoint.observe_end_gap(best_scores[second_len], row, second_len);
    }

    FillOutcome {
        local_endpoint,
        last_column_endpoint,
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_bounded_scalar_for_mode<const LOCAL: bool, A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    first_start_is_free: bool,
    second_start_is_free: bool,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) -> FillOutcome
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    if !second_start_is_free && columns > 1 {
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

    let second_len = second.len();
    let mut local_endpoint = Endpoint::ORIGIN;
    let mut last_column_endpoint = Endpoint::new(best_scores[second_len], 0, second_len);

    let Some((first_symbol, remaining_first)) = first.split_first() else {
        return FillOutcome {
            local_endpoint,
            last_column_endpoint,
        };
    };

    let mut diagonal_base = best_scores[0];
    if first_start_is_free {
        best_scores[0] = 0;
        deletion_scores[0] = 0;
        traceback[columns] = TraceCell::default();
    } else {
        best_scores[0] += open_score;
        deletion_scores[0] = best_scores[0];
        traceback[columns] = TraceCell::default()
            .with_best(BestTrace::Deletion)
            .with_deletion(GapTrace::Open);
    }
    fill_bounded_row::<false, LOCAL, A, S>(
        *first_symbol,
        second,
        substitution,
        open_score,
        extension_score,
        1,
        diagonal_base,
        best_scores[0],
        &mut best_scores[1..],
        &mut deletion_scores[1..],
        &mut traceback[columns + 1..columns * 2],
        &mut local_endpoint,
    );
    last_column_endpoint.observe_end_gap(best_scores[second_len], 1, second_len);

    for (first_index, first_symbol) in remaining_first.iter().copied().enumerate() {
        let row = first_index + 2;
        let row_start = row * columns;
        diagonal_base = best_scores[0];
        if first_start_is_free {
            best_scores[0] = 0;
            deletion_scores[0] = 0;
            traceback[row_start] = TraceCell::default();
        } else {
            best_scores[0] += extension_score;
            deletion_scores[0] = best_scores[0];
            traceback[row_start] = TraceCell::default()
                .with_best(BestTrace::Deletion)
                .with_deletion(GapTrace::Extend);
        }
        fill_bounded_row::<true, LOCAL, A, S>(
            first_symbol,
            second,
            substitution,
            open_score,
            extension_score,
            row,
            diagonal_base,
            best_scores[0],
            &mut best_scores[1..],
            &mut deletion_scores[1..],
            &mut traceback[row_start + 1..row_start + columns],
            &mut local_endpoint,
        );
        last_column_endpoint.observe_end_gap(best_scores[second_len], row, second_len);
    }

    FillOutcome {
        local_endpoint,
        last_column_endpoint,
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_bounded<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: BoundaryMode,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    traceback: &mut [TraceCell],
    columns: usize,
) -> FillOutcome
where
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
        match mode {
            BoundaryMode::Local => fill_bounded_cached_for_mode::<true, A, S>(
                first,
                second,
                substitution,
                substitution_row,
                gap,
                true,
                true,
                best_scores,
                deletion_scores,
                traceback,
                columns,
            ),
            BoundaryMode::EndGap(free_ends) => fill_bounded_cached_for_mode::<false, A, S>(
                first,
                second,
                substitution,
                substitution_row,
                gap,
                free_ends.first_start(),
                free_ends.second_start(),
                best_scores,
                deletion_scores,
                traceback,
                columns,
            ),
        }
    } else {
        match mode {
            BoundaryMode::Local => fill_bounded_scalar_for_mode::<true, A, S>(
                first,
                second,
                substitution,
                gap,
                true,
                true,
                best_scores,
                deletion_scores,
                traceback,
                columns,
            ),
            BoundaryMode::EndGap(free_ends) => fill_bounded_scalar_for_mode::<false, A, S>(
                first,
                second,
                substitution,
                gap,
                free_ends.first_start(),
                free_ends.second_start(),
                best_scores,
                deletion_scores,
                traceback,
                columns,
            ),
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_checked<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: BoundaryMode,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    deletion_representable: &mut [u8],
    traceback: &mut [TraceCell],
    columns: usize,
) -> Result<FillOutcome, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    if !mode.second_start_is_free() {
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
    }

    let second_len = second.len();
    let mut local_endpoint = Endpoint::ORIGIN;
    let mut last_column_endpoint = Endpoint::new(best_scores[second_len], 0, second_len);
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = best_scores[0];
        if mode.first_start_is_free() {
            best_scores[0] = 0;
            deletion_scores[0] = 0;
            deletion_representable[0] = 0;
            traceback[row_start] = TraceCell::default();
        } else {
            let increment = if row == 1 {
                open_score
            } else {
                extension_score
            };
            best_scores[0] = best_scores[0]
                .checked_add(increment)
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
        }

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

            let diagonal = checked_candidate(
                diagonal_base,
                substitution.score(first_symbol, second_symbol),
                row,
                column,
            )?;
            let (best, best_trace) = if mode.is_local() {
                checked_local_best(
                    diagonal,
                    deletion_candidate,
                    insertion_candidate,
                    row,
                    column,
                )?
            } else {
                select_checked_best(
                    diagonal,
                    deletion_candidate,
                    insertion_candidate,
                    row,
                    column,
                )?
            };
            best_scores[column] = best;
            traceback[row_start + column] = TraceCell::default()
                .with_best(best_trace)
                .with_deletion(deletion_trace)
                .with_insertion(insertion_trace);
            if mode.is_local() {
                local_endpoint.observe_local(best, row, column);
            }
            diagonal_base = old_best;
        }
        last_column_endpoint.observe_end_gap(best_scores[second_len], row, second_len);
    }

    Ok(FillOutcome {
        local_endpoint,
        last_column_endpoint,
    })
}

fn select_endpoint(
    mode: BoundaryMode,
    first_len: usize,
    second_len: usize,
    best_scores: &[Score],
    outcome: FillOutcome,
) -> Endpoint {
    match mode {
        BoundaryMode::Local => outcome.local_endpoint,
        BoundaryMode::EndGap(free_ends) => {
            let mut endpoint = Endpoint::new(best_scores[second_len], first_len, second_len);
            if free_ends.first_end() {
                endpoint.observe_end_gap(
                    outcome.last_column_endpoint.score,
                    outcome.last_column_endpoint.row,
                    outcome.last_column_endpoint.column,
                );
            }
            if free_ends.second_end() {
                for (column, score) in best_scores.iter().copied().enumerate() {
                    endpoint.observe_end_gap(score, first_len, column);
                }
            }
            endpoint
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn traceback<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    mode: BoundaryMode,
    endpoint: Endpoint,
    traces: TraceStorage<'_>,
    layout: TraceLayout,
    rows: usize,
    columns: usize,
    maximum_operations: usize,
) -> Result<PairwiseAlignment, AlignmentError> {
    let mut operations = Vec::new();
    operations
        .try_reserve_exact(maximum_operations)
        .map_err(|_| AlignmentError::AllocationFailed {
            buffer: "affine-mode traceback operations",
            elements: maximum_operations,
        })?;
    let end_row = endpoint.row;
    let end_column = endpoint.column;
    let mut row = end_row;
    let mut column = end_column;
    let mut state = TracebackState::Best;

    loop {
        if matches!(state, TracebackState::Best) {
            match mode {
                BoundaryMode::EndGap(free_ends)
                    if (column == 0 && free_ends.first_start())
                        || (row == 0 && free_ends.second_start()) =>
                {
                    break;
                }
                BoundaryMode::EndGap(_) if row == 0 && column == 0 => break,
                BoundaryMode::Local | BoundaryMode::EndGap(_) => {}
            }
        }

        let index = layout
            .index(row, column, rows, columns)
            .ok_or(AlignmentError::TracebackInvariant { row, column })?;
        let trace = traces
            .get(index)
            .ok_or(AlignmentError::TracebackInvariant { row, column })?;
        match state {
            TracebackState::Best => match trace.best() {
                Some(BestTrace::Start) if mode.is_local() => break,
                Some(BestTrace::Diagonal) if row != 0 && column != 0 => {
                    operations.push(
                        if first.as_slice()[row - 1] == second.as_slice()[column - 1] {
                            AlignmentOperation::Match
                        } else {
                            AlignmentOperation::Substitution
                        },
                    );
                    row -= 1;
                    column -= 1;
                }
                Some(BestTrace::Deletion) if row != 0 => {
                    state = TracebackState::Deletion;
                }
                Some(BestTrace::Insertion) if column != 0 => {
                    state = TracebackState::Insertion;
                }
                _ => {
                    return Err(AlignmentError::TracebackInvariant { row, column });
                }
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

    PairwiseAlignment::try_new(endpoint.score, row..end_row, column..end_column, operations)
        .map_err(AlignmentError::InvalidResult)
}

#[allow(clippy::too_many_lines)]
fn align_with_mode<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
    mode: BoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<PairwiseAlignment, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_len = first.len();
    let second_len = second.len();
    let estimate = Gotoh::<S>::memory_requirement(first_len, second_len)?;
    let rows = first_len + 1;
    let columns = second_len + 1;
    let maximum_operations =
        first_len
            .checked_add(second_len)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;

    let bounded = every_path_score_fits::<A, S>(substitution, gap, maximum_operations);
    let simd_mode = match mode {
        BoundaryMode::Local => SimdBoundaryMode::Local,
        BoundaryMode::EndGap(free_ends) => SimdBoundaryMode::EndGap {
            first_start_is_free: free_ends.first_start(),
            first_end_is_free: free_ends.first_end(),
            second_start_is_free: free_ends.second_start(),
            second_end_is_free: free_ends.second_end(),
        },
    };
    let (outcome, trace_layout) = if let Some(simd_outcome) = try_fill_simd(
        first.as_slice(),
        second.as_slice(),
        substitution,
        gap,
        simd_mode,
        maximum_operations,
        workspace,
    )? {
        (
            FillOutcome {
                local_endpoint: simd_outcome.local_endpoint.into(),
                last_column_endpoint: simd_outcome.last_column_endpoint.into(),
            },
            TraceLayout::ColumnMajor,
        )
    } else {
        workspace.prepare(columns, estimate.cells(), !bounded)?;
        if matches!(mode, BoundaryMode::Local)
            || matches!(mode, BoundaryMode::EndGap(free_ends) if free_ends.second_start())
        {
            workspace.best_scores.fill(0);
        }
        if matches!(mode, BoundaryMode::Local) {
            workspace.traceback[..columns].fill(TraceCell::default());
        }
        let best_scores = &mut workspace.best_scores;
        let deletion_scores = &mut workspace.deletion_scores;
        let traces = &mut workspace.traceback;
        let outcome = if bounded {
            fill_bounded(
                first.as_slice(),
                second.as_slice(),
                substitution,
                gap,
                mode,
                best_scores,
                deletion_scores,
                traces,
                columns,
            )
        } else {
            fill_checked(
                first.as_slice(),
                second.as_slice(),
                substitution,
                gap,
                mode,
                best_scores,
                deletion_scores,
                &mut workspace.deletion_representable,
                traces,
                columns,
            )?
        };
        workspace.set_last_backend(AffineBackend::Scalar);
        (outcome, TraceLayout::RowMajor)
    };
    let endpoint = select_endpoint(mode, first_len, second_len, &workspace.best_scores, outcome);
    traceback(
        first,
        second,
        mode,
        endpoint,
        match trace_layout {
            TraceLayout::RowMajor => TraceStorage::Packed(&workspace.traceback),
            TraceLayout::ColumnMajor => TraceStorage::SimdU8(&workspace.simd_traceback),
        },
        trace_layout,
        rows,
        columns,
        maximum_operations,
    )
}

/// Scalar Gotoh alignment with explicit permissions for all four sequence ends.
#[derive(Clone, Debug)]
pub struct AffineEndGapAligner<S> {
    substitution: S,
    gap: AffineGap,
    free_ends: FreeEndGaps,
}

impl<S> AffineEndGapAligner<S> {
    /// Creates an affine aligner with explicit scoring and four-end permissions.
    #[must_use]
    pub const fn new(substitution: S, gap: AffineGap, free_ends: FreeEndGaps) -> Self {
        Self {
            substitution,
            gap,
            free_ends,
        }
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

    /// Returns the configured free-end permissions.
    #[must_use]
    pub const fn free_ends(&self) -> FreeEndGaps {
        self.free_ends
    }

    /// Returns the exact versioned algorithm profile.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        AFFINE_END_GAP_V1
    }

    /// Computes exact element-buffer requirements for two sequence lengths.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions or byte counts
    /// cannot be represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<GotohMemoryEstimate, AlignmentError> {
        Gotoh::<S>::memory_requirement(first_len, second_len)
    }

    /// Aligns the configured sequence ranges with affine internal and scored-end gaps.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or an invalid traceback/result invariant.
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

    /// Aligns with explicit end-gap permissions while reusing scratch storage.
    ///
    /// Results and deterministic tie-breaking are identical to [`Self::align`]. The
    /// returned alignment owns its operations and does not borrow from `workspace`.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or an invalid traceback/result invariant.
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
        align_with_mode(
            first,
            second,
            &self.substitution,
            self.gap,
            BoundaryMode::EndGap(self.free_ends),
            workspace,
        )
    }
}

/// Scalar Smith-Waterman local alignment with Gotoh affine gap states.
#[derive(Clone, Debug)]
pub struct AffineLocalAligner<S> {
    substitution: S,
    gap: AffineGap,
}

impl<S> AffineLocalAligner<S> {
    /// Creates an affine local aligner with explicit scoring.
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

    /// Returns the exact versioned algorithm profile.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        SMITH_WATERMAN_GOTOH_LOCAL_AFFINE_V1
    }

    /// Computes exact element-buffer requirements for two sequence lengths.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions or byte counts
    /// cannot be represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<GotohMemoryEstimate, AlignmentError> {
        Gotoh::<S>::memory_requirement(first_len, second_len)
    }

    /// Finds one deterministic maximum-scoring affine local alignment.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or an invalid traceback/result invariant.
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

    /// Finds a local alignment while reusing caller-owned scratch storage.
    ///
    /// Results and deterministic tie-breaking are identical to [`Self::align`]. The
    /// returned alignment owns its operations and does not borrow from `workspace`.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or an invalid traceback/result invariant.
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
        align_with_mode(
            first,
            second,
            &self.substitution,
            self.gap,
            BoundaryMode::Local,
            workspace,
        )
    }
}

#[cfg(test)]
mod tests {
    use align_core::{AffineGap, DnaSequence, MatchMismatch};

    use super::{
        AFFINE_END_GAP_V1, AffineEndGapAligner, AffineLocalAligner,
        SMITH_WATERMAN_GOTOH_LOCAL_AFFINE_V1,
    };
    use crate::{FreeEndGaps, Gotoh};

    fn gap() -> AffineGap {
        AffineGap::new(-5, -1).expect("valid affine gap")
    }

    #[test]
    fn profiles_configuration_and_memory_are_public() {
        let end_gap =
            AffineEndGapAligner::new(MatchMismatch::new(2, -1), gap(), FreeEndGaps::semiglobal());
        assert_eq!(end_gap.profile_identifier(), AFFINE_END_GAP_V1);
        assert_eq!(end_gap.substitution(), &MatchMismatch::new(2, -1));
        assert_eq!(end_gap.gap(), gap());
        assert_eq!(end_gap.free_ends(), FreeEndGaps::semiglobal());
        assert!(AffineEndGapAligner::<MatchMismatch>::memory_requirement(usize::MAX, 0).is_err());

        let local = AffineLocalAligner::new(MatchMismatch::new(2, -1), gap());
        assert_eq!(
            local.profile_identifier(),
            SMITH_WATERMAN_GOTOH_LOCAL_AFFINE_V1
        );
        assert_eq!(local.substitution(), &MatchMismatch::new(2, -1));
        assert_eq!(local.gap(), gap());
        assert!(AffineLocalAligner::<MatchMismatch>::memory_requirement(2, 3).is_ok());
    }

    #[test]
    fn global_configuration_reduces_to_gotoh() {
        let first = DnaSequence::from_ascii(b"ACGTT").expect("first");
        let second = DnaSequence::from_ascii(b"ACT").expect("second");
        let scoring = MatchMismatch::new(2, -1);
        let expected = Gotoh::new(scoring, gap())
            .align(&first, &second)
            .expect("global Gotoh");
        let observed = AffineEndGapAligner::new(scoring, gap(), FreeEndGaps::global())
            .align(&first, &second)
            .expect("affine global mode");
        assert_eq!(observed, expected);
    }

    #[test]
    fn local_zero_and_positive_results_have_canonical_ranges() {
        let aligner = AffineLocalAligner::new(MatchMismatch::new(2, -3), gap());
        let zero = aligner
            .align(
                &DnaSequence::from_ascii(b"AAAA").expect("first"),
                &DnaSequence::from_ascii(b"CCCC").expect("second"),
            )
            .expect("zero local");
        assert_eq!(zero.score(), 0);
        assert_eq!(zero.first_range(), 0..0);
        assert_eq!(zero.second_range(), 0..0);

        let positive = aligner
            .align(
                &DnaSequence::from_ascii(b"TTACGTAA").expect("first"),
                &DnaSequence::from_ascii(b"GGACGTCC").expect("second"),
            )
            .expect("positive local");
        assert_eq!(positive.score(), 8);
        assert_eq!(positive.first_range(), 2..6);
        assert_eq!(positive.second_range(), 2..6);
    }
}
