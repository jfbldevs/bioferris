//! Exact score-only execution for affine pairwise objectives.

use core::mem::size_of;

use align_core::{AffineGap, Alphabet, Score, Sequence, SubstitutionScorer, Symbol};

use crate::affine_simd::{AffineBackend, SimdBoundaryMode, SimdFillOutcome, try_score_simd};
use crate::gotoh::{
    AffineWorkspace, Gotoh, checked_candidate, every_path_score_fits, reserve_workspace_buffer,
    select_checked_best, select_checked_gap,
};
use crate::{AffineEndGapAligner, AffineLocalAligner, AlignmentError, FreeEndGaps};

/// Versioned exact affine score-only execution profile.
pub const AFFINE_SCORE_ONLY_V1: &str = "AFFINE-SCORE-ONLY-v1";

/// Exact active-buffer memory model for scalar affine score-only execution.
///
/// The checked-path peak includes one representability byte per score column. The
/// bounded peak excludes those bytes because a whole-path proof makes them inactive.
/// Allocation metadata, allocator rounding, input sequences, and retained excess
/// workspace capacity are not included.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AffineScoreMemoryEstimate {
    columns: usize,
    score_buffer_bytes: usize,
    state_validity_buffer_bytes: usize,
    bounded_peak_buffer_bytes: usize,
    checked_peak_buffer_bytes: usize,
}

impl AffineScoreMemoryEstimate {
    /// Returns the number of entries in each rolling score row.
    #[must_use]
    pub const fn columns(self) -> usize {
        self.columns
    }

    /// Returns bytes occupied by both rolling signed-64-bit score rows.
    #[must_use]
    pub const fn score_buffer_bytes(self) -> usize {
        self.score_buffer_bytes
    }

    /// Returns bytes occupied by checked-path deletion-state validity flags.
    #[must_use]
    pub const fn state_validity_buffer_bytes(self) -> usize {
        self.state_validity_buffer_bytes
    }

    /// Returns active bytes when scorer bounds prove every path representable.
    #[must_use]
    pub const fn bounded_peak_buffer_bytes(self) -> usize {
        self.bounded_peak_buffer_bytes
    }

    /// Returns active bytes for the fully checked arithmetic path.
    #[must_use]
    pub const fn checked_peak_buffer_bytes(self) -> usize {
        self.checked_peak_buffer_bytes
    }
}

fn score_memory_requirement(
    first_len: usize,
    second_len: usize,
) -> Result<AffineScoreMemoryEstimate, AlignmentError> {
    first_len
        .checked_add(second_len)
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
    let checked_peak_buffer_bytes = score_buffer_bytes
        .checked_add(state_validity_buffer_bytes)
        .ok_or(AlignmentError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    Ok(AffineScoreMemoryEstimate {
        columns,
        score_buffer_bytes,
        state_validity_buffer_bytes,
        bounded_peak_buffer_bytes: score_buffer_bytes,
        checked_peak_buffer_bytes,
    })
}

impl AffineWorkspace {
    /// Fallibly reserves scalar score-only scratch for the supplied dimensions.
    ///
    /// No traceback or operation storage is reserved. Existing allocations and
    /// contents are retained, including inactive capacity from earlier complete
    /// alignments.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension cannot be
    /// represented, or [`AlignmentError::AllocationFailed`] if reservation fails.
    pub fn reserve_score(
        &mut self,
        first_len: usize,
        second_len: usize,
    ) -> Result<(), AlignmentError> {
        let estimate = score_memory_requirement(first_len, second_len)?;
        reserve_workspace_buffer(
            &mut self.best_scores,
            estimate.columns,
            "affine score-only best-score row",
        )?;
        reserve_workspace_buffer(
            &mut self.deletion_scores,
            estimate.columns,
            "affine score-only deletion-score row",
        )?;
        reserve_workspace_buffer(
            &mut self.deletion_representable,
            estimate.columns,
            "affine score-only state-validity row",
        )
    }

    pub(crate) fn prepare_score(
        &mut self,
        columns: usize,
        checked_arithmetic: bool,
    ) -> Result<(), AlignmentError> {
        reserve_workspace_buffer(
            &mut self.best_scores,
            columns,
            "affine score-only best-score row",
        )?;
        self.best_scores.resize(columns, 0);
        reserve_workspace_buffer(
            &mut self.deletion_scores,
            columns,
            "affine score-only deletion-score row",
        )?;
        self.deletion_scores.resize(columns, 0);
        if checked_arithmetic {
            reserve_workspace_buffer(
                &mut self.deletion_representable,
                columns,
                "affine score-only state-validity row",
            )?;
            self.deletion_representable.resize(columns, 0);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ScoreBoundaryMode {
    EndGap(FreeEndGaps),
    Local,
}

impl ScoreBoundaryMode {
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

    const fn simd_mode(self) -> SimdBoundaryMode {
        match self {
            Self::EndGap(free_ends) => SimdBoundaryMode::EndGap {
                first_start_is_free: free_ends.first_start(),
                first_end_is_free: free_ends.first_end(),
                second_start_is_free: free_ends.second_start(),
                second_end_is_free: free_ends.second_end(),
            },
            Self::Local => SimdBoundaryMode::Local,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ScoreFillOutcome {
    local_best: Score,
    last_column_best: Score,
}

impl From<SimdFillOutcome> for ScoreFillOutcome {
    fn from(outcome: SimdFillOutcome) -> Self {
        Self {
            local_best: outcome.local_endpoint.score,
            last_column_best: outcome.last_column_endpoint.score,
        }
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn fill_bounded_score_row<
    const DELETION_CAN_EXTEND: bool,
    const LOCAL: bool,
    const CACHED: bool,
    A,
    S,
>(
    first_symbol: Symbol<A>,
    second: &[Symbol<A>],
    substitution: &S,
    substitution_row: &[Score],
    open_score: Score,
    extension_score: Score,
    mut diagonal_base: Score,
    mut left_best: Score,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    local_best: &mut Score,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    debug_assert_eq!(best_scores.len(), second.len());
    debug_assert_eq!(deletion_scores.len(), second.len());

    let Some((second_symbol, remaining_second)) = second.split_first() else {
        return;
    };
    let (best_cell, remaining_best) = best_scores
        .split_first_mut()
        .expect("score row matches non-empty sequence");
    let (deletion_cell, remaining_deletions) = deletion_scores
        .split_first_mut()
        .expect("deletion row matches non-empty sequence");

    let old_best = *best_cell;
    let deletion_open = old_best + open_score;
    let deletion = if DELETION_CAN_EXTEND {
        deletion_open.max(*deletion_cell + extension_score)
    } else {
        deletion_open
    };
    *deletion_cell = deletion;
    let mut insertion = left_best + open_score;
    let substitution_score = if CACHED {
        substitution_row[second_symbol.alphabet_index()]
    } else {
        substitution.score(first_symbol, *second_symbol)
    };
    let mut best = (diagonal_base + substitution_score)
        .max(deletion)
        .max(insertion);
    if LOCAL {
        best = best.max(0);
        *local_best = (*local_best).max(best);
    }
    *best_cell = best;
    diagonal_base = old_best;
    left_best = best;

    for ((second_symbol, best_cell), deletion_cell) in remaining_second
        .iter()
        .copied()
        .zip(remaining_best)
        .zip(remaining_deletions)
    {
        let old_best = *best_cell;
        let deletion_open = old_best + open_score;
        let deletion = if DELETION_CAN_EXTEND {
            deletion_open.max(*deletion_cell + extension_score)
        } else {
            deletion_open
        };
        *deletion_cell = deletion;
        insertion = (left_best + open_score).max(insertion + extension_score);
        let substitution_score = if CACHED {
            substitution_row[second_symbol.alphabet_index()]
        } else {
            substitution.score(first_symbol, second_symbol)
        };
        let mut best = (diagonal_base + substitution_score)
            .max(deletion)
            .max(insertion);
        if LOCAL {
            best = best.max(0);
            *local_best = (*local_best).max(best);
        }
        *best_cell = best;
        diagonal_base = old_best;
        left_best = best;
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_one_bounded_score_row<const DELETION_CAN_EXTEND: bool, const LOCAL: bool, A, S>(
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
    local_best: &mut Score,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    if width <= substitution_row.len()
        && substitution.fill_score_row(first_symbol, &mut substitution_row[..width])
    {
        fill_bounded_score_row::<DELETION_CAN_EXTEND, LOCAL, true, A, S>(
            first_symbol,
            second,
            substitution,
            substitution_row,
            open_score,
            extension_score,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            local_best,
        );
    } else {
        fill_bounded_score_row::<DELETION_CAN_EXTEND, LOCAL, false, A, S>(
            first_symbol,
            second,
            substitution,
            &[],
            open_score,
            extension_score,
            diagonal_base,
            left_best,
            best_scores,
            deletion_scores,
            local_best,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_bounded_score_for_mode<const LOCAL: bool, A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    first_start_is_free: bool,
    second_start_is_free: bool,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
) -> ScoreFillOutcome
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    best_scores.fill(0);
    deletion_scores.fill(0);
    if !second_start_is_free && !second.is_empty() {
        let mut boundary = open_score;
        for (column, score) in best_scores.iter_mut().enumerate().skip(1) {
            *score = boundary;
            if column < second.len() {
                boundary += extension_score;
            }
        }
    }

    let mut local_best = 0;
    let mut last_column_best = best_scores[second.len()];
    let mut substitution_row = [0; 256];
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let diagonal_base = best_scores[0];
        if first_start_is_free {
            best_scores[0] = 0;
        } else {
            best_scores[0] += if row == 1 {
                open_score
            } else {
                extension_score
            };
        }
        deletion_scores[0] = best_scores[0];
        if row == 1 {
            fill_one_bounded_score_row::<false, LOCAL, A, S>(
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
                &mut local_best,
            );
        } else {
            fill_one_bounded_score_row::<true, LOCAL, A, S>(
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
                &mut local_best,
            );
        }
        last_column_best = last_column_best.max(best_scores[second.len()]);
    }
    ScoreFillOutcome {
        local_best,
        last_column_best,
    }
}

fn fill_bounded_score<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: ScoreBoundaryMode,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
) -> ScoreFillOutcome
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    match mode {
        ScoreBoundaryMode::Local => fill_bounded_score_for_mode::<true, A, S>(
            first,
            second,
            substitution,
            gap,
            true,
            true,
            best_scores,
            deletion_scores,
        ),
        ScoreBoundaryMode::EndGap(free_ends) => fill_bounded_score_for_mode::<false, A, S>(
            first,
            second,
            substitution,
            gap,
            free_ends.first_start(),
            free_ends.second_start(),
            best_scores,
            deletion_scores,
        ),
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_checked_score<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: ScoreBoundaryMode,
    best_scores: &mut [Score],
    deletion_scores: &mut [Score],
    deletion_representable: &mut [u8],
) -> Result<ScoreFillOutcome, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let open_score = gap.open_score();
    let extension_score = gap.extension_score();
    best_scores.fill(0);
    deletion_scores.fill(0);
    deletion_representable.fill(0);
    if !mode.second_start_is_free() {
        for column in 1..best_scores.len() {
            let increment = if column == 1 {
                open_score
            } else {
                extension_score
            };
            best_scores[column] = best_scores[column - 1]
                .checked_add(increment)
                .ok_or(AlignmentError::ScoreOverflow { row: 0, column })?;
        }
    }

    let mut local_best = 0;
    let mut last_column_best = best_scores[second.len()];
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let mut diagonal_base = best_scores[0];
        if mode.first_start_is_free() {
            best_scores[0] = 0;
            deletion_scores[0] = 0;
            deletion_representable[0] = 0;
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
            let best = if mode.is_local()
                && diagonal.is_none()
                && deletion_candidate.is_none()
                && insertion_candidate.is_none()
            {
                0
            } else {
                let selected = select_checked_best(
                    diagonal,
                    deletion_candidate,
                    insertion_candidate,
                    row,
                    column,
                )?
                .0;
                if mode.is_local() {
                    selected.max(0)
                } else {
                    selected
                }
            };
            best_scores[column] = best;
            if mode.is_local() {
                local_best = local_best.max(best);
            }
            diagonal_base = old_best;
        }
        last_column_best = last_column_best.max(best_scores[second.len()]);
    }
    Ok(ScoreFillOutcome {
        local_best,
        last_column_best,
    })
}

fn select_score(
    mode: ScoreBoundaryMode,
    second_len: usize,
    best_scores: &[Score],
    outcome: ScoreFillOutcome,
) -> Score {
    match mode {
        ScoreBoundaryMode::Local => outcome.local_best,
        ScoreBoundaryMode::EndGap(free_ends) => {
            let mut score = best_scores[second_len];
            if free_ends.first_end() {
                score = score.max(outcome.last_column_best);
            }
            if free_ends.second_end() {
                score = score.max(
                    best_scores
                        .iter()
                        .copied()
                        .max()
                        .expect("the terminal score row includes the origin"),
                );
            }
            score
        }
    }
}

fn score_with_mode<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
    mode: ScoreBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<Score, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_len = first.len();
    let second_len = second.len();
    let maximum_operations =
        first_len
            .checked_add(second_len)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
    let bounded = every_path_score_fits::<A, S>(substitution, gap, maximum_operations);
    let outcome = if let Some(outcome) = try_score_simd(
        first.as_slice(),
        second.as_slice(),
        substitution,
        gap,
        mode.simd_mode(),
        maximum_operations,
        workspace,
    )? {
        outcome.into()
    } else {
        let columns = second_len + 1;
        workspace.prepare_score(columns, !bounded)?;
        let outcome = if bounded {
            fill_bounded_score(
                first.as_slice(),
                second.as_slice(),
                substitution,
                gap,
                mode,
                &mut workspace.best_scores,
                &mut workspace.deletion_scores,
            )
        } else {
            fill_checked_score(
                first.as_slice(),
                second.as_slice(),
                substitution,
                gap,
                mode,
                &mut workspace.best_scores,
                &mut workspace.deletion_scores,
                &mut workspace.deletion_representable,
            )?
        };
        workspace.set_last_backend(AffineBackend::Scalar);
        outcome
    };
    Ok(select_score(
        mode,
        second_len,
        &workspace.best_scores,
        outcome,
    ))
}

impl<S> Gotoh<S> {
    /// Computes the exact global affine optimum without constructing a traceback.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let mut workspace = AffineWorkspace::new();
        self.score_with_workspace(first, second, &mut workspace)
    }

    /// Computes the exact global affine optimum with caller-owned scratch storage.
    ///
    /// This method does not construct traceback cells, operations, or aligned
    /// ranges. Its result is identical to [`Self::align_with_workspace`]'s score.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score_with_workspace<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AffineWorkspace,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        score_with_mode(
            first,
            second,
            self.substitution(),
            self.gap(),
            ScoreBoundaryMode::EndGap(FreeEndGaps::global()),
            workspace,
        )
    }

    /// Computes active scalar score-only scratch requirements.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension or byte count
    /// cannot be represented by `usize`.
    pub fn score_memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<AffineScoreMemoryEstimate, AlignmentError> {
        score_memory_requirement(first_len, second_len)
    }
}

impl<S> AffineEndGapAligner<S> {
    /// Computes the exact configured free-end affine optimum without traceback.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let mut workspace = AffineWorkspace::new();
        self.score_with_workspace(first, second, &mut workspace)
    }

    /// Computes the configured free-end affine optimum with reusable scratch.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score_with_workspace<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AffineWorkspace,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        score_with_mode(
            first,
            second,
            self.substitution(),
            self.gap(),
            ScoreBoundaryMode::EndGap(self.free_ends()),
            workspace,
        )
    }

    /// Computes active scalar score-only scratch requirements.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension or byte count
    /// cannot be represented by `usize`.
    pub fn score_memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<AffineScoreMemoryEstimate, AlignmentError> {
        score_memory_requirement(first_len, second_len)
    }
}

impl<S> AffineLocalAligner<S> {
    /// Computes the exact local affine optimum without constructing a traceback.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let mut workspace = AffineWorkspace::new();
        self.score_with_workspace(first, second, &mut workspace)
    }

    /// Computes the exact local affine optimum with reusable scratch storage.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or strict SIMD unavailability.
    pub fn score_with_workspace<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AffineWorkspace,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        score_with_mode(
            first,
            second,
            self.substitution(),
            self.gap(),
            ScoreBoundaryMode::Local,
            workspace,
        )
    }

    /// Computes active scalar score-only scratch requirements.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension or byte count
    /// cannot be represented by `usize`.
    pub fn score_memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<AffineScoreMemoryEstimate, AlignmentError> {
        score_memory_requirement(first_len, second_len)
    }
}

#[cfg(test)]
mod tests {
    use align_core::{AffineGap, DnaSequence, MatchMismatch};

    use super::{AFFINE_SCORE_ONLY_V1, AffineScoreMemoryEstimate};
    use crate::{AffineWorkspace, Gotoh};

    #[test]
    fn profile_and_scalar_memory_model_are_exact() {
        assert_eq!(AFFINE_SCORE_ONLY_V1, "AFFINE-SCORE-ONLY-v1");
        let estimate = Gotoh::<MatchMismatch>::score_memory_requirement(3, 4)
            .expect("small dimensions are representable");
        assert_eq!(
            estimate,
            AffineScoreMemoryEstimate {
                columns: 5,
                score_buffer_bytes: 80,
                state_validity_buffer_bytes: 5,
                bounded_peak_buffer_bytes: 80,
                checked_peak_buffer_bytes: 85,
            }
        );
        assert!(Gotoh::<MatchMismatch>::score_memory_requirement(usize::MAX, 1).is_err());
    }

    #[test]
    fn manual_global_score_uses_no_traceback_storage() {
        let first = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"ACCT").expect("valid DNA");
        let gap = AffineGap::new(-2, -1).expect("valid gap");
        let aligner = Gotoh::new(MatchMismatch::new(1, -1), gap);
        let mut workspace = AffineWorkspace::new();
        let score = aligner
            .score_with_workspace(&first, &second, &mut workspace)
            .expect("score succeeds");
        assert_eq!(score, 2);
        assert_eq!(workspace.traceback_capacity(), 0);
        assert_eq!(workspace.simd_traceback_capacity(), 0);
    }
}
