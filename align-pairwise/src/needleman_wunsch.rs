use core::mem::size_of;

use align_core::{
    AlignmentOperation, Alphabet, LinearGap, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer, Symbol,
};

use crate::AlignmentError;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum Trace {
    Start,
    Diagonal,
    Up,
    Left,
}

#[cold]
#[inline(never)]
fn select_after_candidate_overflow(
    candidates: [(Score, bool, Score, Trace); 3],
    row: usize,
    column: usize,
) -> Result<(Score, Trace), AlignmentError> {
    if candidates
        .iter()
        .any(|(_, overflowed, increment, _)| *overflowed && *increment > 0)
    {
        // A positively overflowing candidate is greater than every representable
        // score, so the cell's mathematical maximum cannot be stored.
        return Err(AlignmentError::ScoreOverflow { row, column });
    }

    let mut selected = None;
    for (candidate, overflowed, _, predecessor) in candidates {
        if !overflowed {
            match selected {
                None => selected = Some((candidate, predecessor)),
                Some((best, _)) if candidate > best => {
                    selected = Some((candidate, predecessor));
                }
                Some(_) => {}
            }
        }
    }

    // Missing candidates can only be negative overflows at this point. If all
    // three are missing, even the mathematical maximum is below Score::MIN.
    selected.ok_or(AlignmentError::ScoreOverflow { row, column })
}

fn every_path_score_fits<A, S>(
    substitution: &S,
    gap_score: Score,
    maximum_operations: usize,
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
    let minimum_increment = minimum_substitution.min(gap_score);
    let maximum_increment = maximum_substitution.max(gap_score);
    let minimum_path_score = i128::from(minimum_increment) * maximum_operations;
    let maximum_path_score = i128::from(maximum_increment) * maximum_operations;
    minimum_path_score >= i128::from(Score::MIN) && maximum_path_score <= i128::from(Score::MAX)
}

#[inline]
fn select_best(diagonal: Score, up: Score, left: Score) -> (Score, Trace) {
    let (best_gap, gap_trace) = if up >= left {
        (up, Trace::Up)
    } else {
        (left, Trace::Left)
    };
    if diagonal >= best_gap {
        (diagonal, Trace::Diagonal)
    } else {
        (best_gap, gap_trace)
    }
}

#[inline]
fn fill_bounded<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
    traceback: &mut [Trace],
    columns: usize,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut boundary_score = 0;
    for (score, trace) in scores[1..].iter_mut().zip(traceback[1..columns].iter_mut()) {
        boundary_score += gap_score;
        *score = boundary_score;
        *trace = Trace::Left;
    }

    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = scores[0];
        let mut left_score = diagonal_base + gap_score;
        scores[0] = left_score;
        let trace_row = &mut traceback[row_start..row_start + columns];
        trace_row[0] = Trace::Up;
        for ((second_symbol, score), trace) in second
            .iter()
            .copied()
            .zip(scores[1..].iter_mut())
            .zip(trace_row[1..].iter_mut())
        {
            let up_base = *score;
            let diagonal = diagonal_base + substitution.score(first_symbol, second_symbol);
            let up = up_base + gap_score;
            let left = left_score + gap_score;
            let (best_score, predecessor) = select_best(diagonal, up, left);
            *score = best_score;
            *trace = predecessor;
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn fill_checked<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
    traceback: &mut [Trace],
    columns: usize,
) -> Result<(), AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    for column in 1..columns {
        scores[column] = scores[column - 1]
            .checked_add(gap_score)
            .ok_or(AlignmentError::ScoreOverflow { row: 0, column })?;
        traceback[column] = Trace::Left;
    }

    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = scores[0];
        scores[0] = diagonal_base
            .checked_add(gap_score)
            .ok_or(AlignmentError::ScoreOverflow { row, column: 0 })?;
        traceback[row_start] = Trace::Up;
        for (second_index, second_symbol) in second.iter().copied().enumerate() {
            let column = second_index + 1;
            let substitution_score = substitution.score(first_symbol, second_symbol);
            let up_base = scores[column];
            let (diagonal, diagonal_overflowed) = diagonal_base.overflowing_add(substitution_score);
            let (up, up_overflowed) = up_base.overflowing_add(gap_score);
            let (left, left_overflowed) = scores[column - 1].overflowing_add(gap_score);

            let any_overflowed = diagonal_overflowed | up_overflowed | left_overflowed;
            let (best_score, predecessor) = if any_overflowed {
                select_after_candidate_overflow(
                    [
                        (
                            diagonal,
                            diagonal_overflowed,
                            substitution_score,
                            Trace::Diagonal,
                        ),
                        (up, up_overflowed, gap_score, Trace::Up),
                        (left, left_overflowed, gap_score, Trace::Left),
                    ],
                    row,
                    column,
                )?
            } else {
                select_best(diagonal, up, left)
            };
            scores[column] = best_score;
            traceback[row_start + column] = predecessor;
            diagonal_base = up_base;
        }
    }
    Ok(())
}

// Score-only twins of `fill_bounded` and `fill_checked`: the same recurrence,
// boundary values and overflow checks, without traceback writes.
fn score_bounded<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut boundary_score = 0;
    for score in &mut scores[1..] {
        boundary_score += gap_score;
        *score = boundary_score;
    }

    for first_symbol in first.iter().copied() {
        let mut diagonal_base = scores[0];
        let mut left_score = diagonal_base + gap_score;
        scores[0] = left_score;
        for (second_symbol, score) in second.iter().copied().zip(scores[1..].iter_mut()) {
            let up_base = *score;
            let diagonal = diagonal_base + substitution.score(first_symbol, second_symbol);
            let up = up_base + gap_score;
            let left = left_score + gap_score;
            let best_score = diagonal.max(up).max(left);
            *score = best_score;
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
}

fn score_checked<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
) -> Result<(), AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    for column in 1..scores.len() {
        scores[column] = scores[column - 1]
            .checked_add(gap_score)
            .ok_or(AlignmentError::ScoreOverflow { row: 0, column })?;
    }

    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let mut diagonal_base = scores[0];
        scores[0] = diagonal_base
            .checked_add(gap_score)
            .ok_or(AlignmentError::ScoreOverflow { row, column: 0 })?;
        for (second_index, second_symbol) in second.iter().copied().enumerate() {
            let column = second_index + 1;
            let substitution_score = substitution.score(first_symbol, second_symbol);
            let up_base = scores[column];
            let (diagonal, diagonal_overflowed) = diagonal_base.overflowing_add(substitution_score);
            let (up, up_overflowed) = up_base.overflowing_add(gap_score);
            let (left, left_overflowed) = scores[column - 1].overflowing_add(gap_score);

            let best_score = if diagonal_overflowed | up_overflowed | left_overflowed {
                select_after_candidate_overflow(
                    [
                        (
                            diagonal,
                            diagonal_overflowed,
                            substitution_score,
                            Trace::Diagonal,
                        ),
                        (up, up_overflowed, gap_score, Trace::Up),
                        (left, left_overflowed, gap_score, Trace::Left),
                    ],
                    row,
                    column,
                )?
                .0
            } else {
                diagonal.max(up).max(left)
            };
            scores[column] = best_score;
            diagonal_base = up_base;
        }
    }
    Ok(())
}

/// Exact memory model for the algorithm-owned element buffers.
///
/// Allocation metadata, allocator rounding, and input sequences are not included.
/// `peak_buffer_bytes` covers the rolling score buffer, traceback matrix, and the
/// maximum operation-vector capacity retained by the returned result.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MemoryEstimate {
    cells: usize,
    score_buffer_bytes: usize,
    traceback_matrix_bytes: usize,
    maximum_operation_bytes: usize,
    peak_buffer_bytes: usize,
}

impl MemoryEstimate {
    /// Returns the number of dynamic-programming cells, including boundaries.
    #[must_use]
    pub const fn cells(self) -> usize {
        self.cells
    }

    /// Returns bytes occupied by the rolling score buffer.
    #[must_use]
    pub const fn score_buffer_bytes(self) -> usize {
        self.score_buffer_bytes
    }

    /// Returns bytes occupied by score elements.
    ///
    /// This compatibility accessor is equivalent to [`Self::score_buffer_bytes`].
    #[must_use]
    pub const fn score_matrix_bytes(self) -> usize {
        self.score_buffer_bytes
    }

    /// Returns bytes occupied by traceback elements.
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

/// Scalar global Needleman–Wunsch alignment with a linear gap score.
#[derive(Clone, Debug)]
pub struct NeedlemanWunsch<S> {
    substitution: S,
    gap: LinearGap,
}

impl<S> NeedlemanWunsch<S> {
    /// Creates an aligner with an explicit substitution scorer and linear gap score.
    #[must_use]
    pub const fn new(substitution: S, gap: LinearGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured linear gap score.
    #[must_use]
    pub const fn gap(&self) -> LinearGap {
        self.gap
    }

    /// Computes exact element-buffer requirements for two sequence lengths.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if the matrix dimensions or
    /// any byte count cannot be represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<MemoryEstimate, AlignmentError> {
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
        let score_buffer_bytes =
            columns
                .checked_mul(size_of::<Score>())
                .ok_or(AlignmentError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let traceback_matrix_bytes =
            cells
                .checked_mul(size_of::<Trace>())
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
            .checked_add(traceback_matrix_bytes)
            .and_then(|bytes| bytes.checked_add(maximum_operation_bytes))
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;

        Ok(MemoryEstimate {
            cells,
            score_buffer_bytes,
            traceback_matrix_bytes,
            maximum_operation_bytes,
            peak_buffer_bytes,
        })
    }

    /// Aligns two complete sequences and returns one deterministic optimal traceback.
    ///
    /// Ties are resolved in this order: diagonal, deletion (`Up`), insertion
    /// (`Left`). No alternative optimal tracebacks are retained.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, allocation
    /// failure, or an internal traceback/result invariant violation.
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
        let first_len = first.len();
        let second_len = second.len();
        let estimate = Self::memory_requirement(first_len, second_len)?;
        let columns = second_len + 1;

        let mut scores: Vec<Score> = Vec::new();
        scores
            .try_reserve_exact(columns)
            .map_err(|_| AlignmentError::AllocationFailed {
                buffer: "score buffer",
                elements: columns,
            })?;
        scores.resize(columns, 0);

        let mut traceback = Vec::new();
        traceback.try_reserve_exact(estimate.cells).map_err(|_| {
            AlignmentError::AllocationFailed {
                buffer: "traceback matrix",
                elements: estimate.cells,
            }
        })?;
        traceback.resize(estimate.cells, Trace::Start);

        let maximum_operations = first_len + second_len;
        let gap_score = self.gap.score();
        if every_path_score_fits::<A, S>(&self.substitution, gap_score, maximum_operations) {
            fill_bounded(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
                &mut traceback,
                columns,
            );
        } else {
            fill_checked(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
                &mut traceback,
                columns,
            )?;
        }

        let mut operations = Vec::new();
        operations
            .try_reserve_exact(maximum_operations)
            .map_err(|_| AlignmentError::AllocationFailed {
                buffer: "traceback operations",
                elements: maximum_operations,
            })?;
        let mut row = first_len;
        let mut column = second_len;
        while row != 0 || column != 0 {
            let index = row * columns + column;
            let predecessor = traceback
                .get(index)
                .copied()
                .ok_or(AlignmentError::TracebackInvariant { row, column })?;
            match predecessor {
                Trace::Diagonal if row != 0 && column != 0 => {
                    let operation = if first.as_slice()[row - 1] == second.as_slice()[column - 1] {
                        AlignmentOperation::Match
                    } else {
                        AlignmentOperation::Substitution
                    };
                    operations.push(operation);
                    row -= 1;
                    column -= 1;
                }
                Trace::Up if row != 0 => {
                    operations.push(AlignmentOperation::Deletion);
                    row -= 1;
                }
                Trace::Left if column != 0 => {
                    operations.push(AlignmentOperation::Insertion);
                    column -= 1;
                }
                _ => return Err(AlignmentError::TracebackInvariant { row, column }),
            }
        }
        operations.reverse();

        PairwiseAlignment::try_new(scores[second_len], 0..first_len, 0..second_len, operations)
            .map_err(AlignmentError::InvalidResult)
    }

    /// Computes the optimal global score without a traceback.
    ///
    /// The recurrence, boundary values and overflow checks are those of
    /// [`Self::align`], evaluated over one rolling row of `second.len() + 1`
    /// scores, so the result equals `align(first, second)?.score()` while the
    /// algorithm-owned memory stays `8 * (second.len() + 1)` bytes.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] on dimension or score overflow, or when the
    /// score row cannot be allocated. A score overflow names the same cell as in
    /// [`Self::align`].
    pub fn score<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<Score, AlignmentError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let first_len = first.len();
        let second_len = second.len();
        let dimension_overflow = || AlignmentError::DimensionOverflow {
            first_len,
            second_len,
        };
        let columns = second_len.checked_add(1).ok_or_else(dimension_overflow)?;
        let maximum_operations = first_len
            .checked_add(second_len)
            .ok_or_else(dimension_overflow)?;
        let mut scores: Vec<Score> = Vec::new();
        scores
            .try_reserve_exact(columns)
            .map_err(|_| AlignmentError::AllocationFailed {
                buffer: "score buffer",
                elements: columns,
            })?;
        scores.resize(columns, 0);

        let gap_score = self.gap.score();
        if every_path_score_fits::<A, S>(&self.substitution, gap_score, maximum_operations) {
            score_bounded(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
            );
        } else {
            score_checked(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
            )?;
        }
        Ok(scores[second_len])
    }
}

#[cfg(test)]
mod tests {
    use align_core::{AlignmentOperation, Dna, LinearGap, MatchMismatch};

    use super::{NeedlemanWunsch, every_path_score_fits};

    #[test]
    fn memory_model_matches_element_layout() {
        let estimate = NeedlemanWunsch::<MatchMismatch>::memory_requirement(2, 3)
            .expect("small dimensions fit");
        assert_eq!(estimate.cells(), 12);
        assert_eq!(estimate.score_buffer_bytes(), 4 * 8);
        assert_eq!(estimate.score_matrix_bytes(), estimate.score_buffer_bytes());
        assert_eq!(estimate.traceback_matrix_bytes(), 12);
        assert_eq!(
            estimate.maximum_operation_bytes(),
            5 * core::mem::size_of::<AlignmentOperation>()
        );
        assert_eq!(
            estimate.peak_buffer_bytes(),
            estimate.score_buffer_bytes()
                + estimate.traceback_matrix_bytes()
                + estimate.maximum_operation_bytes()
        );
    }

    #[test]
    fn memory_model_rejects_dimension_overflow() {
        assert!(NeedlemanWunsch::<MatchMismatch>::memory_requirement(usize::MAX, 0).is_err());
    }

    #[test]
    fn configuration_accessors_are_exact() {
        let scoring = MatchMismatch::new(2, -1);
        let gap = LinearGap::new(-2).expect("negative gap is valid");
        let aligner = NeedlemanWunsch::new(scoring, gap);
        assert_eq!(*aligner.substitution(), scoring);
        assert_eq!(aligner.gap(), gap);
    }

    #[test]
    fn arithmetic_preflight_accepts_biological_ranges_and_rejects_extremes() {
        let ordinary = MatchMismatch::new(2, -1);
        assert!(every_path_score_fits::<Dna, _>(&ordinary, -2, 2_048));

        let positive_extreme = MatchMismatch::new(i64::MAX, 0);
        assert!(!every_path_score_fits::<Dna, _>(&positive_extreme, 0, 2));

        let negative_extreme = MatchMismatch::new(0, i64::MIN);
        assert!(!every_path_score_fits::<Dna, _>(&negative_extreme, -1, 2));
    }
}
