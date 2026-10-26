use core::mem::size_of;

use align_core::{
    AlignmentOperation, Alphabet, LinearGap, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer, Symbol,
};

use crate::AlignmentError;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum Trace {
    Stop,
    Diagonal,
    Up,
    Left,
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

    fn observe(&mut self, score: Score, row: usize, column: usize) {
        if score > self.score {
            *self = Self { score, row, column };
        }
    }
}

fn every_candidate_score_fits<A, S>(substitution: &S, maximum_operations: usize) -> bool
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

    // Local matrix values are never negative. Adding a non-positive gap or any
    // negative substitution therefore cannot underflow i64. Only accumulated
    // positive substitutions can exceed Score::MAX.
    let maximum_increment = maximum_substitution.max(0);
    i128::from(maximum_increment) * maximum_operations <= i128::from(Score::MAX)
}

#[inline]
fn select_best(diagonal: Score, up: Score, left: Score) -> (Score, Trace) {
    let (best_gap, gap_trace) = if up >= left {
        (up, Trace::Up)
    } else {
        (left, Trace::Left)
    };
    let (best, trace) = if diagonal >= best_gap {
        (diagonal, Trace::Diagonal)
    } else {
        (best_gap, gap_trace)
    };
    if best > 0 {
        (best, trace)
    } else {
        (0, Trace::Stop)
    }
}

#[inline]
fn checked_candidate(
    base: Score,
    increment: Score,
    row: usize,
    column: usize,
) -> Result<Score, AlignmentError> {
    base.checked_add(increment)
        .ok_or(AlignmentError::ScoreOverflow { row, column })
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
) -> Endpoint
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut endpoint = Endpoint::ORIGIN;
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = scores[0];
        let mut left_score = 0;
        scores[0] = 0;
        let trace_row = &mut traceback[row_start..row_start + columns];

        for (second_index, ((second_symbol, score), trace)) in second
            .iter()
            .copied()
            .zip(scores[1..].iter_mut())
            .zip(trace_row[1..].iter_mut())
            .enumerate()
        {
            let column = second_index + 1;
            let up_base = *score;
            let diagonal = diagonal_base + substitution.score(first_symbol, second_symbol);
            let up = up_base + gap_score;
            let left = left_score + gap_score;
            let (best_score, predecessor) = select_best(diagonal, up, left);
            *score = best_score;
            *trace = predecessor;
            endpoint.observe(best_score, row, column);
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
    endpoint
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
) -> Result<Endpoint, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut endpoint = Endpoint::ORIGIN;
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let row_start = row * columns;
        let mut diagonal_base = scores[0];
        let mut left_score = 0;
        scores[0] = 0;

        for (second_index, second_symbol) in second.iter().copied().enumerate() {
            let column = second_index + 1;
            let up_base = scores[column];
            let diagonal = checked_candidate(
                diagonal_base,
                substitution.score(first_symbol, second_symbol),
                row,
                column,
            )?;
            let up = checked_candidate(up_base, gap_score, row, column)?;
            let left = checked_candidate(left_score, gap_score, row, column)?;
            let (best_score, predecessor) = select_best(diagonal, up, left);
            scores[column] = best_score;
            traceback[row_start + column] = predecessor;
            endpoint.observe(best_score, row, column);
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
    Ok(endpoint)
}

// Score-only twins of `fill_bounded` and `fill_checked`: the same recurrence,
// zero floor and overflow checks, keeping only the best cell's score.
fn score_bounded<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
) -> Score
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut best = 0;
    for first_symbol in first.iter().copied() {
        let mut diagonal_base = scores[0];
        let mut left_score = 0;
        scores[0] = 0;
        for (second_symbol, score) in second.iter().copied().zip(scores[1..].iter_mut()) {
            let up_base = *score;
            let diagonal = diagonal_base + substitution.score(first_symbol, second_symbol);
            let up = up_base + gap_score;
            let left = left_score + gap_score;
            let best_score = diagonal.max(up).max(left).max(0);
            *score = best_score;
            best = best.max(best_score);
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
    best
}

fn score_checked<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap_score: Score,
    scores: &mut [Score],
) -> Result<Score, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut best = 0;
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let row = first_index + 1;
        let mut diagonal_base = scores[0];
        let mut left_score = 0;
        scores[0] = 0;
        for (second_index, second_symbol) in second.iter().copied().enumerate() {
            let column = second_index + 1;
            let up_base = scores[column];
            let diagonal = checked_candidate(
                diagonal_base,
                substitution.score(first_symbol, second_symbol),
                row,
                column,
            )?;
            let up = checked_candidate(up_base, gap_score, row, column)?;
            let left = checked_candidate(left_score, gap_score, row, column)?;
            let best_score = diagonal.max(up).max(left).max(0);
            scores[column] = best_score;
            best = best.max(best_score);
            diagonal_base = up_base;
            left_score = best_score;
        }
    }
    Ok(best)
}

/// Exact memory model for Smith–Waterman's algorithm-owned element buffers.
///
/// Allocation metadata, allocator rounding, and input sequences are not included.
/// `peak_buffer_bytes` covers the rolling score buffer, traceback matrix, and the
/// maximum operation-vector capacity retained by the returned result.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SmithWatermanMemoryEstimate {
    cells: usize,
    score_buffer_bytes: usize,
    traceback_matrix_bytes: usize,
    maximum_operation_bytes: usize,
    peak_buffer_bytes: usize,
}

impl SmithWatermanMemoryEstimate {
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

/// Scalar local Smith–Waterman alignment with a linear gap score.
#[derive(Clone, Debug)]
pub struct SmithWaterman<S> {
    substitution: S,
    gap: LinearGap,
}

impl<S> SmithWaterman<S> {
    /// Creates a local aligner with explicit substitution and linear gap scores.
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
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions or byte counts
    /// cannot be represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<SmithWatermanMemoryEstimate, AlignmentError> {
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

        Ok(SmithWatermanMemoryEstimate {
            cells,
            score_buffer_bytes,
            traceback_matrix_bytes,
            maximum_operation_bytes,
            peak_buffer_bytes,
        })
    }

    /// Finds one deterministic maximum-scoring local alignment.
    ///
    /// A zero optimum returns an empty alignment at canonical ranges `0..0`. Positive
    /// cell ties prefer diagonal, deletion, then insertion; equal endpoint maxima keep
    /// their first row-major occurrence.
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
        traceback.resize(estimate.cells, Trace::Stop);

        let maximum_operations = first_len + second_len;
        let gap_score = self.gap.score();
        let endpoint = if every_candidate_score_fits::<A, S>(&self.substitution, maximum_operations)
        {
            fill_bounded(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
                &mut traceback,
                columns,
            )
        } else {
            fill_checked(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
                &mut traceback,
                columns,
            )?
        };

        let mut operations = Vec::new();
        operations
            .try_reserve_exact(maximum_operations)
            .map_err(|_| AlignmentError::AllocationFailed {
                buffer: "traceback operations",
                elements: maximum_operations,
            })?;
        let end_row = endpoint.row;
        let end_column = endpoint.column;
        let mut row = end_row;
        let mut column = end_column;
        loop {
            let index = row
                .checked_mul(columns)
                .and_then(|offset| offset.checked_add(column))
                .ok_or(AlignmentError::TracebackInvariant { row, column })?;
            let predecessor = traceback
                .get(index)
                .copied()
                .ok_or(AlignmentError::TracebackInvariant { row, column })?;
            match predecessor {
                Trace::Stop => break,
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

        PairwiseAlignment::try_new(endpoint.score, row..end_row, column..end_column, operations)
            .map_err(AlignmentError::InvalidResult)
    }

    /// Computes the optimal local score without a traceback.
    ///
    /// The recurrence, zero floor and overflow checks are those of
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
        if every_candidate_score_fits::<A, S>(&self.substitution, maximum_operations) {
            Ok(score_bounded(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
            ))
        } else {
            score_checked(
                first.as_slice(),
                second.as_slice(),
                &self.substitution,
                gap_score,
                &mut scores,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use align_core::{
        AlignmentOperation, Dna, DnaSequence, LinearGap, MatchMismatch, Score, SubstitutionScorer,
        Symbol,
    };

    use super::{SmithWaterman, every_candidate_score_fits};

    #[test]
    fn memory_model_matches_element_layout() {
        let estimate =
            SmithWaterman::<MatchMismatch>::memory_requirement(2, 3).expect("small dimensions fit");
        assert_eq!(estimate.cells(), 12);
        assert_eq!(estimate.score_buffer_bytes(), 4 * 8);
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
        assert!(SmithWaterman::<MatchMismatch>::memory_requirement(usize::MAX, 0).is_err());
    }

    #[test]
    fn configuration_accessors_are_exact() {
        let scoring = MatchMismatch::new(2, -1);
        let gap = LinearGap::new(-2).expect("negative gap is valid");
        let aligner = SmithWaterman::new(scoring, gap);
        assert_eq!(*aligner.substitution(), scoring);
        assert_eq!(aligner.gap(), gap);
    }

    #[test]
    fn arithmetic_preflight_accepts_biological_ranges_and_rejects_positive_extremes() {
        let ordinary = MatchMismatch::new(2, -1);
        assert!(every_candidate_score_fits::<Dna, _>(&ordinary, 2_048));

        let positive_extreme = MatchMismatch::new(i64::MAX, 0);
        assert!(!every_candidate_score_fits::<Dna, _>(&positive_extreme, 2));

        let negative_extreme = MatchMismatch::new(0, i64::MIN);
        assert!(every_candidate_score_fits::<Dna, _>(&negative_extreme, 2));
    }

    #[derive(Clone, Copy)]
    struct UnboundedLiteralScorer;

    impl SubstitutionScorer<Dna> for UnboundedLiteralScorer {
        fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
            if first == second { 2 } else { -1 }
        }
    }

    #[test]
    fn scorer_without_declared_bounds_uses_the_checked_path() {
        assert!(!every_candidate_score_fits::<Dna, _>(
            &UnboundedLiteralScorer,
            8
        ));
        let first = DnaSequence::from_ascii(b"TA").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"CA").expect("valid DNA");
        let aligner = SmithWaterman::new(
            UnboundedLiteralScorer,
            LinearGap::new(-2).expect("valid gap"),
        );
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        assert_eq!(alignment.score(), 2);
        assert_eq!(alignment.first_range(), 1..2);
        assert_eq!(alignment.second_range(), 1..2);
    }
}
