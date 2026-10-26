use core::mem::size_of;
use core::ops::Range;

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer, Symbol,
};

use crate::XDropError;

/// Reason a forward X-drop traversal stopped.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum XDropTermination {
    /// Every real matrix anti-diagonal was considered.
    MatrixExhausted,
    /// Neither of the two preceding retained frontiers could generate a successor.
    FrontiersExhausted,
}

/// Direction in which a forced seed is extended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum XDropDirection {
    /// Extend toward lower sequence indices and finish at the seed.
    Upstream,
    /// Start at the seed and extend toward higher sequence indices.
    Downstream,
    /// Extend independently on both sides of the seed.
    Both,
}

impl XDropDirection {
    const fn includes_upstream(self) -> bool {
        matches!(self, Self::Upstream | Self::Both)
    }

    const fn includes_downstream(self) -> bool {
        matches!(self, Self::Downstream | Self::Both)
    }
}

/// One forced pair of zero-based symbol indices for two-way X-drop extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct XDropSeed {
    first_index: usize,
    second_index: usize,
}

impl XDropSeed {
    /// Constructs a seed from one symbol index in each sequence.
    #[must_use]
    pub const fn new(first_index: usize, second_index: usize) -> Self {
        Self {
            first_index,
            second_index,
        }
    }

    /// Returns the first-sequence symbol index.
    #[must_use]
    pub const fn first_index(self) -> usize {
        self.first_index
    }

    /// Returns the second-sequence symbol index.
    #[must_use]
    pub const fn second_index(self) -> usize {
        self.second_index
    }
}

/// Score, endpoint, and traversal metadata for one forward X-drop extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct XDropScore {
    score: Score,
    first_end: usize,
    second_end: usize,
    evaluated_cells: usize,
    retained_cells: usize,
    antidiagonals: usize,
    maximum_frontier_width: usize,
    termination: XDropTermination,
}

impl XDropScore {
    /// Returns the maximum score among retained cells.
    #[must_use]
    pub const fn score(self) -> Score {
        self.score
    }

    /// Returns the exclusive endpoint in the first extension view.
    #[must_use]
    pub const fn first_end(self) -> usize {
        self.first_end
    }

    /// Returns the exclusive endpoint in the second extension view.
    #[must_use]
    pub const fn second_end(self) -> usize {
        self.second_end
    }

    /// Returns the number of real candidate cells evaluated, including the origin.
    #[must_use]
    pub const fn evaluated_cells(self) -> usize {
        self.evaluated_cells
    }

    /// Returns the number of cells retained after X-drop, including the origin.
    #[must_use]
    pub const fn retained_cells(self) -> usize {
        self.retained_cells
    }

    /// Returns the number of evaluated anti-diagonals, including the origin.
    #[must_use]
    pub const fn antidiagonals(self) -> usize {
        self.antidiagonals
    }

    /// Returns the largest contiguous candidate interval evaluated.
    #[must_use]
    pub const fn maximum_frontier_width(self) -> usize {
        self.maximum_frontier_width
    }

    /// Returns why traversal stopped.
    #[must_use]
    pub const fn termination(self) -> XDropTermination {
        self.termination
    }
}

/// One forward X-drop traceback plus traversal metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct XDropAlignment {
    alignment: PairwiseAlignment,
    summary: XDropScore,
}

impl XDropAlignment {
    /// Borrows the operation-based alignment.
    #[must_use]
    pub const fn alignment(&self) -> &PairwiseAlignment {
        &self.alignment
    }

    /// Consumes this result and returns the operation-based alignment.
    #[must_use]
    pub fn into_alignment(self) -> PairwiseAlignment {
        self.alignment
    }

    /// Returns the maximum retained score.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.summary.score
    }

    /// Returns the exclusive first-sequence endpoint.
    #[must_use]
    pub const fn first_end(&self) -> usize {
        self.summary.first_end
    }

    /// Returns the exclusive second-sequence endpoint.
    #[must_use]
    pub const fn second_end(&self) -> usize {
        self.summary.second_end
    }

    /// Returns the score-only traversal summary.
    #[must_use]
    pub const fn summary(&self) -> XDropScore {
        self.summary
    }

    /// Returns the number of evaluated real candidate cells.
    #[must_use]
    pub const fn evaluated_cells(&self) -> usize {
        self.summary.evaluated_cells
    }

    /// Returns the number of retained cells.
    #[must_use]
    pub const fn retained_cells(&self) -> usize {
        self.summary.retained_cells
    }

    /// Returns why traversal stopped.
    #[must_use]
    pub const fn termination(&self) -> XDropTermination {
        self.summary.termination
    }
}

/// Score, ranges, and component traversals for one forced seed extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct XDropSeedScore {
    score: Score,
    first_start: usize,
    first_end: usize,
    second_start: usize,
    second_end: usize,
    evaluated_cells: usize,
    retained_cells: usize,
    antidiagonals: usize,
    maximum_frontier_width: usize,
    upstream: Option<XDropScore>,
    downstream: Option<XDropScore>,
}

impl XDropSeedScore {
    /// Returns the combined score, including the forced seed pair exactly once.
    #[must_use]
    pub const fn score(self) -> Score {
        self.score
    }

    /// Returns the consumed range in the first sequence.
    #[must_use]
    pub fn first_range(self) -> Range<usize> {
        self.first_start..self.first_end
    }

    /// Returns the consumed range in the second sequence.
    #[must_use]
    pub fn second_range(self) -> Range<usize> {
        self.second_start..self.second_end
    }

    /// Returns aggregate evaluated cells, counting the seed origin once.
    #[must_use]
    pub const fn evaluated_cells(self) -> usize {
        self.evaluated_cells
    }

    /// Returns aggregate retained cells, counting the seed origin once.
    #[must_use]
    pub const fn retained_cells(self) -> usize {
        self.retained_cells
    }

    /// Returns aggregate directional anti-diagonals, counting the seed origin once.
    #[must_use]
    pub const fn antidiagonals(self) -> usize {
        self.antidiagonals
    }

    /// Returns the largest candidate interval on either side of the seed.
    #[must_use]
    pub const fn maximum_frontier_width(self) -> usize {
        self.maximum_frontier_width
    }

    /// Returns the upstream traversal, when requested.
    #[must_use]
    pub const fn upstream(self) -> Option<XDropScore> {
        self.upstream
    }

    /// Returns the downstream traversal, when requested.
    #[must_use]
    pub const fn downstream(self) -> Option<XDropScore> {
        self.downstream
    }
}

/// One forced-seed X-drop traceback plus directional metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct XDropSeedAlignment {
    alignment: PairwiseAlignment,
    summary: XDropSeedScore,
}

impl XDropSeedAlignment {
    /// Borrows the operation-based seeded alignment.
    #[must_use]
    pub const fn alignment(&self) -> &PairwiseAlignment {
        &self.alignment
    }

    /// Consumes this result and returns the operation-based alignment.
    #[must_use]
    pub fn into_alignment(self) -> PairwiseAlignment {
        self.alignment
    }

    /// Returns the combined seeded score.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.summary.score
    }

    /// Returns the consumed first-sequence range.
    #[must_use]
    pub fn first_range(&self) -> Range<usize> {
        self.summary.first_range()
    }

    /// Returns the consumed second-sequence range.
    #[must_use]
    pub fn second_range(&self) -> Range<usize> {
        self.summary.second_range()
    }

    /// Returns the score-only seeded summary.
    #[must_use]
    pub const fn summary(&self) -> XDropSeedScore {
        self.summary
    }
}

/// Worst-case modeled element-buffer requirements for one forward X-drop call.
///
/// Actual traceback history contains only evaluated candidate cells. Allocator
/// metadata, `Vec` headers, inputs, result headers, and stack scalars are excluded.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct XDropMemoryEstimate {
    maximum_frontier_width: usize,
    working_cells: usize,
    maximum_trace_cells: usize,
    metadata_entries: usize,
    maximum_operations: usize,
    score_peak_bytes: usize,
    alignment_peak_bytes: usize,
}

impl XDropMemoryEstimate {
    /// Returns the maximum number of cells on one real anti-diagonal.
    #[must_use]
    pub const fn maximum_frontier_width(self) -> usize {
        self.maximum_frontier_width
    }

    /// Returns affine cells reserved across three rolling frontiers.
    #[must_use]
    pub const fn working_cells(self) -> usize {
        self.working_cells
    }

    /// Returns the full-matrix upper bound on packed trace cells.
    #[must_use]
    pub const fn maximum_trace_cells(self) -> usize {
        self.maximum_trace_cells
    }

    /// Returns the maximum number of anti-diagonal metadata entries.
    #[must_use]
    pub const fn metadata_entries(self) -> usize {
        self.metadata_entries
    }

    /// Returns the maximum output-operation capacity.
    #[must_use]
    pub const fn maximum_operations(self) -> usize {
        self.maximum_operations
    }

    /// Returns peak modeled score-only bytes.
    #[must_use]
    pub const fn score_peak_bytes(self) -> usize {
        self.score_peak_bytes
    }

    /// Returns the worst-case peak modeled alignment bytes.
    #[must_use]
    pub const fn alignment_peak_bytes(self) -> usize {
        self.alignment_peak_bytes
    }
}

/// Affine dynamic-programming X-drop extension.
///
/// The forward primitive fixes the origin and selects the maximum-scoring retained
/// endpoint. X-drop pruning is heuristic and can miss an unrestricted optimum.
#[derive(Clone, Debug)]
pub struct XDrop<S> {
    substitution: S,
    gap: AffineGap,
    threshold: Score,
}

impl<S> XDrop<S> {
    /// Constructs an immutable X-drop configuration.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError::NegativeThreshold`] when `threshold` is negative.
    pub fn new(substitution: S, gap: AffineGap, threshold: Score) -> Result<Self, XDropError> {
        if threshold < 0 {
            return Err(XDropError::NegativeThreshold { threshold });
        }
        Ok(Self {
            substitution,
            gap,
            threshold,
        })
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

    /// Returns the strict X-drop threshold.
    #[must_use]
    pub const fn threshold(&self) -> Score {
        self.threshold
    }

    /// Models worst-case element-buffer requirements for a forward extension.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError::DimensionOverflow`] if a count or byte total cannot be
    /// represented by `usize`.
    pub fn memory_requirement(
        first_len: usize,
        second_len: usize,
    ) -> Result<XDropMemoryEstimate, XDropError> {
        let dimension_error = || XDropError::DimensionOverflow {
            first_len,
            second_len,
        };
        let maximum_operations = first_len
            .checked_add(second_len)
            .ok_or_else(dimension_error)?;
        let metadata_entries = maximum_operations
            .checked_add(1)
            .ok_or_else(dimension_error)?;
        let maximum_frontier_width = first_len
            .min(second_len)
            .checked_add(1)
            .ok_or_else(dimension_error)?;
        let working_cells = maximum_frontier_width
            .checked_mul(3)
            .ok_or_else(dimension_error)?;
        let maximum_trace_cells = first_len
            .checked_add(1)
            .and_then(|rows| {
                second_len
                    .checked_add(1)
                    .and_then(|columns| rows.checked_mul(columns))
            })
            .ok_or_else(dimension_error)?;

        let score_peak_bytes = working_cells
            .checked_mul(size_of::<AffineCell>())
            .ok_or_else(dimension_error)?;
        let trace_bytes = maximum_trace_cells
            .checked_mul(size_of::<TraceCell>())
            .ok_or_else(dimension_error)?;
        let metadata_bytes = metadata_entries
            .checked_mul(size_of::<TraceFrontier>())
            .ok_or_else(dimension_error)?;
        let operation_bytes = maximum_operations
            .checked_mul(size_of::<AlignmentOperation>())
            .ok_or_else(dimension_error)?;
        let persistent_history = trace_bytes
            .checked_add(metadata_bytes)
            .ok_or_else(dimension_error)?;
        let alignment_peak_bytes = persistent_history
            .checked_add(score_peak_bytes.max(operation_bytes))
            .ok_or_else(dimension_error)?;

        Ok(XDropMemoryEstimate {
            maximum_frontier_width,
            working_cells,
            maximum_trace_cells,
            metadata_entries,
            maximum_operations,
            score_peak_bytes,
            alignment_peak_bytes,
        })
    }

    /// Scores a forward extension from the fixed origin.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError`] for unrepresentable dimensions or scores, or if a
    /// rolling frontier cannot be allocated.
    pub fn score<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<XDropScore, XDropError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        fill_forward::<A, S, false>(
            SymbolView::forward(first.as_slice()),
            SymbolView::forward(second.as_slice()),
            &self.substitution,
            self.gap,
            self.threshold,
            None,
        )
    }

    /// Computes one deterministic forward extension traceback.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError`] for dimension or score overflow, allocation failure,
    /// invalid stored traceback, or invalid result structure.
    pub fn align<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<XDropAlignment, XDropError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let run = self.align_view(
            SymbolView::forward(first.as_slice()),
            SymbolView::forward(second.as_slice()),
        )?;
        let alignment = PairwiseAlignment::try_new(
            run.summary.score,
            0..run.summary.first_end,
            0..run.summary.second_end,
            run.operations,
        )
        .map_err(XDropError::InvalidResult)?;
        Ok(XDropAlignment {
            alignment,
            summary: run.summary,
        })
    }

    /// Scores an upstream, downstream, or two-way extension around a forced seed.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError::SeedOutOfBounds`] unless the seed identifies one symbol
    /// in each sequence, plus the execution errors documented by [`Self::score`].
    pub fn score_from_seed<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<XDropSeedScore, XDropError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let (first_seed, second_seed) = validate_seed(first, second, seed)?;
        let upstream = if direction.includes_upstream() {
            Some(fill_forward::<A, S, false>(
                SymbolView::reversed(&first.as_slice()[..seed.first_index]),
                SymbolView::reversed(&second.as_slice()[..seed.second_index]),
                &self.substitution,
                self.gap,
                self.threshold,
                None,
            )?)
        } else {
            None
        };
        let downstream = if direction.includes_downstream() {
            Some(fill_forward::<A, S, false>(
                SymbolView::forward(&first.as_slice()[seed.first_index + 1..]),
                SymbolView::forward(&second.as_slice()[seed.second_index + 1..]),
                &self.substitution,
                self.gap,
                self.threshold,
                None,
            )?)
        } else {
            None
        };
        build_seed_summary(
            first.len(),
            second.len(),
            seed,
            self.substitution.score(first_seed, second_seed),
            upstream,
            downstream,
        )
    }

    /// Computes one deterministic extension traceback around a forced seed.
    ///
    /// # Errors
    ///
    /// Returns [`XDropError::SeedOutOfBounds`] unless the seed identifies one symbol
    /// in each sequence, plus allocation, arithmetic, traceback, or result errors.
    pub fn align_from_seed<A>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<XDropSeedAlignment, XDropError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let (first_seed, second_seed) = validate_seed(first, second, seed)?;
        let mut upstream_run = if direction.includes_upstream() {
            Some(self.align_view(
                SymbolView::reversed(&first.as_slice()[..seed.first_index]),
                SymbolView::reversed(&second.as_slice()[..seed.second_index]),
            )?)
        } else {
            None
        };
        let downstream_run = if direction.includes_downstream() {
            Some(self.align_view(
                SymbolView::forward(&first.as_slice()[seed.first_index + 1..]),
                SymbolView::forward(&second.as_slice()[seed.second_index + 1..]),
            )?)
        } else {
            None
        };

        let summary = build_seed_summary(
            first.len(),
            second.len(),
            seed,
            self.substitution.score(first_seed, second_seed),
            upstream_run.as_ref().map(|run| run.summary),
            downstream_run.as_ref().map(|run| run.summary),
        )?;
        let upstream_operations = upstream_run.as_ref().map_or(0, |run| run.operations.len());
        let downstream_operations = downstream_run
            .as_ref()
            .map_or(0, |run| run.operations.len());
        let operation_count = upstream_operations
            .checked_add(1)
            .and_then(|count| count.checked_add(downstream_operations))
            .ok_or(XDropError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
        let mut operations = Vec::new();
        operations.try_reserve_exact(operation_count).map_err(|_| {
            XDropError::AllocationFailed {
                buffer: "seeded X-drop alignment operations",
                elements: operation_count,
            }
        })?;
        if let Some(run) = upstream_run.as_mut() {
            run.operations.reverse();
            operations.extend_from_slice(&run.operations);
        }
        operations.push(if first_seed == second_seed {
            AlignmentOperation::Match
        } else {
            AlignmentOperation::Substitution
        });
        if let Some(run) = downstream_run {
            operations.extend_from_slice(&run.operations);
        }

        let alignment = PairwiseAlignment::try_new(
            summary.score,
            summary.first_range(),
            summary.second_range(),
            operations,
        )
        .map_err(XDropError::InvalidResult)?;
        Ok(XDropSeedAlignment { alignment, summary })
    }

    fn align_view<A>(
        &self,
        first: SymbolView<'_, A>,
        second: SymbolView<'_, A>,
    ) -> Result<ForwardAlignment, XDropError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let maximum_steps =
            first
                .len()
                .checked_add(second.len())
                .ok_or(XDropError::DimensionOverflow {
                    first_len: first.len(),
                    second_len: second.len(),
                })?;
        let mut history = TraceHistory::new(maximum_steps, first.len(), second.len())?;
        let summary = fill_forward::<A, S, true>(
            first,
            second,
            &self.substitution,
            self.gap,
            self.threshold,
            Some(&mut history),
        )?;
        let operation_capacity = summary.first_end.checked_add(summary.second_end).ok_or(
            XDropError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            },
        )?;
        let mut operations = Vec::new();
        operations
            .try_reserve_exact(operation_capacity)
            .map_err(|_| XDropError::AllocationFailed {
                buffer: "X-drop alignment operations",
                elements: operation_capacity,
            })?;
        traceback(
            first,
            second,
            summary.first_end,
            summary.second_end,
            &history,
            &mut operations,
        )?;
        Ok(ForwardAlignment {
            summary,
            operations,
        })
    }
}

#[derive(Clone, Copy)]
struct SymbolView<'a, A: Alphabet> {
    symbols: &'a [Symbol<A>],
    reversed: bool,
}

impl<'a, A: Alphabet> SymbolView<'a, A> {
    const fn forward(symbols: &'a [Symbol<A>]) -> Self {
        Self {
            symbols,
            reversed: false,
        }
    }

    const fn reversed(symbols: &'a [Symbol<A>]) -> Self {
        Self {
            symbols,
            reversed: true,
        }
    }

    const fn len(self) -> usize {
        self.symbols.len()
    }

    fn symbol(self, index: usize) -> Symbol<A> {
        if self.reversed {
            self.symbols[self.symbols.len() - 1 - index]
        } else {
            self.symbols[index]
        }
    }
}

struct ForwardAlignment {
    summary: XDropScore,
    operations: Vec<AlignmentOperation>,
}

fn validate_seed<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    seed: XDropSeed,
) -> Result<(Symbol<A>, Symbol<A>), XDropError> {
    match (first.get(seed.first_index), second.get(seed.second_index)) {
        (Some(first_seed), Some(second_seed)) => Ok((first_seed, second_seed)),
        _ => Err(XDropError::SeedOutOfBounds {
            first_index: seed.first_index,
            second_index: seed.second_index,
            first_len: first.len(),
            second_len: second.len(),
        }),
    }
}

fn build_seed_summary(
    first_len: usize,
    second_len: usize,
    seed: XDropSeed,
    seed_score: Score,
    upstream: Option<XDropScore>,
    downstream: Option<XDropScore>,
) -> Result<XDropSeedScore, XDropError> {
    let dimension_error = || XDropError::DimensionOverflow {
        first_len,
        second_len,
    };
    let upstream_first = upstream.map_or(0, XDropScore::first_end);
    let upstream_second = upstream.map_or(0, XDropScore::second_end);
    let downstream_first = downstream.map_or(0, XDropScore::first_end);
    let downstream_second = downstream.map_or(0, XDropScore::second_end);
    let first_start = seed
        .first_index
        .checked_sub(upstream_first)
        .ok_or_else(dimension_error)?;
    let second_start = seed
        .second_index
        .checked_sub(upstream_second)
        .ok_or_else(dimension_error)?;
    let first_end = seed
        .first_index
        .checked_add(1)
        .and_then(|value| value.checked_add(downstream_first))
        .ok_or_else(dimension_error)?;
    let second_end = seed
        .second_index
        .checked_add(1)
        .and_then(|value| value.checked_add(downstream_second))
        .ok_or_else(dimension_error)?;

    let mut score = seed_score;
    for component in [upstream, downstream].into_iter().flatten() {
        score = score
            .checked_add(component.score)
            .ok_or(XDropError::ScoreOverflow {
                first_prefix: seed.first_index,
                second_prefix: seed.second_index,
            })?;
    }

    let mut evaluated_cells = 1_usize;
    let mut retained_cells = 1_usize;
    let mut antidiagonals = 1_usize;
    let mut maximum_frontier_width = 1_usize;
    for component in [upstream, downstream].into_iter().flatten() {
        evaluated_cells = evaluated_cells
            .checked_add(component.evaluated_cells - 1)
            .ok_or_else(dimension_error)?;
        retained_cells = retained_cells
            .checked_add(component.retained_cells - 1)
            .ok_or_else(dimension_error)?;
        antidiagonals = antidiagonals
            .checked_add(component.antidiagonals - 1)
            .ok_or_else(dimension_error)?;
        maximum_frontier_width = maximum_frontier_width.max(component.maximum_frontier_width);
    }

    Ok(XDropSeedScore {
        score,
        first_start,
        first_end,
        second_start,
        second_end,
        evaluated_cells,
        retained_cells,
        antidiagonals,
        maximum_frontier_width,
        upstream,
        downstream,
    })
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct ReachableScore {
    value: Score,
    reachable: bool,
}

impl ReachableScore {
    const UNREACHABLE: Self = Self {
        value: 0,
        reachable: false,
    };

    const fn new(value: Score) -> Self {
        Self {
            value,
            reachable: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct AffineCell {
    best: ReachableScore,
    deletion: ReachableScore,
    insertion: ReachableScore,
}

impl AffineCell {
    const UNREACHABLE: Self = Self {
        best: ReachableScore::UNREACHABLE,
        deletion: ReachableScore::UNREACHABLE,
        insertion: ReachableScore::UNREACHABLE,
    };

    const ORIGIN: Self = Self {
        best: ReachableScore::new(0),
        deletion: ReachableScore::UNREACHABLE,
        insertion: ReachableScore::UNREACHABLE,
    };
}

struct Frontier {
    lower: usize,
    cells: Vec<AffineCell>,
    active: Option<(usize, usize)>,
    buffer: &'static str,
}

impl Frontier {
    fn new(_maximum_capacity: usize, buffer: &'static str) -> Result<Self, XDropError> {
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(1)
            .map_err(|_| XDropError::AllocationFailed {
                buffer,
                elements: 1,
            })?;
        Ok(Self {
            lower: 0,
            cells,
            active: None,
            buffer,
        })
    }

    fn set_origin(&mut self) {
        self.lower = 0;
        self.cells.clear();
        self.cells.push(AffineCell::ORIGIN);
        self.active = Some((0, 0));
    }

    fn prepare(&mut self, lower: usize, length: usize) -> Result<(), XDropError> {
        self.lower = lower;
        self.cells.clear();
        if self.cells.capacity() < length {
            self.cells
                .try_reserve_exact(length)
                .map_err(|_| XDropError::AllocationFailed {
                    buffer: self.buffer,
                    elements: length,
                })?;
        }
        self.cells.resize(length, AffineCell::UNREACHABLE);
        self.active = None;
        Ok(())
    }

    fn get(&self, first_prefix: usize) -> AffineCell {
        first_prefix
            .checked_sub(self.lower)
            .and_then(|index| self.cells.get(index))
            .copied()
            .unwrap_or(AffineCell::UNREACHABLE)
    }

    fn get_mut(&mut self, first_prefix: usize) -> Option<&mut AffineCell> {
        first_prefix
            .checked_sub(self.lower)
            .and_then(|index| self.cells.get_mut(index))
    }

    fn finish(&mut self) {
        let first = self.cells.iter().position(|cell| cell.best.reachable);
        let last = self.cells.iter().rposition(|cell| cell.best.reachable);
        self.active = first
            .zip(last)
            .map(|(first, last)| (self.lower + first, self.lower + last));
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum BestTrace {
    Start = 0,
    Diagonal = 1,
    Deletion = 2,
    Insertion = 3,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GapTrace {
    Open,
    Extend,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(transparent)]
struct TraceCell(u8);

impl TraceCell {
    const BEST_MASK: u8 = 0b0000_0011;
    const DELETION_EXTENDS: u8 = 0b0000_0100;
    const INSERTION_EXTENDS: u8 = 0b0000_1000;
    const RETAINED: u8 = 0b1000_0000;

    const fn retained(mut self) -> Self {
        self.0 |= Self::RETAINED;
        self
    }

    const fn is_retained(self) -> bool {
        self.0 & Self::RETAINED != 0
    }

    const fn with_best(mut self, trace: BestTrace) -> Self {
        self.0 = (self.0 & !Self::BEST_MASK) | trace as u8;
        self
    }

    const fn with_deletion(mut self, trace: GapTrace) -> Self {
        if matches!(trace, GapTrace::Extend) {
            self.0 |= Self::DELETION_EXTENDS;
        }
        self
    }

    const fn with_insertion(mut self, trace: GapTrace) -> Self {
        if matches!(trace, GapTrace::Extend) {
            self.0 |= Self::INSERTION_EXTENDS;
        }
        self
    }

    const fn best(self) -> BestTrace {
        match self.0 & Self::BEST_MASK {
            0 => BestTrace::Start,
            1 => BestTrace::Diagonal,
            2 => BestTrace::Deletion,
            3 => BestTrace::Insertion,
            _ => unreachable!(),
        }
    }

    const fn deletion(self) -> GapTrace {
        if self.0 & Self::DELETION_EXTENDS == 0 {
            GapTrace::Open
        } else {
            GapTrace::Extend
        }
    }

    const fn insertion(self) -> GapTrace {
        if self.0 & Self::INSERTION_EXTENDS == 0 {
            GapTrace::Open
        } else {
            GapTrace::Extend
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TracebackState {
    Best,
    Deletion,
    Insertion,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct TraceFrontier {
    offset: usize,
    lower: usize,
    length: usize,
}

struct TraceHistory {
    cells: Vec<TraceCell>,
    frontiers: Vec<TraceFrontier>,
    first_len: usize,
    second_len: usize,
}

impl TraceHistory {
    fn new(maximum_steps: usize, first_len: usize, second_len: usize) -> Result<Self, XDropError> {
        let entries = maximum_steps
            .checked_add(1)
            .ok_or(XDropError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(1)
            .map_err(|_| XDropError::AllocationFailed {
                buffer: "X-drop traceback cells",
                elements: 1,
            })?;
        cells.push(TraceCell::default().retained().with_best(BestTrace::Start));
        let mut frontiers = Vec::new();
        frontiers
            .try_reserve_exact(entries)
            .map_err(|_| XDropError::AllocationFailed {
                buffer: "X-drop traceback metadata",
                elements: entries,
            })?;
        frontiers.push(TraceFrontier {
            offset: 0,
            lower: 0,
            length: 1,
        });
        Ok(Self {
            cells,
            frontiers,
            first_len,
            second_len,
        })
    }

    fn start_frontier(&mut self, lower: usize, length: usize) -> Result<usize, XDropError> {
        let offset = self.cells.len();
        let new_len = offset
            .checked_add(length)
            .ok_or(XDropError::DimensionOverflow {
                first_len: self.first_len,
                second_len: self.second_len,
            })?;
        self.cells
            .try_reserve(length)
            .map_err(|_| XDropError::AllocationFailed {
                buffer: "X-drop traceback cells",
                elements: length,
            })?;
        self.cells.resize(new_len, TraceCell::default());
        self.frontiers.push(TraceFrontier {
            offset,
            lower,
            length,
        });
        Ok(offset)
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_forward<A, S, const TRACE: bool>(
    first: SymbolView<'_, A>,
    second: SymbolView<'_, A>,
    substitution: &S,
    gap: AffineGap,
    threshold: Score,
    mut history: Option<&mut TraceHistory>,
) -> Result<XDropScore, XDropError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_len = first.len();
    let second_len = second.len();
    let maximum_steps = first_len
        .checked_add(second_len)
        .ok_or(XDropError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    let maximum_width =
        first_len
            .min(second_len)
            .checked_add(1)
            .ok_or(XDropError::DimensionOverflow {
                first_len,
                second_len,
            })?;
    first_len
        .checked_add(1)
        .and_then(|rows| {
            second_len
                .checked_add(1)
                .and_then(|columns| rows.checked_mul(columns))
        })
        .ok_or(XDropError::DimensionOverflow {
            first_len,
            second_len,
        })?;

    let mut previous2 = Frontier::new(maximum_width, "X-drop second-prior frontier")?;
    let mut previous1 = Frontier::new(maximum_width, "X-drop prior frontier")?;
    let mut current = Frontier::new(maximum_width, "X-drop current frontier")?;
    previous1.set_origin();

    let mut endpoint = Endpoint {
        score: 0,
        first: 0,
        second: 0,
    };
    let mut global_best = 0;
    let mut evaluated_cells = 1_usize;
    let mut retained_cells = 1_usize;
    let mut antidiagonals = 1_usize;
    let mut maximum_frontier_width = 1_usize;
    let mut termination = XDropTermination::MatrixExhausted;

    for step in 1..=maximum_steps {
        let Some((lower, upper)) =
            candidate_bounds(&previous1, &previous2, step, first_len, second_len)?
        else {
            termination = XDropTermination::FrontiersExhausted;
            break;
        };
        let length = upper - lower + 1;
        current.prepare(lower, length)?;
        maximum_frontier_width = maximum_frontier_width.max(length);
        evaluated_cells =
            evaluated_cells
                .checked_add(length)
                .ok_or(XDropError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        antidiagonals += 1;

        let trace_offset = if TRACE {
            history
                .as_deref_mut()
                .expect("trace-specialized X-drop fill requires history")
                .start_frontier(lower, length)?
        } else {
            0
        };
        let mut next_global_best = global_best;

        for first_prefix in lower..=upper {
            let second_prefix = step - first_prefix;
            let diagonal_base = if first_prefix > 0 {
                previous2.get(first_prefix - 1).best
            } else {
                ReachableScore::UNREACHABLE
            };
            let upper_cell = if first_prefix > 0 {
                previous1.get(first_prefix - 1)
            } else {
                AffineCell::UNREACHABLE
            };
            let left_cell = previous1.get(first_prefix);
            let (candidate, trace) = recurrence_cell(
                first,
                second,
                substitution,
                gap,
                first_prefix,
                second_prefix,
                diagonal_base,
                upper_cell,
                left_cell,
            )?;
            if candidate.best.reachable
                && within_threshold(candidate.best.value, global_best, threshold)
            {
                let slot = current
                    .get_mut(first_prefix)
                    .expect("candidate coordinate lies in prepared frontier");
                *slot = candidate;
                retained_cells =
                    retained_cells
                        .checked_add(1)
                        .ok_or(XDropError::DimensionOverflow {
                            first_len,
                            second_len,
                        })?;
                endpoint.consider(candidate.best.value, first_prefix, second_prefix);
                next_global_best = next_global_best.max(candidate.best.value);
                if TRACE {
                    let index = trace_offset.checked_add(first_prefix - lower).ok_or(
                        XDropError::DimensionOverflow {
                            first_len,
                            second_len,
                        },
                    )?;
                    history
                        .as_deref_mut()
                        .expect("trace-specialized X-drop fill requires history")
                        .cells[index] = trace;
                }
            }
        }
        current.finish();
        global_best = next_global_best;
        core::mem::swap(&mut previous2, &mut previous1);
        core::mem::swap(&mut previous1, &mut current);
    }

    Ok(XDropScore {
        score: endpoint.score,
        first_end: endpoint.first,
        second_end: endpoint.second,
        evaluated_cells,
        retained_cells,
        antidiagonals,
        maximum_frontier_width,
        termination,
    })
}

fn candidate_bounds(
    previous1: &Frontier,
    previous2: &Frontier,
    step: usize,
    first_len: usize,
    second_len: usize,
) -> Result<Option<(usize, usize)>, XDropError> {
    let dimension_error = || XDropError::DimensionOverflow {
        first_len,
        second_len,
    };
    let legal_lower = step.saturating_sub(second_len);
    let legal_upper = first_len.min(step);
    let mut bounds: Option<(usize, usize)> = None;
    if let Some((lower, upper)) = previous1.active {
        merge_legal_interval(
            &mut bounds,
            lower,
            upper.checked_add(1).ok_or_else(dimension_error)?,
            legal_lower,
            legal_upper,
        );
    }
    if let Some((lower, upper)) = previous2.active {
        let diagonal_lower = lower.checked_add(1).ok_or_else(dimension_error)?;
        let diagonal_upper = upper.checked_add(1).ok_or_else(dimension_error)?;
        merge_legal_interval(
            &mut bounds,
            diagonal_lower,
            diagonal_upper,
            legal_lower,
            legal_upper,
        );
    }
    Ok(bounds)
}

fn merge_legal_interval(
    bounds: &mut Option<(usize, usize)>,
    lower: usize,
    upper: usize,
    legal_lower: usize,
    legal_upper: usize,
) {
    let clipped_lower = lower.max(legal_lower);
    let clipped_upper = upper.min(legal_upper);
    if clipped_lower > clipped_upper {
        return;
    }
    *bounds = Some(match *bounds {
        Some((current_lower, current_upper)) => (
            current_lower.min(clipped_lower),
            current_upper.max(clipped_upper),
        ),
        None => (clipped_lower, clipped_upper),
    });
}

#[allow(clippy::too_many_arguments)]
fn recurrence_cell<A, S>(
    first: SymbolView<'_, A>,
    second: SymbolView<'_, A>,
    substitution: &S,
    gap: AffineGap,
    first_prefix: usize,
    second_prefix: usize,
    diagonal_base: ReachableScore,
    upper_cell: AffineCell,
    left_cell: AffineCell,
) -> Result<(AffineCell, TraceCell), XDropError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let (deletion, deletion_trace) = if first_prefix > 0 {
        select_gap_state(
            checked_candidate(
                upper_cell.best,
                gap.open_score(),
                first_prefix,
                second_prefix,
            )?,
            checked_candidate(
                upper_cell.deletion,
                gap.extension_score(),
                first_prefix,
                second_prefix,
            )?,
        )
    } else {
        (ReachableScore::UNREACHABLE, None)
    };
    let (insertion, insertion_trace) = if second_prefix > 0 {
        select_gap_state(
            checked_candidate(
                left_cell.best,
                gap.open_score(),
                first_prefix,
                second_prefix,
            )?,
            checked_candidate(
                left_cell.insertion,
                gap.extension_score(),
                first_prefix,
                second_prefix,
            )?,
        )
    } else {
        (ReachableScore::UNREACHABLE, None)
    };
    let diagonal = if first_prefix > 0 && second_prefix > 0 {
        checked_candidate(
            diagonal_base,
            substitution.score(
                first.symbol(first_prefix - 1),
                second.symbol(second_prefix - 1),
            ),
            first_prefix,
            second_prefix,
        )?
    } else {
        ReachableScore::UNREACHABLE
    };
    let (best, best_trace) = select_best_state(diagonal, deletion, insertion);
    let mut trace = TraceCell::default();
    if best.reachable {
        trace = trace.retained().with_best(best_trace);
        if let Some(gap_trace) = deletion_trace {
            trace = trace.with_deletion(gap_trace);
        }
        if let Some(gap_trace) = insertion_trace {
            trace = trace.with_insertion(gap_trace);
        }
    }
    Ok((
        AffineCell {
            best,
            deletion,
            insertion,
        },
        trace,
    ))
}

fn checked_candidate(
    base: ReachableScore,
    increment: Score,
    first_prefix: usize,
    second_prefix: usize,
) -> Result<ReachableScore, XDropError> {
    if !base.reachable {
        return Ok(ReachableScore::UNREACHABLE);
    }
    if let Some(value) = base.value.checked_add(increment) {
        return Ok(ReachableScore::new(value));
    }
    if increment > 0 {
        Err(XDropError::ScoreOverflow {
            first_prefix,
            second_prefix,
        })
    } else {
        Ok(ReachableScore::UNREACHABLE)
    }
}

fn select_gap_state(
    open: ReachableScore,
    extension: ReachableScore,
) -> (ReachableScore, Option<GapTrace>) {
    if extension.reachable && (!open.reachable || extension.value >= open.value) {
        (extension, Some(GapTrace::Extend))
    } else if open.reachable {
        (open, Some(GapTrace::Open))
    } else {
        (ReachableScore::UNREACHABLE, None)
    }
}

fn select_best_state(
    diagonal: ReachableScore,
    deletion: ReachableScore,
    insertion: ReachableScore,
) -> (ReachableScore, BestTrace) {
    let mut selected = (diagonal, BestTrace::Diagonal);
    if deletion.reachable && (!selected.0.reachable || deletion.value >= selected.0.value) {
        selected = (deletion, BestTrace::Deletion);
    }
    if insertion.reachable && (!selected.0.reachable || insertion.value >= selected.0.value) {
        selected = (insertion, BestTrace::Insertion);
    }
    selected
}

fn within_threshold(score: Score, global_best: Score, threshold: Score) -> bool {
    i128::from(global_best) - i128::from(score) <= i128::from(threshold)
}

#[derive(Clone, Copy, Debug)]
struct Endpoint {
    score: Score,
    first: usize,
    second: usize,
}

impl Endpoint {
    fn consider(&mut self, score: Score, first: usize, second: usize) {
        let total = first + second;
        let selected_total = self.first + self.second;
        let ordered = total == selected_total
            && (second > self.second || (second == self.second && first > self.first));
        if score > self.score || (score == self.score && (total > selected_total || ordered)) {
            *self = Self {
                score,
                first,
                second,
            };
        }
    }
}

#[allow(clippy::too_many_lines)]
fn traceback<A: Alphabet>(
    first: SymbolView<'_, A>,
    second: SymbolView<'_, A>,
    mut first_prefix: usize,
    mut second_prefix: usize,
    history: &TraceHistory,
    operations: &mut Vec<AlignmentOperation>,
) -> Result<(), XDropError> {
    let mut state = TracebackState::Best;
    while first_prefix != 0 || second_prefix != 0 {
        let step =
            first_prefix
                .checked_add(second_prefix)
                .ok_or(XDropError::TracebackInvariant {
                    first_prefix,
                    second_prefix,
                })?;
        let frontier =
            history
                .frontiers
                .get(step)
                .copied()
                .ok_or(XDropError::TracebackInvariant {
                    first_prefix,
                    second_prefix,
                })?;
        let local = first_prefix
            .checked_sub(frontier.lower)
            .filter(|index| *index < frontier.length)
            .ok_or(XDropError::TracebackInvariant {
                first_prefix,
                second_prefix,
            })?;
        let index = frontier
            .offset
            .checked_add(local)
            .ok_or(XDropError::TracebackInvariant {
                first_prefix,
                second_prefix,
            })?;
        let trace = history
            .cells
            .get(index)
            .copied()
            .filter(|cell| cell.is_retained())
            .ok_or(XDropError::TracebackInvariant {
                first_prefix,
                second_prefix,
            })?;

        match state {
            TracebackState::Best => match trace.best() {
                BestTrace::Diagonal if first_prefix > 0 && second_prefix > 0 => {
                    operations.push(
                        if first.symbol(first_prefix - 1) == second.symbol(second_prefix - 1) {
                            AlignmentOperation::Match
                        } else {
                            AlignmentOperation::Substitution
                        },
                    );
                    first_prefix -= 1;
                    second_prefix -= 1;
                }
                BestTrace::Deletion if first_prefix > 0 => {
                    state = TracebackState::Deletion;
                }
                BestTrace::Insertion if second_prefix > 0 => {
                    state = TracebackState::Insertion;
                }
                BestTrace::Start
                | BestTrace::Diagonal
                | BestTrace::Deletion
                | BestTrace::Insertion => {
                    return Err(XDropError::TracebackInvariant {
                        first_prefix,
                        second_prefix,
                    });
                }
            },
            TracebackState::Deletion if first_prefix > 0 => {
                operations.push(AlignmentOperation::Deletion);
                first_prefix -= 1;
                state = if matches!(trace.deletion(), GapTrace::Extend) {
                    TracebackState::Deletion
                } else {
                    TracebackState::Best
                };
            }
            TracebackState::Insertion if second_prefix > 0 => {
                operations.push(AlignmentOperation::Insertion);
                second_prefix -= 1;
                state = if matches!(trace.insertion(), GapTrace::Extend) {
                    TracebackState::Insertion
                } else {
                    TracebackState::Best
                };
            }
            TracebackState::Deletion | TracebackState::Insertion => {
                return Err(XDropError::TracebackInvariant {
                    first_prefix,
                    second_prefix,
                });
            }
        }
    }
    if !matches!(state, TracebackState::Best) {
        return Err(XDropError::TracebackInvariant {
            first_prefix,
            second_prefix,
        });
    }
    operations.reverse();
    Ok(())
}

#[cfg(test)]
mod tests {
    use align_core::{AffineGap, DnaSequence, MatchMismatch};

    use super::{
        AffineCell, Frontier, XDrop, XDropDirection, XDropSeed, candidate_bounds, within_threshold,
    };

    #[test]
    fn threshold_equality_survives_and_strict_excess_fails() {
        assert!(within_threshold(7, 10, 3));
        assert!(!within_threshold(6, 10, 3));
        assert!(within_threshold(i64::MIN, i64::MIN, 0));
    }

    #[test]
    fn candidate_union_keeps_second_prior_diagonal_successors() {
        let mut previous1 = Frontier::new(8, "test").expect("allocation succeeds");
        let mut previous2 = Frontier::new(8, "test").expect("allocation succeeds");
        previous1.lower = 2;
        previous1.cells = vec![AffineCell::ORIGIN, AffineCell::ORIGIN];
        previous1.active = Some((2, 3));
        previous2.lower = 0;
        previous2.cells = vec![AffineCell::ORIGIN];
        previous2.active = Some((0, 0));
        assert_eq!(
            candidate_bounds(&previous1, &previous2, 5, 8, 8).expect("dimensions fit"),
            Some((1, 4))
        );

        previous1.active = None;
        assert_eq!(
            candidate_bounds(&previous1, &previous2, 2, 8, 8).expect("dimensions fit"),
            Some((1, 1))
        );

        previous1.active = Some((7, 7));
        previous2.active = Some((3, 3));
        assert_eq!(
            candidate_bounds(&previous1, &previous2, 13, 7, 8).expect("dimensions fit"),
            Some((7, 7)),
            "an entirely out-of-matrix diagonal interval must not widen a legal gap interval"
        );
    }

    #[test]
    fn configuration_seed_validation_and_memory_are_explicit() {
        let gap = AffineGap::new(-5, -1).expect("valid affine gap");
        assert!(XDrop::new(MatchMismatch::new(2, -3), gap, -1).is_err());
        let aligner =
            XDrop::new(MatchMismatch::new(2, -3), gap, 10).expect("valid X-drop configuration");
        assert_eq!(aligner.threshold(), 10);
        assert_eq!(aligner.gap(), gap);

        let first = DnaSequence::from_ascii(b"AC").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"AC").expect("valid DNA");
        assert!(
            aligner
                .score_from_seed(&first, &second, XDropSeed::new(2, 0), XDropDirection::Both,)
                .is_err()
        );

        let memory =
            XDrop::<MatchMismatch>::memory_requirement(7, 11).expect("small dimensions fit");
        assert_eq!(memory.maximum_frontier_width(), 8);
        assert_eq!(memory.maximum_trace_cells(), 96);
        assert_eq!(memory.metadata_entries(), 19);
        assert_eq!(memory.maximum_operations(), 18);
    }
}
