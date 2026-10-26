use core::cmp::Reverse;
use core::mem::size_of;
use std::collections::{BinaryHeap, VecDeque};

use align_core::{AlignmentOperation, Alphabet, PairwiseAlignment, Score, Sequence};

use crate::{
    WavefrontError,
    compact_wavefront::{
        CompactScoreWorkspace, CompactTracebackSummary, compact_affine_penalty,
        compact_affine_traceback,
    },
};

const UNREACHABLE: isize = -1;
const MATCH_EXTENSION_WORD: usize = size_of::<u64>();
const SMALL_SCHEDULER_MAX_INCREMENT: Score = 63;
const SMALL_HISTORY_SLOTS: usize = 64;
const DENSE_SCHEDULER_MAX_INCREMENT: Score = 4_096;
const EXACT_BOUND_SCAN_MIN: usize = 256;
const EXACT_BOUND_SCAN_MAX: usize = 8_192;
const UNEQUAL_COMPACT_SCAN_MAX: usize = 2_048;
const UNEQUAL_COMPACT_MINIMUM_GAP_DENOMINATOR: usize = 4;
// Main-state completion pruning is profitable for sparse exact bounds and proven
// single-gap paths. On broad substitution fronts it removes only a thin boundary
// while adding a branch to every evaluated diagonal.
const MAIN_PRUNE_SUBSTITUTION_DENOMINATOR: usize = 64;

/// Positive penalties for the exact gap-affine wavefront profile.
///
/// A length-`k` gap costs `gap_opening + k * gap_extension`, following the
/// original WFA publication rather than the score-side `AffineGap` convention.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontAffinePenalties {
    mismatch: Score,
    gap_opening: Score,
    gap_extension: Score,
    gap_open_transition: Score,
}

impl WavefrontAffinePenalties {
    /// Validates and constructs WFA gap-affine penalties.
    ///
    /// # Errors
    ///
    /// The mismatch and extension penalties must be strictly positive, the opening
    /// penalty must be nonnegative, and `opening + extension` must fit [`Score`].
    pub fn new(
        mismatch: Score,
        gap_opening: Score,
        gap_extension: Score,
    ) -> Result<Self, WavefrontError> {
        if mismatch <= 0 {
            return Err(WavefrontError::NonPositiveMismatchPenalty { penalty: mismatch });
        }
        if gap_opening < 0 {
            return Err(WavefrontError::NegativeGapOpeningPenalty {
                penalty: gap_opening,
            });
        }
        if gap_extension <= 0 {
            return Err(WavefrontError::NonPositiveGapExtensionPenalty {
                penalty: gap_extension,
            });
        }
        let gap_open_transition =
            gap_opening
                .checked_add(gap_extension)
                .ok_or(WavefrontError::PenaltyOverflow {
                    penalty: gap_opening,
                })?;
        Ok(Self {
            mismatch,
            gap_opening,
            gap_extension,
            gap_open_transition,
        })
    }

    /// Returns the positive mismatch penalty.
    #[must_use]
    pub const fn mismatch(self) -> Score {
        self.mismatch
    }

    /// Returns the nonnegative gap-opening penalty.
    #[must_use]
    pub const fn gap_opening(self) -> Score {
        self.gap_opening
    }

    /// Returns the positive per-symbol gap-extension penalty.
    #[must_use]
    pub const fn gap_extension(self) -> Score {
        self.gap_extension
    }

    /// Computes the penalty of one gap run.
    ///
    /// A zero length returns zero. Positive lengths use `o + k*e`.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError::PenaltyOverflow`] if the penalty is not
    /// representable.
    pub fn gap_penalty(self, length: usize) -> Result<Score, WavefrontError> {
        if length == 0 {
            return Ok(0);
        }
        let length = Score::try_from(length).map_err(|_| WavefrontError::PenaltyOverflow {
            penalty: self.gap_opening,
        })?;
        self.gap_extension
            .checked_mul(length)
            .and_then(|extension| self.gap_opening.checked_add(extension))
            .ok_or(WavefrontError::PenaltyOverflow {
                penalty: self.gap_opening,
            })
    }
}

/// Execution summary for one exact wavefront calculation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontScore {
    penalty: Score,
    wavefronts: usize,
    computed_offsets: usize,
    source_wavefront_lookups: usize,
    source_offset_lookups: usize,
    match_extensions: usize,
    maximum_wavefront_width: usize,
    modeled_peak_bytes: usize,
}

impl WavefrontScore {
    /// Returns the exact nonnegative minimizing penalty.
    #[must_use]
    pub const fn penalty(self) -> Score {
        self.penalty
    }

    /// Returns the conventional maximizing alignment score, `-penalty`.
    #[must_use]
    pub const fn alignment_score(self) -> Score {
        -self.penalty
    }

    /// Returns the number of nonempty wavefronts computed.
    #[must_use]
    pub const fn wavefronts(self) -> usize {
        self.wavefronts
    }

    /// Returns the number of common diagonal slots evaluated across wavefronts.
    ///
    /// The affine profile computes three state offsets per common slot.
    #[must_use]
    pub const fn computed_offsets(self) -> usize {
        self.computed_offsets
    }

    /// Returns predecessor-wavefront lookups across scheduled penalty attempts.
    ///
    /// The affine recurrence resolves mismatch, gap-open, and gap-extension
    /// predecessors for each attempt, including attempts that produce no wavefront.
    #[must_use]
    pub const fn source_wavefront_lookups(self) -> usize {
        self.source_wavefront_lookups
    }

    /// Returns predecessor state-offset reads in the affine recurrence.
    ///
    /// Each computed diagonal cell reads five predecessor offsets: one mismatch,
    /// two gap-open, and two gap-extension candidates.
    #[must_use]
    pub const fn source_offset_lookups(self) -> usize {
        self.source_offset_lookups
    }

    /// Returns the number of zero-cost equal-symbol extensions performed.
    #[must_use]
    pub const fn match_extensions(self) -> usize {
        self.match_extensions
    }

    /// Returns the largest common diagonal interval computed.
    #[must_use]
    pub const fn maximum_wavefront_width(self) -> usize {
        self.maximum_wavefront_width
    }

    /// Returns modeled peak element-buffer bytes for this execution mode.
    ///
    /// Allocator metadata, collection headers, inputs, and result headers are
    /// excluded.
    #[must_use]
    pub const fn modeled_peak_bytes(self) -> usize {
        self.modeled_peak_bytes
    }
}

/// Checked worst-case element-buffer model up to a supplied terminal penalty.
///
/// The bound assumes every numeric penalty through the terminal penalty has a
/// nonempty wavefront at the maximum possible width. Weighted affine execution is
/// often substantially sparser. Allocator metadata, collection headers, inputs,
/// result headers, and stack scalars are excluded.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontMemoryEstimate {
    terminal_penalty: Score,
    legal_diagonals: usize,
    maximum_wavefronts: usize,
    maximum_score_wavefronts: usize,
    maximum_wavefront_width: usize,
    maximum_common_offsets: usize,
    maximum_state_offsets: usize,
    maximum_trace_cells: usize,
    maximum_operations: usize,
    score_peak_bytes: usize,
    alignment_peak_bytes: usize,
}

/// One forward logical operation run in an exact affine-wavefront traceback.
///
/// A run is equivalent to `count` consecutive per-column operations, but keeps
/// the traceback proportional to the number of operation changes rather than
/// the number of alignment columns.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontTraceRun {
    operation: AlignmentOperation,
    count: usize,
}

impl WavefrontTraceRun {
    pub(crate) const fn from_parts(operation: AlignmentOperation, count: usize) -> Self {
        Self { operation, count }
    }

    /// Returns the repeated alignment operation.
    #[must_use]
    pub const fn operation(self) -> AlignmentOperation {
        self.operation
    }

    /// Returns the positive number of consecutive operations in this run.
    #[must_use]
    pub const fn count(self) -> usize {
        self.count
    }
}

/// Reusable allocation arena for exact affine-wavefront executions.
///
/// A workspace retains wavefront state and traceback buffers between calls but
/// never retains input references or output alignments. It can therefore be
/// reused across alphabets, sequence lengths, and penalty configurations.
#[derive(Debug, Default)]
pub struct WavefrontWorkspace {
    compact_score: CompactScoreWorkspace,
    rolling_history: VecDeque<AffineWavefront>,
    trace_history: Vec<AffineWavefront>,
    recycled: Vec<AffineWavefront>,
    small_history_index: SmallPenaltyHistoryIndex,
    trace_runs: Vec<WavefrontTraceRun>,
    operations: Vec<AlignmentOperation>,
    pending_score_bound: Option<(Score, bool)>,
}

impl WavefrontWorkspace {
    /// Creates an empty workspace that allocates on first use.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            compact_score: CompactScoreWorkspace::new(),
            rolling_history: VecDeque::new(),
            trace_history: Vec::new(),
            recycled: Vec::new(),
            small_history_index: SmallPenaltyHistoryIndex::new(),
            trace_runs: Vec::new(),
            operations: Vec::new(),
            pending_score_bound: None,
        }
    }

    /// Releases all wavefront buffers retained by this workspace.
    pub fn release(&mut self) {
        self.compact_score.release();
        self.rolling_history = VecDeque::new();
        self.trace_history = Vec::new();
        self.recycled = Vec::new();
        self.small_history_index = SmallPenaltyHistoryIndex::new();
        self.trace_runs = Vec::new();
        self.operations = Vec::new();
        self.pending_score_bound = None;
    }

    /// Returns bytes retained by collection element buffers.
    ///
    /// Allocator metadata and collection headers are excluded.
    #[must_use]
    pub fn retained_bytes(&self) -> usize {
        let wavefront_storage = self
            .rolling_history
            .iter()
            .chain(&self.trace_history)
            .chain(&self.recycled)
            .map(|wavefront| {
                wavefront.offsets.capacity() * size_of::<isize>()
                    + wavefront.traces.capacity() * size_of::<AffineTrace>()
            })
            .sum::<usize>();
        wavefront_storage
            + self.compact_score.retained_bytes()
            + self.rolling_history.capacity() * size_of::<AffineWavefront>()
            + self.trace_history.capacity() * size_of::<AffineWavefront>()
            + self.recycled.capacity() * size_of::<AffineWavefront>()
            + self.trace_runs.capacity() * size_of::<WavefrontTraceRun>()
            + self.operations.capacity() * size_of::<AlignmentOperation>()
    }

    fn recycle_trace_history(&mut self, mut history: Vec<AffineWavefront>) {
        self.recycled.append(&mut history);
        self.trace_history = history;
    }
}

impl WavefrontMemoryEstimate {
    /// Returns the modeled terminal penalty.
    #[must_use]
    pub const fn terminal_penalty(self) -> Score {
        self.terminal_penalty
    }

    /// Returns the number of legal matrix diagonals.
    #[must_use]
    pub const fn legal_diagonals(self) -> usize {
        self.legal_diagonals
    }

    /// Returns the maximum number of numeric penalty wavefronts.
    #[must_use]
    pub const fn maximum_wavefronts(self) -> usize {
        self.maximum_wavefronts
    }

    /// Returns the maximum simultaneously retained score-only wavefronts.
    #[must_use]
    pub const fn maximum_score_wavefronts(self) -> usize {
        self.maximum_score_wavefronts
    }

    /// Returns the maximum common diagonal width of one wavefront.
    #[must_use]
    pub const fn maximum_wavefront_width(self) -> usize {
        self.maximum_wavefront_width
    }

    /// Returns the conservative sum of common diagonal slots.
    #[must_use]
    pub const fn maximum_common_offsets(self) -> usize {
        self.maximum_common_offsets
    }

    /// Returns the maximum stored scalar state offsets.
    #[must_use]
    pub const fn maximum_state_offsets(self) -> usize {
        self.maximum_state_offsets
    }

    /// Returns the maximum packed traceback cells.
    #[must_use]
    pub const fn maximum_trace_cells(self) -> usize {
        self.maximum_trace_cells
    }

    /// Returns the maximum output-operation capacity.
    #[must_use]
    pub const fn maximum_operations(self) -> usize {
        self.maximum_operations
    }

    /// Returns worst-case modeled score-only element-buffer bytes.
    #[must_use]
    pub const fn score_peak_bytes(self) -> usize {
        self.score_peak_bytes
    }

    /// Returns worst-case modeled traceback element-buffer bytes.
    #[must_use]
    pub const fn alignment_peak_bytes(self) -> usize {
        self.alignment_peak_bytes
    }
}

/// One exact global wavefront traceback and its execution summary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct WavefrontAlignment {
    alignment: PairwiseAlignment,
    summary: WavefrontScore,
}

/// Exact affine traceback borrowed from a [`WavefrontWorkspace`].
///
/// The operation slice remains valid until the workspace is mutably reused. This
/// mirrors allocator-owning native aligners while [`WavefrontAlignment`] remains
/// the independently owned result type.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontAlignmentView<'a> {
    first_len: usize,
    second_len: usize,
    operations: &'a [AlignmentOperation],
    summary: WavefrontScore,
}

/// Exact compact affine traceback borrowed from a [`WavefrontWorkspace`].
///
/// Runs are exposed in forward alignment order without expanding one value per
/// alignment column. The view remains valid until the workspace is mutably
/// reused.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontTracebackView<'a> {
    first_len: usize,
    second_len: usize,
    columns: usize,
    reverse_runs: &'a [WavefrontTraceRun],
    penalty: Score,
}

impl WavefrontTracebackView<'_> {
    /// Returns the exact nonnegative affine penalty.
    #[must_use]
    pub const fn penalty(&self) -> Score {
        self.penalty
    }

    /// Returns the conventional maximizing score, `-penalty`.
    #[must_use]
    pub const fn alignment_score(&self) -> Score {
        -self.penalty
    }

    /// Iterates compact runs in forward alignment order.
    #[must_use]
    pub fn runs(
        &self,
    ) -> impl DoubleEndedIterator<Item = WavefrontTraceRun> + ExactSizeIterator + '_ {
        self.reverse_runs.iter().rev().copied()
    }

    /// Returns the number of compact operation runs.
    #[must_use]
    pub const fn run_count(&self) -> usize {
        self.reverse_runs.len()
    }

    /// Returns the number of columns represented by all compact runs.
    #[must_use]
    pub const fn columns(&self) -> usize {
        self.columns
    }

    /// Returns the number of consumed symbols in the first sequence.
    #[must_use]
    pub const fn first_len(&self) -> usize {
        self.first_len
    }

    /// Returns the number of consumed symbols in the second sequence.
    #[must_use]
    pub const fn second_len(&self) -> usize {
        self.second_len
    }
}

impl WavefrontAlignmentView<'_> {
    /// Returns the exact nonnegative affine penalty.
    #[must_use]
    pub const fn penalty(&self) -> Score {
        self.summary.penalty
    }

    /// Returns the conventional maximizing score, `-penalty`.
    #[must_use]
    pub const fn alignment_score(&self) -> Score {
        self.summary.alignment_score()
    }

    /// Returns the full per-column operation path.
    #[must_use]
    pub const fn operations(&self) -> &[AlignmentOperation] {
        self.operations
    }

    /// Returns the number of consumed symbols in the first sequence.
    #[must_use]
    pub const fn first_len(&self) -> usize {
        self.first_len
    }

    /// Returns the number of consumed symbols in the second sequence.
    #[must_use]
    pub const fn second_len(&self) -> usize {
        self.second_len
    }

    /// Returns execution statistics for the alignment.
    #[must_use]
    pub const fn summary(&self) -> WavefrontScore {
        self.summary
    }
}

impl WavefrontAlignment {
    /// Borrows the complete global operation path.
    #[must_use]
    pub const fn alignment(&self) -> &PairwiseAlignment {
        &self.alignment
    }

    /// Consumes this result and returns its alignment.
    #[must_use]
    pub fn into_alignment(self) -> PairwiseAlignment {
        self.alignment
    }

    /// Returns the exact nonnegative penalty.
    #[must_use]
    pub const fn penalty(&self) -> Score {
        self.summary.penalty
    }

    /// Returns the score-only execution summary.
    #[must_use]
    pub const fn summary(&self) -> WavefrontScore {
        self.summary
    }
}

/// Exact global gap-affine wavefront alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontAffine {
    penalties: WavefrontAffinePenalties,
}

impl WavefrontAffine {
    /// Constructs an immutable aligner from validated WFA penalties.
    #[must_use]
    pub const fn new(penalties: WavefrontAffinePenalties) -> Self {
        Self { penalties }
    }

    /// Returns the configured penalties.
    #[must_use]
    pub const fn penalties(self) -> WavefrontAffinePenalties {
        self.penalties
    }

    /// Models worst-case buffers through a supplied terminal affine penalty.
    ///
    /// # Errors
    ///
    /// The terminal penalty must be nonnegative and all counts and byte totals must
    /// be representable.
    pub fn memory_requirement(
        self,
        first_len: usize,
        second_len: usize,
        terminal_penalty: Score,
    ) -> Result<WavefrontMemoryEstimate, WavefrontError> {
        let minimum_increment = self
            .penalties
            .mismatch
            .min(self.penalties.gap_open_transition)
            .min(self.penalties.gap_extension);
        let maximum_increment = self
            .penalties
            .mismatch
            .max(self.penalties.gap_open_transition)
            .max(self.penalties.gap_extension);
        let retained_score_wavefronts = usize::try_from(maximum_increment / minimum_increment)
            .ok()
            .and_then(|wavefronts| wavefronts.checked_add(2))
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        memory_requirement(
            first_len,
            second_len,
            terminal_penalty,
            minimum_increment,
            self.penalties.gap_extension,
            3,
            retained_score_wavefronts,
        )
    }

    /// Computes only the exact global gap-affine penalty.
    ///
    /// This lean path omits execution telemetry that is populated by [`Self::score`].
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for unrepresentable dimensions or penalties,
    /// failed allocations, or an internal exact-search invariant failure.
    pub fn penalty<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<Score, WavefrontError> {
        let mut workspace = WavefrontWorkspace::new();
        self.penalty_with_workspace(first, second, &mut workspace)
    }

    /// Computes only the exact affine penalty using a reusable allocation arena.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for unrepresentable dimensions or penalties,
    /// failed allocations, or an internal exact-search invariant failure.
    pub fn penalty_with_workspace<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut WavefrontWorkspace,
    ) -> Result<Score, WavefrontError> {
        run_affine_penalty(
            first.as_indices(),
            second.as_indices(),
            self.penalties,
            workspace,
        )
    }

    pub(crate) fn penalty_with_workspace_bound<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut WavefrontWorkspace,
        upper_bound: Score,
        prune_main: bool,
    ) -> Result<Score, WavefrontError> {
        run_affine_penalty_with_bound(
            first.as_indices(),
            second.as_indices(),
            self.penalties,
            workspace,
            upper_bound,
            prune_main,
        )
    }

    /// Computes the exact global gap-affine penalty and wavefront statistics.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for unrepresentable dimensions or penalties,
    /// failed allocations, or an internal exact-search invariant failure.
    pub fn score<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<WavefrontScore, WavefrontError> {
        let mut workspace = WavefrontWorkspace::new();
        self.score_with_workspace(first, second, &mut workspace)
    }

    /// Computes the exact affine penalty using a reusable allocation arena.
    ///
    /// This method is numerically identical to [`Self::score`].
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for unrepresentable dimensions or penalties,
    /// failed allocations, or an internal exact-search invariant failure.
    pub fn score_with_workspace<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut WavefrontWorkspace,
    ) -> Result<WavefrontScore, WavefrontError> {
        run_affine_score(
            first.as_indices(),
            second.as_indices(),
            self.penalties,
            workspace,
        )
    }

    /// Computes one deterministic exact global gap-affine alignment.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for configuration, dimension, allocation,
    /// traceback, score, or result-construction failures.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<WavefrontAlignment, WavefrontError> {
        let mut workspace = WavefrontWorkspace::new();
        self.align_with_workspace(first, second, &mut workspace)
    }

    /// Computes an exact affine alignment using a reusable allocation arena.
    ///
    /// The returned alignment owns its operations independently of `workspace`.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for configuration, dimension, allocation,
    /// traceback, score, or result-construction failures.
    pub fn align_with_workspace<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut WavefrontWorkspace,
    ) -> Result<WavefrontAlignment, WavefrontError> {
        let (upper_bound, _) =
            exact_alignment_upper_bound(first.as_indices(), second.as_indices(), self.penalties);
        self.align_with_workspace_bound(first, second, workspace, upper_bound)
    }

    pub(crate) fn align_with_workspace_bound<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut WavefrontWorkspace,
        upper_bound: Score,
    ) -> Result<WavefrontAlignment, WavefrontError> {
        let first_indices = first.as_indices();
        let second_indices = second.as_indices();
        let run = run_affine_with_bound::<true>(
            first_indices,
            second_indices,
            self.penalties,
            workspace,
            upper_bound,
        )?;
        if run.summary.penalty == 0 {
            workspace.recycle_trace_history(run.history);
            return build_identity_alignment(first_indices, 3);
        }
        build_affine_alignment(
            first_indices,
            second_indices,
            self.penalties,
            run,
            workspace,
        )
    }

    /// Computes an exact compact affine traceback retained inside `workspace`.
    ///
    /// The returned view yields run-length encoded operations in forward
    /// alignment order and never creates a per-column operation vector. It is
    /// invalidated by the next mutable use of `workspace`.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for configuration, dimension, allocation,
    /// traceback, score, or compact-result validation failures.
    pub fn traceback_in_workspace<'a, A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &'a mut WavefrontWorkspace,
    ) -> Result<WavefrontTracebackView<'a>, WavefrontError> {
        let prepared = prepare_affine_traceback(
            first.as_indices(),
            second.as_indices(),
            self.penalties,
            workspace,
        )?;
        Ok(WavefrontTracebackView {
            first_len: first.len(),
            second_len: second.len(),
            columns: prepared.columns,
            reverse_runs: &workspace.trace_runs,
            penalty: prepared.summary.penalty,
        })
    }

    /// Computes an exact affine traceback retained inside `workspace`.
    ///
    /// Unlike [`Self::align_with_workspace`], this method does not allocate a new
    /// independently owned operation vector on every call. The returned view
    /// borrows the workspace and is invalidated by its next mutable use.
    ///
    /// # Errors
    ///
    /// Returns [`WavefrontError`] for configuration, dimension, allocation,
    /// traceback, score, or result-construction failures.
    pub fn align_in_workspace<'a, A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &'a mut WavefrontWorkspace,
    ) -> Result<WavefrontAlignmentView<'a>, WavefrontError> {
        let (upper_bound, _) =
            exact_alignment_upper_bound(first.as_indices(), second.as_indices(), self.penalties);
        self.align_in_workspace_bound(first, second, workspace, upper_bound)
    }

    pub(crate) fn align_in_workspace_bound<'a, A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &'a mut WavefrontWorkspace,
        upper_bound: Score,
    ) -> Result<WavefrontAlignmentView<'a>, WavefrontError> {
        let prepared = prepare_affine_traceback_with_bound(
            first.as_indices(),
            second.as_indices(),
            self.penalties,
            workspace,
            upper_bound,
        )?;
        workspace.operations.clear();
        if workspace.operations.capacity() < prepared.dimensions.maximum_operations {
            workspace
                .operations
                .try_reserve_exact(prepared.dimensions.maximum_operations)
                .map_err(|_| WavefrontError::AllocationFailed {
                    buffer: "workspace WFA alignment operations",
                    elements: prepared.dimensions.maximum_operations,
                })?;
        }
        expand_trace_runs(
            &workspace.trace_runs,
            &mut workspace.operations,
            prepared.dimensions.maximum_operations,
            first.len(),
            second.len(),
        )?;
        let operation_bytes = checked_bytes(
            prepared.dimensions.maximum_operations,
            size_of::<AlignmentOperation>(),
            first.len(),
            second.len(),
        )?;
        let mut summary = prepared.summary;
        summary.modeled_peak_bytes = summary
            .modeled_peak_bytes
            .checked_add(operation_bytes)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
        Ok(WavefrontAlignmentView {
            first_len: first.len(),
            second_len: second.len(),
            operations: &workspace.operations,
            summary,
        })
    }
}

struct PreparedAffineTraceback {
    dimensions: Dimensions,
    columns: usize,
    summary: WavefrontScore,
}

fn prepare_affine_traceback(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
) -> Result<PreparedAffineTraceback, WavefrontError> {
    let (upper_bound, _) = exact_alignment_upper_bound(first, second, penalties);
    prepare_affine_traceback_with_bound(first, second, penalties, workspace, upper_bound)
}

fn prepare_identity_traceback(
    sequence_len: usize,
    dimensions: Dimensions,
    trace_runs: &mut Vec<WavefrontTraceRun>,
) -> Result<PreparedAffineTraceback, WavefrontError> {
    push_trace_run(
        trace_runs,
        AlignmentOperation::Match,
        sequence_len,
        sequence_len,
        sequence_len,
    )?;
    let state_bytes = checked_bytes(3, size_of::<isize>(), sequence_len, sequence_len)?;
    let trace_bytes = checked_bytes(1, size_of::<AffineTrace>(), sequence_len, sequence_len)?;
    let trace_run_bytes = checked_bytes(
        trace_runs.len(),
        size_of::<WavefrontTraceRun>(),
        sequence_len,
        sequence_len,
    )?;
    let modeled_peak_bytes = state_bytes
        .checked_add(trace_bytes)
        .and_then(|bytes| bytes.checked_add(trace_run_bytes))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len: sequence_len,
            second_len: sequence_len,
        })?;
    Ok(PreparedAffineTraceback {
        dimensions,
        columns: sequence_len,
        summary: WavefrontScore {
            penalty: 0,
            wavefronts: 1,
            computed_offsets: 1,
            source_wavefront_lookups: 0,
            source_offset_lookups: 0,
            match_extensions: sequence_len,
            maximum_wavefront_width: 1,
            modeled_peak_bytes,
        },
    })
}

fn prepare_compact_traceback(
    first_len: usize,
    second_len: usize,
    dimensions: Dimensions,
    trace_runs: &[WavefrontTraceRun],
    compact: CompactTracebackSummary,
) -> Result<PreparedAffineTraceback, WavefrontError> {
    let columns = trace_run_columns(
        trace_runs,
        first_len,
        second_len,
        dimensions.maximum_operations,
    )?;
    let trace_run_bytes = checked_bytes(
        trace_runs.len(),
        size_of::<WavefrontTraceRun>(),
        first_len,
        second_len,
    )?;
    let modeled_peak_bytes = compact
        .modeled_peak_bytes
        .checked_add(trace_run_bytes)
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    Ok(PreparedAffineTraceback {
        dimensions,
        columns,
        summary: WavefrontScore {
            penalty: compact.penalty,
            wavefronts: compact.wavefronts,
            computed_offsets: compact.computed_offsets,
            source_wavefront_lookups: compact.source_wavefront_lookups,
            source_offset_lookups: compact.source_offset_lookups,
            match_extensions: compact.match_extensions,
            maximum_wavefront_width: compact.maximum_wavefront_width,
            modeled_peak_bytes,
        },
    })
}

fn prepare_affine_traceback_with_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
    upper_bound: Score,
) -> Result<PreparedAffineTraceback, WavefrontError> {
    let dimensions = dimensions(first.len(), second.len())?;
    workspace.trace_runs.clear();
    if upper_bound == 0 && first == second {
        return prepare_identity_traceback(first.len(), dimensions, &mut workspace.trace_runs);
    }
    let compact = compact_affine_traceback(
        first,
        second,
        penalties,
        &mut workspace.compact_score,
        upper_bound,
        false,
        &mut workspace.trace_runs,
    )?;
    if let Some(compact) = compact {
        return prepare_compact_traceback(
            first.len(),
            second.len(),
            dimensions,
            &workspace.trace_runs,
            compact,
        );
    }

    let mut run = run_affine_with_bound::<true>(first, second, penalties, workspace, upper_bound)?;
    let dimensions = run.dimensions;
    if run.summary.penalty == 0 {
        push_trace_run(
            &mut workspace.trace_runs,
            AlignmentOperation::Match,
            first.len(),
            first.len(),
            second.len(),
        )?;
    } else {
        traceback_affine(
            first,
            second,
            dimensions,
            penalties,
            run.summary.penalty,
            &run.history,
            &mut workspace.trace_runs,
        )?;
    }
    let columns = trace_run_columns(
        &workspace.trace_runs,
        first.len(),
        second.len(),
        dimensions.maximum_operations,
    )?;
    workspace.recycle_trace_history(core::mem::take(&mut run.history));
    let trace_run_bytes = checked_bytes(
        workspace.trace_runs.len(),
        size_of::<WavefrontTraceRun>(),
        first.len(),
        second.len(),
    )?;
    run.summary.modeled_peak_bytes = run
        .summary
        .modeled_peak_bytes
        .checked_add(trace_run_bytes)
        .ok_or(WavefrontError::DimensionOverflow {
            first_len: first.len(),
            second_len: second.len(),
        })?;
    Ok(PreparedAffineTraceback {
        dimensions,
        columns,
        summary: run.summary,
    })
}

fn build_identity_alignment(
    sequence: &[u8],
    state_components: usize,
) -> Result<WavefrontAlignment, WavefrontError> {
    let operation_capacity =
        sequence
            .len()
            .checked_mul(2)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len: sequence.len(),
                second_len: sequence.len(),
            })?;
    let mut operations = Vec::new();
    operations
        .try_reserve_exact(sequence.len())
        .map_err(|_| WavefrontError::AllocationFailed {
            buffer: "identity WFA alignment operations",
            elements: sequence.len(),
        })?;
    operations.resize(sequence.len(), AlignmentOperation::Match);
    let alignment = PairwiseAlignment::try_new(0, 0..sequence.len(), 0..sequence.len(), operations)
        .map_err(WavefrontError::InvalidResult)?;
    let state_bytes = checked_bytes(
        state_components,
        size_of::<isize>(),
        sequence.len(),
        sequence.len(),
    )?;
    let trace_bytes = size_of::<AffineTrace>().max(1);
    let operation_bytes = checked_bytes(
        operation_capacity,
        size_of::<AlignmentOperation>(),
        sequence.len(),
        sequence.len(),
    )?;
    let modeled_peak_bytes = state_bytes
        .checked_add(trace_bytes)
        .and_then(|bytes| bytes.checked_add(operation_bytes))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len: sequence.len(),
            second_len: sequence.len(),
        })?;
    Ok(WavefrontAlignment {
        alignment,
        summary: WavefrontScore {
            penalty: 0,
            wavefronts: 1,
            computed_offsets: 1,
            source_wavefront_lookups: 0,
            source_offset_lookups: 0,
            match_extensions: sequence.len(),
            maximum_wavefront_width: 1,
            modeled_peak_bytes,
        },
    })
}

#[allow(clippy::too_many_lines)]
fn memory_requirement(
    first_len: usize,
    second_len: usize,
    terminal_penalty: Score,
    minimum_increment: Score,
    diagonal_step_penalty: Score,
    components: usize,
    score_wavefront_bound: usize,
) -> Result<WavefrontMemoryEstimate, WavefrontError> {
    if terminal_penalty < 0 {
        return Err(WavefrontError::NegativeTerminalPenalty {
            penalty: terminal_penalty,
        });
    }
    let dimensions = dimensions(first_len, second_len)?;
    let legal_diagonals =
        dimensions
            .maximum_operations
            .checked_add(1)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
    let maximum_wavefronts = usize::try_from(terminal_penalty)
        .ok()
        .and_then(|penalty| penalty.checked_add(1))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    let radius = usize::try_from(terminal_penalty / diagonal_step_penalty).map_err(|_| {
        WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        }
    })?;
    let maximum_wavefront_width = if radius >= legal_diagonals / 2 {
        legal_diagonals
    } else {
        radius
            .checked_mul(2)
            .and_then(|width| width.checked_add(1))
            .map(|width| width.min(legal_diagonals))
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?
    };
    let maximum_common_offsets = maximum_wavefronts
        .checked_mul(maximum_wavefront_width)
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    let maximum_state_offsets = maximum_common_offsets.checked_mul(components).ok_or(
        WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        },
    )?;
    let maximum_score_wavefronts = maximum_wavefronts.min(score_wavefront_bound);
    let score_state_offsets = maximum_wavefront_width
        .checked_mul(components)
        .and_then(|offsets| offsets.checked_mul(maximum_score_wavefronts))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    let score_peak_bytes = checked_bytes(
        score_state_offsets,
        size_of::<isize>(),
        first_len,
        second_len,
    )?;
    let history_bytes = checked_bytes(
        maximum_state_offsets,
        size_of::<isize>(),
        first_len,
        second_len,
    )?;
    let trace_bytes = checked_bytes(
        maximum_common_offsets,
        size_of::<AffineTrace>(),
        first_len,
        second_len,
    )?;
    let operation_bytes = checked_bytes(
        dimensions.maximum_operations,
        size_of::<AlignmentOperation>(),
        first_len,
        second_len,
    )?;
    let maximum_edits = usize::try_from(terminal_penalty / minimum_increment)
        .map_err(|_| WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?
        .min(dimensions.maximum_operations);
    let maximum_trace_runs = if maximum_edits == 0 {
        0
    } else {
        maximum_edits
            .checked_mul(2)
            .and_then(|runs| runs.checked_add(1))
            .map(|runs| runs.min(dimensions.maximum_operations))
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?
    };
    let trace_run_bytes = checked_bytes(
        maximum_trace_runs,
        size_of::<WavefrontTraceRun>(),
        first_len,
        second_len,
    )?;
    let alignment_peak_bytes = history_bytes
        .checked_add(trace_bytes)
        .and_then(|bytes| bytes.checked_add(trace_run_bytes))
        .and_then(|bytes| bytes.checked_add(operation_bytes))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    Ok(WavefrontMemoryEstimate {
        terminal_penalty,
        legal_diagonals,
        maximum_wavefronts,
        maximum_score_wavefronts,
        maximum_wavefront_width,
        maximum_common_offsets,
        maximum_state_offsets,
        maximum_trace_cells: maximum_common_offsets,
        maximum_operations: dimensions.maximum_operations,
        score_peak_bytes,
        alignment_peak_bytes,
    })
}

#[derive(Clone, Copy, Debug)]
struct Dimensions {
    first: isize,
    second: isize,
    target_diagonal: isize,
    maximum_operations: usize,
}

fn dimensions(first_len: usize, second_len: usize) -> Result<Dimensions, WavefrontError> {
    let maximum_operations =
        first_len
            .checked_add(second_len)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
    if maximum_operations > isize::MAX as usize {
        return Err(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        });
    }
    let first = isize::try_from(first_len).map_err(|_| WavefrontError::DimensionOverflow {
        first_len,
        second_len,
    })?;
    let second = isize::try_from(second_len).map_err(|_| WavefrontError::DimensionOverflow {
        first_len,
        second_len,
    })?;
    Ok(Dimensions {
        first,
        second,
        target_diagonal: second - first,
        maximum_operations,
    })
}

pub(crate) fn count_byte_mismatches(first: &[u8], second: &[u8]) -> usize {
    const MISMATCH_COUNT_WORD: usize = size_of::<u64>();
    const LOW_SEVEN_BITS: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    const HIGH_BITS: u64 = 0x8080_8080_8080_8080;

    let word_bytes = first.len() - first.len() % MISMATCH_COUNT_WORD;
    let first_words = first[..word_bytes].chunks_exact(MISMATCH_COUNT_WORD);
    let second_words = second[..word_bytes].chunks_exact(MISMATCH_COUNT_WORD);
    let mut mismatches = 0_usize;
    for (first_word, second_word) in first_words.zip(second_words) {
        let first_word = u64::from_le_bytes(
            first_word
                .try_into()
                .expect("an exact WFA bound word has eight bytes"),
        );
        let second_word = u64::from_le_bytes(
            second_word
                .try_into()
                .expect("an exact WFA bound word has eight bytes"),
        );
        let difference = first_word ^ second_word;
        if difference != 0 {
            let nonzero_bytes =
                (((difference & LOW_SEVEN_BITS) + LOW_SEVEN_BITS) | difference) & HIGH_BITS;
            mismatches += usize::try_from(nonzero_bytes.count_ones())
                .expect("a word mismatch count fits usize");
        }
    }
    mismatches
        + first[word_bytes..]
            .iter()
            .zip(&second[word_bytes..])
            .filter(|(left, right)| left != right)
            .count()
}

fn exact_gap_or_terminal_substitutions(shorter: &[u8], longer: &[u8]) -> usize {
    debug_assert!(shorter.len() < longer.len());
    let gap_length = longer.len() - shorter.len();
    let (common_prefix, _) = extend_equal_coordinates(shorter, longer, 0, 0);
    if shorter[common_prefix..] == longer[common_prefix + gap_length..] {
        return 0;
    }

    count_byte_mismatches(shorter, &longer[..shorter.len()])
}

fn exact_alignment_upper_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
) -> (Score, bool) {
    let common = first.len().min(second.len());
    if !(EXACT_BOUND_SCAN_MIN..=EXACT_BOUND_SCAN_MAX).contains(&common) {
        return (Score::MAX, false);
    }
    exact_alignment_upper_bound_unbounded(first, second, penalties)
}

pub(crate) fn exact_alignment_upper_bound_unbounded(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
) -> (Score, bool) {
    let common = first.len().min(second.len());
    let gap_length = first.len().abs_diff(second.len());
    let substitutions = if gap_length == 0 {
        count_byte_mismatches(first, second)
    } else if first.len() < second.len() {
        exact_gap_or_terminal_substitutions(first, second)
    } else {
        exact_gap_or_terminal_substitutions(second, first)
    };
    let Some(substitutions) = Score::try_from(substitutions).ok() else {
        return (Score::MAX, false);
    };
    let Some(substitution_penalty) = substitutions.checked_mul(penalties.mismatch) else {
        return (Score::MAX, false);
    };
    let prune_main = usize::try_from(substitutions)
        .ok()
        .is_some_and(|substitutions| {
            substitutions
                .checked_mul(MAIN_PRUNE_SUBSTITUTION_DENOMINATOR)
                .is_some_and(|scaled| scaled <= common)
        });
    if gap_length == 0 {
        return (substitution_penalty, prune_main);
    }
    let Some(gap_length) = Score::try_from(gap_length).ok() else {
        return (Score::MAX, false);
    };
    let penalty = gap_length
        .checked_mul(penalties.gap_extension)
        .and_then(|extension| extension.checked_add(penalties.gap_opening))
        .and_then(|gap| substitution_penalty.checked_add(gap))
        .unwrap_or(Score::MAX);
    (penalty, prune_main && penalty != Score::MAX)
}

fn unequal_compact_upper_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
) -> Option<(Score, bool)> {
    let common = first.len().min(second.len());
    let gap = first.len().abs_diff(second.len());
    if first.len() == second.len()
        || !(EXACT_BOUND_SCAN_MIN..=UNEQUAL_COMPACT_SCAN_MAX).contains(&common)
        || gap < common.div_ceil(UNEQUAL_COMPACT_MINIMUM_GAP_DENOMINATOR)
    {
        return None;
    }
    Some(exact_alignment_upper_bound_unbounded(
        first, second, penalties,
    ))
}

fn constrain_diagonal_bounds(
    lower: &mut isize,
    upper: &mut isize,
    penalty: Score,
    upper_bound: Score,
    gap_extension: Score,
    dimensions: Dimensions,
) {
    if upper_bound == Score::MAX {
        return;
    }
    let remaining = upper_bound - penalty;
    if remaining < 0 {
        *lower = 1;
        *upper = 0;
        return;
    }
    let maximum_shift = isize::try_from(remaining / gap_extension).unwrap_or(isize::MAX);
    *lower = (*lower).max(dimensions.target_diagonal.saturating_sub(maximum_shift));
    *upper = (*upper).min(dimensions.target_diagonal.saturating_add(maximum_shift));
}

fn legal_coordinate(
    offset: isize,
    diagonal: isize,
    dimensions: Dimensions,
) -> Option<(usize, usize)> {
    if offset < 0 || offset > dimensions.second {
        return None;
    }
    let vertical = offset.checked_sub(diagonal)?;
    if vertical < 0 || vertical > dimensions.first {
        return None;
    }
    Some((
        usize::try_from(vertical).ok()?,
        usize::try_from(offset).ok()?,
    ))
}

fn extend_equal(
    first: &[u8],
    second: &[u8],
    diagonal: isize,
    offset: isize,
    dimensions: Dimensions,
) -> Option<(isize, usize)> {
    let (vertical, horizontal) = legal_coordinate(offset, diagonal, dimensions)?;
    let (horizontal, matches) = extend_equal_coordinates(first, second, vertical, horizontal);
    Some((isize::try_from(horizontal).ok()?, matches))
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn extend_equal_coordinates(
    first: &[u8],
    second: &[u8],
    mut vertical: usize,
    mut horizontal: usize,
) -> (usize, usize) {
    let mut matches = 0_usize;
    let remaining = (first.len() - vertical).min(second.len() - horizontal);
    let word_bytes = remaining - remaining % MATCH_EXTENSION_WORD;
    let first_words = first[vertical..vertical + word_bytes].chunks_exact(MATCH_EXTENSION_WORD);
    let second_words =
        second[horizontal..horizontal + word_bytes].chunks_exact(MATCH_EXTENSION_WORD);
    for (word_index, (first_word, second_word)) in first_words.zip(second_words).enumerate() {
        let first_word = u64::from_le_bytes(
            first_word
                .try_into()
                .expect("an exact WFA word chunk has eight bytes"),
        );
        let second_word = u64::from_le_bytes(
            second_word
                .try_into()
                .expect("an exact WFA word chunk has eight bytes"),
        );
        let difference = first_word ^ second_word;
        if difference != 0 {
            let matching_bytes = usize::try_from(difference.trailing_zeros() / u8::BITS)
                .expect("a word prefix length fits usize");
            let matched_words = word_index * MATCH_EXTENSION_WORD;
            return (
                horizontal + matched_words + matching_bytes,
                matches + matched_words + matching_bytes,
            );
        }
    }
    vertical += word_bytes;
    horizontal += word_bytes;
    matches += word_bytes;
    while vertical < first.len()
        && horizontal < second.len()
        && first[vertical] == second[horizontal]
    {
        vertical += 1;
        horizontal += 1;
        matches += 1;
    }
    (horizontal, matches)
}

#[inline]
fn verified_match_length(
    dimensions: Dimensions,
    diagonal: isize,
    source_offset: isize,
    offset: &mut isize,
    penalty: Score,
) -> Result<usize, WavefrontError> {
    if source_offset < 0 || source_offset > *offset {
        return Err(WavefrontError::TracebackInvariant {
            penalty,
            diagonal,
            offset: *offset,
        });
    }
    let Some((source_vertical, source_horizontal)) =
        legal_coordinate(source_offset, diagonal, dimensions)
    else {
        return Err(WavefrontError::TracebackInvariant {
            penalty,
            diagonal,
            offset: *offset,
        });
    };
    let Some((target_vertical, target_horizontal)) =
        legal_coordinate(*offset, diagonal, dimensions)
    else {
        return Err(WavefrontError::TracebackInvariant {
            penalty,
            diagonal,
            offset: *offset,
        });
    };
    if source_vertical > target_vertical
        || source_horizontal > target_horizontal
        || target_vertical - source_vertical != target_horizontal - source_horizontal
    {
        return Err(WavefrontError::TracebackInvariant {
            penalty,
            diagonal,
            offset: *offset,
        });
    }
    let match_length = target_horizontal - source_horizontal;
    *offset = source_offset;
    Ok(match_length)
}

#[inline]
fn push_trace_run(
    runs: &mut Vec<WavefrontTraceRun>,
    operation: AlignmentOperation,
    count: usize,
    first_len: usize,
    second_len: usize,
) -> Result<(), WavefrontError> {
    if count == 0 {
        return Ok(());
    }
    if let Some(last) = runs.last_mut().filter(|run| run.operation == operation) {
        last.count = last
            .count
            .checked_add(count)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        return Ok(());
    }
    if runs.len() == runs.capacity() {
        runs.try_reserve(1)
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "affine WFA traceback runs",
                elements: runs.len().saturating_add(1),
            })?;
    }
    runs.push(WavefrontTraceRun { operation, count });
    Ok(())
}

#[inline]
fn trace_run_columns(
    runs: &[WavefrontTraceRun],
    first_len: usize,
    second_len: usize,
    maximum_operations: usize,
) -> Result<usize, WavefrontError> {
    runs.iter().try_fold(0_usize, |columns, run| {
        columns
            .checked_add(run.count)
            .filter(|columns| *columns <= maximum_operations)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })
    })
}

fn expand_trace_runs(
    runs: &[WavefrontTraceRun],
    operations: &mut Vec<AlignmentOperation>,
    maximum_operations: usize,
    first_len: usize,
    second_len: usize,
) -> Result<(), WavefrontError> {
    for run in runs.iter().rev() {
        let new_length = operations
            .len()
            .checked_add(run.count)
            .filter(|length| *length <= maximum_operations && *length <= operations.capacity())
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        operations.resize(new_length, run.operation);
    }
    Ok(())
}

fn checked_width(
    lower: isize,
    upper: isize,
    first_len: usize,
    second_len: usize,
) -> Result<usize, WavefrontError> {
    if lower > upper {
        return Ok(0);
    }
    let span = upper
        .checked_sub(lower)
        .and_then(|value| value.checked_add(1))
        .and_then(|value| usize::try_from(value).ok())
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })?;
    Ok(span)
}

fn filled_vec<T: Clone>(
    length: usize,
    value: T,
    buffer: &'static str,
) -> Result<Vec<T>, WavefrontError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| WavefrontError::AllocationFailed {
            buffer,
            elements: length,
        })?;
    values.resize(length, value);
    Ok(values)
}

fn checked_bytes(
    elements: usize,
    element_size: usize,
    first_len: usize,
    second_len: usize,
) -> Result<usize, WavefrontError> {
    elements
        .checked_mul(element_size)
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })
}

#[derive(Clone, Copy, Debug, Default)]
struct Statistics {
    wavefronts: usize,
    computed_offsets: usize,
    source_wavefront_lookups: usize,
    source_offset_lookups: usize,
    match_extensions: usize,
    maximum_wavefront_width: usize,
    overflowed: bool,
}

impl Statistics {
    fn add_source_wavefront_lookups(&mut self) {
        let (lookups, overflowed) = self.source_wavefront_lookups.overflowing_add(3);
        self.source_wavefront_lookups = lookups;
        self.overflowed |= overflowed;
    }

    fn add_source_offset_lookups(&mut self, width: usize) {
        let (lookups, multiplication_overflowed) = width.overflowing_mul(5);
        let (total, addition_overflowed) = self.source_offset_lookups.overflowing_add(lookups);
        self.source_offset_lookups = total;
        self.overflowed |= multiplication_overflowed | addition_overflowed;
    }

    fn add_wavefront(&mut self, width: usize, matches: usize) {
        let (wavefronts, wavefronts_overflowed) = self.wavefronts.overflowing_add(1);
        let (computed_offsets, offsets_overflowed) = self.computed_offsets.overflowing_add(width);
        let (match_extensions, matches_overflowed) = self.match_extensions.overflowing_add(matches);
        self.wavefronts = wavefronts;
        self.computed_offsets = computed_offsets;
        self.match_extensions = match_extensions;
        self.maximum_wavefront_width = self.maximum_wavefront_width.max(width);
        self.overflowed |= wavefronts_overflowed | offsets_overflowed | matches_overflowed;
    }

    fn summary(
        self,
        penalty: Score,
        modeled_peak_bytes: usize,
        first_len: usize,
        second_len: usize,
    ) -> Result<WavefrontScore, WavefrontError> {
        if self.overflowed {
            return Err(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            });
        }
        Ok(WavefrontScore {
            penalty,
            wavefronts: self.wavefronts,
            computed_offsets: self.computed_offsets,
            source_wavefront_lookups: self.source_wavefront_lookups,
            source_offset_lookups: self.source_offset_lookups,
            match_extensions: self.match_extensions,
            maximum_wavefront_width: self.maximum_wavefront_width,
            modeled_peak_bytes,
        })
    }
}

fn increment_offset(
    offset: isize,
    first_len: usize,
    second_len: usize,
) -> Result<isize, WavefrontError> {
    if offset == UNREACHABLE {
        return Ok(UNREACHABLE);
    }
    offset
        .checked_add(1)
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })
}

#[allow(clippy::cast_sign_loss)]
const fn wrapping_usize(value: isize) -> usize {
    value as usize
}

#[inline]
fn insertion_candidate(offset: isize, dimensions: Dimensions) -> isize {
    if wrapping_usize(offset) < wrapping_usize(dimensions.second) {
        offset + 1
    } else {
        UNREACHABLE
    }
}

#[inline]
fn deletion_candidate(offset: isize, diagonal: isize, dimensions: Dimensions) -> isize {
    if wrapping_usize(offset.wrapping_sub(diagonal)) <= wrapping_usize(dimensions.first) {
        offset
    } else {
        UNREACHABLE
    }
}

#[inline]
fn substitution_candidate(offset: isize, diagonal: isize, dimensions: Dimensions) -> isize {
    if wrapping_usize(offset) >= wrapping_usize(dimensions.second) {
        return UNREACHABLE;
    }
    if wrapping_usize(offset.wrapping_sub(diagonal)) >= wrapping_usize(dimensions.first) {
        UNREACHABLE
    } else {
        offset + 1
    }
}

#[inline]
fn select_insertion_candidate(
    opening: isize,
    extension: isize,
    dimensions: Dimensions,
) -> (isize, bool) {
    let (selected, extends) = select_gap(opening, extension);
    let candidate = insertion_candidate(selected, dimensions);
    if candidate != UNREACHABLE {
        return (candidate, extends);
    }
    let fallback = insertion_candidate(if extends { opening } else { extension }, dimensions);
    if fallback == UNREACHABLE {
        (UNREACHABLE, false)
    } else {
        (fallback, !extends)
    }
}

#[inline]
fn select_deletion_candidate(
    opening: isize,
    extension: isize,
    diagonal: isize,
    dimensions: Dimensions,
) -> (isize, bool) {
    let (selected, extends) = select_gap(opening, extension);
    let candidate = deletion_candidate(selected, diagonal, dimensions);
    if candidate != UNREACHABLE {
        return (candidate, extends);
    }
    let fallback = deletion_candidate(
        if extends { opening } else { extension },
        diagonal,
        dimensions,
    );
    if fallback == UNREACHABLE {
        (UNREACHABLE, false)
    } else {
        (fallback, !extends)
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn extend_equal_legal(
    first: &[u8],
    second: &[u8],
    diagonal: isize,
    offset: isize,
    dimensions: Dimensions,
) -> (isize, usize) {
    debug_assert!(legal_coordinate(offset, diagonal, dimensions).is_some());
    let vertical = usize::try_from(offset - diagonal)
        .expect("a validated WFA vertical coordinate is nonnegative");
    let horizontal =
        usize::try_from(offset).expect("a validated WFA horizontal coordinate is nonnegative");
    let (horizontal, matches) = extend_equal_coordinates(first, second, vertical, horizontal);
    (
        isize::try_from(horizontal).expect("validated WFA dimensions fit signed coordinates"),
        matches,
    )
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum MainTrace {
    Start = 0,
    Substitution = 1,
    Deletion = 2,
    Insertion = 3,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(transparent)]
struct AffineTrace(u8);

impl AffineTrace {
    const MAIN_MASK: u8 = 0b0000_0011;
    const INSERTION_EXTENDS: u8 = 0b0000_0100;
    const DELETION_EXTENDS: u8 = 0b0000_1000;

    const fn with_main(mut self, trace: MainTrace) -> Self {
        self.0 = (self.0 & !Self::MAIN_MASK) | trace as u8;
        self
    }

    const fn with_insertion_extension(mut self, extends: bool) -> Self {
        if extends {
            self.0 |= Self::INSERTION_EXTENDS;
        }
        self
    }

    const fn with_deletion_extension(mut self, extends: bool) -> Self {
        if extends {
            self.0 |= Self::DELETION_EXTENDS;
        }
        self
    }

    const fn main(self) -> MainTrace {
        match self.0 & Self::MAIN_MASK {
            0 => MainTrace::Start,
            1 => MainTrace::Substitution,
            2 => MainTrace::Deletion,
            3 => MainTrace::Insertion,
            _ => unreachable!(),
        }
    }

    const fn insertion_extends(self) -> bool {
        self.0 & Self::INSERTION_EXTENDS != 0
    }

    const fn deletion_extends(self) -> bool {
        self.0 & Self::DELETION_EXTENDS != 0
    }
}

#[derive(Debug)]
struct AffineWavefront {
    penalty: Score,
    lower: isize,
    width: usize,
    offsets: Vec<isize>,
    traces: Vec<AffineTrace>,
    main_active: Option<(isize, isize)>,
    insertion_active: Option<(isize, isize)>,
    deletion_active: Option<(isize, isize)>,
}

impl AffineWavefront {
    fn new<const TRACE: bool>(
        penalty: Score,
        lower: isize,
        length: usize,
        first_len: usize,
        second_len: usize,
    ) -> Result<Self, WavefrontError> {
        let state_length = length
            .checked_mul(3)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        Ok(Self {
            penalty,
            lower,
            width: length,
            offsets: filled_vec(state_length, UNREACHABLE, "affine WFA state offsets")?,
            traces: if TRACE {
                filled_vec(length, AffineTrace::default(), "affine WFA traces")?
            } else {
                Vec::new()
            },
            main_active: None,
            insertion_active: None,
            deletion_active: None,
        })
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn reset<const TRACE: bool>(
        &mut self,
        penalty: Score,
        lower: isize,
        length: usize,
        first_len: usize,
        second_len: usize,
    ) -> Result<(), WavefrontError> {
        let state_length = length
            .checked_mul(3)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        self.penalty = penalty;
        self.lower = lower;
        self.width = length;
        if self.offsets.len() >= state_length {
            self.offsets.truncate(state_length);
        } else {
            if self.offsets.capacity() < state_length {
                self.offsets
                    .try_reserve(state_length - self.offsets.len())
                    .map_err(|_| WavefrontError::AllocationFailed {
                        buffer: "affine WFA state offsets",
                        elements: state_length,
                    })?;
            }
            self.offsets.resize(state_length, UNREACHABLE);
        }
        if TRACE {
            if self.traces.len() >= length {
                self.traces.truncate(length);
            } else {
                if self.traces.capacity() < length {
                    self.traces
                        .try_reserve(length - self.traces.len())
                        .map_err(|_| WavefrontError::AllocationFailed {
                            buffer: "affine WFA traces",
                            elements: length,
                        })?;
                }
                self.traces.resize(length, AffineTrace::default());
            }
        } else {
            self.traces.clear();
        }
        Ok(())
    }

    fn index(&self, diagonal: isize) -> Option<usize> {
        let index = wrapping_usize(diagonal.wrapping_sub(self.lower));
        (index < self.width).then_some(index)
    }

    fn main(&self, diagonal: isize) -> isize {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    fn insertion(&self, diagonal: isize) -> isize {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(self.width + index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    fn deletion(&self, diagonal: isize) -> isize {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(2 * self.width + index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    #[inline]
    fn main_with_trace(&self, diagonal: isize) -> Option<(isize, AffineTrace)> {
        let index = self.index(diagonal)?;
        Some((*self.offsets.get(index)?, *self.traces.get(index)?))
    }

    #[inline]
    fn insertion_with_trace(&self, diagonal: isize) -> Option<(isize, AffineTrace)> {
        let index = self.index(diagonal)?;
        Some((
            *self.offsets.get(self.width + index)?,
            *self.traces.get(index)?,
        ))
    }

    #[inline]
    fn deletion_with_trace(&self, diagonal: isize) -> Option<(isize, AffineTrace)> {
        let index = self.index(diagonal)?;
        Some((
            *self.offsets.get(2 * self.width + index)?,
            *self.traces.get(index)?,
        ))
    }

    fn reaches_target(&self, dimensions: Dimensions) -> bool {
        self.main(dimensions.target_diagonal) >= dimensions.second
    }
}

struct AffineRun {
    summary: WavefrontScore,
    history: Vec<AffineWavefront>,
    dimensions: Dimensions,
}

#[derive(Clone, Copy)]
struct SourceWavefronts<'a> {
    mismatch: Option<&'a AffineWavefront>,
    opening: Option<&'a AffineWavefront>,
    extension: Option<&'a AffineWavefront>,
}

#[derive(Debug)]
struct SmallPenaltyHistoryIndex {
    indices: [usize; SMALL_HISTORY_SLOTS],
}

impl Default for SmallPenaltyHistoryIndex {
    fn default() -> Self {
        Self::new()
    }
}

impl SmallPenaltyHistoryIndex {
    const fn new() -> Self {
        Self {
            indices: [usize::MAX; SMALL_HISTORY_SLOTS],
        }
    }

    // Truncation cannot change the six low bits used by this 64-slot ring.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::inline_always
    )]
    #[inline(always)]
    fn insert(&mut self, penalty: Score, index: usize) {
        debug_assert!(penalty >= 0);
        self.indices[penalty as usize & (SMALL_HISTORY_SLOTS - 1)] = index;
    }

    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::inline_always
    )]
    #[inline(always)]
    fn get(&self, penalty: Score) -> Option<usize> {
        if penalty < 0 {
            None
        } else {
            Some(self.indices[penalty as usize & (SMALL_HISTORY_SLOTS - 1)])
        }
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn indexed_source_wavefronts<'a>(
    history: &'a [AffineWavefront],
    index: &SmallPenaltyHistoryIndex,
    penalty: Score,
    penalties: WavefrontAffinePenalties,
) -> SourceWavefronts<'a> {
    let find = |source_penalty| {
        index
            .get(source_penalty)
            .and_then(|source_index| history.get(source_index))
            .filter(|wavefront| wavefront.penalty == source_penalty)
    };
    SourceWavefronts {
        mismatch: find(penalty - penalties.mismatch),
        opening: find(penalty - penalties.gap_open_transition),
        extension: find(penalty - penalties.gap_extension),
    }
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn indexed_rolling_source_wavefronts<'a>(
    history: &'a VecDeque<AffineWavefront>,
    history_base: usize,
    index: &SmallPenaltyHistoryIndex,
    penalty: Score,
    penalties: WavefrontAffinePenalties,
) -> SourceWavefronts<'a> {
    let find = |source_penalty| {
        index
            .get(source_penalty)
            .and_then(|source_index| source_index.checked_sub(history_base))
            .and_then(|source_index| history.get(source_index))
            .filter(|wavefront| wavefront.penalty == source_penalty)
    };
    SourceWavefronts {
        mismatch: find(penalty - penalties.mismatch),
        opening: find(penalty - penalties.gap_open_transition),
        extension: find(penalty - penalties.gap_extension),
    }
}

#[inline]
fn find_affine_wavefront_index(
    history: &[AffineWavefront],
    mut index: usize,
    penalty: Score,
) -> Option<usize> {
    loop {
        let wavefront = history.get(index)?;
        if wavefront.penalty == penalty {
            return Some(index);
        }
        if wavefront.penalty < penalty || index == 0 {
            return None;
        }
        index -= 1;
    }
}

fn source_wavefronts(
    history: &[AffineWavefront],
    penalty: Score,
    penalties: WavefrontAffinePenalties,
) -> SourceWavefronts<'_> {
    collect_source_wavefronts(history.iter().rev(), penalty, penalties)
}

fn collect_source_wavefronts<'a>(
    history: impl Iterator<Item = &'a AffineWavefront>,
    penalty: Score,
    penalties: WavefrontAffinePenalties,
) -> SourceWavefronts<'a> {
    let mismatch_penalty = penalty - penalties.mismatch;
    let opening_penalty = penalty - penalties.gap_open_transition;
    let extension_penalty = penalty - penalties.gap_extension;
    let mut sources = SourceWavefronts {
        mismatch: None,
        opening: None,
        extension: None,
    };
    let mut unresolved = 0_u8;
    let mut minimum_penalty = Score::MAX;
    if mismatch_penalty >= 0 {
        unresolved |= 0b001;
        minimum_penalty = minimum_penalty.min(mismatch_penalty);
    }
    if opening_penalty >= 0 {
        unresolved |= 0b010;
        minimum_penalty = minimum_penalty.min(opening_penalty);
    }
    if extension_penalty >= 0 {
        unresolved |= 0b100;
        minimum_penalty = minimum_penalty.min(extension_penalty);
    }
    if unresolved == 0 {
        return sources;
    }
    for wavefront in history {
        if wavefront.penalty < minimum_penalty {
            break;
        }
        if mismatch_penalty == wavefront.penalty {
            sources.mismatch = Some(wavefront);
            unresolved &= !0b001;
        }
        if opening_penalty == wavefront.penalty {
            sources.opening = Some(wavefront);
            unresolved &= !0b010;
        }
        if extension_penalty == wavefront.penalty {
            sources.extension = Some(wavefront);
            unresolved &= !0b100;
        }
        if unresolved == 0 || wavefront.penalty == minimum_penalty {
            break;
        }
    }
    sources
}

fn rolling_source_wavefronts(
    history: &VecDeque<AffineWavefront>,
    penalty: Score,
    penalties: WavefrontAffinePenalties,
) -> SourceWavefronts<'_> {
    collect_source_wavefronts(history.iter().rev(), penalty, penalties)
}

fn merge_shifted_interval(
    bounds_lower: &mut isize,
    bounds_upper: &mut isize,
    interval: Option<(isize, isize)>,
    shift: isize,
    dimensions: Dimensions,
) {
    let Some((lower, upper)) = interval else {
        return;
    };
    let lower = lower.saturating_add(shift);
    let upper = upper.saturating_add(shift);
    let legal_lower = -dimensions.first;
    let legal_upper = dimensions.second;
    let clipped_lower = lower.max(legal_lower);
    let clipped_upper = upper.min(legal_upper);
    if clipped_lower > clipped_upper {
        return;
    }
    *bounds_lower = (*bounds_lower).min(clipped_lower);
    *bounds_upper = (*bounds_upper).max(clipped_upper);
}

#[derive(Clone, Copy)]
struct AffineCursor<'a> {
    main: &'a [isize],
    insertion: &'a [isize],
    deletion: &'a [isize],
    index: isize,
}

impl<'a> AffineCursor<'a> {
    fn new(wavefront: Option<&'a AffineWavefront>, diagonal: isize) -> Self {
        if let Some(source) = wavefront {
            let (main, gaps) = source.offsets.split_at(source.width);
            let (insertion, deletion) = gaps.split_at(source.width);
            Self {
                main,
                insertion,
                deletion,
                index: diagonal - source.lower,
            }
        } else {
            Self {
                main: &[],
                insertion: &[],
                deletion: &[],
                index: 0,
            }
        }
    }

    // This cursor lookup is the innermost affine-wavefront hot path.
    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn get(self, values: &[isize], shift: isize) -> isize {
        let source_index = wrapping_usize(self.index.wrapping_add(shift));
        values.get(source_index).copied().unwrap_or(UNREACHABLE)
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn main(self, shift: isize) -> isize {
        self.get(self.main, shift)
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn insertion(self, shift: isize) -> isize {
        self.get(self.insertion, shift)
    }

    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn deletion(self, shift: isize) -> isize {
        self.get(self.deletion, shift)
    }

    // Advancing is paired with every hot-path cursor lookup.
    #[allow(clippy::inline_always)]
    #[inline(always)]
    fn advance(&mut self) {
        self.index += 1;
    }
}

#[inline]
fn select_gap(opening: isize, extension: isize) -> (isize, bool) {
    if opening >= extension {
        (opening, false)
    } else {
        (extension, true)
    }
}

#[inline]
fn select_affine_main(
    substitution: isize,
    deletion: isize,
    insertion: isize,
) -> (isize, MainTrace) {
    let mut best = substitution;
    let mut trace = MainTrace::Substitution;
    if deletion > best {
        best = deletion;
        trace = MainTrace::Deletion;
    }
    if insertion > best {
        best = insertion;
        trace = MainTrace::Insertion;
    }
    (best, trace)
}

#[derive(Clone, Copy)]
struct AffineWavefrontLimits {
    dimensions: Dimensions,
    upper_bound: Score,
    prune_main: bool,
    gap_opening: Score,
    gap_extension: Score,
}

#[allow(clippy::inline_always, clippy::too_many_lines)]
#[inline(always)]
fn compute_affine_wavefront<const TRACE: bool>(
    first: &[u8],
    second: &[u8],
    penalty: Score,
    limits: AffineWavefrontLimits,
    sources: SourceWavefronts<'_>,
    recycled: Option<AffineWavefront>,
) -> Result<Option<(AffineWavefront, usize)>, WavefrontError> {
    let AffineWavefrontLimits {
        dimensions,
        upper_bound,
        prune_main,
        gap_opening,
        gap_extension,
    } = limits;
    let mut lower = isize::MAX;
    let mut upper = isize::MIN;
    if let Some(mismatch) = sources.mismatch {
        merge_shifted_interval(&mut lower, &mut upper, mismatch.main_active, 0, dimensions);
    }
    if let Some(opening) = sources.opening {
        merge_shifted_interval(&mut lower, &mut upper, opening.main_active, 1, dimensions);
        merge_shifted_interval(&mut lower, &mut upper, opening.main_active, -1, dimensions);
    }
    if let Some(extension) = sources.extension {
        merge_shifted_interval(
            &mut lower,
            &mut upper,
            extension.insertion_active,
            1,
            dimensions,
        );
        merge_shifted_interval(
            &mut lower,
            &mut upper,
            extension.deletion_active,
            -1,
            dimensions,
        );
    }
    constrain_diagonal_bounds(
        &mut lower,
        &mut upper,
        penalty,
        upper_bound,
        gap_extension,
        dimensions,
    );
    if lower > upper {
        return Ok(None);
    }
    let main_completion = if !prune_main || upper_bound == Score::MAX {
        None
    } else {
        let remaining = upper_bound - penalty;
        let shift = if remaining >= gap_opening {
            isize::try_from((remaining - gap_opening) / gap_extension).unwrap_or(isize::MAX)
        } else {
            0
        };
        Some((
            dimensions.target_diagonal.saturating_sub(shift),
            dimensions.target_diagonal.saturating_add(shift),
        ))
    };
    let length = checked_width(lower, upper, first.len(), second.len())?;
    let mut current = if let Some(mut recycled) = recycled {
        recycled.reset::<TRACE>(penalty, lower, length, first.len(), second.len())?;
        recycled
    } else {
        AffineWavefront::new::<TRACE>(penalty, lower, length, first.len(), second.len())?
    };
    let mut matches = 0_usize;
    let mut matches_overflowed = false;
    let mut mismatch = AffineCursor::new(sources.mismatch, lower);
    let mut opening = AffineCursor::new(sources.opening, lower);
    let mut extension = AffineCursor::new(sources.extension, lower);
    let mut diagonal = lower;
    let (main_offsets, gap_offsets) = current.offsets.split_at_mut(length);
    let (insertion_offsets, deletion_offsets) = gap_offsets.split_at_mut(length);
    let traces = &mut current.traces;
    let mut main_lower = isize::MAX;
    let mut main_upper = isize::MIN;
    let mut insertion_lower = isize::MAX;
    let mut insertion_upper = isize::MIN;
    let mut deletion_lower = isize::MAX;
    let mut deletion_upper = isize::MIN;

    for index in 0..length {
        let (insertion, insertion_extends) =
            select_insertion_candidate(opening.main(-1), extension.insertion(-1), dimensions);
        let (deletion, deletion_extends) =
            select_deletion_candidate(opening.main(1), extension.deletion(1), diagonal, dimensions);
        let substitution = substitution_candidate(mismatch.main(0), diagonal, dimensions);
        let (mut candidate, main_trace) = select_affine_main(substitution, deletion, insertion);
        if main_completion
            .is_some_and(|(main_lower, main_upper)| diagonal < main_lower || diagonal > main_upper)
        {
            candidate = UNREACHABLE;
        }

        insertion_offsets[index] = insertion;
        deletion_offsets[index] = deletion;
        if insertion != UNREACHABLE {
            insertion_lower = insertion_lower.min(diagonal);
            insertion_upper = diagonal;
        }
        if deletion != UNREACHABLE {
            deletion_lower = deletion_lower.min(diagonal);
            deletion_upper = diagonal;
        }
        if candidate == UNREACHABLE {
            main_offsets[index] = UNREACHABLE;
            if TRACE {
                traces[index] = AffineTrace::default();
            }
        } else {
            let (offset, match_extension) =
                extend_equal_legal(first, second, diagonal, candidate, dimensions);
            main_offsets[index] = offset;
            main_lower = main_lower.min(diagonal);
            main_upper = diagonal;
            if TRACE {
                traces[index] = AffineTrace::default()
                    .with_main(main_trace)
                    .with_insertion_extension(insertion_extends)
                    .with_deletion_extension(deletion_extends);
            }
            let (next_matches, overflowed) = matches.overflowing_add(match_extension);
            matches = next_matches;
            matches_overflowed |= overflowed;
        }
        if index + 1 < length {
            diagonal += 1;
            mismatch.advance();
            opening.advance();
            extension.advance();
        }
    }

    current.main_active = (main_lower <= main_upper).then_some((main_lower, main_upper));
    current.insertion_active =
        (insertion_lower <= insertion_upper).then_some((insertion_lower, insertion_upper));
    current.deletion_active =
        (deletion_lower <= deletion_upper).then_some((deletion_lower, deletion_upper));

    if matches_overflowed {
        Err(WavefrontError::DimensionOverflow {
            first_len: first.len(),
            second_len: second.len(),
        })
    } else if current.main_active.is_none()
        && (TRACE || (current.insertion_active.is_none() && current.deletion_active.is_none()))
    {
        Ok(None)
    } else {
        Ok(Some((current, matches)))
    }
}

trait AffinePenaltySchedule {
    const INDEXED_HISTORY: bool;

    fn schedule(&mut self, penalty: Score, increment: Score);
    fn pop(&mut self) -> Option<Score>;
    fn overflowed(&self) -> bool;
}

struct SmallPenaltyScheduler {
    slots: u64,
    cursor: Score,
    overflowed: bool,
}

impl SmallPenaltyScheduler {
    const fn new() -> Self {
        Self {
            slots: 0,
            cursor: 0,
            overflowed: false,
        }
    }
}

impl AffinePenaltySchedule for SmallPenaltyScheduler {
    const INDEXED_HISTORY: bool = true;

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn schedule(&mut self, penalty: Score, increment: Score) {
        let Some(next) = penalty.checked_add(increment) else {
            self.overflowed = true;
            return;
        };
        debug_assert!(next > self.cursor);
        debug_assert!(next - self.cursor <= SMALL_SCHEDULER_MAX_INCREMENT);
        self.slots |= 1_u64 << (next as usize & (SMALL_HISTORY_SLOTS - 1));
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn pop(&mut self) -> Option<Score> {
        if self.slots == 0 {
            return None;
        }
        let cursor_slot = self.cursor as usize & (SMALL_HISTORY_SLOTS - 1);
        let shift = u32::try_from((cursor_slot + 1) & (SMALL_HISTORY_SLOTS - 1))
            .expect("a small scheduler rotation fits u32");
        let delta = Score::from(self.slots.rotate_right(shift).trailing_zeros()) + 1;
        let Some(next) = self.cursor.checked_add(delta) else {
            self.overflowed = true;
            return None;
        };
        debug_assert!(delta <= SMALL_SCHEDULER_MAX_INCREMENT);
        self.cursor = next;
        self.slots &= !(1_u64 << (next as usize & (SMALL_HISTORY_SLOTS - 1)));
        Some(next)
    }

    fn overflowed(&self) -> bool {
        self.overflowed
    }
}

struct DensePenaltyScheduler {
    slots: Vec<u8>,
    cursor: Score,
    pending: usize,
    overflowed: bool,
}

impl DensePenaltyScheduler {
    fn new(maximum_increment: Score) -> Result<Self, WavefrontError> {
        let length = usize::try_from(maximum_increment)
            .expect("the bounded affine increment fits usize")
            + 1;
        Ok(Self {
            slots: filled_vec(length, 0_u8, "affine WFA penalty scheduler")?,
            cursor: 0,
            pending: 0,
            overflowed: false,
        })
    }
}

impl AffinePenaltySchedule for DensePenaltyScheduler {
    const INDEXED_HISTORY: bool = false;

    fn schedule(&mut self, penalty: Score, increment: Score) {
        let Some(next) = penalty.checked_add(increment) else {
            self.overflowed = true;
            return;
        };
        debug_assert!(next > self.cursor);
        let index = usize::try_from(next).expect("a nonnegative bounded affine penalty fits usize")
            % self.slots.len();
        if self.slots[index] == 0 {
            self.slots[index] = 1;
            self.pending += 1;
        }
    }

    fn pop(&mut self) -> Option<Score> {
        while self.pending > 0 {
            let Some(next) = self.cursor.checked_add(1) else {
                self.overflowed = true;
                return None;
            };
            self.cursor = next;
            let index = usize::try_from(next)
                .expect("a nonnegative bounded affine penalty fits usize")
                % self.slots.len();
            if self.slots[index] != 0 {
                self.slots[index] = 0;
                self.pending -= 1;
                return Some(next);
            }
        }
        None
    }

    fn overflowed(&self) -> bool {
        self.overflowed
    }
}

struct SparsePenaltyScheduler {
    heap: BinaryHeap<Reverse<Score>>,
    cursor: Score,
    overflowed: bool,
}

impl SparsePenaltyScheduler {
    const fn new() -> Self {
        Self {
            heap: BinaryHeap::new(),
            cursor: 0,
            overflowed: false,
        }
    }
}

impl AffinePenaltySchedule for SparsePenaltyScheduler {
    const INDEXED_HISTORY: bool = false;

    fn schedule(&mut self, penalty: Score, increment: Score) {
        let Some(next) = penalty.checked_add(increment) else {
            self.overflowed = true;
            return;
        };
        self.heap.push(Reverse(next));
    }

    fn pop(&mut self) -> Option<Score> {
        while let Some(Reverse(next)) = self.heap.pop() {
            if next > self.cursor {
                self.cursor = next;
                return Some(next);
            }
        }
        None
    }

    fn overflowed(&self) -> bool {
        self.overflowed
    }
}

enum PenaltyScheduler {
    Small(SmallPenaltyScheduler),
    Dense(DensePenaltyScheduler),
    Sparse(SparsePenaltyScheduler),
}

impl PenaltyScheduler {
    fn new(maximum_increment: Score) -> Result<Self, WavefrontError> {
        if maximum_increment <= SMALL_SCHEDULER_MAX_INCREMENT {
            Ok(Self::Small(SmallPenaltyScheduler::new()))
        } else if maximum_increment <= DENSE_SCHEDULER_MAX_INCREMENT {
            Ok(Self::Dense(DensePenaltyScheduler::new(maximum_increment)?))
        } else {
            Ok(Self::Sparse(SparsePenaltyScheduler::new()))
        }
    }
}

fn initial_affine_wavefront<const TRACE: bool>(
    first: &[u8],
    second: &[u8],
    dimensions: Dimensions,
    recycled: Option<AffineWavefront>,
) -> Result<(AffineWavefront, usize), WavefrontError> {
    let mut wavefront = if let Some(mut recycled) = recycled {
        recycled.reset::<TRACE>(0, 0, 1, first.len(), second.len())?;
        recycled.main_active = None;
        recycled.insertion_active = None;
        recycled.deletion_active = None;
        recycled
    } else {
        AffineWavefront::new::<TRACE>(0, 0, 1, first.len(), second.len())?
    };
    let (offset, matches) = extend_equal(first, second, 0, 0, dimensions).ok_or(
        WavefrontError::TracebackInvariant {
            penalty: 0,
            diagonal: 0,
            offset: 0,
        },
    )?;
    wavefront.offsets[0] = offset;
    wavefront.offsets[1] = UNREACHABLE;
    wavefront.offsets[2] = UNREACHABLE;
    wavefront.main_active = Some((0, 0));
    if TRACE {
        wavefront.traces[0] = AffineTrace::default().with_main(MainTrace::Start);
    }
    Ok((wavefront, matches))
}

struct AffineScoreRun {
    penalty: Score,
    statistics: Statistics,
    peak_state_offsets: usize,
}

fn discard_expired_affine_wavefronts<const METRICS: bool>(
    history: &mut VecDeque<AffineWavefront>,
    current_penalty: Score,
    maximum_increment: Score,
    live_state_offsets: &mut usize,
    recycled: &mut Vec<AffineWavefront>,
) -> usize {
    let mut discarded = 0_usize;
    while history.front().is_some_and(|wavefront| {
        wavefront
            .penalty
            .checked_add(maximum_increment)
            .is_some_and(|last_use| last_use <= current_penalty)
    }) {
        if let Some(expired) = history.pop_front() {
            if METRICS {
                *live_state_offsets -= expired.offsets.len();
            }
            recycled.push(expired);
            discarded += 1;
        }
    }
    discarded
}

#[allow(clippy::too_many_lines)]
fn run_affine_score(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
) -> Result<WavefrontScore, WavefrontError> {
    let run = run_affine_score_core::<true>(first, second, penalties, workspace)?;
    let modeled_peak_bytes = checked_bytes(
        run.peak_state_offsets,
        size_of::<isize>(),
        first.len(),
        second.len(),
    )?;
    run.statistics
        .summary(run.penalty, modeled_peak_bytes, first.len(), second.len())
}

fn run_affine_penalty(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
) -> Result<Score, WavefrontError> {
    if first.len() != second.len() && workspace.pending_score_bound.is_none() {
        if let Some((upper_bound, prune_main)) =
            unequal_compact_upper_bound(first, second, penalties)
        {
            if let Some(penalty) = compact_affine_penalty(
                first,
                second,
                penalties,
                &mut workspace.compact_score,
                upper_bound,
                prune_main,
            )? {
                return Ok(penalty);
            }
        }
        return Ok(run_affine_score_core::<false>(first, second, penalties, workspace)?.penalty);
    }
    let (upper_bound, prune_main) = workspace
        .pending_score_bound
        .unwrap_or_else(|| exact_alignment_upper_bound(first, second, penalties));
    if first.len() == second.len() {
        if upper_bound == 0 {
            return Ok(0);
        }
        if let Some(penalty) = compact_affine_penalty(
            first,
            second,
            penalties,
            &mut workspace.compact_score,
            upper_bound,
            prune_main,
        )? {
            return Ok(penalty);
        }
    }
    Ok(run_affine_score_core::<false>(first, second, penalties, workspace)?.penalty)
}

fn run_affine_penalty_with_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
    upper_bound: Score,
    prune_main: bool,
) -> Result<Score, WavefrontError> {
    debug_assert!(workspace.pending_score_bound.is_none());
    workspace.pending_score_bound = Some((upper_bound, prune_main));
    let result = run_affine_penalty(first, second, penalties, workspace);
    workspace.pending_score_bound = None;
    result
}

#[allow(clippy::too_many_lines)]
fn run_affine_score_core<const METRICS: bool>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
) -> Result<AffineScoreRun, WavefrontError> {
    let (upper_bound, prune_main) = workspace
        .pending_score_bound
        .unwrap_or_else(|| exact_alignment_upper_bound(first, second, penalties));
    let dimensions = dimensions(first.len(), second.len())?;
    let maximum_increment = penalties
        .mismatch
        .max(penalties.gap_open_transition)
        .max(penalties.gap_extension);
    let WavefrontWorkspace {
        rolling_history: history,
        recycled,
        small_history_index: source_index,
        ..
    } = workspace;
    while let Some(expired) = history.pop_front() {
        recycled.push(expired);
    }
    let (initial, initial_matches) =
        initial_affine_wavefront::<false>(first, second, dimensions, recycled.pop())?;
    let initial_reaches_target = initial.reaches_target(dimensions);
    let live_state_offsets = if METRICS { initial.offsets.len() } else { 0 };
    let peak_state_offsets = live_state_offsets;
    if history.capacity() == 0 {
        history
            .try_reserve(1)
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "rolling affine WFA history",
                elements: 1,
            })?;
    }
    history.push_back(initial);
    let mut statistics = Statistics::default();
    if METRICS {
        statistics.add_wavefront(1, initial_matches);
    }
    if initial_reaches_target {
        while let Some(expired) = history.pop_front() {
            recycled.push(expired);
        }
        return Ok(AffineScoreRun {
            penalty: 0,
            statistics,
            peak_state_offsets,
        });
    }

    let scheduler = PenaltyScheduler::new(maximum_increment)?;
    match scheduler {
        PenaltyScheduler::Small(scheduler) => run_affine_score_scheduled::<METRICS, _>(
            first,
            second,
            penalties,
            dimensions,
            history,
            recycled,
            statistics,
            live_state_offsets,
            peak_state_offsets,
            maximum_increment,
            upper_bound,
            prune_main,
            scheduler,
            source_index,
        ),
        PenaltyScheduler::Dense(scheduler) => run_affine_score_scheduled::<METRICS, _>(
            first,
            second,
            penalties,
            dimensions,
            history,
            recycled,
            statistics,
            live_state_offsets,
            peak_state_offsets,
            maximum_increment,
            upper_bound,
            prune_main,
            scheduler,
            source_index,
        ),
        PenaltyScheduler::Sparse(scheduler) => run_affine_score_scheduled::<METRICS, _>(
            first,
            second,
            penalties,
            dimensions,
            history,
            recycled,
            statistics,
            live_state_offsets,
            peak_state_offsets,
            maximum_increment,
            upper_bound,
            prune_main,
            scheduler,
            source_index,
        ),
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run_affine_score_scheduled<const METRICS: bool, S: AffinePenaltySchedule>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    dimensions: Dimensions,
    history: &mut VecDeque<AffineWavefront>,
    recycled: &mut Vec<AffineWavefront>,
    mut statistics: Statistics,
    mut live_state_offsets: usize,
    mut peak_state_offsets: usize,
    maximum_increment: Score,
    upper_bound: Score,
    prune_main: bool,
    mut scheduler: S,
    source_index: &mut SmallPenaltyHistoryIndex,
) -> Result<AffineScoreRun, WavefrontError> {
    let limits = AffineWavefrontLimits {
        dimensions,
        upper_bound,
        prune_main: !METRICS && prune_main,
        gap_opening: penalties.gap_opening,
        gap_extension: penalties.gap_extension,
    };
    let mut history_base = 0_usize;
    if S::INDEXED_HISTORY {
        source_index.insert(0, 0);
    }
    scheduler.schedule(0, penalties.mismatch);
    scheduler.schedule(0, penalties.gap_open_transition);

    while let Some(penalty) = scheduler.pop() {
        if METRICS {
            statistics.add_source_wavefront_lookups();
        }
        let sources = if S::INDEXED_HISTORY {
            indexed_rolling_source_wavefronts(
                history,
                history_base,
                source_index,
                penalty,
                penalties,
            )
        } else {
            rolling_source_wavefronts(history, penalty, penalties)
        };
        let Some((current, matches)) = compute_affine_wavefront::<false>(
            first,
            second,
            penalty,
            limits,
            sources,
            recycled.pop(),
        )?
        else {
            history_base += discard_expired_affine_wavefronts::<METRICS>(
                history,
                penalty,
                maximum_increment,
                &mut live_state_offsets,
                recycled,
            );
            continue;
        };
        if METRICS {
            statistics.add_wavefront(current.width, matches);
            statistics.add_source_offset_lookups(current.width);
        }
        let reaches_target = current.reaches_target(dimensions);
        let has_gap = current.insertion_active.is_some() || current.deletion_active.is_some();
        if METRICS {
            live_state_offsets = live_state_offsets
                .checked_add(current.offsets.len())
                .ok_or(WavefrontError::DimensionOverflow {
                    first_len: first.len(),
                    second_len: second.len(),
                })?;
            peak_state_offsets = peak_state_offsets.max(live_state_offsets);
        }
        if history.len() == history.capacity() {
            history
                .try_reserve(1)
                .map_err(|_| WavefrontError::AllocationFailed {
                    buffer: "rolling affine WFA history",
                    elements: history.len().saturating_add(1),
                })?;
        }
        history.push_back(current);
        if S::INDEXED_HISTORY {
            source_index.insert(penalty, history_base + history.len() - 1);
        }
        if reaches_target {
            while let Some(expired) = history.pop_front() {
                recycled.push(expired);
            }
            return Ok(AffineScoreRun {
                penalty,
                statistics,
                peak_state_offsets,
            });
        }
        scheduler.schedule(penalty, penalties.mismatch);
        scheduler.schedule(penalty, penalties.gap_open_transition);
        if has_gap {
            scheduler.schedule(penalty, penalties.gap_extension);
        }
        history_base += discard_expired_affine_wavefronts::<METRICS>(
            history,
            penalty,
            maximum_increment,
            &mut live_state_offsets,
            recycled,
        );
    }

    let last_penalty = history.back().map_or(0, |wavefront| wavefront.penalty);
    let last_offset = history.back().map_or(UNREACHABLE, |wavefront| {
        wavefront.main(dimensions.target_diagonal)
    });
    while let Some(expired) = history.pop_front() {
        recycled.push(expired);
    }
    if scheduler.overflowed() {
        Err(WavefrontError::PenaltyOverflow {
            penalty: last_penalty,
        })
    } else {
        Err(WavefrontError::TracebackInvariant {
            penalty: last_penalty,
            diagonal: dimensions.target_diagonal,
            offset: last_offset,
        })
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run_affine_with_bound<const TRACE: bool>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut WavefrontWorkspace,
    upper_bound: Score,
) -> Result<AffineRun, WavefrontError> {
    let dimensions = dimensions(first.len(), second.len())?;
    let maximum_increment = penalties
        .mismatch
        .max(penalties.gap_open_transition)
        .max(penalties.gap_extension);
    let mut history = core::mem::take(&mut workspace.trace_history);
    workspace.recycled.append(&mut history);
    if history.capacity() == 0 {
        history
            .try_reserve_exact(1)
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "affine WFA history",
                elements: 1,
            })?;
    }
    let (initial, initial_matches) =
        initial_affine_wavefront::<TRACE>(first, second, dimensions, workspace.recycled.pop())?;
    let initial_reaches_target = initial.reaches_target(dimensions);
    history.push(initial);
    let mut statistics = Statistics::default();
    statistics.add_wavefront(1, initial_matches);
    let terminal_penalty = if initial_reaches_target {
        0
    } else {
        let scheduler = PenaltyScheduler::new(maximum_increment)?;
        match scheduler {
            PenaltyScheduler::Small(scheduler) => run_affine_scheduled::<TRACE, _>(
                first,
                second,
                penalties,
                dimensions,
                &mut history,
                &mut workspace.recycled,
                &mut statistics,
                upper_bound,
                scheduler,
                &mut workspace.small_history_index,
            ),
            PenaltyScheduler::Dense(scheduler) => run_affine_scheduled::<TRACE, _>(
                first,
                second,
                penalties,
                dimensions,
                &mut history,
                &mut workspace.recycled,
                &mut statistics,
                upper_bound,
                scheduler,
                &mut workspace.small_history_index,
            ),
            PenaltyScheduler::Sparse(scheduler) => run_affine_scheduled::<TRACE, _>(
                first,
                second,
                penalties,
                dimensions,
                &mut history,
                &mut workspace.recycled,
                &mut statistics,
                upper_bound,
                scheduler,
                &mut workspace.small_history_index,
            ),
        }?
    };

    let state_slots =
        statistics
            .computed_offsets
            .checked_mul(3)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
    let offset_bytes = checked_bytes(state_slots, size_of::<isize>(), first.len(), second.len())?;
    let trace_bytes = if TRACE {
        checked_bytes(
            statistics.computed_offsets,
            size_of::<AffineTrace>(),
            first.len(),
            second.len(),
        )?
    } else {
        0
    };
    let modeled_peak_bytes =
        offset_bytes
            .checked_add(trace_bytes)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
    Ok(AffineRun {
        summary: statistics.summary(
            terminal_penalty,
            modeled_peak_bytes,
            first.len(),
            second.len(),
        )?,
        history,
        dimensions,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_affine_scheduled<const TRACE: bool, S: AffinePenaltySchedule>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    dimensions: Dimensions,
    history: &mut Vec<AffineWavefront>,
    recycled: &mut Vec<AffineWavefront>,
    statistics: &mut Statistics,
    upper_bound: Score,
    mut scheduler: S,
    source_index: &mut SmallPenaltyHistoryIndex,
) -> Result<Score, WavefrontError> {
    let limits = AffineWavefrontLimits {
        dimensions,
        upper_bound,
        prune_main: false,
        gap_opening: penalties.gap_opening,
        gap_extension: penalties.gap_extension,
    };
    if S::INDEXED_HISTORY {
        source_index.insert(0, 0);
    }
    scheduler.schedule(0, penalties.mismatch);
    scheduler.schedule(0, penalties.gap_open_transition);

    while let Some(penalty) = scheduler.pop() {
        statistics.add_source_wavefront_lookups();
        let sources = if S::INDEXED_HISTORY {
            indexed_source_wavefronts(history, source_index, penalty, penalties)
        } else {
            source_wavefronts(history, penalty, penalties)
        };
        let Some((current, matches)) = compute_affine_wavefront::<TRACE>(
            first,
            second,
            penalty,
            limits,
            sources,
            recycled.pop(),
        )?
        else {
            continue;
        };
        statistics.add_wavefront(current.width, matches);
        statistics.add_source_offset_lookups(current.width);
        let reaches_target = current.reaches_target(dimensions);
        let has_gap = current.insertion_active.is_some() || current.deletion_active.is_some();
        if history.len() == history.capacity() {
            history
                .try_reserve(1)
                .map_err(|_| WavefrontError::AllocationFailed {
                    buffer: "affine WFA history",
                    elements: history.len().saturating_add(1),
                })?;
        }
        history.push(current);
        if S::INDEXED_HISTORY {
            source_index.insert(penalty, history.len() - 1);
        }
        if reaches_target {
            return Ok(penalty);
        }
        scheduler.schedule(penalty, penalties.mismatch);
        scheduler.schedule(penalty, penalties.gap_open_transition);
        if has_gap {
            scheduler.schedule(penalty, penalties.gap_extension);
        }
    }
    let last_penalty = history.last().map_or(0, |wavefront| wavefront.penalty);
    if scheduler.overflowed() {
        Err(WavefrontError::PenaltyOverflow {
            penalty: last_penalty,
        })
    } else {
        Err(WavefrontError::TracebackInvariant {
            penalty: last_penalty,
            diagonal: dimensions.target_diagonal,
            offset: history.last().map_or(UNREACHABLE, |wavefront| {
                wavefront.main(dimensions.target_diagonal)
            }),
        })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AffineTracebackState {
    Main,
    Insertion,
    Deletion,
}

fn build_affine_alignment(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    mut run: AffineRun,
    workspace: &mut WavefrontWorkspace,
) -> Result<WavefrontAlignment, WavefrontError> {
    let dimensions = run.dimensions;
    let mut operations = Vec::new();
    operations
        .try_reserve_exact(dimensions.maximum_operations)
        .map_err(|_| WavefrontError::AllocationFailed {
            buffer: "affine WFA alignment operations",
            elements: dimensions.maximum_operations,
        })?;
    workspace.trace_runs.clear();
    traceback_affine(
        first,
        second,
        dimensions,
        penalties,
        run.summary.penalty,
        &run.history,
        &mut workspace.trace_runs,
    )?;
    expand_trace_runs(
        &workspace.trace_runs,
        &mut operations,
        dimensions.maximum_operations,
        first.len(),
        second.len(),
    )?;
    workspace.recycle_trace_history(core::mem::take(&mut run.history));
    let alignment = PairwiseAlignment::try_new(
        run.summary.alignment_score(),
        0..first.len(),
        0..second.len(),
        operations,
    )
    .map_err(WavefrontError::InvalidResult)?;
    let operation_bytes = checked_bytes(
        dimensions.maximum_operations,
        size_of::<AlignmentOperation>(),
        first.len(),
        second.len(),
    )?;
    let trace_run_bytes = checked_bytes(
        workspace.trace_runs.len(),
        size_of::<WavefrontTraceRun>(),
        first.len(),
        second.len(),
    )?;
    run.summary.modeled_peak_bytes = run
        .summary
        .modeled_peak_bytes
        .checked_add(operation_bytes)
        .and_then(|bytes| bytes.checked_add(trace_run_bytes))
        .ok_or(WavefrontError::DimensionOverflow {
            first_len: first.len(),
            second_len: second.len(),
        })?;
    Ok(WavefrontAlignment {
        alignment,
        summary: run.summary,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn traceback_affine(
    first: &[u8],
    second: &[u8],
    dimensions: Dimensions,
    penalties: WavefrontAffinePenalties,
    terminal_penalty: Score,
    history: &[AffineWavefront],
    runs: &mut Vec<WavefrontTraceRun>,
) -> Result<(), WavefrontError> {
    let mut penalty = terminal_penalty;
    let mut diagonal = dimensions.target_diagonal;
    let mut offset = dimensions.second;
    let mut state = AffineTracebackState::Main;
    let mut history_index = history.len().saturating_sub(1);

    loop {
        history_index = find_affine_wavefront_index(history, history_index, penalty).ok_or(
            WavefrontError::TracebackInvariant {
                penalty,
                diagonal,
                offset,
            },
        )?;
        let wavefront = &history[history_index];
        match state {
            AffineTracebackState::Main => {
                let (main, trace) = wavefront.main_with_trace(diagonal).ok_or(
                    WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    },
                )?;
                if main != offset {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    });
                }
                let main_trace = trace.main();
                let mut next_history_index = history_index;
                let source_offset = match main_trace {
                    MainTrace::Start => 0,
                    MainTrace::Substitution => {
                        let source_penalty = penalty.checked_sub(penalties.mismatch).ok_or(
                            WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal,
                                offset,
                            },
                        )?;
                        next_history_index =
                            find_affine_wavefront_index(history, history_index, source_penalty)
                                .ok_or(WavefrontError::TracebackInvariant {
                                    penalty,
                                    diagonal,
                                    offset,
                                })?;
                        increment_offset(
                            history[next_history_index].main(diagonal),
                            first.len(),
                            second.len(),
                        )?
                    }
                    MainTrace::Deletion => wavefront.deletion(diagonal),
                    MainTrace::Insertion => wavefront.insertion(diagonal),
                };
                let match_length = verified_match_length(
                    dimensions,
                    diagonal,
                    source_offset,
                    &mut offset,
                    penalty,
                )?;
                push_trace_run(
                    runs,
                    AlignmentOperation::Match,
                    match_length,
                    first.len(),
                    second.len(),
                )?;
                match main_trace {
                    MainTrace::Start => {
                        if penalty != 0 || diagonal != 0 || offset != 0 {
                            return Err(WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal,
                                offset,
                            });
                        }
                        break;
                    }
                    MainTrace::Substitution => {
                        push_trace_run(
                            runs,
                            AlignmentOperation::Substitution,
                            1,
                            first.len(),
                            second.len(),
                        )?;
                        history_index = next_history_index;
                        offset -= 1;
                        penalty = penalty.checked_sub(penalties.mismatch).ok_or(
                            WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal,
                                offset,
                            },
                        )?;
                    }
                    MainTrace::Deletion => state = AffineTracebackState::Deletion,
                    MainTrace::Insertion => state = AffineTracebackState::Insertion,
                }
            }
            AffineTracebackState::Deletion => {
                let (deletion, trace) = wavefront.deletion_with_trace(diagonal).ok_or(
                    WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    },
                )?;
                if deletion != offset {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    });
                }
                let extends = trace.deletion_extends();
                push_trace_run(
                    runs,
                    AlignmentOperation::Deletion,
                    1,
                    first.len(),
                    second.len(),
                )?;
                diagonal += 1;
                let decrement = if extends {
                    penalties.gap_extension
                } else {
                    penalties.gap_open_transition
                };
                penalty =
                    penalty
                        .checked_sub(decrement)
                        .ok_or(WavefrontError::TracebackInvariant {
                            penalty,
                            diagonal,
                            offset,
                        })?;
                state = if extends {
                    AffineTracebackState::Deletion
                } else {
                    AffineTracebackState::Main
                };
            }
            AffineTracebackState::Insertion => {
                let (insertion, trace) = wavefront.insertion_with_trace(diagonal).ok_or(
                    WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    },
                )?;
                if insertion != offset {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal,
                        offset,
                    });
                }
                let extends = trace.insertion_extends();
                push_trace_run(
                    runs,
                    AlignmentOperation::Insertion,
                    1,
                    first.len(),
                    second.len(),
                )?;
                offset -= 1;
                diagonal -= 1;
                let decrement = if extends {
                    penalties.gap_extension
                } else {
                    penalties.gap_open_transition
                };
                penalty =
                    penalty
                        .checked_sub(decrement)
                        .ok_or(WavefrontError::TracebackInvariant {
                            penalty,
                            diagonal,
                            offset,
                        })?;
                state = if extends {
                    AffineTracebackState::Insertion
                } else {
                    AffineTracebackState::Main
                };
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        MainTrace, UNREACHABLE, WavefrontAffinePenalties, WavefrontWorkspace,
        count_byte_mismatches, exact_alignment_upper_bound, extend_equal_coordinates,
        run_affine_penalty, run_affine_score_core, select_affine_main, select_gap,
        unequal_compact_upper_bound,
    };

    #[test]
    fn predecessor_selectors_implement_the_documented_priority() {
        assert_eq!(select_gap(4, 4), (4, false), "opening precedes extension");
        assert_eq!(select_affine_main(4, 4, 4), (4, MainTrace::Substitution));
        assert_eq!(
            select_affine_main(UNREACHABLE, 4, 4),
            (4, MainTrace::Deletion)
        );
    }

    #[test]
    fn match_extension_finds_word_and_large_chunk_boundaries() {
        let first = vec![b'A'; 4_097];
        let mut second = first.clone();

        for mismatch in [0, 1, 7, 8, 9, 1_023, 1_024, 1_025, 4_095, 4_096] {
            second[mismatch] = b'C';
            assert_eq!(
                extend_equal_coordinates(&first, &second, 0, 0),
                (mismatch, mismatch),
                "first mismatch at {mismatch}"
            );
            second[mismatch] = b'A';
        }

        assert_eq!(
            extend_equal_coordinates(&first, &second, 0, 0),
            (first.len(), first.len())
        );
        assert_eq!(
            extend_equal_coordinates(&first, &second, 13, 17),
            (first.len(), first.len() - 17)
        );
    }

    #[test]
    fn packed_mismatch_count_matches_scalar_bytes_at_every_word_boundary() {
        for length in (0..=257).chain([1_023, 1_024, 1_025, 2_047, 2_048, 2_049, 4_097]) {
            let first: Vec<u8> = (0..length)
                .map(|index| u8::try_from((index * 37 + index / 5) % 256).unwrap())
                .collect();
            let mut second = first.clone();
            for index in (1..length).step_by(7) {
                second[index] = second[index].wrapping_add(1);
            }
            let expected = first
                .iter()
                .zip(&second)
                .filter(|(left, right)| left != right)
                .count();
            assert_eq!(count_byte_mismatches(&first, &second), expected);
        }
    }

    #[test]
    fn main_completion_pruning_is_limited_to_sparse_exact_bounds() {
        let penalties = WavefrontAffinePenalties::new(4, 6, 2).unwrap();
        let short = vec![b'A'; 255];
        assert_eq!(
            exact_alignment_upper_bound(&short, &short, penalties),
            (i64::MAX, false),
            "small matrices avoid a redundant exact-bound scan"
        );
        let first = vec![b'A'; 1_024];
        let mut sparse = first.clone();
        for index in (31..1_024).step_by(97).take(10) {
            sparse[index] = b'C';
        }
        let (sparse_bound, sparse_pruning) =
            exact_alignment_upper_bound(&first, &sparse, penalties);
        assert_eq!(sparse_bound, 40);
        assert!(sparse_pruning);

        let mut dense = first.clone();
        for index in (17..1_024).step_by(37).take(17) {
            dense[index] = b'C';
        }
        let (dense_bound, dense_pruning) = exact_alignment_upper_bound(&first, &dense, penalties);
        assert_eq!(dense_bound, 68);
        assert!(!dense_pruning);

        let mut one_gap = first.clone();
        one_gap.insert(512, b'C');
        let (gap_bound, gap_pruning) = exact_alignment_upper_bound(&first, &one_gap, penalties);
        assert_eq!(gap_bound, 8);
        assert!(gap_pruning);
    }

    #[test]
    fn unequal_compact_gate_preserves_scores_for_broad_and_single_gap_inputs() {
        let penalties = WavefrontAffinePenalties::new(2, 4, 1).unwrap();
        let first: Vec<u8> = (0..384)
            .map(|index| b"ACGT"[(index * 13 + index / 7) % 4])
            .collect();
        let second: Vec<u8> = (0..256)
            .map(|index| b"TGCA"[(index * 11 + index / 5) % 4])
            .collect();
        assert!(unequal_compact_upper_bound(&first, &second, penalties).is_some());
        assert!(unequal_compact_upper_bound(&first[..256], &second, penalties).is_none());
        assert!(
            unequal_compact_upper_bound(&first, &first[..320], penalties).is_none(),
            "a small relative length difference retains the legacy route"
        );

        let mut legacy_workspace = WavefrontWorkspace::new();
        let expected =
            run_affine_score_core::<false>(&first, &second, penalties, &mut legacy_workspace)
                .unwrap()
                .penalty;
        let mut compact_workspace = WavefrontWorkspace::new();
        assert_eq!(
            run_affine_penalty(&first, &second, penalties, &mut compact_workspace).unwrap(),
            expected
        );

        let gap_first = vec![b'A'; 256];
        let mut gap_second = gap_first.clone();
        gap_second.splice(128..128, [b'C'; 128]);
        let mut gap_legacy_workspace = WavefrontWorkspace::new();
        let expected_gap = run_affine_score_core::<false>(
            &gap_first,
            &gap_second,
            penalties,
            &mut gap_legacy_workspace,
        )
        .unwrap()
        .penalty;
        let mut gap_compact_workspace = WavefrontWorkspace::new();
        assert_eq!(
            run_affine_penalty(
                &gap_first,
                &gap_second,
                penalties,
                &mut gap_compact_workspace,
            )
            .unwrap(),
            expected_gap
        );
    }
}
