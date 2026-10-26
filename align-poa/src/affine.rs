use core::mem::size_of;
use core::ops::Range;

use align_core::{AffineGap, Alphabet, Score, Sequence, SubstitutionScorer, Symbol};

use crate::{PoaAlignedPair, PoaError, PoaGraph, PoaNodeId, PoaOperation};

type WideScore = i128;

/// Published non-negative POASTA costs with an exact score-side conversion.
///
/// A gap of positive length `k` costs `opening + k * extension`. The converted
/// [`AffineGap`] therefore scores it as the exact negation of that cost under the
/// workspace length-one-open convention.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoastaAffineCosts {
    mismatch: Score,
    gap_opening: Score,
    gap_extension: Score,
    gap: AffineGap,
}

impl PoastaAffineCosts {
    /// Validates the affine costs required by the POASTA paper.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] when a cost is negative, opening is not strictly
    /// greater than extension, or `opening + extension` is unrepresentable.
    pub fn new(
        mismatch: Score,
        gap_opening: Score,
        gap_extension: Score,
    ) -> Result<Self, PoaError> {
        for (parameter, cost) in [
            ("mismatch", mismatch),
            ("gap-opening", gap_opening),
            ("gap-extension", gap_extension),
        ] {
            if cost < 0 {
                return Err(PoaError::NegativeAffineCost { parameter, cost });
            }
        }
        if gap_opening <= gap_extension {
            return Err(PoaError::AffineGapOpeningNotGreaterThanExtension {
                opening: gap_opening,
                extension: gap_extension,
            });
        }
        let first_symbol_cost =
            gap_opening
                .checked_add(gap_extension)
                .ok_or(PoaError::AffineCostOverflow {
                    opening: gap_opening,
                    extension: gap_extension,
                })?;
        let open_score = first_symbol_cost
            .checked_neg()
            .ok_or(PoaError::AffineCostOverflow {
                opening: gap_opening,
                extension: gap_extension,
            })?;
        let extension_score = gap_extension
            .checked_neg()
            .ok_or(PoaError::AffineCostOverflow {
                opening: gap_opening,
                extension: gap_extension,
            })?;
        let gap = AffineGap::new(open_score, extension_score).map_err(|_| {
            PoaError::AffineCostOverflow {
                opening: gap_opening,
                extension: gap_extension,
            }
        })?;
        Ok(Self {
            mismatch,
            gap_opening,
            gap_extension,
            gap,
        })
    }

    /// Returns the non-negative mismatch cost.
    #[must_use]
    pub const fn mismatch(self) -> Score {
        self.mismatch
    }

    /// Returns the non-negative cost charged once when a gap is opened.
    #[must_use]
    pub const fn gap_opening(self) -> Score {
        self.gap_opening
    }

    /// Returns the non-negative cost charged for every symbol in a gap.
    #[must_use]
    pub const fn gap_extension(self) -> Score {
        self.gap_extension
    }

    /// Returns the exactly converted workspace affine score.
    #[must_use]
    pub const fn affine_gap(self) -> AffineGap {
        self.gap
    }

    /// Computes `opening + length * extension` exactly in a wide integer.
    #[must_use]
    pub fn gap_cost(self, length: usize) -> WideScore {
        if length == 0 {
            return 0;
        }
        WideScore::from(self.gap_opening)
            + (length as WideScore) * WideScore::from(self.gap_extension)
    }
}

impl<A: Alphabet> SubstitutionScorer<A> for PoastaAffineCosts {
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score {
        if first == second { 0 } else { -self.mismatch }
    }

    fn score_bounds(&self) -> Option<(Score, Score)> {
        Some((-self.mismatch, 0))
    }
}

/// Checked element-buffer model for scalar global affine sequence-to-DAG alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AffinePoaAlignmentMemoryEstimate {
    cells: usize,
    wide_score_matrix_bytes: usize,
    traceback_matrix_bytes: usize,
    node_index_bytes: usize,
    maximum_operation_bytes: usize,
    maximum_aligned_pair_bytes: usize,
    peak_buffer_bytes: usize,
}

impl AffinePoaAlignmentMemoryEstimate {
    /// Returns matrix cells including the virtual-source row and query column zero.
    #[must_use]
    pub const fn cells(self) -> usize {
        self.cells
    }

    /// Returns combined bytes for the three full `i128` state matrices.
    #[must_use]
    pub const fn wide_score_matrix_bytes(self) -> usize {
        self.wide_score_matrix_bytes
    }

    /// Returns combined bytes for match, deletion, and insertion traces.
    #[must_use]
    pub const fn traceback_matrix_bytes(self) -> usize {
        self.traceback_matrix_bytes
    }

    /// Returns bytes for topological order and the inverse node-to-row map.
    #[must_use]
    pub const fn node_index_bytes(self) -> usize {
        self.node_index_bytes
    }

    /// Returns bytes reserved for the maximum forward operation count.
    #[must_use]
    pub const fn maximum_operation_bytes(self) -> usize {
        self.maximum_operation_bytes
    }

    /// Returns bytes reserved for the maximum aligned-pair projection.
    #[must_use]
    pub const fn maximum_aligned_pair_bytes(self) -> usize {
        self.maximum_aligned_pair_bytes
    }

    /// Returns the modeled simultaneous element-buffer peak.
    #[must_use]
    pub const fn peak_buffer_bytes(self) -> usize {
        self.peak_buffer_bytes
    }
}

/// One deterministic exact global affine alignment of a query to a DAG path.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GlobalAffineSequenceGraphAlignment {
    score: Score,
    query_len: usize,
    operations: Box<[PoaOperation]>,
    aligned_pairs: Box<[PoaAlignedPair]>,
    memory: AffinePoaAlignmentMemoryEstimate,
}

impl GlobalAffineSequenceGraphAlignment {
    /// Returns the exact maximizing alignment score.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Returns the complete consumed query range.
    #[must_use]
    pub const fn query_range(&self) -> Range<usize> {
        0..self.query_len
    }

    /// Returns canonical forward traceback operations.
    #[must_use]
    pub const fn operations(&self) -> &[PoaOperation] {
        &self.operations
    }

    /// Returns the diagonal graph-node/query-residue projection.
    #[must_use]
    pub const fn aligned_pairs(&self) -> &[PoaAlignedPair] {
        &self.aligned_pairs
    }

    /// Returns exact algorithm-owned buffer accounting for this call.
    #[must_use]
    pub const fn memory_estimate(&self) -> AffinePoaAlignmentMemoryEstimate {
        self.memory
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        Score,
        usize,
        Box<[PoaOperation]>,
        Box<[PoaAlignedPair]>,
        AffinePoaAlignmentMemoryEstimate,
    ) {
        (
            self.score,
            self.query_len,
            self.operations,
            self.aligned_pairs,
            self.memory,
        )
    }

    pub(crate) fn integration_targets(
        &self,
        query_len: usize,
    ) -> Result<Vec<Option<PoaNodeId>>, PoaError> {
        if query_len != self.query_len {
            return verification("affine integration target length differs from the query");
        }
        let mut targets = try_filled(query_len, None, "affine POA integration targets")?;
        for pair in &self.aligned_pairs {
            if pair.query_index() >= query_len
                || targets[pair.query_index()].replace(pair.node()).is_some()
            {
                return verification(
                    "affine aligned-pair mapping is not one-to-one with the query",
                );
            }
        }
        Ok(targets)
    }
}

/// Independent replay evidence for a global affine sequence-to-DAG alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AffinePoaVerificationReport {
    recomputed_score: Score,
    aligned_operations: usize,
    mismatch_operations: usize,
    graph_only_operations: usize,
    query_only_operations: usize,
    gap_runs: usize,
    selected_graph_nodes: usize,
}

impl AffinePoaVerificationReport {
    /// Returns the score recomputed without using DP or traceback buffers.
    #[must_use]
    pub const fn recomputed_score(self) -> Score {
        self.recomputed_score
    }

    /// Returns the number of aligned graph/query pairs.
    #[must_use]
    pub const fn aligned_operations(self) -> usize {
        self.aligned_operations
    }

    /// Returns the number of aligned pairs whose encoded symbols differ.
    #[must_use]
    pub const fn mismatch_operations(self) -> usize {
        self.mismatch_operations
    }

    /// Returns the number of graph-only symbols.
    #[must_use]
    pub const fn graph_only_operations(self) -> usize {
        self.graph_only_operations
    }

    /// Returns the number of query-only symbols.
    #[must_use]
    pub const fn query_only_operations(self) -> usize {
        self.query_only_operations
    }

    /// Returns the number of maximal affine gap runs.
    #[must_use]
    pub const fn gap_runs(self) -> usize {
        self.gap_runs
    }

    /// Returns the number of real graph nodes on the selected source-to-sink path.
    #[must_use]
    pub const fn selected_graph_nodes(self) -> usize {
        self.selected_graph_nodes
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
enum MatchTraceKind {
    Unreachable = 0,
    Start = 1,
    Diagonal = 2,
    CloseDeletion = 3,
    CloseInsertion = 4,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(transparent)]
struct MatchTrace(usize);

impl MatchTrace {
    const KIND_SHIFT: u32 = usize::BITS - 3;
    const PREDECESSOR_MASK: usize = usize::MAX >> 3;
    const UNREACHABLE: Self = Self(0);
    const START: Self = Self((MatchTraceKind::Start as usize) << Self::KIND_SHIFT);
    const CLOSE_DELETION: Self = Self((MatchTraceKind::CloseDeletion as usize) << Self::KIND_SHIFT);
    const CLOSE_INSERTION: Self =
        Self((MatchTraceKind::CloseInsertion as usize) << Self::KIND_SHIFT);

    const fn diagonal(predecessor_row: usize) -> Self {
        Self(((MatchTraceKind::Diagonal as usize) << Self::KIND_SHIFT) | (predecessor_row + 1))
    }

    const fn kind(self) -> MatchTraceKind {
        match self.0 >> Self::KIND_SHIFT {
            0 => MatchTraceKind::Unreachable,
            1 => MatchTraceKind::Start,
            2 => MatchTraceKind::Diagonal,
            3 => MatchTraceKind::CloseDeletion,
            _ => MatchTraceKind::CloseInsertion,
        }
    }

    const fn predecessor_row(self) -> Option<usize> {
        let encoded = self.0 & Self::PREDECESSOR_MASK;
        if encoded == 0 {
            None
        } else {
            Some(encoded - 1)
        }
    }

    const fn is_reachable(self) -> bool {
        !matches!(self.kind(), MatchTraceKind::Unreachable)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
enum DeletionTraceKind {
    Unreachable = 0,
    Open = 1,
    Extend = 2,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(transparent)]
struct DeletionTrace(usize);

impl DeletionTrace {
    const KIND_SHIFT: u32 = usize::BITS - 2;
    const PREDECESSOR_MASK: usize = usize::MAX >> 2;
    const UNREACHABLE: Self = Self(0);

    const fn open(predecessor_row: usize) -> Self {
        Self(((DeletionTraceKind::Open as usize) << Self::KIND_SHIFT) | (predecessor_row + 1))
    }

    const fn extend(predecessor_row: usize) -> Self {
        Self(((DeletionTraceKind::Extend as usize) << Self::KIND_SHIFT) | (predecessor_row + 1))
    }

    const fn kind(self) -> DeletionTraceKind {
        match self.0 >> Self::KIND_SHIFT {
            0 => DeletionTraceKind::Unreachable,
            1 => DeletionTraceKind::Open,
            _ => DeletionTraceKind::Extend,
        }
    }

    const fn predecessor_row(self) -> Option<usize> {
        let encoded = self.0 & Self::PREDECESSOR_MASK;
        if encoded == 0 {
            None
        } else {
            Some(encoded - 1)
        }
    }

    const fn is_reachable(self) -> bool {
        !matches!(self.kind(), DeletionTraceKind::Unreachable)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum InsertionTrace {
    Unreachable = 0,
    Open = 1,
    Extend = 2,
}

impl InsertionTrace {
    const fn is_reachable(self) -> bool {
        !matches!(self, Self::Unreachable)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TracebackState {
    Match,
    Deletion,
    Insertion,
}

impl TracebackState {
    const fn name(self) -> &'static str {
        match self {
            Self::Match => "M",
            Self::Deletion => "D",
            Self::Insertion => "I",
        }
    }
}

/// Scalar exact global affine alignment of one query to a path through a POA DAG.
#[derive(Clone, Debug)]
pub struct GlobalAffinePoa<S> {
    substitution: S,
    gap: AffineGap,
}

impl<S> GlobalAffinePoa<S> {
    /// Creates the generic score-side affine sequence-to-DAG oracle.
    #[must_use]
    pub const fn new(substitution: S, gap: AffineGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured length-one-open affine gap score.
    #[must_use]
    pub const fn gap(&self) -> AffineGap {
        self.gap
    }

    /// Computes exact element-buffer requirements without allocating them.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError::DimensionOverflow`] if a dimension, packed predecessor,
    /// operation count, byte count, or matrix cell count is unrepresentable.
    pub fn memory_requirement(
        graph_nodes: usize,
        query_len: usize,
    ) -> Result<AffinePoaAlignmentMemoryEstimate, PoaError> {
        let rows = graph_nodes
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        if rows > MatchTrace::PREDECESSOR_MASK || rows > DeletionTrace::PREDECESSOR_MASK {
            return Err(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            });
        }
        let columns = query_len
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        let cells = rows
            .checked_mul(columns)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        let maximum_operations =
            graph_nodes
                .checked_add(query_len)
                .ok_or(PoaError::DimensionOverflow {
                    graph_nodes,
                    query_len,
                })?;
        let one_score_matrix = checked_bytes::<WideScore>(cells, graph_nodes, query_len)?;
        let wide_score_matrix_bytes =
            one_score_matrix
                .checked_mul(3)
                .ok_or(PoaError::DimensionOverflow {
                    graph_nodes,
                    query_len,
                })?;
        let traceback_matrix_bytes = checked_bytes::<MatchTrace>(cells, graph_nodes, query_len)?
            .checked_add(checked_bytes::<DeletionTrace>(
                cells,
                graph_nodes,
                query_len,
            )?)
            .and_then(|bytes| {
                checked_bytes::<InsertionTrace>(cells, graph_nodes, query_len)
                    .ok()
                    .and_then(|insertion| bytes.checked_add(insertion))
            })
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        let node_index_bytes = checked_bytes::<usize>(graph_nodes, graph_nodes, query_len)?
            .checked_mul(2)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        let maximum_operation_bytes =
            checked_bytes::<PoaOperation>(maximum_operations, graph_nodes, query_len)?;
        let maximum_aligned_pair_bytes =
            checked_bytes::<PoaAlignedPair>(graph_nodes.min(query_len), graph_nodes, query_len)?;
        let peak_buffer_bytes = wide_score_matrix_bytes
            .checked_add(traceback_matrix_bytes)
            .and_then(|bytes| bytes.checked_add(node_index_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_operation_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_aligned_pair_bytes))
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        Ok(AffinePoaAlignmentMemoryEstimate {
            cells,
            wide_score_matrix_bytes,
            traceback_matrix_bytes,
            node_index_bytes,
            maximum_operation_bytes,
            maximum_aligned_pair_bytes,
            peak_buffer_bytes,
        })
    }

    /// Aligns a complete query globally to one complete source-to-sink DAG path.
    ///
    /// Ties follow the frozen profile: gap extension before opening, smaller
    /// predecessor IDs within a transition, diagonal before deletion before
    /// insertion, and the smallest sink ID at the endpoint.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for an invalid graph, unrepresentable dimension or
    /// optimum, allocation failure, invalid traceback, or failed replay.
    #[allow(clippy::too_many_lines)]
    pub fn align<A>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
    ) -> Result<GlobalAffineSequenceGraphAlignment, PoaError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let (order, _) = graph.validate_with_orders()?;
        let graph_nodes = graph.node_count();
        let query_len = query.len();
        let memory = Self::memory_requirement(graph_nodes, query_len)?;
        let columns = query_len + 1;
        let maximum_operations = graph_nodes + query_len;

        let mut node_to_row = try_filled(graph_nodes, usize::MAX, "affine POA node-to-row map")?;
        for (offset, node) in order.iter().copied().enumerate() {
            node_to_row[node.index()] = offset + 1;
        }

        let mut match_scores = try_filled(memory.cells, 0, "affine POA M scores")?;
        let mut deletion_scores = try_filled(memory.cells, 0, "affine POA D scores")?;
        let mut insertion_scores = try_filled(memory.cells, 0, "affine POA I scores")?;
        let mut match_trace =
            try_filled(memory.cells, MatchTrace::UNREACHABLE, "affine POA M trace")?;
        let mut deletion_trace = try_filled(
            memory.cells,
            DeletionTrace::UNREACHABLE,
            "affine POA D trace",
        )?;
        let mut insertion_trace = try_filled(
            memory.cells,
            InsertionTrace::Unreachable,
            "affine POA I trace",
        )?;

        match_trace[0] = MatchTrace::START;
        let open_score = WideScore::from(self.gap.open_score());
        let extension_score = WideScore::from(self.gap.extension_score());

        for column in 1..=query_len {
            let index = column;
            let previous = index - 1;
            let mut best: Option<(WideScore, InsertionTrace)> = None;
            if insertion_trace[previous].is_reachable() {
                choose_strict_max(
                    &mut best,
                    insertion_scores[previous] + extension_score,
                    InsertionTrace::Extend,
                );
            }
            if match_trace[previous].is_reachable() {
                choose_strict_max(
                    &mut best,
                    match_scores[previous] + open_score,
                    InsertionTrace::Open,
                );
            }
            let Some((score, trace)) = best else {
                return affine_traceback(None, column, TracebackState::Insertion);
            };
            insertion_scores[index] = score;
            insertion_trace[index] = trace;
            match_scores[index] = score;
            match_trace[index] = MatchTrace::CLOSE_INSERTION;
        }

        for (offset, node_id) in order.iter().copied().enumerate() {
            let row = offset + 1;
            let node = graph
                .node(node_id)
                .ok_or(PoaError::UnknownNode { node: node_id })?;

            for column in 0..=query_len {
                let index = row * columns + column;

                if column > 0 {
                    let previous = index - 1;
                    let mut best: Option<(WideScore, InsertionTrace)> = None;
                    if insertion_trace[previous].is_reachable() {
                        choose_strict_max(
                            &mut best,
                            insertion_scores[previous] + extension_score,
                            InsertionTrace::Extend,
                        );
                    }
                    if match_trace[previous].is_reachable() {
                        choose_strict_max(
                            &mut best,
                            match_scores[previous] + open_score,
                            InsertionTrace::Open,
                        );
                    }
                    if let Some((score, trace)) = best {
                        insertion_scores[index] = score;
                        insertion_trace[index] = trace;
                    }
                }

                let mut deletion_best: Option<(WideScore, DeletionTrace)> = None;
                if node.predecessors().is_empty() {
                    choose_deletion_extension(
                        &mut deletion_best,
                        0,
                        column,
                        columns,
                        extension_score,
                        &deletion_scores,
                        &deletion_trace,
                    );
                    choose_deletion_open(
                        &mut deletion_best,
                        0,
                        column,
                        columns,
                        open_score,
                        &match_scores,
                        &match_trace,
                    );
                } else {
                    for predecessor in node.predecessors() {
                        let predecessor_row = node_to_row[predecessor.index()];
                        choose_deletion_extension(
                            &mut deletion_best,
                            predecessor_row,
                            column,
                            columns,
                            extension_score,
                            &deletion_scores,
                            &deletion_trace,
                        );
                    }
                    for predecessor in node.predecessors() {
                        let predecessor_row = node_to_row[predecessor.index()];
                        choose_deletion_open(
                            &mut deletion_best,
                            predecessor_row,
                            column,
                            columns,
                            open_score,
                            &match_scores,
                            &match_trace,
                        );
                    }
                }
                if let Some((score, trace)) = deletion_best {
                    deletion_scores[index] = score;
                    deletion_trace[index] = trace;
                }

                let mut match_best: Option<(WideScore, MatchTrace)> = None;
                if column > 0 {
                    let query_symbol =
                        query.get(column - 1).ok_or(PoaError::DimensionOverflow {
                            graph_nodes,
                            query_len,
                        })?;
                    let substitution =
                        WideScore::from(self.substitution.score(node.symbol(), query_symbol));
                    if node.predecessors().is_empty() {
                        choose_diagonal(
                            &mut match_best,
                            0,
                            column,
                            columns,
                            substitution,
                            &match_scores,
                            &match_trace,
                        );
                    } else {
                        for predecessor in node.predecessors() {
                            choose_diagonal(
                                &mut match_best,
                                node_to_row[predecessor.index()],
                                column,
                                columns,
                                substitution,
                                &match_scores,
                                &match_trace,
                            );
                        }
                    }
                }
                if deletion_trace[index].is_reachable() {
                    choose_strict_max(
                        &mut match_best,
                        deletion_scores[index],
                        MatchTrace::CLOSE_DELETION,
                    );
                }
                if insertion_trace[index].is_reachable() {
                    choose_strict_max(
                        &mut match_best,
                        insertion_scores[index],
                        MatchTrace::CLOSE_INSERTION,
                    );
                }
                let Some((score, trace)) = match_best else {
                    return affine_traceback(Some(node_id), column, TracebackState::Match);
                };
                match_scores[index] = score;
                match_trace[index] = trace;
            }
        }

        let (mut row, endpoint_score) = if graph_nodes == 0 {
            (0, match_scores[query_len])
        } else {
            let mut endpoint: Option<(WideScore, usize)> = None;
            for node in graph.nodes() {
                if node.successors().is_empty() {
                    let sink_row = node_to_row[node.id().index()];
                    let index = sink_row * columns + query_len;
                    if match_trace[index].is_reachable()
                        && endpoint
                            .as_ref()
                            .is_none_or(|(best, _)| match_scores[index] > *best)
                    {
                        endpoint = Some((match_scores[index], sink_row));
                    }
                }
            }
            let Some((score, sink_row)) = endpoint else {
                return affine_traceback(None, query_len, TracebackState::Match);
            };
            (sink_row, score)
        };
        let score =
            Score::try_from(endpoint_score).map_err(|_| PoaError::AffineScoreOutOfRange {
                score: endpoint_score,
            })?;

        let mut operations = Vec::new();
        operations
            .try_reserve_exact(maximum_operations)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "affine POA traceback operations",
                elements: maximum_operations,
            })?;
        let mut column = query_len;
        let mut state = TracebackState::Match;
        loop {
            let index = row * columns + column;
            match state {
                TracebackState::Match => match match_trace[index].kind() {
                    MatchTraceKind::Start if row == 0 && column == 0 => break,
                    MatchTraceKind::Diagonal if row > 0 && column > 0 => {
                        let Some(predecessor_row) = match_trace[index].predecessor_row() else {
                            return affine_traceback(order.get(row - 1).copied(), column, state);
                        };
                        if predecessor_row >= row {
                            return affine_traceback(order.get(row - 1).copied(), column, state);
                        }
                        operations.push(PoaOperation::Aligned {
                            node: order[row - 1],
                            query_index: column - 1,
                        });
                        row = predecessor_row;
                        column -= 1;
                    }
                    MatchTraceKind::CloseDeletion if row > 0 => {
                        state = TracebackState::Deletion;
                    }
                    MatchTraceKind::CloseInsertion if column > 0 => {
                        state = TracebackState::Insertion;
                    }
                    _ => {
                        return affine_traceback(
                            row.checked_sub(1)
                                .and_then(|offset| order.get(offset))
                                .copied(),
                            column,
                            state,
                        );
                    }
                },
                TracebackState::Deletion if row > 0 => {
                    let trace = deletion_trace[index];
                    let Some(predecessor_row) = trace.predecessor_row() else {
                        return affine_traceback(order.get(row - 1).copied(), column, state);
                    };
                    if predecessor_row >= row {
                        return affine_traceback(order.get(row - 1).copied(), column, state);
                    }
                    operations.push(PoaOperation::GraphOnly {
                        node: order[row - 1],
                    });
                    row = predecessor_row;
                    state = match trace.kind() {
                        DeletionTraceKind::Open => TracebackState::Match,
                        DeletionTraceKind::Extend => TracebackState::Deletion,
                        DeletionTraceKind::Unreachable => {
                            return affine_traceback(
                                order.get(row.saturating_sub(1)).copied(),
                                column,
                                TracebackState::Deletion,
                            );
                        }
                    };
                }
                TracebackState::Insertion if column > 0 => {
                    let trace = insertion_trace[index];
                    operations.push(PoaOperation::QueryOnly {
                        query_index: column - 1,
                    });
                    column -= 1;
                    state = match trace {
                        InsertionTrace::Open => TracebackState::Match,
                        InsertionTrace::Extend => TracebackState::Insertion,
                        InsertionTrace::Unreachable => {
                            return affine_traceback(
                                row.checked_sub(1)
                                    .and_then(|offset| order.get(offset))
                                    .copied(),
                                column,
                                TracebackState::Insertion,
                            );
                        }
                    };
                }
                _ => {
                    return affine_traceback(
                        row.checked_sub(1)
                            .and_then(|offset| order.get(offset))
                            .copied(),
                        column,
                        state,
                    );
                }
            }
            if operations.len() > maximum_operations {
                return affine_traceback(
                    row.checked_sub(1)
                        .and_then(|offset| order.get(offset))
                        .copied(),
                    column,
                    state,
                );
            }
        }
        operations.reverse();

        let aligned_pair_count = operations
            .iter()
            .filter(|operation| matches!(operation, PoaOperation::Aligned { .. }))
            .count();
        let mut aligned_pairs = Vec::new();
        aligned_pairs
            .try_reserve_exact(aligned_pair_count)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "affine POA aligned pairs",
                elements: aligned_pair_count,
            })?;
        aligned_pairs.extend(operations.iter().filter_map(|operation| match *operation {
            PoaOperation::Aligned { node, query_index } => {
                Some(PoaAlignedPair::new(node, query_index))
            }
            PoaOperation::GraphOnly { .. } | PoaOperation::QueryOnly { .. } => None,
        }));

        let alignment = GlobalAffineSequenceGraphAlignment {
            score,
            query_len,
            operations: operations.into_boxed_slice(),
            aligned_pairs: aligned_pairs.into_boxed_slice(),
            memory,
        };
        verify_global_affine_sequence_graph_alignment_validated(
            &alignment,
            graph,
            query,
            &self.substitution,
            self.gap,
        )?;
        Ok(alignment)
    }
}

impl GlobalAffinePoa<PoastaAffineCosts> {
    /// Creates the exact paper-cost profile with no silent gap-convention change.
    #[must_use]
    pub const fn from_poasta_costs(costs: PoastaAffineCosts) -> Self {
        Self::new(costs, costs.affine_gap())
    }
}

/// Independently replays a global affine sequence-to-DAG result.
///
/// # Errors
///
/// Returns [`PoaError`] when the graph path, complete query consumption, affine
/// run score, aligned-pair projection, or reported score is invalid.
pub fn verify_global_affine_sequence_graph_alignment<A, S>(
    alignment: &GlobalAffineSequenceGraphAlignment,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
) -> Result<AffinePoaVerificationReport, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    graph.validate()?;
    verify_global_affine_sequence_graph_alignment_validated(
        alignment,
        graph,
        query,
        substitution,
        gap,
    )
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GapDirection {
    GraphOnly,
    QueryOnly,
}

#[allow(clippy::too_many_lines)]
pub(crate) fn verify_global_affine_sequence_graph_alignment_validated<A, S>(
    alignment: &GlobalAffineSequenceGraphAlignment,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
) -> Result<AffinePoaVerificationReport, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    verify_global_affine_operations_validated(
        alignment.score,
        alignment.query_len,
        &alignment.operations,
        &alignment.aligned_pairs,
        graph,
        query,
        substitution,
        gap,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn verify_global_affine_operations_validated<A, S>(
    reported_score: Score,
    reported_query_len: usize,
    operations: &[PoaOperation],
    aligned_pairs: &[PoaAlignedPair],
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
) -> Result<AffinePoaVerificationReport, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    if reported_query_len != query.len() {
        return verification("global affine result does not expose the complete query range");
    }
    let maximum_operations =
        graph
            .node_count()
            .checked_add(query.len())
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
    if operations.len() > maximum_operations {
        return verification("global affine traceback exceeds graph-plus-query length");
    }

    let mut score = 0_i128;
    let mut previous_graph_node = None;
    let mut next_query_index = 0_usize;
    let mut aligned_pair_index = 0_usize;
    let mut aligned_operations = 0_usize;
    let mut mismatch_operations = 0_usize;
    let mut graph_only_operations = 0_usize;
    let mut query_only_operations = 0_usize;
    let mut selected_graph_nodes = 0_usize;
    let mut active_gap: Option<(GapDirection, usize)> = None;

    for operation in operations {
        match *operation {
            PoaOperation::Aligned { node, query_index } => {
                flush_gap(&mut score, &mut active_gap, gap)?;
                verify_global_graph_step(graph, previous_graph_node, node)?;
                if query_index != next_query_index {
                    return verification(
                        "aligned operation breaks complete consecutive query consumption",
                    );
                }
                let query_symbol =
                    query
                        .get(query_index)
                        .ok_or_else(|| PoaError::VerificationFailed {
                            message: format!("aligned query index {query_index} is out of range"),
                        })?;
                let graph_symbol = graph
                    .node(node)
                    .ok_or(PoaError::UnknownNode { node })?
                    .symbol();
                score += WideScore::from(substitution.score(graph_symbol, query_symbol));
                if graph_symbol != query_symbol {
                    mismatch_operations =
                        checked_increment(mismatch_operations, "affine mismatch operation count")?;
                }
                previous_graph_node = Some(node);
                next_query_index = checked_increment(next_query_index, "affine query position")?;
                aligned_operations =
                    checked_increment(aligned_operations, "affine aligned operation count")?;
                selected_graph_nodes =
                    checked_increment(selected_graph_nodes, "affine selected graph nodes")?;
                if aligned_pairs.get(aligned_pair_index)
                    != Some(&PoaAlignedPair::new(node, query_index))
                {
                    return verification("affine aligned-pair projection disagrees with traceback");
                }
                aligned_pair_index =
                    checked_increment(aligned_pair_index, "affine aligned-pair index")?;
            }
            PoaOperation::GraphOnly { node } => {
                verify_global_graph_step(graph, previous_graph_node, node)?;
                continue_gap(&mut score, &mut active_gap, GapDirection::GraphOnly, gap)?;
                previous_graph_node = Some(node);
                graph_only_operations =
                    checked_increment(graph_only_operations, "affine graph-only operation count")?;
                selected_graph_nodes =
                    checked_increment(selected_graph_nodes, "affine selected graph nodes")?;
            }
            PoaOperation::QueryOnly { query_index } => {
                if query_index != next_query_index || query.get(query_index).is_none() {
                    return verification(
                        "query-only operation breaks complete consecutive query consumption",
                    );
                }
                continue_gap(&mut score, &mut active_gap, GapDirection::QueryOnly, gap)?;
                next_query_index = checked_increment(next_query_index, "affine query position")?;
                query_only_operations =
                    checked_increment(query_only_operations, "affine query-only operation count")?;
            }
        }
    }
    flush_gap(&mut score, &mut active_gap, gap)?;

    if next_query_index != query.len() {
        return verification("global affine traceback does not consume the complete query");
    }
    if aligned_pair_index != aligned_pairs.len() {
        return verification("affine aligned-pair projection has trailing entries");
    }
    if graph.node_count() == 0 {
        if previous_graph_node.is_some() || selected_graph_nodes != 0 {
            return verification("an empty graph result consumes a real graph node");
        }
    } else {
        let Some(last_node) = previous_graph_node else {
            return verification("a non-empty graph result selects no source-to-sink path");
        };
        if !graph
            .node(last_node)
            .ok_or(PoaError::UnknownNode { node: last_node })?
            .successors()
            .is_empty()
        {
            return verification("global affine graph path does not terminate at a sink");
        }
    }
    if selected_graph_nodes > graph.node_count() {
        return verification("global affine graph path repeats more nodes than the DAG contains");
    }

    let recomputed_score =
        Score::try_from(score).map_err(|_| PoaError::AffineScoreOutOfRange { score })?;
    if recomputed_score != reported_score {
        return verification(format!(
            "global affine replay score {recomputed_score} differs from reported score {reported_score}"
        ));
    }

    let gap_runs = operations
        .iter()
        .fold((None, 0_usize), |(previous, runs), operation| {
            let direction = match operation {
                PoaOperation::GraphOnly { .. } => Some(GapDirection::GraphOnly),
                PoaOperation::QueryOnly { .. } => Some(GapDirection::QueryOnly),
                PoaOperation::Aligned { .. } => None,
            };
            let next_runs = if direction.is_some() && direction != previous {
                runs + 1
            } else {
                runs
            };
            (direction, next_runs)
        })
        .1;

    Ok(AffinePoaVerificationReport {
        recomputed_score,
        aligned_operations,
        mismatch_operations,
        graph_only_operations,
        query_only_operations,
        gap_runs,
        selected_graph_nodes,
    })
}

fn choose_deletion_extension(
    best: &mut Option<(WideScore, DeletionTrace)>,
    predecessor_row: usize,
    column: usize,
    columns: usize,
    extension_score: WideScore,
    scores: &[WideScore],
    trace: &[DeletionTrace],
) {
    let predecessor = predecessor_row * columns + column;
    if trace[predecessor].is_reachable() {
        choose_strict_max(
            best,
            scores[predecessor] + extension_score,
            DeletionTrace::extend(predecessor_row),
        );
    }
}

fn choose_deletion_open(
    best: &mut Option<(WideScore, DeletionTrace)>,
    predecessor_row: usize,
    column: usize,
    columns: usize,
    open_score: WideScore,
    scores: &[WideScore],
    trace: &[MatchTrace],
) {
    let predecessor = predecessor_row * columns + column;
    if trace[predecessor].is_reachable() {
        choose_strict_max(
            best,
            scores[predecessor] + open_score,
            DeletionTrace::open(predecessor_row),
        );
    }
}

fn choose_diagonal(
    best: &mut Option<(WideScore, MatchTrace)>,
    predecessor_row: usize,
    column: usize,
    columns: usize,
    substitution: WideScore,
    scores: &[WideScore],
    trace: &[MatchTrace],
) {
    let predecessor = predecessor_row * columns + column - 1;
    if trace[predecessor].is_reachable() {
        choose_strict_max(
            best,
            scores[predecessor] + substitution,
            MatchTrace::diagonal(predecessor_row),
        );
    }
}

fn choose_strict_max<T: Copy>(best: &mut Option<(WideScore, T)>, score: WideScore, value: T) {
    if best.as_ref().is_none_or(|(current, _)| score > *current) {
        *best = Some((score, value));
    }
}

fn verify_global_graph_step<A: Alphabet>(
    graph: &PoaGraph<A>,
    previous: Option<PoaNodeId>,
    node: PoaNodeId,
) -> Result<(), PoaError> {
    let current = graph.node(node).ok_or(PoaError::UnknownNode { node })?;
    match previous {
        Some(previous) if !graph.has_edge(previous, node) => verification(format!(
            "global affine traceback has no graph edge {} -> {}",
            previous.index(),
            node.index()
        )),
        None if !current.predecessors().is_empty() => {
            verification("global affine graph path does not begin at a source")
        }
        Some(_) | None => Ok(()),
    }
}

fn continue_gap(
    score: &mut WideScore,
    active: &mut Option<(GapDirection, usize)>,
    direction: GapDirection,
    gap: AffineGap,
) -> Result<(), PoaError> {
    match active {
        Some((current, length)) if *current == direction => {
            *length = checked_increment(*length, "affine gap-run length")?;
        }
        Some(_) => {
            flush_gap(score, active, gap)?;
            *active = Some((direction, 1));
        }
        None => *active = Some((direction, 1)),
    }
    Ok(())
}

fn flush_gap(
    score: &mut WideScore,
    active: &mut Option<(GapDirection, usize)>,
    gap: AffineGap,
) -> Result<(), PoaError> {
    let Some((_, length)) = active.take() else {
        return Ok(());
    };
    let extension_count = length.checked_sub(1).ok_or(PoaError::ArithmeticOverflow {
        context: "affine gap extension count",
    })?;
    *score += WideScore::from(gap.open_score())
        + (extension_count as WideScore) * WideScore::from(gap.extension_score());
    Ok(())
}

fn checked_increment(value: usize, context: &'static str) -> Result<usize, PoaError> {
    value
        .checked_add(1)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_bytes<T>(
    elements: usize,
    graph_nodes: usize,
    query_len: usize,
) -> Result<usize, PoaError> {
    elements
        .checked_mul(size_of::<T>())
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes,
            query_len,
        })
}

fn try_filled<T: Clone>(
    elements: usize,
    value: T,
    buffer: &'static str,
) -> Result<Vec<T>, PoaError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| PoaError::AllocationFailed { buffer, elements })?;
    values.resize(elements, value);
    Ok(values)
}

fn affine_traceback<T>(
    node: Option<PoaNodeId>,
    column: usize,
    state: TracebackState,
) -> Result<T, PoaError> {
    Err(PoaError::AffineTracebackInvariant {
        node,
        column,
        state: state.name(),
    })
}

fn verification<T>(message: impl Into<String>) -> Result<T, PoaError> {
    Err(PoaError::VerificationFailed {
        message: message.into(),
    })
}

#[cfg(test)]
mod tests {
    use core::mem::size_of;

    use align_core::{
        AffineGap, Dna, DnaSequence, MatchMismatch, Protein, ProteinSequence, Rna, RnaSequence,
        Score, Sequence,
    };
    use align_msa::SequenceRecord;

    use super::{
        AffinePoaVerificationReport, DeletionTrace, GlobalAffinePoa, InsertionTrace, MatchTrace,
        PoastaAffineCosts, verify_global_affine_sequence_graph_alignment,
    };
    use crate::{PoaError, PoaGraph, PoaNodeId, PoaOperation};

    fn record<A: align_core::Alphabet>(identifier: &str, sequence: &str) -> SequenceRecord<A> {
        SequenceRecord::new(identifier, Sequence::<A>::from_ascii(sequence).unwrap()).unwrap()
    }

    fn linear_graph<A: align_core::Alphabet>(sequence: &str) -> PoaGraph<A> {
        let mut graph = PoaGraph::new();
        let source = record("source", sequence);
        graph
            .append_record(&source, &vec![None; source.sequence().len()])
            .unwrap();
        graph.validate().unwrap();
        graph
    }

    fn paper_costs() -> PoastaAffineCosts {
        PoastaAffineCosts::new(4, 6, 2).unwrap()
    }

    fn paper_aligner() -> GlobalAffinePoa<PoastaAffineCosts> {
        GlobalAffinePoa::from_poasta_costs(paper_costs())
    }

    #[test]
    fn published_costs_convert_without_changing_gap_convention() {
        let costs = paper_costs();
        assert_eq!(costs.mismatch(), 4);
        assert_eq!(costs.gap_opening(), 6);
        assert_eq!(costs.gap_extension(), 2);
        assert_eq!(costs.affine_gap().open_score(), -8);
        assert_eq!(costs.affine_gap().extension_score(), -2);
        assert_eq!(costs.gap_cost(0), 0);
        assert_eq!(costs.gap_cost(1), 8);
        assert_eq!(costs.gap_cost(4), 14);
        assert_eq!(costs.affine_gap().score_for_length(1), Some(-8));
        assert_eq!(costs.affine_gap().score_for_length(4), Some(-14));
    }

    #[test]
    fn invalid_published_costs_are_typed_errors() {
        assert_eq!(
            PoastaAffineCosts::new(-1, 6, 2),
            Err(PoaError::NegativeAffineCost {
                parameter: "mismatch",
                cost: -1,
            })
        );
        assert_eq!(
            PoastaAffineCosts::new(4, 2, 2),
            Err(PoaError::AffineGapOpeningNotGreaterThanExtension {
                opening: 2,
                extension: 2,
            })
        );
        assert_eq!(
            PoastaAffineCosts::new(4, Score::MAX, Score::MAX - 1),
            Err(PoaError::AffineCostOverflow {
                opening: Score::MAX,
                extension: Score::MAX - 1,
            })
        );
    }

    #[test]
    fn identical_linear_dna_uses_only_diagonals() {
        let graph = linear_graph::<Dna>("ACGT");
        let query = DnaSequence::from_ascii("ACGT").unwrap();
        let alignment = paper_aligner().align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), 0);
        assert_eq!(alignment.query_range(), 0..4);
        assert_eq!(
            alignment.operations(),
            &[
                PoaOperation::Aligned {
                    node: PoaNodeId::new(0),
                    query_index: 0,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(1),
                    query_index: 1,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(2),
                    query_index: 2,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(3),
                    query_index: 3,
                },
            ]
        );
        let report = verify_global_affine_sequence_graph_alignment(
            &alignment,
            &graph,
            &query,
            &paper_costs(),
            paper_costs().affine_gap(),
        )
        .unwrap();
        assert_eq!(
            report,
            AffinePoaVerificationReport {
                recomputed_score: 0,
                aligned_operations: 4,
                mismatch_operations: 0,
                graph_only_operations: 0,
                query_only_operations: 0,
                gap_runs: 0,
                selected_graph_nodes: 4,
            }
        );
    }

    #[test]
    fn internal_deletion_charges_one_open_plus_one_extension() {
        let graph = linear_graph::<Dna>("ACGT");
        let query = DnaSequence::from_ascii("AGT").unwrap();
        let alignment = paper_aligner().align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), -8);
        assert_eq!(
            alignment.operations(),
            &[
                PoaOperation::Aligned {
                    node: PoaNodeId::new(0),
                    query_index: 0,
                },
                PoaOperation::GraphOnly {
                    node: PoaNodeId::new(1),
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(2),
                    query_index: 1,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(3),
                    query_index: 2,
                },
            ]
        );
    }

    #[test]
    fn internal_insertion_charges_one_open_plus_one_extension() {
        let graph = linear_graph::<Dna>("ACGT");
        let query = DnaSequence::from_ascii("ACAGT").unwrap();
        let alignment = paper_aligner().align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), -8);
        assert!(matches!(
            alignment.operations()[2],
            PoaOperation::QueryOnly { query_index: 2 }
        ));
    }

    #[test]
    fn leading_and_trailing_gap_runs_use_affine_lengths() {
        let graph = linear_graph::<Dna>("ACGT");
        let leading = paper_aligner()
            .align(&graph, &DnaSequence::from_ascii("TTACGT").unwrap())
            .unwrap();
        assert_eq!(leading.score(), -10);
        assert!(matches!(
            leading.operations()[0],
            PoaOperation::QueryOnly { query_index: 0 }
        ));
        assert!(matches!(
            leading.operations()[1],
            PoaOperation::QueryOnly { query_index: 1 }
        ));

        let trailing = paper_aligner()
            .align(&graph, &DnaSequence::from_ascii("ACGTAA").unwrap())
            .unwrap();
        assert_eq!(trailing.score(), -10);
        assert!(matches!(
            trailing.operations()[4],
            PoaOperation::QueryOnly { query_index: 4 }
        ));
        assert!(matches!(
            trailing.operations()[5],
            PoaOperation::QueryOnly { query_index: 5 }
        ));
    }

    #[test]
    fn strict_cost_report_reconstructs_the_published_objective() {
        let graph = linear_graph::<Dna>("ACGT");
        let query = DnaSequence::from_ascii("ATT").unwrap();
        let costs = paper_costs();
        let alignment = GlobalAffinePoa::from_poasta_costs(costs)
            .align(&graph, &query)
            .unwrap();
        let report = verify_global_affine_sequence_graph_alignment(
            &alignment,
            &graph,
            &query,
            &costs,
            costs.affine_gap(),
        )
        .unwrap();
        let gap_symbols = report.graph_only_operations() + report.query_only_operations();
        let published_cost = i128::from(costs.mismatch()) * (report.mismatch_operations() as i128)
            + i128::from(costs.gap_opening()) * (report.gap_runs() as i128)
            + i128::from(costs.gap_extension()) * (gap_symbols as i128);
        assert_eq!(i128::from(alignment.score()), -published_cost);
    }

    #[test]
    fn all_empty_boundaries_have_canonical_global_results() {
        let empty_graph = PoaGraph::<Dna>::new();
        let empty_query = DnaSequence::from_ascii("").unwrap();
        let both_empty = paper_aligner().align(&empty_graph, &empty_query).unwrap();
        assert_eq!(both_empty.score(), 0);
        assert!(both_empty.operations().is_empty());

        let query = DnaSequence::from_ascii("ACG").unwrap();
        let no_graph = paper_aligner().align(&empty_graph, &query).unwrap();
        assert_eq!(no_graph.score(), -12);
        assert_eq!(no_graph.operations().len(), 3);
        assert!(
            no_graph
                .operations()
                .iter()
                .all(|operation| matches!(operation, PoaOperation::QueryOnly { .. }))
        );

        let graph = linear_graph::<Dna>("ACG");
        let no_query = paper_aligner().align(&graph, &empty_query).unwrap();
        assert_eq!(no_query.score(), -12);
        assert_eq!(no_query.operations().len(), 3);
        assert!(
            no_query
                .operations()
                .iter()
                .all(|operation| matches!(operation, PoaOperation::GraphOnly { .. }))
        );
    }

    #[test]
    fn mismatches_can_lose_to_adjacent_opposite_gap_runs() {
        let graph = linear_graph::<Dna>("A");
        let query = DnaSequence::from_ascii("C").unwrap();
        let gap = AffineGap::new(-2, -1).unwrap();
        let aligner = GlobalAffinePoa::new(MatchMismatch::new(0, -100), gap);
        let alignment = aligner.align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), -4);
        assert_eq!(
            alignment.operations(),
            &[
                PoaOperation::QueryOnly { query_index: 0 },
                PoaOperation::GraphOnly {
                    node: PoaNodeId::new(0),
                },
            ]
        );
        let report = verify_global_affine_sequence_graph_alignment(
            &alignment,
            &graph,
            &query,
            &MatchMismatch::new(0, -100),
            gap,
        )
        .unwrap();
        assert_eq!(report.gap_runs(), 2);
    }

    #[test]
    fn branch_alignment_selects_the_best_graph_path() {
        let mut graph = linear_graph::<Dna>("ACG");
        graph
            .append_record(
                &record("allele", "ATG"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    Some(PoaNodeId::new(2)),
                ],
            )
            .unwrap();
        graph.validate().unwrap();
        let query = DnaSequence::from_ascii("ATG").unwrap();
        let alignment = paper_aligner().align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), 0);
        assert_eq!(
            alignment
                .operations()
                .iter()
                .filter_map(|operation| match operation {
                    PoaOperation::Aligned { node, .. } => Some(node.index()),
                    _ => None,
                })
                .collect::<Vec<_>>(),
            vec![0, 3, 2]
        );
    }

    #[test]
    fn equal_disconnected_paths_choose_the_smallest_sink_id() {
        let mut graph = linear_graph::<Dna>("A");
        graph
            .append_record(&record("second", "A"), &[None])
            .unwrap();
        graph.validate().unwrap();
        let query = DnaSequence::from_ascii("A").unwrap();
        let alignment = paper_aligner().align(&graph, &query).unwrap();
        assert_eq!(
            alignment.operations(),
            &[PoaOperation::Aligned {
                node: PoaNodeId::new(0),
                query_index: 0,
            }]
        );
    }

    #[test]
    fn protein_and_rna_share_the_same_exact_kernel() {
        let protein_graph = linear_graph::<Protein>("ARND");
        let protein = ProteinSequence::from_ascii("ARXD").unwrap();
        let protein_alignment = paper_aligner().align(&protein_graph, &protein).unwrap();
        assert_eq!(protein_alignment.score(), -4);

        let rna_graph = linear_graph::<Rna>("ACGU");
        let rna = RnaSequence::from_ascii("AGU").unwrap();
        let rna_alignment = paper_aligner().align(&rna_graph, &rna).unwrap();
        assert_eq!(rna_alignment.score(), -8);
    }

    #[test]
    fn exact_wide_optimum_reports_public_score_overflow() {
        let graph = linear_graph::<Dna>("AA");
        let query = DnaSequence::from_ascii("AA").unwrap();
        let aligner = GlobalAffinePoa::new(
            MatchMismatch::new(Score::MAX, -1),
            AffineGap::new(-1, 0).unwrap(),
        );
        assert_eq!(
            aligner.align(&graph, &query),
            Err(PoaError::AffineScoreOutOfRange {
                score: i128::from(Score::MAX) * 2,
            })
        );
    }

    #[test]
    fn memory_model_matches_every_allocated_matrix_kind() {
        let memory = GlobalAffinePoa::<MatchMismatch>::memory_requirement(3, 4).unwrap();
        assert_eq!(memory.cells(), 20);
        assert_eq!(memory.wide_score_matrix_bytes(), 20 * size_of::<i128>() * 3);
        assert_eq!(
            memory.traceback_matrix_bytes(),
            20 * (size_of::<MatchTrace>()
                + size_of::<DeletionTrace>()
                + size_of::<InsertionTrace>())
        );
        assert_eq!(memory.node_index_bytes(), 6 * size_of::<usize>());
        assert!(GlobalAffinePoa::<MatchMismatch>::memory_requirement(1, usize::MAX).is_err());
    }

    #[test]
    fn independent_replay_rejects_score_and_path_tampering() {
        let graph = linear_graph::<Dna>("AC");
        let query = DnaSequence::from_ascii("AC").unwrap();
        let mut alignment = paper_aligner().align(&graph, &query).unwrap();
        alignment.score = -1;
        assert!(
            verify_global_affine_sequence_graph_alignment(
                &alignment,
                &graph,
                &query,
                &paper_costs(),
                paper_costs().affine_gap(),
            )
            .is_err()
        );

        let mut alignment = paper_aligner().align(&graph, &query).unwrap();
        alignment.operations.swap(0, 1);
        assert!(
            verify_global_affine_sequence_graph_alignment(
                &alignment,
                &graph,
                &query,
                &paper_costs(),
                paper_costs().affine_gap(),
            )
            .is_err()
        );
    }
}
