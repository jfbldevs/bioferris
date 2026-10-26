use core::mem::size_of;
use core::ops::Range;

use align_core::{Alphabet, LinearGap, Score, Sequence, SubstitutionScorer};

use crate::{PoaAlignedPair, PoaError, PoaGraph, PoaNodeId, PoaOperation};

/// Checked element-buffer model for scalar local sequence-to-DAG alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoaAlignmentMemoryEstimate {
    cells: usize,
    score_matrix_bytes: usize,
    traceback_matrix_bytes: usize,
    node_index_bytes: usize,
    maximum_operation_bytes: usize,
    maximum_aligned_pair_bytes: usize,
    peak_buffer_bytes: usize,
}

impl PoaAlignmentMemoryEstimate {
    /// Returns the number of matrix cells, including query column zero.
    #[must_use]
    pub const fn cells(self) -> usize {
        self.cells
    }

    /// Returns bytes occupied by the full score matrix.
    #[must_use]
    pub const fn score_matrix_bytes(self) -> usize {
        self.score_matrix_bytes
    }

    /// Returns bytes occupied by the full traceback matrix.
    #[must_use]
    pub const fn traceback_matrix_bytes(self) -> usize {
        self.traceback_matrix_bytes
    }

    /// Returns bytes for topological order and the node-to-row map.
    #[must_use]
    pub const fn node_index_bytes(self) -> usize {
        self.node_index_bytes
    }

    /// Returns bytes reserved for the maximum possible operation count.
    #[must_use]
    pub const fn maximum_operation_bytes(self) -> usize {
        self.maximum_operation_bytes
    }

    /// Returns bytes reserved for the maximum number of aligned-pair records.
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

/// One deterministic exact local alignment of a sequence to a partial-order DAG.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceGraphAlignment {
    score: Score,
    query_range: Range<usize>,
    operations: Box<[PoaOperation]>,
    aligned_pairs: Box<[PoaAlignedPair]>,
    memory: PoaAlignmentMemoryEstimate,
}

impl SequenceGraphAlignment {
    /// Returns the exact local alignment score.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Returns the half-open query interval consumed by the local traceback.
    #[must_use]
    pub fn query_range(&self) -> Range<usize> {
        self.query_range.clone()
    }

    /// Returns forward traceback operations.
    #[must_use]
    pub const fn operations(&self) -> &[PoaOperation] {
        &self.operations
    }

    /// Returns homologously aligned graph-node/query-residue pairs.
    #[must_use]
    pub const fn aligned_pairs(&self) -> &[PoaAlignedPair] {
        &self.aligned_pairs
    }

    /// Returns the checked memory model used for this call.
    #[must_use]
    pub const fn memory_estimate(&self) -> PoaAlignmentMemoryEstimate {
        self.memory
    }

    pub(crate) fn canonical_empty(memory: PoaAlignmentMemoryEstimate) -> Self {
        Self {
            score: 0,
            query_range: 0..0,
            operations: Box::new([]),
            aligned_pairs: Box::new([]),
            memory,
        }
    }

    pub(crate) fn integration_targets(
        &self,
        query_len: usize,
    ) -> Result<Vec<Option<PoaNodeId>>, PoaError> {
        let mut targets = Vec::new();
        targets
            .try_reserve_exact(query_len)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "POA integration targets",
                elements: query_len,
            })?;
        targets.resize(query_len, None);
        for pair in &self.aligned_pairs {
            if pair.query_index() >= query_len
                || targets[pair.query_index()].replace(pair.node()).is_some()
            {
                return verification("aligned-pair mapping is not one-to-one with the query");
            }
        }
        Ok(targets)
    }
}

/// Independent replay evidence for a sequence-to-graph alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoaVerificationReport {
    recomputed_score: Score,
    aligned_operations: usize,
    graph_only_operations: usize,
    query_only_operations: usize,
}

impl PoaVerificationReport {
    /// Returns the score recomputed without using DP buffers.
    #[must_use]
    pub const fn recomputed_score(self) -> Score {
        self.recomputed_score
    }

    /// Returns the number of homologous diagonal operations.
    #[must_use]
    pub const fn aligned_operations(self) -> usize {
        self.aligned_operations
    }

    /// Returns the number of graph-only gap operations.
    #[must_use]
    pub const fn graph_only_operations(self) -> usize {
        self.graph_only_operations
    }

    /// Returns the number of query-only gap operations.
    #[must_use]
    pub const fn query_only_operations(self) -> usize {
        self.query_only_operations
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
enum TraceKind {
    Stop = 0,
    Diagonal = 1,
    GraphOnly = 2,
    QueryOnly = 3,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(transparent)]
struct Trace(usize);

impl Trace {
    const KIND_SHIFT: u32 = usize::BITS - 2;
    const PREDECESSOR_MASK: usize = usize::MAX >> 2;
    const STOP: Self = Self(0);

    const fn diagonal(predecessor: Option<usize>) -> Self {
        Self::with_predecessor(TraceKind::Diagonal, predecessor)
    }

    const fn graph_only(predecessor: Option<usize>) -> Self {
        Self::with_predecessor(TraceKind::GraphOnly, predecessor)
    }

    const fn query_only() -> Self {
        Self((TraceKind::QueryOnly as usize) << Self::KIND_SHIFT)
    }

    const fn with_predecessor(kind: TraceKind, predecessor: Option<usize>) -> Self {
        let encoded = match predecessor {
            Some(row) => row + 1,
            None => 0,
        };
        Self(((kind as usize) << Self::KIND_SHIFT) | encoded)
    }

    const fn kind(self) -> TraceKind {
        match self.0 >> Self::KIND_SHIFT {
            0 => TraceKind::Stop,
            1 => TraceKind::Diagonal,
            2 => TraceKind::GraphOnly,
            _ => TraceKind::QueryOnly,
        }
    }

    const fn predecessor(self) -> Option<usize> {
        let encoded = self.0 & Self::PREDECESSOR_MASK;
        if encoded == 0 {
            None
        } else {
            Some(encoded - 1)
        }
    }
}

#[derive(Clone, Copy)]
struct Endpoint {
    score: Score,
    row: usize,
    column: usize,
}

/// Scalar exact local linear-gap alignment of a sequence to a POA DAG.
#[derive(Clone, Debug)]
pub struct LocalPoa<S> {
    substitution: S,
    gap: LinearGap,
}

impl<S> LocalPoa<S> {
    /// Creates an aligner with explicit substitution and linear-gap scoring.
    #[must_use]
    pub const fn new(substitution: S, gap: LinearGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured non-positive linear gap score.
    #[must_use]
    pub const fn gap(&self) -> LinearGap {
        self.gap
    }

    /// Computes exact algorithm-owned element-buffer requirements.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError::DimensionOverflow`] if a dimension, cell count, or
    /// byte count cannot fit in `usize`.
    pub fn memory_requirement(
        graph_nodes: usize,
        query_len: usize,
    ) -> Result<PoaAlignmentMemoryEstimate, PoaError> {
        if graph_nodes > Trace::PREDECESSOR_MASK {
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
        let cells = graph_nodes
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
        let score_matrix_bytes = checked_bytes::<Score>(cells, graph_nodes, query_len)?;
        let traceback_matrix_bytes = checked_bytes::<Trace>(cells, graph_nodes, query_len)?;
        let one_node_index = checked_bytes::<usize>(graph_nodes, graph_nodes, query_len)?;
        let node_index_bytes =
            one_node_index
                .checked_mul(2)
                .ok_or(PoaError::DimensionOverflow {
                    graph_nodes,
                    query_len,
                })?;
        let maximum_operation_bytes =
            checked_bytes::<PoaOperation>(maximum_operations, graph_nodes, query_len)?;
        let maximum_aligned_pair_bytes =
            checked_bytes::<PoaAlignedPair>(graph_nodes.min(query_len), graph_nodes, query_len)?;
        let peak_buffer_bytes = score_matrix_bytes
            .checked_add(traceback_matrix_bytes)
            .and_then(|bytes| bytes.checked_add(node_index_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_operation_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_aligned_pair_bytes))
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        Ok(PoaAlignmentMemoryEstimate {
            cells,
            score_matrix_bytes,
            traceback_matrix_bytes,
            node_index_bytes,
            maximum_operation_bytes,
            maximum_aligned_pair_bytes,
            peak_buffer_bytes,
        })
    }

    /// Aligns one sequence locally against every traversal of a validated DAG.
    ///
    /// Ties follow profile `LEE-2002-POA-local-linear-v1`: diagonals before
    /// graph-only moves, stable predecessor ID within either move, query-only
    /// last, and the first strict endpoint maximum in stable topological
    /// row-major order.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for an invalid graph, dimension or positive score
    /// overflow, allocation failure, invalid traceback, or failed independent
    /// replay.
    #[allow(clippy::too_many_lines)]
    pub fn align<A>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
    ) -> Result<SequenceGraphAlignment, PoaError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let (node_order, _) = graph.validate_with_orders()?;
        self.align_with_order(graph, query, &node_order)
    }

    #[allow(clippy::too_many_lines)]
    pub(crate) fn align_with_order<A>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        order: &[PoaNodeId],
    ) -> Result<SequenceGraphAlignment, PoaError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        let graph_nodes = graph.node_count();
        let query_len = query.len();
        let memory = Self::memory_requirement(graph_nodes, query_len)?;
        if graph_nodes == 0 || query_len == 0 {
            return Ok(SequenceGraphAlignment::canonical_empty(memory));
        }
        let columns = query_len + 1;
        let mut node_to_row = try_zeroed(graph_nodes, "POA node-to-row map")?;
        for (row, node_id) in order.iter().copied().enumerate() {
            node_to_row[node_id.index()] = row;
        }
        let mut scores = try_filled(memory.cells, 0, "POA score matrix")?;
        let mut traceback = try_filled(memory.cells, Trace::STOP, "POA traceback matrix")?;
        let gap_score = self.gap.score();
        let maximum_operations = graph_nodes + query_len;
        let endpoint =
            if every_local_score_fits::<A, S>(&self.substitution, gap_score, maximum_operations) {
                fill_matrix::<A, S, true>(
                    graph,
                    query,
                    &self.substitution,
                    gap_score,
                    order,
                    &node_to_row,
                    &mut scores,
                    &mut traceback,
                    columns,
                )?
            } else {
                fill_matrix::<A, S, false>(
                    graph,
                    query,
                    &self.substitution,
                    gap_score,
                    order,
                    &node_to_row,
                    &mut scores,
                    &mut traceback,
                    columns,
                )?
            };

        if endpoint.score == 0 {
            return Ok(SequenceGraphAlignment::canonical_empty(memory));
        }
        let mut operations = Vec::new();
        operations
            .try_reserve_exact(maximum_operations)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "POA traceback operations",
                elements: maximum_operations,
            })?;
        let mut row = endpoint.row;
        let mut column = endpoint.column;
        loop {
            let node_id = order[row];
            let trace = traceback.get(row * columns + column).copied().ok_or(
                PoaError::TracebackInvariant {
                    node: node_id,
                    column,
                },
            )?;
            match trace.kind() {
                TraceKind::Diagonal if column != 0 => {
                    operations.push(PoaOperation::Aligned {
                        node: node_id,
                        query_index: column - 1,
                    });
                    column -= 1;
                    if let Some(predecessor) = trace.predecessor() {
                        row = predecessor;
                    } else {
                        break;
                    }
                }
                TraceKind::GraphOnly => {
                    operations.push(PoaOperation::GraphOnly { node: node_id });
                    if let Some(predecessor) = trace.predecessor() {
                        row = predecessor;
                    } else {
                        break;
                    }
                }
                TraceKind::QueryOnly if column != 0 => {
                    operations.push(PoaOperation::QueryOnly {
                        query_index: column - 1,
                    });
                    column -= 1;
                }
                TraceKind::Stop => break,
                _ => {
                    return Err(PoaError::TracebackInvariant {
                        node: node_id,
                        column,
                    });
                }
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
                buffer: "POA aligned pairs",
                elements: aligned_pair_count,
            })?;
        aligned_pairs.extend(operations.iter().filter_map(|operation| match operation {
            PoaOperation::Aligned { node, query_index } => {
                Some(PoaAlignedPair::new(*node, *query_index))
            }
            PoaOperation::GraphOnly { .. } | PoaOperation::QueryOnly { .. } => None,
        }));
        let query_start = operations
            .iter()
            .find_map(consumed_query_index)
            .ok_or_else(|| PoaError::VerificationFailed {
                message: "positive local alignment consumes no query residue".to_owned(),
            })?;
        let query_end = operations
            .iter()
            .rev()
            .find_map(consumed_query_index)
            .and_then(|index| index.checked_add(1))
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes,
                query_len,
            })?;
        let alignment = SequenceGraphAlignment {
            score: endpoint.score,
            query_range: query_start..query_end,
            operations: operations.into_boxed_slice(),
            aligned_pairs: aligned_pairs.into_boxed_slice(),
            memory,
        };
        verify_sequence_graph_alignment_validated(
            &alignment,
            graph,
            query,
            &self.substitution,
            self.gap,
        )?;
        Ok(alignment)
    }
}

/// Replays and verifies a sequence-to-graph result independently of DP buffers.
///
/// # Errors
///
/// Returns [`PoaError`] if an operation is illegal, graph/query continuity is
/// broken, a local prefix is non-positive, or the recomputed score differs.
#[allow(clippy::too_many_lines)]
pub fn verify_sequence_graph_alignment<A, S>(
    alignment: &SequenceGraphAlignment,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap: LinearGap,
) -> Result<PoaVerificationReport, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    graph.validate()?;
    verify_sequence_graph_alignment_validated(alignment, graph, query, substitution, gap)
}

#[allow(clippy::too_many_lines)]
pub(crate) fn verify_sequence_graph_alignment_validated<A, S>(
    alignment: &SequenceGraphAlignment,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap: LinearGap,
) -> Result<PoaVerificationReport, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    if alignment.score == 0 {
        if alignment.query_range != (0..0)
            || !alignment.operations.is_empty()
            || !alignment.aligned_pairs.is_empty()
        {
            return verification("zero-score result is not the canonical empty alignment");
        }
        return Ok(PoaVerificationReport {
            recomputed_score: 0,
            aligned_operations: 0,
            graph_only_operations: 0,
            query_only_operations: 0,
        });
    }
    if alignment.score < 0 || alignment.operations.is_empty() {
        return verification("a non-empty local alignment must have a positive score");
    }
    if alignment.query_range.start >= alignment.query_range.end
        || alignment.query_range.end > query.len()
    {
        return verification("local query range is empty or outside the query");
    }

    let mut score = 0_i64;
    let mut previous_graph_node = None;
    let mut next_query_index = alignment.query_range.start;
    let mut aligned_pair_index = 0_usize;
    let mut aligned_operations = 0_usize;
    let mut graph_only_operations = 0_usize;
    let mut query_only_operations = 0_usize;
    for operation in &alignment.operations {
        let increment = match *operation {
            PoaOperation::Aligned { node, query_index } => {
                verify_graph_continuity(graph, previous_graph_node, node)?;
                if query_index != next_query_index {
                    return verification("aligned operation breaks consecutive query consumption");
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
                next_query_index =
                    next_query_index
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "replayed query position",
                        })?;
                previous_graph_node = Some(node);
                aligned_operations =
                    aligned_operations
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "aligned operation count",
                        })?;
                if alignment.aligned_pairs.get(aligned_pair_index)
                    != Some(&PoaAlignedPair::new(node, query_index))
                {
                    return verification(
                        "reported aligned pairs disagree with traceback operations",
                    );
                }
                aligned_pair_index =
                    aligned_pair_index
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "replayed aligned-pair index",
                        })?;
                substitution.score(graph_symbol, query_symbol)
            }
            PoaOperation::GraphOnly { node } => {
                verify_graph_continuity(graph, previous_graph_node, node)?;
                previous_graph_node = Some(node);
                graph_only_operations =
                    graph_only_operations
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "graph-only operation count",
                        })?;
                gap.score()
            }
            PoaOperation::QueryOnly { query_index } => {
                if previous_graph_node.is_none() {
                    return verification("query-only operation precedes every graph operation");
                }
                if query_index != next_query_index || query.get(query_index).is_none() {
                    return verification(
                        "query-only operation breaks consecutive query consumption",
                    );
                }
                next_query_index =
                    next_query_index
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "replayed query position",
                        })?;
                query_only_operations =
                    query_only_operations
                        .checked_add(1)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "query-only operation count",
                        })?;
                gap.score()
            }
        };
        score = score
            .checked_add(increment)
            .ok_or_else(|| PoaError::VerificationFailed {
                message: "replayed score overflows the score type".to_owned(),
            })?;
        if score <= 0 {
            return verification("a replayed local traceback prefix is non-positive");
        }
    }
    if next_query_index != alignment.query_range.end {
        return verification("traceback query consumption disagrees with its reported range");
    }
    if aligned_pair_index != alignment.aligned_pairs.len() {
        return verification("reported aligned pairs disagree with traceback operations");
    }
    if score != alignment.score {
        return verification("recomputed traceback score differs from the reported optimum");
    }
    Ok(PoaVerificationReport {
        recomputed_score: score,
        aligned_operations,
        graph_only_operations,
        query_only_operations,
    })
}

fn verify_graph_continuity<A: Alphabet>(
    graph: &PoaGraph<A>,
    previous: Option<PoaNodeId>,
    current: PoaNodeId,
) -> Result<(), PoaError> {
    if graph.node(current).is_none() {
        return Err(PoaError::UnknownNode { node: current });
    }
    if previous.is_some_and(|previous| !graph.has_edge(previous, current)) {
        return verification("consecutive graph-consuming operations do not follow an edge");
    }
    Ok(())
}

fn every_local_score_fits<A, S>(
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
    let maximum_increment = maximum_substitution.max(gap_score).max(0);
    i128::from(maximum_increment) * maximum_operations <= i128::from(Score::MAX)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_matrix<A, S, const BOUNDED: bool>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    substitution: &S,
    gap_score: Score,
    order: &[PoaNodeId],
    node_to_row: &[usize],
    scores: &mut [Score],
    traceback: &mut [Trace],
    columns: usize,
) -> Result<Endpoint, PoaError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut endpoint = Endpoint {
        score: 0,
        row: 0,
        column: 0,
    };
    for (row, node_id) in order.iter().copied().enumerate() {
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        for (query_index, query_symbol) in query.iter().enumerate() {
            let column = query_index + 1;
            let cell = row * columns + column;
            let mut best_score = 0;
            let mut best_trace = Trace::STOP;
            let substitution_score = substitution.score(node.symbol(), query_symbol);

            if node.predecessors().is_empty() {
                consider_candidate::<BOUNDED>(
                    &mut best_score,
                    &mut best_trace,
                    0,
                    substitution_score,
                    Trace::diagonal(None),
                    node_id,
                    column,
                )?;
            } else {
                for predecessor in node.predecessors() {
                    let predecessor_row = node_to_row[predecessor.index()];
                    let base = scores[predecessor_row * columns + query_index];
                    consider_candidate::<BOUNDED>(
                        &mut best_score,
                        &mut best_trace,
                        base,
                        substitution_score,
                        Trace::diagonal(Some(predecessor_row)),
                        node_id,
                        column,
                    )?;
                }
            }

            if node.predecessors().is_empty() {
                consider_candidate::<BOUNDED>(
                    &mut best_score,
                    &mut best_trace,
                    0,
                    gap_score,
                    Trace::graph_only(None),
                    node_id,
                    column,
                )?;
            } else {
                for predecessor in node.predecessors() {
                    let predecessor_row = node_to_row[predecessor.index()];
                    let base = scores[predecessor_row * columns + column];
                    consider_candidate::<BOUNDED>(
                        &mut best_score,
                        &mut best_trace,
                        base,
                        gap_score,
                        Trace::graph_only(Some(predecessor_row)),
                        node_id,
                        column,
                    )?;
                }
            }

            consider_candidate::<BOUNDED>(
                &mut best_score,
                &mut best_trace,
                scores[cell - 1],
                gap_score,
                Trace::query_only(),
                node_id,
                column,
            )?;
            scores[cell] = best_score;
            traceback[cell] = best_trace;
            if best_score > endpoint.score {
                endpoint = Endpoint {
                    score: best_score,
                    row,
                    column,
                };
            }
        }
    }
    Ok(endpoint)
}

#[inline]
fn consider_candidate<const BOUNDED: bool>(
    best_score: &mut Score,
    best_trace: &mut Trace,
    base: Score,
    increment: Score,
    candidate_trace: Trace,
    node: PoaNodeId,
    column: usize,
) -> Result<(), PoaError> {
    if BOUNDED {
        let candidate = base + increment;
        if candidate > *best_score {
            *best_score = candidate;
            *best_trace = candidate_trace;
        }
        return Ok(());
    }
    match base.checked_add(increment) {
        Some(candidate) if candidate > *best_score => {
            *best_score = candidate;
            *best_trace = candidate_trace;
        }
        Some(_) => {}
        None if increment < 0 => {}
        None => return Err(PoaError::ScoreOverflow { node, column }),
    }
    Ok(())
}

fn consumed_query_index(operation: &PoaOperation) -> Option<usize> {
    match operation {
        PoaOperation::Aligned { query_index, .. } | PoaOperation::QueryOnly { query_index } => {
            Some(*query_index)
        }
        PoaOperation::GraphOnly { .. } => None,
    }
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

fn try_zeroed(elements: usize, buffer: &'static str) -> Result<Vec<usize>, PoaError> {
    try_filled(elements, 0, buffer)
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

fn verification<T>(message: impl Into<String>) -> Result<T, PoaError> {
    Err(PoaError::VerificationFailed {
        message: message.into(),
    })
}

#[cfg(test)]
mod tests {
    use align_core::{
        Dna, DnaSequence, LinearGap, MatchMismatch, Score, SubstitutionScorer, Symbol,
    };
    use align_msa::SequenceRecord;

    use super::{LocalPoa, verify_sequence_graph_alignment};
    use crate::{PoaError, PoaGraph, PoaNodeId, PoaOperation};

    #[derive(Clone, Copy)]
    struct UnboundedMatchMismatch;

    impl SubstitutionScorer<Dna> for UnboundedMatchMismatch {
        fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
            if first == second { 2 } else { -1 }
        }
    }

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn graph(sequence: &str) -> PoaGraph<Dna> {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("reference", sequence), &vec![None; sequence.len()])
            .unwrap();
        graph
    }

    #[test]
    fn linear_graph_reduces_to_hand_computed_smith_waterman() {
        let graph = graph("ACGT");
        let query = DnaSequence::from_ascii("TTACGAA").unwrap();
        let aligner = LocalPoa::new(MatchMismatch::new(2, -1), LinearGap::new(-2).unwrap());
        let alignment = aligner.align(&graph, &query).unwrap();
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.query_range(), 2..5);
        assert_eq!(
            alignment.operations(),
            &[
                PoaOperation::Aligned {
                    node: PoaNodeId::new(0),
                    query_index: 2,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(1),
                    query_index: 3,
                },
                PoaOperation::Aligned {
                    node: PoaNodeId::new(2),
                    query_index: 4,
                },
            ]
        );
        assert_eq!(
            verify_sequence_graph_alignment(
                &alignment,
                &graph,
                &query,
                &MatchMismatch::new(2, -1),
                LinearGap::new(-2).unwrap(),
            )
            .unwrap()
            .recomputed_score(),
            6
        );
    }

    #[test]
    fn no_positive_hit_returns_the_canonical_empty_mapping() {
        let graph = graph("AAAA");
        let query = DnaSequence::from_ascii("CCCC").unwrap();
        let alignment = LocalPoa::new(MatchMismatch::new(2, -3), LinearGap::new(-2).unwrap())
            .align(&graph, &query)
            .unwrap();
        assert_eq!(alignment.score(), 0);
        assert_eq!(alignment.query_range(), 0..0);
        assert!(alignment.operations().is_empty());
    }

    #[test]
    fn traceback_contains_an_internal_graph_only_move() {
        let graph = graph("ACGT");
        let query = DnaSequence::from_ascii("AGT").unwrap();
        let alignment = LocalPoa::new(MatchMismatch::new(2, -3), LinearGap::new(-1).unwrap())
            .align(&graph, &query)
            .unwrap();
        assert_eq!(alignment.score(), 5);
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
    fn traceback_contains_an_internal_query_only_move() {
        let graph = graph("ACGT");
        let query = DnaSequence::from_ascii("ACAGT").unwrap();
        let alignment = LocalPoa::new(MatchMismatch::new(2, -3), LinearGap::new(-1).unwrap())
            .align(&graph, &query)
            .unwrap();
        assert_eq!(alignment.score(), 7);
        assert!(matches!(
            alignment.operations()[2],
            PoaOperation::QueryOnly { query_index: 2 }
        ));
    }

    #[test]
    fn endpoint_and_predecessor_ties_use_stable_node_ids() {
        let linear = graph("AA");
        let query = DnaSequence::from_ascii("A").unwrap();
        let scoring = MatchMismatch::new(2, -3);
        let gap = LinearGap::new(-1).unwrap();
        let endpoint_tie = LocalPoa::new(scoring, gap).align(&linear, &query).unwrap();
        assert_eq!(
            endpoint_tie.operations(),
            &[PoaOperation::Aligned {
                node: PoaNodeId::new(0),
                query_index: 0,
            }]
        );

        let mut bubble = graph("ACG");
        bubble
            .append_record(
                &record("allele", "ATG"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    Some(PoaNodeId::new(2)),
                ],
            )
            .unwrap();
        bubble.validate().unwrap();
        let query = DnaSequence::from_ascii("AG").unwrap();
        let predecessor_tie = LocalPoa::new(scoring, gap).align(&bubble, &query).unwrap();
        assert_eq!(predecessor_tie.score(), 3);
        assert_eq!(
            predecessor_tie.operations(),
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
            ]
        );
    }

    #[test]
    fn memory_model_is_checked_and_matches_element_layout() {
        let memory = LocalPoa::<MatchMismatch>::memory_requirement(3, 4).unwrap();
        assert_eq!(memory.cells(), 15);
        assert_eq!(memory.score_matrix_bytes(), 15 * size_of::<i64>());
        assert_eq!(memory.traceback_matrix_bytes(), 15 * size_of::<usize>());
        assert_eq!(memory.node_index_bytes(), 6 * size_of::<usize>());
        assert_eq!(
            memory.maximum_aligned_pair_bytes(),
            3 * size_of::<super::PoaAlignedPair>()
        );
        assert!(LocalPoa::<MatchMismatch>::memory_requirement(1, usize::MAX).is_err());
    }

    #[test]
    fn checked_arithmetic_reports_an_exact_positive_score_overflow() {
        let graph = graph("AA");
        let query = DnaSequence::from_ascii("AA").unwrap();
        let error = LocalPoa::new(
            MatchMismatch::new(Score::MAX, -1),
            LinearGap::new(0).unwrap(),
        )
        .align(&graph, &query)
        .unwrap_err();
        assert_eq!(
            error,
            PoaError::ScoreOverflow {
                node: PoaNodeId::new(1),
                column: 2,
            }
        );
    }

    #[test]
    fn scorers_without_declared_bounds_use_the_checked_path() {
        let graph = graph("ACGT");
        let query = DnaSequence::from_ascii("TTACGAA").unwrap();
        let alignment = LocalPoa::new(UnboundedMatchMismatch, LinearGap::new(-2).unwrap())
            .align(&graph, &query)
            .unwrap();
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.query_range(), 2..5);
    }
}
