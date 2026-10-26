use core::cmp::{max, min};
use core::ops::Index;

use crate::{PoaError, PoaGraph, PoaNodeId};
use align_core::Alphabet;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Superbubble {
    entrance: usize,
    exit: usize,
}

impl Superbubble {
    #[cfg(test)]
    pub(crate) const fn entrance(self) -> usize {
        self.entrance
    }

    pub(crate) const fn exit(self) -> usize {
        self.exit
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct BubbleMembership {
    bubble: usize,
    min_edges_to_exit: usize,
    max_edges_to_exit: usize,
}

impl BubbleMembership {
    pub(crate) const fn bubble(self) -> usize {
        self.bubble
    }

    pub(crate) const fn min_edges_to_exit(self) -> usize {
        self.min_edges_to_exit
    }

    pub(crate) const fn max_edges_to_exit(self) -> usize {
        self.max_edges_to_exit
    }
}

#[derive(Clone, Debug)]
struct SuccessorIndex {
    offsets: Box<[usize]>,
    targets: Box<[usize]>,
}

impl SuccessorIndex {
    fn from_graph<A: Alphabet>(graph: &PoaGraph<A>) -> Result<Self, PoaError> {
        let sink = graph
            .node_count()
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: 0,
            })?;
        let cursor_count = sink.checked_add(1).ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: 0,
        })?;
        let entries = if graph.node_count() == 0 {
            1
        } else {
            graph
                .nodes()
                .iter()
                .try_fold(0_usize, |total, node| {
                    let outgoing = if node.successors().is_empty() {
                        1
                    } else {
                        node.successors().len()
                    };
                    total
                        .checked_add(outgoing)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA compact successor count",
                        })
                })?
                .checked_add(
                    graph
                        .nodes()
                        .iter()
                        .filter(|node| node.predecessors().is_empty())
                        .count(),
                )
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "POASTA compact successor count",
                })?
        };
        let offset_count = cursor_count
            .checked_add(1)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA compact successor offset count",
            })?;
        let mut offsets = Vec::new();
        try_reserve_exact(
            &mut offsets,
            offset_count,
            "POASTA compact successor offsets",
        )?;
        let mut targets = Vec::new();
        try_reserve_exact(&mut targets, entries, "POASTA compact successor targets")?;
        offsets.push(0);
        if graph.node_count() == 0 {
            targets.push(sink);
        } else {
            targets.extend(
                graph
                    .nodes()
                    .iter()
                    .filter(|node| node.predecessors().is_empty())
                    .map(|node| node.id().index() + 1),
            );
        }
        offsets.push(targets.len());
        for node in graph.nodes() {
            if node.successors().is_empty() {
                targets.push(sink);
            } else {
                targets.extend(
                    node.successors()
                        .iter()
                        .map(|successor| successor.index() + 1),
                );
            }
            offsets.push(targets.len());
        }
        offsets.push(targets.len());
        debug_assert_eq!(offsets.len(), offset_count);
        debug_assert_eq!(targets.len(), entries);
        Ok(Self {
            offsets: offsets.into_boxed_slice(),
            targets: targets.into_boxed_slice(),
        })
    }

    #[cfg(test)]
    fn from_nested(rows: Box<[Box<[usize]>]>) -> Result<Self, PoaError> {
        let entries = rows.iter().try_fold(0_usize, |total, row| {
            total
                .checked_add(row.len())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "POASTA compact successor count",
                })
        })?;
        let offset_count = rows
            .len()
            .checked_add(1)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA compact successor offset count",
            })?;
        let mut offsets = Vec::new();
        try_reserve_exact(
            &mut offsets,
            offset_count,
            "POASTA compact successor offsets",
        )?;
        let mut targets = Vec::new();
        try_reserve_exact(&mut targets, entries, "POASTA compact successor targets")?;
        offsets.push(0);
        for row in rows {
            targets.extend_from_slice(&row);
            offsets.push(targets.len());
        }
        Ok(Self {
            offsets: offsets.into_boxed_slice(),
            targets: targets.into_boxed_slice(),
        })
    }

    fn transposed(&self) -> Result<Self, PoaError> {
        let offset_count = self
            .len()
            .checked_add(1)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA transposed-adjacency offset count",
            })?;
        let mut offsets = filled(offset_count, 0_usize, "POASTA transposed-adjacency offsets")?;
        for &target in &self.targets {
            let end = target
                .checked_add(1)
                .and_then(|cursor| offsets.get_mut(cursor))
                .ok_or_else(|| {
                    internal_error("an adjacency target is outside its transposed index")
                })?;
            *end = checked_increment(*end, "POASTA transposed-adjacency degree")?;
        }
        for cursor in 1..offsets.len() {
            offsets[cursor] = checked_add(
                offsets[cursor],
                offsets[cursor - 1],
                "POASTA transposed-adjacency offset",
            )?;
        }
        let mut next = Vec::new();
        try_reserve_exact(&mut next, self.len(), "POASTA transposed-adjacency cursors")?;
        next.extend_from_slice(&offsets[..self.len()]);
        let mut targets = filled(
            self.targets.len(),
            0_usize,
            "POASTA transposed-adjacency targets",
        )?;
        for source in 0..self.len() {
            for &target in self.get(source) {
                let write = next[target];
                targets[write] = source;
                next[target] =
                    checked_increment(next[target], "POASTA transposed-adjacency write cursor")?;
            }
        }
        Ok(Self {
            offsets: offsets.into_boxed_slice(),
            targets: targets.into_boxed_slice(),
        })
    }

    fn get(&self, cursor: usize) -> &[usize] {
        let Some(end_cursor) = cursor.checked_add(1) else {
            return &[];
        };
        let (Some(&start), Some(&end)) = (self.offsets.get(cursor), self.offsets.get(end_cursor))
        else {
            return &[];
        };
        self.targets.get(start..end).unwrap_or(&[])
    }

    fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    fn modeled_payload_bytes(&self) -> Result<usize, PoaError> {
        checked_product(
            self.offsets.len().checked_add(self.targets.len()).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "POASTA compact successor payload",
                },
            )?,
            core::mem::size_of::<usize>(),
            "POASTA compact successor payload",
        )
    }
}

impl Index<usize> for SuccessorIndex {
    type Output = [usize];

    fn index(&self, cursor: usize) -> &Self::Output {
        self.get(cursor)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct PoastaGraphIndex {
    source: usize,
    sink: usize,
    successors: SuccessorIndex,
    min_edges_to_sink: Box<[usize]>,
    max_edges_to_sink: Box<[usize]>,
    bubbles: Box<[Superbubble]>,
    membership_offsets: Box<[usize]>,
    memberships: Box<[BubbleMembership]>,
    exit_offsets: Box<[usize]>,
    exit_bubbles: Box<[usize]>,
}

type SuperbubbleIndex = (
    Box<[Superbubble]>,
    Box<[usize]>,
    Box<[BubbleMembership]>,
    Box<[usize]>,
    Box<[usize]>,
);

impl PoastaGraphIndex {
    pub(crate) fn build<A: Alphabet>(
        graph: &PoaGraph<A>,
        topological_nodes: &[PoaNodeId],
        include_superbubbles: bool,
    ) -> Result<Self, PoaError> {
        if !include_superbubbles {
            return Self::build_compact(graph, topological_nodes);
        }
        let dag = AugmentedDag::build(graph, topological_nodes)?;
        let (min_edges_to_sink, max_edges_to_sink) = terminal_distances(&dag)?;
        let (bubbles, membership_offsets, memberships, exit_offsets, exit_bubbles) =
            build_superbubble_index(&dag, &min_edges_to_sink, &max_edges_to_sink)?;
        let successors = dag.successors;
        Ok(Self {
            source: dag.source,
            sink: dag.sink,
            successors,
            min_edges_to_sink: min_edges_to_sink.into_boxed_slice(),
            max_edges_to_sink: max_edges_to_sink.into_boxed_slice(),
            bubbles,
            membership_offsets,
            memberships,
            exit_offsets,
            exit_bubbles,
        })
    }

    fn build_compact<A: Alphabet>(
        graph: &PoaGraph<A>,
        topological_nodes: &[PoaNodeId],
    ) -> Result<Self, PoaError> {
        let source = 0;
        let sink = graph
            .node_count()
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: 0,
            })?;
        let successors = SuccessorIndex::from_graph(graph)?;
        let (min_edges_to_sink, max_edges_to_sink) =
            compact_terminal_distances(&successors, topological_nodes, source, sink)?;
        Ok(Self {
            source,
            sink,
            successors,
            min_edges_to_sink: min_edges_to_sink.into_boxed_slice(),
            max_edges_to_sink: max_edges_to_sink.into_boxed_slice(),
            bubbles: Box::default(),
            membership_offsets: Box::default(),
            memberships: Box::default(),
            exit_offsets: Box::default(),
            exit_bubbles: Box::default(),
        })
    }

    pub(crate) const fn source(&self) -> usize {
        self.source
    }

    pub(crate) const fn sink(&self) -> usize {
        self.sink
    }

    #[cfg(test)]
    pub(crate) fn cursor_for_node(&self, node: PoaNodeId) -> Option<usize> {
        let cursor = node.index().checked_add(1)?;
        (cursor < self.sink).then_some(cursor)
    }

    pub(crate) fn node_for_cursor(&self, cursor: usize) -> Option<PoaNodeId> {
        (cursor > self.source && cursor < self.sink).then(|| PoaNodeId::new(cursor - 1))
    }

    pub(crate) fn successors(&self, cursor: usize) -> &[usize] {
        self.successors.get(cursor)
    }

    pub(crate) fn min_remaining_nodes(&self, cursor: usize) -> Option<usize> {
        self.min_edges_to_sink
            .get(cursor)
            .and_then(|distance| distance.checked_sub(1))
    }

    pub(crate) fn max_remaining_nodes(&self, cursor: usize) -> Option<usize> {
        self.max_edges_to_sink
            .get(cursor)
            .and_then(|distance| distance.checked_sub(1))
    }

    pub(crate) fn bubbles(&self) -> &[Superbubble] {
        &self.bubbles
    }

    pub(crate) fn memberships(&self, cursor: usize) -> &[BubbleMembership] {
        indexed_slice(&self.membership_offsets, &self.memberships, cursor)
    }

    pub(crate) fn bubbles_ending_at(&self, cursor: usize) -> &[usize] {
        indexed_slice(&self.exit_offsets, &self.exit_bubbles, cursor)
    }

    pub(crate) fn membership_entries(&self) -> usize {
        self.memberships.len()
    }

    pub(crate) fn modeled_payload_bytes(&self) -> Result<usize, PoaError> {
        let cursor_count = self.successors.len();
        let index_offsets = self
            .membership_offsets
            .len()
            .checked_add(self.exit_offsets.len())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA indexed-slice offset count",
            })?;
        let mut bytes = self.successors.modeled_payload_bytes()?;
        bytes = checked_add(
            bytes,
            checked_product(
                index_offsets,
                core::mem::size_of::<usize>(),
                "POASTA indexed-slice offsets",
            )?,
            "POASTA graph-index payload",
        )?;
        bytes = checked_add(
            bytes,
            checked_product(
                cursor_count,
                core::mem::size_of::<usize>() * 2,
                "POASTA terminal distances",
            )?,
            "POASTA graph-index payload",
        )?;
        bytes = checked_add(
            bytes,
            checked_product(
                self.bubbles.len(),
                core::mem::size_of::<Superbubble>(),
                "POASTA superbubbles",
            )?,
            "POASTA graph-index payload",
        )?;
        bytes = checked_add(
            bytes,
            checked_product(
                self.memberships.len(),
                core::mem::size_of::<BubbleMembership>(),
                "POASTA bubble memberships",
            )?,
            "POASTA graph-index payload",
        )?;
        checked_add(
            bytes,
            checked_product(
                self.exit_bubbles.len(),
                core::mem::size_of::<usize>(),
                "POASTA bubble-exit index",
            )?,
            "POASTA graph-index payload",
        )
    }
}

fn indexed_slice<'a, T>(offsets: &[usize], values: &'a [T], cursor: usize) -> &'a [T] {
    let Some(end_cursor) = cursor.checked_add(1) else {
        return &[];
    };
    let (Some(&start), Some(&end)) = (offsets.get(cursor), offsets.get(end_cursor)) else {
        return &[];
    };
    values.get(start..end).unwrap_or(&[])
}

#[derive(Clone, Debug)]
struct AugmentedDag {
    source: usize,
    sink: usize,
    successors: SuccessorIndex,
    predecessors: SuccessorIndex,
    topological: Box<[usize]>,
}

impl AugmentedDag {
    fn build<A: Alphabet>(
        graph: &PoaGraph<A>,
        topological_nodes: &[PoaNodeId],
    ) -> Result<Self, PoaError> {
        let sink = graph
            .node_count()
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: 0,
            })?;
        let successors = SuccessorIndex::from_graph(graph)?;
        let cursor_count = successors.len();
        let source = 0;
        let predecessors = successors.transposed()?;

        let mut topological = Vec::new();
        try_reserve_exact(
            &mut topological,
            cursor_count,
            "POASTA augmented topological order",
        )?;
        topological.push(source);
        topological.extend(topological_nodes.iter().map(|node| node.index() + 1));
        topological.push(sink);
        if topological.len() != cursor_count {
            return internal("the augmented topological order has the wrong length");
        }

        Ok(Self {
            source,
            sink,
            successors,
            predecessors,
            topological: topological.into_boxed_slice(),
        })
    }

    fn len(&self) -> usize {
        self.successors.len()
    }
}

fn terminal_distances(dag: &AugmentedDag) -> Result<(Vec<usize>, Vec<usize>), PoaError> {
    let mut minimum = filled(dag.len(), usize::MAX, "POASTA minimum terminal distances")?;
    let mut maximum = filled(dag.len(), 0, "POASTA maximum terminal distances")?;
    minimum[dag.sink] = 0;
    for &cursor in dag.topological.iter().rev() {
        if cursor == dag.sink {
            continue;
        }
        let successors = &dag.successors[cursor];
        if successors.is_empty() {
            return internal("a non-sink augmented cursor has no successor");
        }
        let mut min_distance = usize::MAX;
        let mut max_distance = 0;
        for &successor in successors {
            let successor_min = *minimum
                .get(successor)
                .ok_or_else(|| internal_error("terminal-distance successor is out of range"))?;
            if successor_min == usize::MAX {
                return internal("terminal distances violate topological order");
            }
            min_distance = min(
                min_distance,
                checked_increment(successor_min, "minimum terminal distance")?,
            );
            max_distance = max(
                max_distance,
                checked_increment(maximum[successor], "maximum terminal distance")?,
            );
        }
        minimum[cursor] = min_distance;
        maximum[cursor] = max_distance;
    }
    Ok((minimum, maximum))
}

fn compact_terminal_distances(
    successors: &SuccessorIndex,
    topological_nodes: &[PoaNodeId],
    source: usize,
    sink: usize,
) -> Result<(Vec<usize>, Vec<usize>), PoaError> {
    let mut minimum = filled(
        successors.len(),
        usize::MAX,
        "POASTA minimum terminal distances",
    )?;
    let mut maximum = filled(successors.len(), 0, "POASTA maximum terminal distances")?;
    minimum[sink] = 0;
    for cursor in topological_nodes
        .iter()
        .rev()
        .map(|node| node.index() + 1)
        .chain(core::iter::once(source))
    {
        let targets = successors.get(cursor);
        if targets.is_empty() {
            return internal("a non-sink augmented cursor has no successor");
        }
        let mut min_distance = usize::MAX;
        let mut max_distance = 0;
        for &successor in targets {
            let successor_min = *minimum
                .get(successor)
                .ok_or_else(|| internal_error("terminal-distance successor is out of range"))?;
            if successor_min == usize::MAX {
                return internal("terminal distances violate topological order");
            }
            min_distance = min(
                min_distance,
                checked_increment(successor_min, "minimum terminal distance")?,
            );
            max_distance = max(
                max_distance,
                checked_increment(maximum[successor], "maximum terminal distance")?,
            );
        }
        minimum[cursor] = min_distance;
        maximum[cursor] = max_distance;
    }
    Ok((minimum, maximum))
}

#[allow(clippy::too_many_lines)]
fn build_superbubble_index(
    dag: &AugmentedDag,
    min_edges_to_sink: &[usize],
    max_edges_to_sink: &[usize],
) -> Result<SuperbubbleIndex, PoaError> {
    let reverse_postorder = reverse_dfs_postorder(dag)?;
    let mut position = filled(
        dag.len(),
        usize::MAX,
        "superbubble reverse-postorder inverse",
    )?;
    for (rank, &cursor) in reverse_postorder.iter().enumerate() {
        position[cursor] = rank;
    }

    let mut out_parent = filled(dag.len(), None, "superbubble OutParent")?;
    let mut out_child = filled(dag.len(), None, "superbubble OutChild")?;
    for cursor in 0..dag.len() {
        out_parent[cursor] = dag.predecessors[cursor]
            .iter()
            .map(|predecessor| position[*predecessor])
            .min();
        out_child[cursor] = dag.successors[cursor]
            .iter()
            .map(|successor| position[*successor])
            .max();
    }

    let mut detected = detect_superbubbles(&reverse_postorder, &position, &out_parent, &out_child)?;
    detected.sort_unstable_by_key(|bubble| (position[bubble.entrance], position[bubble.exit]));

    let mut membership_cursors = filled(dag.len(), 0_usize, "POASTA membership counts")?;
    let mut exit_cursors = filled(dag.len(), 0_usize, "POASTA bubble-exit counts")?;
    for bubble in detected.iter().copied() {
        let entrance_position = position[bubble.entrance];
        let exit_position = position[bubble.exit];
        if entrance_position >= exit_position {
            return internal("a detected superbubble is reversed in DFS order");
        }
        exit_cursors[bubble.exit] =
            checked_increment(exit_cursors[bubble.exit], "POASTA bubbles sharing an exit")?;
        for &cursor in reverse_postorder
            .iter()
            .take(exit_position)
            .skip(entrance_position)
        {
            membership_cursors[cursor] = checked_increment(
                membership_cursors[cursor],
                "POASTA node bubble membership count",
            )?;
        }
    }
    let membership_offsets =
        offsets_from_counts(&membership_cursors, "POASTA bubble-membership offsets")?;
    let exit_offsets = offsets_from_counts(&exit_cursors, "POASTA bubble-exit offsets")?;
    let membership_entries = *membership_offsets.last().unwrap_or(&0);
    let exit_entries = *exit_offsets.last().unwrap_or(&0);
    let mut memberships = filled(
        membership_entries,
        BubbleMembership {
            bubble: 0,
            min_edges_to_exit: 0,
            max_edges_to_exit: 0,
        },
        "POASTA bubble memberships",
    )?;
    let mut exit_bubbles = filled(exit_entries, 0_usize, "POASTA bubble exits")?;
    membership_cursors.copy_from_slice(&membership_offsets[..dag.len()]);
    exit_cursors.copy_from_slice(&exit_offsets[..dag.len()]);

    for (bubble_index, bubble) in detected.iter().copied().enumerate() {
        let exit_write = exit_cursors[bubble.exit];
        exit_bubbles[exit_write] = bubble_index;
        exit_cursors[bubble.exit] =
            checked_increment(exit_cursors[bubble.exit], "POASTA bubble-exit write cursor")?;
        let entrance_position = position[bubble.entrance];
        let exit_position = position[bubble.exit];
        for &cursor in reverse_postorder
            .iter()
            .take(exit_position)
            .skip(entrance_position)
        {
            let min_edges_to_exit = min_edges_to_sink[cursor]
                .checked_sub(min_edges_to_sink[bubble.exit])
                .ok_or_else(|| internal_error("minimum bubble distance underflowed"))?;
            let max_edges_to_exit = max_edges_to_sink[cursor]
                .checked_sub(max_edges_to_sink[bubble.exit])
                .ok_or_else(|| internal_error("maximum bubble distance underflowed"))?;
            let membership_write = membership_cursors[cursor];
            memberships[membership_write] = BubbleMembership {
                bubble: bubble_index,
                min_edges_to_exit,
                max_edges_to_exit,
            };
            membership_cursors[cursor] = checked_increment(
                membership_cursors[cursor],
                "POASTA bubble-membership write cursor",
            )?;
        }
    }

    Ok((
        detected.into_boxed_slice(),
        membership_offsets.into_boxed_slice(),
        memberships.into_boxed_slice(),
        exit_offsets.into_boxed_slice(),
        exit_bubbles.into_boxed_slice(),
    ))
}

fn reverse_dfs_postorder(dag: &AugmentedDag) -> Result<Vec<usize>, PoaError> {
    let mut visited = filled(dag.len(), false, "superbubble DFS visited flags")?;
    let mut stack = Vec::new();
    try_reserve_exact(&mut stack, dag.len(), "superbubble DFS stack")?;
    let mut postorder = Vec::new();
    try_reserve_exact(&mut postorder, dag.len(), "superbubble DFS postorder")?;

    visited[dag.source] = true;
    stack.push((dag.source, 0_usize));
    while let Some((cursor, next_child)) = stack.last_mut() {
        if let Some(&child) = dag.successors[*cursor].get(*next_child) {
            *next_child = checked_increment(*next_child, "superbubble DFS child index")?;
            if !visited[child] {
                visited[child] = true;
                stack.push((child, 0));
            }
        } else {
            postorder.push(*cursor);
            stack.pop();
        }
    }
    if postorder.len() != dag.len() || visited.iter().any(|seen| !seen) {
        return internal("the virtual source does not reach every augmented graph cursor");
    }
    postorder.reverse();
    Ok(postorder)
}

fn detect_superbubbles(
    reverse_postorder: &[usize],
    position: &[usize],
    out_parent: &[Option<usize>],
    out_child: &[Option<usize>],
) -> Result<Vec<Superbubble>, PoaError> {
    let cursor_count = reverse_postorder.len();
    let mut stack: Vec<Option<usize>> = Vec::new();
    try_reserve_exact(&mut stack, cursor_count, "superbubble exit stack")?;
    let mut out_parent_map = filled(cursor_count, None, "superbubble interval parents")?;
    let mut exit: Option<usize> = None;
    let mut detected = Vec::new();
    try_reserve_exact(&mut detected, cursor_count, "detected superbubbles")?;

    for rank in (0..cursor_count).rev() {
        let cursor = reverse_postorder[rank];
        let child = out_child[cursor];
        if child == rank.checked_add(1) {
            stack.push(exit);
            exit = reverse_postorder.get(rank + 1).copied();
        }

        while child_is_after_exit(child, exit, position) {
            let previous_exit = exit;
            exit = stack
                .pop()
                .ok_or_else(|| internal_error("superbubble candidate stack underflow"))?;
            merge_interval_parent(&mut out_parent_map, previous_exit, exit);
        }

        if interval_parent(&out_parent_map, exit) == Some(rank) {
            let reported_exit =
                exit.ok_or_else(|| internal_error("a superbubble candidate has no exit"))?;
            detected.push(Superbubble {
                entrance: cursor,
                exit: reported_exit,
            });
            let previous_exit = exit;
            exit = stack
                .pop()
                .ok_or_else(|| internal_error("superbubble completion stack underflow"))?;
            merge_interval_parent(&mut out_parent_map, previous_exit, exit);
        }

        out_parent_map[cursor] = out_parent[cursor];
        let merged = parent_min(
            interval_parent(&out_parent_map, exit),
            out_parent_map[cursor],
        );
        set_interval_parent(&mut out_parent_map, exit, merged);
    }
    Ok(detected)
}

fn child_is_after_exit(child: Option<usize>, exit: Option<usize>, position: &[usize]) -> bool {
    match (child, exit) {
        (None, _) | (_, None) => child.is_none() && exit.is_some(),
        (Some(child), Some(exit)) => child > position[exit],
    }
}

fn parent_min(left: Option<usize>, right: Option<usize>) -> Option<usize> {
    match (left, right) {
        (None, _) | (_, None) => None,
        (Some(left), Some(right)) => Some(min(left, right)),
    }
}

fn interval_parent(values: &[Option<usize>], cursor: Option<usize>) -> Option<usize> {
    cursor.and_then(|cursor| values.get(cursor).copied().flatten())
}

fn set_interval_parent(values: &mut [Option<usize>], cursor: Option<usize>, value: Option<usize>) {
    if let Some(cursor) = cursor {
        values[cursor] = value;
    }
}

fn merge_interval_parent(
    values: &mut [Option<usize>],
    previous_exit: Option<usize>,
    exit: Option<usize>,
) {
    let merged = parent_min(
        interval_parent(values, previous_exit),
        interval_parent(values, exit),
    );
    set_interval_parent(values, exit, merged);
}

fn filled<T: Clone>(elements: usize, value: T, buffer: &'static str) -> Result<Vec<T>, PoaError> {
    let mut values = Vec::new();
    try_reserve_exact(&mut values, elements, buffer)?;
    values.resize(elements, value);
    Ok(values)
}

fn offsets_from_counts(counts: &[usize], buffer: &'static str) -> Result<Vec<usize>, PoaError> {
    let entries = counts
        .len()
        .checked_add(1)
        .ok_or(PoaError::ArithmeticOverflow {
            context: "POASTA indexed-slice offset count",
        })?;
    let mut offsets = Vec::new();
    try_reserve_exact(&mut offsets, entries, buffer)?;
    offsets.push(0);
    for &count in counts {
        let next = checked_add(
            *offsets.last().unwrap_or(&0),
            count,
            "POASTA indexed-slice offset",
        )?;
        offsets.push(next);
    }
    Ok(offsets)
}

fn try_reserve_exact<T>(
    values: &mut Vec<T>,
    additional: usize,
    buffer: &'static str,
) -> Result<(), PoaError> {
    values
        .try_reserve_exact(additional)
        .map_err(|_| PoaError::AllocationFailed {
            buffer,
            elements: additional,
        })
}

fn checked_increment(value: usize, context: &'static str) -> Result<usize, PoaError> {
    value
        .checked_add(1)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_product(left: usize, right: usize, context: &'static str) -> Result<usize, PoaError> {
    left.checked_mul(right)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_add(left: usize, right: usize, context: &'static str) -> Result<usize, PoaError> {
    left.checked_add(right)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn internal<T>(message: &'static str) -> Result<T, PoaError> {
    Err(internal_error(message))
}

fn internal_error(message: &'static str) -> PoaError {
    PoaError::InvalidGraph {
        message: format!("POASTA preprocessing invariant: {message}"),
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use super::*;

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(
            identifier,
            DnaSequence::from_ascii(sequence.as_bytes()).expect("valid DNA"),
        )
        .expect("valid record")
    }

    fn linear_graph(sequence: &str) -> PoaGraph<Dna> {
        let source = record("backbone", sequence);
        let mut graph = PoaGraph::new();
        graph
            .append_record(&source, &vec![None; source.sequence().len()])
            .expect("seed graph");
        graph
    }

    fn ordered_dag_from_mask(vertex_count: usize, edge_mask: usize) -> Option<AugmentedDag> {
        let source = 0;
        let sink = vertex_count.checked_sub(1)?;
        let mut successors = vec![Vec::new(); vertex_count];
        let mut bit = 0_usize;
        for (from, targets) in successors.iter_mut().enumerate() {
            for to in (from + 1)..vertex_count {
                if edge_mask & (1_usize << bit) != 0 {
                    targets.push(to);
                }
                bit += 1;
            }
        }

        let mut reachable = vec![false; vertex_count];
        reachable[source] = true;
        for cursor in 0..vertex_count {
            if reachable[cursor] {
                for &successor in &successors[cursor] {
                    reachable[successor] = true;
                }
            }
        }
        if reachable.iter().any(|is_reachable| !is_reachable) {
            return None;
        }

        let mut reaches_sink = vec![false; vertex_count];
        reaches_sink[sink] = true;
        for cursor in (0..sink).rev() {
            reaches_sink[cursor] = successors[cursor]
                .iter()
                .any(|successor| reaches_sink[*successor]);
        }
        if reaches_sink.iter().any(|does_reach| !does_reach) {
            return None;
        }

        let successors = SuccessorIndex::from_nested(
            successors
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        )
        .ok()?;
        let predecessors = successors.transposed().ok()?;
        Some(AugmentedDag {
            source,
            sink,
            successors,
            predecessors,
            topological: (0..vertex_count).collect::<Vec<_>>().into_boxed_slice(),
        })
    }

    fn reachable_without_crossing(
        adjacency: &SuccessorIndex,
        start: usize,
        stop: usize,
    ) -> Vec<bool> {
        let mut reached = vec![false; adjacency.len()];
        let mut stack = vec![start];
        reached[start] = true;
        while let Some(cursor) = stack.pop() {
            if cursor == stop {
                continue;
            }
            for &next in adjacency[cursor].iter().rev() {
                if !reached[next] {
                    reached[next] = true;
                    stack.push(next);
                }
            }
        }
        reached
    }

    fn is_superbubbloid(dag: &AugmentedDag, entrance: usize, exit: usize) -> bool {
        if entrance == exit {
            return false;
        }
        let forward = reachable_without_crossing(&dag.successors, entrance, exit);
        if !forward[exit] {
            return false;
        }
        let reverse = reachable_without_crossing(&dag.predecessors, exit, entrance);
        forward == reverse
    }

    fn brute_force_superbubbles(dag: &AugmentedDag) -> Vec<(usize, usize)> {
        let mut bubbles = Vec::new();
        for entrance in 0..dag.len() {
            for exit in 0..dag.len() {
                if !is_superbubbloid(dag, entrance, exit) {
                    continue;
                }
                let vertices = reachable_without_crossing(&dag.successors, entrance, exit);
                let has_nearer_entrance =
                    vertices.iter().enumerate().any(|(candidate, is_inside)| {
                        *is_inside
                            && candidate != entrance
                            && is_superbubbloid(dag, candidate, exit)
                    });
                if !has_nearer_entrance {
                    bubbles.push((entrance, exit));
                }
            }
        }
        bubbles.sort_unstable();
        bubbles
    }

    fn published_superbubbles(dag: &AugmentedDag) -> Vec<(usize, usize)> {
        let reverse_postorder = reverse_dfs_postorder(dag).expect("reverse DFS postorder");
        let mut position = vec![usize::MAX; dag.len()];
        for (rank, &cursor) in reverse_postorder.iter().enumerate() {
            position[cursor] = rank;
        }
        let out_parent = (0..dag.len())
            .map(|cursor| {
                dag.predecessors[cursor]
                    .iter()
                    .map(|predecessor| position[*predecessor])
                    .min()
            })
            .collect::<Vec<_>>();
        let out_child = (0..dag.len())
            .map(|cursor| {
                dag.successors[cursor]
                    .iter()
                    .map(|successor| position[*successor])
                    .max()
            })
            .collect::<Vec<_>>();
        let mut bubbles =
            detect_superbubbles(&reverse_postorder, &position, &out_parent, &out_child)
                .expect("published detector")
                .into_iter()
                .map(|bubble| (bubble.entrance, bubble.exit))
                .collect::<Vec<_>>();
        bubbles.sort_unstable();
        bubbles
    }

    fn local_distances_to_exit(
        dag: &AugmentedDag,
        entrance: usize,
        exit: usize,
    ) -> Vec<Option<(usize, usize)>> {
        let inside = reachable_without_crossing(&dag.successors, entrance, exit);
        let mut minimum = vec![usize::MAX; dag.len()];
        let mut maximum = vec![0; dag.len()];
        minimum[exit] = 0;
        for &cursor in dag.topological.iter().rev() {
            if cursor == exit || !inside[cursor] {
                continue;
            }
            for &successor in &dag.successors[cursor] {
                if !inside[successor] {
                    continue;
                }
                assert_ne!(minimum[successor], usize::MAX);
                minimum[cursor] = min(minimum[cursor], minimum[successor] + 1);
                maximum[cursor] = max(maximum[cursor], maximum[successor] + 1);
            }
            assert_ne!(minimum[cursor], usize::MAX);
        }
        (0..dag.len())
            .map(|cursor| {
                (inside[cursor] && cursor != exit).then_some((minimum[cursor], maximum[cursor]))
            })
            .collect()
    }

    #[test]
    fn terminal_distances_include_virtual_edges_but_report_remaining_nodes() {
        let graph = linear_graph("ACG");
        let order = graph.stable_topological_order().expect("order");
        let index = PoastaGraphIndex::build(&graph, &order, false).expect("index");
        assert_eq!(index.min_remaining_nodes(index.source()), Some(3));
        assert_eq!(index.max_remaining_nodes(index.source()), Some(3));
        for (node, expected) in graph.nodes().iter().zip([2, 1, 0]) {
            let cursor = index.cursor_for_node(node.id()).expect("cursor");
            assert_eq!(index.min_remaining_nodes(cursor), Some(expected));
            assert_eq!(index.max_remaining_nodes(cursor), Some(expected));
        }
    }

    #[test]
    fn empty_graph_has_one_virtual_edge_and_zero_remaining_nodes() {
        let graph = PoaGraph::<Dna>::new();
        let order = graph.stable_topological_order().expect("order");
        let index = PoastaGraphIndex::build(&graph, &order, true).expect("index");
        assert_eq!(index.successors(index.source()), &[index.sink()]);
        assert_eq!(index.min_remaining_nodes(index.source()), Some(0));
        assert_eq!(index.max_remaining_nodes(index.source()), Some(0));
    }

    #[test]
    fn detects_branching_superbubble_and_indexes_both_alleles() {
        let mut graph = linear_graph("ACG");
        let alternate = record("alternate", "ATG");
        graph
            .append_record(
                &alternate,
                &[Some(PoaNodeId::new(0)), None, Some(PoaNodeId::new(2))],
            )
            .expect("add allele");
        let order = graph.stable_topological_order().expect("order");
        let index = PoastaGraphIndex::build(&graph, &order, true).expect("index");

        let entrance = index.cursor_for_node(PoaNodeId::new(0)).expect("entrance");
        let exit = index.cursor_for_node(PoaNodeId::new(2)).expect("exit");
        let bubble_index = index
            .bubbles()
            .iter()
            .position(|bubble| bubble.entrance() == entrance && bubble.exit() == exit)
            .expect("branch bubble");
        for node in [PoaNodeId::new(0), PoaNodeId::new(1), PoaNodeId::new(3)] {
            let cursor = index.cursor_for_node(node).expect("cursor");
            let membership = index
                .memberships(cursor)
                .iter()
                .find(|membership| membership.bubble() == bubble_index)
                .expect("membership");
            assert!(membership.min_edges_to_exit() >= 1);
            assert!(membership.max_edges_to_exit() >= membership.min_edges_to_exit());
        }
        assert!(
            index
                .memberships(exit)
                .iter()
                .all(|membership| membership.bubble() != bubble_index)
        );
    }

    #[test]
    fn retains_nested_superbubbles_and_multiple_memberships() {
        let mut graph = linear_graph("ACGTA");
        graph
            .append_record(
                &record("inner-alternate", "ACCTA"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    None,
                    Some(PoaNodeId::new(3)),
                    Some(PoaNodeId::new(4)),
                ],
            )
            .expect("inner allele");
        graph
            .append_record(
                &record("outer-alternate", "AGA"),
                &[Some(PoaNodeId::new(0)), None, Some(PoaNodeId::new(4))],
            )
            .expect("outer allele");
        let order = graph.stable_topological_order().expect("order");
        let index = PoastaGraphIndex::build(&graph, &order, true).expect("index");

        let outer_entrance = index.cursor_for_node(PoaNodeId::new(0)).unwrap();
        let outer_exit = index.cursor_for_node(PoaNodeId::new(4)).unwrap();
        let inner_entrance = index.cursor_for_node(PoaNodeId::new(1)).unwrap();
        let inner_exit = index.cursor_for_node(PoaNodeId::new(3)).unwrap();
        let outer = index
            .bubbles()
            .iter()
            .position(|bubble| bubble.entrance() == outer_entrance && bubble.exit() == outer_exit)
            .expect("outer bubble");
        let inner = index
            .bubbles()
            .iter()
            .position(|bubble| bubble.entrance() == inner_entrance && bubble.exit() == inner_exit)
            .expect("inner bubble");
        let inner_branch = index.cursor_for_node(PoaNodeId::new(2)).unwrap();
        let memberships = index.memberships(inner_branch);
        assert!(
            memberships
                .iter()
                .any(|membership| membership.bubble() == outer)
        );
        assert!(
            memberships
                .iter()
                .any(|membership| membership.bubble() == inner)
        );
    }

    #[test]
    fn published_detector_matches_the_definition_for_every_ordered_dag_through_six_vertices() {
        let mut checked = 0_usize;
        for vertex_count in 2..=6 {
            let possible_edges = vertex_count * (vertex_count - 1) / 2;
            for edge_mask in 0..(1_usize << possible_edges) {
                let Some(dag) = ordered_dag_from_mask(vertex_count, edge_mask) else {
                    continue;
                };
                let expected = brute_force_superbubbles(&dag);
                let actual = published_superbubbles(&dag);
                assert_eq!(
                    actual, expected,
                    "vertex_count={vertex_count}, edge_mask={edge_mask:#x}"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 3_481);
    }

    #[test]
    fn compact_membership_distances_match_local_oracle_for_every_ordered_dag_through_six_vertices()
    {
        let mut checked = 0_usize;
        for vertex_count in 2..=6 {
            let possible_edges = vertex_count * (vertex_count - 1) / 2;
            for edge_mask in 0..(1_usize << possible_edges) {
                let Some(dag) = ordered_dag_from_mask(vertex_count, edge_mask) else {
                    continue;
                };
                let (minimum, maximum) = terminal_distances(&dag).expect("terminal distances");
                let (bubbles, offsets, memberships, exit_offsets, exit_bubbles) =
                    build_superbubble_index(&dag, &minimum, &maximum)
                        .expect("compact superbubble index");
                for (bubble_index, bubble) in bubbles.iter().copied().enumerate() {
                    let expected = local_distances_to_exit(&dag, bubble.entrance(), bubble.exit());
                    for (cursor, expected_distance) in expected.into_iter().enumerate() {
                        let actual = indexed_slice(&offsets, &memberships, cursor)
                            .iter()
                            .find(|membership| membership.bubble() == bubble_index)
                            .map(|membership| {
                                (
                                    membership.min_edges_to_exit(),
                                    membership.max_edges_to_exit(),
                                )
                            });
                        assert_eq!(
                            actual, expected_distance,
                            "vertex_count={vertex_count}, edge_mask={edge_mask:#x}, \
                             bubble={bubble_index}, cursor={cursor}"
                        );
                    }
                    assert!(
                        indexed_slice(&exit_offsets, &exit_bubbles, bubble.exit())
                            .contains(&bubble_index)
                    );
                }
                checked += 1;
            }
        }
        assert_eq!(checked, 3_481);
    }
}
