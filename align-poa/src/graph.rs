use core::cmp::Reverse;
use core::mem::size_of;
use std::collections::BinaryHeap;

use align_core::{Alphabet, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

use crate::PoaError;

/// Stable identifier of one residue node in a partial-order graph.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct PoaNodeId(usize);

impl PoaNodeId {
    pub(crate) const fn new(index: usize) -> Self {
        Self(index)
    }

    /// Returns the stable zero-based arena index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

/// Stable identifier of one homologous alignment class.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct PoaClassId(usize);

impl PoaClassId {
    const fn new(index: usize) -> Self {
        Self(index)
    }

    /// Returns the stable zero-based arena index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

/// Stable identifier of one named source path.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct PoaPathId(usize);

impl PoaPathId {
    const fn new(index: usize) -> Self {
        Self(index)
    }

    /// Returns the stable zero-based arena index.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

/// One source-sequence residue represented by a graph node.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PoaOccurrence {
    path: PoaPathId,
    source_position: usize,
}

impl PoaOccurrence {
    /// Returns the source path containing this residue.
    #[must_use]
    pub const fn path(self) -> PoaPathId {
        self.path
    }

    /// Returns the zero-based residue position in the ungapped source.
    #[must_use]
    pub const fn source_position(self) -> usize {
        self.source_position
    }
}

/// One typed residue node in a partial-order graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoaNode<A: Alphabet> {
    id: PoaNodeId,
    symbol: Symbol<A>,
    class: PoaClassId,
    occurrences: Vec<PoaOccurrence>,
    predecessors: Vec<PoaNodeId>,
    successors: Vec<PoaNodeId>,
}

impl<A: Alphabet> PoaNode<A> {
    /// Returns this node's stable identifier.
    #[must_use]
    pub const fn id(&self) -> PoaNodeId {
        self.id
    }

    /// Returns the biological symbol represented by this node.
    #[must_use]
    pub const fn symbol(&self) -> Symbol<A> {
        self.symbol
    }

    /// Returns the homologous alignment class containing this allele.
    #[must_use]
    pub const fn class(&self) -> PoaClassId {
        self.class
    }

    /// Returns all source occurrences in stable path-position order.
    #[must_use]
    pub fn occurrences(&self) -> &[PoaOccurrence] {
        self.occurrences.as_slice()
    }

    /// Returns sorted incoming neighbors.
    #[must_use]
    pub fn predecessors(&self) -> &[PoaNodeId] {
        self.predecessors.as_slice()
    }

    /// Returns sorted outgoing neighbors.
    #[must_use]
    pub fn successors(&self) -> &[PoaNodeId] {
        self.successors.as_slice()
    }
}

/// One homologous alignment class containing mutually alternative residue nodes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoaAlignmentClass {
    id: PoaClassId,
    nodes: Vec<PoaNodeId>,
}

impl PoaAlignmentClass {
    /// Returns this class's stable identifier.
    #[must_use]
    pub const fn id(&self) -> PoaClassId {
        self.id
    }

    /// Returns allele nodes in stable identifier order.
    #[must_use]
    pub fn nodes(&self) -> &[PoaNodeId] {
        self.nodes.as_slice()
    }
}

/// One supported directed edge between consecutive source residues.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PoaEdge {
    source: PoaNodeId,
    target: PoaNodeId,
    support: usize,
}

impl PoaEdge {
    /// Returns the edge source.
    #[must_use]
    pub const fn source(self) -> PoaNodeId {
        self.source
    }

    /// Returns the edge target.
    #[must_use]
    pub const fn target(self) -> PoaNodeId {
        self.target
    }

    /// Returns the number of stored source paths traversing this edge directly.
    #[must_use]
    pub const fn support(self) -> usize {
        self.support
    }
}

/// One complete named source path through the partial-order graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoaPath {
    id: PoaPathId,
    identifier: String,
    nodes: Vec<PoaNodeId>,
    source_symbols: Vec<u8>,
}

impl PoaPath {
    /// Returns this path's stable identifier.
    #[must_use]
    pub const fn id(&self) -> PoaPathId {
        self.id
    }

    /// Returns the original sequence identifier.
    #[must_use]
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    /// Returns graph nodes in exact source-residue order.
    #[must_use]
    pub fn nodes(&self) -> &[PoaNodeId] {
        self.nodes.as_slice()
    }

    /// Returns compact alphabet indices in exact source-residue order.
    pub(crate) fn source_symbols(&self) -> &[u8] {
        self.source_symbols.as_slice()
    }
}

/// A validated residue-level partial-order multiple-sequence alignment graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoaGraph<A: Alphabet> {
    nodes: Vec<PoaNode<A>>,
    classes: Vec<PoaAlignmentClass>,
    edges: Vec<PoaEdge>,
    paths: Vec<PoaPath>,
}

impl<A: Alphabet> Default for PoaGraph<A> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: Alphabet> PoaGraph<A> {
    /// Creates an empty, valid graph.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            nodes: Vec::new(),
            classes: Vec::new(),
            edges: Vec::new(),
            paths: Vec::new(),
        }
    }

    /// Returns all residue nodes in stable identifier order.
    #[must_use]
    pub fn nodes(&self) -> &[PoaNode<A>] {
        self.nodes.as_slice()
    }

    /// Returns one node by stable identifier.
    #[must_use]
    pub fn node(&self, id: PoaNodeId) -> Option<&PoaNode<A>> {
        self.nodes.get(id.index())
    }

    /// Returns all homologous alignment classes in stable identifier order.
    #[must_use]
    pub fn classes(&self) -> &[PoaAlignmentClass] {
        self.classes.as_slice()
    }

    /// Returns one class by stable identifier.
    #[must_use]
    pub fn alignment_class(&self, id: PoaClassId) -> Option<&PoaAlignmentClass> {
        self.classes.get(id.index())
    }

    /// Returns all supported edges sorted by source and then target identifier.
    #[must_use]
    pub fn edges(&self) -> &[PoaEdge] {
        self.edges.as_slice()
    }

    /// Returns all source paths in stable input order.
    #[must_use]
    pub fn paths(&self) -> &[PoaPath] {
        self.paths.as_slice()
    }

    /// Returns one source path by stable identifier.
    #[must_use]
    pub fn path(&self, id: PoaPathId) -> Option<&PoaPath> {
        self.paths.get(id.index())
    }

    /// Returns the number of residue nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the number of alignment classes.
    #[must_use]
    pub fn class_count(&self) -> usize {
        self.classes.len()
    }

    /// Returns the number of supported directed edges.
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Returns the number of stored source paths.
    #[must_use]
    pub fn path_count(&self) -> usize {
        self.paths.len()
    }

    /// Returns whether a directed edge exists.
    #[must_use]
    pub fn has_edge(&self, source: PoaNodeId, target: PoaNodeId) -> bool {
        self.edge_position(source, target).is_ok()
    }

    /// Returns edge support, if the directed edge exists.
    #[must_use]
    pub fn edge_support(&self, source: PoaNodeId, target: PoaNodeId) -> Option<usize> {
        self.edge_position(source, target)
            .ok()
            .map(|index| self.edges[index].support)
    }

    /// Computes the canonical stable topological node order.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] if allocation fails or the graph contains a cycle.
    pub fn stable_topological_order(&self) -> Result<Vec<PoaNodeId>, PoaError> {
        self.node_topological_order_unchecked()
    }

    /// Computes the canonical stable topological alignment-class order.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] if allocation fails or the class quotient contains a
    /// cycle.
    pub fn stable_class_topological_order(&self) -> Result<Vec<PoaClassId>, PoaError> {
        self.class_topological_order_unchecked()
    }

    /// Validates every graph, class, edge, occurrence, and path invariant.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError::InvalidGraph`] at the first deterministic invariant
    /// violation, or an allocation/size error while constructing audit buffers.
    pub fn validate(&self) -> Result<(), PoaError> {
        self.validate_with_orders().map(|_| ())
    }

    #[allow(clippy::too_many_lines)]
    pub(crate) fn validate_with_orders(
        &self,
    ) -> Result<(Vec<PoaNodeId>, Vec<PoaClassId>), PoaError> {
        for (index, node) in self.nodes.iter().enumerate() {
            if node.id != PoaNodeId::new(index) {
                return invalid_graph(format!("node arena entry {index} has a mismatched ID"));
            }
            if node.class.index() >= self.classes.len() {
                return invalid_graph(format!(
                    "node {index} refers to missing class {}",
                    node.class.index()
                ));
            }
            if node.occurrences.is_empty() {
                return invalid_graph(format!("node {index} has no source occurrence"));
            }
            if !strictly_sorted(&node.predecessors) {
                return invalid_graph(format!("node {index} predecessors are not strictly sorted"));
            }
            if !strictly_sorted(&node.successors) {
                return invalid_graph(format!("node {index} successors are not strictly sorted"));
            }
            if !strictly_sorted(&node.occurrences) {
                return invalid_graph(format!("node {index} occurrences are not strictly sorted"));
            }
        }

        for (index, class) in self.classes.iter().enumerate() {
            if class.id != PoaClassId::new(index) {
                return invalid_graph(format!("class arena entry {index} has a mismatched ID"));
            }
            if class.nodes.is_empty() {
                return invalid_graph(format!("alignment class {index} is empty"));
            }
            if !strictly_sorted(&class.nodes) {
                return invalid_graph(format!(
                    "alignment class {index} nodes are not strictly sorted"
                ));
            }
            for (allele_index, node_id) in class.nodes.iter().copied().enumerate() {
                let Some(node) = self.node(node_id) else {
                    return invalid_graph(format!(
                        "alignment class {index} refers to missing node {}",
                        node_id.index()
                    ));
                };
                if node.class != class.id {
                    return invalid_graph(format!(
                        "alignment class {index} contains node {} assigned to class {}",
                        node_id.index(),
                        node.class.index()
                    ));
                }
                if class.nodes[..allele_index]
                    .iter()
                    .any(|previous| self.nodes[previous.index()].symbol == node.symbol)
                {
                    return invalid_graph(format!(
                        "alignment class {index} contains duplicate symbol {}",
                        node.symbol
                    ));
                }
            }
        }

        let mut expected_predecessors = empty_nested(self.nodes.len(), "predecessor audit")?;
        let mut expected_successors = empty_nested(self.nodes.len(), "successor audit")?;
        let mut previous_edge = None;
        for edge in &self.edges {
            if edge.source.index() >= self.nodes.len() || edge.target.index() >= self.nodes.len() {
                return invalid_graph("an edge refers to a missing node".to_owned());
            }
            let key = (edge.source, edge.target);
            if previous_edge.is_some_and(|previous| previous >= key) {
                return invalid_graph("edges are not strictly sorted and unique".to_owned());
            }
            previous_edge = Some(key);
            if edge.support == 0 {
                return invalid_graph(format!(
                    "edge {} -> {} has zero support",
                    edge.source.index(),
                    edge.target.index()
                ));
            }
            if self.nodes[edge.source.index()].class == self.nodes[edge.target.index()].class {
                return invalid_graph(format!(
                    "edge {} -> {} stays inside one alignment class",
                    edge.source.index(),
                    edge.target.index()
                ));
            }
            expected_successors[edge.source.index()]
                .try_reserve(1)
                .map_err(|_| allocation("successor audit", 1))?;
            expected_predecessors[edge.target.index()]
                .try_reserve(1)
                .map_err(|_| allocation("predecessor audit", 1))?;
            expected_successors[edge.source.index()].push(edge.target);
            expected_predecessors[edge.target.index()].push(edge.source);
        }
        for (index, node) in self.nodes.iter().enumerate() {
            if node.predecessors != expected_predecessors[index] {
                return invalid_graph(format!("node {index} predecessor index is inconsistent"));
            }
            if node.successors != expected_successors[index] {
                return invalid_graph(format!("node {index} successor index is inconsistent"));
            }
        }

        let mut path_identifiers = Vec::new();
        reserve_exact(
            &mut path_identifiers,
            self.paths.len(),
            "source identifier audit",
        )?;
        for (path_index, path) in self.paths.iter().enumerate() {
            if path.id != PoaPathId::new(path_index) {
                return invalid_graph(format!("path arena entry {path_index} has a mismatched ID"));
            }
            if path.identifier.trim().is_empty()
                || path.identifier.as_bytes().iter().any(u8::is_ascii_control)
            {
                return invalid_graph(format!("path {path_index} has an invalid identifier"));
            }
            path_identifiers.push((path.identifier.as_str(), path_index));
        }
        path_identifiers.sort_unstable();
        if let Some(pair) = path_identifiers
            .windows(2)
            .find(|pair| pair[0].0 == pair[1].0)
        {
            return invalid_graph(format!(
                "path {} repeats identifier {:?}",
                pair[1].1, pair[1].0
            ));
        }

        let mut class_seen = try_zeroed(self.classes.len(), "path-class audit")?;
        let mut observed_edge_support = try_zeroed(self.edges.len(), "edge-support audit")?;
        let mut total_path_residues = 0_usize;
        for (path_index, path) in self.paths.iter().enumerate() {
            if path.source_symbols.len() != path.nodes.len() {
                return invalid_graph(format!(
                    "path {path_index} has {} cached source symbols for {} nodes",
                    path.source_symbols.len(),
                    path.nodes.len()
                ));
            }
            let marker = path_index
                .checked_add(1)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "path validation marker",
                })?;
            for (position, node_id) in path.nodes.iter().copied().enumerate() {
                let Some(node) = self.node(node_id) else {
                    return invalid_graph(format!(
                        "path {path_index} refers to missing node {}",
                        node_id.index()
                    ));
                };
                if usize::from(path.source_symbols[position]) != node.symbol.alphabet_index() {
                    return invalid_graph(format!(
                        "path {path_index} cached source symbol at position {position} disagrees with node {}",
                        node_id.index()
                    ));
                }
                let class_index = node.class.index();
                if class_seen[class_index] == marker {
                    return invalid_graph(format!(
                        "path {path_index} visits alignment class {class_index} twice"
                    ));
                }
                class_seen[class_index] = marker;
                let occurrence = PoaOccurrence {
                    path: path.id,
                    source_position: position,
                };
                if node.occurrences.binary_search(&occurrence).is_err() {
                    return invalid_graph(format!(
                        "path {path_index} position {position} is missing from node occurrences"
                    ));
                }
            }
            for pair in path.nodes.windows(2) {
                let edge_index =
                    self.edge_position(pair[0], pair[1])
                        .map_err(|_| PoaError::InvalidGraph {
                            message: format!(
                                "path {path_index} traverses missing edge {} -> {}",
                                pair[0].index(),
                                pair[1].index()
                            ),
                        })?;
                observed_edge_support[edge_index] = observed_edge_support[edge_index]
                    .checked_add(1)
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "observed edge support",
                    })?;
            }
            total_path_residues = total_path_residues.checked_add(path.nodes.len()).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "source residue count",
                },
            )?;
        }

        let mut total_occurrences = 0_usize;
        for node in &self.nodes {
            total_occurrences = total_occurrences
                .checked_add(node.occurrences.len())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "node occurrence count",
                })?;
            for occurrence in &node.occurrences {
                let Some(path) = self.path(occurrence.path) else {
                    return Err(PoaError::UnknownPath {
                        path: occurrence.path,
                    });
                };
                if path.nodes.get(occurrence.source_position) != Some(&node.id) {
                    return invalid_graph(format!(
                        "node {} occurrence ({},{}) disagrees with its source path",
                        node.id.index(),
                        occurrence.path.index(),
                        occurrence.source_position
                    ));
                }
            }
        }
        if total_occurrences != total_path_residues {
            return invalid_graph(format!(
                "{total_occurrences} node occurrences do not equal {total_path_residues} path residues"
            ));
        }

        for (edge, observed_support) in self.edges.iter().zip(observed_edge_support) {
            if observed_support != edge.support {
                return invalid_graph(format!(
                    "edge {} -> {} support is {}, expected {}",
                    edge.source.index(),
                    edge.target.index(),
                    edge.support,
                    observed_support
                ));
            }
        }

        let node_order = self.node_topological_order_unchecked()?;
        let class_order = self.class_topological_order_unchecked()?;
        Ok((node_order, class_order))
    }

    /// Serializes alignment classes into deterministic rectangular columns.
    ///
    /// The graph remains authoritative: the chosen relative order of incomparable
    /// classes has no biological meaning.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] if graph validation or rectangular MSA construction
    /// fails.
    pub fn to_multiple_alignment(&self) -> Result<MultipleAlignment<A>, PoaError> {
        let (_, class_order) = self.validate_with_orders()?;
        self.to_multiple_alignment_with_class_order(&class_order)
    }

    pub(crate) fn to_multiple_alignment_with_class_order(
        &self,
        order: &[PoaClassId],
    ) -> Result<MultipleAlignment<A>, PoaError> {
        let mut class_to_column = try_zeroed(self.classes.len(), "class-to-column map")?;
        for (column, class_id) in order.iter().copied().enumerate() {
            class_to_column[class_id.index()] = column;
        }
        let mut rows = Vec::new();
        reserve_exact(&mut rows, self.paths.len(), "rectangular POA rows")?;
        for path in &self.paths {
            let mut symbols = Vec::new();
            reserve_exact(&mut symbols, order.len(), "rectangular POA row")?;
            symbols.resize(order.len(), None);
            for node_id in &path.nodes {
                let node = &self.nodes[node_id.index()];
                let column = class_to_column[node.class.index()];
                if symbols[column].replace(node.symbol).is_some() {
                    return invalid_graph(format!(
                        "path {} occupies class column {column} twice",
                        path.id.index()
                    ));
                }
            }
            rows.push(AlignedRow::new(
                checked_string_clone(&path.identifier, "rectangular POA identifier")?,
                symbols,
            )?);
        }
        MultipleAlignment::try_new(rows).map_err(PoaError::from)
    }

    /// Returns modeled bytes occupied by graph-owned element and identifier data.
    ///
    /// Allocator metadata, capacity rounding, and the `PoaGraph` value itself are
    /// excluded.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError::ArithmeticOverflow`] if the byte total cannot fit in
    /// `usize`.
    pub fn modeled_storage_bytes(&self) -> Result<usize, PoaError> {
        let mut bytes = checked_product(self.nodes.len(), size_of::<PoaNode<A>>(), "node arena")?;
        bytes = checked_sum(
            bytes,
            checked_product(
                self.classes.len(),
                size_of::<PoaAlignmentClass>(),
                "class arena",
            )?,
            "graph storage",
        )?;
        bytes = checked_sum(
            bytes,
            checked_product(self.edges.len(), size_of::<PoaEdge>(), "edge arena")?,
            "graph storage",
        )?;
        bytes = checked_sum(
            bytes,
            checked_product(self.paths.len(), size_of::<PoaPath>(), "path arena")?,
            "graph storage",
        )?;
        for node in &self.nodes {
            bytes = checked_sum(
                bytes,
                checked_product(
                    node.occurrences.len(),
                    size_of::<PoaOccurrence>(),
                    "node occurrences",
                )?,
                "graph storage",
            )?;
            bytes = checked_sum(
                bytes,
                checked_product(
                    node.predecessors.len(),
                    size_of::<PoaNodeId>(),
                    "node predecessors",
                )?,
                "graph storage",
            )?;
            bytes = checked_sum(
                bytes,
                checked_product(
                    node.successors.len(),
                    size_of::<PoaNodeId>(),
                    "node successors",
                )?,
                "graph storage",
            )?;
        }
        for class in &self.classes {
            bytes = checked_sum(
                bytes,
                checked_product(class.nodes.len(), size_of::<PoaNodeId>(), "class alleles")?,
                "graph storage",
            )?;
        }
        for path in &self.paths {
            bytes = checked_sum(bytes, path.identifier.len(), "graph storage")?;
            bytes = checked_sum(bytes, path.source_symbols.len(), "graph storage")?;
            bytes = checked_sum(
                bytes,
                checked_product(path.nodes.len(), size_of::<PoaNodeId>(), "path nodes")?,
                "graph storage",
            )?;
        }
        Ok(bytes)
    }

    pub(crate) fn append_record(
        &mut self,
        record: &SequenceRecord<A>,
        aligned_targets: &[Option<PoaNodeId>],
    ) -> Result<PoaPathId, PoaError> {
        if aligned_targets.len() != record.sequence().len() {
            return invalid_graph(format!(
                "integration mapping has length {}, expected {}",
                aligned_targets.len(),
                record.sequence().len()
            ));
        }
        self.nodes
            .try_reserve(aligned_targets.len())
            .map_err(|_| allocation("graph nodes", aligned_targets.len()))?;
        self.classes
            .try_reserve(aligned_targets.len())
            .map_err(|_| allocation("alignment classes", aligned_targets.len()))?;
        self.paths
            .try_reserve(1)
            .map_err(|_| allocation("source paths", 1))?;
        self.edges
            .try_reserve(aligned_targets.len().saturating_sub(1))
            .map_err(|_| allocation("graph edges", aligned_targets.len().saturating_sub(1)))?;

        let path_id = PoaPathId::new(self.paths.len());
        let mut path_nodes = Vec::new();
        reserve_exact(&mut path_nodes, aligned_targets.len(), "source path nodes")?;
        let mut source_symbols = Vec::new();
        reserve_exact(
            &mut source_symbols,
            record.sequence().len(),
            "source path symbols",
        )?;
        source_symbols.extend_from_slice(record.sequence().as_indices());
        for (position, (symbol, target)) in record
            .sequence()
            .iter()
            .zip(aligned_targets.iter().copied())
            .enumerate()
        {
            let node_id = if let Some(target) = target {
                let class_id = self
                    .node(target)
                    .ok_or(PoaError::UnknownNode { node: target })?
                    .class;
                if let Some(existing) = self.classes[class_id.index()]
                    .nodes
                    .iter()
                    .copied()
                    .find(|candidate| self.nodes[candidate.index()].symbol == symbol)
                {
                    existing
                } else {
                    self.add_node(symbol, class_id)?
                }
            } else {
                let class_id = self.add_class();
                self.add_node(symbol, class_id)?
            };
            let occurrence = PoaOccurrence {
                path: path_id,
                source_position: position,
            };
            self.nodes[node_id.index()]
                .occurrences
                .try_reserve(1)
                .map_err(|_| allocation("node occurrences", 1))?;
            self.nodes[node_id.index()].occurrences.push(occurrence);
            path_nodes.push(node_id);
        }
        for pair in path_nodes.windows(2) {
            self.add_supported_edge(pair[0], pair[1])?;
        }
        self.paths.push(PoaPath {
            id: path_id,
            identifier: checked_string_clone(record.identifier(), "source path identifier")?,
            nodes: path_nodes,
            source_symbols,
        });
        Ok(path_id)
    }

    #[cfg(test)]
    pub(crate) fn into_appended_record(
        mut self,
        record: &SequenceRecord<A>,
        aligned_targets: &[Option<PoaNodeId>],
    ) -> Result<(Self, PoaPathId, Vec<PoaNodeId>, Vec<PoaClassId>), PoaError> {
        let path = self.append_record(record, aligned_targets)?;
        let (node_order, class_order) = self.validate_with_orders()?;
        Ok((self, path, node_order, class_order))
    }

    /// Appends independently verified integration targets and recomputes only
    /// the orders needed by the next iterative alignment.
    ///
    /// The caller must have validated the current graph and replayed the
    /// sequence-to-graph alignment that produced `aligned_targets`. Debug
    /// builds audit the completed graph exhaustively; optimized builds retain
    /// the canonical class-order and source-recovery checks without rescanning
    /// construction invariants that these checked mutations preserve.
    pub(crate) fn into_appended_record_trusted(
        mut self,
        record: &SequenceRecord<A>,
        aligned_targets: &[Option<PoaNodeId>],
    ) -> Result<(Self, PoaPathId, Vec<PoaNodeId>), PoaError> {
        let path = self.append_record(record, aligned_targets)?;
        let node_order = self.node_topological_order_unchecked()?;
        Ok((self, path, node_order))
    }

    fn add_class(&mut self) -> PoaClassId {
        let id = PoaClassId::new(self.classes.len());
        self.classes.push(PoaAlignmentClass {
            id,
            nodes: Vec::new(),
        });
        id
    }

    fn add_node(&mut self, symbol: Symbol<A>, class: PoaClassId) -> Result<PoaNodeId, PoaError> {
        let class_index = class.index();
        if class_index >= self.classes.len() {
            return invalid_graph(format!("cannot add a node to missing class {class_index}"));
        }
        let id = PoaNodeId::new(self.nodes.len());
        self.classes[class_index]
            .nodes
            .try_reserve(1)
            .map_err(|_| allocation("class allele nodes", 1))?;
        self.classes[class_index].nodes.push(id);
        self.nodes.push(PoaNode {
            id,
            symbol,
            class,
            occurrences: Vec::new(),
            predecessors: Vec::new(),
            successors: Vec::new(),
        });
        Ok(id)
    }

    fn add_supported_edge(&mut self, source: PoaNodeId, target: PoaNodeId) -> Result<(), PoaError> {
        if source.index() >= self.nodes.len() {
            return Err(PoaError::UnknownNode { node: source });
        }
        if target.index() >= self.nodes.len() {
            return Err(PoaError::UnknownNode { node: target });
        }
        match self.edge_position(source, target) {
            Ok(index) => {
                self.edges[index].support = self.edges[index].support.checked_add(1).ok_or(
                    PoaError::ArithmeticOverflow {
                        context: "edge support",
                    },
                )?;
            }
            Err(index) => {
                self.nodes[source.index()]
                    .successors
                    .try_reserve(1)
                    .map_err(|_| allocation("node successors", 1))?;
                self.nodes[target.index()]
                    .predecessors
                    .try_reserve(1)
                    .map_err(|_| allocation("node predecessors", 1))?;
                let successor_index = self.nodes[source.index()]
                    .successors
                    .binary_search(&target)
                    .unwrap_or_else(core::convert::identity);
                self.nodes[source.index()]
                    .successors
                    .insert(successor_index, target);
                let predecessor_index = self.nodes[target.index()]
                    .predecessors
                    .binary_search(&source)
                    .unwrap_or_else(core::convert::identity);
                self.nodes[target.index()]
                    .predecessors
                    .insert(predecessor_index, source);
                self.edges.insert(
                    index,
                    PoaEdge {
                        source,
                        target,
                        support: 1,
                    },
                );
            }
        }
        Ok(())
    }

    fn edge_position(&self, source: PoaNodeId, target: PoaNodeId) -> Result<usize, usize> {
        self.edges
            .binary_search_by_key(&(source, target), |edge| (edge.source, edge.target))
    }

    fn node_topological_order_unchecked(&self) -> Result<Vec<PoaNodeId>, PoaError> {
        let mut indegrees = Vec::new();
        reserve_exact(&mut indegrees, self.nodes.len(), "node indegrees")?;
        indegrees.extend(self.nodes.iter().map(|node| node.predecessors.len()));
        let mut ready = BinaryHeap::new();
        ready
            .try_reserve(self.nodes.len())
            .map_err(|_| allocation("node topological queue", self.nodes.len()))?;
        for (index, indegree) in indegrees.iter().copied().enumerate() {
            if indegree == 0 {
                ready.push(Reverse(PoaNodeId::new(index)));
            }
        }
        let mut order = Vec::new();
        reserve_exact(&mut order, self.nodes.len(), "node topological order")?;
        while let Some(Reverse(node_id)) = ready.pop() {
            order.push(node_id);
            for successor in &self.nodes[node_id.index()].successors {
                let indegree = indegrees
                    .get_mut(successor.index())
                    .ok_or(PoaError::UnknownNode { node: *successor })?;
                *indegree = indegree
                    .checked_sub(1)
                    .ok_or_else(|| PoaError::InvalidGraph {
                        message: format!(
                            "node {} topological indegree underflowed",
                            successor.index()
                        ),
                    })?;
                if *indegree == 0 {
                    ready.push(Reverse(*successor));
                }
            }
        }
        if order.len() != self.nodes.len() {
            return invalid_graph("node graph contains a directed cycle".to_owned());
        }
        Ok(order)
    }

    fn class_topological_order_unchecked(&self) -> Result<Vec<PoaClassId>, PoaError> {
        let mut quotient_edges = Vec::new();
        reserve_exact(
            &mut quotient_edges,
            self.edges.len(),
            "class quotient edges",
        )?;
        quotient_edges.extend(self.edges.iter().map(|edge| {
            (
                self.nodes[edge.source.index()].class,
                self.nodes[edge.target.index()].class,
            )
        }));
        quotient_edges.sort_unstable();
        quotient_edges.dedup();

        let mut successors = empty_nested(self.classes.len(), "class successors")?;
        let mut indegrees = try_zeroed(self.classes.len(), "class indegrees")?;
        for (source, target) in quotient_edges {
            if source == target {
                return invalid_graph(format!(
                    "alignment class {} has a self edge",
                    source.index()
                ));
            }
            successors[source.index()]
                .try_reserve(1)
                .map_err(|_| allocation("class successors", 1))?;
            successors[source.index()].push(target);
            indegrees[target.index()] =
                indegrees[target.index()]
                    .checked_add(1)
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "alignment-class indegree",
                    })?;
        }

        let mut ready = BinaryHeap::new();
        ready
            .try_reserve(self.classes.len())
            .map_err(|_| allocation("class topological queue", self.classes.len()))?;
        for (index, indegree) in indegrees.iter().copied().enumerate() {
            if indegree == 0 {
                ready.push(Reverse(PoaClassId::new(index)));
            }
        }
        let mut order = Vec::new();
        reserve_exact(&mut order, self.classes.len(), "class topological order")?;
        while let Some(Reverse(class_id)) = ready.pop() {
            order.push(class_id);
            for successor in &successors[class_id.index()] {
                indegrees[successor.index()] = indegrees[successor.index()]
                    .checked_sub(1)
                    .ok_or_else(|| PoaError::InvalidGraph {
                        message: format!(
                            "class {} topological indegree underflowed",
                            successor.index()
                        ),
                    })?;
                if indegrees[successor.index()] == 0 {
                    ready.push(Reverse(*successor));
                }
            }
        }
        if order.len() != self.classes.len() {
            return invalid_graph("alignment-class quotient contains a directed cycle".to_owned());
        }
        Ok(order)
    }
}

fn invalid_graph<T>(message: String) -> Result<T, PoaError> {
    Err(PoaError::InvalidGraph { message })
}

fn allocation(buffer: &'static str, elements: usize) -> PoaError {
    PoaError::AllocationFailed { buffer, elements }
}

fn reserve_exact<T>(
    values: &mut Vec<T>,
    elements: usize,
    buffer: &'static str,
) -> Result<(), PoaError> {
    values
        .try_reserve_exact(elements)
        .map_err(|_| allocation(buffer, elements))
}

fn checked_string_clone(source: &str, buffer: &'static str) -> Result<String, PoaError> {
    let mut clone = String::new();
    clone
        .try_reserve_exact(source.len())
        .map_err(|_| allocation(buffer, source.len()))?;
    clone.push_str(source);
    Ok(clone)
}

fn try_zeroed(elements: usize, buffer: &'static str) -> Result<Vec<usize>, PoaError> {
    let mut values = Vec::new();
    reserve_exact(&mut values, elements, buffer)?;
    values.resize(elements, 0);
    Ok(values)
}

fn empty_nested<T>(elements: usize, buffer: &'static str) -> Result<Vec<Vec<T>>, PoaError> {
    let mut values = Vec::new();
    reserve_exact(&mut values, elements, buffer)?;
    values.resize_with(elements, Vec::new);
    Ok(values)
}

fn strictly_sorted<T: Ord>(values: &[T]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn checked_product(left: usize, right: usize, context: &'static str) -> Result<usize, PoaError> {
    left.checked_mul(right)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_sum(left: usize, right: usize, context: &'static str) -> Result<usize, PoaError> {
    left.checked_add(right)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use super::{PoaGraph, PoaNodeId};

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn one_path_builds_a_valid_linear_graph_and_serializes_exactly() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "ACG"), &[None, None, None])
            .unwrap();
        graph.validate().unwrap();
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.class_count(), 3);
        assert_eq!(graph.edge_count(), 2);
        assert_eq!(
            graph.edge_support(PoaNodeId::new(0), PoaNodeId::new(1)),
            Some(1)
        );
        assert_eq!(
            graph
                .to_multiple_alignment()
                .unwrap()
                .row(0)
                .unwrap()
                .to_ascii_with_gaps(),
            "ACG"
        );
    }

    #[test]
    fn identical_target_nodes_are_fused_and_edge_support_is_counted() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph
            .append_record(
                &record("second", "AC"),
                &[Some(PoaNodeId::new(0)), Some(PoaNodeId::new(1))],
            )
            .unwrap();
        graph.validate().unwrap();
        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.path_count(), 2);
        assert_eq!(
            graph.edge_support(PoaNodeId::new(0), PoaNodeId::new(1)),
            Some(2)
        );
        assert_eq!(
            graph.node(PoaNodeId::new(0)).unwrap().occurrences().len(),
            2
        );
    }

    #[test]
    fn aggressive_fusion_reuses_an_existing_same_symbol_allele() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "ACG"), &[None, None, None])
            .unwrap();
        let targets = [
            Some(PoaNodeId::new(0)),
            Some(PoaNodeId::new(1)),
            Some(PoaNodeId::new(2)),
        ];
        graph
            .append_record(&record("second", "ATG"), &targets)
            .unwrap();
        assert_eq!(graph.node_count(), 4);
        graph
            .append_record(&record("third", "ATG"), &targets)
            .unwrap();
        graph.validate().unwrap();
        assert_eq!(graph.node_count(), 4);
        let thymine = graph.classes()[1]
            .nodes()
            .iter()
            .copied()
            .find(|node| graph.node(*node).unwrap().symbol().as_ascii() == b'T')
            .unwrap();
        assert_eq!(graph.node(thymine).unwrap().occurrences().len(), 2);
        assert_eq!(graph.paths()[1].nodes()[1], thymine);
        assert_eq!(graph.paths()[2].nodes()[1], thymine);
    }

    #[test]
    fn validation_rejects_a_cycle_and_bad_support() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.nodes[0].predecessors.push(PoaNodeId::new(1));
        graph.nodes[1].successors.push(PoaNodeId::new(0));
        graph.edges.push(super::PoaEdge {
            source: PoaNodeId::new(1),
            target: PoaNodeId::new(0),
            support: 1,
        });
        assert!(graph.validate().is_err());

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.edges[0].support = 2;
        assert!(graph.validate().is_err());
    }

    #[test]
    fn validation_rejects_corrupt_ids_occurrences_paths_and_alleles() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.nodes[0].id = PoaNodeId::new(9);
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("mismatched ID")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.nodes[0].occurrences[0].source_position = 7;
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("missing from node occurrences")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.paths[0].nodes[0] = PoaNodeId::new(7);
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("missing node")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.paths[0].source_symbols[0] = graph.paths[0].source_symbols[0].wrapping_add(1);
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("cached source symbol")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "ACG"), &[None, None, None])
            .unwrap();
        graph
            .append_record(
                &record("second", "ATG"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    Some(PoaNodeId::new(2)),
                ],
            )
            .unwrap();
        graph.nodes[3].symbol = graph.nodes[1].symbol;
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("duplicate symbol")
        );
    }

    #[test]
    fn validation_rejects_same_class_unsupported_and_out_of_range_edges() {
        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.nodes[1].class = graph.nodes[0].class;
        graph.classes[0].nodes.push(PoaNodeId::new(1));
        graph.classes.pop();
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("stays inside one alignment class")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "ACG"), &[None, None, None])
            .unwrap();
        graph
            .add_supported_edge(PoaNodeId::new(0), PoaNodeId::new(2))
            .unwrap();
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("support is 1, expected 0")
        );

        let mut graph = PoaGraph::new();
        graph
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        graph.nodes[1].successors.push(PoaNodeId::new(9));
        graph.edges.push(super::PoaEdge {
            source: PoaNodeId::new(1),
            target: PoaNodeId::new(9),
            support: 1,
        });
        assert!(
            graph
                .validate()
                .unwrap_err()
                .to_string()
                .contains("missing node")
        );
    }

    #[test]
    fn validation_rejects_repeated_path_classes_and_class_quotient_cycles() {
        let mut repeated = PoaGraph::new();
        repeated
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        repeated
            .append_record(
                &record("repeated", "ACA"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    Some(PoaNodeId::new(0)),
                ],
            )
            .unwrap();
        assert!(
            repeated
                .validate()
                .unwrap_err()
                .to_string()
                .contains("visits alignment class 0 twice")
        );

        let mut quotient_cycle = PoaGraph::new();
        quotient_cycle
            .append_record(&record("first", "AC"), &[None, None])
            .unwrap();
        quotient_cycle
            .append_record(
                &record("second", "GT"),
                &[Some(PoaNodeId::new(1)), Some(PoaNodeId::new(0))],
            )
            .unwrap();
        assert!(
            quotient_cycle
                .validate()
                .unwrap_err()
                .to_string()
                .contains("class quotient contains a directed cycle")
        );
    }
}
