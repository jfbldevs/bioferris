use crate::PoaNodeId;

/// One forward operation in an exact sequence-to-graph traceback.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PoaOperation {
    /// Consume one graph node and one homologously aligned query residue.
    Aligned {
        /// Stable graph node used by the dynamic program.
        node: PoaNodeId,
        /// Zero-based query residue index.
        query_index: usize,
    },
    /// Consume one graph node while emitting a gap in the query.
    GraphOnly {
        /// Stable graph node consumed by the operation.
        node: PoaNodeId,
    },
    /// Consume one query residue without consuming a real graph node.
    ///
    /// In a global traceback this can occur at the implicit source, at a real
    /// graph node, or after the selected sink.
    QueryOnly {
        /// Zero-based query residue index.
        query_index: usize,
    },
}

/// One homologous graph-node/query-residue pair from a traceback.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PoaAlignedPair {
    node: PoaNodeId,
    query_index: usize,
}

impl PoaAlignedPair {
    pub(crate) const fn new(node: PoaNodeId, query_index: usize) -> Self {
        Self { node, query_index }
    }

    /// Returns the graph node participating in this aligned pair.
    #[must_use]
    pub const fn node(self) -> PoaNodeId {
        self.node
    }

    /// Returns the zero-based query residue index.
    #[must_use]
    pub const fn query_index(self) -> usize {
        self.query_index
    }
}
