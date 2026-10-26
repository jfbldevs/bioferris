use crate::MsaError;

/// One directed edge below a rooted guide-tree node.
#[derive(Clone, Debug, PartialEq)]
pub struct GuideTreeEdge {
    length: f64,
    child: Box<GuideTreeNode>,
}

impl GuideTreeEdge {
    /// Creates an edge with a finite, possibly negative, inferred branch length.
    ///
    /// Negative NJ limbs are retained as observations; algorithms may clamp them
    /// when deriving quantities that require nonnegative contributions.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] if `length` is NaN or infinite.
    pub fn new(length: f64, child: GuideTreeNode) -> Result<Self, MsaError> {
        if !length.is_finite() {
            return Err(MsaError::NonFiniteBranchLength { length });
        }
        Ok(Self {
            length,
            child: Box::new(child),
        })
    }

    /// Returns the signed inferred branch length.
    #[must_use]
    pub const fn length(&self) -> f64 {
        self.length
    }

    /// Returns the child node.
    #[must_use]
    pub const fn child(&self) -> &GuideTreeNode {
        &self.child
    }
}

/// A node in a rooted binary guide tree.
#[derive(Clone, Debug, PartialEq)]
pub enum GuideTreeNode {
    /// One input sequence.
    Leaf {
        /// Zero-based input index.
        sequence_index: usize,
    },
    /// A deterministic binary merge.
    Branch {
        /// Child containing the smaller minimum input index.
        left: GuideTreeEdge,
        /// Other child.
        right: GuideTreeEdge,
    },
}

impl GuideTreeNode {
    /// Creates a leaf node.
    #[must_use]
    pub const fn leaf(sequence_index: usize) -> Self {
        Self::Leaf { sequence_index }
    }

    /// Creates an ordered branch node.
    #[must_use]
    pub const fn branch(left: GuideTreeEdge, right: GuideTreeEdge) -> Self {
        Self::Branch { left, right }
    }

    /// Returns the minimum input index below this node.
    #[must_use]
    pub fn minimum_leaf_index(&self) -> usize {
        match self {
            Self::Leaf { sequence_index } => *sequence_index,
            Self::Branch { left, right } => left
                .child()
                .minimum_leaf_index()
                .min(right.child().minimum_leaf_index()),
        }
    }
}

/// A validated rooted binary guide tree containing each input exactly once.
#[derive(Clone, Debug, PartialEq)]
pub struct GuideTree {
    root: GuideTreeNode,
    leaf_count: usize,
}

impl GuideTree {
    /// Validates a rooted guide tree.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] unless every leaf in `0..leaf_count` occurs exactly
    /// once. A guide tree must contain at least one leaf.
    pub fn try_new(root: GuideTreeNode, leaf_count: usize) -> Result<Self, MsaError> {
        if leaf_count == 0 {
            return Err(MsaError::EmptyCollection);
        }
        let mut seen = vec![false; leaf_count];
        validate_node(&root, leaf_count, &mut seen)?;
        if let Some(index) = seen.iter().position(|present| !present) {
            return Err(MsaError::MissingLeaf { index });
        }
        Ok(Self { root, leaf_count })
    }

    /// Returns the root node.
    #[must_use]
    pub const fn root(&self) -> &GuideTreeNode {
        &self.root
    }

    /// Returns the number of leaves.
    #[must_use]
    pub const fn leaf_count(&self) -> usize {
        self.leaf_count
    }

    /// Returns leaf indices in deterministic postorder.
    #[must_use]
    pub fn postorder_leaves(&self) -> Vec<usize> {
        let mut leaves = Vec::with_capacity(self.leaf_count);
        collect_postorder(&self.root, &mut leaves);
        leaves
    }
}

fn validate_node(
    node: &GuideTreeNode,
    leaf_count: usize,
    seen: &mut [bool],
) -> Result<(), MsaError> {
    match node {
        GuideTreeNode::Leaf { sequence_index } => {
            if *sequence_index >= leaf_count {
                return Err(MsaError::LeafIndexOutOfRange {
                    index: *sequence_index,
                    leaf_count,
                });
            }
            if seen[*sequence_index] {
                return Err(MsaError::DuplicateLeaf {
                    index: *sequence_index,
                });
            }
            seen[*sequence_index] = true;
        }
        GuideTreeNode::Branch { left, right } => {
            validate_node(left.child(), leaf_count, seen)?;
            validate_node(right.child(), leaf_count, seen)?;
        }
    }
    Ok(())
}

fn collect_postorder(node: &GuideTreeNode, leaves: &mut Vec<usize>) {
    match node {
        GuideTreeNode::Leaf { sequence_index } => leaves.push(*sequence_index),
        GuideTreeNode::Branch { left, right } => {
            collect_postorder(left.child(), leaves);
            collect_postorder(right.child(), leaves);
        }
    }
}
