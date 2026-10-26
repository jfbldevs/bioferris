use core::cmp::Ordering;
use core::mem::size_of;
use std::collections::BinaryHeap;

use align_msa::{DistanceMatrix, GuideTree, GuideTreeEdge, GuideTreeNode};

use crate::{GuideTreePermutation, Muscle5Error};

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Muscle5Node {
    Leaf {
        sequence_index: usize,
    },
    Branch {
        left: Box<Muscle5Node>,
        right: Box<Muscle5Node>,
        left_length: f64,
        right_length: f64,
    },
}

impl Muscle5Node {
    pub(crate) fn minimum_leaf(&self) -> usize {
        match self {
            Self::Leaf { sequence_index } => *sequence_index,
            Self::Branch { left, right, .. } => left.minimum_leaf().min(right.minimum_leaf()),
        }
    }

    pub(crate) fn leaf_count(&self) -> usize {
        match self {
            Self::Leaf { .. } => 1,
            Self::Branch { left, right, .. } => left.leaf_count() + right.leaf_count(),
        }
    }

    #[cfg(test)]
    pub(crate) fn members(&self) -> Vec<usize> {
        let mut members = Vec::with_capacity(self.leaf_count());
        self.collect_members(&mut members);
        members.sort_unstable();
        members
    }

    #[cfg(test)]
    fn collect_members(&self, members: &mut Vec<usize>) {
        match self {
            Self::Leaf { sequence_index } => members.push(*sequence_index),
            Self::Branch { left, right, .. } => {
                left.collect_members(members);
                right.collect_members(members);
            }
        }
    }

    fn to_public(&self) -> Result<GuideTreeNode, Muscle5Error> {
        match self {
            Self::Leaf { sequence_index } => Ok(GuideTreeNode::leaf(*sequence_index)),
            Self::Branch {
                left,
                right,
                left_length,
                right_length,
            } => Ok(GuideTreeNode::branch(
                GuideTreeEdge::new(*left_length, left.to_public()?)?,
                GuideTreeEdge::new(*right_length, right.to_public()?)?,
            )),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Muscle5Tree {
    root: Muscle5Node,
    leaf_count: usize,
    workspace_bytes: usize,
}

impl Muscle5Tree {
    pub(crate) const fn root(&self) -> &Muscle5Node {
        &self.root
    }

    pub(crate) const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    pub(crate) fn to_public(&self) -> Result<GuideTree, Muscle5Error> {
        GuideTree::try_new(self.root.to_public()?, self.leaf_count).map_err(Muscle5Error::from)
    }

    pub(crate) fn permuted(&self, permutation: GuideTreePermutation) -> Result<Self, Muscle5Error> {
        if permutation == GuideTreePermutation::None || self.leaf_count < 3 {
            return Ok(self.clone());
        }
        let root = permute_root(&self.root, permutation)?;
        Ok(Self {
            root,
            leaf_count: self.leaf_count,
            workspace_bytes: self.workspace_bytes,
        })
    }
}

#[derive(Clone, Debug)]
struct Cluster {
    node: Muscle5Node,
    members: Box<[usize]>,
    size: usize,
    minimum_leaf: usize,
    height: f64,
}

#[derive(Clone, Copy, Debug)]
struct HeapEntry {
    distance: f64,
    first: usize,
    second: usize,
    first_minimum: usize,
    second_minimum: usize,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.distance.to_bits() == other.distance.to_bits()
            && self.first_minimum == other.first_minimum
            && self.second_minimum == other.second_minimum
            && self.first == other.first
            && self.second == other.second
    }
}

impl Eq for HeapEntry {}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.first_minimum.cmp(&self.first_minimum))
            .then_with(|| other.second_minimum.cmp(&self.second_minimum))
            .then_with(|| other.first.cmp(&self.first))
            .then_with(|| other.second.cmp(&self.second))
    }
}

pub(crate) fn expected_error_upgma(
    distances: &DistanceMatrix,
) -> Result<Muscle5Tree, Muscle5Error> {
    weighted_expected_error_upgma(distances, &vec![1; distances.size()])
}

#[allow(clippy::too_many_lines)]
pub(crate) fn weighted_expected_error_upgma(
    distances: &DistanceMatrix,
    leaf_weights: &[usize],
) -> Result<Muscle5Tree, Muscle5Error> {
    let leaf_count = distances.size();
    if leaf_count == 0 {
        return Err(Muscle5Error::InternalInvariant {
            message: "UPGMA requires at least one sequence",
        });
    }
    if leaf_count == 1 {
        if leaf_weights != [1] && (leaf_weights.len() != 1 || leaf_weights[0] == 0) {
            return Err(Muscle5Error::InvalidParameter {
                parameter: "leaf_weights",
                requirement: "contain one positive value per guide-tree leaf",
            });
        }
        return Ok(Muscle5Tree {
            root: Muscle5Node::Leaf { sequence_index: 0 },
            leaf_count,
            workspace_bytes: 0,
        });
    }
    if leaf_weights.len() != leaf_count || leaf_weights.contains(&0) {
        return Err(Muscle5Error::InvalidParameter {
            parameter: "leaf_weights",
            requirement: "contain one positive value per guide-tree leaf",
        });
    }
    let maximum_nodes = leaf_count
        .checked_mul(2)
        .and_then(|value| value.checked_sub(1))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "UPGMA node count",
        })?;
    let distance_cells =
        maximum_nodes
            .checked_mul(maximum_nodes)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "UPGMA distance workspace",
            })?;
    let mut matrix = vec![f64::INFINITY; distance_cells];
    let mut clusters = vec![None; maximum_nodes];
    let mut active = vec![false; maximum_nodes];
    for leaf in 0..leaf_count {
        clusters[leaf] = Some(Cluster {
            node: Muscle5Node::Leaf {
                sequence_index: leaf,
            },
            members: vec![leaf].into_boxed_slice(),
            size: leaf_weights[leaf],
            minimum_leaf: leaf,
            height: 0.0,
        });
        active[leaf] = true;
    }
    let mut heap = BinaryHeap::new();
    for first in 0..leaf_count {
        for second in first + 1..leaf_count {
            let distance = distances
                .get(first, second)
                .ok_or(Muscle5Error::InternalInvariant {
                    message: "expected-error matrix index missing",
                })?;
            matrix[first * maximum_nodes + second] = distance;
            matrix[second * maximum_nodes + first] = distance;
            heap.push(heap_entry(first, second, distance, &clusters));
        }
    }

    let mut next_id = leaf_count;
    let mut active_count = leaf_count;
    while active_count > 1 {
        let selected = loop {
            let candidate = heap.pop().ok_or(Muscle5Error::InternalInvariant {
                message: "UPGMA heap exhausted before the final merge",
            })?;
            if active[candidate.first] && active[candidate.second] {
                break candidate;
            }
        };
        let first_id = selected.first;
        let second_id = selected.second;
        let first = clusters[first_id]
            .take()
            .ok_or(Muscle5Error::InternalInvariant {
                message: "active UPGMA cluster is absent",
            })?;
        let second = clusters[second_id]
            .take()
            .ok_or(Muscle5Error::InternalInvariant {
                message: "active UPGMA cluster is absent",
            })?;
        active[first_id] = false;
        active[second_id] = false;
        let parent_height = selected.distance / 2.0;
        let first_size = first.size;
        let second_size = second.size;
        let (left, right) = if first.minimum_leaf <= second.minimum_leaf {
            (first, second)
        } else {
            (second, first)
        };
        let mut members = Vec::with_capacity(left.size + right.size);
        members.extend_from_slice(&left.members);
        members.extend_from_slice(&right.members);
        members.sort_unstable();
        let combined_size =
            left.size
                .checked_add(right.size)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "weighted UPGMA cluster size",
                })?;
        let node = Muscle5Node::Branch {
            left: Box::new(left.node),
            right: Box::new(right.node),
            left_length: (parent_height - left.height).max(0.0),
            right_length: (parent_height - right.height).max(0.0),
        };
        let new_cluster = Cluster {
            node,
            members: members.into_boxed_slice(),
            size: combined_size,
            minimum_leaf: left.minimum_leaf.min(right.minimum_leaf),
            height: parent_height,
        };
        clusters[next_id] = Some(new_cluster);
        active[next_id] = true;

        let first_weight =
            f64::from(
                u32::try_from(first_size).map_err(|_| Muscle5Error::ArithmeticOverflow {
                    stage: "weighted UPGMA floating cluster size",
                })?,
            );
        let second_weight = f64::from(u32::try_from(second_size).map_err(|_| {
            Muscle5Error::ArithmeticOverflow {
                stage: "weighted UPGMA floating cluster size",
            }
        })?);
        let combined_weight = f64::from(u32::try_from(combined_size).map_err(|_| {
            Muscle5Error::ArithmeticOverflow {
                stage: "weighted UPGMA floating cluster size",
            }
        })?);

        for other in 0..next_id {
            if !active[other] {
                continue;
            }
            let first_distance = matrix[first_id * maximum_nodes + other];
            let second_distance = matrix[second_id * maximum_nodes + other];
            let distance =
                (first_distance * first_weight + second_distance * second_weight) / combined_weight;
            matrix[next_id * maximum_nodes + other] = distance;
            matrix[other * maximum_nodes + next_id] = distance;
            heap.push(heap_entry(next_id, other, distance, &clusters));
        }
        next_id += 1;
        active_count -= 1;
    }
    let root_id =
        active
            .iter()
            .position(|is_active| *is_active)
            .ok_or(Muscle5Error::InternalInvariant {
                message: "UPGMA final root is absent",
            })?;
    let root = clusters[root_id]
        .take()
        .ok_or(Muscle5Error::InternalInvariant {
            message: "UPGMA final cluster is absent",
        })?
        .node;
    let workspace_bytes = matrix
        .len()
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(active.len() * size_of::<bool>()))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "UPGMA workspace bytes",
        })?;
    Ok(Muscle5Tree {
        root,
        leaf_count,
        workspace_bytes,
    })
}

fn heap_entry(
    first: usize,
    second: usize,
    distance: f64,
    clusters: &[Option<Cluster>],
) -> HeapEntry {
    let first_minimum = clusters[first]
        .as_ref()
        .expect("heap cluster exists")
        .minimum_leaf;
    let second_minimum = clusters[second]
        .as_ref()
        .expect("heap cluster exists")
        .minimum_leaf;
    if first_minimum <= second_minimum {
        HeapEntry {
            distance,
            first,
            second,
            first_minimum,
            second_minimum,
        }
    } else {
        HeapEntry {
            distance,
            first: second,
            second: first,
            first_minimum: second_minimum,
            second_minimum: first_minimum,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct GraphEdge {
    neighbor: usize,
    length: f64,
}

#[derive(Clone, Debug)]
struct GraphNode {
    leaf: Option<usize>,
    edges: Vec<GraphEdge>,
}

#[derive(Clone, Debug)]
struct TreeGraph {
    nodes: Vec<GraphNode>,
}

#[derive(Clone, Debug)]
struct Component {
    nodes: Vec<bool>,
    leaves: Vec<usize>,
}

#[derive(Clone, Debug)]
struct OrientedCut {
    inside: usize,
    outside: usize,
}

#[derive(Clone, Debug)]
struct BuiltComponent {
    node: Muscle5Node,
    stem_length: f64,
}

fn permute_root(
    root: &Muscle5Node,
    permutation: GuideTreePermutation,
) -> Result<Muscle5Node, Muscle5Error> {
    let graph = TreeGraph::from_tree(root);
    let all_nodes = vec![true; graph.nodes.len()];
    let first = select_cut(&graph, &[], &all_nodes, root.leaf_count(), 3, 2)?;
    let first_key = edge_key(first.inside, first.outside);
    let first_blocks = [first_key];
    let remainder = graph.component(first.outside, &first_blocks)?;
    let second = select_cut(
        &graph,
        &first_blocks,
        &remainder.nodes,
        remainder.leaves.len(),
        2,
        1,
    )?;
    let second_key = edge_key(second.inside, second.outside);
    let blocks = [first_key, second_key];
    let a = graph.build_component(first.inside, None, &blocks)?.node;
    let b = graph.build_component(second.inside, None, &blocks)?.node;
    let c = graph.build_component(second.outside, None, &blocks)?.node;
    let result = match permutation {
        GuideTreePermutation::None => root.clone(),
        GuideTreePermutation::Abc => join(join(a, b), c),
        GuideTreePermutation::Acb => join(join(a, c), b),
        GuideTreePermutation::Bca => join(join(b, c), a),
    };
    Ok(result)
}

impl TreeGraph {
    fn from_tree(root: &Muscle5Node) -> Self {
        let mut graph = Self { nodes: Vec::new() };
        graph.add_tree_node(root);
        graph
    }

    fn add_tree_node(&mut self, node: &Muscle5Node) -> usize {
        let identifier = self.nodes.len();
        self.nodes.push(GraphNode {
            leaf: match node {
                Muscle5Node::Leaf { sequence_index } => Some(*sequence_index),
                Muscle5Node::Branch { .. } => None,
            },
            edges: Vec::new(),
        });
        if let Muscle5Node::Branch {
            left,
            right,
            left_length,
            right_length,
        } = node
        {
            let left_id = self.add_tree_node(left);
            let right_id = self.add_tree_node(right);
            self.add_edge(identifier, left_id, *left_length);
            self.add_edge(identifier, right_id, *right_length);
        }
        identifier
    }

    fn add_edge(&mut self, first: usize, second: usize, length: f64) {
        self.nodes[first].edges.push(GraphEdge {
            neighbor: second,
            length,
        });
        self.nodes[second].edges.push(GraphEdge {
            neighbor: first,
            length,
        });
    }

    fn component(
        &self,
        start: usize,
        blocked: &[(usize, usize)],
    ) -> Result<Component, Muscle5Error> {
        let mut nodes = vec![false; self.nodes.len()];
        let mut leaves = Vec::new();
        let mut stack = vec![start];
        nodes[start] = true;
        while let Some(current) = stack.pop() {
            if let Some(leaf) = self.nodes[current].leaf {
                leaves.push(leaf);
            }
            for edge in &self.nodes[current].edges {
                if is_blocked(current, edge.neighbor, blocked) || nodes[edge.neighbor] {
                    continue;
                }
                nodes[edge.neighbor] = true;
                stack.push(edge.neighbor);
            }
        }
        if leaves.is_empty() {
            return Err(Muscle5Error::InternalInvariant {
                message: "guide-tree cut produced a leafless component",
            });
        }
        leaves.sort_unstable();
        Ok(Component { nodes, leaves })
    }

    fn build_component(
        &self,
        current: usize,
        parent: Option<usize>,
        blocked: &[(usize, usize)],
    ) -> Result<BuiltComponent, Muscle5Error> {
        if let Some(sequence_index) = self.nodes[current].leaf {
            return Ok(BuiltComponent {
                node: Muscle5Node::Leaf { sequence_index },
                stem_length: 0.0,
            });
        }
        let mut children = Vec::new();
        for edge in &self.nodes[current].edges {
            if Some(edge.neighbor) == parent || is_blocked(current, edge.neighbor, blocked) {
                continue;
            }
            let child = self.build_component(edge.neighbor, Some(current), blocked)?;
            children.push((child.node, edge.length + child.stem_length));
        }
        match children.len() {
            0 => Err(Muscle5Error::InternalInvariant {
                message: "guide-tree component ended at an internal node",
            }),
            1 => {
                let (node, stem_length) = children.pop().expect("one component child");
                Ok(BuiltComponent { node, stem_length })
            }
            2 => {
                let (first, first_length) = children.remove(0);
                let (second, second_length) = children.remove(0);
                let (left, left_length, right, right_length) =
                    if first.minimum_leaf() <= second.minimum_leaf() {
                        (first, first_length, second, second_length)
                    } else {
                        (second, second_length, first, first_length)
                    };
                Ok(BuiltComponent {
                    node: Muscle5Node::Branch {
                        left: Box::new(left),
                        right: Box::new(right),
                        left_length,
                        right_length,
                    },
                    stem_length: 0.0,
                })
            }
            _ => Err(Muscle5Error::InternalInvariant {
                message: "binary guide-tree component has more than two children",
            }),
        }
    }
}

fn select_cut(
    graph: &TreeGraph,
    blocked: &[(usize, usize)],
    allowed_nodes: &[bool],
    total_leaves: usize,
    balance_factor: usize,
    minimum_complement: usize,
) -> Result<OrientedCut, Muscle5Error> {
    let mut best: Option<(usize, Vec<usize>, usize, usize)> = None;
    for first in 0..graph.nodes.len() {
        if !allowed_nodes[first] {
            continue;
        }
        for edge in &graph.nodes[first].edges {
            let second = edge.neighbor;
            if first >= second || !allowed_nodes[second] || is_blocked(first, second, blocked) {
                continue;
            }
            for (inside, outside) in [(first, second), (second, first)] {
                let mut candidate_blocks = blocked.to_vec();
                candidate_blocks.push(edge_key(inside, outside));
                let component = graph.component(inside, &candidate_blocks)?;
                let size = component.leaves.len();
                if size == 0 || size + minimum_complement > total_leaves {
                    continue;
                }
                let scaled =
                    size.checked_mul(balance_factor)
                        .ok_or(Muscle5Error::ArithmeticOverflow {
                            stage: "guide-tree split balance",
                        })?;
                let score = scaled.abs_diff(total_leaves);
                let key = (score, component.leaves, inside, outside);
                if best.as_ref().is_none_or(|current| key < *current) {
                    best = Some(key);
                }
            }
        }
    }
    let (_, _, inside, outside) = best.ok_or(Muscle5Error::InternalInvariant {
        message: "no eligible guide-tree split was found",
    })?;
    Ok(OrientedCut { inside, outside })
}

fn join(left: Muscle5Node, right: Muscle5Node) -> Muscle5Node {
    Muscle5Node::Branch {
        left: Box::new(left),
        right: Box::new(right),
        left_length: 0.0,
        right_length: 0.0,
    }
}

fn edge_key(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn is_blocked(first: usize, second: usize, blocked: &[(usize, usize)]) -> bool {
    blocked.contains(&edge_key(first, second))
}

#[cfg(test)]
mod tests {
    use align_msa::{DistanceMatrix, GuideTreeNode};

    use super::{Muscle5Node, Muscle5Tree, expected_error_upgma, weighted_expected_error_upgma};
    use crate::GuideTreePermutation;

    fn matrix(size: usize, upper_triangle: &[(usize, usize, f64)]) -> DistanceMatrix {
        let mut values = vec![0.0; size * size];
        for &(first, second, value) in upper_triangle {
            values[first * size + second] = value;
            values[second * size + first] = value;
        }
        DistanceMatrix::try_from_flat(size, values).expect("valid hand matrix")
    }

    fn branch(
        left: Muscle5Node,
        right: Muscle5Node,
        left_length: f64,
        right_length: f64,
    ) -> Muscle5Node {
        Muscle5Node::Branch {
            left: Box::new(left),
            right: Box::new(right),
            left_length,
            right_length,
        }
    }

    fn leaf(sequence_index: usize) -> Muscle5Node {
        Muscle5Node::Leaf { sequence_index }
    }

    fn root_groups(node: &Muscle5Node) -> Vec<Vec<usize>> {
        let Muscle5Node::Branch { left, right, .. } = node else {
            panic!("expected a branch");
        };
        let mut groups = vec![left.members(), right.members()];
        groups.sort();
        groups
    }

    #[test]
    fn upgma_uses_weighted_average_and_ultrametric_branch_lengths() {
        let distances = matrix(3, &[(0, 1, 0.2), (0, 2, 0.4), (1, 2, 0.8)]);
        let tree = expected_error_upgma(&distances).expect("UPGMA succeeds");

        let Muscle5Node::Branch {
            left,
            right,
            left_length,
            right_length,
        } = tree.root()
        else {
            panic!("expected root branch");
        };
        assert_eq!(left.members(), vec![0, 1]);
        assert_eq!(right.members(), vec![2]);
        assert!((*left_length - 0.2).abs() < 1.0e-12);
        assert!((*right_length - 0.3).abs() < 1.0e-12);

        let Muscle5Node::Branch {
            left_length,
            right_length,
            ..
        } = left.as_ref()
        else {
            panic!("expected first merge");
        };
        assert!((*left_length - 0.1).abs() < 1.0e-12);
        assert!((*right_length - 0.1).abs() < 1.0e-12);
        assert!(tree.workspace_bytes() > 0);
    }

    #[test]
    fn biased_upgma_uses_supplied_leaf_population_weights() {
        let distances = matrix(3, &[(0, 1, 0.2), (0, 2, 0.4), (1, 2, 0.8)]);
        let tree = weighted_expected_error_upgma(&distances, &[1, 3, 1]).unwrap();
        let Muscle5Node::Branch {
            left,
            right,
            left_length,
            right_length,
        } = tree.root()
        else {
            panic!("expected root branch");
        };
        assert_eq!(left.members(), vec![0, 1]);
        assert_eq!(right.members(), vec![2]);
        assert!((*left_length - 0.25).abs() < 1.0e-12);
        assert!((*right_length - 0.35).abs() < 1.0e-12);
    }

    #[test]
    fn biased_upgma_rejects_missing_or_zero_weights() {
        let distances = matrix(2, &[(0, 1, 0.2)]);
        assert!(weighted_expected_error_upgma(&distances, &[1]).is_err());
        assert!(weighted_expected_error_upgma(&distances, &[1, 0]).is_err());
    }

    #[test]
    fn upgma_ties_are_resolved_by_smallest_active_leaf_pair() {
        let distances = matrix(
            4,
            &[
                (0, 1, 0.5),
                (0, 2, 0.5),
                (0, 3, 0.5),
                (1, 2, 0.5),
                (1, 3, 0.5),
                (2, 3, 0.5),
            ],
        );
        let first = expected_error_upgma(&distances).expect("UPGMA succeeds");
        let second = expected_error_upgma(&distances).expect("UPGMA succeeds");
        assert_eq!(first, second);

        let Muscle5Node::Branch { left, right, .. } = first.root() else {
            panic!("expected root branch");
        };
        assert_eq!(left.members(), vec![0, 1, 2]);
        assert_eq!(right.members(), vec![3]);
    }

    #[test]
    fn whole_tree_permutations_preserve_leaves_and_select_published_resolutions() {
        let pair_01 = branch(leaf(0), leaf(1), 0.1, 0.1);
        let pair_23 = branch(leaf(2), leaf(3), 0.2, 0.2);
        let pair_45 = branch(leaf(4), leaf(5), 0.3, 0.3);
        let remainder = branch(pair_23, pair_45, 0.4, 0.4);
        let tree = Muscle5Tree {
            root: branch(pair_01, remainder, 0.5, 0.5),
            leaf_count: 6,
            workspace_bytes: 123,
        };

        let abc = tree
            .permuted(GuideTreePermutation::Abc)
            .expect("abc permutation");
        let acb = tree
            .permuted(GuideTreePermutation::Acb)
            .expect("acb permutation");
        let bca = tree
            .permuted(GuideTreePermutation::Bca)
            .expect("bca permutation");

        assert_eq!(root_groups(abc.root()), vec![vec![0, 1, 2, 3], vec![4, 5]]);
        assert_eq!(root_groups(acb.root()), vec![vec![0, 1, 4, 5], vec![2, 3]]);
        assert_eq!(root_groups(bca.root()), vec![vec![0, 1], vec![2, 3, 4, 5]]);
        for permuted in [&abc, &acb, &bca] {
            assert_eq!(permuted.root().members(), vec![0, 1, 2, 3, 4, 5]);
            assert_eq!(permuted.workspace_bytes(), 123);
        }
    }

    #[test]
    fn public_tree_conversion_validates_every_leaf() {
        let distances = matrix(2, &[(0, 1, 0.6)]);
        let tree = expected_error_upgma(&distances).expect("UPGMA succeeds");
        let public = tree.to_public().expect("valid public guide tree");
        assert_eq!(public.leaf_count(), 2);
        assert_eq!(public.postorder_leaves(), vec![0, 1]);
        let GuideTreeNode::Branch { left, right } = public.root() else {
            panic!("expected public branch");
        };
        assert!((left.length() - 0.3).abs() < 1.0e-12);
        assert!((right.length() - 0.3).abs() < 1.0e-12);
    }

    #[test]
    fn singleton_tree_is_preserved_by_every_permutation() {
        let distances = matrix(1, &[]);
        let tree = expected_error_upgma(&distances).expect("singleton tree");
        for permutation in [
            GuideTreePermutation::None,
            GuideTreePermutation::Abc,
            GuideTreePermutation::Acb,
            GuideTreePermutation::Bca,
        ] {
            assert_eq!(tree.permuted(permutation).expect("permutation"), tree);
        }
    }
}
