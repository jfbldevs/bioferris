use core::mem::size_of;

use align_msa::{GuideTree, GuideTreeNode};

use crate::{FftNsiParameters, MafftError};

/// One canonical tree-edge bipartition examined by FFT-NS-i.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FftNsiPartition {
    left_members: Box<[usize]>,
    right_members: Box<[usize]>,
}

impl FftNsiPartition {
    /// Returns the canonical first side as original zero-based input indices.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns the complementary side as original zero-based input indices.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }
}

/// Dense symmetric Gotoh pair weights used by one FFT-NS-i objective.
#[derive(Clone, Debug, PartialEq)]
pub struct FftNsiPairWeights {
    size: usize,
    values: Box<[f64]>,
}

impl FftNsiPairWeights {
    /// Returns the number of sequence rows represented by the matrix.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// Returns one pair weight, or `None` when either index is out of range.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.offset(row, column).map(|offset| self.values[offset])
    }

    /// Returns the row-major dense matrix, whose diagonal is zero.
    #[must_use]
    pub const fn as_flat_slice(&self) -> &[f64] {
        &self.values
    }

    /// Returns the sum over distinct unordered sequence pairs.
    #[must_use]
    pub fn sum_of_pairs(&self) -> f64 {
        let mut total = 0.0;
        for row in 0..self.size {
            for column in 0..row {
                total += self.values[row * self.size + column];
            }
        }
        total
    }

    /// Returns the smallest distinct-pair weight, or `None` for a singleton.
    #[must_use]
    pub fn minimum_pair_weight(&self) -> Option<f64> {
        self.distinct_pair_values().min_by(f64::total_cmp)
    }

    /// Returns the largest distinct-pair weight, or `None` for a singleton.
    #[must_use]
    pub fn maximum_pair_weight(&self) -> Option<f64> {
        self.distinct_pair_values().max_by(f64::total_cmp)
    }

    pub(crate) fn resident_bytes(&self) -> Result<usize, MafftError> {
        self.values
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-i pair-weight memory accounting",
            })
    }

    fn distinct_pair_values(&self) -> impl Iterator<Item = f64> + '_ {
        (0..self.size)
            .flat_map(move |row| (0..row).map(move |column| self.values[row * self.size + column]))
    }

    fn offset(&self, row: usize, column: usize) -> Option<usize> {
        (row < self.size && column < self.size).then_some(row * self.size + column)
    }
}

#[derive(Clone, Debug)]
struct TreeEdge {
    first: usize,
    second: usize,
    length: f64,
    log_factor: f64,
}

impl TreeEdge {
    fn other(&self, node: usize) -> Result<usize, MafftError> {
        if node == self.first {
            Ok(self.second)
        } else if node == self.second {
            Ok(self.first)
        } else {
            Err(MafftError::InternalInvariant {
                message: "unrooted edge does not contain the requested node",
            })
        }
    }
}

#[derive(Clone, Debug)]
struct UnrootedTree {
    leaf_count: usize,
    leaf_at_node: Vec<Option<usize>>,
    adjacency: Vec<Vec<usize>>,
    edges: Vec<TreeEdge>,
}

impl UnrootedTree {
    fn from_guide_tree(tree: &GuideTree) -> Result<Self, MafftError> {
        let leaf_count = tree.leaf_count();
        let mut result = Self {
            leaf_count,
            leaf_at_node: (0..leaf_count).map(Some).collect(),
            adjacency: vec![Vec::new(); leaf_count],
            edges: Vec::with_capacity(leaf_count.saturating_mul(2).saturating_sub(3)),
        };

        match tree.root() {
            GuideTreeNode::Leaf { sequence_index } => {
                if leaf_count != 1 || *sequence_index != 0 {
                    return Err(MafftError::InternalInvariant {
                        message: "singleton guide tree has an invalid leaf",
                    });
                }
            }
            GuideTreeNode::Branch { left, right } => {
                let left_node = result.build_subtree(left.child())?;
                let right_node = result.build_subtree(right.child())?;
                let root_length = left.length() + right.length();
                if !root_length.is_finite() {
                    return Err(MafftError::NonFiniteComputation {
                        stage: "FFT-NS-i root suppression",
                    });
                }
                result.add_edge(left_node, right_node, root_length)?;
            }
        }
        result.validate_shape()?;
        Ok(result)
    }

    fn build_subtree(&mut self, node: &GuideTreeNode) -> Result<usize, MafftError> {
        match node {
            GuideTreeNode::Leaf { sequence_index } => Ok(*sequence_index),
            GuideTreeNode::Branch { left, right } => {
                let internal = self.leaf_at_node.len();
                self.leaf_at_node.push(None);
                self.adjacency.push(Vec::with_capacity(3));
                let left_node = self.build_subtree(left.child())?;
                let right_node = self.build_subtree(right.child())?;
                self.add_edge(internal, left_node, left.length())?;
                self.add_edge(internal, right_node, right.length())?;
                Ok(internal)
            }
        }
    }

    fn add_edge(&mut self, first: usize, second: usize, length: f64) -> Result<(), MafftError> {
        if first >= self.adjacency.len() || second >= self.adjacency.len() || first == second {
            return Err(MafftError::InternalInvariant {
                message: "invalid unrooted tree edge endpoints",
            });
        }
        if !length.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT-NS-i unrooted branch construction",
            });
        }
        let edge_index = self.edges.len();
        self.edges.push(TreeEdge {
            first,
            second,
            length,
            log_factor: 0.0,
        });
        self.adjacency[first].push(edge_index);
        self.adjacency[second].push(edge_index);
        Ok(())
    }

    fn validate_shape(&self) -> Result<(), MafftError> {
        let expected_edges = self.leaf_count.saturating_mul(2).saturating_sub(3);
        if self.edges.len() != expected_edges {
            return Err(MafftError::InternalInvariant {
                message: "suppressed guide tree does not have 2N-3 edges",
            });
        }
        for (node, leaf) in self.leaf_at_node.iter().enumerate() {
            let expected_degree = if leaf.is_some() {
                usize::from(self.leaf_count > 1)
            } else {
                3
            };
            if self.adjacency[node].len() != expected_degree {
                return Err(MafftError::InternalInvariant {
                    message: "unrooted binary tree has an invalid node degree",
                });
            }
        }
        Ok(())
    }

    #[allow(clippy::similar_names)]
    fn assign_gotoh_factors(
        &mut self,
        equalization_factor: f64,
        branch_length_floor: f64,
    ) -> Result<(), MafftError> {
        if self.edges.is_empty() {
            return Ok(());
        }
        if self.edges.iter().all(|edge| edge.length <= 0.0) {
            for edge in &mut self.edges {
                edge.log_factor = 0.0;
            }
            return Ok(());
        }

        let mut local_logs = vec![Vec::<(usize, f64)>::new(); self.adjacency.len()];
        for (node, local_at_node) in local_logs.iter_mut().enumerate() {
            if self.leaf_at_node[node].is_some() {
                continue;
            }
            let edge_indices = &self.adjacency[node];
            if edge_indices.len() != 3 {
                return Err(MafftError::InternalInvariant {
                    message: "Gotoh three-way node does not have three directions",
                });
            }
            let mut effective = [0.0; 3];
            for (slot, edge_index) in edge_indices.iter().copied().enumerate() {
                effective[slot] =
                    self.effective_direction(node, edge_index, branch_length_floor)?;
            }
            let [a, b, c] = effective;
            let sum = a * b + a * c + b * c;
            if !sum.is_finite() || sum <= 0.0 {
                return Err(MafftError::NonFiniteComputation {
                    stage: "FFT-NS-i Gotoh three-way normalization",
                });
            }
            let denominator = 2.0 * sum;
            let pair_ab = c * (a + b) / denominator;
            let pair_ac = b * (a + c) / denominator;
            let pair_bc = a * (b + c) / denominator;
            if [pair_ab, pair_ac, pair_bc]
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            {
                return Err(MafftError::NonFiniteComputation {
                    stage: "FFT-NS-i Gotoh local pair weights",
                });
            }
            let exponent = equalization_factor;
            let log_a = 0.5 * (pair_ab.ln() + pair_ac.ln() - pair_bc.ln()) * exponent;
            let log_b = 0.5 * (pair_ab.ln() + pair_bc.ln() - pair_ac.ln()) * exponent;
            let log_c = 0.5 * (pair_ac.ln() + pair_bc.ln() - pair_ab.ln()) * exponent;
            for (edge_index, log_factor) in edge_indices.iter().copied().zip([log_a, log_b, log_c])
            {
                if !log_factor.is_finite() {
                    return Err(MafftError::NonFiniteComputation {
                        stage: "FFT-NS-i Gotoh local edge factors",
                    });
                }
                local_at_node.push((edge_index, log_factor));
            }
        }

        for (edge_index, edge) in self.edges.iter_mut().enumerate() {
            let mut log_factor = 0.0;
            for endpoint in [edge.first, edge.second] {
                if self.leaf_at_node[endpoint].is_none() {
                    let local = local_logs[endpoint]
                        .iter()
                        .find_map(|(candidate, value)| (*candidate == edge_index).then_some(*value))
                        .ok_or(MafftError::InternalInvariant {
                            message: "missing Gotoh factor at an internal edge endpoint",
                        })?;
                    log_factor += local;
                }
            }
            if !log_factor.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "FFT-NS-i Gotoh edge factorization",
                });
            }
            edge.log_factor = log_factor;
        }
        Ok(())
    }

    fn effective_direction(
        &self,
        from: usize,
        edge_index: usize,
        floor: f64,
    ) -> Result<f64, MafftError> {
        let edge = &self.edges[edge_index];
        let neighbor = edge.other(from)?;
        let length = edge.length.max(floor);
        let effective = if self.leaf_at_node[neighbor].is_some() {
            length
        } else {
            let remaining = self.adjacency[neighbor]
                .iter()
                .copied()
                .filter(|candidate| *candidate != edge_index)
                .collect::<Vec<_>>();
            if remaining.len() != 2 {
                return Err(MafftError::InternalInvariant {
                    message: "effective-resistance recursion found a non-binary node",
                });
            }
            let first = self.effective_direction(neighbor, remaining[0], floor)?;
            let second = self.effective_direction(neighbor, remaining[1], floor)?;
            let smaller = first.min(second);
            let larger = first.max(second);
            let parallel = smaller / (1.0 + smaller / larger);
            length + parallel
        };
        if !effective.is_finite() || effective <= 0.0 {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT-NS-i synthetic resistance",
            });
        }
        Ok(effective)
    }

    fn partitions(&self) -> Result<Vec<TreePartition>, MafftError> {
        let mut partitions = Vec::with_capacity(self.edges.len());
        for (edge_index, edge) in self.edges.iter().enumerate() {
            let mut first_members = Vec::new();
            self.collect_component_leaves(edge.first, edge_index, &mut first_members)?;
            first_members.sort_unstable();
            let mut second_members = Vec::new();
            self.collect_component_leaves(edge.second, edge_index, &mut second_members)?;
            second_members.sort_unstable();

            let keep_first = first_members.len() < second_members.len()
                || (first_members.len() == second_members.len() && first_members < second_members);
            let (left_members, right_members, left_endpoint, right_endpoint) = if keep_first {
                (first_members, second_members, edge.first, edge.second)
            } else {
                (second_members, first_members, edge.second, edge.first)
            };
            partitions.push(TreePartition {
                public: FftNsiPartition {
                    left_members: left_members.into_boxed_slice(),
                    right_members: right_members.into_boxed_slice(),
                },
                edge_index,
                left_endpoint,
                right_endpoint,
            });
        }
        partitions.sort_by(|first, second| {
            first
                .public
                .left_members
                .as_ref()
                .cmp(second.public.left_members.as_ref())
                .then_with(|| {
                    first
                        .public
                        .right_members
                        .as_ref()
                        .cmp(second.public.right_members.as_ref())
                })
        });
        for pair in partitions.windows(2) {
            if pair[0].public == pair[1].public {
                return Err(MafftError::InternalInvariant {
                    message: "unrooted guide tree produced a duplicate bipartition",
                });
            }
        }
        Ok(partitions)
    }

    fn collect_component_leaves(
        &self,
        node: usize,
        blocked_edge: usize,
        output: &mut Vec<usize>,
    ) -> Result<(), MafftError> {
        if let Some(leaf) = self.leaf_at_node[node] {
            output.push(leaf);
            return Ok(());
        }
        for edge_index in self.adjacency[node].iter().copied() {
            if edge_index == blocked_edge {
                continue;
            }
            let next = self.edges[edge_index].other(node)?;
            self.collect_component_leaves_from(next, node, output)?;
        }
        Ok(())
    }

    fn collect_component_leaves_from(
        &self,
        node: usize,
        parent: usize,
        output: &mut Vec<usize>,
    ) -> Result<(), MafftError> {
        if let Some(leaf) = self.leaf_at_node[node] {
            output.push(leaf);
            return Ok(());
        }
        for edge_index in self.adjacency[node].iter().copied() {
            let next = self.edges[edge_index].other(node)?;
            if next != parent {
                self.collect_component_leaves_from(next, node, output)?;
            }
        }
        Ok(())
    }

    fn pair_weights(&self) -> Result<FftNsiPairWeights, MafftError> {
        let matrix_len =
            self.leaf_count
                .checked_mul(self.leaf_count)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "FFT-NS-i pair-weight matrix dimensions",
                })?;
        let mut values = vec![0.0; matrix_len];
        if self.leaf_count < 2 {
            return Ok(FftNsiPairWeights {
                size: self.leaf_count,
                values: values.into_boxed_slice(),
            });
        }

        let mut logarithms = Vec::with_capacity(
            self.leaf_count
                .checked_mul(self.leaf_count - 1)
                .and_then(|value| value.checked_div(2))
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "FFT-NS-i pair-weight count",
                })?,
        );
        let mut maximum_log = f64::NEG_INFINITY;
        for first in 0..self.leaf_count {
            for second in 0..first {
                let log_weight =
                    self.path_log(first, second, None)?
                        .ok_or(MafftError::InternalInvariant {
                            message: "unrooted tree has no path between two leaves",
                        })?;
                maximum_log = maximum_log.max(log_weight);
                logarithms.push((first, second, log_weight));
            }
        }
        let mut shifted_sum = 0.0;
        for (_, _, log_weight) in &logarithms {
            shifted_sum += (*log_weight - maximum_log).exp().max(f64::MIN_POSITIVE);
        }
        #[allow(clippy::cast_precision_loss)]
        let target_sum = logarithms.len() as f64;
        let scale = target_sum / shifted_sum;
        if !scale.is_finite() || scale <= 0.0 {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT-NS-i pair-weight normalization",
            });
        }
        for (first, second, log_weight) in logarithms {
            let weight = (log_weight - maximum_log).exp().max(f64::MIN_POSITIVE) * scale;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(MafftError::NonFiniteComputation {
                    stage: "FFT-NS-i normalized pair weights",
                });
            }
            values[first * self.leaf_count + second] = weight;
            values[second * self.leaf_count + first] = weight;
        }
        Ok(FftNsiPairWeights {
            size: self.leaf_count,
            values: values.into_boxed_slice(),
        })
    }

    fn profile_weights(
        &self,
        members: &[usize],
        endpoint: usize,
        blocked_edge: usize,
    ) -> Result<Vec<f64>, MafftError> {
        if members.is_empty() {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-i partition contains an empty subgroup",
            });
        }
        let mut logarithms = Vec::with_capacity(members.len());
        let mut maximum_log = f64::NEG_INFINITY;
        for member in members.iter().copied() {
            let log_weight = self.path_log(member, endpoint, Some(blocked_edge))?.ok_or(
                MafftError::InternalInvariant {
                    message: "partition leaf cannot reach its cut endpoint",
                },
            )?;
            maximum_log = maximum_log.max(log_weight);
            logarithms.push(log_weight);
        }
        let mut weights = logarithms
            .iter()
            .map(|value| (*value - maximum_log).exp().max(f64::MIN_POSITIVE))
            .collect::<Vec<_>>();
        let total: f64 = weights.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT-NS-i subgroup weight normalization",
            });
        }
        for weight in &mut weights {
            *weight /= total;
        }
        Ok(weights)
    }

    fn path_log(
        &self,
        start: usize,
        target: usize,
        blocked_edge: Option<usize>,
    ) -> Result<Option<f64>, MafftError> {
        self.path_log_from(start, usize::MAX, target, blocked_edge)
    }

    fn path_log_from(
        &self,
        node: usize,
        parent: usize,
        target: usize,
        blocked_edge: Option<usize>,
    ) -> Result<Option<f64>, MafftError> {
        if node == target {
            return Ok(Some(0.0));
        }
        for edge_index in self.adjacency[node].iter().copied() {
            if Some(edge_index) == blocked_edge {
                continue;
            }
            let next = self.edges[edge_index].other(node)?;
            if next == parent {
                continue;
            }
            if let Some(suffix) = self.path_log_from(next, node, target, blocked_edge)? {
                let result = self.edges[edge_index].log_factor + suffix;
                if !result.is_finite() {
                    return Err(MafftError::NonFiniteComputation {
                        stage: "FFT-NS-i logarithmic path weight",
                    });
                }
                return Ok(Some(result));
            }
        }
        Ok(None)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct TreePartition {
    pub(crate) public: FftNsiPartition,
    edge_index: usize,
    left_endpoint: usize,
    right_endpoint: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct GotohWeighting {
    tree: UnrootedTree,
    pair_weights: FftNsiPairWeights,
    partitions: Vec<TreePartition>,
}

impl GotohWeighting {
    pub(crate) fn from_guide_tree(
        guide_tree: &GuideTree,
        parameters: FftNsiParameters,
    ) -> Result<Self, MafftError> {
        Self::from_guide_tree_values(
            guide_tree,
            parameters.equalization_factor(),
            parameters.branch_length_floor(),
        )
    }

    pub(crate) fn from_guide_tree_values(
        guide_tree: &GuideTree,
        equalization_factor: f64,
        branch_length_floor: f64,
    ) -> Result<Self, MafftError> {
        if !equalization_factor.is_finite()
            || equalization_factor <= 0.0
            || !branch_length_floor.is_finite()
            || branch_length_floor <= 0.0
        {
            return Err(MafftError::InternalInvariant {
                message: "Gotoh weighting parameters are invalid",
            });
        }
        let mut tree = UnrootedTree::from_guide_tree(guide_tree)?;
        tree.assign_gotoh_factors(equalization_factor, branch_length_floor)?;
        let pair_weights = tree.pair_weights()?;
        let partitions = tree.partitions()?;
        Ok(Self {
            tree,
            pair_weights,
            partitions,
        })
    }

    pub(crate) const fn pair_weights(&self) -> &FftNsiPairWeights {
        &self.pair_weights
    }

    pub(crate) fn partitions(&self) -> &[TreePartition] {
        &self.partitions
    }

    pub(crate) fn partition_profile_weights(
        &self,
        partition: &TreePartition,
    ) -> Result<(Vec<f64>, Vec<f64>), MafftError> {
        let left = self.tree.profile_weights(
            partition.public.left_members(),
            partition.left_endpoint,
            partition.edge_index,
        )?;
        let right = self.tree.profile_weights(
            partition.public.right_members(),
            partition.right_endpoint,
            partition.edge_index,
        )?;
        Ok((left, right))
    }
}

#[cfg(test)]
mod tests {
    use align_msa::{GuideTree, GuideTreeEdge, GuideTreeNode};

    use super::{GotohWeighting, UnrootedTree};
    use crate::FftNsiParameters;

    fn edge(length: f64, child: GuideTreeNode) -> GuideTreeEdge {
        GuideTreeEdge::new(length, child).unwrap()
    }

    fn four_leaf_tree() -> GuideTree {
        let left = GuideTreeNode::branch(
            edge(0.1, GuideTreeNode::leaf(0)),
            edge(0.2, GuideTreeNode::leaf(1)),
        );
        let right = GuideTreeNode::branch(
            edge(0.3, GuideTreeNode::leaf(2)),
            edge(0.4, GuideTreeNode::leaf(3)),
        );
        GuideTree::try_new(GuideTreeNode::branch(edge(0.5, left), edge(0.6, right)), 4).unwrap()
    }

    #[test]
    fn root_suppression_yields_every_unrooted_edge_once() {
        let weighting =
            GotohWeighting::from_guide_tree(&four_leaf_tree(), FftNsiParameters::paper_2002())
                .unwrap();
        assert_eq!(weighting.partitions().len(), 5);
        let splits = weighting
            .partitions()
            .iter()
            .map(|partition| {
                (
                    partition.public.left_members().to_vec(),
                    partition.public.right_members().to_vec(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            splits,
            vec![
                (vec![0], vec![1, 2, 3]),
                (vec![0, 1], vec![2, 3]),
                (vec![1], vec![0, 2, 3]),
                (vec![2], vec![0, 1, 3]),
                (vec![3], vec![0, 1, 2]),
            ]
        );
    }

    #[test]
    fn gotoh_pair_weights_are_positive_symmetric_and_mean_one() {
        let weighting =
            GotohWeighting::from_guide_tree(&four_leaf_tree(), FftNsiParameters::paper_2002())
                .unwrap();
        let weights = weighting.pair_weights();
        assert_eq!(weights.size(), 4);
        assert!((weights.sum_of_pairs() - 6.0).abs() < 1.0e-12);
        for first in 0..4 {
            assert_eq!(weights.get(first, first), Some(0.0));
            for second in 0..first {
                let forward = weights.get(first, second).unwrap();
                assert!(forward > 0.0 && forward.is_finite());
                assert_eq!(
                    forward.to_bits(),
                    weights.get(second, first).unwrap().to_bits()
                );
            }
        }
    }

    #[test]
    fn subgroup_weights_factor_every_crossing_pair() {
        let weighting =
            GotohWeighting::from_guide_tree(&four_leaf_tree(), FftNsiParameters::paper_2002())
                .unwrap();
        let partition = weighting
            .partitions()
            .iter()
            .find(|partition| partition.public.left_members() == [0, 1])
            .unwrap();
        let (left, right) = weighting.partition_profile_weights(partition).unwrap();
        let pair = weighting.pair_weights();
        let mut ratios = Vec::new();
        for (left_slot, first) in partition.public.left_members().iter().copied().enumerate() {
            for (right_slot, second) in partition.public.right_members().iter().copied().enumerate()
            {
                ratios
                    .push(pair.get(first, second).unwrap() / (left[left_slot] * right[right_slot]));
            }
        }
        let reference = ratios[0];
        assert!(
            ratios
                .iter()
                .all(|ratio| (*ratio - reference).abs() <= 1.0e-11 * reference.abs())
        );
    }

    #[test]
    fn entirely_zero_tree_uses_equal_pair_weights() {
        let tree = GuideTree::try_new(
            GuideTreeNode::branch(
                edge(0.0, GuideTreeNode::leaf(0)),
                edge(
                    0.0,
                    GuideTreeNode::branch(
                        edge(0.0, GuideTreeNode::leaf(1)),
                        edge(0.0, GuideTreeNode::leaf(2)),
                    ),
                ),
            ),
            3,
        )
        .unwrap();
        let weighting =
            GotohWeighting::from_guide_tree(&tree, FftNsiParameters::paper_2002()).unwrap();
        assert!((weighting.pair_weights().sum_of_pairs() - 3.0).abs() < 1.0e-12);
        assert_eq!(weighting.pair_weights().get(0, 1), Some(1.0));
        assert_eq!(weighting.pair_weights().get(0, 2), Some(1.0));
        assert_eq!(weighting.pair_weights().get(1, 2), Some(1.0));
    }

    #[test]
    fn singleton_unrooted_tree_has_no_partitions() {
        let guide = GuideTree::try_new(GuideTreeNode::leaf(0), 1).unwrap();
        let tree = UnrootedTree::from_guide_tree(&guide).unwrap();
        assert!(tree.edges.is_empty());
        assert!(tree.partitions().unwrap().is_empty());
        assert!(tree.pair_weights().unwrap().sum_of_pairs().abs() < f64::EPSILON);
    }
}
