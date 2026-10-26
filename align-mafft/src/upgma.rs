use align_msa::{DistanceMatrix, GuideTree, GuideTreeEdge, GuideTreeNode};

use crate::MafftError;

#[derive(Clone, Debug)]
struct Cluster {
    node: GuideTreeNode,
    size: usize,
    minimum_leaf: usize,
    height: f64,
}

#[derive(Clone, Copy, Debug)]
enum Linkage {
    SizeWeightedAverage,
    MafftMixed { arithmetic_fraction: f64 },
}

pub(crate) fn build_upgma(distances: &DistanceMatrix) -> Result<GuideTree, MafftError> {
    build_upgma_from_values(distances.size(), distances.as_flat_slice())
}

pub(crate) fn build_upgma_from_values(
    count: usize,
    distances: &[f64],
) -> Result<GuideTree, MafftError> {
    build_linkage_tree(count, distances, Linkage::SizeWeightedAverage)
}

pub(crate) fn build_mixed_linkage_from_values(
    count: usize,
    distances: &[f64],
    arithmetic_fraction: f64,
) -> Result<GuideTree, MafftError> {
    if !arithmetic_fraction.is_finite() || !(0.0..=1.0).contains(&arithmetic_fraction) {
        return Err(MafftError::InternalInvariant {
            message: "MAFFT mixed-linkage fraction is outside [0, 1]",
        });
    }
    build_linkage_tree(
        count,
        distances,
        Linkage::MafftMixed {
            arithmetic_fraction,
        },
    )
}

fn build_linkage_tree(
    count: usize,
    distances: &[f64],
    linkage: Linkage,
) -> Result<GuideTree, MafftError> {
    if count == 0 {
        return Err(align_msa::MsaError::EmptyCollection.into());
    }
    let expected = count
        .checked_mul(count)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "UPGMA distance matrix dimensions",
        })?;
    if distances.len() != expected {
        return Err(MafftError::InternalInvariant {
            message: "UPGMA distance matrix shape mismatch",
        });
    }
    validate_distance_values(count, distances)?;
    let mut matrix = distances.to_vec();
    let mut clusters = (0..count)
        .map(|index| {
            Some(Cluster {
                node: GuideTreeNode::leaf(index),
                size: 1,
                minimum_leaf: index,
                height: 0.0,
            })
        })
        .collect::<Vec<_>>();
    let mut active = count;

    while active > 1 {
        let (first_index, second_index, merge_distance) = best_pair(&clusters, &matrix, count)
            .ok_or(MafftError::InternalInvariant {
                message: "UPGMA has multiple active clusters but no pair",
            })?;

        let first = clusters[first_index]
            .take()
            .ok_or(MafftError::InternalInvariant {
                message: "UPGMA selected an inactive first cluster",
            })?;
        let second = clusters[second_index]
            .take()
            .ok_or(MafftError::InternalInvariant {
                message: "UPGMA selected an inactive second cluster",
            })?;

        let combined_size =
            first
                .size
                .checked_add(second.size)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "UPGMA cluster size",
                })?;
        for other in 0..count {
            if clusters[other].is_none() {
                continue;
            }
            let first_distance = matrix[first_index * count + other];
            let second_distance = matrix[second_index * count + other];
            #[allow(clippy::cast_precision_loss)]
            let updated = match linkage {
                Linkage::SizeWeightedAverage => {
                    (first.size as f64 * first_distance + second.size as f64 * second_distance)
                        / combined_size as f64
                }
                Linkage::MafftMixed {
                    arithmetic_fraction,
                } => {
                    arithmetic_fraction * (first_distance + second_distance) * 0.5
                        + (1.0 - arithmetic_fraction) * first_distance.min(second_distance)
                }
            };
            matrix[first_index * count + other] = updated;
            matrix[other * count + first_index] = updated;
        }
        matrix[first_index * count + first_index] = 0.0;

        let parent_height = merge_distance * 0.5;
        let (left, right) = if first.minimum_leaf <= second.minimum_leaf {
            (first, second)
        } else {
            (second, first)
        };
        let branch_length = |child_height: f64| match linkage {
            Linkage::SizeWeightedAverage => parent_height - child_height,
            Linkage::MafftMixed { .. } => (parent_height - child_height).max(0.0),
        };
        let left_edge = GuideTreeEdge::new(branch_length(left.height), left.node)?;
        let right_edge = GuideTreeEdge::new(branch_length(right.height), right.node)?;
        clusters[first_index] = Some(Cluster {
            node: GuideTreeNode::branch(left_edge, right_edge),
            size: combined_size,
            minimum_leaf: left.minimum_leaf.min(right.minimum_leaf),
            height: parent_height,
        });
        active -= 1;
    }

    let root = clusters
        .into_iter()
        .flatten()
        .next()
        .ok_or(MafftError::InternalInvariant {
            message: "UPGMA lost its root cluster",
        })?
        .node;
    GuideTree::try_new(root, count).map_err(Into::into)
}

fn validate_distance_values(count: usize, distances: &[f64]) -> Result<(), MafftError> {
    for row in 0..count {
        for column in 0..count {
            let value = distances[row * count + column];
            if !value.is_finite() || value < 0.0 {
                return Err(MafftError::InternalInvariant {
                    message: "UPGMA distance matrix contains an invalid value",
                });
            }
            if row == column && value != 0.0 {
                return Err(MafftError::InternalInvariant {
                    message: "UPGMA distance matrix has a nonzero diagonal",
                });
            }
            if column < row && value.total_cmp(&distances[column * count + row]).is_ne() {
                return Err(MafftError::InternalInvariant {
                    message: "UPGMA distance matrix is asymmetric",
                });
            }
        }
    }
    Ok(())
}

fn best_pair(
    clusters: &[Option<Cluster>],
    matrix: &[f64],
    count: usize,
) -> Option<(usize, usize, f64)> {
    let mut best: Option<(usize, usize, f64, (usize, usize))> = None;
    for first in 0..count {
        let Some(first_cluster) = &clusters[first] else {
            continue;
        };
        for second in (first + 1)..count {
            let Some(second_cluster) = &clusters[second] else {
                continue;
            };
            let distance = matrix[first * count + second];
            let key = if first_cluster.minimum_leaf <= second_cluster.minimum_leaf {
                (first_cluster.minimum_leaf, second_cluster.minimum_leaf)
            } else {
                (second_cluster.minimum_leaf, first_cluster.minimum_leaf)
            };
            let replace = match best {
                None => true,
                Some((_, _, best_distance, best_key)) => {
                    distance.total_cmp(&best_distance).is_lt()
                        || (distance.total_cmp(&best_distance).is_eq() && key < best_key)
                }
            };
            if replace {
                best = Some((first, second, distance, key));
            }
        }
    }
    best.map(|(first, second, distance, _)| (first, second, distance))
}

pub(crate) fn clustalw_sequence_weights(tree: &GuideTree) -> Vec<f64> {
    let mut weights = vec![0.0; tree.leaf_count()];
    add_edge_weights(tree.root(), &mut weights);
    let total: f64 = weights.iter().sum();
    if !total.is_finite() || total <= f64::EPSILON {
        #[allow(clippy::cast_precision_loss)]
        let equal = 1.0 / tree.leaf_count() as f64;
        weights.fill(equal);
        return weights;
    }
    let floor = total * 1.0e-12;
    for weight in &mut weights {
        *weight = weight.max(floor);
    }
    let adjusted_total: f64 = weights.iter().sum();
    for weight in &mut weights {
        *weight /= adjusted_total;
    }
    weights
}

fn add_edge_weights(node: &GuideTreeNode, weights: &mut [f64]) {
    match node {
        GuideTreeNode::Leaf { .. } => {}
        GuideTreeNode::Branch { left, right } => {
            for edge in [left, right] {
                let leaf_count = descendant_count(edge.child());
                #[allow(clippy::cast_precision_loss)]
                let contribution = edge.length().max(0.0) / leaf_count as f64;
                add_to_descendants(edge.child(), contribution, weights);
                add_edge_weights(edge.child(), weights);
            }
        }
    }
}

fn descendant_count(node: &GuideTreeNode) -> usize {
    match node {
        GuideTreeNode::Leaf { .. } => 1,
        GuideTreeNode::Branch { left, right } => {
            descendant_count(left.child()) + descendant_count(right.child())
        }
    }
}

fn add_to_descendants(node: &GuideTreeNode, value: f64, weights: &mut [f64]) {
    match node {
        GuideTreeNode::Leaf { sequence_index } => weights[*sequence_index] += value,
        GuideTreeNode::Branch { left, right } => {
            add_to_descendants(left.child(), value, weights);
            add_to_descendants(right.child(), value, weights);
        }
    }
}

#[cfg(test)]
mod tests {
    use align_msa::{DistanceMatrix, GuideTreeNode};

    use super::{build_mixed_linkage_from_values, build_upgma, clustalw_sequence_weights};

    #[test]
    fn hand_calculated_three_leaf_tree_and_weights() {
        let matrix =
            DistanceMatrix::try_from_flat(3, vec![0.0, 0.2, 0.8, 0.2, 0.0, 0.6, 0.8, 0.6, 0.0])
                .unwrap();
        let tree = build_upgma(&matrix).unwrap();
        let GuideTreeNode::Branch { left, right } = tree.root() else {
            panic!("three leaves require a branch");
        };
        assert_eq!(right.child().minimum_leaf_index(), 2);
        assert!((right.length() - 0.35).abs() < 1.0e-12);
        assert!((left.length() - 0.25).abs() < 1.0e-12);
        let weights = clustalw_sequence_weights(&tree);
        assert!((weights[0] - 0.28125).abs() < 1.0e-12);
        assert!((weights[1] - 0.28125).abs() < 1.0e-12);
        assert!((weights[2] - 0.4375).abs() < 1.0e-12);
    }

    #[test]
    fn zero_length_tree_uses_equal_positive_weights() {
        let matrix = DistanceMatrix::try_from_flat(2, vec![0.0, 0.0, 0.0, 0.0]).unwrap();
        let tree = build_upgma(&matrix).unwrap();
        assert_eq!(clustalw_sequence_weights(&tree), vec![0.5, 0.5]);
    }

    #[test]
    fn ties_use_minimum_input_indices() {
        let matrix = DistanceMatrix::try_from_flat(
            4,
            vec![
                0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0,
            ],
        )
        .unwrap();
        let tree = build_upgma(&matrix).unwrap();
        assert_eq!(tree.postorder_leaves(), vec![0, 1, 2, 3]);
    }

    #[test]
    fn mafft_mixed_linkage_uses_point_one_arithmetic_component() {
        let values = vec![0.0, 0.2, 0.8, 0.2, 0.0, 0.6, 0.8, 0.6, 0.0];
        let tree = build_mixed_linkage_from_values(3, &values, 0.1).unwrap();
        let GuideTreeNode::Branch { left, right } = tree.root() else {
            panic!("three leaves require a branch");
        };
        assert_eq!(right.child().minimum_leaf_index(), 2);
        assert!((right.length() - 0.305).abs() < 1.0e-12);
        assert!((left.length() - 0.205).abs() < 1.0e-12);
    }
}
