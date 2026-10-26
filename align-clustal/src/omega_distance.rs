use core::cmp::Ordering;
use core::mem::size_of;

use align_core::{Alphabet, Sequence};
use align_msa::{DistanceMatrix, GuideTree, GuideTreeEdge, GuideTreeNode, SequenceRecord};

use crate::ClustalOmegaError;

pub(crate) const DISTANCE_SCALE: u64 = 1_000_000_000;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum FixedTreeNode {
    Leaf {
        sequence_index: usize,
    },
    Branch {
        left_length: u64,
        left: Box<Self>,
        right_length: u64,
        right: Box<Self>,
    },
}

impl FixedTreeNode {
    fn into_public(self) -> Result<GuideTreeNode, ClustalOmegaError> {
        match self {
            Self::Leaf { sequence_index } => Ok(GuideTreeNode::leaf(sequence_index)),
            Self::Branch {
                left_length,
                left,
                right_length,
                right,
            } => {
                let left = GuideTreeEdge::new(fixed_to_f64(left_length), left.into_public()?)?;
                let right = GuideTreeEdge::new(fixed_to_f64(right_length), right.into_public()?)?;
                Ok(GuideTreeNode::branch(left, right))
            }
        }
    }

    pub(crate) fn replace_cluster_leaves(
        self,
        subtrees: &mut [Option<FixedTreeNode>],
    ) -> Result<Self, ClustalOmegaError> {
        match self {
            Self::Leaf { sequence_index } => subtrees
                .get_mut(sequence_index)
                .and_then(Option::take)
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "mBed supertree referred to a missing cluster subtree",
                }),
            Self::Branch {
                left_length,
                left,
                right_length,
                right,
            } => {
                let left = left.replace_cluster_leaves(subtrees)?;
                let right = right.replace_cluster_leaves(subtrees)?;
                Ok(Self::Branch {
                    left_length,
                    left: Box::new(left),
                    right_length,
                    right: Box::new(right),
                })
            }
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FixedDistanceMatrix {
    size: usize,
    values: Vec<u64>,
}

impl FixedDistanceMatrix {
    fn zeroed(size: usize) -> Result<Self, ClustalOmegaError> {
        let cells = checked_square(size, "one-tuple distance matrix dimensions")?;
        Ok(Self {
            size,
            values: vec![0; cells],
        })
    }

    pub(crate) const fn size(&self) -> usize {
        self.size
    }

    pub(crate) fn get(&self, row: usize, column: usize) -> u64 {
        self.values[row * self.size + column]
    }

    fn set_symmetric(&mut self, row: usize, column: usize, value: u64) {
        self.values[row * self.size + column] = value;
        self.values[column * self.size + row] = value;
    }

    pub(crate) fn to_public(&self) -> Result<DistanceMatrix, ClustalOmegaError> {
        let values = self
            .values
            .iter()
            .map(|value| fixed_to_f64(*value))
            .collect();
        Ok(DistanceMatrix::try_from_flat(self.size, values)?)
    }

    pub(crate) fn memory_bytes(&self) -> Result<usize, ClustalOmegaError> {
        self.values.len().checked_mul(size_of::<u64>()).ok_or(
            ClustalOmegaError::ArithmeticOverflow {
                stage: "one-tuple distance matrix memory accounting",
            },
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct UpgmaResult {
    pub(crate) root: FixedTreeNode,
    pub(crate) workspace_bytes: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct MbedResult {
    pub(crate) root: FixedTreeNode,
    pub(crate) seed_indices: Vec<usize>,
    pub(crate) clusters: Vec<Vec<usize>>,
    pub(crate) workspace_bytes: usize,
}

#[derive(Clone, Debug)]
struct ActiveCluster {
    node: FixedTreeNode,
    height: u64,
    minimum_leaf: usize,
}

#[derive(Clone, Debug)]
struct VectorCluster {
    members: Vec<usize>,
    barycentre: Vec<u64>,
}

pub(crate) fn one_tuple_distance<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    gap_penalty: u32,
) -> Result<u64, ClustalOmegaError> {
    let denominator = first.len().min(second.len());
    if denominator == 0 {
        return Ok(if first.is_empty() == second.is_empty() {
            0
        } else {
            DISTANCE_SCALE
        });
    }

    let score = one_tuple_similarity_score(first, second, gap_penalty)?;
    let similarity = rounded_ratio(score, denominator, DISTANCE_SCALE)?;
    DISTANCE_SCALE
        .checked_sub(similarity)
        .ok_or(ClustalOmegaError::InternalInvariant {
            message: "one-tuple similarity exceeded its shorter sequence",
        })
}

pub(crate) fn full_distance_matrix<A: Alphabet>(
    records: &[SequenceRecord<A>],
    gap_penalty: u32,
    thread_count: usize,
) -> Result<(FixedDistanceMatrix, usize), ClustalOmegaError> {
    let mut matrix = FixedDistanceMatrix::zeroed(records.len())?;
    let pair_count = records
        .len()
        .checked_mul(records.len().saturating_sub(1))
        .map(|ordered| ordered / 2)
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "one-tuple pair count",
        })?;
    let worker_count = thread_count.min(records.len()).min(pair_count);
    let parallel = worker_count > 1
        && pair_count >= 16
        && records.iter().any(|record| record.sequence().len() >= 64);
    if !parallel {
        for first in 0..records.len() {
            for second in first + 1..records.len() {
                let distance = one_tuple_distance(
                    records[first].sequence(),
                    records[second].sequence(),
                    gap_penalty,
                )?;
                matrix.set_symmetric(first, second, distance);
            }
        }
        return Ok((matrix, 0));
    }
    let pairs = (0..records.len())
        .flat_map(|first| (first + 1..records.len()).map(move |second| (first, second)))
        .collect::<Vec<_>>();
    let chunk_size = pairs.len().div_ceil(worker_count);
    let chunk_count = pairs.len().div_ceil(chunk_size);
    let observations = std::thread::scope(|scope| {
        let calling_start = (chunk_count - 1) * chunk_size;
        let handles = (0..chunk_count - 1)
            .map(|worker| {
                let start = worker * chunk_size;
                let end = (start + chunk_size).min(pairs.len());
                let chunk = &pairs[start..end];
                scope.spawn(move || evaluate_distance_pairs(records, chunk, gap_penalty))
            })
            .collect::<Vec<_>>();
        let mut calling_observations =
            evaluate_distance_pairs(records, &pairs[calling_start..], gap_penalty)?;
        let mut observations = Vec::with_capacity(pairs.len());
        for handle in handles {
            let mut chunk = handle
                .join()
                .map_err(|_| ClustalOmegaError::InternalInvariant {
                    message: "a deterministic guide-distance worker panicked",
                })??;
            observations.append(&mut chunk);
        }
        observations.append(&mut calling_observations);
        Ok::<_, ClustalOmegaError>(observations)
    })?;
    for (first, second, distance) in observations {
        matrix.set_symmetric(first, second, distance);
    }
    let workspace_bytes = pairs
        .len()
        .checked_mul(
            size_of::<(usize, usize)>()
                .checked_add(size_of::<(usize, usize, u64)>())
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "parallel guide-distance element accounting",
                })?,
        )
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "parallel guide-distance workspace accounting",
        })?;
    Ok((matrix, workspace_bytes))
}

fn evaluate_distance_pairs<A: Alphabet>(
    records: &[SequenceRecord<A>],
    pairs: &[(usize, usize)],
    gap_penalty: u32,
) -> Result<Vec<(usize, usize, u64)>, ClustalOmegaError> {
    pairs
        .iter()
        .map(|(first, second)| {
            Ok((
                *first,
                *second,
                one_tuple_distance(
                    records[*first].sequence(),
                    records[*second].sequence(),
                    gap_penalty,
                )?,
            ))
        })
        .collect()
}

fn subset_distance_matrix<A: Alphabet>(
    records: &[SequenceRecord<A>],
    members: &[usize],
    gap_penalty: u32,
) -> Result<FixedDistanceMatrix, ClustalOmegaError> {
    let mut matrix = FixedDistanceMatrix::zeroed(members.len())?;
    for first in 0..members.len() {
        for second in first + 1..members.len() {
            let distance = one_tuple_distance(
                records[members[first]].sequence(),
                records[members[second]].sequence(),
                gap_penalty,
            )?;
            matrix.set_symmetric(first, second, distance);
        }
    }
    Ok(matrix)
}

#[allow(clippy::too_many_lines)]
pub(crate) fn upgma(
    matrix: &FixedDistanceMatrix,
    leaf_indices: &[usize],
) -> Result<UpgmaResult, ClustalOmegaError> {
    if matrix.size() == 0 || matrix.size() != leaf_indices.len() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "UPGMA matrix and leaf collection differ in size",
        });
    }
    let capacity = matrix
        .size()
        .checked_mul(2)
        .and_then(|value| value.checked_sub(1))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "UPGMA cluster capacity",
        })?;
    let cells = checked_square(capacity, "UPGMA workspace dimensions")?;
    let mut distances = vec![0_u64; cells];
    for row in 0..matrix.size() {
        for column in 0..matrix.size() {
            distances[row * capacity + column] = matrix.get(row, column);
        }
    }

    let mut active = vec![false; capacity];
    let mut clusters = vec![None; capacity];
    for (slot, sequence_index) in leaf_indices.iter().copied().enumerate() {
        active[slot] = true;
        clusters[slot] = Some(ActiveCluster {
            node: FixedTreeNode::Leaf { sequence_index },
            height: 0,
            minimum_leaf: sequence_index,
        });
    }
    let mut active_count = matrix.size();
    let mut next_slot = matrix.size();
    while active_count > 1 {
        let (first, second, selected_distance) =
            select_closest(&active, &clusters, &distances, capacity)?;
        let first_cluster = clusters[first]
            .take()
            .ok_or(ClustalOmegaError::InternalInvariant {
                message: "UPGMA selected a missing first cluster",
            })?;
        let second_cluster =
            clusters[second]
                .take()
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "UPGMA selected a missing second cluster",
                })?;
        let height = selected_distance
            .checked_add(1)
            .map(|value| value / 2)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "UPGMA node height",
            })?
            .max(first_cluster.height)
            .max(second_cluster.height);
        let first_length = height - first_cluster.height;
        let second_length = height - second_cluster.height;
        let minimum_leaf = first_cluster.minimum_leaf.min(second_cluster.minimum_leaf);

        for other in 0..next_slot {
            if !active[other] || other == first || other == second {
                continue;
            }
            // Clustal Omega's observable guide-tree update gives both child
            // nodes equal weight, irrespective of their leaf counts.
            let numerator = u128::from(distances[first * capacity + other])
                + u128::from(distances[second * capacity + other]);
            let updated =
                rounded_div_u128(numerator, 2, "Clustal Omega guide-tree node-mean update")?;
            distances[next_slot * capacity + other] = updated;
            distances[other * capacity + next_slot] = updated;
        }

        active[first] = false;
        active[second] = false;
        active[next_slot] = true;
        clusters[next_slot] = Some(ActiveCluster {
            node: FixedTreeNode::Branch {
                left_length: first_length,
                left: Box::new(first_cluster.node),
                right_length: second_length,
                right: Box::new(second_cluster.node),
            },
            height,
            minimum_leaf,
        });
        active_count -= 1;
        next_slot += 1;
    }
    let root_index = active.iter().position(|is_active| *is_active).ok_or(
        ClustalOmegaError::InternalInvariant {
            message: "UPGMA finished without a root cluster",
        },
    )?;
    let root = clusters[root_index]
        .take()
        .ok_or(ClustalOmegaError::InternalInvariant {
            message: "UPGMA root cluster is absent",
        })?
        .node;
    let workspace_bytes = distances
        .len()
        .checked_mul(size_of::<u64>())
        .and_then(|bytes| bytes.checked_add(active.len() * size_of::<bool>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "UPGMA workspace accounting",
        })?;
    Ok(UpgmaResult {
        root,
        workspace_bytes,
    })
}

pub(crate) fn to_public_tree(
    root: FixedTreeNode,
    leaf_count: usize,
) -> Result<GuideTree, ClustalOmegaError> {
    Ok(GuideTree::try_new(root.into_public()?, leaf_count)?)
}

pub(crate) fn mbed<A: Alphabet>(
    records: &[SequenceRecord<A>],
    maximum_cluster_size: usize,
    seed: u64,
    gap_penalty: u32,
) -> Result<MbedResult, ClustalOmegaError> {
    if records.is_empty() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "mBed requires at least one input record",
        });
    }
    let seed_indices = systematic_seed_indices(records);
    let vectors = embedding_vectors(records, &seed_indices, gap_penalty)?;
    let clusters = bisecting_clusters(&vectors, maximum_cluster_size, seed)?;

    let mut subtrees = Vec::with_capacity(clusters.len());
    let mut peak_workspace = vectors
        .len()
        .checked_mul(seed_indices.len())
        .and_then(|cells| cells.checked_mul(size_of::<u64>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "mBed embedding memory accounting",
        })?;
    let mut barycentres = Vec::with_capacity(clusters.len());
    for members in &clusters {
        let distances = subset_distance_matrix(records, members, gap_penalty)?;
        let tree = upgma(&distances, members)?;
        peak_workspace = peak_workspace.max(
            distances
                .memory_bytes()?
                .checked_add(tree.workspace_bytes)
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "mBed leaf-cluster workspace accounting",
                })?,
        );
        subtrees.push(Some(tree.root));
        barycentres.push(integer_barycentre(members, &vectors)?);
    }

    let cluster_matrix = barycentre_distance_matrix(&barycentres)?;
    let cluster_leaves = (0..clusters.len()).collect::<Vec<_>>();
    let supertree = upgma(&cluster_matrix, &cluster_leaves)?;
    peak_workspace = peak_workspace.max(
        cluster_matrix
            .memory_bytes()?
            .checked_add(supertree.workspace_bytes)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "mBed supertree workspace accounting",
            })?,
    );
    let root = supertree.root.replace_cluster_leaves(&mut subtrees)?;
    if subtrees.iter().any(Option::is_some) {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "mBed supertree did not consume every cluster subtree",
        });
    }
    Ok(MbedResult {
        root,
        seed_indices,
        clusters,
        workspace_bytes: peak_workspace,
    })
}

fn select_closest(
    active: &[bool],
    clusters: &[Option<ActiveCluster>],
    distances: &[u64],
    stride: usize,
) -> Result<(usize, usize, u64), ClustalOmegaError> {
    let mut selected: Option<(u64, usize, usize, usize, usize)> = None;
    for first in 0..active.len() {
        if !active[first] {
            continue;
        }
        for second in first + 1..active.len() {
            if !active[second] {
                continue;
            }
            let first_minimum = clusters[first]
                .as_ref()
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "active UPGMA slot has no cluster",
                })?
                .minimum_leaf;
            let second_minimum = clusters[second]
                .as_ref()
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "active UPGMA slot has no cluster",
                })?
                .minimum_leaf;
            let (minimum, maximum) = if first_minimum <= second_minimum {
                (first_minimum, second_minimum)
            } else {
                (second_minimum, first_minimum)
            };
            let candidate = (
                distances[first * stride + second],
                minimum,
                maximum,
                first,
                second,
            );
            if selected.as_ref().is_none_or(|current| candidate < *current) {
                selected = Some(candidate);
            }
        }
    }
    selected
        .map(|(distance, _, _, first, second)| (first, second, distance))
        .ok_or(ClustalOmegaError::InternalInvariant {
            message: "UPGMA could not select two active clusters",
        })
}

fn systematic_seed_indices<A: Alphabet>(records: &[SequenceRecord<A>]) -> Vec<usize> {
    let mut ordered = (0..records.len()).collect::<Vec<_>>();
    ordered.sort_by_key(|index| (records[*index].sequence().len(), *index));
    let requested = ceil_log2(records.len()).max(1).min(records.len());
    let mut seeds: Vec<usize> = Vec::with_capacity(requested);
    for position in 0..requested {
        let ordered_index = if requested == 1 {
            0
        } else {
            position * (ordered.len() - 1) / (requested - 1)
        };
        let candidate = ordered[ordered_index];
        if !seeds.iter().any(|existing| {
            records[*existing].sequence().as_slice() == records[candidate].sequence().as_slice()
        }) {
            seeds.push(candidate);
        }
    }
    if seeds.is_empty() {
        seeds.push(ordered[0]);
    }
    seeds
}

fn embedding_vectors<A: Alphabet>(
    records: &[SequenceRecord<A>],
    seeds: &[usize],
    gap_penalty: u32,
) -> Result<Vec<Vec<u64>>, ClustalOmegaError> {
    records
        .iter()
        .map(|record| {
            seeds
                .iter()
                .map(|seed| {
                    one_tuple_distance(record.sequence(), records[*seed].sequence(), gap_penalty)
                })
                .collect()
        })
        .collect()
}

fn bisecting_clusters(
    vectors: &[Vec<u64>],
    maximum_size: usize,
    seed: u64,
) -> Result<Vec<Vec<usize>>, ClustalOmegaError> {
    let members = (0..vectors.len()).collect::<Vec<_>>();
    let barycentre = integer_barycentre(&members, vectors)?;
    let mut clusters = vec![VectorCluster {
        members,
        barycentre,
    }];
    while clusters
        .iter()
        .any(|cluster| cluster.members.len() > maximum_size)
    {
        let mut selected = None;
        for (index, cluster) in clusters.iter().enumerate() {
            if cluster.members.len() <= maximum_size {
                continue;
            }
            let error = cluster_error(cluster, vectors)?;
            let minimum_member =
                *cluster
                    .members
                    .iter()
                    .min()
                    .ok_or(ClustalOmegaError::InternalInvariant {
                        message: "mBed contains an empty candidate cluster",
                    })?;
            let candidate = (error, core::cmp::Reverse(minimum_member), index);
            if selected.as_ref().is_none_or(|current| candidate > *current) {
                selected = Some(candidate);
            }
        }
        let selected =
            selected
                .map(|(_, _, index)| index)
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "mBed could not select an oversized cluster",
                })?;
        let cluster = clusters.swap_remove(selected);
        let split_seed = derive_cluster_seed(seed, &cluster.members);
        let (left, right) = split_cluster(cluster.members, vectors, split_seed)?;
        clusters.push(VectorCluster {
            barycentre: integer_barycentre(&left, vectors)?,
            members: left,
        });
        clusters.push(VectorCluster {
            barycentre: integer_barycentre(&right, vectors)?,
            members: right,
        });
    }
    let mut result = clusters
        .into_iter()
        .map(|mut cluster| {
            cluster.members.sort_unstable();
            cluster.members
        })
        .collect::<Vec<_>>();
    result.sort_by_key(|members| members[0]);
    Ok(result)
}

fn split_cluster(
    mut members: Vec<usize>,
    vectors: &[Vec<u64>],
    seed: u64,
) -> Result<(Vec<usize>, Vec<usize>), ClustalOmegaError> {
    members.sort_unstable();
    if members.len() < 2 {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "bisecting k-means received fewer than two members",
        });
    }
    let all_equal = members
        .iter()
        .skip(1)
        .all(|member| vectors[*member] == vectors[members[0]]);
    if all_equal {
        let midpoint = members.len() / 2;
        return Ok((members[..midpoint].to_vec(), members[midpoint..].to_vec()));
    }

    let mut random = SplitMix64::new(seed);
    let first_position = usize::try_from(random.next() % members.len() as u64)
        .expect("random member position fits usize");
    let mut centroids = [vectors[members[first_position]].clone(), Vec::new()];
    let weights = members
        .iter()
        .map(|member| squared_distance(&vectors[*member], &centroids[0]))
        .collect::<Result<Vec<_>, _>>()?;
    let total = weights
        .iter()
        .try_fold(0_u128, |sum, weight| sum.checked_add(*weight))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "mBed k-means++ sampling weights",
        })?;
    if total == 0 {
        let midpoint = members.len() / 2;
        return Ok((members[..midpoint].to_vec(), members[midpoint..].to_vec()));
    }
    let draw = random.next_u128() % total;
    let mut cumulative = 0_u128;
    let mut second_position = members.len() - 1;
    for (position, weight) in weights.iter().copied().enumerate() {
        cumulative =
            cumulative
                .checked_add(weight)
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "mBed k-means++ cumulative weights",
                })?;
        if draw < cumulative {
            second_position = position;
            break;
        }
    }
    centroids[1].clone_from(&vectors[members[second_position]]);

    let mut assignments = vec![usize::MAX; members.len()];
    for _ in 0..100 {
        let mut changed = false;
        for (position, member) in members.iter().copied().enumerate() {
            let first = squared_distance(&vectors[member], &centroids[0])?;
            let second = squared_distance(&vectors[member], &centroids[1])?;
            let assignment = usize::from(second < first);
            changed |= assignments[position] != assignment;
            assignments[position] = assignment;
        }
        repair_empty_side(&members, vectors, &centroids, &mut assignments)?;
        let left = members
            .iter()
            .zip(&assignments)
            .filter_map(|(member, side)| (*side == 0).then_some(*member))
            .collect::<Vec<_>>();
        let right = members
            .iter()
            .zip(&assignments)
            .filter_map(|(member, side)| (*side == 1).then_some(*member))
            .collect::<Vec<_>>();
        let next = [
            integer_barycentre(&left, vectors)?,
            integer_barycentre(&right, vectors)?,
        ];
        if !changed || next == centroids {
            return Ok((left, right));
        }
        centroids = next;
    }
    let left = members
        .iter()
        .zip(&assignments)
        .filter_map(|(member, side)| (*side == 0).then_some(*member))
        .collect();
    let right = members
        .iter()
        .zip(&assignments)
        .filter_map(|(member, side)| (*side == 1).then_some(*member))
        .collect();
    Ok((left, right))
}

fn repair_empty_side(
    members: &[usize],
    vectors: &[Vec<u64>],
    centroids: &[Vec<u64>; 2],
    assignments: &mut [usize],
) -> Result<(), ClustalOmegaError> {
    let first_count = assignments.iter().filter(|side| **side == 0).count();
    if first_count != 0 && first_count != assignments.len() {
        return Ok(());
    }
    let occupied = usize::from(first_count == 0);
    let empty = 1 - occupied;
    let mut selected: Option<(u128, core::cmp::Reverse<usize>, usize)> = None;
    for (position, member) in members.iter().copied().enumerate() {
        let distance = squared_distance(&vectors[member], &centroids[occupied])?;
        let candidate = (distance, core::cmp::Reverse(member), position);
        if selected.as_ref().is_none_or(|current| candidate > *current) {
            selected = Some(candidate);
        }
    }
    let position =
        selected
            .map(|(_, _, position)| position)
            .ok_or(ClustalOmegaError::InternalInvariant {
                message: "mBed could not repair an empty k-means side",
            })?;
    assignments[position] = empty;
    Ok(())
}

#[allow(clippy::needless_range_loop)] // Each iteration reduces one column across selected rows.
fn integer_barycentre(
    members: &[usize],
    vectors: &[Vec<u64>],
) -> Result<Vec<u64>, ClustalOmegaError> {
    let dimensions = vectors.first().map_or(0, Vec::len);
    if members.is_empty() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "cannot calculate an empty mBed barycentre",
        });
    }
    let mut result = Vec::with_capacity(dimensions);
    for dimension in 0..dimensions {
        let total = members
            .iter()
            .try_fold(0_u128, |sum, member| {
                sum.checked_add(u128::from(vectors[*member][dimension]))
            })
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "mBed integer barycentre",
            })?;
        result.push(rounded_div_u128(
            total,
            members.len() as u128,
            "mBed integer barycentre",
        )?);
    }
    Ok(result)
}

fn cluster_error(cluster: &VectorCluster, vectors: &[Vec<u64>]) -> Result<u128, ClustalOmegaError> {
    cluster.members.iter().try_fold(0_u128, |sum, member| {
        sum.checked_add(squared_distance(&vectors[*member], &cluster.barycentre)?)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "mBed cluster error",
            })
    })
}

fn squared_distance(first: &[u64], second: &[u64]) -> Result<u128, ClustalOmegaError> {
    if first.len() != second.len() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "mBed vectors have different dimensions",
        });
    }
    first
        .iter()
        .zip(second)
        .try_fold(0_u128, |sum, (first, second)| {
            let difference = first.abs_diff(*second);
            let square = u128::from(difference)
                .checked_mul(u128::from(difference))
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "mBed squared Euclidean distance",
                })?;
            sum.checked_add(square)
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "mBed squared Euclidean distance",
                })
        })
}

fn barycentre_distance_matrix(
    barycentres: &[Vec<u64>],
) -> Result<FixedDistanceMatrix, ClustalOmegaError> {
    let mut matrix = FixedDistanceMatrix::zeroed(barycentres.len())?;
    for first in 0..barycentres.len() {
        for second in first + 1..barycentres.len() {
            let squared = squared_distance(&barycentres[first], &barycentres[second])?;
            let dimensions = barycentres[first].len().max(1) as u128;
            let mean_square =
                rounded_div_u128(squared, dimensions, "mBed normalized Euclidean distance")?;
            let distance = integer_sqrt_nearest(u128::from(mean_square));
            matrix.set_symmetric(first, second, distance.min(DISTANCE_SCALE));
        }
    }
    Ok(matrix)
}

fn integer_sqrt_nearest(value: u128) -> u64 {
    if value == 0 {
        return 0;
    }
    let mut low = 1_u128;
    let mut high = value.min(u128::from(u64::MAX));
    while low <= high {
        let middle = low + (high - low) / 2;
        match middle
            .checked_mul(middle)
            .map_or(Ordering::Greater, |square| square.cmp(&value))
        {
            Ordering::Equal => return u64::try_from(middle).expect("root bounded by u64"),
            Ordering::Less => low = middle + 1,
            Ordering::Greater => high = middle - 1,
        }
    }
    let floor = high;
    let ceiling = low;
    let floor_error = value - floor * floor;
    let ceiling_error = ceiling
        .checked_mul(ceiling)
        .map_or(u128::MAX, |square| square - value);
    u64::try_from(if ceiling_error <= floor_error {
        ceiling
    } else {
        floor
    })
    .expect("root bounded by u64")
}

fn derive_cluster_seed(seed: u64, members: &[usize]) -> u64 {
    let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
    for member in members {
        state ^= (*member as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
        state = splitmix_value(state);
    }
    state
}

#[derive(Clone, Copy, Debug)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        splitmix_value(self.state)
    }

    fn next_u128(&mut self) -> u128 {
        (u128::from(self.next()) << 64) | u128::from(self.next())
    }
}

const fn splitmix_value(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

fn checked_square(value: usize, stage: &'static str) -> Result<usize, ClustalOmegaError> {
    value
        .checked_mul(value)
        .ok_or(ClustalOmegaError::ArithmeticOverflow { stage })
}

fn rounded_ratio(
    numerator: usize,
    denominator: usize,
    scale: u64,
) -> Result<u64, ClustalOmegaError> {
    let numerator = (numerator as u128).checked_mul(u128::from(scale)).ok_or(
        ClustalOmegaError::ArithmeticOverflow {
            stage: "fixed one-tuple ratio",
        },
    )?;
    rounded_div_u128(numerator, denominator as u128, "fixed one-tuple ratio")
}

fn one_tuple_similarity_score<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    gap_penalty: u32,
) -> Result<usize, ClustalOmegaError> {
    let (row_sequence, column_sequence) = if first.len() >= second.len() {
        (first.as_slice(), second.as_slice())
    } else {
        (second.as_slice(), first.as_slice())
    };
    if column_sequence.is_empty() {
        return Ok(0);
    }

    let diagonal_count = row_sequence
        .len()
        .checked_add(column_sequence.len())
        .and_then(|value| value.checked_sub(1))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "sparse one-tuple diagonal dimensions",
        })?;
    let mut positions = vec![Vec::new(); A::SYMBOLS.len()];
    for (column, symbol) in column_sequence.iter().copied().enumerate() {
        positions[symbol.alphabet_index()].push(column);
    }

    let mut diagonal_best = vec![0_usize; diagonal_count];
    let mut dominance_best = vec![0_usize; column_sequence.len() + 1];
    let mut pending = Vec::new();
    let penalty =
        usize::try_from(gap_penalty).map_err(|_| ClustalOmegaError::ArithmeticOverflow {
            stage: "sparse one-tuple gap penalty conversion",
        })?;
    let mut best = 0_usize;

    for (row, symbol) in row_sequence.iter().copied().enumerate() {
        pending.clear();
        for column in positions[symbol.alphabet_index()].iter().copied() {
            let diagonal = row + column_sequence.len() - 1 - column;
            let changed_diagonal =
                fenwick_prefix_max(&dominance_best, column).saturating_sub(penalty);
            let score = diagonal_best[diagonal]
                .max(changed_diagonal)
                .checked_add(1)
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "sparse one-tuple chain score",
                })?;
            best = best.max(score);
            pending.push((column, diagonal, score));
        }
        for (column, diagonal, score) in pending.iter().copied() {
            diagonal_best[diagonal] = diagonal_best[diagonal].max(score);
            fenwick_update_max(&mut dominance_best, column + 1, score)?;
        }
    }
    Ok(best)
}

fn fenwick_prefix_max(tree: &[usize], mut end: usize) -> usize {
    let mut best = 0;
    while end > 0 {
        best = best.max(tree[end]);
        end &= end - 1;
    }
    best
}

fn fenwick_update_max(
    tree: &mut [usize],
    mut position: usize,
    value: usize,
) -> Result<(), ClustalOmegaError> {
    while position < tree.len() {
        tree[position] = tree[position].max(value);
        let step = position & (!position + 1);
        position = position
            .checked_add(step)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "sparse one-tuple Fenwick update",
            })?;
    }
    Ok(())
}

#[cfg(test)]
fn dense_one_tuple_similarity_score<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    gap_penalty: u32,
) -> Result<usize, ClustalOmegaError> {
    const UNREACHABLE: i64 = i64::MIN / 4;
    let (row_sequence, column_sequence) = if first.len() >= second.len() {
        (first.as_slice(), second.as_slice())
    } else {
        (second.as_slice(), first.as_slice())
    };
    let columns =
        column_sequence
            .len()
            .checked_add(1)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "one-tuple score row dimensions",
            })?;
    let penalty = i64::from(gap_penalty);
    let mut previous_match = vec![0_i64; columns];
    let mut previous_gap_in_right = vec![UNREACHABLE; columns];
    let mut previous_gap_in_left = vec![UNREACHABLE; columns];
    let mut current_match = vec![0_i64; columns];
    let mut current_gap_in_right = vec![UNREACHABLE; columns];
    let mut current_gap_in_left = vec![UNREACHABLE; columns];
    let mut best_endpoint = 0_i64;

    for (row_index, row_symbol) in row_sequence.iter().copied().enumerate() {
        current_match[0] = 0;
        current_gap_in_right[0] = UNREACHABLE;
        current_gap_in_left[0] = UNREACHABLE;
        for (column_index, column_symbol) in column_sequence.iter().copied().enumerate() {
            let column = column_index + 1;
            let identity = i64::from(row_symbol == column_symbol);
            current_match[column] = previous_match[column - 1]
                .max(previous_gap_in_right[column - 1])
                .max(previous_gap_in_left[column - 1])
                .checked_add(identity)
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "one-tuple identity score",
                })?;
            let opened_in_right = previous_match[column]
                .max(previous_gap_in_left[column])
                .checked_sub(penalty)
                .unwrap_or(UNREACHABLE);
            current_gap_in_right[column] = opened_in_right.max(previous_gap_in_right[column]);
            let opened_in_left = current_match[column - 1]
                .max(current_gap_in_right[column - 1])
                .checked_sub(penalty)
                .unwrap_or(UNREACHABLE);
            current_gap_in_left[column] = opened_in_left.max(current_gap_in_left[column - 1]);
            if column == columns - 1 || row_index == row_sequence.len() - 1 {
                best_endpoint = best_endpoint
                    .max(current_match[column])
                    .max(current_gap_in_right[column])
                    .max(current_gap_in_left[column]);
            }
        }
        core::mem::swap(&mut previous_match, &mut current_match);
        core::mem::swap(&mut previous_gap_in_right, &mut current_gap_in_right);
        core::mem::swap(&mut previous_gap_in_left, &mut current_gap_in_left);
    }
    usize::try_from(best_endpoint.max(0)).map_err(|_| ClustalOmegaError::ArithmeticOverflow {
        stage: "one-tuple similarity conversion",
    })
}

fn rounded_div_u128(
    numerator: u128,
    denominator: u128,
    stage: &'static str,
) -> Result<u64, ClustalOmegaError> {
    if denominator == 0 {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "rounded division received a zero denominator",
        });
    }
    let rounded = numerator
        .checked_add(denominator / 2)
        .ok_or(ClustalOmegaError::ArithmeticOverflow { stage })?
        / denominator;
    u64::try_from(rounded).map_err(|_| ClustalOmegaError::ArithmeticOverflow { stage })
}

#[allow(clippy::cast_precision_loss)]
fn fixed_to_f64(value: u64) -> f64 {
    value as f64 / DISTANCE_SCALE as f64
}

const fn ceil_log2(value: usize) -> usize {
    if value <= 1 {
        0
    } else {
        usize::BITS as usize - (value - 1).leading_zeros() as usize
    }
}

#[cfg(test)]
mod tests {
    use align_core::{DnaSequence, ProteinSequence};
    use align_msa::{GuideTreeNode, SequenceRecord};

    use super::{
        DISTANCE_SCALE, FixedDistanceMatrix, dense_one_tuple_similarity_score,
        integer_sqrt_nearest, mbed, one_tuple_distance, one_tuple_similarity_score,
        systematic_seed_indices, to_public_tree, upgma,
    };

    fn binary_sequences(maximum_length: usize) -> Vec<DnaSequence> {
        let mut sequences = Vec::new();
        for length in 0..=maximum_length {
            for mask in 0..(1_usize << length) {
                let symbols = (0..length)
                    .map(|position| {
                        if mask & (1 << position) == 0 {
                            b'A'
                        } else {
                            b'C'
                        }
                    })
                    .collect::<Vec<_>>();
                sequences.push(DnaSequence::from_ascii(symbols).unwrap());
            }
        }
        sequences
    }

    #[test]
    fn sparse_one_tuple_score_is_identical_to_the_dense_scalar_oracle() {
        let sequences = binary_sequences(4);
        for gap_penalty in [0, 1, 3, 5] {
            for first in &sequences {
                for second in &sequences {
                    assert_eq!(
                        one_tuple_similarity_score(first, second, gap_penalty).unwrap(),
                        dense_one_tuple_similarity_score(first, second, gap_penalty).unwrap(),
                        "gap={gap_penalty}, first={first}, second={second}"
                    );
                }
            }
        }
    }

    #[test]
    fn one_tuple_distance_matches_the_published_gap_event_score() {
        let first = ProteinSequence::from_ascii(b"AABC").unwrap();
        let second = ProteinSequence::from_ascii(b"ABBD").unwrap();
        assert_eq!(
            one_tuple_distance(&first, &second, 3).unwrap(),
            DISTANCE_SCALE / 2
        );

        let reference = DnaSequence::from_ascii(b"ACGTACGT").unwrap();
        let deletion = DnaSequence::from_ascii(b"ACGACGT").unwrap();
        assert_eq!(
            one_tuple_distance(&reference, &deletion, 3).unwrap(),
            428_571_429
        );

        let contained = DnaSequence::from_ascii(b"CGTACG").unwrap();
        assert_eq!(one_tuple_distance(&reference, &contained, 3).unwrap(), 0);
    }

    #[test]
    fn empty_one_tuple_boundary_is_total_and_symmetric() {
        let empty = DnaSequence::from_ascii(b"").unwrap();
        let nonempty = DnaSequence::from_ascii(b"A").unwrap();
        assert_eq!(one_tuple_distance(&empty, &empty, 3).unwrap(), 0);
        assert_eq!(
            one_tuple_distance(&empty, &nonempty, 3).unwrap(),
            DISTANCE_SCALE
        );
        assert_eq!(
            one_tuple_distance(&nonempty, &empty, 3).unwrap(),
            DISTANCE_SCALE
        );
    }

    #[test]
    fn upgma_matches_a_hand_calculated_three_leaf_tree() {
        let matrix = FixedDistanceMatrix {
            size: 3,
            values: vec![
                0,
                200_000_000,
                800_000_000,
                200_000_000,
                0,
                600_000_000,
                800_000_000,
                600_000_000,
                0,
            ],
        };
        let result = upgma(&matrix, &[0, 1, 2]).unwrap();
        let tree = to_public_tree(result.root, 3).unwrap();
        match tree.root() {
            GuideTreeNode::Branch { left, right } => {
                assert!((left.length() - 0.35).abs() < f64::EPSILON);
                assert!((right.length() - 0.25).abs() < f64::EPSILON);
                assert_eq!(left.child(), &GuideTreeNode::leaf(2));
                match right.child() {
                    GuideTreeNode::Branch { left, right } => {
                        assert!((left.length() - 0.1).abs() < f64::EPSILON);
                        assert!((right.length() - 0.1).abs() < f64::EPSILON);
                    }
                    GuideTreeNode::Leaf { .. } => panic!("expected the (0,1) branch"),
                }
            }
            GuideTreeNode::Leaf { .. } => panic!("expected a three-leaf tree"),
        }
    }

    #[test]
    fn guide_tree_update_weights_child_nodes_not_their_leaf_counts() {
        let matrix = FixedDistanceMatrix {
            size: 4,
            values: vec![
                0,
                100_000_000,
                200_000_000,
                800_000_000,
                100_000_000,
                0,
                400_000_000,
                600_000_000,
                200_000_000,
                400_000_000,
                0,
                900_000_000,
                800_000_000,
                600_000_000,
                900_000_000,
                0,
            ],
        };
        let result = upgma(&matrix, &[0, 1, 2, 3]).unwrap();
        let tree = to_public_tree(result.root, 4).unwrap();
        match tree.root() {
            GuideTreeNode::Branch { left, right } => {
                assert!((left.length() - 0.4).abs() < f64::EPSILON);
                assert!((right.length() - 0.25).abs() < f64::EPSILON);
                assert_eq!(left.child(), &GuideTreeNode::leaf(3));
            }
            GuideTreeNode::Leaf { .. } => panic!("expected a four-leaf tree"),
        }
    }

    #[test]
    fn systematic_seeds_are_length_ordered_and_deduplicated() {
        let records = [
            SequenceRecord::new("long", DnaSequence::from_ascii(b"AAAA").unwrap()).unwrap(),
            SequenceRecord::new("short", DnaSequence::from_ascii(b"A").unwrap()).unwrap(),
            SequenceRecord::new("duplicate", DnaSequence::from_ascii(b"AAAA").unwrap()).unwrap(),
            SequenceRecord::new("middle", DnaSequence::from_ascii(b"AC").unwrap()).unwrap(),
        ];
        assert_eq!(systematic_seed_indices(&records), vec![1, 2]);
    }

    #[test]
    fn mbed_enforces_the_leaf_cluster_cap_even_for_tied_vectors() {
        let records = (0..11)
            .map(|index| {
                SequenceRecord::new(
                    format!("s{index}"),
                    DnaSequence::from_ascii(b"AAAA").unwrap(),
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let result = mbed(&records, 3, 7, 3).unwrap();
        assert!(result.clusters.iter().all(|cluster| cluster.len() <= 3));
        let mut members = result
            .clusters
            .iter()
            .flatten()
            .copied()
            .collect::<Vec<_>>();
        members.sort_unstable();
        assert_eq!(members, (0..11).collect::<Vec<_>>());
        let tree = to_public_tree(result.root, records.len()).unwrap();
        assert_eq!(tree.leaf_count(), records.len());
    }

    #[test]
    fn integer_square_root_rounds_to_nearest_with_upward_ties() {
        assert_eq!(integer_sqrt_nearest(0), 0);
        assert_eq!(integer_sqrt_nearest(15), 4);
        assert_eq!(integer_sqrt_nearest(20), 4);
        assert_eq!(integer_sqrt_nearest(21), 5);
    }
}
