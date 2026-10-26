use core::mem::size_of;
use std::collections::{BTreeMap, BTreeSet};

use align_core::Sequence;
use align_msa::SequenceRecord;

use crate::muscle5::{
    Muscle5, Muscle5Alignment, Muscle5MemoryEstimate, Muscle5ProgressiveStep,
    Muscle5RefinementAttempt, estimate_returned_alignment, unordered_pair_count, validate_records,
    verify_ungapped_replay,
};
use crate::muscle5_emissions::{Muscle5Alphabet, canonical_word_index, super5_word_length};
use crate::muscle5_hmm::{HmmParameters, PreparedHmm};
use crate::muscle5_posterior::{
    PairOperation, PairPosteriorCache, PosteriorLibrary, PosteriorProvider, SparsePosterior,
    expected_error_matrix, maximum_expected_accuracy, pair_posterior,
};
use crate::muscle5_profile::{
    PairSelection, PositionedSymbol, WorkingProfile, WorkingRow, align_profiles,
    quantized_objective,
};
use crate::muscle5_rng::{CLUSTER_DOMAIN, SAMPLE_DOMAIN, SplitMix64, derived_seed};
use crate::muscle5_tree::{Muscle5Node, weighted_expected_error_upgma};
use crate::{Muscle5Error, Muscle5Mode};

const FNV_OFFSET: u64 = 0xCBF2_9CE4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;

/// Profile-profile posterior approximation selected by Super5.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Muscle5ProfileAlignmentMethod {
    /// Evaluate every cross-profile sequence pair (`PPP-pa`).
    PppPa,
    /// Evaluate the deterministic random pair sample (`PPP-ps`).
    PppPs,
}

/// A sequence removed by the Super5 redundancy pass.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5RedundancyAssignment {
    sequence_index: usize,
    representative_index: usize,
    expected_error: f64,
    pair_path_operations: usize,
}

impl Muscle5RedundancyAssignment {
    /// Returns the deferred stable input index.
    #[must_use]
    pub const fn sequence_index(&self) -> usize {
        self.sequence_index
    }

    /// Returns the representative stable input index.
    #[must_use]
    pub const fn representative_index(&self) -> usize {
        self.representative_index
    }

    /// Returns posterior-decoding expected error for the accepted pair.
    #[must_use]
    pub const fn expected_error(&self) -> f64 {
        self.expected_error
    }

    /// Returns the number of operations in the stored representative pair path.
    #[must_use]
    pub const fn pair_path_operation_count(&self) -> usize {
        self.pair_path_operations
    }
}

/// One bounded representative cluster aligned internally by PPP.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5Super5Cluster {
    source_coarse_cluster: usize,
    representative_members: Box<[usize]>,
    randomly_subdivided: bool,
    consensus: String,
    alignment_width: usize,
    raw_posterior_checksum: u64,
    consistent_posterior_checksum: u64,
    progressive_merges: usize,
    refinement_attempts: usize,
}

impl Muscle5Super5Cluster {
    /// Returns the pre-subdivision coarse-cluster index.
    #[must_use]
    pub const fn source_coarse_cluster(&self) -> usize {
        self.source_coarse_cluster
    }

    /// Returns stable indices of representative sequences in this cluster.
    #[must_use]
    pub const fn representative_members(&self) -> &[usize] {
        &self.representative_members
    }

    /// Returns whether the final size bound required random subdivision.
    #[must_use]
    pub const fn randomly_subdivided(&self) -> bool {
        self.randomly_subdivided
    }

    /// Returns the typed consensus rendered as normalized ASCII.
    #[must_use]
    pub fn consensus(&self) -> &str {
        &self.consensus
    }

    /// Returns the intra-cluster PPP alignment width.
    #[must_use]
    pub const fn alignment_width(&self) -> usize {
        self.alignment_width
    }

    /// Returns the raw intra-cluster posterior checksum.
    #[must_use]
    pub const fn raw_posterior_checksum(&self) -> u64 {
        self.raw_posterior_checksum
    }

    /// Returns the consistency-transformed intra-cluster posterior checksum.
    #[must_use]
    pub const fn consistent_posterior_checksum(&self) -> u64 {
        self.consistent_posterior_checksum
    }

    /// Returns the number of intra-cluster progressive merges.
    #[must_use]
    pub const fn progressive_merge_count(&self) -> usize {
        self.progressive_merges
    }

    /// Returns the number of intra-cluster randomized refinement attempts.
    #[must_use]
    pub const fn refinement_attempt_count(&self) -> usize {
        self.refinement_attempts
    }
}

/// One guide-tree merge of already aligned Super5 clusters.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5Super5Merge {
    left_members: Box<[usize]>,
    right_members: Box<[usize]>,
    method: Muscle5ProfileAlignmentMethod,
    evaluated_pairs: Box<[(usize, usize)]>,
    expected_accuracy: f64,
    result_width: usize,
}

impl Muscle5Super5Merge {
    /// Returns stable representative indices in the left profile.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns stable representative indices in the right profile.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }

    /// Returns whether this merge used PPP-pa or PPP-ps.
    #[must_use]
    pub const fn method(&self) -> Muscle5ProfileAlignmentMethod {
        self.method
    }

    /// Returns every canonical input pair whose posterior contributed.
    #[must_use]
    pub const fn evaluated_pairs(&self) -> &[(usize, usize)] {
        &self.evaluated_pairs
    }

    /// Returns the selected profile path's summed posterior score.
    #[must_use]
    pub const fn expected_accuracy(&self) -> f64 {
        self.expected_accuracy
    }

    /// Returns the merged representative MSA width.
    #[must_use]
    pub const fn result_width(&self) -> usize {
        self.result_width
    }
}

/// Auditable evidence for every characteristic Super5 divide-and-conquer stage.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5Super5Evidence {
    decreasing_length_order: Box<[usize]>,
    representatives: Box<[usize]>,
    redundancy_assignments: Box<[Muscle5RedundancyAssignment]>,
    clusters: Box<[Muscle5Super5Cluster]>,
    profile_merges: Box<[Muscle5Super5Merge]>,
    uclust_ee_comparisons: usize,
    direct_pair_posteriors: usize,
}

impl Muscle5Super5Evidence {
    /// Returns stable input indices after decreasing-length U-sort ordering.
    #[must_use]
    pub const fn decreasing_length_order(&self) -> &[usize] {
        &self.decreasing_length_order
    }

    /// Returns representatives retained after redundancy reduction.
    #[must_use]
    pub const fn representatives(&self) -> &[usize] {
        &self.representatives
    }

    /// Returns all deferred-to-representative assignments.
    #[must_use]
    pub const fn redundancy_assignments(&self) -> &[Muscle5RedundancyAssignment] {
        &self.redundancy_assignments
    }

    /// Returns final bounded coarse clusters in progressive leaf order.
    #[must_use]
    pub const fn clusters(&self) -> &[Muscle5Super5Cluster] {
        &self.clusters
    }

    /// Returns all inter-cluster PPP-pa/PPP-ps merges.
    #[must_use]
    pub const fn profile_merges(&self) -> &[Muscle5Super5Merge] {
        &self.profile_merges
    }

    /// Returns the number of expensive UCLUST-EE candidate comparisons.
    #[must_use]
    pub const fn uclust_ee_comparison_count(&self) -> usize {
        self.uclust_ee_comparisons
    }

    /// Returns all directly calculated pair posterior matrices across Super5 stages.
    #[must_use]
    pub const fn direct_pair_posterior_count(&self) -> usize {
        self.direct_pair_posteriors
    }
}

#[derive(Clone, Debug)]
struct DeferredSequence {
    sequence_index: usize,
    representative_index: usize,
    expected_error: f64,
    expected_accuracy: f64,
    operations: Box<[PairOperation]>,
}

#[derive(Clone, Debug)]
struct PairComparison {
    expected_error: f64,
    expected_accuracy: f64,
    operations: Box<[PairOperation]>,
}

#[derive(Clone, Debug)]
struct UclustResult {
    representatives: Vec<usize>,
    clusters: Vec<Vec<usize>>,
    assignments: Vec<DeferredSequence>,
    comparisons: usize,
}

#[derive(Clone, Debug)]
struct BoundedCluster {
    source_coarse_cluster: usize,
    members: Vec<usize>,
    randomly_subdivided: bool,
}

#[derive(Clone, Debug)]
struct PosteriorAudit {
    raw_checksum: u64,
    consistent_checksum: u64,
    raw_entries: usize,
    consistent_entries: usize,
    direct_pair_posteriors: usize,
    posterior_library_peak_bytes: usize,
    pair_hmm_peak_bytes: usize,
    consistency_peak_bytes: usize,
}

impl PosteriorAudit {
    const fn new() -> Self {
        Self {
            raw_checksum: FNV_OFFSET,
            consistent_checksum: FNV_OFFSET,
            raw_entries: 0,
            consistent_entries: 0,
            direct_pair_posteriors: 0,
            posterior_library_peak_bytes: 0,
            pair_hmm_peak_bytes: 0,
            consistency_peak_bytes: 0,
        }
    }

    fn absorb_alignment<A: Muscle5Alphabet>(
        &mut self,
        alignment: &Muscle5Alignment<A>,
    ) -> Result<(), Muscle5Error> {
        self.raw_checksum = hash_u64(self.raw_checksum, alignment.raw_posterior_checksum());
        self.consistent_checksum = hash_u64(
            self.consistent_checksum,
            alignment.consistent_posterior_checksum(),
        );
        self.raw_entries = self
            .raw_entries
            .checked_add(alignment.raw_posterior_entry_count())
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 raw posterior entry audit",
            })?;
        self.consistent_entries = self
            .consistent_entries
            .checked_add(alignment.consistent_posterior_entry_count())
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 consistent posterior entry audit",
            })?;
        let memory = alignment.memory_estimate();
        self.posterior_library_peak_bytes = self
            .posterior_library_peak_bytes
            .max(memory.posterior_library_bytes());
        self.pair_hmm_peak_bytes = self.pair_hmm_peak_bytes.max(memory.pair_hmm_peak_bytes());
        self.consistency_peak_bytes = self
            .consistency_peak_bytes
            .max(memory.consistency_peak_bytes());
        Ok(())
    }

    fn absorb_direct(
        &mut self,
        first: usize,
        second: usize,
        posterior: &SparsePosterior,
        peak_scratch_bytes: usize,
    ) -> Result<(), Muscle5Error> {
        let first_hash = u64::try_from(first).map_err(|_| Muscle5Error::ArithmeticOverflow {
            stage: "Super5 posterior index checksum",
        })?;
        let second_hash = u64::try_from(second).map_err(|_| Muscle5Error::ArithmeticOverflow {
            stage: "Super5 posterior index checksum",
        })?;
        self.raw_checksum = hash_u64(self.raw_checksum, first_hash);
        self.raw_checksum = hash_u64(self.raw_checksum, second_hash);
        self.raw_checksum = hash_u64(self.raw_checksum, posterior.checksum());
        self.consistent_checksum = hash_u64(self.consistent_checksum, first_hash);
        self.consistent_checksum = hash_u64(self.consistent_checksum, second_hash);
        self.consistent_checksum = hash_u64(self.consistent_checksum, posterior.checksum());
        self.raw_entries = self
            .raw_entries
            .checked_add(posterior.entry_count())
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 direct raw posterior entries",
            })?;
        self.consistent_entries = self
            .consistent_entries
            .checked_add(posterior.entry_count())
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 direct consistent posterior entries",
            })?;
        self.direct_pair_posteriors =
            self.direct_pair_posteriors
                .checked_add(1)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 direct posterior count",
                })?;
        self.pair_hmm_peak_bytes = self.pair_hmm_peak_bytes.max(peak_scratch_bytes);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn absorb_library(
        &mut self,
        raw_checksum: u64,
        raw_entries: usize,
        consistent_checksum: u64,
        consistent_entries: usize,
        storage_bytes: usize,
        pair_hmm_peak_bytes: usize,
        consistency_peak_bytes: usize,
    ) -> Result<(), Muscle5Error> {
        self.raw_checksum = hash_u64(self.raw_checksum, raw_checksum);
        self.consistent_checksum = hash_u64(self.consistent_checksum, consistent_checksum);
        self.raw_entries =
            self.raw_entries
                .checked_add(raw_entries)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 raw posterior library entries",
                })?;
        self.consistent_entries = self
            .consistent_entries
            .checked_add(consistent_entries)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 consistent posterior library entries",
            })?;
        self.posterior_library_peak_bytes = self.posterior_library_peak_bytes.max(storage_bytes);
        self.pair_hmm_peak_bytes = self.pair_hmm_peak_bytes.max(pair_hmm_peak_bytes);
        self.consistency_peak_bytes = self.consistency_peak_bytes.max(consistency_peak_bytes);
        Ok(())
    }
}

fn decreasing_length_order<A: Muscle5Alphabet>(records: &[SequenceRecord<A>]) -> Vec<usize> {
    let mut order = (0..records.len()).collect::<Vec<_>>();
    order.sort_by(|first, second| {
        records[*second]
            .sequence()
            .len()
            .cmp(&records[*first].sequence().len())
            .then_with(|| first.cmp(second))
    });
    order
}

fn all_kmer_counts<A: Muscle5Alphabet>(records: &[SequenceRecord<A>]) -> Vec<BTreeMap<u64, usize>> {
    records
        .iter()
        .map(|record| kmer_counts(record.sequence()))
        .collect()
}

#[allow(clippy::cast_possible_truncation)]
fn kmer_counts<A: Muscle5Alphabet>(sequence: &Sequence<A>) -> BTreeMap<u64, usize> {
    let word_length = super5_word_length::<A>();
    let base = if word_length == 3 { 20_u64 } else { 4_u64 };
    let mut counts = BTreeMap::new();
    if sequence.len() < word_length {
        return counts;
    }
    for start in 0..=sequence.len() - word_length {
        let mut word = 0_u64;
        let mut canonical = true;
        for symbol in &sequence.as_slice()[start..start + word_length] {
            let Some(index) = canonical_word_index(*symbol) else {
                canonical = false;
                break;
            };
            word = word * base + index as u64;
        }
        if canonical {
            *counts.entry(word).or_insert(0) += 1;
        }
    }
    counts
}

fn shared_kmer_count(first: &BTreeMap<u64, usize>, second: &BTreeMap<u64, usize>) -> usize {
    let (smaller, larger) = if first.len() <= second.len() {
        (first, second)
    } else {
        (second, first)
    };
    smaller
        .iter()
        .map(|(word, count)| count.min(larger.get(word).unwrap_or(&0)))
        .sum()
}

fn compare_pair<A: Muscle5Alphabet>(
    representative: usize,
    query: usize,
    records: &[SequenceRecord<A>],
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    audit: &mut PosteriorAudit,
    cache: &mut PairPosteriorCache<A>,
) -> Result<PairComparison, Muscle5Error> {
    let calculation = cache.pair_posterior(
        records[representative].sequence(),
        records[query].sequence(),
        hmm,
        cutoff,
    )?;
    if !calculation.log_partition.is_finite() {
        return Err(Muscle5Error::NonFiniteComputation {
            stage: "Super5 UCLUST pair partition",
        });
    }
    audit.absorb_direct(
        representative,
        query,
        &calculation.posterior,
        calculation.peak_scratch_bytes,
    )?;
    let decoded = maximum_expected_accuracy(&calculation.posterior)?;
    Ok(PairComparison {
        expected_error: decoded.expected_error(),
        expected_accuracy: decoded.expected_accuracy(),
        operations: decoded.operations().to_vec().into_boxed_slice(),
    })
}

#[allow(clippy::too_many_arguments)]
fn uclust_ee<A: Muscle5Alphabet>(
    ordered_indices: &[usize],
    records: &[SequenceRecord<A>],
    kmers: &[BTreeMap<u64, usize>],
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    expected_error_threshold: f64,
    maximum_rejections: usize,
    audit: &mut PosteriorAudit,
    cache: &mut PairPosteriorCache<A>,
) -> Result<UclustResult, Muscle5Error> {
    let mut representatives: Vec<usize> = Vec::new();
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut assignments = Vec::new();
    let mut comparisons = 0_usize;
    for &query in ordered_indices {
        let mut ranked = representatives
            .iter()
            .copied()
            .enumerate()
            .map(|(cluster, representative)| {
                (
                    representative,
                    cluster,
                    shared_kmer_count(&kmers[query], &kmers[representative]),
                )
            })
            .collect::<Vec<_>>();
        ranked.sort_by(|first, second| second.2.cmp(&first.2).then_with(|| first.0.cmp(&second.0)));
        let mut matched = None;
        for (representative, cluster, _) in ranked.into_iter().take(maximum_rejections) {
            let comparison =
                compare_pair(representative, query, records, hmm, cutoff, audit, cache)?;
            comparisons = comparisons
                .checked_add(1)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "UCLUST-EE comparison count",
                })?;
            if comparison.expected_error < expected_error_threshold {
                matched = Some((representative, cluster, comparison));
                break;
            }
        }
        if let Some((representative, cluster, comparison)) = matched {
            clusters[cluster].push(query);
            assignments.push(DeferredSequence {
                sequence_index: query,
                representative_index: representative,
                expected_error: comparison.expected_error,
                expected_accuracy: comparison.expected_accuracy,
                operations: comparison.operations,
            });
        } else {
            representatives.push(query);
            clusters.push(vec![query]);
        }
    }
    Ok(UclustResult {
        representatives,
        clusters,
        assignments,
        comparisons,
    })
}

#[allow(clippy::too_many_arguments)]
fn bound_coarse_clusters<A: Muscle5Alphabet>(
    coarse_clusters: &[Vec<usize>],
    records: &[SequenceRecord<A>],
    kmers: &[BTreeMap<u64, usize>],
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    subdivision_threshold: f64,
    maximum_size: usize,
    maximum_rejections: usize,
    rng: &mut SplitMix64,
    audit: &mut PosteriorAudit,
    cache: &mut PairPosteriorCache<A>,
) -> Result<(Vec<BoundedCluster>, usize), Muscle5Error> {
    let mut bounded = Vec::new();
    let mut comparisons = 0_usize;
    for (source_coarse_cluster, cluster) in coarse_clusters.iter().enumerate() {
        let subdivisions = if cluster.len() > maximum_size {
            let mut ordered = cluster.clone();
            ordered.sort_by(|first, second| {
                records[*second]
                    .sequence()
                    .len()
                    .cmp(&records[*first].sequence().len())
                    .then_with(|| first.cmp(second))
            });
            let result = uclust_ee(
                &ordered,
                records,
                kmers,
                hmm,
                cutoff,
                subdivision_threshold,
                maximum_rejections,
                audit,
                cache,
            )?;
            comparisons = comparisons.checked_add(result.comparisons).ok_or(
                Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 subdivision comparison count",
                },
            )?;
            result.clusters
        } else {
            vec![cluster.clone()]
        };
        for subdivision in subdivisions {
            if subdivision.len() <= maximum_size {
                bounded.push(BoundedCluster {
                    source_coarse_cluster,
                    members: subdivision,
                    randomly_subdivided: false,
                });
                continue;
            }
            let mut shuffled = subdivision;
            for index in (1..shuffled.len()).rev() {
                let selected = rng.index(index + 1);
                shuffled.swap(index, selected);
            }
            for chunk in shuffled.chunks(maximum_size) {
                bounded.push(BoundedCluster {
                    source_coarse_cluster,
                    members: chunk.to_vec(),
                    randomly_subdivided: true,
                });
            }
        }
    }
    Ok((bounded, comparisons))
}

#[derive(Clone, Debug)]
struct ClusterRun<A: Muscle5Alphabet> {
    profile: WorkingProfile<A>,
    consensus: Sequence<A>,
    evidence: Muscle5Super5Cluster,
    progressive_steps: Vec<Muscle5ProgressiveStep>,
    refinement_attempts: Vec<Muscle5RefinementAttempt>,
    fixed_objective: i128,
    objective_pairs: usize,
    distance_and_tree_bytes: usize,
    profile_dp_peak_bytes: usize,
}

fn align_bounded_cluster<A: Muscle5Alphabet>(
    cluster: &BoundedCluster,
    records: &[SequenceRecord<A>],
    aligner: &Muscle5<A>,
    audit: &mut PosteriorAudit,
    cache: &PairPosteriorCache<A>,
) -> Result<ClusterRun<A>, Muscle5Error> {
    let local_records = cluster
        .members
        .iter()
        .map(|index| records[*index].clone())
        .collect::<Vec<_>>();
    let local_parameters = aligner.parameters().clone().with_mode(Muscle5Mode::Align);
    let local =
        Muscle5::<A>::with_parameters(local_parameters)?.align_reusing(&local_records, cache)?;
    audit.absorb_alignment(&local)?;
    let profile = WorkingProfile::from_alignment(local.alignment(), &cluster.members)?;
    let consensus = consensus_sequence(&profile)?;
    let progressive_steps = local
        .progressive_steps()
        .iter()
        .map(|step| Muscle5ProgressiveStep {
            left_members: map_local_members(step.left_members(), &cluster.members),
            right_members: map_local_members(step.right_members(), &cluster.members),
            result_width: step.result_width(),
            expected_accuracy: step.expected_accuracy(),
            operation_count: step.operation_count(),
        })
        .collect::<Vec<_>>();
    let refinement_attempts = local
        .refinement_attempts()
        .iter()
        .map(|attempt| Muscle5RefinementAttempt {
            iteration: attempt.iteration(),
            first_partition: map_local_members(attempt.first_partition(), &cluster.members),
            fixed_score_before: attempt.fixed_score_before,
            fixed_score_after: attempt.fixed_score_after,
            accepted: attempt.accepted(),
            candidate_width: attempt.candidate_width(),
            profile_operations: attempt.profile_operation_count(),
        })
        .collect::<Vec<_>>();
    let evidence = Muscle5Super5Cluster {
        source_coarse_cluster: cluster.source_coarse_cluster,
        representative_members: cluster.members.clone().into_boxed_slice(),
        randomly_subdivided: cluster.randomly_subdivided,
        consensus: consensus.to_ascii_string(),
        alignment_width: local.alignment().width(),
        raw_posterior_checksum: local.raw_posterior_checksum(),
        consistent_posterior_checksum: local.consistent_posterior_checksum(),
        progressive_merges: local.progressive_steps().len(),
        refinement_attempts: local.refinement_attempts().len(),
    };
    Ok(ClusterRun {
        profile,
        consensus,
        evidence,
        progressive_steps,
        refinement_attempts,
        fixed_objective: quantized_objective(local.final_objective_score())?,
        objective_pairs: unordered_pair_count(cluster.members.len())?,
        distance_and_tree_bytes: local.memory_estimate().distance_and_tree_bytes(),
        profile_dp_peak_bytes: local.memory_estimate().profile_dp_peak_bytes(),
    })
}

fn map_local_members(local: &[usize], global: &[usize]) -> Box<[usize]> {
    let mut mapped = local.iter().map(|index| global[*index]).collect::<Vec<_>>();
    mapped.sort_unstable();
    mapped.into_boxed_slice()
}

fn consensus_sequence<A: Muscle5Alphabet>(
    profile: &WorkingProfile<A>,
) -> Result<Sequence<A>, Muscle5Error> {
    let width = profile.width();
    let mut consensus = Vec::with_capacity(width);
    for column in 0..width {
        let mut counts = vec![0_usize; A::SYMBOLS.len()];
        let mut gaps = 0_usize;
        for row in &profile.rows {
            if let Some(cell) = row.cells[column] {
                counts[cell.symbol.alphabet_index()] += 1;
            } else {
                gaps += 1;
            }
        }
        let selected = counts
            .iter()
            .enumerate()
            .max_by(|first, second| first.1.cmp(second.1).then_with(|| second.0.cmp(&first.0)));
        let Some((symbol_index, count)) = selected else {
            return Err(Muscle5Error::InternalInvariant {
                message: "Super5 consensus alphabet is empty",
            });
        };
        if *count <= gaps {
            continue;
        }
        consensus.push(A::SYMBOLS[symbol_index]);
    }
    Sequence::from_ascii(consensus).map_err(|_| Muscle5Error::InternalInvariant {
        message: "Super5 consensus emitted an invalid typed symbol",
    })
}

#[derive(Clone, Debug)]
struct SelectedPosteriorStore {
    lengths: Box<[usize]>,
    matrices: BTreeMap<(usize, usize), SparsePosterior>,
}

impl SelectedPosteriorStore {
    fn build<A: Muscle5Alphabet>(
        records: &[SequenceRecord<A>],
        pairs: &[(usize, usize)],
        hmm: &PreparedHmm<A>,
        cutoff: f64,
        audit: &mut PosteriorAudit,
        cache: &PairPosteriorCache<A>,
    ) -> Result<Self, Muscle5Error> {
        let lengths = records
            .iter()
            .map(|record| record.sequence().len())
            .collect::<Vec<_>>()
            .into_boxed_slice();
        let mut matrices = BTreeMap::new();
        let mut storage_bytes = lengths.len().checked_mul(size_of::<usize>()).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "Super5 selected posterior lengths",
            },
        )?;
        for &(first, second) in pairs {
            if first >= records.len() || second >= records.len() || first >= second {
                return Err(Muscle5Error::InternalInvariant {
                    message: "Super5 selected posterior pair is not canonical",
                });
            }
            let (first_sequence, second_sequence) =
                (records[first].sequence(), records[second].sequence());
            let calculation = match cache.stored(first_sequence, second_sequence, hmm, cutoff) {
                Some(calculation) => calculation,
                None => pair_posterior(first_sequence, second_sequence, hmm, cutoff)?,
            };
            if !calculation.log_partition.is_finite() {
                return Err(Muscle5Error::NonFiniteComputation {
                    stage: "Super5 profile-pair partition",
                });
            }
            audit.absorb_direct(
                first,
                second,
                &calculation.posterior,
                calculation.peak_scratch_bytes,
            )?;
            let reverse = calculation.posterior.transpose()?;
            storage_bytes = storage_bytes
                .checked_add(calculation.posterior.storage_bytes())
                .and_then(|bytes| bytes.checked_add(reverse.storage_bytes()))
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 selected posterior storage",
                })?;
            matrices.insert((first, second), calculation.posterior);
            matrices.insert((second, first), reverse);
        }
        audit.posterior_library_peak_bytes = audit.posterior_library_peak_bytes.max(storage_bytes);
        Ok(Self { lengths, matrices })
    }
}

impl PosteriorProvider for SelectedPosteriorStore {
    fn posterior(&self, first: usize, second: usize) -> Option<&SparsePosterior> {
        self.matrices.get(&(first, second))
    }

    fn sequence_length(&self, sequence_index: usize) -> Option<usize> {
        self.lengths.get(sequence_index).copied()
    }
}

fn selected_cross_pairs<A: Muscle5Alphabet>(
    left: &WorkingProfile<A>,
    right: &WorkingProfile<A>,
    limit: usize,
    rng: &mut SplitMix64,
) -> Result<(Muscle5ProfileAlignmentMethod, Vec<(usize, usize)>), Muscle5Error> {
    let left_members = left.members();
    let right_members = right.members();
    let product = left_members.len().checked_mul(right_members.len()).ok_or(
        Muscle5Error::ArithmeticOverflow {
            stage: "Super5 cross-profile pair count",
        },
    )?;
    let selected_indices = if product <= limit {
        None
    } else {
        let mut selected = BTreeSet::new();
        while selected.len() < limit {
            selected.insert(rng.index(product));
        }
        Some(selected)
    };
    let mut pairs = Vec::with_capacity(product.min(limit));
    if let Some(selected) = selected_indices {
        for pair_index in selected {
            let left_index = pair_index / right_members.len();
            let right_index = pair_index % right_members.len();
            pairs.push(canonical_pair(
                left_members[left_index],
                right_members[right_index],
            ));
        }
    } else {
        for &first in &left_members {
            for &second in &right_members {
                pairs.push(canonical_pair(first, second));
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    if pairs.len() != product.min(limit) {
        return Err(Muscle5Error::InternalInvariant {
            message: "Super5 cross-profile sampling did not produce the requested pair count",
        });
    }
    let method = if product <= limit {
        Muscle5ProfileAlignmentMethod::PppPa
    } else {
        Muscle5ProfileAlignmentMethod::PppPs
    };
    Ok((method, pairs))
}

const fn canonical_pair(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

#[allow(clippy::too_many_arguments)]
fn merge_cluster_tree<A: Muscle5Alphabet>(
    node: &Muscle5Node,
    cluster_profiles: &[WorkingProfile<A>],
    records: &[SequenceRecord<A>],
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    pair_limit: usize,
    rng: &mut SplitMix64,
    audit: &mut PosteriorAudit,
    merges: &mut Vec<Muscle5Super5Merge>,
    progressive_steps: &mut Vec<Muscle5ProgressiveStep>,
    profile_dp_peak_bytes: &mut usize,
    fixed_objective: &mut i128,
    objective_pairs: &mut usize,
    cache: &PairPosteriorCache<A>,
) -> Result<WorkingProfile<A>, Muscle5Error> {
    match node {
        Muscle5Node::Leaf { sequence_index } => cluster_profiles
            .get(*sequence_index)
            .cloned()
            .ok_or(Muscle5Error::InternalInvariant {
                message: "Super5 guide-tree leaf is outside cluster profiles",
            }),
        Muscle5Node::Branch { left, right, .. } => {
            let left_profile = merge_cluster_tree(
                left,
                cluster_profiles,
                records,
                hmm,
                cutoff,
                pair_limit,
                rng,
                audit,
                merges,
                progressive_steps,
                profile_dp_peak_bytes,
                fixed_objective,
                objective_pairs,
                cache,
            )?;
            let right_profile = merge_cluster_tree(
                right,
                cluster_profiles,
                records,
                hmm,
                cutoff,
                pair_limit,
                rng,
                audit,
                merges,
                progressive_steps,
                profile_dp_peak_bytes,
                fixed_objective,
                objective_pairs,
                cache,
            )?;
            let left_members = left_profile.members().into_boxed_slice();
            let right_members = right_profile.members().into_boxed_slice();
            let (method, pairs) =
                selected_cross_pairs(&left_profile, &right_profile, pair_limit, rng)?;
            let store = SelectedPosteriorStore::build(records, &pairs, hmm, cutoff, audit, cache)?;
            let profile_merge = align_profiles(
                left_profile,
                right_profile,
                &store,
                PairSelection::Selected(&pairs),
            )?;
            *profile_dp_peak_bytes = (*profile_dp_peak_bytes).max(profile_merge.workspace_bytes);
            *fixed_objective = fixed_objective
                .checked_add(quantized_objective(profile_merge.expected_accuracy)?)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 evaluated objective",
                })?;
            *objective_pairs = objective_pairs.checked_add(pairs.len()).ok_or(
                Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 evaluated objective pair count",
                },
            )?;
            let result_width = profile_merge.profile.width();
            let operation_count = profile_merge.operations.len();
            progressive_steps.push(Muscle5ProgressiveStep {
                left_members: left_members.clone(),
                right_members: right_members.clone(),
                result_width,
                expected_accuracy: profile_merge.expected_accuracy,
                operation_count,
            });
            merges.push(Muscle5Super5Merge {
                left_members,
                right_members,
                method,
                evaluated_pairs: pairs.into_boxed_slice(),
                expected_accuracy: profile_merge.expected_accuracy,
                result_width,
            });
            Ok(profile_merge.profile)
        }
    }
}

#[allow(clippy::too_many_lines)]
fn reinsert_deferred<A: Muscle5Alphabet>(
    profile: &mut WorkingProfile<A>,
    deferred: &DeferredSequence,
    records: &[SequenceRecord<A>],
) -> Result<(), Muscle5Error> {
    if profile
        .rows
        .iter()
        .any(|row| row.input_index == deferred.sequence_index)
    {
        return Err(Muscle5Error::InternalInvariant {
            message: "Super5 attempted to insert a row twice",
        });
    }
    let representative_length = records[deferred.representative_index].sequence().len();
    let sequence = records[deferred.sequence_index].sequence();
    let mut anchored = vec![None; representative_length];
    let mut insertions = vec![Vec::new(); representative_length + 1];
    let mut representative_position = 0_usize;
    let mut sequence_position = 0_usize;
    for operation in &deferred.operations {
        match operation {
            PairOperation::Match => {
                let symbol =
                    sequence
                        .get(sequence_position)
                        .ok_or(Muscle5Error::InternalInvariant {
                            message: "redundant pair path exceeded the deferred sequence",
                        })?;
                let slot = anchored.get_mut(representative_position).ok_or(
                    Muscle5Error::InternalInvariant {
                        message: "redundant pair path exceeded its representative",
                    },
                )?;
                *slot = Some(PositionedSymbol {
                    symbol,
                    position: sequence_position,
                });
                representative_position += 1;
                sequence_position += 1;
            }
            PairOperation::ConsumeFirst => representative_position += 1,
            PairOperation::ConsumeSecond => {
                let symbol =
                    sequence
                        .get(sequence_position)
                        .ok_or(Muscle5Error::InternalInvariant {
                            message: "redundant insertion exceeded the deferred sequence",
                        })?;
                let interval = insertions.get_mut(representative_position).ok_or(
                    Muscle5Error::InternalInvariant {
                        message: "redundant insertion exceeded its representative",
                    },
                )?;
                interval.push(PositionedSymbol {
                    symbol,
                    position: sequence_position,
                });
                sequence_position += 1;
            }
        }
    }
    if representative_position != representative_length || sequence_position != sequence.len() {
        return Err(Muscle5Error::InternalInvariant {
            message: "redundant pair path did not consume both sequences",
        });
    }

    let representative_row = profile
        .rows
        .iter()
        .position(|row| row.input_index == deferred.representative_index)
        .ok_or(Muscle5Error::InternalInvariant {
            message: "Super5 representative row is absent during reinsertion",
        })?;
    let mut representative_columns = vec![None; representative_length];
    for (column, cell) in profile.rows[representative_row].cells.iter().enumerate() {
        if let Some(cell) = cell {
            let slot = representative_columns.get_mut(cell.position).ok_or(
                Muscle5Error::InternalInvariant {
                    message: "representative MSA position exceeds its raw sequence",
                },
            )?;
            *slot = Some(column);
        }
    }
    if representative_columns.iter().any(Option::is_none) {
        return Err(Muscle5Error::InternalInvariant {
            message: "representative MSA omitted a biological residue",
        });
    }
    let mut representative_columns = representative_columns
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(Muscle5Error::InternalInvariant {
            message: "representative MSA omitted a biological residue",
        })?;
    let mut placed_insertions = vec![Vec::new(); representative_length + 1];
    for interval in 0..=representative_length {
        let start = if interval == 0 {
            0
        } else {
            representative_columns[interval - 1] + 1
        };
        let end = if interval == representative_length {
            profile.width()
        } else {
            representative_columns[interval]
        };
        let required = insertions[interval].len();
        let available = end.saturating_sub(start);
        let insertion_end = end.checked_add(required.saturating_sub(available)).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "Super5 insertion-column placement",
            },
        )?;
        for insertion_column in end..insertion_end {
            for row in &mut profile.rows {
                row.cells.insert(insertion_column, None);
            }
            for column in &mut representative_columns {
                if *column >= insertion_column {
                    *column += 1;
                }
            }
        }
        placed_insertions[interval] = (start..start + required).collect();
    }
    let mut cells = vec![None; profile.width()];
    for (position, cell) in anchored.into_iter().enumerate() {
        if let Some(cell) = cell {
            cells[representative_columns[position]] = Some(cell);
        }
    }
    for (interval, symbols) in insertions.into_iter().enumerate() {
        for (cell, column) in symbols
            .into_iter()
            .zip(placed_insertions[interval].iter().copied())
        {
            cells[column] = Some(cell);
        }
    }
    profile.rows.push(WorkingRow {
        input_index: deferred.sequence_index,
        cells,
    });
    profile.rows.sort_by_key(|row| row.input_index);
    Ok(())
}

pub(crate) fn run_super5<A: Muscle5Alphabet>(
    aligner: &Muscle5<A>,
    records: &[SequenceRecord<A>],
) -> Result<Muscle5Alignment<A>, Muscle5Error> {
    validate_records(records)?;
    aligner.parameters().validate()?;
    let parameters = aligner.parameters();
    let hmm_parameters = HmmParameters::published::<A>(
        parameters.perturbation_seed(),
        parameters.perturbation_amplitude(),
    )?;
    let hmm = hmm_parameters.prepare::<A>();
    // Later stages request many pairs that UCLUST already compared; reusing
    // those posteriors is exact, because each is a pure function of its pair.
    let mut cache = PairPosteriorCache::new(&hmm, parameters.posterior_cutoff());
    run_super5_using(aligner, records, &hmm, &mut cache)
}

/// [`run_super5`] after validation, with the run's prepared HMM and the cache
/// that UCLUST fills and later stages read.
#[allow(clippy::too_many_lines)]
fn run_super5_using<A: Muscle5Alphabet>(
    aligner: &Muscle5<A>,
    records: &[SequenceRecord<A>],
    hmm: &PreparedHmm<A>,
    cache: &mut PairPosteriorCache<A>,
) -> Result<Muscle5Alignment<A>, Muscle5Error> {
    let parameters = aligner.parameters();
    let decreasing_length_order = decreasing_length_order(records);
    let kmers = all_kmer_counts(records);
    let mut audit = PosteriorAudit::new();

    let redundancy = uclust_ee(
        &decreasing_length_order,
        records,
        &kmers,
        hmm,
        parameters.posterior_cutoff(),
        parameters.redundancy_ee_threshold(),
        parameters.maximum_uclust_rejections(),
        &mut audit,
        cache,
    )?;
    let coarse = uclust_ee(
        &redundancy.representatives,
        records,
        &kmers,
        hmm,
        parameters.posterior_cutoff(),
        parameters.coarse_ee_threshold(),
        parameters.maximum_uclust_rejections(),
        &mut audit,
        cache,
    )?;
    let mut cluster_rng =
        SplitMix64::new(derived_seed(parameters.perturbation_seed(), CLUSTER_DOMAIN));
    let (bounded_clusters, subdivision_comparisons) = bound_coarse_clusters(
        &coarse.clusters,
        records,
        &kmers,
        hmm,
        parameters.posterior_cutoff(),
        parameters.subdivision_ee_threshold(),
        parameters.maximum_cluster_size(),
        parameters.maximum_uclust_rejections(),
        &mut cluster_rng,
        &mut audit,
        cache,
    )?;
    if bounded_clusters.is_empty() {
        return Err(Muscle5Error::InternalInvariant {
            message: "Super5 clustering produced no bounded clusters",
        });
    }

    let mut cluster_profiles = Vec::with_capacity(bounded_clusters.len());
    let mut consensus_records = Vec::with_capacity(bounded_clusters.len());
    let mut cluster_evidence = Vec::with_capacity(bounded_clusters.len());
    let mut progressive_steps = Vec::new();
    let mut refinement_attempts = Vec::new();
    let mut fixed_objective = 0_i128;
    let mut objective_pairs = 0_usize;
    let mut profile_dp_peak_bytes = 0_usize;
    let mut distance_and_tree_bytes = 0_usize;
    for (cluster_index, cluster) in bounded_clusters.iter().enumerate() {
        let run = align_bounded_cluster(cluster, records, aligner, &mut audit, cache)?;
        fixed_objective = fixed_objective.checked_add(run.fixed_objective).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "Super5 intra-cluster objective",
            },
        )?;
        objective_pairs = objective_pairs.checked_add(run.objective_pairs).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "Super5 intra-cluster objective pairs",
            },
        )?;
        profile_dp_peak_bytes = profile_dp_peak_bytes.max(run.profile_dp_peak_bytes);
        distance_and_tree_bytes = distance_and_tree_bytes.max(run.distance_and_tree_bytes);
        progressive_steps.extend(run.progressive_steps);
        refinement_attempts.extend(run.refinement_attempts);
        consensus_records.push(
            SequenceRecord::new(
                format!("super5_cluster_{cluster_index}"),
                run.consensus.clone(),
            )
            .map_err(Muscle5Error::from)?,
        );
        cluster_profiles.push(run.profile);
        cluster_evidence.push(run.evidence);
    }

    let mut consensus_library = PosteriorLibrary::build(
        &consensus_records,
        hmm,
        parameters.posterior_cutoff(),
        Some(cache),
    )?;
    let consensus_raw_checksum = consensus_library.checksum();
    let consensus_raw_entries = consensus_library.entry_count();
    let consensus_raw_storage = consensus_library.storage_bytes();
    consensus_library.transform_consistency(
        parameters.consistency_iterations(),
        parameters.posterior_cutoff(),
    )?;
    audit.absorb_library(
        consensus_raw_checksum,
        consensus_raw_entries,
        consensus_library.checksum(),
        consensus_library.entry_count(),
        consensus_raw_storage.max(consensus_library.storage_bytes()),
        consensus_library.pair_hmm_peak_bytes(),
        consensus_library.consistency_peak_bytes(),
    )?;
    let consensus_distances = expected_error_matrix(&consensus_library)?;
    let mut representative_population = vec![1_usize; records.len()];
    for assignment in &redundancy.assignments {
        representative_population[assignment.representative_index] = representative_population
            [assignment.representative_index]
            .checked_add(1)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 representative population",
            })?;
    }
    let cluster_weights = bounded_clusters
        .iter()
        .map(|cluster| {
            cluster.members.iter().try_fold(0_usize, |total, member| {
                total.checked_add(representative_population[*member]).ok_or(
                    Muscle5Error::ArithmeticOverflow {
                        stage: "Super5 biased guide-tree leaf weight",
                    },
                )
            })
        })
        .collect::<Result<Vec<_>, Muscle5Error>>()?;
    let original_tree = weighted_expected_error_upgma(&consensus_distances, &cluster_weights)?;
    let cluster_tree = original_tree.permuted(parameters.guide_tree_permutation())?;
    let guide_tree = cluster_tree.to_public()?;
    distance_and_tree_bytes = distance_and_tree_bytes.max(
        consensus_distances
            .as_flat_slice()
            .len()
            .checked_mul(size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(original_tree.workspace_bytes()))
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 consensus distance and tree storage",
            })?,
    );

    let mut sampling_rng =
        SplitMix64::new(derived_seed(parameters.perturbation_seed(), SAMPLE_DOMAIN));
    let mut profile_merges = Vec::new();
    let mut representative_profile = merge_cluster_tree(
        cluster_tree.root(),
        &cluster_profiles,
        records,
        hmm,
        parameters.posterior_cutoff(),
        parameters.ppp_pair_limit(),
        &mut sampling_rng,
        &mut audit,
        &mut profile_merges,
        &mut progressive_steps,
        &mut profile_dp_peak_bytes,
        &mut fixed_objective,
        &mut objective_pairs,
        cache,
    )?;

    let mut deferred = redundancy.assignments;
    deferred.sort_by_key(|assignment| assignment.sequence_index);
    for assignment in &deferred {
        reinsert_deferred(&mut representative_profile, assignment, records)?;
        fixed_objective = fixed_objective
            .checked_add(quantized_objective(assignment.expected_accuracy)?)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 redundancy objective",
            })?;
        objective_pairs =
            objective_pairs
                .checked_add(1)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "Super5 redundancy objective pairs",
                })?;
    }
    let alignment = representative_profile.to_alignment(records)?;
    verify_ungapped_replay(&alignment, records)?;
    let returned_alignment_bytes = estimate_returned_alignment(&alignment)?;
    let super5_auxiliary_bytes = estimate_super5_auxiliary(
        &kmers,
        &bounded_clusters,
        &cluster_profiles,
        &profile_merges,
    )?;
    let phase_peak = audit
        .pair_hmm_peak_bytes
        .max(audit.consistency_peak_bytes)
        .max(profile_dp_peak_bytes);
    let estimated_peak_bytes = hmm
        .storage_bytes()
        .checked_add(audit.posterior_library_peak_bytes)
        .and_then(|bytes| bytes.checked_add(distance_and_tree_bytes))
        .and_then(|bytes| bytes.checked_add(super5_auxiliary_bytes))
        .and_then(|bytes| bytes.checked_add(returned_alignment_bytes))
        .and_then(|bytes| bytes.checked_add(phase_peak))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "Super5 peak memory estimate",
        })?;
    let uclust_ee_comparisons = redundancy
        .comparisons
        .checked_add(coarse.comparisons)
        .and_then(|count| count.checked_add(subdivision_comparisons))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "Super5 UCLUST-EE comparisons",
        })?;
    let public_assignments = deferred
        .iter()
        .map(|assignment| Muscle5RedundancyAssignment {
            sequence_index: assignment.sequence_index,
            representative_index: assignment.representative_index,
            expected_error: assignment.expected_error,
            pair_path_operations: assignment.operations.len(),
        })
        .collect::<Vec<_>>()
        .into_boxed_slice();
    let evidence = Muscle5Super5Evidence {
        decreasing_length_order: decreasing_length_order.into_boxed_slice(),
        representatives: redundancy.representatives.into_boxed_slice(),
        redundancy_assignments: public_assignments,
        clusters: cluster_evidence.into_boxed_slice(),
        profile_merges: profile_merges.into_boxed_slice(),
        uclust_ee_comparisons,
        direct_pair_posteriors: audit.direct_pair_posteriors,
    };
    let possible_objective_pairs = unordered_pair_count(records.len())?;
    if objective_pairs > possible_objective_pairs {
        return Err(Muscle5Error::InternalInvariant {
            message: "Super5 evaluated a final sequence pair more than once",
        });
    }
    Ok(Muscle5Alignment {
        alignment,
        mode: Muscle5Mode::Super5,
        perturbation_seed: parameters.perturbation_seed(),
        guide_tree_permutation: parameters.guide_tree_permutation(),
        expected_error_distances: consensus_distances,
        guide_tree,
        raw_posterior_checksum: audit.raw_checksum,
        consistent_posterior_checksum: audit.consistent_checksum,
        raw_posterior_entries: audit.raw_entries,
        consistent_posterior_entries: audit.consistent_entries,
        progressive_steps: progressive_steps.into_boxed_slice(),
        refinement_attempts: refinement_attempts.into_boxed_slice(),
        fixed_final_objective: fixed_objective,
        evaluated_objective_pairs: objective_pairs,
        possible_objective_pairs,
        super5_evidence: Some(evidence),
        memory: Muscle5MemoryEstimate {
            hmm_bytes: hmm.storage_bytes(),
            posterior_library_bytes: audit.posterior_library_peak_bytes,
            pair_hmm_peak_bytes: audit.pair_hmm_peak_bytes,
            consistency_peak_bytes: audit.consistency_peak_bytes,
            distance_and_tree_bytes,
            profile_dp_peak_bytes,
            super5_auxiliary_bytes,
            returned_alignment_bytes,
            estimated_peak_bytes,
        },
    })
}

fn estimate_super5_auxiliary<A: Muscle5Alphabet>(
    kmers: &[BTreeMap<u64, usize>],
    clusters: &[BoundedCluster],
    profiles: &[WorkingProfile<A>],
    merges: &[Muscle5Super5Merge],
) -> Result<usize, Muscle5Error> {
    let kmer_entries = kmers.iter().map(BTreeMap::len).sum::<usize>();
    let kmer_bytes = kmer_entries.checked_mul(size_of::<(u64, usize)>()).ok_or(
        Muscle5Error::ArithmeticOverflow {
            stage: "Super5 k-mer storage estimate",
        },
    )?;
    let cluster_members = clusters
        .iter()
        .map(|cluster| cluster.members.len())
        .sum::<usize>();
    let cluster_bytes = cluster_members.checked_mul(size_of::<usize>()).ok_or(
        Muscle5Error::ArithmeticOverflow {
            stage: "Super5 cluster storage estimate",
        },
    )?;
    let profile_cells = profiles.iter().try_fold(0_usize, |total, profile| {
        profile
            .rows
            .len()
            .checked_mul(profile.width())
            .and_then(|cells| total.checked_add(cells))
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "Super5 retained profile cell count",
            })
    })?;
    let profile_bytes = profile_cells
        .checked_mul(size_of::<Option<PositionedSymbol<A>>>())
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "Super5 retained profile storage estimate",
        })?;
    let sampled_pairs = merges
        .iter()
        .map(|merge| merge.evaluated_pairs.len())
        .sum::<usize>();
    let sampled_bytes = sampled_pairs
        .checked_mul(size_of::<(usize, usize)>())
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "Super5 sampled-pair storage estimate",
        })?;
    kmer_bytes
        .checked_add(cluster_bytes)
        .and_then(|bytes| bytes.checked_add(profile_bytes))
        .and_then(|bytes| bytes.checked_add(sampled_bytes))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "Super5 auxiliary memory estimate",
        })
}

fn hash_u64(mut hash: u64, value: u64) -> u64 {
    for byte in value.to_le_bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::SequenceRecord;

    use super::{
        DeferredSequence, Muscle5ProfileAlignmentMethod, PosteriorAudit, SplitMix64,
        WorkingProfile, WorkingRow, all_kmer_counts, bound_coarse_clusters, consensus_sequence,
        decreasing_length_order, kmer_counts, reinsert_deferred, run_super5_using,
        selected_cross_pairs, uclust_ee,
    };
    use crate::muscle5_hmm::HmmParameters;
    use crate::muscle5_posterior::{PairOperation, PairPosteriorCache};
    use crate::muscle5_profile::{PositionedSymbol, WorkingProfile as Profile};
    use crate::{Muscle5, Muscle5Parameters};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn fast_super5() -> Muscle5<Dna> {
        let parameters = Muscle5Parameters::published_super5()
            .with_iterations(0, 0)
            .with_super5_limits(0.0, 0.30, 0.10, 2, 16, 2)
            .unwrap();
        Muscle5::with_parameters(parameters).unwrap()
    }

    #[test]
    fn canonical_kmers_preserve_multiplicity_and_ambiguity_breaks_words() {
        let repeated = DnaSequence::from_ascii("AAAAAAA").unwrap();
        let counts = kmer_counts(&repeated);
        assert_eq!(counts.values().copied().collect::<Vec<_>>(), vec![2]);
        let ambiguous = DnaSequence::from_ascii("AAANAAA").unwrap();
        assert!(kmer_counts(&ambiguous).is_empty());
        let protein = ProteinSequence::from_ascii("AAAA").unwrap();
        assert_eq!(kmer_counts(&protein).values().next(), Some(&2));
    }

    #[test]
    fn decreasing_length_order_uses_stable_input_index_ties() {
        let records = [dna("a", "A"), dna("b", "AAA"), dna("c", "AAA")];
        assert_eq!(decreasing_length_order(&records), vec![1, 2, 0]);
    }

    #[test]
    fn uclust_ee_defers_an_identical_sequence_to_the_first_representative() {
        let records = [dna("first", "ACGTAC"), dna("same", "ACGTAC")];
        let kmers = all_kmer_counts(&records);
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let mut audit = PosteriorAudit::new();
        let mut cache = PairPosteriorCache::new(&hmm, 0.01);
        let result = uclust_ee(
            &[0, 1],
            &records,
            &kmers,
            &hmm,
            0.01,
            0.01,
            16,
            &mut audit,
            &mut cache,
        )
        .unwrap();
        assert_eq!(result.representatives, vec![0]);
        assert_eq!(result.clusters, vec![vec![0, 1]]);
        assert_eq!(result.assignments[0].representative_index, 0);
        assert!(result.assignments[0].expected_error < 0.01);
        assert_eq!(result.comparisons, 1);
    }

    #[test]
    fn unresolved_oversized_cluster_is_randomly_bounded() {
        let records = (0..5)
            .map(|index| dna(&format!("s{index}"), ""))
            .collect::<Vec<_>>();
        let kmers = all_kmer_counts(&records);
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let mut rng = SplitMix64::new(9);
        let mut audit = PosteriorAudit::new();
        let mut cache = PairPosteriorCache::new(&hmm, 0.01);
        let (clusters, _) = bound_coarse_clusters(
            &[(0..5).collect()],
            &records,
            &kmers,
            &hmm,
            0.01,
            1.0,
            2,
            16,
            &mut rng,
            &mut audit,
            &mut cache,
        )
        .unwrap();
        assert_eq!(
            clusters
                .iter()
                .map(|cluster| cluster.members.len())
                .sum::<usize>(),
            5
        );
        assert!(clusters.iter().all(|cluster| cluster.members.len() <= 2));
        assert!(clusters.iter().all(|cluster| cluster.randomly_subdivided));
    }

    #[test]
    fn consensus_deletes_gap_majority_and_breaks_symbol_ties_canonically() {
        let a = DnaSequence::from_ascii("A").unwrap().get(0).unwrap();
        let c = DnaSequence::from_ascii("C").unwrap().get(0).unwrap();
        let g = DnaSequence::from_ascii("G").unwrap().get(0).unwrap();
        let profile = Profile {
            rows: vec![
                WorkingRow {
                    input_index: 0,
                    cells: vec![
                        Some(PositionedSymbol {
                            symbol: a,
                            position: 0,
                        }),
                        Some(PositionedSymbol {
                            symbol: c,
                            position: 1,
                        }),
                        None,
                    ],
                },
                WorkingRow {
                    input_index: 1,
                    cells: vec![
                        Some(PositionedSymbol {
                            symbol: a,
                            position: 0,
                        }),
                        None,
                        Some(PositionedSymbol {
                            symbol: g,
                            position: 1,
                        }),
                    ],
                },
                WorkingRow {
                    input_index: 2,
                    cells: vec![
                        None,
                        None,
                        Some(PositionedSymbol {
                            symbol: g,
                            position: 0,
                        }),
                    ],
                },
            ],
        };
        assert_eq!(consensus_sequence(&profile).unwrap().to_string(), "AG");
    }

    #[test]
    fn ppp_ps_selects_exactly_the_limit_without_replacement() {
        let symbol = DnaSequence::from_ascii("A").unwrap().get(0).unwrap();
        let make_profile = |start: usize| WorkingProfile {
            rows: (start..start + 45)
                .map(|input_index| WorkingRow {
                    input_index,
                    cells: vec![Some(PositionedSymbol {
                        symbol,
                        position: 0,
                    })],
                })
                .collect(),
        };
        let left = make_profile(0);
        let right = make_profile(45);
        let mut first_rng = SplitMix64::new(7);
        let mut second_rng = SplitMix64::new(7);
        let first = selected_cross_pairs(&left, &right, 2_000, &mut first_rng).unwrap();
        let second = selected_cross_pairs(&left, &right, 2_000, &mut second_rng).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.0, Muscle5ProfileAlignmentMethod::PppPs);
        assert_eq!(first.1.len(), 2_000);
        assert!(first.1.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn redundant_pair_path_inserts_internal_and_terminal_columns() {
        let records = [dna("representative", "AC"), dna("member", "TAGC")];
        let mut profile = WorkingProfile::single(&records[0], 0);
        let deferred = DeferredSequence {
            sequence_index: 1,
            representative_index: 0,
            expected_error: 0.0,
            expected_accuracy: 2.0,
            operations: vec![
                PairOperation::ConsumeSecond,
                PairOperation::Match,
                PairOperation::ConsumeSecond,
                PairOperation::Match,
            ]
            .into_boxed_slice(),
        };
        reinsert_deferred(&mut profile, &deferred, &records).unwrap();
        let alignment = profile.to_alignment(&records).unwrap();
        assert_eq!(alignment.row(0).unwrap().to_ascii_with_gaps(), "-A-C");
        assert_eq!(alignment.row(1).unwrap().to_ascii_with_gaps(), "TAGC");
    }

    /// Mutated copies of one random sequence, so that clustering compares pairs
    /// that later stages request again.
    fn related_dna_family(count: usize, length: usize, seed: u64) -> Vec<SequenceRecord<Dna>> {
        let mut rng = SplitMix64::new(seed);
        let base = (0..length)
            .map(|_| b"ACGT"[rng.index(4)])
            .collect::<Vec<_>>();
        (0..count)
            .map(|member| {
                let mut sequence = base.clone();
                for _ in 0..length / 8 {
                    let at = rng.index(sequence.len());
                    sequence[at] = b"ACGT"[rng.index(4)];
                }
                sequence.truncate(length - rng.index(4));
                dna(
                    &format!("s{member}"),
                    std::str::from_utf8(&sequence).unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn reused_posteriors_leave_every_super5_result_unchanged() {
        let aligner = fast_super5();
        let records = related_dna_family(10, 60, 4);
        let parameters = aligner.parameters();
        let hmm = HmmParameters::published::<Dna>(
            parameters.perturbation_seed(),
            parameters.perturbation_amplitude(),
        )
        .unwrap()
        .prepare::<Dna>();
        let cutoff = parameters.posterior_cutoff();
        let mut reusing = PairPosteriorCache::new(&hmm, cutoff);
        let reused = run_super5_using(&aligner, &records, &hmm, &mut reusing).unwrap();
        let mut recomputing = PairPosteriorCache::without_reuse(&hmm, cutoff);
        let recomputed = run_super5_using(&aligner, &records, &hmm, &mut recomputing).unwrap();
        assert!(reusing.hits() > 0);
        assert_eq!(recomputing.hits(), 0);
        let reused = format!("{reused:?}");
        assert_eq!(reused, format!("{recomputed:?}"));
        assert_eq!(reused, format!("{:?}", aligner.align(&records).unwrap()));
    }

    #[test]
    fn super5_pipeline_is_deterministic_and_preserves_every_input() {
        let records = [
            dna("a", "AAAAAA"),
            dna("b", "CCCCCC"),
            dna("c", "GGGGGG"),
            dna("d", "TTTTTT"),
            dna("e", "ACGTAC"),
        ];
        let aligner = fast_super5();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(first.super5_evidence(), second.super5_evidence());
        let evidence = first.super5_evidence().unwrap();
        assert!(
            evidence
                .clusters()
                .iter()
                .all(|cluster| cluster.representative_members().len() <= 2)
        );
        assert_eq!(first.guide_tree().leaf_count(), evidence.clusters().len());
        assert!(first.evaluated_objective_pair_count() <= first.possible_objective_pair_count());
        assert!(first.memory_estimate().super5_auxiliary_bytes() > 0);
        for (row, record) in first.alignment().rows().iter().zip(&records) {
            assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
        }
    }

    #[test]
    fn super5_common_api_supports_protein_rna_and_redundant_empty_dna() {
        let protein_records = [
            SequenceRecord::new("p1", ProteinSequence::from_ascii("ARND").unwrap()).unwrap(),
            SequenceRecord::new("p2", ProteinSequence::from_ascii("CEQG").unwrap()).unwrap(),
        ];
        let protein_parameters = Muscle5Parameters::published_super5()
            .with_iterations(0, 0)
            .with_super5_limits(0.0, 0.30, 0.10, 2, 16, 2)
            .unwrap();
        let protein = Muscle5::<Protein>::with_parameters(protein_parameters)
            .unwrap()
            .align(&protein_records)
            .unwrap();
        assert_eq!(protein.alignment().sequence_count(), 2);

        let rna_records = [
            SequenceRecord::new("r1", RnaSequence::from_ascii("ACGU").unwrap()).unwrap(),
            SequenceRecord::new("r2", RnaSequence::from_ascii("UGCA").unwrap()).unwrap(),
        ];
        let rna_parameters = Muscle5Parameters::published_super5()
            .with_iterations(0, 0)
            .with_super5_limits(0.0, 0.30, 0.10, 2, 16, 2)
            .unwrap();
        let rna = Muscle5::<Rna>::with_parameters(rna_parameters)
            .unwrap()
            .align(&rna_records)
            .unwrap();
        assert_eq!(rna.alignment().row(0).unwrap().to_ungapped_ascii(), "ACGU");

        let empty_records = [dna("e1", ""), dna("e2", "")];
        let empty_parameters = Muscle5Parameters::published_super5().with_iterations(0, 0);
        let empty = Muscle5::<Dna>::with_parameters(empty_parameters)
            .unwrap()
            .align(&empty_records)
            .unwrap();
        assert_eq!(empty.alignment().width(), 0);
        assert_eq!(
            empty
                .super5_evidence()
                .unwrap()
                .redundancy_assignments()
                .len(),
            1
        );
    }
}
