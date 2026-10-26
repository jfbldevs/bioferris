use core::marker::PhantomData;
use core::mem::size_of;
use std::collections::{BTreeMap, HashSet};

use align_core::Symbol;
use align_msa::{DistanceMatrix, GuideTree, MsaError, MultipleAlignment, SequenceRecord};

use crate::muscle5_hmm::HmmParameters;
use crate::muscle5_posterior::{PairPosteriorCache, PosteriorLibrary, expected_error_matrix};
use crate::muscle5_profile::{
    ProgressiveMergeEvidence, progressive_alignment, quantized_objective, realign_partition,
    total_objective,
};
use crate::muscle5_rng::{REFINE_DOMAIN, SplitMix64, derived_seed};
use crate::muscle5_super5::{Muscle5Super5Evidence, run_super5};
use crate::muscle5_tree::expected_error_upgma;
use crate::{GuideTreePermutation, Muscle5Alphabet, Muscle5Error, Muscle5Mode, Muscle5Parameters};

/// One postorder probabilistic-profile merge in MUSCLE 5.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5ProgressiveStep {
    pub(crate) left_members: Box<[usize]>,
    pub(crate) right_members: Box<[usize]>,
    pub(crate) result_width: usize,
    pub(crate) expected_accuracy: f64,
    pub(crate) operation_count: usize,
}

impl Muscle5ProgressiveStep {
    /// Returns stable input indices in the left child profile.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns stable input indices in the right child profile.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }

    /// Returns the MSA width produced by this merge.
    #[must_use]
    pub const fn result_width(&self) -> usize {
        self.result_width
    }

    /// Returns the summed cross-profile posterior score on the selected path.
    #[must_use]
    pub const fn expected_accuracy(&self) -> f64 {
        self.expected_accuracy
    }

    /// Returns the number of profile-column traceback operations.
    #[must_use]
    pub const fn operation_count(&self) -> usize {
        self.operation_count
    }
}

impl From<ProgressiveMergeEvidence> for Muscle5ProgressiveStep {
    fn from(evidence: ProgressiveMergeEvidence) -> Self {
        Self {
            left_members: evidence.left_members,
            right_members: evidence.right_members,
            result_width: evidence.result_width,
            expected_accuracy: evidence.expected_accuracy,
            operation_count: evidence.operation_count,
        }
    }
}

/// One deterministic random-partition refinement decision.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Muscle5RefinementAttempt {
    pub(crate) iteration: usize,
    pub(crate) first_partition: Box<[usize]>,
    pub(crate) fixed_score_before: i128,
    pub(crate) fixed_score_after: i128,
    pub(crate) accepted: bool,
    pub(crate) candidate_width: usize,
    pub(crate) profile_operations: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RejectedRefinement {
    fixed_score_before: i128,
    fixed_score_after: i128,
    candidate_width: usize,
    profile_operations: usize,
}

impl Muscle5RefinementAttempt {
    /// Returns the zero-based attempt index.
    #[must_use]
    pub const fn iteration(&self) -> usize {
        self.iteration
    }

    /// Returns stable row indices assigned to the first partition.
    #[must_use]
    pub const fn first_partition(&self) -> &[usize] {
        &self.first_partition
    }

    /// Returns the independently replayed cross-partition score before realignment.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn score_before(&self) -> f64 {
        self.fixed_score_before as f64 / 1_000_000_000.0
    }

    /// Returns the independently replayed cross-partition candidate score.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn score_after(&self) -> f64 {
        self.fixed_score_after as f64 / 1_000_000_000.0
    }

    /// Returns whether strict fixed-precision improvement retained the candidate.
    #[must_use]
    pub const fn accepted(&self) -> bool {
        self.accepted
    }

    /// Returns the candidate MSA width.
    #[must_use]
    pub const fn candidate_width(&self) -> usize {
        self.candidate_width
    }

    /// Returns the number of operations in the candidate profile path.
    #[must_use]
    pub const fn profile_operation_count(&self) -> usize {
        self.profile_operations
    }
}

/// Algorithm-owned storage and high-water estimates for one MUSCLE 5 call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Muscle5MemoryEstimate {
    pub(crate) hmm_bytes: usize,
    pub(crate) posterior_library_bytes: usize,
    pub(crate) pair_hmm_peak_bytes: usize,
    pub(crate) consistency_peak_bytes: usize,
    pub(crate) distance_and_tree_bytes: usize,
    pub(crate) profile_dp_peak_bytes: usize,
    pub(crate) super5_auxiliary_bytes: usize,
    pub(crate) returned_alignment_bytes: usize,
    pub(crate) estimated_peak_bytes: usize,
}

impl Muscle5MemoryEstimate {
    /// Returns prepared five-state HMM storage.
    #[must_use]
    pub const fn hmm_bytes(self) -> usize {
        self.hmm_bytes
    }

    /// Returns final sparse posterior-library storage.
    #[must_use]
    pub const fn posterior_library_bytes(self) -> usize {
        self.posterior_library_bytes
    }

    /// Returns the largest forward/backward pair-HMM scratch allocation.
    #[must_use]
    pub const fn pair_hmm_peak_bytes(self) -> usize {
        self.pair_hmm_peak_bytes
    }

    /// Returns the largest consistency-accumulator allocation.
    #[must_use]
    pub const fn consistency_peak_bytes(self) -> usize {
        self.consistency_peak_bytes
    }

    /// Returns expected-error matrix plus UPGMA workspace storage.
    #[must_use]
    pub const fn distance_and_tree_bytes(self) -> usize {
        self.distance_and_tree_bytes
    }

    /// Returns the largest progressive or refinement profile-DP workspace.
    #[must_use]
    pub const fn profile_dp_peak_bytes(self) -> usize {
        self.profile_dp_peak_bytes
    }

    /// Returns Super5 k-mer, clustering, sampling, and retained-profile storage.
    #[must_use]
    pub const fn super5_auxiliary_bytes(self) -> usize {
        self.super5_auxiliary_bytes
    }

    /// Returns aligned cells and identifiers owned by the returned MSA.
    #[must_use]
    pub const fn returned_alignment_bytes(self) -> usize {
        self.returned_alignment_bytes
    }

    /// Returns a conservative algorithm-owned high-water estimate.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// Complete deterministic PPP result and scientific audit evidence.
#[derive(Clone, Debug)]
pub struct Muscle5Alignment<A: Muscle5Alphabet> {
    pub(crate) alignment: MultipleAlignment<A>,
    pub(crate) mode: Muscle5Mode,
    pub(crate) perturbation_seed: u64,
    pub(crate) guide_tree_permutation: GuideTreePermutation,
    pub(crate) expected_error_distances: DistanceMatrix,
    pub(crate) guide_tree: GuideTree,
    pub(crate) raw_posterior_checksum: u64,
    pub(crate) consistent_posterior_checksum: u64,
    pub(crate) raw_posterior_entries: usize,
    pub(crate) consistent_posterior_entries: usize,
    pub(crate) progressive_steps: Box<[Muscle5ProgressiveStep]>,
    pub(crate) refinement_attempts: Box<[Muscle5RefinementAttempt]>,
    pub(crate) fixed_final_objective: i128,
    pub(crate) evaluated_objective_pairs: usize,
    pub(crate) possible_objective_pairs: usize,
    pub(crate) super5_evidence: Option<Muscle5Super5Evidence>,
    pub(crate) memory: Muscle5MemoryEstimate,
}

impl<A: Muscle5Alphabet> Muscle5Alignment<A> {
    /// Returns the validated MSA in stable input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns whether `Align` or `Super5` produced this result.
    #[must_use]
    pub const fn mode(&self) -> Muscle5Mode {
        self.mode
    }

    /// Returns the HMM perturbation seed.
    #[must_use]
    pub const fn perturbation_seed(&self) -> u64 {
        self.perturbation_seed
    }

    /// Returns the guide-tree resolution used for progressive alignment.
    #[must_use]
    pub const fn guide_tree_permutation(&self) -> GuideTreePermutation {
        self.guide_tree_permutation
    }

    /// Returns pairwise posterior expected-error distances.
    #[must_use]
    pub const fn expected_error_distances(&self) -> &DistanceMatrix {
        &self.expected_error_distances
    }

    /// Returns the deterministic, possibly permuted guide tree.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        &self.guide_tree
    }

    /// Returns the checksum before probabilistic consistency.
    #[must_use]
    pub const fn raw_posterior_checksum(&self) -> u64 {
        self.raw_posterior_checksum
    }

    /// Returns the checksum after all configured consistency rounds.
    #[must_use]
    pub const fn consistent_posterior_checksum(&self) -> u64 {
        self.consistent_posterior_checksum
    }

    /// Returns retained posterior entries before consistency, including both orientations.
    #[must_use]
    pub const fn raw_posterior_entry_count(&self) -> usize {
        self.raw_posterior_entries
    }

    /// Returns retained posterior entries after consistency, including both orientations.
    #[must_use]
    pub const fn consistent_posterior_entry_count(&self) -> usize {
        self.consistent_posterior_entries
    }

    /// Returns postorder progressive-profile merge evidence.
    #[must_use]
    pub const fn progressive_steps(&self) -> &[Muscle5ProgressiveStep] {
        &self.progressive_steps
    }

    /// Returns every configured random-partition refinement decision.
    #[must_use]
    pub const fn refinement_attempts(&self) -> &[Muscle5RefinementAttempt] {
        &self.refinement_attempts
    }

    /// Returns the final all-pairs posterior objective at `1e-9` resolution.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn final_objective_score(&self) -> f64 {
        self.fixed_final_objective as f64 / 1_000_000_000.0
    }

    /// Returns the number of sequence pairs represented in the final objective.
    #[must_use]
    pub const fn evaluated_objective_pair_count(&self) -> usize {
        self.evaluated_objective_pairs
    }

    /// Returns the number of possible unordered input pairs.
    #[must_use]
    pub const fn possible_objective_pair_count(&self) -> usize {
        self.possible_objective_pairs
    }

    /// Returns characteristic Super5 evidence, or `None` for PPP `Align`.
    #[must_use]
    pub const fn super5_evidence(&self) -> Option<&Muscle5Super5Evidence> {
        self.super5_evidence.as_ref()
    }

    /// Returns explicit algorithm-owned memory estimates.
    #[must_use]
    pub const fn memory_estimate(&self) -> Muscle5MemoryEstimate {
        self.memory
    }
}

/// One named member of a MUSCLE 5 stratified or diversified ensemble.
#[derive(Clone, Debug)]
pub struct Muscle5EnsembleMember<A: Muscle5Alphabet> {
    name: String,
    alignment: Muscle5Alignment<A>,
}

impl<A: Muscle5Alphabet> Muscle5EnsembleMember<A> {
    /// Returns the stable `<permutation>.<seed>` replicate name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the complete independently owned alignment result.
    #[must_use]
    pub const fn result(&self) -> &Muscle5Alignment<A> {
        &self.alignment
    }
}

/// Ordered collection of independently owned MUSCLE 5 replicate results.
#[derive(Clone, Debug)]
pub struct Muscle5Ensemble<A: Muscle5Alphabet> {
    members: Box<[Muscle5EnsembleMember<A>]>,
}

impl<A: Muscle5Alphabet> Muscle5Ensemble<A> {
    /// Returns ensemble members in deterministic generation order.
    #[must_use]
    pub const fn members(&self) -> &[Muscle5EnsembleMember<A>] {
        &self.members
    }

    /// Returns the number of generated replicates.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.members.len()
    }

    /// Returns whether the ensemble contains no replicates.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.members.is_empty()
    }
}

/// Scalar reconstruction of `MUSCLE-5.1.0-Nature-2022-default-v2`.
///
/// Version 2 changes only how refinement scores the current alignment before
/// comparing it with a candidate: both scores now come from the candidate's own
/// column-score matrix, so rounding can no longer report a decrease.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5<A: Muscle5Alphabet> {
    parameters: Muscle5Parameters,
    marker: PhantomData<A>,
}

impl<A: Muscle5Alphabet> Muscle5<A> {
    /// Creates the publication-default PPP `Align` profile.
    #[must_use]
    pub const fn published_align() -> Self {
        Self {
            parameters: Muscle5Parameters::published_align(),
            marker: PhantomData,
        }
    }

    /// Creates the publication-default `Super5` profile.
    #[must_use]
    pub const fn published_super5() -> Self {
        Self {
            parameters: Muscle5Parameters::published_super5(),
            marker: PhantomData,
        }
    }

    /// Creates an aligner from a validated versioned parameter profile.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] for an invalid probability, threshold, or hard limit.
    pub fn with_parameters(parameters: Muscle5Parameters) -> Result<Self, Muscle5Error> {
        parameters.validate()?;
        Ok(Self {
            parameters,
            marker: PhantomData,
        })
    }

    /// Returns the immutable versioned parameter profile.
    #[must_use]
    pub const fn parameters(&self) -> &Muscle5Parameters {
        &self.parameters
    }

    /// Aligns one or more typed, named, ungapped sequences.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] for invalid input, arithmetic failure, or a
    /// structural/replay contradiction.
    pub fn align(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Muscle5Alignment<A>, Muscle5Error> {
        match self.parameters.mode() {
            Muscle5Mode::Align => self.align_ppp(records, None),
            Muscle5Mode::Super5 => run_super5(self, records),
        }
    }

    /// Builds the published stratified schedule of four permutations per seed.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] if `seed_count` is zero or a replicate fails.
    pub fn stratified_ensemble(
        &self,
        records: &[SequenceRecord<A>],
        seed_count: usize,
    ) -> Result<Muscle5Ensemble<A>, Muscle5Error> {
        if seed_count == 0 {
            return Err(Muscle5Error::InvalidParameter {
                parameter: "seed_count",
                requirement: "be at least one",
            });
        }
        let capacity = seed_count
            .checked_mul(GuideTreePermutation::ALL.len())
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "stratified ensemble size",
            })?;
        let mut members = Vec::with_capacity(capacity);
        for seed in 0..seed_count {
            let seed = u64::try_from(seed).map_err(|_| Muscle5Error::ArithmeticOverflow {
                stage: "stratified ensemble seed",
            })?;
            for permutation in GuideTreePermutation::ALL {
                members.push(self.ensemble_member(records, seed, permutation)?);
            }
        }
        Ok(Muscle5Ensemble {
            members: members.into_boxed_slice(),
        })
    }

    /// Builds a diversified schedule cycling the four guide-tree permutations.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] if `replicate_count` is zero or a replicate fails.
    pub fn diversified_ensemble(
        &self,
        records: &[SequenceRecord<A>],
        replicate_count: usize,
    ) -> Result<Muscle5Ensemble<A>, Muscle5Error> {
        if replicate_count == 0 {
            return Err(Muscle5Error::InvalidParameter {
                parameter: "replicate_count",
                requirement: "be at least one",
            });
        }
        let mut members = Vec::with_capacity(replicate_count);
        for replicate in 0..replicate_count {
            let seed = u64::try_from(replicate).map_err(|_| Muscle5Error::ArithmeticOverflow {
                stage: "diversified ensemble seed",
            })?;
            let permutation = GuideTreePermutation::ALL[replicate % 4];
            members.push(self.ensemble_member(records, seed, permutation)?);
        }
        Ok(Muscle5Ensemble {
            members: members.into_boxed_slice(),
        })
    }

    fn ensemble_member(
        &self,
        records: &[SequenceRecord<A>],
        seed: u64,
        permutation: GuideTreePermutation,
    ) -> Result<Muscle5EnsembleMember<A>, Muscle5Error> {
        let parameters = self
            .parameters
            .clone()
            .with_perturbation_seed(seed)
            .with_guide_tree_permutation(permutation);
        let aligner = Self::with_parameters(parameters)?;
        let alignment = aligner.align(records)?;
        Ok(Muscle5EnsembleMember {
            name: format!("{}.{seed}", permutation.label()),
            alignment,
        })
    }

    /// Aligns in [`Muscle5Mode::Align`], taking pair posteriors already in
    /// `cache` from it; Super5 aligns its clusters this way.
    pub(crate) fn align_reusing(
        &self,
        records: &[SequenceRecord<A>],
        cache: &PairPosteriorCache<A>,
    ) -> Result<Muscle5Alignment<A>, Muscle5Error> {
        self.align_ppp(records, Some(cache))
    }

    #[allow(clippy::too_many_lines)]
    fn align_ppp(
        &self,
        records: &[SequenceRecord<A>],
        cache: Option<&PairPosteriorCache<A>>,
    ) -> Result<Muscle5Alignment<A>, Muscle5Error> {
        validate_records(records)?;
        self.parameters.validate()?;
        let hmm_parameters = HmmParameters::published::<A>(
            self.parameters.perturbation_seed(),
            self.parameters.perturbation_amplitude(),
        )?;
        let hmm = hmm_parameters.prepare::<A>();
        let mut library =
            PosteriorLibrary::build(records, &hmm, self.parameters.posterior_cutoff(), cache)?;
        let raw_posterior_checksum = library.checksum();
        let raw_posterior_entries = library.entry_count();
        library.transform_consistency(
            self.parameters.consistency_iterations(),
            self.parameters.posterior_cutoff(),
        )?;
        let consistent_posterior_checksum = library.checksum();
        let consistent_posterior_entries = library.entry_count();
        let expected_error_distances = expected_error_matrix(&library)?;
        let original_tree = expected_error_upgma(&expected_error_distances)?;
        let tree = original_tree.permuted(self.parameters.guide_tree_permutation())?;
        let guide_tree = tree.to_public()?;
        let progressive = progressive_alignment(tree.root(), records, &library)?;
        let mut profile = progressive.profile;
        let mut profile_dp_peak_bytes = progressive.peak_workspace_bytes;
        let mut refinement_attempts = Vec::with_capacity(
            self.parameters
                .refinement_iterations()
                .min(usize::from(records.len() > 1) * self.parameters.refinement_iterations()),
        );
        if records.len() > 1 {
            let mut rng = SplitMix64::new(derived_seed(
                self.parameters.perturbation_seed(),
                REFINE_DOMAIN,
            ));
            // A rejected partition is a pure function of the current profile and
            // posterior library. Replaying it cannot become acceptable until an
            // accepted candidate changes the profile, at which point the cache is
            // invalidated. The published 100-attempt schedule and evidence remain
            // intact while duplicate dynamic programs are avoided exactly.
            let mut rejected_partitions = BTreeMap::<Vec<bool>, RejectedRefinement>::new();
            for iteration in 0..self.parameters.refinement_iterations() {
                let mut first_partition = (0..records.len())
                    .map(|_| rng.next_u64() & 1 == 1)
                    .collect::<Vec<_>>();
                if first_partition.iter().all(|selected| !selected) {
                    first_partition[0] = true;
                } else if first_partition.iter().all(|selected| *selected) {
                    first_partition[0] = false;
                }
                if let Some(cached) = rejected_partitions.get(&first_partition).copied() {
                    let members = selected_members(&first_partition);
                    refinement_attempts.push(Muscle5RefinementAttempt {
                        iteration,
                        first_partition: members,
                        fixed_score_before: cached.fixed_score_before,
                        fixed_score_after: cached.fixed_score_after,
                        accepted: false,
                        candidate_width: cached.candidate_width,
                        profile_operations: cached.profile_operations,
                    });
                    continue;
                }
                let realignment = realign_partition(&profile, &first_partition, &library)?;
                let candidate = realignment.candidate;
                profile_dp_peak_bytes = profile_dp_peak_bytes.max(candidate.workspace_bytes);
                // Both scores come from one column-score matrix, each added along
                // its own path in path order: the optimum is never below the
                // current arrangement, and an identical candidate scores exactly
                // the same. Version 1 summed the current objective row pair by row
                // pair, which could round it across a fixed-precision step.
                let score_before = realignment.current_score;
                let score_after = candidate.expected_accuracy;
                // Checked builds keep the independent row-pair objectives as a
                // replay of both scores; optimized builds skip the extra scans.
                #[cfg(debug_assertions)]
                {
                    use crate::muscle5_profile::cross_partition_objective;
                    let replayed_before =
                        cross_partition_objective(&profile, &first_partition, &library)?;
                    if (replayed_before - score_before).abs() > 1.0e-9 * score_before.abs().max(1.0)
                    {
                        return Err(Muscle5Error::InternalInvariant {
                            message: "PPP refinement current objective disagrees with its path replay",
                        });
                    }
                    let replayed_after =
                        cross_partition_objective(&candidate.profile, &first_partition, &library)?;
                    if (replayed_after - score_after).abs() > 1.0e-9 * score_after.abs().max(1.0) {
                        return Err(Muscle5Error::InternalInvariant {
                            message: "PPP refinement traceback objective disagrees with its DP optimum",
                        });
                    }
                }
                let fixed_score_before = quantized_objective(score_before)?;
                let fixed_score_after = quantized_objective(score_after)?;
                if fixed_score_after < fixed_score_before {
                    return Err(Muscle5Error::InternalInvariant {
                        message: "PPP refinement candidate decreased its partition objective",
                    });
                }
                let accepted = fixed_score_after > fixed_score_before;
                let candidate_width = candidate.profile.width();
                let profile_operations = candidate.operations.len();
                let members = selected_members(&first_partition);
                refinement_attempts.push(Muscle5RefinementAttempt {
                    iteration,
                    first_partition: members,
                    fixed_score_before,
                    fixed_score_after,
                    accepted,
                    candidate_width,
                    profile_operations,
                });
                if accepted {
                    profile = candidate.profile;
                    rejected_partitions.clear();
                } else {
                    rejected_partitions.insert(
                        first_partition,
                        RejectedRefinement {
                            fixed_score_before,
                            fixed_score_after,
                            candidate_width,
                            profile_operations,
                        },
                    );
                }
            }
        }
        let final_objective = total_objective(&profile, &library)?;
        let fixed_final_objective = quantized_objective(final_objective)?;
        let alignment = profile.to_alignment(records)?;
        verify_ungapped_replay(&alignment, records)?;
        let returned_alignment_bytes = estimate_returned_alignment(&alignment)?;
        let distance_bytes = expected_error_distances
            .as_flat_slice()
            .len()
            .checked_mul(size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(original_tree.workspace_bytes()))
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "MUSCLE 5 distance and tree memory estimate",
            })?;
        let memory = estimate_memory(
            hmm.storage_bytes(),
            &library,
            distance_bytes,
            profile_dp_peak_bytes,
            returned_alignment_bytes,
        )?;
        Ok(Muscle5Alignment {
            alignment,
            mode: Muscle5Mode::Align,
            perturbation_seed: self.parameters.perturbation_seed(),
            guide_tree_permutation: self.parameters.guide_tree_permutation(),
            expected_error_distances,
            guide_tree,
            raw_posterior_checksum,
            consistent_posterior_checksum,
            raw_posterior_entries,
            consistent_posterior_entries,
            progressive_steps: progressive
                .merges
                .into_iter()
                .map(Muscle5ProgressiveStep::from)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            refinement_attempts: refinement_attempts.into_boxed_slice(),
            fixed_final_objective,
            evaluated_objective_pairs: unordered_pair_count(records.len())?,
            possible_objective_pairs: unordered_pair_count(records.len())?,
            super5_evidence: None,
            memory,
        })
    }
}

fn selected_members(partition: &[bool]) -> Box<[usize]> {
    partition
        .iter()
        .enumerate()
        .filter_map(|(index, selected)| selected.then_some(index))
        .collect::<Vec<_>>()
        .into_boxed_slice()
}

pub(crate) fn validate_records<A: Muscle5Alphabet>(
    records: &[SequenceRecord<A>],
) -> Result<(), Muscle5Error> {
    if records.is_empty() {
        return Err(MsaError::EmptyCollection.into());
    }
    let mut identifiers = HashSet::with_capacity(records.len());
    for record in records {
        if !identifiers.insert(record.identifier()) {
            return Err(MsaError::DuplicateIdentifier {
                identifier: record.identifier().to_owned(),
            }
            .into());
        }
    }
    Ok(())
}

pub(crate) fn verify_ungapped_replay<A: Muscle5Alphabet>(
    alignment: &MultipleAlignment<A>,
    records: &[SequenceRecord<A>],
) -> Result<(), Muscle5Error> {
    if alignment.sequence_count() != records.len() {
        return Err(Muscle5Error::InternalInvariant {
            message: "final MUSCLE 5 alignment has the wrong row count",
        });
    }
    for (row, record) in alignment.rows().iter().zip(records) {
        let replay = row
            .symbols()
            .iter()
            .filter_map(|cell| *cell)
            .collect::<Vec<_>>();
        if row.identifier() != record.identifier() || replay != record.sequence().as_slice() {
            return Err(Muscle5Error::InternalInvariant {
                message: "final MUSCLE 5 alignment changed an input biological sequence",
            });
        }
    }
    Ok(())
}

pub(crate) fn estimate_returned_alignment<A: Muscle5Alphabet>(
    alignment: &MultipleAlignment<A>,
) -> Result<usize, Muscle5Error> {
    let cells = alignment
        .sequence_count()
        .checked_mul(alignment.width())
        .and_then(|value| value.checked_mul(size_of::<Option<Symbol<A>>>()))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "returned MUSCLE 5 alignment memory estimate",
        })?;
    let identifiers = alignment.rows().iter().try_fold(0_usize, |total, row| {
        total.checked_add(row.identifier().len())
    });
    cells
        .checked_add(identifiers.ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "returned MUSCLE 5 alignment memory estimate",
        })?)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "returned MUSCLE 5 alignment memory estimate",
        })
}

fn estimate_memory(
    hmm_bytes: usize,
    library: &PosteriorLibrary,
    distance_and_tree_bytes: usize,
    profile_dp_peak_bytes: usize,
    returned_alignment_bytes: usize,
) -> Result<Muscle5MemoryEstimate, Muscle5Error> {
    let posterior_library_bytes = library.storage_bytes();
    let pair_hmm_peak_bytes = library.pair_hmm_peak_bytes();
    let consistency_peak_bytes = library.consistency_peak_bytes();
    let phase_peak = pair_hmm_peak_bytes
        .max(consistency_peak_bytes)
        .max(profile_dp_peak_bytes);
    let estimated_peak_bytes = hmm_bytes
        .checked_add(posterior_library_bytes)
        .and_then(|bytes| bytes.checked_add(distance_and_tree_bytes))
        .and_then(|bytes| bytes.checked_add(returned_alignment_bytes))
        .and_then(|bytes| bytes.checked_add(phase_peak))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "MUSCLE 5 peak memory estimate",
        })?;
    Ok(Muscle5MemoryEstimate {
        hmm_bytes,
        posterior_library_bytes,
        pair_hmm_peak_bytes,
        consistency_peak_bytes,
        distance_and_tree_bytes,
        profile_dp_peak_bytes,
        super5_auxiliary_bytes: 0,
        returned_alignment_bytes,
        estimated_peak_bytes,
    })
}

pub(crate) fn unordered_pair_count(count: usize) -> Result<usize, Muscle5Error> {
    count
        .checked_mul(count.saturating_sub(1))
        .and_then(|value| value.checked_div(2))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "unordered sequence-pair count",
        })
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::SequenceRecord;

    use super::Muscle5;
    use crate::{GuideTreePermutation, Muscle5Parameters};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn align_pipeline_is_deterministic_and_preserves_dna() {
        let records = [
            dna("a", "ACGT"),
            dna("b", "AGT"),
            dna("c", "ACCT"),
            dna("d", "TACGT"),
        ];
        let parameters = Muscle5Parameters::published_align().with_iterations(1, 5);
        let aligner = Muscle5::<Dna>::with_parameters(parameters).unwrap();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(
            first.consistent_posterior_checksum(),
            second.consistent_posterior_checksum()
        );
        assert_eq!(first.progressive_steps().len(), records.len() - 1);
        assert_eq!(first.refinement_attempts().len(), 5);
        assert!(
            first.refinement_attempts().iter().all(
                |attempt| !attempt.accepted() || attempt.score_after() > attempt.score_before()
            )
        );
        for (row, record) in first.alignment().rows().iter().zip(&records) {
            assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
        }
        assert!(first.memory_estimate().estimated_peak_bytes() > 0);
    }

    #[test]
    fn protein_rna_and_empty_singleton_use_the_common_api() {
        let proteins = [
            SequenceRecord::new("p1", ProteinSequence::from_ascii("ARND").unwrap()).unwrap(),
            SequenceRecord::new("p2", ProteinSequence::from_ascii("ARD").unwrap()).unwrap(),
        ];
        let protein = Muscle5::<Protein>::with_parameters(
            Muscle5Parameters::published_align().with_iterations(0, 0),
        )
        .unwrap()
        .align(&proteins)
        .unwrap();
        assert_eq!(protein.alignment().sequence_count(), 2);

        let rna = [
            SequenceRecord::new("r1", RnaSequence::from_ascii("ACGU").unwrap()).unwrap(),
            SequenceRecord::new("r2", RnaSequence::from_ascii("AGU").unwrap()).unwrap(),
        ];
        let rna = Muscle5::<Rna>::with_parameters(
            Muscle5Parameters::published_align().with_iterations(0, 0),
        )
        .unwrap()
        .align(&rna)
        .unwrap();
        assert_eq!(rna.alignment().row(0).unwrap().to_ungapped_ascii(), "ACGU");

        let empty = [dna("empty", "")];
        let singleton = Muscle5::<Dna>::published_align().align(&empty).unwrap();
        assert_eq!(singleton.alignment().width(), 0);
        assert!(singleton.progressive_steps().is_empty());
        assert!(singleton.refinement_attempts().is_empty());
    }

    #[test]
    fn duplicate_identifiers_and_empty_collection_are_rejected() {
        assert!(Muscle5::<Dna>::published_align().align(&[]).is_err());
        let duplicate = [dna("same", "A"), dna("same", "C")];
        assert!(Muscle5::<Dna>::published_align().align(&duplicate).is_err());
    }

    #[test]
    fn ensemble_schedules_follow_the_published_order() {
        let records = [dna("a", "AC"), dna("b", "AG"), dna("c", "AT")];
        let aligner = Muscle5::<Dna>::with_parameters(
            Muscle5Parameters::published_align().with_iterations(0, 0),
        )
        .unwrap();
        let stratified = aligner.stratified_ensemble(&records, 1).unwrap();
        assert_eq!(stratified.len(), 4);
        assert_eq!(
            stratified
                .members()
                .iter()
                .map(super::Muscle5EnsembleMember::name)
                .collect::<Vec<_>>(),
            ["none.0", "abc.0", "acb.0", "bca.0"]
        );
        let diversified = aligner.diversified_ensemble(&records, 5).unwrap();
        assert_eq!(
            diversified
                .members()
                .iter()
                .map(|member| member.result().guide_tree_permutation())
                .collect::<Vec<_>>(),
            [
                GuideTreePermutation::None,
                GuideTreePermutation::Abc,
                GuideTreePermutation::Acb,
                GuideTreePermutation::Bca,
                GuideTreePermutation::None,
            ]
        );
    }
}
