use core::marker::PhantomData;
use core::mem::size_of;
use std::time::{Duration, Instant};

use align_core::{PhaseTimings, PhaseTimingsBuilder, Profiled, Symbol};
use align_msa::{GuideTree, MultipleAlignment, SequenceRecord};

use crate::aligned_distance::FftNs2DistanceMatrix;
use crate::fft_ns1::{
    FftNs1Merge, MergePhaseTimings, ProgressiveMemory, merge_pair, merge_pair_profiled,
};
use crate::matrix::ScientificScorer;
use crate::profile::{WorkingProfile, WorkingRow};
use crate::refinement::{FftNsiPairWeights, FftNsiPartition, GotohWeighting, TreePartition};
use crate::upgma::build_upgma_from_values;
use crate::wsp::{WspPairScoreCache, weighted_sum_of_pairs};
use crate::{FftNsiParameters, MafftAlphabet, MafftError, MafftFftNs2};

/// Stable identifier for the native-compatible 2002 FFT-NS-i scientific profile.
pub const FFT_NSI_PROFILE_IDENTIFIER: &str = "MAFFT-2002-FFT-NS-i-v2";

/// Why an FFT-NS-i call stopped refining its initial FFT-NS-2 alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FftNsiTermination {
    /// A complete best-first cycle found no strictly better WSP alignment.
    Converged,
    /// An accepted candidate score repeated an earlier best-first score.
    Oscillating,
    /// The versioned maximum number of accepted cycles was reached.
    CycleLimit,
}

/// Aggregate opt-in timings for work performed while evaluating FFT-NS-i candidates.
///
/// Durations are summed over individual candidates. With multiple candidate
/// workers this is aggregate worker work rather than elapsed batch wall time;
/// [`PhaseTimings`] retains the complete-call and batch wall-clock view.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FftNsiCandidatePhaseTimings {
    candidate_count: usize,
    total: Duration,
    projection: Duration,
    profile_preparation: Duration,
    fft_anchor_search: Duration,
    profile_dp: Duration,
    profile_dp_block_count: usize,
    profile_dp_cell_count: usize,
    profile_dp_homology_term_evaluations: usize,
    profile_dp_homology_term_cell_counts: [usize; 6],
    profile_dp_workspace_initialization: Duration,
    profile_dp_recurrence_fill: Duration,
    profile_dp_endpoint_traceback: Duration,
    profile_dp_block_accounting: Duration,
    profile_dp_score_replay: Duration,
    profile_assembly: Duration,
    alignment_materialization: Duration,
    validation: Duration,
    wsp_scoring: Duration,
    accounting: Duration,
}

impl FftNsiCandidatePhaseTimings {
    /// Returns the number of tree-edge candidates included in the aggregates.
    #[must_use]
    pub const fn candidate_count(self) -> usize {
        self.candidate_count
    }

    /// Returns aggregate wall time across all individual candidate evaluations.
    #[must_use]
    pub const fn total(self) -> Duration {
        self.total
    }

    /// Returns the sum of the named, non-overlapping candidate subphases.
    #[must_use]
    pub fn attributed(self) -> Duration {
        self.phases()
            .into_iter()
            .fold(Duration::ZERO, |total, (_, elapsed)| {
                total.saturating_add(elapsed)
            })
    }

    /// Returns candidate work not attributed to a named subphase.
    #[must_use]
    pub fn unattributed(self) -> Duration {
        self.total.saturating_sub(self.attributed())
    }

    /// Returns the stable candidate subphase names and aggregate durations.
    #[must_use]
    pub const fn phases(self) -> [(&'static str, Duration); 9] {
        [
            ("projection", self.projection),
            ("profile_preparation", self.profile_preparation),
            ("fft_anchor_search", self.fft_anchor_search),
            ("profile_dp", self.profile_dp),
            ("profile_assembly", self.profile_assembly),
            ("alignment_materialization", self.alignment_materialization),
            ("validation", self.validation),
            ("wsp_scoring", self.wsp_scoring),
            ("accounting", self.accounting),
        ]
    }

    /// Returns the number of profile-DP blocks evaluated across candidates.
    #[must_use]
    pub const fn profile_dp_block_count(self) -> usize {
        self.profile_dp_block_count
    }

    /// Returns the number of non-empty profile-DP cells evaluated.
    #[must_use]
    pub const fn profile_dp_cell_count(self) -> usize {
        self.profile_dp_cell_count
    }

    /// Returns the number of sparse homology terms evaluated across profile-DP cells.
    #[must_use]
    pub const fn profile_dp_homology_term_evaluations(self) -> usize {
        self.profile_dp_homology_term_evaluations
    }

    /// Returns cell counts for 0, 1, 2, 3, 4, and at least 5 homology terms.
    #[must_use]
    pub const fn profile_dp_homology_term_cell_counts(self) -> [usize; 6] {
        self.profile_dp_homology_term_cell_counts
    }

    /// Returns profile-DP work including orchestration outside its subphases.
    #[must_use]
    pub const fn profile_dp_total(self) -> Duration {
        self.profile_dp
    }

    /// Returns stable nested profile-DP phase names and durations.
    #[must_use]
    pub const fn profile_dp_phases(self) -> [(&'static str, Duration); 5] {
        [
            (
                "workspace_initialization",
                self.profile_dp_workspace_initialization,
            ),
            ("recurrence_fill", self.profile_dp_recurrence_fill),
            ("endpoint_traceback", self.profile_dp_endpoint_traceback),
            ("block_accounting", self.profile_dp_block_accounting),
            ("score_replay", self.profile_dp_score_replay),
        ]
    }

    /// Returns profile-DP work attributed to nested phases.
    #[must_use]
    pub fn profile_dp_attributed(self) -> Duration {
        self.profile_dp_phases()
            .into_iter()
            .fold(Duration::ZERO, |total, (_, elapsed)| {
                total.saturating_add(elapsed)
            })
    }

    /// Returns profile-DP orchestration not attributed to a nested phase.
    #[must_use]
    pub fn profile_dp_unattributed(self) -> Duration {
        self.profile_dp.saturating_sub(self.profile_dp_attributed())
    }

    fn observe(&mut self, other: Self) {
        self.candidate_count = self.candidate_count.saturating_add(other.candidate_count);
        self.total = self.total.saturating_add(other.total);
        self.projection = self.projection.saturating_add(other.projection);
        self.profile_preparation = self
            .profile_preparation
            .saturating_add(other.profile_preparation);
        self.fft_anchor_search = self
            .fft_anchor_search
            .saturating_add(other.fft_anchor_search);
        self.profile_dp = self.profile_dp.saturating_add(other.profile_dp);
        self.profile_dp_block_count = self
            .profile_dp_block_count
            .saturating_add(other.profile_dp_block_count);
        self.profile_dp_cell_count = self
            .profile_dp_cell_count
            .saturating_add(other.profile_dp_cell_count);
        self.profile_dp_homology_term_evaluations = self
            .profile_dp_homology_term_evaluations
            .saturating_add(other.profile_dp_homology_term_evaluations);
        for (total, observed) in self
            .profile_dp_homology_term_cell_counts
            .iter_mut()
            .zip(other.profile_dp_homology_term_cell_counts)
        {
            *total = total.saturating_add(observed);
        }
        self.profile_dp_workspace_initialization = self
            .profile_dp_workspace_initialization
            .saturating_add(other.profile_dp_workspace_initialization);
        self.profile_dp_recurrence_fill = self
            .profile_dp_recurrence_fill
            .saturating_add(other.profile_dp_recurrence_fill);
        self.profile_dp_endpoint_traceback = self
            .profile_dp_endpoint_traceback
            .saturating_add(other.profile_dp_endpoint_traceback);
        self.profile_dp_block_accounting = self
            .profile_dp_block_accounting
            .saturating_add(other.profile_dp_block_accounting);
        self.profile_dp_score_replay = self
            .profile_dp_score_replay
            .saturating_add(other.profile_dp_score_replay);
        self.profile_assembly = self.profile_assembly.saturating_add(other.profile_assembly);
        self.alignment_materialization = self
            .alignment_materialization
            .saturating_add(other.alignment_materialization);
        self.validation = self.validation.saturating_add(other.validation);
        self.wsp_scoring = self.wsp_scoring.saturating_add(other.wsp_scoring);
        self.accounting = self.accounting.saturating_add(other.accounting);
    }
}

/// One FFT-NS-i result with coarse complete-call and detailed candidate timings.
#[derive(Clone, Debug)]
pub struct FftNsiDetailedProfile<A: MafftAlphabet> {
    result: FftNsiAlignment<A>,
    timings: PhaseTimings,
    candidate_timings: FftNsiCandidatePhaseTimings,
}

impl<A: MafftAlphabet> FftNsiDetailedProfile<A> {
    /// Returns the unchanged scientific alignment result.
    #[must_use]
    pub const fn result(&self) -> &FftNsiAlignment<A> {
        &self.result
    }

    /// Returns complete-call phase timings.
    #[must_use]
    pub const fn timings(&self) -> &PhaseTimings {
        &self.timings
    }

    /// Returns aggregate detailed candidate-evaluation timings.
    #[must_use]
    pub const fn candidate_timings(&self) -> FftNsiCandidatePhaseTimings {
        self.candidate_timings
    }

    /// Separates the scientific result and both timing views.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        FftNsiAlignment<A>,
        PhaseTimings,
        FftNsiCandidatePhaseTimings,
    ) {
        (self.result, self.timings, self.candidate_timings)
    }
}

/// One auditable best-first refinement attempt.
#[derive(Clone, Debug)]
pub struct FftNsiCycle {
    cycle_index: usize,
    guide_tree: GuideTree,
    partition_count: usize,
    score_before: f64,
    best_candidate_score: Option<f64>,
    score_after: f64,
    best_partition: Option<FftNsiPartition>,
    best_candidate_merge: Option<FftNs1Merge>,
    accepted: bool,
    minimum_pair_weight: Option<f64>,
    maximum_pair_weight: Option<f64>,
}

impl FftNsiCycle {
    /// Returns the zero-based attempted-cycle index.
    #[must_use]
    pub const fn cycle_index(&self) -> usize {
        self.cycle_index
    }

    /// Returns the guide tree reconstructed from this cycle's current MSA.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        &self.guide_tree
    }

    /// Returns the number of unique tree-edge candidates examined.
    #[must_use]
    pub const fn partition_count(&self) -> usize {
        self.partition_count
    }

    /// Returns WSP for the current alignment under this cycle's fixed weights.
    #[must_use]
    pub const fn score_before(&self) -> f64 {
        self.score_before
    }

    /// Returns the highest candidate WSP, or `None` when no partition exists.
    #[must_use]
    pub const fn best_candidate_score(&self) -> Option<f64> {
        self.best_candidate_score
    }

    /// Returns the retained WSP after this attempt.
    #[must_use]
    pub const fn score_after(&self) -> f64 {
        self.score_after
    }

    /// Returns the best examined partition, whether accepted or rejected.
    #[must_use]
    pub const fn best_partition(&self) -> Option<&FftNsiPartition> {
        self.best_partition.as_ref()
    }

    /// Returns the accepted partition, or `None` for a converged attempt.
    #[must_use]
    pub fn accepted_partition(&self) -> Option<&FftNsiPartition> {
        self.accepted
            .then_some(self.best_partition.as_ref())
            .flatten()
    }

    /// Returns the best candidate's independently replayed FFT profile merge.
    #[must_use]
    pub const fn best_candidate_merge(&self) -> Option<&FftNs1Merge> {
        self.best_candidate_merge.as_ref()
    }

    /// Returns whether this attempt strictly improved its fixed WSP objective.
    #[must_use]
    pub const fn accepted(&self) -> bool {
        self.accepted
    }

    /// Returns the smallest distinct Gotoh pair weight for this cycle.
    #[must_use]
    pub const fn minimum_pair_weight(&self) -> Option<f64> {
        self.minimum_pair_weight
    }

    /// Returns the largest distinct Gotoh pair weight for this cycle.
    #[must_use]
    pub const fn maximum_pair_weight(&self) -> Option<f64> {
        self.maximum_pair_weight
    }
}

/// Major-buffer and high-water memory accounting for an FFT-NS-i call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FftNsiMemoryEstimate {
    initial_estimated_peak_bytes: usize,
    retained_initial_major_bytes: usize,
    cycle_distance_matrix_peak_bytes: usize,
    cycle_pair_weights_peak_bytes: usize,
    cycle_pair_score_cache_peak_bytes: usize,
    candidate_profile_score_peak_bytes: usize,
    candidate_fft_workspace_peak_bytes: usize,
    candidate_profile_dp_peak_bytes: usize,
    candidate_alignment_peak_bytes: usize,
    candidate_parallel_workspace_peak_bytes: usize,
    returned_final_alignment_bytes: usize,
    estimated_peak_bytes: usize,
}

impl FftNsiMemoryEstimate {
    /// Returns the complete modeled peak of the nested FFT-NS-2 call.
    #[must_use]
    pub const fn initial_estimated_peak_bytes(self) -> usize {
        self.initial_estimated_peak_bytes
    }

    /// Returns major FFT-NS-2 buffers retained while refinement runs.
    #[must_use]
    pub const fn retained_initial_major_bytes(self) -> usize {
        self.retained_initial_major_bytes
    }

    /// Returns the largest aligned-distance matrix resident in one cycle.
    #[must_use]
    pub const fn cycle_distance_matrix_peak_bytes(self) -> usize {
        self.cycle_distance_matrix_peak_bytes
    }

    /// Returns the largest dense Gotoh pair-weight matrix.
    #[must_use]
    pub const fn cycle_pair_weights_peak_bytes(self) -> usize {
        self.cycle_pair_weights_peak_bytes
    }

    /// Returns the largest cache of fixed weighted pair contributions.
    #[must_use]
    pub const fn cycle_pair_score_cache_peak_bytes(self) -> usize {
        self.cycle_pair_score_cache_peak_bytes
    }

    /// Returns the largest transformed profile-score workspace.
    #[must_use]
    pub const fn candidate_profile_score_peak_bytes(self) -> usize {
        self.candidate_profile_score_peak_bytes
    }

    /// Returns the largest transform-buffer and implementation-scratch workspace.
    #[must_use]
    pub const fn candidate_fft_workspace_peak_bytes(self) -> usize {
        self.candidate_fft_workspace_peak_bytes
    }

    /// Returns the largest candidate score/traceback workspace.
    #[must_use]
    pub const fn candidate_profile_dp_peak_bytes(self) -> usize {
        self.candidate_profile_dp_peak_bytes
    }

    /// Returns the largest typed candidate alignment buffer.
    #[must_use]
    pub const fn candidate_alignment_peak_bytes(self) -> usize {
        self.candidate_alignment_peak_bytes
    }

    /// Returns the conservative concurrent candidate-evaluation workspace.
    #[must_use]
    pub const fn candidate_parallel_workspace_peak_bytes(self) -> usize {
        self.candidate_parallel_workspace_peak_bytes
    }

    /// Returns typed cell storage in the final alignment.
    #[must_use]
    pub const fn returned_final_alignment_bytes(self) -> usize {
        self.returned_final_alignment_bytes
    }

    /// Returns the conservative major-buffer high-water estimate.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// A completed `MAFFT-2002-FFT-NS-i-v2` alignment and refinement evidence.
#[derive(Clone, Debug)]
pub struct FftNsiAlignment<A: MafftAlphabet> {
    initial: crate::FftNs2Alignment<A>,
    cycles: Box<[FftNsiCycle]>,
    termination: FftNsiTermination,
    final_distance_matrix: FftNs2DistanceMatrix,
    final_guide_tree: GuideTree,
    final_pair_weights: FftNsiPairWeights,
    final_score: f64,
    alignment: MultipleAlignment<A>,
    memory: FftNsiMemoryEstimate,
}

impl<A: MafftAlphabet> FftNsiAlignment<A> {
    /// Returns the final refined alignment in original input-row order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the complete FFT-NS-2 initial-alignment audit.
    #[must_use]
    pub const fn initial(&self) -> &crate::FftNs2Alignment<A> {
        &self.initial
    }

    /// Returns every accepted attempt and the terminal convergence attempt.
    #[must_use]
    pub const fn cycles(&self) -> &[FftNsiCycle] {
        &self.cycles
    }

    /// Returns the number of strictly improving best-first cycles.
    #[must_use]
    pub fn accepted_cycle_count(&self) -> usize {
        self.cycles.iter().filter(|cycle| cycle.accepted()).count()
    }

    /// Returns whether convergence or the cycle limit stopped refinement.
    #[must_use]
    pub const fn termination(&self) -> FftNsiTermination {
        self.termination
    }

    /// Returns distances reconstructed from the final alignment.
    #[must_use]
    pub const fn final_distance_matrix(&self) -> &FftNs2DistanceMatrix {
        &self.final_distance_matrix
    }

    /// Returns the deterministic UPGMA tree reconstructed from the final MSA.
    #[must_use]
    pub const fn final_guide_tree(&self) -> &GuideTree {
        &self.final_guide_tree
    }

    /// Returns normalized Gotoh weights for every final sequence pair.
    #[must_use]
    pub const fn final_pair_weights(&self) -> &FftNsiPairWeights {
        &self.final_pair_weights
    }

    /// Returns final WSP evaluated with the final tree-derived weights.
    #[must_use]
    pub const fn final_score(&self) -> f64 {
        self.final_score
    }

    /// Returns whether refinement changed any aligned cell from FFT-NS-2.
    #[must_use]
    pub fn changed_initial_alignment(&self) -> bool {
        self.initial.alignment() != &self.alignment
    }

    /// Returns explicit major-buffer accounting for the complete pipeline.
    #[must_use]
    pub const fn memory_estimate(&self) -> FftNsiMemoryEstimate {
        self.memory
    }

    /// Consumes the audit record and returns only the final alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.alignment
    }
}

/// Deterministic implementation of the 2002 FFT-NS-i scientific profile.
#[derive(Clone, Copy, Debug)]
pub struct MafftFftNsi<A: MafftAlphabet> {
    parameters: FftNsiParameters,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftFftNsi<A> {
    /// Creates the native-compatible `MAFFT-2002-FFT-NS-i-v2` scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: FftNsiParameters::paper_2002(),
            marker: PhantomData,
        }
    }

    /// Returns the stable identifier of the reproduced scientific profile.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        FFT_NSI_PROFILE_IDENTIFIER
    }

    /// Returns the immutable, versioned iterative parameters.
    #[must_use]
    pub const fn parameters(self) -> FftNsiParameters {
        self.parameters
    }

    /// Replaces the deterministic candidate-evaluation worker count.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError::InvalidThreadCount`] when `count` is zero.
    pub fn with_thread_count(mut self, count: usize) -> Result<Self, MafftError> {
        self.parameters = self.parameters.with_thread_count(count)?;
        Ok(self)
    }

    /// Aligns named protein, DNA, or RNA records with FFT-NS-i.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError`] for an empty collection, duplicate identifiers,
    /// checked-size overflow, non-finite computation, or an internal structural
    /// validation failure.
    #[allow(clippy::too_many_lines)]
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<FftNsiAlignment<A>, MafftError> {
        self.align_with_phase_timings::<false>(records, None, None)
    }

    /// Aligns records and reports monotonic wall time by FFT-NS-i phase.
    ///
    /// This diagnostic entry point returns the same scientific result as
    /// [`Self::align`]. Repeated refinement-cycle phase names are aggregated,
    /// and clock reads occur only on this opt-in path.
    ///
    /// # Errors
    ///
    /// Returns the same [`MafftError`] values as [`Self::align`].
    pub fn align_profiled(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Profiled<FftNsiAlignment<A>>, MafftError> {
        let mut timings = PhaseTimingsBuilder::new();
        let result = self.align_with_phase_timings::<false>(records, Some(&mut timings), None)?;
        Ok(Profiled::new(result, timings.finish()))
    }

    /// Aligns records and adds detailed candidate-evaluation phase timings.
    ///
    /// The ordinary [`Self::align`] and coarse [`Self::align_profiled`] paths do
    /// not execute these inner clocks. Candidate durations are aggregate worker
    /// work when more than one candidate worker is selected.
    ///
    /// # Errors
    ///
    /// Returns the same [`MafftError`] values as [`Self::align`].
    pub fn align_detailed_profiled(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<FftNsiDetailedProfile<A>, MafftError> {
        let mut timings = PhaseTimingsBuilder::new();
        let mut candidate_timings = FftNsiCandidatePhaseTimings::default();
        let result = self.align_with_phase_timings::<true>(
            records,
            Some(&mut timings),
            Some(&mut candidate_timings),
        )?;
        Ok(FftNsiDetailedProfile {
            result,
            timings: timings.finish(),
            candidate_timings,
        })
    }

    #[allow(clippy::too_many_lines)]
    fn align_with_phase_timings<const DETAILED: bool>(
        &self,
        records: &[SequenceRecord<A>],
        mut timings: Option<&mut PhaseTimingsBuilder>,
        mut candidate_timings: Option<&mut FftNsiCandidatePhaseTimings>,
    ) -> Result<FftNsiAlignment<A>, MafftError> {
        let initial = measure_phase(&mut timings, "initial_alignment", || {
            MafftFftNs2::<A>::new().align(records)
        })?;
        let (mut current, scorer) = measure_phase(
            &mut timings,
            "refinement_setup",
            || -> Result<_, MafftError> {
                Ok((
                    initial.alignment().clone(),
                    ScientificScorer::new::<A>(self.parameters.initial().progressive())?,
                ))
            },
        )?;
        let mut cycles = Vec::new();
        let mut accepted_score_history = Vec::new();
        let mut termination = FftNsiTermination::CycleLimit;
        let mut distance_peak = 0usize;
        let mut pair_weights_peak = 0usize;
        let mut pair_score_cache_peak = 0usize;
        let mut profile_score_peak = 0usize;
        let mut fft_peak = 0usize;
        let mut profile_dp_peak = 0usize;
        let mut candidate_alignment_peak = 0usize;
        let mut candidate_parallel_peak = 0usize;

        for cycle_index in 0..self.parameters.maximum_cycles() {
            let (
                guide_tree,
                weighting,
                pair_score_cache,
                score_before,
                minimum_pair_weight,
                maximum_pair_weight,
            ) = measure_phase(&mut timings, "cycle_setup", || -> Result<_, MafftError> {
                let distances = FftNs2DistanceMatrix::from_alignment(
                    &current,
                    self.parameters.initial().poisson_identity_floor(),
                )?;
                distance_peak = distance_peak.max(distances.resident_bytes()?);
                let guide_tree =
                    build_upgma_from_values(distances.size(), distances.as_flat_slice())?;
                let weighting = GotohWeighting::from_guide_tree(&guide_tree, self.parameters)?;
                pair_weights_peak =
                    pair_weights_peak.max(weighting.pair_weights().resident_bytes()?);
                let pair_score_cache = WspPairScoreCache::from_alignment(
                    &current,
                    weighting.pair_weights(),
                    &scorer,
                    self.parameters.initial().progressive().gap_open(),
                )?;
                pair_score_cache_peak =
                    pair_score_cache_peak.max(pair_score_cache.resident_bytes()?);
                let score_before = pair_score_cache.total();
                let minimum_pair_weight = weighting.pair_weights().minimum_pair_weight();
                let maximum_pair_weight = weighting.pair_weights().maximum_pair_weight();
                Ok((
                    guide_tree,
                    weighting,
                    pair_score_cache,
                    score_before,
                    minimum_pair_weight,
                    maximum_pair_weight,
                ))
            })?;
            let batch = measure_phase(&mut timings, "candidate_evaluation", || {
                evaluate_partitions::<A, DETAILED>(
                    weighting.partitions(),
                    &weighting,
                    records,
                    &current,
                    &pair_score_cache,
                    &scorer,
                    self.parameters,
                )
            })?;
            if DETAILED {
                candidate_timings
                    .as_deref_mut()
                    .expect("detailed FFT-NS-i profiling requires a collector")
                    .observe(batch.timings);
            }
            let should_stop = measure_phase(
                &mut timings,
                "cycle_commit",
                || -> Result<bool, MafftError> {
                    profile_score_peak = profile_score_peak.max(batch.memory.profile_score);
                    fft_peak = fft_peak.max(batch.memory.fft_workspace);
                    profile_dp_peak = profile_dp_peak.max(batch.memory.profile_dp);
                    candidate_alignment_peak =
                        candidate_alignment_peak.max(batch.memory.candidate_alignment);
                    candidate_parallel_peak = candidate_parallel_peak.max(batch.concurrent_peak);
                    let best = batch.best;
                    let partition_count = weighting.partitions().len();
                    if let Some(candidate) = best {
                        let candidate_score = candidate.score;
                        let accepted = candidate_score > score_before;
                        let oscillating = accepted
                            && repeats_native_best_first_score(
                                &accepted_score_history,
                                candidate_score,
                            );
                        let score_after = if accepted {
                            candidate_score
                        } else {
                            score_before
                        };
                        let best_candidate_score = Some(candidate_score);
                        let best_partition = Some(candidate.partition);
                        let best_candidate_merge = Some(candidate.merge);
                        if accepted {
                            current = candidate.alignment;
                            accepted_score_history.push(candidate_score);
                        }
                        cycles.push(FftNsiCycle {
                            cycle_index,
                            guide_tree,
                            partition_count,
                            score_before,
                            best_candidate_score,
                            score_after,
                            best_partition,
                            best_candidate_merge,
                            accepted,
                            minimum_pair_weight,
                            maximum_pair_weight,
                        });
                        if oscillating {
                            termination = FftNsiTermination::Oscillating;
                        } else if !accepted {
                            termination = FftNsiTermination::Converged;
                        }
                        Ok(oscillating || !accepted)
                    } else {
                        cycles.push(FftNsiCycle {
                            cycle_index,
                            guide_tree,
                            partition_count,
                            score_before,
                            best_candidate_score: None,
                            score_after: score_before,
                            best_partition: None,
                            best_candidate_merge: None,
                            accepted: false,
                            minimum_pair_weight,
                            maximum_pair_weight,
                        });
                        termination = FftNsiTermination::Converged;
                        Ok(true)
                    }
                },
            )?;
            if should_stop {
                break;
            }
        }

        let (final_distance_matrix, final_guide_tree, final_pair_weights, final_score) =
            measure_phase(
                &mut timings,
                "final_scoring",
                || -> Result<_, MafftError> {
                    let final_distance_matrix = FftNs2DistanceMatrix::from_alignment(
                        &current,
                        self.parameters.initial().poisson_identity_floor(),
                    )?;
                    distance_peak = distance_peak.max(final_distance_matrix.resident_bytes()?);
                    let final_guide_tree = build_upgma_from_values(
                        final_distance_matrix.size(),
                        final_distance_matrix.as_flat_slice(),
                    )?;
                    let final_weighting =
                        GotohWeighting::from_guide_tree(&final_guide_tree, self.parameters)?;
                    pair_weights_peak =
                        pair_weights_peak.max(final_weighting.pair_weights().resident_bytes()?);
                    let final_pair_weights = final_weighting.pair_weights().clone();
                    let final_score = weighted_sum_of_pairs(
                        &current,
                        &final_pair_weights,
                        &scorer,
                        self.parameters.initial().progressive().gap_open(),
                    )?;
                    Ok((
                        final_distance_matrix,
                        final_guide_tree,
                        final_pair_weights,
                        final_score,
                    ))
                },
            )?;

        measure_phase(
            &mut timings,
            "result_accounting",
            || -> Result<_, MafftError> {
                let initial_memory = initial.memory_estimate();
                let retained_initial_major_bytes = initial_memory
                    .retained_first_pass_major_bytes()
                    .checked_add(initial_memory.second_distance_matrix_bytes())
                    .and_then(|bytes| {
                        bytes.checked_add(initial_memory.second_sequence_weights_bytes())
                    })
                    .and_then(|bytes| {
                        bytes.checked_add(initial_memory.returned_final_alignment_bytes())
                    })
                    .ok_or(MafftError::ArithmeticOverflow {
                        stage: "FFT-NS-i retained initial memory accounting",
                    })?;
                let returned_final_alignment_bytes =
                    alignment_bytes::<A>(current.sequence_count(), current.width())?;
                let refinement_peak = retained_initial_major_bytes
                    .checked_add(returned_final_alignment_bytes)
                    .and_then(|bytes| bytes.checked_add(distance_peak))
                    .and_then(|bytes| bytes.checked_add(pair_weights_peak))
                    .and_then(|bytes| bytes.checked_add(pair_score_cache_peak))
                    .and_then(|bytes| bytes.checked_add(candidate_parallel_peak))
                    .ok_or(MafftError::ArithmeticOverflow {
                        stage: "FFT-NS-i peak memory accounting",
                    })?;
                let memory = FftNsiMemoryEstimate {
                    initial_estimated_peak_bytes: initial_memory.estimated_peak_bytes(),
                    retained_initial_major_bytes,
                    cycle_distance_matrix_peak_bytes: distance_peak,
                    cycle_pair_weights_peak_bytes: pair_weights_peak,
                    cycle_pair_score_cache_peak_bytes: pair_score_cache_peak,
                    candidate_profile_score_peak_bytes: profile_score_peak,
                    candidate_fft_workspace_peak_bytes: fft_peak,
                    candidate_profile_dp_peak_bytes: profile_dp_peak,
                    candidate_alignment_peak_bytes: candidate_alignment_peak,
                    candidate_parallel_workspace_peak_bytes: candidate_parallel_peak,
                    returned_final_alignment_bytes,
                    estimated_peak_bytes: initial_memory
                        .estimated_peak_bytes()
                        .max(refinement_peak),
                };

                Ok(FftNsiAlignment {
                    initial,
                    cycles: cycles.into_boxed_slice(),
                    termination,
                    final_distance_matrix,
                    final_guide_tree,
                    final_pair_weights,
                    final_score,
                    alignment: current,
                    memory,
                })
            },
        )
    }
}

fn measure_phase<T>(
    timings: &mut Option<&mut PhaseTimingsBuilder>,
    name: &'static str,
    operation: impl FnOnce() -> T,
) -> T {
    if let Some(timings) = timings.as_deref_mut() {
        timings.measure(name, operation)
    } else {
        operation()
    }
}

#[allow(clippy::float_cmp)]
fn repeats_native_best_first_score(history: &[f64], candidate_score: f64) -> bool {
    // MAFFT's released BESTFIRST controller deliberately starts at history[1],
    // compares doubles exactly, and retains the repeated candidate before it
    // stops. Preserve those details rather than introducing a tolerance or an
    // alignment-state policy that would define a different algorithm.
    history
        .get(1..)
        .is_some_and(|previous| previous.iter().rev().any(|score| *score == candidate_score))
}

impl<A: MafftAlphabet> Default for MafftFftNsi<A> {
    fn default() -> Self {
        Self::new()
    }
}

struct Candidate<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    partition: FftNsiPartition,
    merge: FftNs1Merge,
    score: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct CandidateMemory {
    profile_score: usize,
    fft_workspace: usize,
    profile_dp: usize,
    candidate_alignment: usize,
}

impl CandidateMemory {
    fn observe(&mut self, other: Self) {
        self.profile_score = self.profile_score.max(other.profile_score);
        self.fft_workspace = self.fft_workspace.max(other.fft_workspace);
        self.profile_dp = self.profile_dp.max(other.profile_dp);
        self.candidate_alignment = self.candidate_alignment.max(other.candidate_alignment);
    }

    fn conservative_worker_peak(self) -> Result<usize, MafftError> {
        self.profile_score
            .checked_add(self.fft_workspace)
            .and_then(|bytes| bytes.checked_add(self.profile_dp))
            .and_then(|bytes| {
                self.candidate_alignment
                    .checked_mul(2)
                    .and_then(|alignment| bytes.checked_add(alignment))
            })
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-i candidate worker memory accounting",
            })
    }
}

struct CandidateEvaluation<A: MafftAlphabet> {
    candidate: Candidate<A>,
    memory: CandidateMemory,
    timings: FftNsiCandidatePhaseTimings,
}

struct CandidateWorker<A: MafftAlphabet> {
    best: Option<Candidate<A>>,
    memory: CandidateMemory,
    timings: FftNsiCandidatePhaseTimings,
}

struct CandidateBatch<A: MafftAlphabet> {
    best: Option<Candidate<A>>,
    memory: CandidateMemory,
    concurrent_peak: usize,
    timings: FftNsiCandidatePhaseTimings,
}

#[allow(clippy::too_many_arguments)]
fn evaluate_partitions<A: MafftAlphabet, const DETAILED: bool>(
    partitions: &[TreePartition],
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    pair_score_cache: &WspPairScoreCache,
    scorer: &ScientificScorer,
    parameters: FftNsiParameters,
) -> Result<CandidateBatch<A>, MafftError> {
    if partitions.is_empty() {
        return Ok(CandidateBatch {
            best: None,
            memory: CandidateMemory::default(),
            concurrent_peak: 0,
            timings: FftNsiCandidatePhaseTimings::default(),
        });
    }
    let worker_count = parameters.thread_count().min(partitions.len());
    let workers = if worker_count == 1 {
        vec![evaluate_candidate_worker::<A, DETAILED>(
            partitions,
            weighting,
            records,
            current,
            pair_score_cache,
            scorer,
            parameters,
        )?]
    } else {
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(worker_count);
            let base_len = partitions.len() / worker_count;
            let longer_count = partitions.len() % worker_count;
            let mut start = 0usize;
            for worker_index in 0..worker_count {
                let chunk_len = base_len + usize::from(worker_index < longer_count);
                let end = start + chunk_len;
                let chunk = &partitions[start..end];
                start = end;
                handles.push(scope.spawn(move || {
                    evaluate_candidate_worker::<A, DETAILED>(
                        chunk,
                        weighting,
                        records,
                        current,
                        pair_score_cache,
                        scorer,
                        parameters,
                    )
                }));
            }
            debug_assert_eq!(start, partitions.len());
            let mut workers = Vec::with_capacity(handles.len());
            for handle in handles {
                workers.push(handle.join().map_err(|_| MafftError::InternalInvariant {
                    message: "an FFT-NS-i candidate worker panicked",
                })??);
            }
            Ok::<_, MafftError>(workers)
        })?
    };

    let mut best = None;
    let mut memory = CandidateMemory::default();
    let mut concurrent_peak = 0usize;
    let mut timings = FftNsiCandidatePhaseTimings::default();
    for worker in workers {
        memory.observe(worker.memory);
        if DETAILED {
            timings.observe(worker.timings);
        }
        concurrent_peak = concurrent_peak
            .checked_add(worker.memory.conservative_worker_peak()?)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-i parallel candidate memory accounting",
            })?;
        if let Some(candidate) = worker.best {
            let replace = best
                .as_ref()
                .is_none_or(|current_best: &Candidate<A>| candidate.score > current_best.score);
            if replace {
                best = Some(candidate);
            }
        }
    }
    Ok(CandidateBatch {
        best,
        memory,
        concurrent_peak,
        timings,
    })
}

#[allow(clippy::too_many_arguments)]
fn evaluate_candidate_worker<A: MafftAlphabet, const DETAILED: bool>(
    partitions: &[TreePartition],
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    pair_score_cache: &WspPairScoreCache,
    scorer: &ScientificScorer,
    parameters: FftNsiParameters,
) -> Result<CandidateWorker<A>, MafftError> {
    let mut best = None;
    let mut memory = CandidateMemory::default();
    let mut timings = FftNsiCandidatePhaseTimings::default();
    for partition in partitions {
        let evaluation = evaluate_candidate::<A, DETAILED>(
            partition,
            weighting,
            records,
            current,
            pair_score_cache,
            scorer,
            parameters,
        )?;
        memory.observe(evaluation.memory);
        if DETAILED {
            timings.observe(evaluation.timings);
        }
        let replace = best.as_ref().is_none_or(|current_best: &Candidate<A>| {
            evaluation.candidate.score > current_best.score
        });
        if replace {
            best = Some(evaluation.candidate);
        }
    }
    Ok(CandidateWorker {
        best,
        memory,
        timings,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn evaluate_candidate<A: MafftAlphabet, const DETAILED: bool>(
    partition: &TreePartition,
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    pair_score_cache: &WspPairScoreCache,
    scorer: &ScientificScorer,
    parameters: FftNsiParameters,
) -> Result<CandidateEvaluation<A>, MafftError> {
    let total_started = DETAILED.then(Instant::now);
    let mut timings = FftNsiCandidatePhaseTimings {
        candidate_count: 1,
        ..FftNsiCandidatePhaseTimings::default()
    };

    let projection_started = DETAILED.then(Instant::now);
    let (left_weights, right_weights) = weighting.partition_profile_weights(partition)?;
    let left = project_profile(current, partition.public.left_members(), &left_weights)?;
    let right = project_profile(current, partition.public.right_members(), &right_weights)?;
    finish_candidate_phase(projection_started, &mut timings.projection);

    let mut merges = Vec::with_capacity(1);
    let mut progressive_memory = ProgressiveMemory::default();
    let mut merge_timings = MergePhaseTimings::default();
    let candidate_profile = if DETAILED {
        merge_pair_profiled(
            left,
            right,
            scorer,
            parameters.initial().progressive(),
            &mut merges,
            &mut progressive_memory,
            &mut merge_timings,
        )?
    } else {
        merge_pair(
            left,
            right,
            scorer,
            parameters.initial().progressive(),
            &mut merges,
            &mut progressive_memory,
        )?
    };
    timings.profile_preparation = merge_timings.profile_preparation;
    timings.fft_anchor_search = merge_timings.fft_anchor_search;
    timings.profile_dp = merge_timings.profile_dp;
    timings.profile_dp_block_count = merge_timings.profile_dp_details.block_count;
    timings.profile_dp_cell_count = merge_timings.profile_dp_details.cell_count;
    timings.profile_dp_homology_term_evaluations =
        merge_timings.profile_dp_details.homology_term_evaluations;
    timings.profile_dp_homology_term_cell_counts =
        merge_timings.profile_dp_details.homology_term_cell_counts;
    timings.profile_dp_workspace_initialization =
        merge_timings.profile_dp_details.workspace_initialization;
    timings.profile_dp_recurrence_fill = merge_timings.profile_dp_details.recurrence_fill;
    timings.profile_dp_endpoint_traceback = merge_timings.profile_dp_details.endpoint_traceback;
    timings.profile_dp_block_accounting = merge_timings.profile_dp_details.block_accounting;
    timings.profile_dp_score_replay = merge_timings.profile_dp_score_replay;
    timings.profile_assembly = merge_timings.profile_assembly;

    let materialization_started = DETAILED.then(Instant::now);
    let candidate_alignment = candidate_profile.into_alignment()?;
    finish_candidate_phase(
        materialization_started,
        &mut timings.alignment_materialization,
    );

    let validation_started = DETAILED.then(Instant::now);
    validate_candidate(records, current, &candidate_alignment, &partition.public)?;
    finish_candidate_phase(validation_started, &mut timings.validation);

    let wsp_started = DETAILED.then(Instant::now);
    let candidate_score = pair_score_cache.score_candidate(
        &candidate_alignment,
        &partition.public,
        weighting.pair_weights(),
        scorer,
        parameters.initial().progressive().gap_open(),
    )?;
    finish_candidate_phase(wsp_started, &mut timings.wsp_scoring);
    #[cfg(debug_assertions)]
    {
        let direct_score = weighted_sum_of_pairs(
            &candidate_alignment,
            weighting.pair_weights(),
            scorer,
            parameters.initial().progressive().gap_open(),
        )?;
        debug_assert_eq!(candidate_score.to_bits(), direct_score.to_bits());
    }
    let accounting_started = DETAILED.then(Instant::now);
    let merge = merges.pop().ok_or(MafftError::InternalInvariant {
        message: "FFT-NS-i candidate did not emit one profile merge",
    })?;
    if !merges.is_empty() {
        return Err(MafftError::InternalInvariant {
            message: "FFT-NS-i candidate emitted multiple profile merges",
        });
    }
    let candidate_alignment_bytes = alignment_bytes::<A>(
        candidate_alignment.sequence_count(),
        candidate_alignment.width(),
    )?;
    let mut evaluation = CandidateEvaluation {
        candidate: Candidate {
            alignment: candidate_alignment,
            partition: partition.public.clone(),
            merge,
            score: candidate_score,
        },
        memory: CandidateMemory {
            profile_score: progressive_memory.profile_score,
            fft_workspace: progressive_memory.fft_workspace,
            profile_dp: progressive_memory.profile_dp,
            candidate_alignment: candidate_alignment_bytes,
        },
        timings,
    };
    finish_candidate_phase(accounting_started, &mut evaluation.timings.accounting);
    finish_candidate_phase(total_started, &mut evaluation.timings.total);
    Ok(evaluation)
}

#[inline]
fn finish_candidate_phase(started: Option<Instant>, elapsed: &mut Duration) {
    if let Some(started) = started {
        *elapsed = elapsed.saturating_add(started.elapsed());
    }
}

fn project_profile<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    members: &[usize],
    weights: &[f64],
) -> Result<WorkingProfile<A>, MafftError> {
    if members.is_empty() || members.len() != weights.len() {
        return Err(MafftError::InternalInvariant {
            message: "FFT-NS-i projected-profile dimensions disagree",
        });
    }
    let mut keep = vec![false; alignment.width()];
    for member in members.iter().copied() {
        let row = alignment.row(member).ok_or(MafftError::InternalInvariant {
            message: "FFT-NS-i partition member is out of range",
        })?;
        for (column, cell) in row.symbols().iter().enumerate() {
            keep[column] |= cell.is_some();
        }
    }
    let width = keep.iter().filter(|retain| **retain).count();
    let mut rows = Vec::with_capacity(members.len());
    for (slot, member) in members.iter().copied().enumerate() {
        let row = &alignment.rows()[member];
        let weight = weights[slot];
        if !weight.is_finite() || weight <= 0.0 {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT-NS-i projected profile weights",
            });
        }
        let cells = row
            .symbols()
            .iter()
            .copied()
            .zip(keep.iter().copied())
            .filter_map(|(cell, retain)| retain.then_some(cell))
            .collect();
        rows.push(WorkingRow {
            input_index: member,
            identifier: row.identifier().to_owned(),
            cells,
            weight,
        });
    }
    Ok(WorkingProfile { rows, width })
}

fn validate_candidate<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    previous: &MultipleAlignment<A>,
    candidate: &MultipleAlignment<A>,
    partition: &FftNsiPartition,
) -> Result<(), MafftError> {
    if candidate.sequence_count() != records.len() || previous.sequence_count() != records.len() {
        return Err(MafftError::InternalInvariant {
            message: "FFT-NS-i candidate sequence count changed",
        });
    }
    for (index, record) in records.iter().enumerate() {
        let row = &candidate.rows()[index];
        if row.identifier() != record.identifier()
            || !row
                .symbols()
                .iter()
                .filter_map(|cell| *cell)
                .eq(record.sequence().iter())
        {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-i candidate changed an input sequence",
            });
        }
    }
    for members in [partition.left_members(), partition.right_members()] {
        if !same_projected_subalignment(previous, candidate, members) {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-i candidate changed a fixed subgroup alignment",
            });
        }
    }
    Ok(())
}

fn same_projected_subalignment<A: MafftAlphabet>(
    previous: &MultipleAlignment<A>,
    candidate: &MultipleAlignment<A>,
    members: &[usize],
) -> bool {
    let mut previous_start = 0usize;
    let mut candidate_start = 0usize;
    loop {
        let previous_column = next_subgroup_column(previous, members, previous_start);
        let candidate_column = next_subgroup_column(candidate, members, candidate_start);
        match (previous_column, candidate_column) {
            (None, None) => return true,
            (Some(previous_column), Some(candidate_column)) => {
                if members.iter().copied().any(|member| {
                    previous.rows()[member].symbols()[previous_column]
                        != candidate.rows()[member].symbols()[candidate_column]
                }) {
                    return false;
                }
                previous_start = previous_column + 1;
                candidate_start = candidate_column + 1;
            }
            (None, Some(_)) | (Some(_), None) => return false,
        }
    }
}

fn next_subgroup_column<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    members: &[usize],
    start: usize,
) -> Option<usize> {
    (start..alignment.width()).find(|column| {
        members
            .iter()
            .copied()
            .any(|member| alignment.rows()[member].symbols()[*column].is_some())
    })
}

fn alignment_bytes<A: MafftAlphabet>(
    sequence_count: usize,
    width: usize,
) -> Result<usize, MafftError> {
    sequence_count
        .checked_mul(width)
        .and_then(|cells| cells.checked_mul(size_of::<Option<Symbol<A>>>()))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "FFT-NS-i alignment memory accounting",
        })
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence, Symbol};
    use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

    use super::{
        FftNsiTermination, MafftFftNsi, repeats_native_best_first_score,
        same_projected_subalignment,
    };
    use crate::{FFT_NSI_PROFILE_IDENTIFIER, MafftError, MafftFftNs2};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn aligned(rows: &[(&str, &str)]) -> MultipleAlignment<Dna> {
        MultipleAlignment::try_new(
            rows.iter()
                .map(|(identifier, text)| {
                    AlignedRow::new(
                        (*identifier).to_owned(),
                        text.bytes()
                            .map(|byte| {
                                (byte != b'-').then(|| Symbol::<Dna>::from_ascii(byte).unwrap())
                            })
                            .collect(),
                    )
                    .unwrap()
                })
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn native_best_first_history_skips_first_score_and_detects_later_repeats() {
        assert!(!repeats_native_best_first_score(&[], 4.0));
        assert!(!repeats_native_best_first_score(&[4.0], 4.0));
        assert!(!repeats_native_best_first_score(&[4.0, 5.0], 4.0));
        assert!(repeats_native_best_first_score(&[4.0, 5.0], 5.0));
        assert!(repeats_native_best_first_score(&[4.0, 5.0, 6.0], 5.0));
    }

    #[test]
    fn singleton_has_one_zero_candidate_convergence_attempt() {
        let result = MafftFftNsi::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(result.termination(), FftNsiTermination::Converged);
        assert_eq!(result.accepted_cycle_count(), 0);
        assert_eq!(result.cycles().len(), 1);
        assert_eq!(result.cycles()[0].partition_count(), 0);
        assert_eq!(result.cycles()[0].best_candidate_score(), None);
        assert!(result.final_pair_weights().sum_of_pairs().abs() < f64::EPSILON);
        assert!(result.final_score().abs() < f64::EPSILON);
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
    }

    #[test]
    fn initial_audit_is_exactly_a_direct_fft_ns2_call() {
        let records = [
            dna("a", "ACGTACGT"),
            dna("b", "ACGTCGT"),
            dna("c", "ACGTTCGT"),
            dna("d", "ACTTCGT"),
        ];
        let direct = MafftFftNs2::<Dna>::new().align(&records).unwrap();
        let iterative = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(iterative.initial().alignment(), direct.alignment());
        assert_eq!(
            iterative.initial().second_distance_matrix(),
            direct.second_distance_matrix()
        );
        assert_eq!(
            iterative.initial().second_guide_tree(),
            direct.second_guide_tree()
        );
        assert_eq!(
            iterative.initial().second_sequence_weights(),
            direct.second_sequence_weights()
        );
    }

    #[test]
    fn protein_dna_and_rna_share_the_iterative_pipeline() {
        let proteins = [
            SequenceRecord::<Protein>::new("p1", ProteinSequence::from_ascii("MKTAYIAK").unwrap())
                .unwrap(),
            SequenceRecord::<Protein>::new("p2", ProteinSequence::from_ascii("MKTIAK").unwrap())
                .unwrap(),
        ];
        let dna_records = [dna("d1", "ACGTT"), dna("d2", "ACTT")];
        let rna_records = [
            SequenceRecord::<Rna>::new("r1", RnaSequence::from_ascii("ACGUU").unwrap()).unwrap(),
            SequenceRecord::<Rna>::new("r2", RnaSequence::from_ascii("ACUU").unwrap()).unwrap(),
        ];
        assert_eq!(
            MafftFftNsi::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNsi::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNsi::<Rna>::new()
                .align(&rna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
    }

    #[test]
    fn cycles_examine_every_tree_edge_and_accept_only_strict_improvement() {
        let records = [
            dna("a", "AAAACCCCGGGG"),
            dna("b", "AAAAGGGG"),
            dna("c", "CCCCGGGG"),
            dna("d", "AAAACCCC"),
        ];
        let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        for cycle in result.cycles() {
            assert_eq!(cycle.partition_count(), 2 * records.len() - 3);
            if cycle.accepted() {
                assert!(cycle.score_after() > cycle.score_before());
                assert_eq!(
                    cycle.score_after().to_bits(),
                    cycle.best_candidate_score().unwrap().to_bits()
                );
                assert!(cycle.accepted_partition().is_some());
            } else {
                assert_eq!(
                    cycle.score_after().to_bits(),
                    cycle.score_before().to_bits()
                );
                assert!(
                    cycle
                        .best_candidate_score()
                        .is_none_or(|score| score <= cycle.score_before())
                );
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let expected_weight_sum = (records.len() * (records.len() - 1) / 2) as f64;
        assert!((result.final_pair_weights().sum_of_pairs() - expected_weight_sum).abs() < 1.0e-12);
    }

    #[test]
    fn profile_identifier_determinism_and_memory_are_explicit() {
        let records = [dna("a", "ACGTACGT"), dna("b", "ACGTCGT"), dna("c", "ACGT")];
        let aligner = MafftFftNsi::<Dna>::new();
        assert_eq!(aligner.profile_identifier(), FFT_NSI_PROFILE_IDENTIFIER);
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(first.final_pair_weights(), second.final_pair_weights());
        assert_eq!(
            first.final_score().to_bits(),
            second.final_score().to_bits()
        );
        assert_eq!(first.cycles().len(), second.cycles().len());
        for (left, right) in first.cycles().iter().zip(second.cycles()) {
            assert_eq!(
                left.score_before().to_bits(),
                right.score_before().to_bits()
            );
            assert_eq!(left.score_after().to_bits(), right.score_after().to_bits());
            assert_eq!(left.best_partition(), right.best_partition());
        }
        let memory = first.memory_estimate();
        assert!(memory.initial_estimated_peak_bytes() > 0);
        assert!(memory.cycle_distance_matrix_peak_bytes() > 0);
        assert!(memory.cycle_pair_weights_peak_bytes() > 0);
        assert!(memory.cycle_pair_score_cache_peak_bytes() > 0);
        assert!(memory.returned_final_alignment_bytes() > 0);
        assert!(memory.candidate_parallel_workspace_peak_bytes() > 0);
        assert!(memory.estimated_peak_bytes() >= memory.initial_estimated_peak_bytes());
    }

    #[test]
    fn fixed_divergent_case_strictly_improves_then_converges() {
        let records = [
            dna("s0", "TTGCCCTTTGTC"),
            dna("s1", "ATCATACA"),
            dna("s2", "TAGTGATA"),
            dna("s3", "TCGCCCCGGCAT"),
            dna("s4", "ATCGGAGCTC"),
        ];
        let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(result.accepted_cycle_count(), 1);
        assert_eq!(result.cycles().len(), 2);
        assert!(result.cycles()[0].accepted());
        assert!(result.cycles()[0].score_after() > result.cycles()[0].score_before());
        assert!(!result.cycles()[1].accepted());
        assert_eq!(result.termination(), FftNsiTermination::Converged);
        assert!(result.changed_initial_alignment());
    }

    #[test]
    fn fixed_native_best_first_repeat_stops_as_oscillation() {
        let records = [
            dna(
                "a",
                "TTCGCATTCGTGAATCGCATAGCTCCCTGCAAATAGCCTCACATGGGATCAGTACTCTCGGTAGCATGTGCCC",
            ),
            dna(
                "b",
                "TTCGCATTCGTAAATCGCATGGCTCCTGCAAATAGCCTCACATTGAATCAGTACTCTCGGTAGCATGGGCGCT",
            ),
            dna(
                "c",
                "TTCGCATTCGTGAATCGCATGGCTCCCTGCAACTAGCCTCACATGGAATCAGTACTCTGGTAGAATAGGCGCC",
            ),
            dna(
                "d",
                "TTCGCATTCTTGAATCGCATGGCTCCCTGCTAATAGCCTCACATCGAATCAATACTCTCGGTAGCATGGGCTC",
            ),
            dna(
                "e",
                "TTCGCATTCGTGAATCGCATGGCTCCCTGCAAATAGCCTCGACATGGAATCAGTACTCTCGGGAGCATCGGCA",
            ),
        ];
        let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(result.termination(), FftNsiTermination::Oscillating);
        assert_eq!(result.cycles().len(), 5);
        assert_eq!(result.accepted_cycle_count(), 5);
        let repeated = result.cycles().last().unwrap().score_after();
        assert!(
            result.cycles()[1..result.cycles().len() - 1]
                .iter()
                .any(|cycle| cycle.score_after().to_bits() == repeated.to_bits())
        );
    }

    #[test]
    fn subgroup_validation_ignores_only_inserted_all_gap_columns() {
        let previous = aligned(&[("a", "A-CG"), ("b", "ATCG"), ("c", "ACG-"), ("d", "A-GT")]);
        let preserved = aligned(&[
            ("a", "A--CG"),
            ("b", "AT-CG"),
            ("c", "A-CG-"),
            ("d", "A--GT"),
        ]);
        let changed = aligned(&[("a", "AC-G"), ("b", "ATCG"), ("c", "ACG-"), ("d", "A-GT")]);
        assert!(same_projected_subalignment(&previous, &preserved, &[0, 1]));
        assert!(!same_projected_subalignment(&previous, &changed, &[0, 1]));
    }

    #[test]
    fn parallel_candidates_are_bit_exact_and_reduce_in_partition_order() {
        let records = [
            dna("s0", "TTGCCCTTTGTC"),
            dna("s1", "ATCATACA"),
            dna("s2", "TAGTGATA"),
            dna("s3", "TCGCCCCGGCAT"),
            dna("s4", "ATCGGAGCTC"),
        ];
        let sequential = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        let parallel_aligner = MafftFftNsi::<Dna>::new().with_thread_count(4).unwrap();
        assert_eq!(parallel_aligner.parameters().thread_count(), 4);
        let parallel = parallel_aligner.align(&records).unwrap();
        assert_eq!(sequential.alignment(), parallel.alignment());
        assert_eq!(
            sequential.final_score().to_bits(),
            parallel.final_score().to_bits()
        );
        assert_eq!(
            sequential.final_pair_weights(),
            parallel.final_pair_weights()
        );
        assert_eq!(sequential.cycles().len(), parallel.cycles().len());
        for (sequential, parallel) in sequential.cycles().iter().zip(parallel.cycles()) {
            assert_eq!(sequential.best_partition(), parallel.best_partition());
            assert_eq!(
                sequential.score_before().to_bits(),
                parallel.score_before().to_bits()
            );
            assert_eq!(
                sequential.score_after().to_bits(),
                parallel.score_after().to_bits()
            );
            assert_eq!(
                sequential.best_candidate_score().map(f64::to_bits),
                parallel.best_candidate_score().map(f64::to_bits)
            );
            assert_eq!(
                sequential.best_candidate_merge(),
                parallel.best_candidate_merge()
            );
        }
        assert!(matches!(
            MafftFftNsi::<Dna>::new().with_thread_count(0),
            Err(MafftError::InvalidThreadCount { value: 0 })
        ));
    }
}
