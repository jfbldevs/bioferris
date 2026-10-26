use core::marker::PhantomData;
use core::mem::size_of;
use std::collections::BTreeSet;
use std::time::{Duration, Instant};

use align_core::Symbol;
use align_msa::{
    DistanceMatrix, GuideTree, GuideTreeNode, MsaError, MultipleAlignment, SequenceRecord,
};

use crate::fft::{best_lags, correlate_channels};
use crate::matrix::ScientificScorer;
use crate::profile::{
    PreparedProfile, ProfileDpPhaseTimings, ProfileOperation, ProfilePairScores, WorkingProfile,
    align_block_profiled_without_bonus, align_block_without_bonus, merge_profiles,
    replay_operations_score_without_bonus,
};
use crate::segment::{find_segments, select_monotone_chain};
use crate::tuple::TupleCollection;
use crate::upgma::{build_upgma, clustalw_sequence_weights};
use crate::{FftNs1Parameters, MafftAlphabet, MafftError};

/// Stable identifier for the immutable 2002 FFT-NS-1 scientific profile.
pub const FFT_NS1_PROFILE_IDENTIFIER: &str = "MAFFT-2002-FFT-NS-1-v1";

/// One homologous FFT segment retained by the segment-level dynamic program.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FftNs1Anchor {
    left_start: usize,
    right_start: usize,
    length: usize,
    score: f64,
}

impl FftNs1Anchor {
    /// Returns the zero-based start column in the left child profile.
    #[must_use]
    pub const fn left_start(self) -> usize {
        self.left_start
    }

    /// Returns the zero-based start column in the right child profile.
    #[must_use]
    pub const fn right_start(self) -> usize {
        self.right_start
    }

    /// Returns the retained homologous segment length.
    #[must_use]
    pub const fn len(self) -> usize {
        self.length
    }

    /// Returns whether the retained segment is empty.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.length == 0
    }

    /// Returns the segment's accumulated normalized homology score.
    #[must_use]
    pub const fn score(self) -> f64 {
        self.score
    }

    /// Returns the fixed center-column pair used to divide the profile matrix.
    #[must_use]
    pub const fn center(self) -> (usize, usize) {
        (
            self.left_start + self.length / 2,
            self.right_start + self.length / 2,
        )
    }
}

/// Auditable evidence for one progressive FFT-NS-1 profile merge.
#[derive(Clone, Debug, PartialEq)]
pub struct FftNs1Merge {
    left_members: Box<[usize]>,
    right_members: Box<[usize]>,
    left_width: usize,
    right_width: usize,
    correlation_peak_count: usize,
    candidate_segment_count: usize,
    selected_anchors: Box<[FftNs1Anchor]>,
    operation_count: usize,
    result_width: usize,
    score: f64,
    fft_workspace_bytes: usize,
    profile_dp_peak_bytes: usize,
}

impl FftNs1Merge {
    /// Returns sorted original input indices in the left child profile.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns sorted original input indices in the right child profile.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }

    /// Returns the left child profile width before this merge.
    #[must_use]
    pub const fn left_width(&self) -> usize {
        self.left_width
    }

    /// Returns the right child profile width before this merge.
    #[must_use]
    pub const fn right_width(&self) -> usize {
        self.right_width
    }

    /// Returns the number of highest FFT lags inspected.
    #[must_use]
    pub const fn correlation_peak_count(&self) -> usize {
        self.correlation_peak_count
    }

    /// Returns the number of homologous segments before monotone chaining.
    #[must_use]
    pub const fn candidate_segment_count(&self) -> usize {
        self.candidate_segment_count
    }

    /// Returns the monotone homologous segments used to divide the DP matrix.
    #[must_use]
    pub const fn selected_anchors(&self) -> &[FftNs1Anchor] {
        &self.selected_anchors
    }

    /// Returns whether FFT evidence divided this profile matrix.
    #[must_use]
    pub const fn used_fft_partition(&self) -> bool {
        !self.selected_anchors.is_empty()
    }

    /// Returns the number of output path operations.
    #[must_use]
    pub const fn operation_count(&self) -> usize {
        self.operation_count
    }

    /// Returns the merged profile width.
    #[must_use]
    pub const fn result_width(&self) -> usize {
        self.result_width
    }

    /// Returns the independently replayed profile-DP score.
    #[must_use]
    pub const fn score(&self) -> f64 {
        self.score
    }

    /// Returns transform buffers and implementation scratch at this merge.
    #[must_use]
    pub const fn fft_workspace_bytes(&self) -> usize {
        self.fft_workspace_bytes
    }

    /// Returns the largest partition traceback workspace at this merge.
    #[must_use]
    pub const fn profile_dp_peak_bytes(&self) -> usize {
        self.profile_dp_peak_bytes
    }
}

/// Explicit major-buffer and high-water memory accounting for an FFT-NS-1 call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FftNs1MemoryEstimate {
    tuple_storage_bytes: usize,
    distance_matrix_bytes: usize,
    profile_score_peak_bytes: usize,
    fft_workspace_peak_bytes: usize,
    profile_dp_peak_bytes: usize,
    returned_alignment_bytes: usize,
    estimated_peak_bytes: usize,
}

impl FftNs1MemoryEstimate {
    /// Returns sparse sorted compressed-tuple storage.
    #[must_use]
    pub const fn tuple_storage_bytes(self) -> usize {
        self.tuple_storage_bytes
    }

    /// Returns dense rough-distance matrix storage.
    #[must_use]
    pub const fn distance_matrix_bytes(self) -> usize {
        self.distance_matrix_bytes
    }

    /// Returns the largest transformed profile-score workspace.
    #[must_use]
    pub const fn profile_score_peak_bytes(self) -> usize {
        self.profile_score_peak_bytes
    }

    /// Returns the largest transform-buffer and implementation-scratch workspace.
    #[must_use]
    pub const fn fft_workspace_peak_bytes(self) -> usize {
        self.fft_workspace_peak_bytes
    }

    /// Returns the largest partition score/traceback workspace.
    #[must_use]
    pub const fn profile_dp_peak_bytes(self) -> usize {
        self.profile_dp_peak_bytes
    }

    /// Returns typed cell storage in the returned alignment, excluding allocator overhead.
    #[must_use]
    pub const fn returned_alignment_bytes(self) -> usize {
        self.returned_alignment_bytes
    }

    /// Returns the sum of the modeled persistent arrays and temporary high-water marks.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// A completed `MAFFT-2002-FFT-NS-1-v1` alignment and its audit evidence.
#[derive(Clone, Debug)]
pub struct FftNs1Alignment<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    distance_matrix: DistanceMatrix,
    guide_tree: GuideTree,
    sequence_weights: Box<[f64]>,
    effective_tuple_length: usize,
    merges: Box<[FftNs1Merge]>,
    root_merge_score: Option<f64>,
    memory: FftNs1MemoryEstimate,
}

impl<A: MafftAlphabet> FftNs1Alignment<A> {
    /// Returns the final alignment in original input-row order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the publication's rough shared-tuple distance matrix.
    #[must_use]
    pub const fn distance_matrix(&self) -> &DistanceMatrix {
        &self.distance_matrix
    }

    /// Returns the single UPGMA guide tree.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        &self.guide_tree
    }

    /// Returns normalized CLUSTALW path weights in original input order.
    #[must_use]
    pub const fn sequence_weights(&self) -> &[f64] {
        &self.sequence_weights
    }

    /// Returns the compressed tuple length actually used by the input collection.
    #[must_use]
    pub const fn effective_tuple_length(&self) -> usize {
        self.effective_tuple_length
    }

    /// Returns progressive merges in deterministic postorder.
    #[must_use]
    pub const fn merges(&self) -> &[FftNs1Merge] {
        &self.merges
    }

    /// Returns the last profile merge score, or `None` for a singleton input.
    #[must_use]
    pub const fn root_merge_score(&self) -> Option<f64> {
        self.root_merge_score
    }

    /// Returns explicit memory accounting for this call.
    #[must_use]
    pub const fn memory_estimate(&self) -> FftNs1MemoryEstimate {
        self.memory
    }

    /// Consumes the audit record and returns only the alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.alignment
    }
}

/// Deterministic implementation of the 2002 single-pass FFT-NS-1 pipeline.
#[derive(Clone, Copy, Debug)]
pub struct MafftFftNs1<A: MafftAlphabet> {
    parameters: FftNs1Parameters,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftFftNs1<A> {
    /// Creates the exact `MAFFT-2002-FFT-NS-1-v1` scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: FftNs1Parameters::paper_2002(),
            marker: PhantomData,
        }
    }

    /// Returns the stable identifier of the reproduced scientific profile.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        FFT_NS1_PROFILE_IDENTIFIER
    }

    /// Returns the immutable versioned scientific parameters.
    #[must_use]
    pub const fn parameters(self) -> FftNs1Parameters {
        self.parameters
    }

    /// Aligns named protein, DNA, or RNA records with one FFT-NS-1 pass.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError`] for an empty collection, duplicate identifiers,
    /// checked-size overflow, non-finite computation, or an internal structural
    /// validation failure.
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<FftNs1Alignment<A>, MafftError> {
        validate_records(records)?;
        let scorer = ScientificScorer::new::<A>(self.parameters)?;
        let tuples = TupleCollection::from_records(records, self.parameters)?;
        let effective_tuple_length = tuples.effective_length();
        let distance_matrix = tuples.distance_matrix(records)?;
        let guide_tree = build_upgma(&distance_matrix)?;
        let sequence_weights = clustalw_sequence_weights(&guide_tree);
        let progressive = run_progressive_pass(
            records,
            &guide_tree,
            &sequence_weights,
            &scorer,
            self.parameters,
        )?;
        let distance_matrix_bytes = records
            .len()
            .checked_mul(records.len())
            .and_then(|cells| cells.checked_mul(size_of::<f64>()))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-1 distance memory accounting",
            })?;
        let persistent = tuples
            .stored_bytes()
            .checked_add(distance_matrix_bytes)
            .and_then(|bytes| bytes.checked_add(progressive.returned_alignment_bytes))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-1 persistent memory accounting",
            })?;
        let temporary = progressive
            .memory
            .profile_score
            .checked_add(progressive.memory.fft_workspace)
            .and_then(|bytes| bytes.checked_add(progressive.memory.profile_dp))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-1 temporary memory accounting",
            })?;
        let estimated_peak_bytes =
            persistent
                .checked_add(temporary)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "FFT-NS-1 peak memory accounting",
                })?;
        let memory = FftNs1MemoryEstimate {
            tuple_storage_bytes: tuples.stored_bytes(),
            distance_matrix_bytes,
            profile_score_peak_bytes: progressive.memory.profile_score,
            fft_workspace_peak_bytes: progressive.memory.fft_workspace,
            profile_dp_peak_bytes: progressive.memory.profile_dp,
            returned_alignment_bytes: progressive.returned_alignment_bytes,
            estimated_peak_bytes,
        };
        Ok(FftNs1Alignment {
            alignment: progressive.alignment,
            distance_matrix,
            guide_tree,
            sequence_weights: sequence_weights.into_boxed_slice(),
            effective_tuple_length,
            merges: progressive.merges.into_boxed_slice(),
            root_merge_score: progressive.root_merge_score,
            memory,
        })
    }
}

impl<A: MafftAlphabet> Default for MafftFftNs1<A> {
    fn default() -> Self {
        Self::new()
    }
}

pub(crate) struct ProgressivePass<A: MafftAlphabet> {
    pub(crate) alignment: MultipleAlignment<A>,
    pub(crate) merges: Vec<FftNs1Merge>,
    pub(crate) root_merge_score: Option<f64>,
    pub(crate) memory: ProgressiveMemory,
    pub(crate) returned_alignment_bytes: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ProgressiveMemory {
    pub(crate) profile_score: usize,
    pub(crate) fft_workspace: usize,
    pub(crate) profile_dp: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct MergePhaseTimings {
    pub(crate) profile_preparation: Duration,
    pub(crate) fft_anchor_search: Duration,
    pub(crate) profile_dp: Duration,
    pub(crate) profile_dp_details: ProfileDpPhaseTimings,
    pub(crate) profile_dp_score_replay: Duration,
    pub(crate) profile_assembly: Duration,
}

pub(crate) fn run_progressive_pass<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    guide_tree: &GuideTree,
    sequence_weights: &[f64],
    scorer: &ScientificScorer,
    parameters: FftNs1Parameters,
) -> Result<ProgressivePass<A>, MafftError> {
    if guide_tree.leaf_count() != records.len() || sequence_weights.len() != records.len() {
        return Err(MafftError::InternalInvariant {
            message: "progressive pass dimensions disagree",
        });
    }
    let mut merges = Vec::with_capacity(records.len().saturating_sub(1));
    let mut memory = ProgressiveMemory::default();
    let profile = align_node(
        guide_tree.root(),
        records,
        sequence_weights,
        scorer,
        parameters,
        &mut merges,
        &mut memory,
    )?;
    let root_merge_score = merges.last().map(FftNs1Merge::score);
    let returned_alignment_bytes = profile
        .rows
        .len()
        .checked_mul(profile.width)
        .and_then(|cells| cells.checked_mul(size_of::<Option<Symbol<A>>>()))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "returned progressive alignment accounting",
        })?;
    let alignment = profile.into_alignment()?;
    Ok(ProgressivePass {
        alignment,
        merges,
        root_merge_score,
        memory,
        returned_alignment_bytes,
    })
}

#[allow(clippy::too_many_arguments)]
fn align_node<A: MafftAlphabet>(
    node: &GuideTreeNode,
    records: &[SequenceRecord<A>],
    sequence_weights: &[f64],
    scorer: &ScientificScorer,
    parameters: FftNs1Parameters,
    merges: &mut Vec<FftNs1Merge>,
    memory: &mut ProgressiveMemory,
) -> Result<WorkingProfile<A>, MafftError> {
    match node {
        GuideTreeNode::Leaf { sequence_index } => Ok(WorkingProfile::leaf(
            *sequence_index,
            &records[*sequence_index],
            sequence_weights[*sequence_index],
        )),
        GuideTreeNode::Branch { left, right } => {
            let left_profile = align_node(
                left.child(),
                records,
                sequence_weights,
                scorer,
                parameters,
                merges,
                memory,
            )?;
            let right_profile = align_node(
                right.child(),
                records,
                sequence_weights,
                scorer,
                parameters,
                merges,
                memory,
            )?;
            merge_pair(
                left_profile,
                right_profile,
                scorer,
                parameters,
                merges,
                memory,
            )
        }
    }
}

#[allow(clippy::too_many_lines)]
pub(crate) fn merge_pair<A: MafftAlphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    scorer: &ScientificScorer,
    parameters: FftNs1Parameters,
    merges: &mut Vec<FftNs1Merge>,
    memory: &mut ProgressiveMemory,
) -> Result<WorkingProfile<A>, MafftError> {
    merge_pair_with_phase_timings::<A, false>(
        left,
        right,
        scorer,
        parameters,
        merges,
        memory,
        &mut MergePhaseTimings::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn merge_pair_profiled<A: MafftAlphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    scorer: &ScientificScorer,
    parameters: FftNs1Parameters,
    merges: &mut Vec<FftNs1Merge>,
    memory: &mut ProgressiveMemory,
    timings: &mut MergePhaseTimings,
) -> Result<WorkingProfile<A>, MafftError> {
    merge_pair_with_phase_timings::<A, true>(
        left, right, scorer, parameters, merges, memory, timings,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn merge_pair_with_phase_timings<A: MafftAlphabet, const PROFILED: bool>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    scorer: &ScientificScorer,
    parameters: FftNs1Parameters,
    merges: &mut Vec<FftNs1Merge>,
    memory: &mut ProgressiveMemory,
    timings: &mut MergePhaseTimings,
) -> Result<WorkingProfile<A>, MafftError> {
    let preparation_started = PROFILED.then(Instant::now);
    let left_members = left.member_indices();
    let right_members = right.member_indices();
    let left_width = left.width;
    let right_width = right.width;
    let prepared_left = PreparedProfile::new(&left, scorer)?;
    let prepared_right = PreparedProfile::new(&right, scorer)?;
    let pair_scores = ProfilePairScores::new(&prepared_left, &prepared_right, scorer)?;
    memory.profile_score = memory.profile_score.max(pair_scores.resident_bytes()?);
    finish_merge_phase(preparation_started, &mut timings.profile_preparation);

    let fft_started = PROFILED.then(Instant::now);
    let correlations = correlate_channels(prepared_left.signals(), prepared_right.signals())?;
    memory.fft_workspace = memory.fft_workspace.max(correlations.workspace_bytes);
    let peaks = best_lags(&correlations.values, parameters.peak_count());
    let candidates = find_segments(
        &prepared_left,
        &prepared_right,
        &pair_scores,
        &peaks,
        parameters,
    );
    let selected = select_monotone_chain(&candidates);
    let public_anchors = selected
        .iter()
        .map(|segment| FftNs1Anchor {
            left_start: segment.left_start,
            right_start: segment.right_start,
            length: segment.length,
            score: segment.score,
        })
        .collect::<Vec<_>>();
    finish_merge_phase(fft_started, &mut timings.fft_anchor_search);

    let dp_started = PROFILED.then(Instant::now);
    let operation_capacity =
        left_width
            .checked_add(right_width)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile operation capacity",
            })?;
    let mut operations = Vec::with_capacity(operation_capacity);
    let mut score = 0.0;
    let mut previous_left = 0usize;
    let mut previous_right = 0usize;
    let mut dp_peak_bytes = 0usize;
    for segment in &selected {
        let (center_left, center_right) = segment.center();
        if center_left < previous_left || center_right < previous_right {
            return Err(MafftError::InternalInvariant {
                message: "selected FFT segment centers are not monotone",
            });
        }
        let block = if PROFILED {
            align_block_profiled_without_bonus(
                &prepared_left,
                &prepared_right,
                &pair_scores,
                previous_left..center_left,
                previous_right..center_right,
                parameters.gap_open(),
                &mut timings.profile_dp_details,
            )?
        } else {
            align_block_without_bonus(
                &prepared_left,
                &prepared_right,
                &pair_scores,
                previous_left..center_left,
                previous_right..center_right,
                parameters.gap_open(),
            )?
        };
        score += block.score;
        dp_peak_bytes = dp_peak_bytes.max(block.peak_bytes);
        operations.extend(block.operations);
        operations.push(ProfileOperation::Diagonal);
        score += pair_scores.homology_without_bonus(center_left, center_right);
        previous_left = center_left + 1;
        previous_right = center_right + 1;
    }
    let trailing = if PROFILED {
        align_block_profiled_without_bonus(
            &prepared_left,
            &prepared_right,
            &pair_scores,
            previous_left..left_width,
            previous_right..right_width,
            parameters.gap_open(),
            &mut timings.profile_dp_details,
        )?
    } else {
        align_block_without_bonus(
            &prepared_left,
            &prepared_right,
            &pair_scores,
            previous_left..left_width,
            previous_right..right_width,
            parameters.gap_open(),
        )?
    };
    score += trailing.score;
    dp_peak_bytes = dp_peak_bytes.max(trailing.peak_bytes);
    operations.extend(trailing.operations);
    let replay_started = PROFILED.then(Instant::now);
    let replayed = replay_operations_score_without_bonus(
        &prepared_left,
        &prepared_right,
        &pair_scores,
        0..left_width,
        0..right_width,
        &operations,
        parameters.gap_open(),
    )?;
    finish_merge_phase(replay_started, &mut timings.profile_dp_score_replay);
    let tolerance = 1.0e-9 * (1.0 + score.abs().max(replayed.abs()));
    if (score - replayed).abs() > tolerance {
        return Err(MafftError::InternalInvariant {
            message: "independent profile score replay disagrees with dynamic programming",
        });
    }
    score = replayed;
    memory.profile_dp = memory.profile_dp.max(dp_peak_bytes);
    finish_merge_phase(dp_started, &mut timings.profile_dp);

    let assembly_started = PROFILED.then(Instant::now);
    let operation_count = operations.len();
    let merged_profile = merge_profiles(left, right, &operations)?;
    merges.push(FftNs1Merge {
        left_members: left_members.into_boxed_slice(),
        right_members: right_members.into_boxed_slice(),
        left_width,
        right_width,
        correlation_peak_count: peaks.len(),
        candidate_segment_count: candidates.len(),
        selected_anchors: public_anchors.into_boxed_slice(),
        operation_count,
        result_width: merged_profile.width,
        score,
        fft_workspace_bytes: correlations.workspace_bytes,
        profile_dp_peak_bytes: dp_peak_bytes,
    });
    finish_merge_phase(assembly_started, &mut timings.profile_assembly);
    Ok(merged_profile)
}

#[inline]
fn finish_merge_phase(started: Option<Instant>, elapsed: &mut Duration) {
    if let Some(started) = started {
        *elapsed = elapsed.saturating_add(started.elapsed());
    }
}

pub(crate) fn validate_records<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
) -> Result<(), MafftError> {
    if records.is_empty() {
        return Err(MsaError::EmptyCollection.into());
    }
    let mut identifiers = BTreeSet::new();
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

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::{MsaError, SequenceRecord};

    use super::MafftFftNs1;
    use crate::MafftError;

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn empty_collection_is_rejected() {
        let error = MafftFftNs1::<Dna>::new().align(&[]).unwrap_err();
        assert_eq!(error, MafftError::Msa(MsaError::EmptyCollection));
    }

    #[test]
    fn singleton_is_returned_without_a_merge() {
        let result = MafftFftNs1::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
        assert!(result.merges().is_empty());
        assert_eq!(result.root_merge_score(), None);
        assert_eq!(result.sequence_weights(), &[1.0]);
    }

    #[test]
    fn duplicate_identifiers_are_rejected_before_alignment() {
        let error = MafftFftNs1::<Dna>::new()
            .align(&[dna("same", "A"), dna("same", "C")])
            .unwrap_err();
        assert!(matches!(
            error,
            MafftError::Msa(MsaError::DuplicateIdentifier { .. })
        ));
    }

    #[test]
    fn protein_dna_and_rna_share_the_pipeline() {
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
            MafftFftNs1::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNs1::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNs1::<Rna>::new()
                .align(&rna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
    }

    #[test]
    fn output_is_deterministic_and_preserves_input_row_order() {
        let records = [dna("z", "ACGTTAC"), dna("a", "ACGTAC"), dna("m", "AGTTAC")];
        let aligner = MafftFftNs1::<Dna>::new();
        assert_eq!(aligner.profile_identifier(), "MAFFT-2002-FFT-NS-1-v1");
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        let first_rows = first
            .alignment()
            .rows()
            .iter()
            .map(|row| (row.identifier(), row.to_ascii_with_gaps()))
            .collect::<Vec<_>>();
        let second_rows = second
            .alignment()
            .rows()
            .iter()
            .map(|row| (row.identifier(), row.to_ascii_with_gaps()))
            .collect::<Vec<_>>();
        assert_eq!(first_rows, second_rows);
        assert_eq!(
            first
                .alignment()
                .rows()
                .iter()
                .map(align_msa::AlignedRow::identifier)
                .collect::<Vec<_>>(),
            vec!["z", "a", "m"]
        );
    }

    #[test]
    fn long_identical_sequences_activate_fft_partitioning() {
        let sequence = "ACGT".repeat(20);
        let records = [dna("a", &sequence), dna("b", &sequence)];
        let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
        assert!(result.merges()[0].candidate_segment_count() > 0);
        assert!(result.merges()[0].used_fft_partition());
        assert_eq!(result.alignment().width(), sequence.len());
    }
}
