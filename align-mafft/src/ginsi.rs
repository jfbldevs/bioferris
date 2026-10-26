use core::marker::PhantomData;
use core::mem::size_of;
use core::ops::Range;

use align_core::Symbol;
use align_msa::{DistanceMatrix, GuideTree, GuideTreeNode, MultipleAlignment, SequenceRecord};

use crate::fft_ns1::{ProgressiveMemory, merge_pair, validate_records};
use crate::generalized_affine::{GeneralizedParameters, align_generalized_affine};
use crate::local_affine::align_local_affine;
use crate::matrix::ScientificScorer;
use crate::profile::{
    PreparedProfile, ProfilePairScores, WorkingProfile, WorkingRow, align_block, merge_profiles,
    replay_operations_score,
};
use crate::refinement::{FftNsiPairWeights, FftNsiPartition, GotohWeighting, TreePartition};
use crate::upgma::{build_mixed_linkage_from_values, clustalw_sequence_weights};
use crate::wsp::{pair_alignment_score, weighted_sum_of_pairs};
use crate::{GinsiParameters, MafftAlphabet, MafftError};

/// Stable identifier for the ancestral MAFFT 5.3 G-INS-i profile.
pub const GINSI_PROFILE_IDENTIFIER: &str = "MAFFT-2005-G-INS-i-v1";

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum InsiPairwiseStage {
    GlobalFft,
    LocalAffine {
        gap_open: f64,
        offset_penalty: f64,
        gap_extension: f64,
    },
    GeneralizedAffine {
        ordinary_gap_open: f64,
        offset_penalty: f64,
        ordinary_gap_extension: f64,
        skip_gap_open: f64,
        paired_skip_extension: f64,
    },
}

/// Public name for a canonical fixed-tree G-INS-i partition.
pub type GinsiPartition = FftNsiPartition;

/// Public name for the dense fixed Gotoh pair-weight matrix.
pub type GinsiPairWeights = FftNsiPairWeights;

/// Why G-INS-i stopped iterative refinement.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GinsiTermination {
    /// A complete best-first attempt found no strict WSP+I improvement.
    Converged,
    /// The versioned maximum number of accepted attempts was reached.
    CycleLimit,
}

/// Auditable evidence for one all-pairs G/L/E-INS-i alignment.
#[derive(Clone, Debug, PartialEq)]
pub struct GinsiPairEvidence {
    first_sequence: usize,
    second_sequence: usize,
    first_aligned_start: usize,
    first_aligned_end: usize,
    second_aligned_start: usize,
    second_aligned_end: usize,
    alignment_score: f64,
    alignment_column_count: usize,
    normalization_score: f64,
    distance: f64,
    fft_anchor_count: usize,
    segment_count: usize,
    library_entry_count: usize,
    generalized_skip_opening_count: usize,
    generalized_skipped_pair_count: usize,
    generalized_unpaired_residue_count: usize,
}

impl GinsiPairEvidence {
    /// Returns the first original input index.
    #[must_use]
    pub const fn first_sequence(&self) -> usize {
        self.first_sequence
    }

    /// Returns the second original input index.
    #[must_use]
    pub const fn second_sequence(&self) -> usize {
        self.second_sequence
    }

    /// Returns the half-open first-sequence range selected by the pair stage.
    ///
    /// G-INS-i reports the complete source range. L-INS-i and E-INS-i report
    /// their selected local or generalized-local range through public aliases.
    #[must_use]
    pub fn first_aligned_range(&self) -> core::ops::Range<usize> {
        self.first_aligned_start..self.first_aligned_end
    }

    /// Returns the half-open second-sequence range selected by the pair stage.
    #[must_use]
    pub fn second_aligned_range(&self) -> core::ops::Range<usize> {
        self.second_aligned_start..self.second_aligned_end
    }

    /// Returns the independently replayed pair-alignment score.
    #[must_use]
    pub const fn alignment_score(&self) -> f64 {
        self.alignment_score
    }

    /// Returns the serialized column count in the selected pair traceback.
    ///
    /// One E-INS-i paired skip occupies two one-sided columns so its residues
    /// cannot be mistaken for a homologous pair.
    #[must_use]
    pub const fn alignment_column_count(&self) -> usize {
        self.alignment_column_count
    }

    /// Returns `min(t_ss, t_tt)`, before explicit distance totalization.
    #[must_use]
    pub const fn normalization_score(&self) -> f64 {
        self.normalization_score
    }

    /// Returns the finite normalized score distance in [0, 1].
    #[must_use]
    pub const fn distance(&self) -> f64 {
        self.distance
    }

    /// Returns the number of FFT anchors retained in this pair alignment.
    #[must_use]
    pub const fn fft_anchor_count(&self) -> usize {
        self.fft_anchor_count
    }

    /// Returns the number of maximal residue-residue segments.
    #[must_use]
    pub const fn segment_count(&self) -> usize {
        self.segment_count
    }

    /// Returns the number of sparse primary-library residue pairs.
    #[must_use]
    pub const fn library_entry_count(&self) -> usize {
        self.library_entry_count
    }

    /// Returns the number of generalized unaligned regions in the pair path.
    ///
    /// This is zero for G-INS-i and L-INS-i.
    #[must_use]
    pub const fn generalized_skip_opening_count(&self) -> usize {
        self.generalized_skip_opening_count
    }

    /// Returns residue pairs consumed without asserting homology.
    ///
    /// This is zero for G-INS-i and L-INS-i.
    #[must_use]
    pub const fn generalized_skipped_pair_count(&self) -> usize {
        self.generalized_skipped_pair_count
    }

    /// Returns one-sided residues consumed inside generalized skips.
    ///
    /// This is zero for G-INS-i and L-INS-i.
    #[must_use]
    pub const fn generalized_unpaired_residue_count(&self) -> usize {
        self.generalized_unpaired_residue_count
    }
}

/// Aggregate statistics for the primary pairwise consistency library.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GinsiLibraryStatistics {
    pair_count: usize,
    segment_count: usize,
    entry_count: usize,
    maximum_pair_entries: usize,
    resident_bytes: usize,
}

impl GinsiLibraryStatistics {
    /// Returns the number of unordered sequence-pair libraries.
    #[must_use]
    pub const fn pair_count(self) -> usize {
        self.pair_count
    }

    /// Returns the number of maximal gap-free segments.
    #[must_use]
    pub const fn segment_count(self) -> usize {
        self.segment_count
    }

    /// Returns the number of retained sparse residue-pair entries.
    #[must_use]
    pub const fn entry_count(self) -> usize {
        self.entry_count
    }

    /// Returns the largest entry count of one unordered sequence pair.
    #[must_use]
    pub const fn maximum_pair_entries(self) -> usize {
        self.maximum_pair_entries
    }

    /// Returns modeled payload bytes retained by the sparse library.
    #[must_use]
    pub const fn resident_bytes(self) -> usize {
        self.resident_bytes
    }
}

/// Evidence for one full-rectangle consistency-aware profile merge.
#[derive(Clone, Debug, PartialEq)]
pub struct GinsiMerge {
    left_members: Box<[usize]>,
    right_members: Box<[usize]>,
    left_width: usize,
    right_width: usize,
    result_width: usize,
    score: f64,
    profile_score_bytes: usize,
    profile_dp_peak_bytes: usize,
}

impl GinsiMerge {
    /// Returns sorted input indices in the left child profile.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns sorted input indices in the right child profile.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }

    /// Returns the left child width.
    #[must_use]
    pub const fn left_width(&self) -> usize {
        self.left_width
    }

    /// Returns the right child width.
    #[must_use]
    pub const fn right_width(&self) -> usize {
        self.right_width
    }

    /// Returns the merged alignment width.
    #[must_use]
    pub const fn result_width(&self) -> usize {
        self.result_width
    }

    /// Returns the independently replayed base-plus-importance DP score.
    #[must_use]
    pub const fn score(&self) -> f64 {
        self.score
    }

    /// Returns transformed base scores plus the dense temporary bonus bytes.
    #[must_use]
    pub const fn profile_score_bytes(&self) -> usize {
        self.profile_score_bytes
    }

    /// Returns the scalar full-rectangle DP high-water bytes.
    #[must_use]
    pub const fn profile_dp_peak_bytes(&self) -> usize {
        self.profile_dp_peak_bytes
    }
}

/// One auditable best-first fixed-tree refinement attempt.
#[derive(Clone, Debug)]
pub struct GinsiCycle {
    cycle_index: usize,
    partition_count: usize,
    score_before: f64,
    best_candidate_score: Option<f64>,
    score_after: f64,
    best_partition: Option<GinsiPartition>,
    best_candidate_merge: Option<GinsiMerge>,
    accepted: bool,
}

impl GinsiCycle {
    /// Returns the zero-based attempted-cycle index.
    #[must_use]
    pub const fn cycle_index(&self) -> usize {
        self.cycle_index
    }

    /// Returns the number of fixed tree-edge candidates examined.
    #[must_use]
    pub const fn partition_count(&self) -> usize {
        self.partition_count
    }

    /// Returns WSP+I before this attempt.
    #[must_use]
    pub const fn score_before(&self) -> f64 {
        self.score_before
    }

    /// Returns the best candidate score, or None for a singleton.
    #[must_use]
    pub const fn best_candidate_score(&self) -> Option<f64> {
        self.best_candidate_score
    }

    /// Returns the retained WSP+I after this attempt.
    #[must_use]
    pub const fn score_after(&self) -> f64 {
        self.score_after
    }

    /// Returns the best examined fixed-tree partition.
    #[must_use]
    pub const fn best_partition(&self) -> Option<&GinsiPartition> {
        self.best_partition.as_ref()
    }

    /// Returns the best candidate's independently replayed profile merge.
    #[must_use]
    pub const fn best_candidate_merge(&self) -> Option<&GinsiMerge> {
        self.best_candidate_merge.as_ref()
    }

    /// Returns whether this attempt strictly improved WSP+I.
    #[must_use]
    pub const fn accepted(&self) -> bool {
        self.accepted
    }
}

/// Major-buffer and high-water memory accounting for one G-INS-i call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GinsiMemoryEstimate {
    pair_alignment_peak_bytes: usize,
    pair_score_rows_peak_bytes: usize,
    pair_traceback_peak_bytes: usize,
    pair_operation_peak_bytes: usize,
    pair_profile_peak_bytes: usize,
    distance_matrix_bytes: usize,
    raw_library_peak_bytes: usize,
    primary_library_bytes: usize,
    progressive_profile_score_peak_bytes: usize,
    progressive_profile_dp_peak_bytes: usize,
    pair_weights_bytes: usize,
    objective_cache_bytes: usize,
    candidate_profile_score_peak_bytes: usize,
    candidate_profile_dp_peak_bytes: usize,
    candidate_alignment_peak_bytes: usize,
    candidate_parallel_workspace_peak_bytes: usize,
    returned_final_alignment_bytes: usize,
    estimated_peak_bytes: usize,
}

impl GinsiMemoryEstimate {
    /// Returns the largest temporary all-pairs FFT alignment workspace.
    #[must_use]
    pub const fn pair_alignment_peak_bytes(self) -> usize {
        self.pair_alignment_peak_bytes
    }

    /// Returns rolling affine-score row bytes for the largest local pair.
    ///
    /// This is zero for G-INS-i and populated through the L-INS-i alias.
    #[must_use]
    pub const fn pair_score_rows_peak_bytes(self) -> usize {
        self.pair_score_rows_peak_bytes
    }

    /// Returns packed traceback bytes for the largest local pair.
    #[must_use]
    pub const fn pair_traceback_peak_bytes(self) -> usize {
        self.pair_traceback_peak_bytes
    }

    /// Returns temporary operation-buffer bytes for the largest local pair.
    #[must_use]
    pub const fn pair_operation_peak_bytes(self) -> usize {
        self.pair_operation_peak_bytes
    }

    /// Returns padded typed pair-profile bytes for the largest local pair.
    #[must_use]
    pub const fn pair_profile_peak_bytes(self) -> usize {
        self.pair_profile_peak_bytes
    }

    /// Returns bytes in the dense normalized-score distance matrix.
    #[must_use]
    pub const fn distance_matrix_bytes(self) -> usize {
        self.distance_matrix_bytes
    }

    /// Returns the estimated peak raw segment payload before library reduction.
    #[must_use]
    pub const fn raw_library_peak_bytes(self) -> usize {
        self.raw_library_peak_bytes
    }

    /// Returns retained sparse primary-library payload bytes.
    #[must_use]
    pub const fn primary_library_bytes(self) -> usize {
        self.primary_library_bytes
    }

    /// Returns the largest progressive transformed-score workspace.
    #[must_use]
    pub const fn progressive_profile_score_peak_bytes(self) -> usize {
        self.progressive_profile_score_peak_bytes
    }

    /// Returns the largest progressive scalar DP workspace.
    #[must_use]
    pub const fn progressive_profile_dp_peak_bytes(self) -> usize {
        self.progressive_profile_dp_peak_bytes
    }

    /// Returns bytes in the dense fixed Gotoh pair weights.
    #[must_use]
    pub const fn pair_weights_bytes(self) -> usize {
        self.pair_weights_bytes
    }

    /// Returns bytes in the exact fixed-pair objective cache.
    #[must_use]
    pub const fn objective_cache_bytes(self) -> usize {
        self.objective_cache_bytes
    }

    /// Returns the largest refinement transformed-score workspace.
    #[must_use]
    pub const fn candidate_profile_score_peak_bytes(self) -> usize {
        self.candidate_profile_score_peak_bytes
    }

    /// Returns the largest refinement scalar DP workspace.
    #[must_use]
    pub const fn candidate_profile_dp_peak_bytes(self) -> usize {
        self.candidate_profile_dp_peak_bytes
    }

    /// Returns the largest typed candidate alignment buffer.
    #[must_use]
    pub const fn candidate_alignment_peak_bytes(self) -> usize {
        self.candidate_alignment_peak_bytes
    }

    /// Returns a conservative sum of concurrently resident candidate workers.
    #[must_use]
    pub const fn candidate_parallel_workspace_peak_bytes(self) -> usize {
        self.candidate_parallel_workspace_peak_bytes
    }

    /// Returns typed cell bytes in the returned final alignment.
    #[must_use]
    pub const fn returned_final_alignment_bytes(self) -> usize {
        self.returned_final_alignment_bytes
    }

    /// Returns the conservative modeled major-buffer high-water mark.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// A completed MAFFT-2005-G-INS-i-v1 alignment and its audit evidence.
#[derive(Clone, Debug)]
pub struct GinsiAlignment<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    initial_alignment: MultipleAlignment<A>,
    distance_matrix: DistanceMatrix,
    guide_tree: GuideTree,
    progressive_sequence_weights: Box<[f64]>,
    pair_evidence: Box<[GinsiPairEvidence]>,
    library_statistics: GinsiLibraryStatistics,
    progressive_merges: Box<[GinsiMerge]>,
    cycles: Box<[GinsiCycle]>,
    termination: GinsiTermination,
    pair_weights: GinsiPairWeights,
    final_score: f64,
    memory: GinsiMemoryEstimate,
}

impl<A: MafftAlphabet> GinsiAlignment<A> {
    /// Returns the final refined alignment in original input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the G-INS-1 progressive alignment before refinement.
    #[must_use]
    pub const fn initial_alignment(&self) -> &MultipleAlignment<A> {
        &self.initial_alignment
    }

    /// Returns the all-pairs normalized-score distance matrix.
    #[must_use]
    pub const fn distance_matrix(&self) -> &DistanceMatrix {
        &self.distance_matrix
    }

    /// Returns the one fixed mixed-linkage guide tree.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        &self.guide_tree
    }

    /// Returns normalized Thompson/CLUSTALW progressive weights.
    #[must_use]
    pub const fn progressive_sequence_weights(&self) -> &[f64] {
        &self.progressive_sequence_weights
    }

    /// Returns all unordered pair-alignment evidence in lexical index order.
    #[must_use]
    pub const fn pair_evidence(&self) -> &[GinsiPairEvidence] {
        &self.pair_evidence
    }

    /// Returns aggregate primary-library statistics.
    #[must_use]
    pub const fn library_statistics(&self) -> GinsiLibraryStatistics {
        self.library_statistics
    }

    /// Returns progressive consistency-aware merges in tree postorder.
    #[must_use]
    pub const fn progressive_merges(&self) -> &[GinsiMerge] {
        &self.progressive_merges
    }

    /// Returns accepted attempts and the terminal convergence attempt.
    #[must_use]
    pub const fn cycles(&self) -> &[GinsiCycle] {
        &self.cycles
    }

    /// Returns the number of strictly improving attempts.
    #[must_use]
    pub fn accepted_cycle_count(&self) -> usize {
        self.cycles.iter().filter(|cycle| cycle.accepted()).count()
    }

    /// Returns whether convergence or the cycle limit stopped refinement.
    #[must_use]
    pub const fn termination(&self) -> GinsiTermination {
        self.termination
    }

    /// Returns the fixed dense Gotoh pair weights.
    #[must_use]
    pub const fn pair_weights(&self) -> &GinsiPairWeights {
        &self.pair_weights
    }

    /// Returns the independently verified final WSP+I objective.
    #[must_use]
    pub const fn final_score(&self) -> f64 {
        self.final_score
    }

    /// Returns explicit major-buffer memory accounting.
    #[must_use]
    pub const fn memory_estimate(&self) -> GinsiMemoryEstimate {
        self.memory
    }

    /// Consumes the audit record and returns only the final alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.alignment
    }
}

/// Deterministic scalar reconstruction of ancestral MAFFT 5.3 G-INS-i.
#[derive(Clone, Copy, Debug)]
pub struct MafftGinsi<A: MafftAlphabet> {
    parameters: GinsiParameters,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftGinsi<A> {
    /// Creates the immutable MAFFT-2005-G-INS-i-v1 scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: GinsiParameters::paper_2005(),
            marker: PhantomData,
        }
    }

    /// Returns the stable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        GINSI_PROFILE_IDENTIFIER
    }

    /// Returns the versioned scientific parameters.
    #[must_use]
    pub const fn parameters(self) -> GinsiParameters {
        self.parameters
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence, Symbol};
    use align_msa::SequenceRecord;

    use super::{
        GinsiPairEvidence, GinsiTermination, ImportanceEntry, ImportanceLibrary, PairLibrary,
        RawPair, RawSegment, ResiduePair, WorkingProfile, WorkingRow, build_importance_library,
        consistency_bonus, extract_segments, normalized_score_distance,
    };
    use crate::{GINSI_PROFILE_IDENTIFIER, MafftGinsi};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn cells(text: &str) -> Vec<Option<Symbol<Dna>>> {
        text.bytes()
            .map(|byte| {
                if byte == b'-' {
                    None
                } else {
                    Some(Symbol::from_ascii(byte).unwrap())
                }
            })
            .collect()
    }

    fn evidence(first: usize, second: usize) -> GinsiPairEvidence {
        GinsiPairEvidence {
            first_sequence: first,
            second_sequence: second,
            first_aligned_start: 0,
            first_aligned_end: 1,
            second_aligned_start: 0,
            second_aligned_end: 1,
            alignment_score: 0.0,
            alignment_column_count: 1,
            normalization_score: 1.0,
            distance: 1.0,
            fft_anchor_count: 0,
            segment_count: 1,
            library_entry_count: 0,
            generalized_skip_opening_count: 0,
            generalized_skipped_pair_count: 0,
            generalized_unpaired_residue_count: 0,
        }
    }

    #[test]
    fn score_distance_totalizes_empty_and_nonpositive_cases() {
        assert_eq!(
            normalized_score_distance(0.0, 0.0, 0, 0).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            normalized_score_distance(0.0, 0.0, 0, 1).to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            normalized_score_distance(-5.0, -1.0, 1, 1).to_bits(),
            1.0_f64.to_bits()
        );
        assert!((normalized_score_distance(7.5, 10.0, 4, 4) - 0.25).abs() < f64::EPSILON);
        assert_eq!(
            normalized_score_distance(11.0, 10.0, 4, 4).to_bits(),
            0.0_f64.to_bits()
        );
    }

    #[test]
    fn maximal_segments_break_at_every_gapped_column() {
        let profile = WorkingProfile {
            rows: vec![
                WorkingRow {
                    input_index: 0,
                    identifier: "a".to_owned(),
                    cells: cells("AC-GT"),
                    weight: 1.0,
                },
                WorkingRow {
                    input_index: 1,
                    identifier: "b".to_owned(),
                    cells: cells("A-CGT"),
                    weight: 1.0,
                },
            ],
            width: 5,
        };
        let scorer = crate::matrix::ScientificScorer::ginsi_2005::<Dna>(0.123).unwrap();
        let segments = extract_segments(&profile, 0, 1, &scorer).unwrap();
        assert_eq!(segments.len(), 2);
        assert_eq!(segments[0].pairs.len(), 1);
        assert_eq!(segments[0].pairs[0].first, 0);
        assert_eq!(segments[0].pairs[0].second, 0);
        assert_eq!(
            segments[1]
                .pairs
                .iter()
                .map(|pair| (pair.first, pair.second))
                .collect::<Vec<_>>(),
            vec![(2, 2), (3, 3)]
        );
    }

    #[test]
    fn hand_computed_frequencies_and_segment_importance_match_equations() {
        let records = [dna("a", "A"), dna("b", "A"), dna("c", "A")];
        let raw = [(0, 1), (0, 2), (1, 2)]
            .into_iter()
            .map(|(first, second)| RawPair {
                first_sequence: first,
                second_sequence: second,
                segments: vec![RawSegment {
                    substitution_score: 1.0,
                    pairs: vec![ResiduePair {
                        first: 0,
                        second: 0,
                    }],
                }],
            })
            .collect();
        let mut pair_evidence = vec![evidence(0, 1), evidence(0, 2), evidence(1, 2)];
        let library =
            build_importance_library(&records, raw, &[0.2, 0.3, 0.5], &mut pair_evidence).unwrap();
        assert!((library.importance(0, 1, 0, 0) - 0.75).abs() < 1.0e-12);
        assert!((library.importance(0, 2, 0, 0) - 0.65).abs() < 1.0e-12);
        assert!((library.importance(1, 2, 0, 0) - 0.60).abs() < 1.0e-12);
        assert_eq!(library.statistics.entry_count(), 3);
        assert!(
            pair_evidence
                .iter()
                .all(|pair| pair.library_entry_count() == 1)
        );
    }

    #[test]
    fn one_cell_profile_bonus_is_weight_times_sparse_importance() {
        let mut pairs = vec![None; 4];
        pairs[1] = Some(PairLibrary {
            entries: vec![ImportanceEntry {
                first_residue: 0,
                second_residue: 0,
                value: 0.75,
            }]
            .into_boxed_slice(),
        });
        let library = ImportanceLibrary {
            size: 2,
            pairs,
            statistics: super::GinsiLibraryStatistics {
                pair_count: 1,
                segment_count: 1,
                entry_count: 1,
                maximum_pair_entries: 1,
                resident_bytes: size_of::<ImportanceEntry>(),
            },
        };
        let left = WorkingProfile::leaf(0, &dna("a", "A"), 1.0);
        let right = WorkingProfile::leaf(1, &dna("b", "A"), 1.0);
        let bonus = consistency_bonus(&left, &right, &library, 2.7).unwrap();
        assert_eq!(bonus.len(), 1);
        assert!((bonus[0] - 2.025).abs() < 1.0e-12);
    }

    #[test]
    fn singleton_is_total_and_converges_without_candidates() {
        let result = MafftGinsi::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(
            MafftGinsi::<Dna>::new().profile_identifier(),
            GINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
        assert_eq!(result.termination(), GinsiTermination::Converged);
        assert_eq!(result.cycles().len(), 1);
        assert_eq!(result.cycles()[0].partition_count(), 0);
        assert!(result.final_score().abs() < f64::EPSILON);
    }

    #[test]
    fn protein_dna_and_rna_share_the_complete_pipeline() {
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
            MafftGinsi::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftGinsi::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftGinsi::<Rna>::new()
                .align(&rna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
    }
}

impl<A: MafftAlphabet> Default for MafftGinsi<A> {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Copy, Debug)]
struct ResiduePair {
    first: usize,
    second: usize,
}

#[derive(Clone, Debug)]
struct RawSegment {
    substitution_score: f64,
    pairs: Vec<ResiduePair>,
}

#[derive(Clone, Debug)]
struct RawPair {
    first_sequence: usize,
    second_sequence: usize,
    segments: Vec<RawSegment>,
}

#[derive(Clone, Copy, Debug)]
struct ImportanceEntry {
    first_residue: usize,
    second_residue: usize,
    value: f64,
}

#[derive(Clone, Debug)]
struct PairLibrary {
    entries: Box<[ImportanceEntry]>,
}

#[derive(Clone, Debug)]
struct ImportanceLibrary {
    size: usize,
    pairs: Vec<Option<PairLibrary>>,
    statistics: GinsiLibraryStatistics,
}

impl ImportanceLibrary {
    fn pair(&self, first: usize, second: usize) -> Option<(&PairLibrary, bool)> {
        if first >= self.size || second >= self.size || first == second {
            return None;
        }
        let (low, high, reversed) = if first < second {
            (first, second, false)
        } else {
            (second, first, true)
        };
        self.pairs[low * self.size + high]
            .as_ref()
            .map(|pair| (pair, reversed))
    }

    fn importance(
        &self,
        first_sequence: usize,
        second_sequence: usize,
        first_residue: usize,
        second_residue: usize,
    ) -> f64 {
        let Some((pair, reversed)) = self.pair(first_sequence, second_sequence) else {
            return 0.0;
        };
        let key = if reversed {
            (second_residue, first_residue)
        } else {
            (first_residue, second_residue)
        };
        pair.entries
            .binary_search_by_key(&key, |entry| (entry.first_residue, entry.second_residue))
            .map_or(0.0, |index| pair.entries[index].value)
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct PairStageMemory {
    alignment_peak: usize,
    score_rows_peak: usize,
    traceback_peak: usize,
    operation_peak: usize,
    profile_peak: usize,
    raw_library: usize,
}

struct PairAlignment<A: MafftAlphabet> {
    profile: WorkingProfile<A>,
    score: f64,
    fft_anchor_count: usize,
    first_range: Range<usize>,
    second_range: Range<usize>,
    column_count: usize,
    temporary_bytes: usize,
    generalized_skip_opening_count: usize,
    generalized_skipped_pair_count: usize,
    generalized_unpaired_residue_count: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct MergeStageMemory {
    profile_score: usize,
    profile_dp: usize,
}

impl MergeStageMemory {
    fn observe(&mut self, merge: &GinsiMerge) {
        self.profile_score = self.profile_score.max(merge.profile_score_bytes);
        self.profile_dp = self.profile_dp.max(merge.profile_dp_peak_bytes);
    }
}

struct Candidate<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    partition: GinsiPartition,
    partition_index: usize,
    merge: GinsiMerge,
    score: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct CandidateMemory {
    profile_score: usize,
    profile_dp: usize,
    candidate_alignment: usize,
}

impl CandidateMemory {
    fn observe<A: MafftAlphabet>(&mut self, candidate: &Candidate<A>) -> Result<(), MafftError> {
        self.profile_score = self.profile_score.max(candidate.merge.profile_score_bytes);
        self.profile_dp = self.profile_dp.max(candidate.merge.profile_dp_peak_bytes);
        self.candidate_alignment = self.candidate_alignment.max(alignment_bytes::<A>(
            candidate.alignment.sequence_count(),
            candidate.alignment.width(),
        )?);
        Ok(())
    }

    fn merge(&mut self, other: Self) {
        self.profile_score = self.profile_score.max(other.profile_score);
        self.profile_dp = self.profile_dp.max(other.profile_dp);
        self.candidate_alignment = self.candidate_alignment.max(other.candidate_alignment);
    }

    fn conservative_worker_peak(self) -> Result<usize, MafftError> {
        self.profile_score
            .checked_add(self.profile_dp)
            .and_then(|bytes| {
                self.candidate_alignment
                    .checked_mul(2)
                    .and_then(|alignments| bytes.checked_add(alignments))
            })
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i candidate worker memory accounting",
            })
    }
}

struct CandidateWorker<A: MafftAlphabet> {
    best: Option<Candidate<A>>,
    memory: CandidateMemory,
}

struct CandidateBatch<A: MafftAlphabet> {
    best: Option<Candidate<A>>,
    memory: CandidateMemory,
    concurrent_peak: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct RefinementMemory {
    objective_cache: usize,
    profile_score: usize,
    profile_dp: usize,
    candidate_alignment: usize,
    parallel_workspace: usize,
}

struct RefinementResult<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    cycles: Vec<GinsiCycle>,
    termination: GinsiTermination,
    final_score: f64,
    memory: RefinementMemory,
}

struct ObjectiveCache {
    size: usize,
    contributions: Box<[f64]>,
    total: f64,
}

type AllPairsResult = (
    Vec<RawPair>,
    DistanceMatrix,
    Vec<GinsiPairEvidence>,
    PairStageMemory,
);

fn raw_pair_resident_bytes(segments: &[RawSegment]) -> Result<usize, MafftError> {
    let segment_bytes = segments.len().checked_mul(size_of::<RawSegment>()).ok_or(
        MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i raw segment memory accounting",
        },
    )?;
    let residue_pair_count = segments
        .iter()
        .try_fold(0usize, |total, segment| {
            total.checked_add(segment.pairs.len())
        })
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i raw residue-pair count",
        })?;
    let residue_pair_bytes = residue_pair_count
        .checked_mul(size_of::<ResiduePair>())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i raw residue-pair memory accounting",
        })?;
    size_of::<RawPair>()
        .checked_add(segment_bytes)
        .and_then(|bytes| bytes.checked_add(residue_pair_bytes))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i raw pair memory accounting",
        })
}

fn align_global_pair<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
) -> Result<PairAlignment<A>, MafftError> {
    let left = WorkingProfile::leaf(first, &records[first], 1.0);
    let right = WorkingProfile::leaf(second, &records[second], 1.0);
    let mut merges = Vec::with_capacity(1);
    let mut memory = ProgressiveMemory::default();
    let profile = merge_pair(
        left,
        right,
        scorer,
        crate::FftNs1Parameters::ginsi_2005_pairwise(),
        &mut merges,
        &mut memory,
    )?;
    let merge = merges.pop().ok_or(MafftError::InternalInvariant {
        message: "[GLE]-INS-i pair alignment produced no merge evidence",
    })?;
    let temporary_bytes = memory
        .profile_score
        .checked_add(memory.fft_workspace)
        .and_then(|bytes| bytes.checked_add(memory.profile_dp))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i pair alignment memory accounting",
        })?;
    let column_count = profile.width;
    Ok(PairAlignment {
        profile,
        score: merge.score(),
        fft_anchor_count: merge.selected_anchors().len(),
        first_range: 0..records[first].sequence().len(),
        second_range: 0..records[second].sequence().len(),
        column_count,
        temporary_bytes,
        generalized_skip_opening_count: 0,
        generalized_skipped_pair_count: 0,
        generalized_unpaired_residue_count: 0,
    })
}

fn align_local_pair<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    gap_open: f64,
    gap_extension: f64,
    memory: &mut PairStageMemory,
) -> Result<PairAlignment<A>, MafftError> {
    let local = align_local_affine(records, first, second, scorer, gap_open, gap_extension)?;
    memory.score_rows_peak = memory.score_rows_peak.max(local.memory.score_rows_bytes);
    memory.traceback_peak = memory.traceback_peak.max(local.memory.traceback_bytes);
    memory.operation_peak = memory.operation_peak.max(local.memory.operation_bytes);
    memory.profile_peak = memory.profile_peak.max(local.memory.profile_bytes);
    Ok(PairAlignment {
        profile: local.profile,
        score: local.score,
        fft_anchor_count: 0,
        first_range: local.first_range,
        second_range: local.second_range,
        column_count: local.operations.len(),
        temporary_bytes: local.memory.estimated_peak_bytes,
        generalized_skip_opening_count: 0,
        generalized_skipped_pair_count: 0,
        generalized_unpaired_residue_count: 0,
    })
}

#[allow(clippy::too_many_arguments)]
fn align_generalized_pair<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    ordinary_gap_open: f64,
    ordinary_gap_extension: f64,
    skip_gap_open: f64,
    paired_skip_extension: f64,
    memory: &mut PairStageMemory,
) -> Result<PairAlignment<A>, MafftError> {
    let generalized = align_generalized_affine(
        records,
        first,
        second,
        scorer,
        GeneralizedParameters {
            ordinary_open: ordinary_gap_open,
            ordinary_extension: ordinary_gap_extension,
            skip_open: skip_gap_open,
            paired_skip_extension,
        },
    )?;
    memory.score_rows_peak = memory
        .score_rows_peak
        .max(generalized.memory.score_rows_bytes);
    memory.traceback_peak = memory
        .traceback_peak
        .max(generalized.memory.traceback_bytes);
    memory.operation_peak = memory
        .operation_peak
        .max(generalized.memory.operation_bytes);
    memory.profile_peak = memory.profile_peak.max(generalized.memory.profile_bytes);
    let column_count = generalized
        .operations
        .len()
        .checked_add(generalized.skipped_pair_count)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "E-INS-i serialized pair column count",
        })?;
    Ok(PairAlignment {
        profile: generalized.profile,
        score: generalized.score,
        fft_anchor_count: 0,
        first_range: generalized.first_range,
        second_range: generalized.second_range,
        column_count,
        temporary_bytes: generalized.memory.estimated_peak_bytes,
        generalized_skip_opening_count: generalized.skip_opening_count,
        generalized_skipped_pair_count: generalized.skipped_pair_count,
        generalized_unpaired_residue_count: generalized.skipped_unpaired_residue_count,
    })
}

fn align_pair<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    stage: InsiPairwiseStage,
    memory: &mut PairStageMemory,
) -> Result<PairAlignment<A>, MafftError> {
    match stage {
        InsiPairwiseStage::GlobalFft => align_global_pair(records, first, second, scorer),
        InsiPairwiseStage::LocalAffine {
            gap_open,
            gap_extension,
            ..
        } => align_local_pair(
            records,
            first,
            second,
            scorer,
            gap_open,
            gap_extension,
            memory,
        ),
        InsiPairwiseStage::GeneralizedAffine {
            ordinary_gap_open,
            ordinary_gap_extension,
            skip_gap_open,
            paired_skip_extension,
            ..
        } => align_generalized_pair(
            records,
            first,
            second,
            scorer,
            ordinary_gap_open,
            ordinary_gap_extension,
            skip_gap_open,
            paired_skip_extension,
            memory,
        ),
    }
}

fn compute_all_pairs<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    group_scorer: &ScientificScorer,
    pairwise_stage: InsiPairwiseStage,
) -> Result<AllPairsResult, MafftError> {
    let count = records.len();
    let matrix_len = count
        .checked_mul(count)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i distance matrix dimensions",
        })?;
    let mut values = vec![0.0; matrix_len];
    let pair_capacity =
        count
            .checked_mul(count.saturating_sub(1))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i pair evidence capacity",
            })?
            / 2;
    let mut raw_pairs = Vec::with_capacity(pair_capacity);
    let mut evidence = Vec::with_capacity(raw_pairs.capacity());
    let local_scorer = match pairwise_stage {
        InsiPairwiseStage::GlobalFft => None,
        InsiPairwiseStage::LocalAffine { offset_penalty, .. }
        | InsiPairwiseStage::GeneralizedAffine { offset_penalty, .. } => {
            Some(ScientificScorer::ginsi_2005::<A>(-offset_penalty)?)
        }
    };
    let pair_scorer = local_scorer.as_ref().unwrap_or(group_scorer);
    let self_scores = records
        .iter()
        .map(|record| {
            record
                .sequence()
                .iter()
                .map(|symbol| {
                    pair_scorer.score_indices(symbol.alphabet_index(), symbol.alphabet_index())
                })
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let mut memory = PairStageMemory::default();

    for first in 0..count {
        for second in (first + 1)..count {
            let pair = align_pair(
                records,
                first,
                second,
                pair_scorer,
                pairwise_stage,
                &mut memory,
            )?;
            memory.alignment_peak = memory.alignment_peak.max(pair.temporary_bytes);

            validate_pair_profile(&pair.profile, records, first, second)?;
            let segments = extract_segments(&pair.profile, first, second, pair_scorer)?;
            let normalization_score = self_scores[first].min(self_scores[second]);
            let distance = normalized_score_distance(
                pair.score,
                normalization_score,
                records[first].sequence().len(),
                records[second].sequence().len(),
            );
            values[first * count + second] = distance;
            values[second * count + first] = distance;
            let raw_pair_bytes = raw_pair_resident_bytes(&segments)?;
            memory.raw_library = memory.raw_library.checked_add(raw_pair_bytes).ok_or(
                MafftError::ArithmeticOverflow {
                    stage: "[GLE]-INS-i raw library memory accounting",
                },
            )?;
            evidence.push(GinsiPairEvidence {
                first_sequence: first,
                second_sequence: second,
                first_aligned_start: pair.first_range.start,
                first_aligned_end: pair.first_range.end,
                second_aligned_start: pair.second_range.start,
                second_aligned_end: pair.second_range.end,
                alignment_score: pair.score,
                alignment_column_count: pair.column_count,
                normalization_score,
                distance,
                fft_anchor_count: pair.fft_anchor_count,
                segment_count: segments.len(),
                library_entry_count: 0,
                generalized_skip_opening_count: pair.generalized_skip_opening_count,
                generalized_skipped_pair_count: pair.generalized_skipped_pair_count,
                generalized_unpaired_residue_count: pair.generalized_unpaired_residue_count,
            });
            raw_pairs.push(RawPair {
                first_sequence: first,
                second_sequence: second,
                segments,
            });
        }
    }
    let distance_matrix = DistanceMatrix::try_from_flat(count, values)?;
    Ok((raw_pairs, distance_matrix, evidence, memory))
}

fn normalized_score_distance(
    pair_score: f64,
    normalization_score: f64,
    first_length: usize,
    second_length: usize,
) -> f64 {
    if first_length == 0 && second_length == 0 {
        0.0
    } else if first_length == 0 || second_length == 0 || normalization_score <= 0.0 {
        1.0
    } else {
        (1.0 - pair_score / normalization_score).clamp(0.0, 1.0)
    }
}

fn validate_pair_profile<A: MafftAlphabet>(
    profile: &WorkingProfile<A>,
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
) -> Result<(), MafftError> {
    if profile.rows.len() != 2
        || profile
            .rows
            .iter()
            .any(|row| row.cells.len() != profile.width)
    {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i pair profile is not rectangular",
        });
    }
    for index in [first, second] {
        let row = profile
            .rows
            .iter()
            .find(|row| row.input_index == index)
            .ok_or(MafftError::InternalInvariant {
                message: "[GLE]-INS-i pair profile lost an input row",
            })?;
        if !row
            .cells
            .iter()
            .filter_map(|cell| *cell)
            .eq(records[index].sequence().iter())
        {
            return Err(MafftError::InternalInvariant {
                message: "[GLE]-INS-i pair profile changed an input sequence",
            });
        }
    }
    Ok(())
}

fn extract_segments<A: MafftAlphabet>(
    profile: &WorkingProfile<A>,
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
) -> Result<Vec<RawSegment>, MafftError> {
    let left = profile
        .rows
        .iter()
        .find(|row| row.input_index == first)
        .ok_or(MafftError::InternalInvariant {
            message: "[GLE]-INS-i segment extraction lost the first row",
        })?;
    let right = profile
        .rows
        .iter()
        .find(|row| row.input_index == second)
        .ok_or(MafftError::InternalInvariant {
            message: "[GLE]-INS-i segment extraction lost the second row",
        })?;
    let mut segments = Vec::new();
    let mut current: Option<RawSegment> = None;
    let mut first_residue = 0usize;
    let mut second_residue = 0usize;
    for (left_cell, right_cell) in left.cells.iter().copied().zip(right.cells.iter().copied()) {
        match (left_cell, right_cell) {
            (Some(left_symbol), Some(right_symbol)) => {
                let segment = current.get_or_insert_with(|| RawSegment {
                    substitution_score: 0.0,
                    pairs: Vec::new(),
                });
                segment.substitution_score += scorer
                    .score_indices(left_symbol.alphabet_index(), right_symbol.alphabet_index());
                segment.pairs.push(ResiduePair {
                    first: first_residue,
                    second: second_residue,
                });
                first_residue += 1;
                second_residue += 1;
            }
            (left_cell, right_cell) => {
                if let Some(segment) = current.take() {
                    segments.push(segment);
                }
                first_residue += usize::from(left_cell.is_some());
                second_residue += usize::from(right_cell.is_some());
            }
        }
    }
    if let Some(segment) = current {
        segments.push(segment);
    }
    if segments
        .iter()
        .any(|segment| segment.pairs.is_empty() || !segment.substitution_score.is_finite())
    {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i gap-free segment construction",
        });
    }
    Ok(segments)
}

fn validate_unique_importance_entries(entries: &[ImportanceEntry]) -> Result<(), MafftError> {
    if entries.windows(2).any(|window| {
        window[0].first_residue == window[1].first_residue
            && window[0].second_residue == window[1].second_residue
    }) {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i primary library contains a duplicate residue pair",
        });
    }
    Ok(())
}

fn build_importance_library<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    raw_pairs: Vec<RawPair>,
    sequence_weights: &[f64],
    evidence: &mut [GinsiPairEvidence],
) -> Result<ImportanceLibrary, MafftError> {
    let size = records.len();
    if sequence_weights.len() != size || evidence.len() != raw_pairs.len() {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i library dimensions disagree",
        });
    }
    let mut frequencies = records
        .iter()
        .map(|record| vec![0.0; record.sequence().len()])
        .collect::<Vec<_>>();
    for raw in &raw_pairs {
        for segment in &raw.segments {
            for pair in &segment.pairs {
                frequencies[raw.first_sequence][pair.first] +=
                    sequence_weights[raw.second_sequence];
                frequencies[raw.second_sequence][pair.second] +=
                    sequence_weights[raw.first_sequence];
            }
        }
    }

    let pair_slots = size
        .checked_mul(size)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i primary library dimensions",
        })?;
    let mut pairs = vec![None; pair_slots];
    let mut total_segments = 0usize;
    let mut total_entries = 0usize;
    let mut maximum_pair_entries = 0usize;
    for (pair_index, raw) in raw_pairs.into_iter().enumerate() {
        let mut entries = Vec::new();
        total_segments = total_segments.checked_add(raw.segments.len()).ok_or(
            MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i segment count",
            },
        )?;
        for segment in raw.segments {
            #[allow(clippy::cast_precision_loss)]
            let denominator = 2.0 * segment.pairs.len() as f64;
            let endpoint_frequency = segment
                .pairs
                .iter()
                .map(|pair| {
                    frequencies[raw.first_sequence][pair.first]
                        + frequencies[raw.second_sequence][pair.second]
                })
                .sum::<f64>();
            let importance = segment.substitution_score * endpoint_frequency / denominator;
            if !importance.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "[GLE]-INS-i segment importance",
                });
            }
            entries.extend(segment.pairs.into_iter().map(|pair| ImportanceEntry {
                first_residue: pair.first,
                second_residue: pair.second,
                value: importance,
            }));
        }
        entries.sort_by_key(|entry| (entry.first_residue, entry.second_residue));
        validate_unique_importance_entries(&entries)?;
        maximum_pair_entries = maximum_pair_entries.max(entries.len());
        total_entries =
            total_entries
                .checked_add(entries.len())
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "[GLE]-INS-i primary library entry count",
                })?;
        evidence[pair_index].library_entry_count = entries.len();
        pairs[raw.first_sequence * size + raw.second_sequence] = Some(PairLibrary {
            entries: entries.into_boxed_slice(),
        });
    }
    let pair_index_bytes = pairs
        .len()
        .checked_mul(size_of::<Option<PairLibrary>>())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i primary library pair-index memory accounting",
        })?;
    let resident_bytes = total_entries
        .checked_mul(size_of::<ImportanceEntry>())
        .and_then(|bytes| bytes.checked_add(pair_index_bytes))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i primary library memory accounting",
        })?;
    let statistics = GinsiLibraryStatistics {
        pair_count: evidence.len(),
        segment_count: total_segments,
        entry_count: total_entries,
        maximum_pair_entries,
        resident_bytes,
    };
    Ok(ImportanceLibrary {
        size,
        pairs,
        statistics,
    })
}

fn normalized_profile_weights<A: MafftAlphabet>(
    profile: &WorkingProfile<A>,
) -> Result<Vec<f64>, MafftError> {
    let total = profile.rows.iter().map(|row| row.weight).sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i profile weight normalization",
        });
    }
    let weights = profile
        .rows
        .iter()
        .map(|row| row.weight / total)
        .collect::<Vec<_>>();
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight <= 0.0)
    {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i normalized profile weights",
        });
    }
    Ok(weights)
}

fn residue_columns<A: MafftAlphabet>(row: &WorkingRow<A>) -> Vec<usize> {
    row.cells
        .iter()
        .enumerate()
        .filter_map(|(column, cell)| cell.is_some().then_some(column))
        .collect()
}

fn consistency_bonus<A: MafftAlphabet>(
    left: &WorkingProfile<A>,
    right: &WorkingProfile<A>,
    library: &ImportanceLibrary,
    importance_weight: f64,
) -> Result<Vec<f64>, MafftError> {
    let cell_count = left
        .width
        .checked_mul(right.width)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i dense consistency bonus dimensions",
        })?;
    let mut bonus = vec![0.0; cell_count];
    let left_weights = normalized_profile_weights(left)?;
    let right_weights = normalized_profile_weights(right)?;
    let left_columns = left.rows.iter().map(residue_columns).collect::<Vec<_>>();
    let right_columns = right.rows.iter().map(residue_columns).collect::<Vec<_>>();

    for (left_slot, left_row) in left.rows.iter().enumerate() {
        for (right_slot, right_row) in right.rows.iter().enumerate() {
            let Some((pair, reversed)) = library.pair(left_row.input_index, right_row.input_index)
            else {
                return Err(MafftError::InternalInvariant {
                    message: "[GLE]-INS-i profile merge lacks a pair library",
                });
            };
            let factor = importance_weight * left_weights[left_slot] * right_weights[right_slot];
            for entry in &pair.entries {
                let (left_residue, right_residue) = if reversed {
                    (entry.second_residue, entry.first_residue)
                } else {
                    (entry.first_residue, entry.second_residue)
                };
                let left_column = *left_columns[left_slot].get(left_residue).ok_or(
                    MafftError::InternalInvariant {
                        message: "[GLE]-INS-i library left residue is out of range",
                    },
                )?;
                let right_column = *right_columns[right_slot].get(right_residue).ok_or(
                    MafftError::InternalInvariant {
                        message: "[GLE]-INS-i library right residue is out of range",
                    },
                )?;
                bonus[left_column * right.width + right_column] += factor * entry.value;
            }
        }
    }
    if bonus.iter().any(|value| !value.is_finite()) {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i dense consistency bonus",
        });
    }
    Ok(bonus)
}

fn merge_with_library<A: MafftAlphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<(WorkingProfile<A>, GinsiMerge), MafftError> {
    let left_members = left.member_indices();
    let right_members = right.member_indices();
    let left_width = left.width;
    let right_width = right.width;
    let prepared_left = PreparedProfile::new(&left, scorer)?;
    let prepared_right = PreparedProfile::new(&right, scorer)?;
    let bonus = consistency_bonus(&left, &right, library, parameters.importance_weight())?;
    let pair_scores =
        ProfilePairScores::new(&prepared_left, &prepared_right, scorer)?.with_dense_bonus(bonus)?;
    let profile_score_bytes = pair_scores.resident_bytes()?;
    let block = align_block(
        &prepared_left,
        &prepared_right,
        &pair_scores,
        0..left_width,
        0..right_width,
        parameters.gap_open(),
    )?;
    let replayed = replay_operations_score(
        &prepared_left,
        &prepared_right,
        &pair_scores,
        0..left_width,
        0..right_width,
        &block.operations,
        parameters.gap_open(),
    )?;
    let tolerance = 1.0e-9 * (1.0 + block.score.abs().max(replayed.abs()));
    if (block.score - replayed).abs() > tolerance {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i DP and independent profile replay disagree",
        });
    }
    let profile_dp_peak_bytes = block.peak_bytes;
    let merged = merge_profiles(left, right, &block.operations)?;
    let evidence = GinsiMerge {
        left_members: left_members.into_boxed_slice(),
        right_members: right_members.into_boxed_slice(),
        left_width,
        right_width,
        result_width: merged.width,
        score: replayed,
        profile_score_bytes,
        profile_dp_peak_bytes,
    };
    Ok((merged, evidence))
}

fn progressive_alignment<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    guide_tree: &GuideTree,
    sequence_weights: &[f64],
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<(MultipleAlignment<A>, Vec<GinsiMerge>, MergeStageMemory), MafftError> {
    if guide_tree.leaf_count() != records.len() || sequence_weights.len() != records.len() {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i progressive dimensions disagree",
        });
    }
    let mut merges = Vec::with_capacity(records.len().saturating_sub(1));
    let mut memory = MergeStageMemory::default();
    let profile = progressive_node(
        guide_tree.root(),
        records,
        sequence_weights,
        scorer,
        library,
        parameters,
        &mut merges,
        &mut memory,
    )?;
    let alignment = profile.into_alignment()?;
    Ok((alignment, merges, memory))
}

#[allow(clippy::too_many_arguments)]
fn progressive_node<A: MafftAlphabet>(
    node: &GuideTreeNode,
    records: &[SequenceRecord<A>],
    sequence_weights: &[f64],
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
    merges: &mut Vec<GinsiMerge>,
    memory: &mut MergeStageMemory,
) -> Result<WorkingProfile<A>, MafftError> {
    match node {
        GuideTreeNode::Leaf { sequence_index } => Ok(WorkingProfile::leaf(
            *sequence_index,
            &records[*sequence_index],
            sequence_weights[*sequence_index],
        )),
        GuideTreeNode::Branch { left, right } => {
            let left_profile = progressive_node(
                left.child(),
                records,
                sequence_weights,
                scorer,
                library,
                parameters,
                merges,
                memory,
            )?;
            let right_profile = progressive_node(
                right.child(),
                records,
                sequence_weights,
                scorer,
                library,
                parameters,
                merges,
                memory,
            )?;
            let (profile, merge) =
                merge_with_library(left_profile, right_profile, scorer, library, parameters)?;
            memory.observe(&merge);
            merges.push(merge);
            Ok(profile)
        }
    }
}

fn project_profile<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    members: &[usize],
    weights: &[f64],
) -> Result<WorkingProfile<A>, MafftError> {
    if members.is_empty() || members.len() != weights.len() {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i projected-profile dimensions disagree",
        });
    }
    let mut keep = vec![false; alignment.width()];
    for member in members.iter().copied() {
        let row = alignment.row(member).ok_or(MafftError::InternalInvariant {
            message: "[GLE]-INS-i partition member is out of range",
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
                stage: "[GLE]-INS-i projected profile weights",
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
    partition: &GinsiPartition,
) -> Result<(), MafftError> {
    if candidate.sequence_count() != records.len() || previous.sequence_count() != records.len() {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i candidate sequence count changed",
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
                message: "[GLE]-INS-i candidate changed an input sequence",
            });
        }
    }
    for members in [partition.left_members(), partition.right_members()] {
        if !same_projected_subalignment(previous, candidate, members) {
            return Err(MafftError::InternalInvariant {
                message: "[GLE]-INS-i candidate changed a fixed subgroup alignment",
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
            stage: "[GLE]-INS-i alignment memory accounting",
        })
}

fn aligned_importance<A: MafftAlphabet>(
    first_sequence: usize,
    second_sequence: usize,
    first: &[Option<Symbol<A>>],
    second: &[Option<Symbol<A>>],
    library: &ImportanceLibrary,
) -> Result<f64, MafftError> {
    if first.len() != second.len() {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i induced importance pair is not rectangular",
        });
    }
    let mut first_residue = 0usize;
    let mut second_residue = 0usize;
    let mut total = 0.0;
    for (first_cell, second_cell) in first.iter().copied().zip(second.iter().copied()) {
        if first_cell.is_some() && second_cell.is_some() {
            total += library.importance(
                first_sequence,
                second_sequence,
                first_residue,
                second_residue,
            );
        }
        first_residue += usize::from(first_cell.is_some());
        second_residue += usize::from(second_cell.is_some());
    }
    if !total.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i aligned importance score",
        });
    }
    Ok(total)
}

fn pair_objective_contribution<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    first: usize,
    second: usize,
    weights: &GinsiPairWeights,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<f64, MafftError> {
    let weight = weights
        .get(first, second)
        .ok_or(MafftError::InternalInvariant {
            message: "[GLE]-INS-i objective pair weight is out of range",
        })?;
    let first_row = alignment.rows()[first].symbols();
    let second_row = alignment.rows()[second].symbols();
    let base = pair_alignment_score(first_row, second_row, scorer, parameters.gap_open())?;
    let importance = aligned_importance(first, second, first_row, second_row, library)?;
    let contribution = weight * (base + parameters.importance_weight() * importance);
    if !contribution.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i weighted pair objective",
        });
    }
    Ok(contribution)
}

fn direct_wsp_plus_importance<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    weights: &GinsiPairWeights,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<f64, MafftError> {
    let base = weighted_sum_of_pairs(alignment, weights, scorer, parameters.gap_open())?;
    let mut importance = 0.0;
    for first in 0..alignment.sequence_count() {
        for second in 0..first {
            let weight = weights
                .get(first, second)
                .ok_or(MafftError::InternalInvariant {
                    message: "[GLE]-INS-i direct objective weight is out of range",
                })?;
            importance += weight
                * aligned_importance(
                    first,
                    second,
                    alignment.rows()[first].symbols(),
                    alignment.rows()[second].symbols(),
                    library,
                )?;
        }
    }
    let total = base + parameters.importance_weight() * importance;
    if !total.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "[GLE]-INS-i direct WSP+I objective",
        });
    }
    Ok(total)
}

impl ObjectiveCache {
    fn from_alignment<A: MafftAlphabet>(
        alignment: &MultipleAlignment<A>,
        weights: &GinsiPairWeights,
        scorer: &ScientificScorer,
        library: &ImportanceLibrary,
        parameters: GinsiParameters,
    ) -> Result<Self, MafftError> {
        let size = alignment.sequence_count();
        if weights.size() != size || library.size != size {
            return Err(MafftError::InternalInvariant {
                message: "[GLE]-INS-i objective-cache dimensions disagree",
            });
        }
        let matrix_len = size
            .checked_mul(size)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i objective-cache dimensions",
            })?;
        let mut contributions = vec![0.0; matrix_len];
        let mut total = 0.0;
        for first in 0..size {
            for second in 0..first {
                let contribution = pair_objective_contribution(
                    alignment, first, second, weights, scorer, library, parameters,
                )?;
                contributions[first * size + second] = contribution;
                total += contribution;
            }
        }
        if !total.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "[GLE]-INS-i cached WSP+I objective",
            });
        }
        Ok(Self {
            size,
            contributions: contributions.into_boxed_slice(),
            total,
        })
    }

    const fn total(&self) -> f64 {
        self.total
    }

    fn score_candidate<A: MafftAlphabet>(
        &self,
        alignment: &MultipleAlignment<A>,
        partition: &GinsiPartition,
        weights: &GinsiPairWeights,
        scorer: &ScientificScorer,
        library: &ImportanceLibrary,
        parameters: GinsiParameters,
    ) -> Result<f64, MafftError> {
        if alignment.sequence_count() != self.size {
            return Err(MafftError::InternalInvariant {
                message: "[GLE]-INS-i candidate and objective-cache dimensions disagree",
            });
        }
        let sides = partition_sides(partition, self.size)?;
        let mut total = 0.0;
        for first in 0..self.size {
            for second in 0..first {
                if sides[first] == sides[second] {
                    total += self.contributions[first * self.size + second];
                } else {
                    total += pair_objective_contribution(
                        alignment, first, second, weights, scorer, library, parameters,
                    )?;
                }
            }
        }
        if !total.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "[GLE]-INS-i cached candidate WSP+I objective",
            });
        }
        Ok(total)
    }

    fn resident_bytes(&self) -> Result<usize, MafftError> {
        self.contributions
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i objective-cache memory accounting",
            })
    }
}

fn partition_sides(partition: &GinsiPartition, size: usize) -> Result<Vec<u8>, MafftError> {
    let mut sides = vec![u8::MAX; size];
    for (side, members) in [partition.left_members(), partition.right_members()]
        .into_iter()
        .enumerate()
    {
        for member in members.iter().copied() {
            let label = sides.get_mut(member).ok_or(MafftError::InternalInvariant {
                message: "[GLE]-INS-i objective partition member is out of range",
            })?;
            if *label != u8::MAX {
                return Err(MafftError::InternalInvariant {
                    message: "[GLE]-INS-i objective partition repeats a member",
                });
            }
            *label = u8::try_from(side).expect("two partition sides fit u8");
        }
    }
    if sides.contains(&u8::MAX) {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i objective partition does not cover every row",
        });
    }
    Ok(sides)
}

#[allow(clippy::too_many_arguments)]
fn evaluate_candidate<A: MafftAlphabet>(
    partition_index: usize,
    partition: &TreePartition,
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    cache: &ObjectiveCache,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<Candidate<A>, MafftError> {
    let (left_weights, right_weights) = weighting.partition_profile_weights(partition)?;
    let left = project_profile(current, partition.public.left_members(), &left_weights)?;
    let right = project_profile(current, partition.public.right_members(), &right_weights)?;
    let (candidate_profile, merge) = merge_with_library(left, right, scorer, library, parameters)?;
    let alignment = candidate_profile.into_alignment()?;
    validate_candidate(records, current, &alignment, &partition.public)?;
    let score = cache.score_candidate(
        &alignment,
        &partition.public,
        weighting.pair_weights(),
        scorer,
        library,
        parameters,
    )?;
    Ok(Candidate {
        alignment,
        partition: partition.public.clone(),
        partition_index,
        merge,
        score,
    })
}

#[allow(clippy::too_many_arguments)]
fn evaluate_candidate_worker<A: MafftAlphabet>(
    partitions: &[TreePartition],
    first_partition: usize,
    partition_stride: usize,
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    cache: &ObjectiveCache,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
) -> Result<CandidateWorker<A>, MafftError> {
    let mut best = None;
    let mut memory = CandidateMemory::default();
    for partition_index in (first_partition..partitions.len()).step_by(partition_stride) {
        let partition = &partitions[partition_index];
        let candidate = evaluate_candidate(
            partition_index,
            partition,
            weighting,
            records,
            current,
            cache,
            scorer,
            library,
            parameters,
        )?;
        memory.observe(&candidate)?;
        let replace = best
            .as_ref()
            .is_none_or(|current_best| candidate_precedes(&candidate, current_best));
        if replace {
            best = Some(candidate);
        }
    }
    Ok(CandidateWorker { best, memory })
}

fn candidate_precedes<A: MafftAlphabet>(candidate: &Candidate<A>, current: &Candidate<A>) -> bool {
    match candidate.score.partial_cmp(&current.score) {
        Some(core::cmp::Ordering::Greater) => true,
        Some(core::cmp::Ordering::Equal) => candidate.partition_index < current.partition_index,
        Some(core::cmp::Ordering::Less) | None => false,
    }
}

#[allow(clippy::too_many_arguments)]
fn evaluate_partitions<A: MafftAlphabet>(
    partitions: &[TreePartition],
    weighting: &GotohWeighting,
    records: &[SequenceRecord<A>],
    current: &MultipleAlignment<A>,
    cache: &ObjectiveCache,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
    thread_count: usize,
) -> Result<CandidateBatch<A>, MafftError> {
    if partitions.is_empty() {
        return Ok(CandidateBatch {
            best: None,
            memory: CandidateMemory::default(),
            concurrent_peak: 0,
        });
    }
    let worker_count = thread_count.min(partitions.len());
    let workers = if worker_count == 1 {
        vec![evaluate_candidate_worker(
            partitions, 0, 1, weighting, records, current, cache, scorer, library, parameters,
        )?]
    } else {
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(worker_count);
            for worker_index in 0..worker_count {
                handles.push(scope.spawn(move || {
                    evaluate_candidate_worker(
                        partitions,
                        worker_index,
                        worker_count,
                        weighting,
                        records,
                        current,
                        cache,
                        scorer,
                        library,
                        parameters,
                    )
                }));
            }
            let mut workers = Vec::with_capacity(handles.len());
            for handle in handles {
                workers.push(handle.join().map_err(|_| MafftError::InternalInvariant {
                    message: "a [GLE]-INS-i candidate worker panicked",
                })??);
            }
            Ok::<_, MafftError>(workers)
        })?
    };

    let mut best = None;
    let mut memory = CandidateMemory::default();
    let mut concurrent_peak = 0usize;
    // Interleaving tree-partition indices balances differently sized profile
    // rectangles. The original index explicitly preserves the scalar tie rule.
    for worker in workers {
        memory.merge(worker.memory);
        concurrent_peak = concurrent_peak
            .checked_add(worker.memory.conservative_worker_peak()?)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "[GLE]-INS-i parallel candidate memory accounting",
            })?;
        if let Some(candidate) = worker.best {
            let replace = best
                .as_ref()
                .is_none_or(|current_best| candidate_precedes(&candidate, current_best));
            if replace {
                best = Some(candidate);
            }
        }
    }
    Ok(CandidateBatch {
        best,
        memory,
        concurrent_peak,
    })
}

fn scores_agree(first: f64, second: f64) -> bool {
    let tolerance = 1.0e-9 * (1.0 + first.abs().max(second.abs()));
    (first - second).abs() <= tolerance
}

fn checked_sum(stage: &'static str, values: &[usize]) -> Result<usize, MafftError> {
    values
        .iter()
        .try_fold(0usize, |total, value| total.checked_add(*value))
        .ok_or(MafftError::ArithmeticOverflow { stage })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run_refinement<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    initial: MultipleAlignment<A>,
    weighting: &GotohWeighting,
    scorer: &ScientificScorer,
    library: &ImportanceLibrary,
    parameters: GinsiParameters,
    thread_count: usize,
) -> Result<RefinementResult<A>, MafftError> {
    let mut current = initial;
    let mut cache = ObjectiveCache::from_alignment(
        &current,
        weighting.pair_weights(),
        scorer,
        library,
        parameters,
    )?;
    let direct_initial = direct_wsp_plus_importance(
        &current,
        weighting.pair_weights(),
        scorer,
        library,
        parameters,
    )?;
    if !scores_agree(cache.total(), direct_initial) {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i initial cached and direct WSP+I disagree",
        });
    }
    let mut memory = RefinementMemory {
        objective_cache: cache.resident_bytes()?,
        ..RefinementMemory::default()
    };
    let mut cycles = Vec::new();
    let mut termination = GinsiTermination::CycleLimit;

    for cycle_index in 0..parameters.maximum_cycles() {
        let score_before = cache.total();
        let batch = evaluate_partitions(
            weighting.partitions(),
            weighting,
            records,
            &current,
            &cache,
            scorer,
            library,
            parameters,
            thread_count,
        )?;
        memory.profile_score = memory.profile_score.max(batch.memory.profile_score);
        memory.profile_dp = memory.profile_dp.max(batch.memory.profile_dp);
        memory.candidate_alignment = memory
            .candidate_alignment
            .max(batch.memory.candidate_alignment);
        memory.parallel_workspace = memory.parallel_workspace.max(batch.concurrent_peak);
        let best = batch.best;
        let partition_count = weighting.partitions().len();
        let Some(best) = best else {
            cycles.push(GinsiCycle {
                cycle_index,
                partition_count,
                score_before,
                best_candidate_score: None,
                score_after: score_before,
                best_partition: None,
                best_candidate_merge: None,
                accepted: false,
            });
            termination = GinsiTermination::Converged;
            break;
        };
        let Candidate {
            alignment,
            partition,
            partition_index: _,
            merge,
            score,
        } = best;
        let accepted = score > score_before;
        let score_after = if accepted { score } else { score_before };
        cycles.push(GinsiCycle {
            cycle_index,
            partition_count,
            score_before,
            best_candidate_score: Some(score),
            score_after,
            best_partition: Some(partition),
            best_candidate_merge: Some(merge),
            accepted,
        });
        if !accepted {
            termination = GinsiTermination::Converged;
            break;
        }
        current = alignment;
        cache = ObjectiveCache::from_alignment(
            &current,
            weighting.pair_weights(),
            scorer,
            library,
            parameters,
        )?;
        memory.objective_cache = memory.objective_cache.max(cache.resident_bytes()?);
        if !scores_agree(cache.total(), score) {
            return Err(MafftError::InternalInvariant {
                message: "[GLE]-INS-i accepted candidate cache changed its WSP+I score",
            });
        }
    }
    let direct_final = direct_wsp_plus_importance(
        &current,
        weighting.pair_weights(),
        scorer,
        library,
        parameters,
    )?;
    if !scores_agree(cache.total(), direct_final) {
        return Err(MafftError::InternalInvariant {
            message: "[GLE]-INS-i final cached and direct WSP+I disagree",
        });
    }
    Ok(RefinementResult {
        alignment: current,
        cycles,
        termination,
        final_score: direct_final,
        memory,
    })
}

#[allow(clippy::too_many_lines)]
pub(crate) fn align_insi<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    parameters: GinsiParameters,
    pairwise_stage: InsiPairwiseStage,
    thread_count: usize,
) -> Result<GinsiAlignment<A>, MafftError> {
    if thread_count == 0 {
        return Err(MafftError::InvalidThreadCount {
            value: thread_count,
        });
    }
    validate_records(records)?;
    let scorer = ScientificScorer::ginsi_2005::<A>(parameters.random_match_offset())?;
    let (raw_pairs, distance_matrix, mut pair_evidence, pair_memory) =
        compute_all_pairs(records, &scorer, pairwise_stage)?;
    let guide_tree = build_mixed_linkage_from_values(
        distance_matrix.size(),
        distance_matrix.as_flat_slice(),
        parameters.mixed_linkage_fraction(),
    )?;
    let progressive_sequence_weights = clustalw_sequence_weights(&guide_tree);
    let library = build_importance_library(
        records,
        raw_pairs,
        &progressive_sequence_weights,
        &mut pair_evidence,
    )?;
    let library_statistics = library.statistics;
    let (initial_alignment, progressive_merges, progressive_memory) = progressive_alignment(
        records,
        &guide_tree,
        &progressive_sequence_weights,
        &scorer,
        &library,
        parameters,
    )?;
    let initial_alignment_bytes = alignment_bytes::<A>(
        initial_alignment.sequence_count(),
        initial_alignment.width(),
    )?;
    let weighting = GotohWeighting::from_guide_tree_values(
        &guide_tree,
        parameters.equalization_factor(),
        parameters.branch_length_floor(),
    )?;
    let pair_weights_bytes = weighting.pair_weights().resident_bytes()?;
    let pair_weights = weighting.pair_weights().clone();
    let refinement = run_refinement(
        records,
        initial_alignment.clone(),
        &weighting,
        &scorer,
        &library,
        parameters,
        thread_count,
    )?;
    let returned_final_alignment_bytes = alignment_bytes::<A>(
        refinement.alignment.sequence_count(),
        refinement.alignment.width(),
    )?;
    let distance_matrix_bytes = distance_matrix
        .as_flat_slice()
        .len()
        .checked_mul(size_of::<f64>())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "[GLE]-INS-i distance memory accounting",
        })?;
    let pair_phase = checked_sum(
        "[GLE]-INS-i pair phase memory accounting",
        &[
            distance_matrix_bytes,
            pair_memory.raw_library,
            pair_memory.alignment_peak,
        ],
    )?;
    let progressive_phase = checked_sum(
        "[GLE]-INS-i progressive phase memory accounting",
        &[
            distance_matrix_bytes,
            library_statistics.resident_bytes(),
            initial_alignment_bytes,
            progressive_memory.profile_score,
            progressive_memory.profile_dp,
        ],
    )?;
    let refinement_phase = checked_sum(
        "[GLE]-INS-i refinement phase memory accounting",
        &[
            distance_matrix_bytes,
            library_statistics.resident_bytes(),
            pair_weights_bytes,
            refinement.memory.objective_cache,
            initial_alignment_bytes,
            returned_final_alignment_bytes,
            refinement.memory.parallel_workspace,
        ],
    )?;
    let memory = GinsiMemoryEstimate {
        pair_alignment_peak_bytes: pair_memory.alignment_peak,
        pair_score_rows_peak_bytes: pair_memory.score_rows_peak,
        pair_traceback_peak_bytes: pair_memory.traceback_peak,
        pair_operation_peak_bytes: pair_memory.operation_peak,
        pair_profile_peak_bytes: pair_memory.profile_peak,
        distance_matrix_bytes,
        raw_library_peak_bytes: pair_memory.raw_library,
        primary_library_bytes: library_statistics.resident_bytes(),
        progressive_profile_score_peak_bytes: progressive_memory.profile_score,
        progressive_profile_dp_peak_bytes: progressive_memory.profile_dp,
        pair_weights_bytes,
        objective_cache_bytes: refinement.memory.objective_cache,
        candidate_profile_score_peak_bytes: refinement.memory.profile_score,
        candidate_profile_dp_peak_bytes: refinement.memory.profile_dp,
        candidate_alignment_peak_bytes: refinement.memory.candidate_alignment,
        candidate_parallel_workspace_peak_bytes: refinement.memory.parallel_workspace,
        returned_final_alignment_bytes,
        estimated_peak_bytes: pair_phase.max(progressive_phase).max(refinement_phase),
    };
    Ok(GinsiAlignment {
        alignment: refinement.alignment,
        initial_alignment,
        distance_matrix,
        guide_tree,
        progressive_sequence_weights: progressive_sequence_weights.into_boxed_slice(),
        pair_evidence: pair_evidence.into_boxed_slice(),
        library_statistics,
        progressive_merges: progressive_merges.into_boxed_slice(),
        cycles: refinement.cycles.into_boxed_slice(),
        termination: refinement.termination,
        pair_weights,
        final_score: refinement.final_score,
        memory,
    })
}

impl<A: MafftAlphabet> MafftGinsi<A> {
    /// Aligns named protein, DNA, or RNA records with ancestral G-INS-i.
    ///
    /// # Errors
    ///
    /// Returns `MafftError` for invalid input, checked-size overflow, non-finite
    /// scientific calculations, or a violated alignment invariant.
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<GinsiAlignment<A>, MafftError> {
        align_insi(records, self.parameters, InsiPairwiseStage::GlobalFft, 1)
    }
}
