use core::marker::PhantomData;
use core::mem::size_of;

use align_msa::{GuideTree, MultipleAlignment, SequenceRecord};

use crate::aligned_distance::FftNs2DistanceMatrix;
use crate::fft_ns1::{FftNs1Alignment, FftNs1Anchor, FftNs1Merge, run_progressive_pass};
use crate::matrix::ScientificScorer;
use crate::upgma::{build_upgma_from_values, clustalw_sequence_weights};
use crate::{FftNs2Parameters, MafftAlphabet, MafftError, MafftFftNs1};

/// Stable identifier for the immutable 2002 FFT-NS-2 scientific profile.
pub const FFT_NS2_PROFILE_IDENTIFIER: &str = "MAFFT-2002-FFT-NS-2-v1";

/// FFT anchor evidence shared by either progressive pass of FFT-NS-2.
pub type FftNs2Anchor = FftNs1Anchor;

/// One auditable profile merge from either progressive pass of FFT-NS-2.
pub type FftNs2Merge = FftNs1Merge;

/// Major-buffer and high-water memory accounting for an FFT-NS-2 call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FftNs2MemoryEstimate {
    first_pass_estimated_peak_bytes: usize,
    retained_first_pass_major_bytes: usize,
    second_distance_matrix_bytes: usize,
    second_sequence_weights_bytes: usize,
    second_profile_score_peak_bytes: usize,
    second_fft_workspace_peak_bytes: usize,
    second_profile_dp_peak_bytes: usize,
    returned_final_alignment_bytes: usize,
    estimated_peak_bytes: usize,
}

impl FftNs2MemoryEstimate {
    /// Returns the complete major-buffer peak modeled by the nested first pass.
    #[must_use]
    pub const fn first_pass_estimated_peak_bytes(self) -> usize {
        self.first_pass_estimated_peak_bytes
    }

    /// Returns retained first-alignment, first-distance, and first-weight storage.
    #[must_use]
    pub const fn retained_first_pass_major_bytes(self) -> usize {
        self.retained_first_pass_major_bytes
    }

    /// Returns second-distance values plus comparable and identical count arrays.
    #[must_use]
    pub const fn second_distance_matrix_bytes(self) -> usize {
        self.second_distance_matrix_bytes
    }

    /// Returns storage for second-tree sequence weights.
    #[must_use]
    pub const fn second_sequence_weights_bytes(self) -> usize {
        self.second_sequence_weights_bytes
    }

    /// Returns the largest transformed profile-score workspace in pass two.
    #[must_use]
    pub const fn second_profile_score_peak_bytes(self) -> usize {
        self.second_profile_score_peak_bytes
    }

    /// Returns the largest transform-buffer and implementation-scratch workspace in pass two.
    #[must_use]
    pub const fn second_fft_workspace_peak_bytes(self) -> usize {
        self.second_fft_workspace_peak_bytes
    }

    /// Returns the largest partition score/traceback workspace in pass two.
    #[must_use]
    pub const fn second_profile_dp_peak_bytes(self) -> usize {
        self.second_profile_dp_peak_bytes
    }

    /// Returns typed cell storage in the final alignment.
    #[must_use]
    pub const fn returned_final_alignment_bytes(self) -> usize {
        self.returned_final_alignment_bytes
    }

    /// Returns the larger modeled high-water mark across the two phases.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// A completed `MAFFT-2002-FFT-NS-2-v1` alignment and both passes' evidence.
#[derive(Clone, Debug)]
pub struct FftNs2Alignment<A: MafftAlphabet> {
    first_pass: FftNs1Alignment<A>,
    second_distance_matrix: FftNs2DistanceMatrix,
    second_guide_tree: GuideTree,
    second_sequence_weights: Box<[f64]>,
    second_merges: Box<[FftNs2Merge]>,
    second_root_merge_score: Option<f64>,
    alignment: MultipleAlignment<A>,
    memory: FftNs2MemoryEstimate,
}

impl<A: MafftAlphabet> FftNs2Alignment<A> {
    /// Returns the final second-pass alignment in original input-row order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the complete, independently inspectable FFT-NS-1 first pass.
    #[must_use]
    pub const fn first_pass(&self) -> &FftNs1Alignment<A> {
        &self.first_pass
    }

    /// Returns Poisson-corrected distances inferred from the first alignment.
    #[must_use]
    pub const fn second_distance_matrix(&self) -> &FftNs2DistanceMatrix {
        &self.second_distance_matrix
    }

    /// Returns the guide tree reconstructed from the first alignment.
    #[must_use]
    pub const fn second_guide_tree(&self) -> &GuideTree {
        &self.second_guide_tree
    }

    /// Returns normalized path weights recomputed from the second guide tree.
    #[must_use]
    pub const fn second_sequence_weights(&self) -> &[f64] {
        &self.second_sequence_weights
    }

    /// Returns second-pass progressive merges in deterministic postorder.
    #[must_use]
    pub const fn second_merges(&self) -> &[FftNs2Merge] {
        &self.second_merges
    }

    /// Returns the last second-pass profile merge score, or `None` for a singleton.
    #[must_use]
    pub const fn second_root_merge_score(&self) -> Option<f64> {
        self.second_root_merge_score
    }

    /// Returns whether guide-tree topology or branch lengths changed after pass one.
    #[must_use]
    pub fn guide_tree_changed(&self) -> bool {
        self.first_pass.guide_tree() != &self.second_guide_tree
    }

    /// Returns explicit major-buffer accounting across both passes.
    #[must_use]
    pub const fn memory_estimate(&self) -> FftNs2MemoryEstimate {
        self.memory
    }

    /// Consumes the audit record and returns only the final alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.alignment
    }
}

/// Deterministic implementation of the 2002 two-pass FFT-NS-2 pipeline.
#[derive(Clone, Copy, Debug)]
pub struct MafftFftNs2<A: MafftAlphabet> {
    parameters: FftNs2Parameters,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftFftNs2<A> {
    /// Creates the exact `MAFFT-2002-FFT-NS-2-v1` scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: FftNs2Parameters::paper_2002(),
            marker: PhantomData,
        }
    }

    /// Returns the stable identifier of the reproduced scientific profile.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        FFT_NS2_PROFILE_IDENTIFIER
    }

    /// Returns the immutable, versioned two-pass parameters.
    #[must_use]
    pub const fn parameters(self) -> FftNs2Parameters {
        self.parameters
    }

    /// Aligns named protein, DNA, or RNA records with two progressive passes.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError`] for an empty collection, duplicate identifiers,
    /// checked-size overflow, non-finite computation, or an internal structural
    /// validation failure.
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<FftNs2Alignment<A>, MafftError> {
        let first_pass = MafftFftNs1::<A>::new().align(records)?;
        let second_distance_matrix = FftNs2DistanceMatrix::from_alignment(
            first_pass.alignment(),
            self.parameters.poisson_identity_floor(),
        )?;
        let second_guide_tree = build_upgma_from_values(
            second_distance_matrix.size(),
            second_distance_matrix.as_flat_slice(),
        )?;
        let second_sequence_weights = clustalw_sequence_weights(&second_guide_tree);
        let scorer = ScientificScorer::new::<A>(self.parameters.progressive())?;
        let second = run_progressive_pass(
            records,
            &second_guide_tree,
            &second_sequence_weights,
            &scorer,
            self.parameters.progressive(),
        )?;

        let first_memory = first_pass.memory_estimate();
        let first_weight_bytes =
            records
                .len()
                .checked_mul(size_of::<f64>())
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "FFT-NS-2 first weight memory accounting",
                })?;
        let retained_first_pass_major_bytes = first_memory
            .distance_matrix_bytes()
            .checked_add(first_memory.returned_alignment_bytes())
            .and_then(|bytes| bytes.checked_add(first_weight_bytes))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-2 retained first-pass memory accounting",
            })?;
        let second_distance_matrix_bytes = second_distance_matrix.resident_bytes()?;
        let second_sequence_weights_bytes =
            records
                .len()
                .checked_mul(size_of::<f64>())
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "FFT-NS-2 second weight memory accounting",
                })?;
        let second_temporary_bytes = second
            .memory
            .profile_score
            .checked_add(second.memory.fft_workspace)
            .and_then(|bytes| bytes.checked_add(second.memory.profile_dp))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-2 second-pass temporary memory accounting",
            })?;
        let second_phase_peak = retained_first_pass_major_bytes
            .checked_add(second_distance_matrix_bytes)
            .and_then(|bytes| bytes.checked_add(second_sequence_weights_bytes))
            .and_then(|bytes| bytes.checked_add(second.returned_alignment_bytes))
            .and_then(|bytes| bytes.checked_add(second_temporary_bytes))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-2 second-phase peak memory accounting",
            })?;
        let memory = FftNs2MemoryEstimate {
            first_pass_estimated_peak_bytes: first_memory.estimated_peak_bytes(),
            retained_first_pass_major_bytes,
            second_distance_matrix_bytes,
            second_sequence_weights_bytes,
            second_profile_score_peak_bytes: second.memory.profile_score,
            second_fft_workspace_peak_bytes: second.memory.fft_workspace,
            second_profile_dp_peak_bytes: second.memory.profile_dp,
            returned_final_alignment_bytes: second.returned_alignment_bytes,
            estimated_peak_bytes: first_memory.estimated_peak_bytes().max(second_phase_peak),
        };

        Ok(FftNs2Alignment {
            first_pass,
            second_distance_matrix,
            second_guide_tree,
            second_sequence_weights: second_sequence_weights.into_boxed_slice(),
            second_merges: second.merges.into_boxed_slice(),
            second_root_merge_score: second.root_merge_score,
            alignment: second.alignment,
            memory,
        })
    }
}

impl<A: MafftAlphabet> Default for MafftFftNs2<A> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::SequenceRecord;

    use super::MafftFftNs2;

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn singleton_runs_two_zero_merge_passes() {
        let result = MafftFftNs2::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(result.first_pass().merges().len(), 0);
        assert_eq!(result.second_merges().len(), 0);
        assert_eq!(result.second_root_merge_score(), None);
        assert_eq!(result.second_distance_matrix().get(0, 0), Some(0.0));
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
    }

    #[test]
    fn protein_dna_and_rna_share_both_passes() {
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
            MafftFftNs2::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNs2::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftFftNs2::<Rna>::new()
                .align(&rna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
    }

    #[test]
    fn profile_identifier_and_memory_are_explicit() {
        let records = [dna("a", "ACGTACGT"), dna("b", "ACGTCGT"), dna("c", "ACGT")];
        let aligner = MafftFftNs2::<Dna>::new();
        assert_eq!(aligner.profile_identifier(), "MAFFT-2002-FFT-NS-2-v1");
        let result = aligner.align(&records).unwrap();
        let memory = result.memory_estimate();
        assert!(memory.first_pass_estimated_peak_bytes() > 0);
        assert!(memory.second_distance_matrix_bytes() >= 9 * size_of::<f64>());
        assert!(memory.returned_final_alignment_bytes() > 0);
        assert!(memory.estimated_peak_bytes() >= memory.first_pass_estimated_peak_bytes());
    }
}
