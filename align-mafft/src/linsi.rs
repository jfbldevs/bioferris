use core::marker::PhantomData;

use align_msa::{DistanceMatrix, GuideTree, MultipleAlignment, SequenceRecord};

use crate::ginsi::{InsiPairwiseStage, align_insi};
use crate::{
    GinsiAlignment, GinsiCycle, GinsiLibraryStatistics, GinsiMemoryEstimate, GinsiMerge,
    GinsiPairEvidence, GinsiPairWeights, GinsiPartition, GinsiTermination, LinsiParameters,
    MafftAlphabet, MafftError,
};

/// Stable identifier for the ancestral scalar L-INS-i scientific profile.
pub const LINSI_PROFILE_IDENTIFIER: &str = "MAFFT-2005-L-INS-i-v1";

/// Auditable local-pair evidence used by L-INS-i.
pub type LinsiPairEvidence = GinsiPairEvidence;

/// Aggregate statistics for the primary local consistency library.
pub type LinsiLibraryStatistics = GinsiLibraryStatistics;

/// Evidence for one complete-profile consistency-aware merge.
pub type LinsiMerge = GinsiMerge;

/// Evidence for one best-first fixed-tree refinement attempt.
pub type LinsiCycle = GinsiCycle;

/// Major-buffer memory accounting for one L-INS-i call.
pub type LinsiMemoryEstimate = GinsiMemoryEstimate;

/// Why L-INS-i stopped fixed-tree iterative refinement.
pub type LinsiTermination = GinsiTermination;

/// Canonical fixed-tree partition used during L-INS-i refinement.
pub type LinsiPartition = GinsiPartition;

/// Dense fixed Gotoh pair-weight matrix used by L-INS-i.
pub type LinsiPairWeights = GinsiPairWeights;

/// A completed `MAFFT-2005-L-INS-i-v1` alignment and its audit evidence.
#[derive(Clone, Debug)]
pub struct LinsiAlignment<A: MafftAlphabet> {
    inner: GinsiAlignment<A>,
}

impl<A: MafftAlphabet> LinsiAlignment<A> {
    /// Returns the final full-length refined alignment in input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        self.inner.alignment()
    }

    /// Returns the full-length L-INS-1 alignment before refinement.
    #[must_use]
    pub const fn initial_alignment(&self) -> &MultipleAlignment<A> {
        self.inner.initial_alignment()
    }

    /// Returns the all-pairs normalized local-score distance matrix.
    #[must_use]
    pub const fn distance_matrix(&self) -> &DistanceMatrix {
        self.inner.distance_matrix()
    }

    /// Returns the one fixed mixed-linkage guide tree.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        self.inner.guide_tree()
    }

    /// Returns normalized Thompson/CLUSTALW progressive weights.
    #[must_use]
    pub const fn progressive_sequence_weights(&self) -> &[f64] {
        self.inner.progressive_sequence_weights()
    }

    /// Returns local pair scores, source ranges, distances, and segment counts.
    #[must_use]
    pub const fn pair_evidence(&self) -> &[LinsiPairEvidence] {
        self.inner.pair_evidence()
    }

    /// Returns aggregate primary local-library statistics.
    #[must_use]
    pub const fn library_statistics(&self) -> LinsiLibraryStatistics {
        self.inner.library_statistics()
    }

    /// Returns complete-profile progressive merges in tree postorder.
    #[must_use]
    pub const fn progressive_merges(&self) -> &[LinsiMerge] {
        self.inner.progressive_merges()
    }

    /// Returns accepted attempts and the terminal convergence attempt.
    #[must_use]
    pub const fn cycles(&self) -> &[LinsiCycle] {
        self.inner.cycles()
    }

    /// Returns the number of strictly improving refinement attempts.
    #[must_use]
    pub fn accepted_cycle_count(&self) -> usize {
        self.inner.accepted_cycle_count()
    }

    /// Returns whether convergence or the cycle limit stopped refinement.
    #[must_use]
    pub const fn termination(&self) -> LinsiTermination {
        self.inner.termination()
    }

    /// Returns the fixed dense Gotoh pair weights.
    #[must_use]
    pub const fn pair_weights(&self) -> &LinsiPairWeights {
        self.inner.pair_weights()
    }

    /// Returns the independently verified final complete-MSA WSP+I objective.
    #[must_use]
    pub const fn final_score(&self) -> f64 {
        self.inner.final_score()
    }

    /// Returns checked major-buffer memory accounting.
    #[must_use]
    pub const fn memory_estimate(&self) -> LinsiMemoryEstimate {
        self.inner.memory_estimate()
    }

    /// Consumes the audit record and returns only the final full-length MSA.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.inner.into_alignment()
    }
}

/// Deterministic reconstruction of ancestral MAFFT L-INS-i.
#[derive(Clone, Copy, Debug)]
pub struct MafftLinsi<A: MafftAlphabet> {
    parameters: LinsiParameters,
    thread_count: usize,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftLinsi<A> {
    /// Creates the immutable `MAFFT-2005-L-INS-i-v1` scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: LinsiParameters::paper_2005(),
            thread_count: 1,
            marker: PhantomData,
        }
    }

    /// Returns the stable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        LINSI_PROFILE_IDENTIFIER
    }

    /// Returns the immutable published and versioned parameters.
    #[must_use]
    pub const fn parameters(self) -> LinsiParameters {
        self.parameters
    }

    /// Returns the deterministic refinement worker count.
    #[must_use]
    pub const fn thread_count(self) -> usize {
        self.thread_count
    }

    /// Replaces the deterministic refinement worker count.
    ///
    /// Candidate partitions are evaluated concurrently and reduced in canonical
    /// tree-partition order, so this execution setting does not change the
    /// scientific result or its tie rules.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError::InvalidThreadCount`] when `count` is zero.
    pub const fn with_thread_count(mut self, count: usize) -> Result<Self, MafftError> {
        if count == 0 {
            return Err(MafftError::InvalidThreadCount { value: count });
        }
        self.thread_count = count;
        Ok(self)
    }

    /// Aligns named protein, DNA, or RNA records with ancestral L-INS-i.
    ///
    /// The local affine pair stage supplies distances and consistency evidence;
    /// the returned progressive and refined MSAs retain every source residue.
    ///
    /// # Errors
    ///
    /// Returns `MafftError` for invalid input, checked-size overflow, non-finite
    /// scientific calculations, or a violated alignment invariant.
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<LinsiAlignment<A>, MafftError> {
        let inner = align_insi(
            records,
            self.parameters.common(),
            InsiPairwiseStage::LocalAffine {
                gap_open: self.parameters.local_gap_open(),
                offset_penalty: self.parameters.local_offset_penalty(),
                gap_extension: self.parameters.local_gap_extension(),
            },
            self.thread_count,
        )?;
        Ok(LinsiAlignment { inner })
    }
}

impl<A: MafftAlphabet> Default for MafftLinsi<A> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::SequenceRecord;

    use super::{LINSI_PROFILE_IDENTIFIER, LinsiTermination, MafftLinsi};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn singleton_is_full_length_and_converges_without_candidates() {
        let result = MafftLinsi::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(
            MafftLinsi::<Dna>::new().profile_identifier(),
            LINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
        assert_eq!(result.termination(), LinsiTermination::Converged);
        assert_eq!(result.cycles().len(), 1);
        assert_eq!(result.cycles()[0].partition_count(), 0);
    }

    #[test]
    fn local_pair_ranges_exclude_flanks_but_final_rows_recover_full_sources() {
        let records = [dna("a", "TTACGTAA"), dna("b", "GGACGTCC")];
        let result = MafftLinsi::<Dna>::new().align(&records).unwrap();
        let pair = &result.pair_evidence()[0];
        assert_eq!(pair.first_aligned_range(), 2..6);
        assert_eq!(pair.second_aligned_range(), 2..6);
        assert_eq!(pair.fft_anchor_count(), 0);
        assert!(pair.library_entry_count() >= 4);
        for (row, expected) in result
            .alignment()
            .rows()
            .iter()
            .zip(["TTACGTAA", "GGACGTCC"])
        {
            assert_eq!(row.to_ungapped_ascii(), expected);
        }
    }

    #[test]
    fn protein_dna_and_rna_use_the_same_complete_pipeline() {
        let proteins = [
            SequenceRecord::<Protein>::new(
                "p1",
                ProteinSequence::from_ascii("XXMKTAYIAKXX").unwrap(),
            )
            .unwrap(),
            SequenceRecord::<Protein>::new(
                "p2",
                ProteinSequence::from_ascii("GGMKTIAKGG").unwrap(),
            )
            .unwrap(),
        ];
        let dna_records = [dna("d1", "TTACGTTAA"), dna("d2", "GGACTTCC")];
        let rna_records = [
            SequenceRecord::<Rna>::new("r1", RnaSequence::from_ascii("UUACGUUAA").unwrap())
                .unwrap(),
            SequenceRecord::<Rna>::new("r2", RnaSequence::from_ascii("GGACUUCC").unwrap()).unwrap(),
        ];
        assert_eq!(
            MafftLinsi::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftLinsi::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftLinsi::<Rna>::new()
                .align(&rna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
    }

    #[test]
    fn repeated_complete_runs_are_bitwise_deterministic() {
        let records = [
            dna("a", "TTACGTAA"),
            dna("b", "GGACGTCC"),
            dna("c", "CCACGTTT"),
        ];
        let first = MafftLinsi::<Dna>::new().align(&records).unwrap();
        let second = MafftLinsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(
            first.final_score().to_bits(),
            second.final_score().to_bits()
        );
        assert_eq!(
            first
                .pair_evidence()
                .iter()
                .map(|pair| (
                    pair.first_aligned_range(),
                    pair.second_aligned_range(),
                    pair.alignment_score().to_bits(),
                ))
                .collect::<Vec<_>>(),
            second
                .pair_evidence()
                .iter()
                .map(|pair| (
                    pair.first_aligned_range(),
                    pair.second_aligned_range(),
                    pair.alignment_score().to_bits(),
                ))
                .collect::<Vec<_>>()
        );
    }
}
