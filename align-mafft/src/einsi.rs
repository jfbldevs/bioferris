use core::marker::PhantomData;

use align_msa::{DistanceMatrix, GuideTree, MultipleAlignment, SequenceRecord};

use crate::ginsi::{InsiPairwiseStage, align_insi};
use crate::{
    EinsiParameters, GinsiAlignment, GinsiCycle, GinsiLibraryStatistics, GinsiMemoryEstimate,
    GinsiMerge, GinsiPairEvidence, GinsiPairWeights, GinsiPartition, GinsiTermination,
    MafftAlphabet, MafftError,
};

/// Stable identifier for the ancestral scalar E-INS-i scientific profile.
pub const EINSI_PROFILE_IDENTIFIER: &str = "MAFFT-2008-E-INS-i-v1";

/// Auditable generalized-local pair evidence used by E-INS-i.
pub type EinsiPairEvidence = GinsiPairEvidence;

/// Aggregate statistics for the primary generalized-local consistency library.
pub type EinsiLibraryStatistics = GinsiLibraryStatistics;

/// Evidence for one complete-profile consistency-aware merge.
pub type EinsiMerge = GinsiMerge;

/// Evidence for one best-first fixed-tree refinement attempt.
pub type EinsiCycle = GinsiCycle;

/// Major-buffer memory accounting for one E-INS-i call.
pub type EinsiMemoryEstimate = GinsiMemoryEstimate;

/// Why E-INS-i stopped fixed-tree iterative refinement.
pub type EinsiTermination = GinsiTermination;

/// Canonical fixed-tree partition used during E-INS-i refinement.
pub type EinsiPartition = GinsiPartition;

/// Dense fixed Gotoh pair-weight matrix used by E-INS-i.
pub type EinsiPairWeights = GinsiPairWeights;

/// A completed `MAFFT-2008-E-INS-i-v1` alignment and its audit evidence.
#[derive(Clone, Debug)]
pub struct EinsiAlignment<A: MafftAlphabet> {
    inner: GinsiAlignment<A>,
}

impl<A: MafftAlphabet> EinsiAlignment<A> {
    /// Returns the final full-length refined alignment in input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        self.inner.alignment()
    }

    /// Returns the full-length E-INS-1 alignment before refinement.
    #[must_use]
    pub const fn initial_alignment(&self) -> &MultipleAlignment<A> {
        self.inner.initial_alignment()
    }

    /// Returns the all-pairs normalized generalized-local-score distance matrix.
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

    /// Returns generalized-local scores, ranges, distances, and skip evidence.
    #[must_use]
    pub const fn pair_evidence(&self) -> &[EinsiPairEvidence] {
        self.inner.pair_evidence()
    }

    /// Returns aggregate primary generalized-local-library statistics.
    #[must_use]
    pub const fn library_statistics(&self) -> EinsiLibraryStatistics {
        self.inner.library_statistics()
    }

    /// Returns complete-profile progressive merges in tree postorder.
    #[must_use]
    pub const fn progressive_merges(&self) -> &[EinsiMerge] {
        self.inner.progressive_merges()
    }

    /// Returns accepted attempts and the terminal convergence attempt.
    #[must_use]
    pub const fn cycles(&self) -> &[EinsiCycle] {
        self.inner.cycles()
    }

    /// Returns the number of strictly improving refinement attempts.
    #[must_use]
    pub fn accepted_cycle_count(&self) -> usize {
        self.inner.accepted_cycle_count()
    }

    /// Returns whether convergence or the cycle limit stopped refinement.
    #[must_use]
    pub const fn termination(&self) -> EinsiTermination {
        self.inner.termination()
    }

    /// Returns the fixed dense Gotoh pair weights.
    #[must_use]
    pub const fn pair_weights(&self) -> &EinsiPairWeights {
        self.inner.pair_weights()
    }

    /// Returns the independently verified final complete-MSA WSP+I objective.
    #[must_use]
    pub const fn final_score(&self) -> f64 {
        self.inner.final_score()
    }

    /// Returns checked major-buffer memory accounting.
    #[must_use]
    pub const fn memory_estimate(&self) -> EinsiMemoryEstimate {
        self.inner.memory_estimate()
    }

    /// Consumes the audit record and returns only the final full-length MSA.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.inner.into_alignment()
    }
}

/// Deterministic scalar reconstruction of ancestral MAFFT E-INS-i.
#[derive(Clone, Copy, Debug)]
pub struct MafftEinsi<A: MafftAlphabet> {
    parameters: EinsiParameters,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> MafftEinsi<A> {
    /// Creates the immutable `MAFFT-2008-E-INS-i-v1` scientific profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            parameters: EinsiParameters::paper_2008(),
            marker: PhantomData,
        }
    }

    /// Returns the stable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        EINSI_PROFILE_IDENTIFIER
    }

    /// Returns the immutable published and versioned parameters.
    #[must_use]
    pub const fn parameters(self) -> EinsiParameters {
        self.parameters
    }

    /// Aligns named protein, DNA, or RNA records with ancestral E-INS-i.
    ///
    /// The generalized-local pair stage connects conserved blocks while
    /// excluding skipped residues from the primary homology library. The
    /// progressive and refined MSAs retain every source residue.
    ///
    /// # Errors
    ///
    /// Returns `MafftError` for invalid input, checked-size overflow, non-finite
    /// scientific calculations, or a violated alignment invariant.
    pub fn align(&self, records: &[SequenceRecord<A>]) -> Result<EinsiAlignment<A>, MafftError> {
        let inner = align_insi(
            records,
            self.parameters.common(),
            InsiPairwiseStage::GeneralizedAffine {
                ordinary_gap_open: self.parameters.local_gap_open(),
                offset_penalty: self.parameters.local_offset_penalty(),
                ordinary_gap_extension: self.parameters.local_gap_extension(),
                skip_gap_open: self.parameters.skip_gap_open(),
                paired_skip_extension: self.parameters.paired_skip_extension(),
            },
            1,
        )?;
        Ok(EinsiAlignment { inner })
    }
}

impl<A: MafftAlphabet> Default for MafftEinsi<A> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
    use align_msa::SequenceRecord;

    use super::{EINSI_PROFILE_IDENTIFIER, EinsiTermination, MafftEinsi};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn singleton_is_full_length_and_converges_without_candidates() {
        let result = MafftEinsi::<Dna>::new()
            .align(&[dna("only", "ACGT")])
            .unwrap();
        assert_eq!(
            MafftEinsi::<Dna>::new().profile_identifier(),
            EINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
        assert_eq!(result.termination(), EinsiTermination::Converged);
        assert_eq!(result.cycles().len(), 1);
        assert_eq!(result.cycles()[0].partition_count(), 0);
    }

    #[test]
    fn generalized_pair_retains_full_sources_and_reports_skips() {
        let records = [
            dna("a", "TTACGTACGTACGTAAAAAAAAAAAAACGTACGTACGTAA"),
            dna("b", "GGACGTACGTACGTCCCCCCCCCCCCCACGTACGTACGTCC"),
        ];
        let result = MafftEinsi::<Dna>::new().align(&records).unwrap();
        let pair = &result.pair_evidence()[0];
        assert!(pair.generalized_skip_opening_count() >= 1);
        assert!(pair.generalized_skipped_pair_count() >= 1);
        assert_eq!(pair.fft_anchor_count(), 0);
        assert!(pair.library_entry_count() >= 20);
        for (row, expected) in result.alignment().rows().iter().zip(
            records
                .iter()
                .map(|record| record.sequence().to_ascii_string()),
        ) {
            assert_eq!(row.to_ungapped_ascii(), expected);
        }
    }

    #[test]
    fn protein_dna_and_rna_use_the_same_complete_pipeline() {
        let proteins = [
            SequenceRecord::<Protein>::new(
                "p1",
                ProteinSequence::from_ascii("MKTAYIAKAAAAAAAAAAAAGHPEPTIDE").unwrap(),
            )
            .unwrap(),
            SequenceRecord::<Protein>::new(
                "p2",
                ProteinSequence::from_ascii("MKTAYIAKCCCCCCCCCCCCGHPEPTIDE").unwrap(),
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
            MafftEinsi::<Protein>::new()
                .align(&proteins)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftEinsi::<Dna>::new()
                .align(&dna_records)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            MafftEinsi::<Rna>::new()
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
            dna("a", "ACGTACGTAAAAAAAAACGTACGT"),
            dna("b", "ACGTACGTCCCCCCCCCACGTACGT"),
            dna("c", "ACGTACGTGGGGGGGGGACGTACGT"),
        ];
        let first = MafftEinsi::<Dna>::new().align(&records).unwrap();
        let second = MafftEinsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(
            first.final_score().to_bits(),
            second.final_score().to_bits()
        );
        assert_eq!(first.pair_evidence(), second.pair_evidence());
    }
}
