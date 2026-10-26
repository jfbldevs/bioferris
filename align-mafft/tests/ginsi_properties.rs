//! Property checks for G-INS-i determinism, invariants, and fixed objective monotonicity.

use align_core::{Dna, DnaSequence};
use align_mafft::{GinsiTermination, MafftGinsi};
use align_msa::SequenceRecord;
use proptest::prelude::*;

fn records(sequences: &[Vec<u8>]) -> Vec<SequenceRecord<Dna>> {
    sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| {
            SequenceRecord::new(
                format!("s{index}"),
                DnaSequence::from_ascii(sequence).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn arbitrary_dna_is_deterministic_rectangular_and_consumption_preserving(
        sequences in prop::collection::vec(
            prop::collection::vec(
                prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']),
                0..15,
            ),
            1..7,
        ),
    ) {
        let records = records(&sequences);
        let aligner = MafftGinsi::<Dna>::new();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.initial_alignment(), second.initial_alignment());
        prop_assert_eq!(first.distance_matrix(), second.distance_matrix());
        prop_assert_eq!(first.guide_tree(), second.guide_tree());
        prop_assert_eq!(first.progressive_sequence_weights(), second.progressive_sequence_weights());
        prop_assert_eq!(first.pair_evidence(), second.pair_evidence());
        prop_assert_eq!(first.progressive_merges(), second.progressive_merges());
        prop_assert_eq!(first.pair_weights(), second.pair_weights());
        prop_assert_eq!(first.final_score().to_bits(), second.final_score().to_bits());
        prop_assert_eq!(first.memory_estimate(), second.memory_estimate());
        prop_assert_eq!(first.termination(), second.termination());
        prop_assert_eq!(first.cycles().len(), second.cycles().len());

        for (record, row) in records.iter().zip(first.alignment().rows()) {
            prop_assert_eq!(row.identifier(), record.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
            prop_assert_eq!(row.len(), first.alignment().width());
        }
        for column in 0..first.alignment().width() {
            prop_assert!(first.alignment().column(column).unwrap().any(|cell| cell.is_some()));
        }
        for row in 0..records.len() {
            prop_assert_eq!(first.distance_matrix().get(row, row), Some(0.0));
            for column in 0..row {
                let distance = first.distance_matrix().get(row, column).unwrap();
                prop_assert!((0.0..=1.0).contains(&distance));
                prop_assert_eq!(
                    distance.to_bits(),
                    first.distance_matrix().get(column, row).unwrap().to_bits()
                );
            }
        }
        let expected_pairs = records.len() * records.len().saturating_sub(1) / 2;
        prop_assert_eq!(first.pair_evidence().len(), expected_pairs);
        prop_assert_eq!(first.library_statistics().pair_count(), expected_pairs);
        prop_assert_eq!(
            first
                .pair_evidence()
                .iter()
                .map(align_mafft::GinsiPairEvidence::library_entry_count)
                .sum::<usize>(),
            first.library_statistics().entry_count()
        );
        for pair in first.pair_evidence() {
            prop_assert_eq!(pair.generalized_skip_opening_count(), 0);
            prop_assert_eq!(pair.generalized_skipped_pair_count(), 0);
            prop_assert_eq!(pair.generalized_unpaired_residue_count(), 0);
        }
        for cycle in first.cycles() {
            let expected_partitions = records.len().saturating_mul(2).saturating_sub(3);
            prop_assert_eq!(cycle.partition_count(), expected_partitions);
            if cycle.accepted() {
                prop_assert!(cycle.score_after() > cycle.score_before());
            } else {
                prop_assert_eq!(cycle.score_after().to_bits(), cycle.score_before().to_bits());
            }
        }
        if first.termination() == GinsiTermination::Converged {
            prop_assert!(!first.cycles().last().unwrap().accepted());
        }
    }

    #[test]
    fn gotoh_weights_are_symmetric_positive_and_have_unit_pair_mean(
        sequences in prop::collection::vec(
            prop::collection::vec(
                prop::sample::select(vec![b'A', b'C', b'G', b'T']),
                0..16,
            ),
            1..8,
        ),
    ) {
        let records = records(&sequences);
        let result = MafftGinsi::<Dna>::new().align(&records).unwrap();
        let weights = result.pair_weights();
        for first in 0..records.len() {
            prop_assert_eq!(weights.get(first, first), Some(0.0));
            for second in 0..first {
                let value = weights.get(first, second).unwrap();
                prop_assert!(value.is_finite() && value > 0.0);
                prop_assert_eq!(value.to_bits(), weights.get(second, first).unwrap().to_bits());
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let expected = (records.len() * records.len().saturating_sub(1) / 2) as f64;
        prop_assert!((weights.sum_of_pairs() - expected).abs() < 1.0e-10);
        prop_assert!(result.final_score().is_finite());
        prop_assert!(
            result.memory_estimate().estimated_peak_bytes()
                >= result.memory_estimate().returned_final_alignment_bytes()
        );
    }
}
