//! Property checks for iterative refinement, pair weights, and determinism.

use align_core::{Dna, DnaSequence};
use align_mafft::{FftNsiTermination, MafftFftNs2, MafftFftNsi};
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
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn arbitrary_dna_is_deterministic_rectangular_and_consumption_preserving(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']), 0..16),
            1..7,
        ),
    ) {
        let records = records(&sequences);
        let aligner = MafftFftNsi::<Dna>::new();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        let parallel = aligner.with_thread_count(4).unwrap().align(&records).unwrap();
        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.alignment(), parallel.alignment());
        prop_assert_eq!(first.final_distance_matrix(), second.final_distance_matrix());
        prop_assert_eq!(first.final_distance_matrix(), parallel.final_distance_matrix());
        prop_assert_eq!(first.final_guide_tree(), second.final_guide_tree());
        prop_assert_eq!(first.final_guide_tree(), parallel.final_guide_tree());
        prop_assert_eq!(first.final_pair_weights(), second.final_pair_weights());
        prop_assert_eq!(first.final_pair_weights(), parallel.final_pair_weights());
        prop_assert_eq!(first.final_score().to_bits(), second.final_score().to_bits());
        prop_assert_eq!(first.final_score().to_bits(), parallel.final_score().to_bits());
        prop_assert_eq!(first.termination(), second.termination());
        prop_assert_eq!(first.termination(), parallel.termination());
        prop_assert_eq!(first.cycles().len(), second.cycles().len());
        prop_assert_eq!(first.cycles().len(), parallel.cycles().len());
        for (sequential_cycle, parallel_cycle) in first.cycles().iter().zip(parallel.cycles()) {
            prop_assert_eq!(
                sequential_cycle.best_partition(),
                parallel_cycle.best_partition()
            );
            prop_assert_eq!(
                sequential_cycle.score_before().to_bits(),
                parallel_cycle.score_before().to_bits()
            );
            prop_assert_eq!(
                sequential_cycle.best_candidate_score().map(f64::to_bits),
                parallel_cycle.best_candidate_score().map(f64::to_bits)
            );
            prop_assert_eq!(
                sequential_cycle.score_after().to_bits(),
                parallel_cycle.score_after().to_bits()
            );
            prop_assert_eq!(
                sequential_cycle.best_candidate_merge(),
                parallel_cycle.best_candidate_merge()
            );
        }
        for (record, row) in records.iter().zip(first.alignment().rows()) {
            prop_assert_eq!(row.identifier(), record.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
            prop_assert_eq!(row.len(), first.alignment().width());
        }
        for cycle in first.cycles() {
            let expected_partitions = records.len().saturating_mul(2).saturating_sub(3);
            prop_assert_eq!(cycle.partition_count(), expected_partitions);
            if cycle.accepted() {
                prop_assert!(cycle.score_after() > cycle.score_before());
                prop_assert!(cycle.accepted_partition().is_some());
            } else {
                prop_assert_eq!(cycle.score_after().to_bits(), cycle.score_before().to_bits());
            }
        }
        match first.termination() {
            FftNsiTermination::Converged => {
                prop_assert!(!first.cycles().last().unwrap().accepted());
            }
            FftNsiTermination::Oscillating => {
                let last = first.cycles().last().unwrap();
                prop_assert!(last.accepted());
                let repeated = last.score_after();
                prop_assert!(
                    first.cycles()[1..first.cycles().len() - 1]
                        .iter()
                        .any(|cycle| cycle.score_after().to_bits() == repeated.to_bits())
                );
            }
            FftNsiTermination::CycleLimit => {
                prop_assert_eq!(
                    first.accepted_cycle_count(),
                    aligner.parameters().maximum_cycles()
                );
            }
        }
    }

    #[test]
    fn gotoh_matrix_is_symmetric_positive_and_has_unit_pair_mean(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 0..18),
            1..8,
        ),
    ) {
        let records = records(&sequences);
        let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        let weights = result.final_pair_weights();
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
    }

    #[test]
    fn nested_initial_stage_equals_direct_fft_ns2(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 0..14),
            1..6,
        ),
    ) {
        let records = records(&sequences);
        let direct = MafftFftNs2::<Dna>::new().align(&records).unwrap();
        let iterative = MafftFftNsi::<Dna>::new().align(&records).unwrap();
        prop_assert_eq!(iterative.initial().alignment(), direct.alignment());
        prop_assert_eq!(
            iterative.initial().second_distance_matrix(),
            direct.second_distance_matrix()
        );
        prop_assert_eq!(iterative.initial().second_guide_tree(), direct.second_guide_tree());
        prop_assert_eq!(
            iterative.initial().second_sequence_weights(),
            direct.second_sequence_weights()
        );
        prop_assert_eq!(iterative.initial().second_merges(), direct.second_merges());
    }
}
