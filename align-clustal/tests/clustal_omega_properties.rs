//! Property tests for public Clustal Omega invariants and determinism.

use align_clustal::{ClustalOmega, ClustalOmegaParameters, OmegaEngineMode, OmegaGuideTreeMode};
use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;
use proptest::prelude::*;

fn records(raw: &[Vec<u8>]) -> Vec<SequenceRecord<Dna>> {
    raw.iter()
        .enumerate()
        .map(|(index, sequence)| {
            SequenceRecord::new(
                format!("seq{index}"),
                DnaSequence::from_ascii(sequence).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

fn dna_sequence() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
        0..11,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn full_viterbi_pipeline_preserves_every_public_structure(
        raw in prop::collection::vec(dna_sequence(), 1..6),
    ) {
        let records = records(&raw);
        let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
            .with_guide_tree_mode(OmegaGuideTreeMode::Full)
            .with_engine_mode(OmegaEngineMode::Viterbi);
        let result = ClustalOmega::<Dna>::with_parameters(parameters)
            .unwrap()
            .align(&records)
            .unwrap();
        let alignment = result.alignment();

        prop_assert_eq!(alignment.sequence_count(), raw.len());
        prop_assert_eq!(result.guide_tree().leaf_count(), raw.len());
        prop_assert_eq!(result.progressive_steps().len(), raw.len().saturating_sub(1));
        for (index, row) in alignment.rows().iter().enumerate() {
            prop_assert_eq!(row.identifier(), format!("seq{index}"));
            let ungapped = row.to_ungapped_ascii();
            prop_assert_eq!(ungapped.as_bytes(), raw[index].as_slice());
            prop_assert_eq!(row.len(), alignment.width());
        }
        for column in 0..alignment.width() {
            prop_assert!(alignment.column(column).unwrap().any(|cell| cell.is_some()));
        }

        let distances = result.guide_evidence().distances().unwrap();
        prop_assert_eq!(distances.size(), raw.len());
        for first in 0..raw.len() {
            prop_assert_eq!(distances.get(first, first), Some(0.0));
            for second in 0..raw.len() {
                let forward = distances.get(first, second).unwrap();
                let reverse = distances.get(second, first).unwrap();
                prop_assert!((0.0..=1.0).contains(&forward));
                prop_assert_eq!(forward.to_bits(), reverse.to_bits());
            }
        }
        let progressive_steps_are_valid = result.progressive_steps().iter().all(|step| {
            step.score().is_finite()
                && step.result_width() <= alignment.width()
                && step.left_core().start <= step.left_core().end
                && step.right_core().start <= step.right_core().end
        });
        prop_assert!(progressive_steps_are_valid);
        let memory = result.memory_estimate();
        prop_assert!(memory.estimated_peak_bytes() >= memory.returned_alignment_bytes());
        prop_assert!(memory.estimated_peak_bytes() >= memory.retained_evidence_bytes());
    }

    #[test]
    fn repeated_mac_calls_are_bit_deterministic(
        raw in prop::collection::vec(dna_sequence(), 1..5),
    ) {
        let records = records(&raw);
        let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
            .with_guide_tree_mode(OmegaGuideTreeMode::Full)
            .with_engine_mode(OmegaEngineMode::Mac);
        let aligner = ClustalOmega::<Dna>::with_parameters(parameters).unwrap();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();

        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.guide_tree(), second.guide_tree());
        prop_assert_eq!(first.progressive_steps(), second.progressive_steps());
        prop_assert_eq!(first.memory_estimate(), second.memory_estimate());
        prop_assert_eq!(
            first.guide_evidence().distances().unwrap().as_flat_slice(),
            second.guide_evidence().distances().unwrap().as_flat_slice(),
        );
    }

    #[test]
    fn one_tuple_distance_is_symmetric_and_bounded(
        first in dna_sequence(),
        second in dna_sequence(),
    ) {
        let first = DnaSequence::from_ascii(&first).unwrap();
        let second = DnaSequence::from_ascii(&second).unwrap();
        let aligner = ClustalOmega::<Dna>::published();
        let forward = aligner.pair_distance(&first, &second).unwrap();
        let reverse = aligner.pair_distance(&second, &first).unwrap();
        prop_assert!((0.0..=1.0).contains(&forward));
        prop_assert_eq!(forward.to_bits(), reverse.to_bits());
        if first == second {
            prop_assert_eq!(forward.to_bits(), 0.0_f64.to_bits());
        }
    }
}
