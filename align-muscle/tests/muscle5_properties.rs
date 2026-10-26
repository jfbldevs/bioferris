//! Property tests for public MUSCLE 5 determinism and structural guarantees.

use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;
use align_muscle::{Muscle5, Muscle5Parameters};
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

fn assert_public_invariants(
    result: &align_muscle::Muscle5Alignment<Dna>,
    raw: &[Vec<u8>],
) -> Result<(), TestCaseError> {
    let alignment = result.alignment();
    prop_assert_eq!(alignment.sequence_count(), raw.len());
    for (index, row) in alignment.rows().iter().enumerate() {
        prop_assert_eq!(row.identifier(), format!("seq{index}"));
        let ungapped = row.to_ungapped_ascii();
        prop_assert_eq!(ungapped.as_bytes(), raw[index].as_slice());
        prop_assert_eq!(row.len(), alignment.width());
    }
    for column in 0..alignment.width() {
        prop_assert!(alignment.column(column).unwrap().any(|cell| cell.is_some()));
    }
    for attempt in result.refinement_attempts() {
        prop_assert!(attempt.score_after() + 1.0e-9 >= attempt.score_before());
        if attempt.accepted() {
            prop_assert!(attempt.score_after() > attempt.score_before());
        }
    }
    prop_assert!(result.evaluated_objective_pair_count() <= result.possible_objective_pair_count());
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn align_preserves_rows_distances_and_bit_determinism(
        raw in prop::collection::vec(
            prop::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')],
                0..9,
            ),
            1..6,
        ),
        seed in 0_u64..16,
    ) {
        let records = records(&raw);
        let parameters = Muscle5Parameters::published_align()
            .with_iterations(1, 4)
            .with_perturbation_seed(seed);
        let aligner = Muscle5::<Dna>::with_parameters(parameters).unwrap();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_public_invariants(&first, &raw)?;
        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.raw_posterior_checksum(), second.raw_posterior_checksum());
        prop_assert_eq!(first.consistent_posterior_checksum(), second.consistent_posterior_checksum());
        prop_assert_eq!(first.expected_error_distances(), second.expected_error_distances());
        prop_assert_eq!(first.progressive_steps(), second.progressive_steps());
        prop_assert_eq!(first.refinement_attempts(), second.refinement_attempts());
        prop_assert_eq!(first.memory_estimate(), second.memory_estimate());
        for row in 0..raw.len() {
            for column in 0..raw.len() {
                let forward = first.expected_error_distances().get(row, column).unwrap();
                let reverse = first.expected_error_distances().get(column, row).unwrap();
                prop_assert!((0.0..=1.0).contains(&forward));
                prop_assert!((forward - reverse).abs() <= 1.0e-12);
            }
        }
    }

    #[test]
    fn super5_preserves_rows_and_characteristic_evidence(
        raw in prop::collection::vec(
            prop::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                0..8,
            ),
            1..8,
        ),
        seed in 0_u64..8,
    ) {
        let records = records(&raw);
        let parameters = Muscle5Parameters::published_super5()
            .with_iterations(0, 2)
            .with_perturbation_seed(seed)
            .with_super5_limits(0.0, 0.30, 0.10, 3, 4, 3)
            .unwrap();
        let aligner = Muscle5::<Dna>::with_parameters(parameters).unwrap();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_public_invariants(&first, &raw)?;
        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.super5_evidence(), second.super5_evidence());
        let evidence = first.super5_evidence().unwrap();
        prop_assert!(evidence.clusters().iter().all(|cluster| cluster.representative_members().len() <= 3));
        prop_assert_eq!(first.guide_tree().leaf_count(), evidence.clusters().len());
    }
}
