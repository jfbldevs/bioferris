//! Property checks for the two-pass pipeline and its independently replayed distances.

use align_core::{Dna, DnaSequence};
use align_mafft::{FFT_NS2_POISSON_IDENTITY_FLOOR, MafftFftNs1, MafftFftNs2};
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
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn arbitrary_dna_is_rectangular_deterministic_and_consumption_preserving(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']), 0..20),
            1..7,
        ),
    ) {
        let records = records(&sequences);
        let aligner = MafftFftNs2::<Dna>::new();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        prop_assert_eq!(first.alignment(), second.alignment());
        prop_assert_eq!(first.second_distance_matrix(), second.second_distance_matrix());
        prop_assert_eq!(first.second_guide_tree(), second.second_guide_tree());
        prop_assert_eq!(first.second_sequence_weights(), second.second_sequence_weights());
        prop_assert_eq!(first.second_merges(), second.second_merges());
        for (record, row) in records.iter().zip(first.alignment().rows()) {
            prop_assert_eq!(row.identifier(), record.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
            prop_assert_eq!(row.len(), first.alignment().width());
        }
    }

    #[test]
    fn exposed_poisson_matrix_matches_an_independent_column_oracle(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']), 0..24),
            1..8,
        ),
    ) {
        let records = records(&sequences);
        let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
        let alignment = result.first_pass().alignment();
        let matrix = result.second_distance_matrix();
        for row in 0..records.len() {
            prop_assert_eq!(matrix.get(row, row), Some(0.0));
            for column in 0..records.len() {
                let left = alignment.rows()[row].symbols();
                let right = alignment.rows()[column].symbols();
                let mut comparable = 0usize;
                let mut identical = 0usize;
                for (left_cell, right_cell) in left.iter().zip(right) {
                    if let (Some(left_symbol), Some(right_symbol)) = (left_cell, right_cell) {
                        comparable += 1;
                        identical += usize::from(left_symbol == right_symbol);
                    }
                }
                let left_empty = left.iter().all(Option::is_none);
                let right_empty = right.iter().all(Option::is_none);
                let expected = if row == column || (comparable == 0 && left_empty && right_empty) {
                    0.0
                } else if comparable == 0 {
                    -FFT_NS2_POISSON_IDENTITY_FLOOR.ln()
                } else {
                    #[allow(clippy::cast_precision_loss)]
                    let unchanged = identical as f64 / comparable as f64;
                    -unchanged.max(FFT_NS2_POISSON_IDENTITY_FLOOR).ln()
                };
                let observed = matrix.get(row, column).unwrap();
                prop_assert!(observed.is_finite() && observed >= 0.0);
                prop_assert!((observed - expected).abs() < 1.0e-12);
                prop_assert_eq!(matrix.get(column, row), Some(observed));
                prop_assert_eq!(matrix.comparable_columns(row, column), Some(comparable));
                prop_assert_eq!(matrix.identical_columns(row, column), Some(identical));
            }
        }
    }

    #[test]
    fn nested_first_pass_equals_a_direct_fft_ns1_call(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 0..18),
            1..6,
        ),
    ) {
        let records = records(&sequences);
        let direct = MafftFftNs1::<Dna>::new().align(&records).unwrap();
        let nested = MafftFftNs2::<Dna>::new().align(&records).unwrap();
        prop_assert_eq!(nested.first_pass().alignment(), direct.alignment());
        prop_assert_eq!(nested.first_pass().distance_matrix(), direct.distance_matrix());
        prop_assert_eq!(nested.first_pass().guide_tree(), direct.guide_tree());
        prop_assert_eq!(nested.first_pass().sequence_weights(), direct.sequence_weights());
        prop_assert_eq!(nested.first_pass().merges(), direct.merges());
    }
}
