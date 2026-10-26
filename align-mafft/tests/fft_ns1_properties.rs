//! Property checks for deterministic, rectangular, consumption-preserving output.

use align_core::{Dna, DnaSequence};
use align_mafft::MafftFftNs1;
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
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 0..20),
            1..7,
        ),
    ) {
        let records = records(&sequences);
        let aligner = MafftFftNs1::<Dna>::new();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        prop_assert_eq!(first.alignment().width(), second.alignment().width());
        for ((record, first_row), second_row) in records
            .iter()
            .zip(first.alignment().rows())
            .zip(second.alignment().rows())
        {
            prop_assert_eq!(first_row.to_ungapped_ascii(), record.sequence().to_ascii_string());
            prop_assert_eq!(first_row.to_ascii_with_gaps(), second_row.to_ascii_with_gaps());
            prop_assert_eq!(first_row.len(), first.alignment().width());
        }
    }

    #[test]
    fn rough_distance_is_symmetric_and_bounded(
        sequences in prop::collection::vec(
            prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']), 0..24),
            1..8,
        ),
    ) {
        let records = records(&sequences);
        let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
        let distances = result.distance_matrix();
        for row in 0..records.len() {
            prop_assert_eq!(distances.get(row, row), Some(0.0));
            for column in 0..records.len() {
                let value = distances.get(row, column).unwrap();
                prop_assert!((0.0..=1.0).contains(&value));
                prop_assert_eq!(Some(value), distances.get(column, row));
            }
        }
    }
}
