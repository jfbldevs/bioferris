//! Property tests for shared alignment invariants and metrics.

use align_core::{Dna, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, compare_alignments};
use proptest::prelude::*;

fn make_alignment(raw: &[Vec<u8>], slots: &[Vec<usize>]) -> MultipleAlignment<Dna> {
    let width = raw
        .iter()
        .zip(slots)
        .map(|(sequence, positions)| sequence.len() + positions.len())
        .max()
        .unwrap_or(0);
    let rows = raw
        .iter()
        .enumerate()
        .map(|(row_index, sequence)| {
            let mut cells = sequence
                .iter()
                .map(|byte| Some(Symbol::<Dna>::from_ascii(*byte).unwrap()))
                .collect::<Vec<_>>();
            for position in &slots[row_index] {
                let position = (*position).min(cells.len());
                cells.insert(position, None);
            }
            cells.resize(width, None);
            AlignedRow::new(format!("s{row_index}"), cells).unwrap()
        })
        .collect::<Vec<_>>();

    // Remove columns that became all gaps after right-padding.
    let keep = (0..width)
        .filter(|column| rows.iter().any(|row| row.symbols()[*column].is_some()))
        .collect::<Vec<_>>();
    MultipleAlignment::try_new(
        rows.into_iter()
            .map(|row| {
                AlignedRow::new(
                    row.identifier(),
                    keep.iter().map(|column| row.symbols()[*column]).collect(),
                )
                .unwrap()
            })
            .collect(),
    )
    .unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn self_comparison_is_exact_and_bounded(
        raw in prop::collection::vec(
            prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 0..16),
            1..7,
        ),
        gap_seed in prop::collection::vec(prop::collection::vec(0_usize..20, 0..5), 1..7),
    ) {
        let slots = (0..raw.len())
            .map(|index| gap_seed.get(index).cloned().unwrap_or_default())
            .collect::<Vec<_>>();
        let alignment = make_alignment(&raw, &slots);
        let metrics = compare_alignments(&alignment, &alignment).unwrap();
        prop_assert!((metrics.sp_score() - 1.0).abs() < f64::EPSILON);
        prop_assert!((metrics.tc_score() - 1.0).abs() < f64::EPSILON);
        prop_assert!(metrics.correct_pairs() <= metrics.reference_pairs());
        prop_assert!(metrics.correct_columns() <= metrics.reference_columns());
    }
}
