//! Property tests against an independent row dynamic program.

use align_core::DnaSequence;
use align_fast::{MyersPattern, verify_unit_edit_alignment};
use proptest::prelude::*;

fn dna(indices: &[u8]) -> DnaSequence {
    let bytes: Vec<u8> = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index) % 4])
        .collect();
    DnaSequence::from_ascii(bytes).expect("generator emits canonical DNA")
}

fn reference_global(pattern: &[u8], text: &[u8]) -> usize {
    let mut previous: Vec<usize> = (0..=pattern.len()).collect();
    let mut current = vec![0; pattern.len() + 1];
    for (column, text_symbol) in text.iter().copied().enumerate() {
        current[0] = column + 1;
        for row in 1..=pattern.len() {
            current[row] = (previous[row - 1] + usize::from(pattern[row - 1] != text_symbol))
                .min(previous[row] + 1)
                .min(current[row - 1] + 1);
        }
        core::mem::swap(&mut previous, &mut current);
    }
    previous[pattern.len()]
}

fn reference_search(pattern: &[u8], text: &[u8]) -> Vec<usize> {
    let mut previous: Vec<usize> = (0..=pattern.len()).collect();
    let mut current = vec![0; pattern.len() + 1];
    let mut scores = Vec::with_capacity(text.len());
    for text_symbol in text.iter().copied() {
        current[0] = 0;
        for row in 1..=pattern.len() {
            current[row] = (previous[row - 1] + usize::from(pattern[row - 1] != text_symbol))
                .min(previous[row] + 1)
                .min(current[row - 1] + 1);
        }
        core::mem::swap(&mut previous, &mut current);
        scores.push(previous[pattern.len()]);
    }
    scores
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(160))]

    #[test]
    fn global_distance_and_traceback_match_independent_oracle(
        pattern_indices in prop::collection::vec(0_u8..4, 0..140),
        text_indices in prop::collection::vec(0_u8..4, 0..140),
    ) {
        let pattern = dna(&pattern_indices);
        let text = dna(&text_indices);
        let pattern_bytes: Vec<u8> = pattern.iter().map(align_core::Symbol::as_ascii).collect();
        let text_bytes: Vec<u8> = text.iter().map(align_core::Symbol::as_ascii).collect();
        let expected = reference_global(&pattern_bytes, &text_bytes);
        let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
        let distance = myers.global_distance(&text).expect("distance succeeds");
        prop_assert_eq!(distance, expected);
        prop_assert_eq!(myers.global_distance(&pattern).expect("identity succeeds"), 0);
        prop_assert!(distance >= pattern.len().abs_diff(text.len()));
        prop_assert!(distance <= pattern.len().max(text.len()));

        let alignment = myers.align_global(&text).expect("traceback succeeds");
        let report = verify_unit_edit_alignment(&pattern, &text, &alignment)
            .expect("independent replay succeeds");
        prop_assert_eq!(report.distance(), expected);
        prop_assert_eq!(alignment.score(), -i64::try_from(expected).expect("generated distance fits"));
        prop_assert_eq!(myers.align_global(&text).expect("repeat succeeds"), alignment);

        let reverse = MyersPattern::new(&text)
            .expect("reverse preprocessing succeeds")
            .global_distance(&pattern)
            .expect("reverse distance succeeds");
        prop_assert_eq!(reverse, distance);
    }

    #[test]
    fn approximate_search_matches_independent_boundary_matrix(
        pattern_indices in prop::collection::vec(0_u8..4, 0..150),
        text_indices in prop::collection::vec(0_u8..4, 0..150),
        threshold in 0_usize..12,
    ) {
        let pattern = dna(&pattern_indices);
        let text = dna(&text_indices);
        let pattern_bytes: Vec<u8> = pattern.iter().map(align_core::Symbol::as_ascii).collect();
        let text_bytes: Vec<u8> = text.iter().map(align_core::Symbol::as_ascii).collect();
        let expected_scores = reference_search(&pattern_bytes, &text_bytes);
        let expected: Vec<(usize, usize)> = expected_scores
            .into_iter()
            .enumerate()
            .filter(|(_, distance)| *distance <= threshold)
            .map(|(column, distance)| (column + 1, distance))
            .collect();
        let observed: Vec<(usize, usize)> = MyersPattern::new(&pattern)
            .expect("preprocessing succeeds")
            .find_all(&text, threshold)
            .expect("search succeeds")
            .into_iter()
            .map(|hit| (hit.end(), hit.distance()))
            .collect();
        prop_assert_eq!(observed, expected);
    }

    #[test]
    fn global_distance_satisfies_triangle_inequality(
        first_indices in prop::collection::vec(0_u8..4, 0..45),
        second_indices in prop::collection::vec(0_u8..4, 0..45),
        third_indices in prop::collection::vec(0_u8..4, 0..45),
    ) {
        let first = dna(&first_indices);
        let second = dna(&second_indices);
        let third = dna(&third_indices);
        let first_pattern = MyersPattern::new(&first).expect("preprocessing succeeds");
        let second_pattern = MyersPattern::new(&second).expect("preprocessing succeeds");
        let direct = first_pattern.global_distance(&third).expect("distance succeeds");
        let via_second = first_pattern.global_distance(&second).expect("distance succeeds")
            + second_pattern.global_distance(&third).expect("distance succeeds");
        prop_assert!(direct <= via_second);
    }
}
