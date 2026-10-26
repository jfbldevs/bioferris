//! Property tests for affine X-drop traversal, traceback, symmetry, and forced seeds.

mod support;

use align_core::{AffineGap, DnaSequence, MatchMismatch, SubstitutionScorer};
use align_fast::{XDrop, XDropDirection, XDropSeed, XDropTermination};
use proptest::prelude::*;
use support::xdrop::{OracleTermination, oracle, replay};

fn dna(indices: &[u8]) -> DnaSequence {
    let bytes: Vec<u8> = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index) % 4])
        .collect();
    DnaSequence::from_ascii(bytes).expect("generator emits canonical DNA")
}

fn direction(index: u8) -> XDropDirection {
    match index % 3 {
        0 => XDropDirection::Upstream,
        1 => XDropDirection::Downstream,
        _ => XDropDirection::Both,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn rolling_frontiers_match_dense_oracle_and_replay(
        first_indices in prop::collection::vec(0_u8..4, 0..28),
        second_indices in prop::collection::vec(0_u8..4, 0..28),
        match_score in -2_i64..7,
        mismatch_score in -7_i64..3,
        extension_score in -4_i64..1,
        open_extra in 0_i64..8,
        threshold in 0_i64..24,
    ) {
        let first = dna(&first_indices);
        let second = dna(&second_indices);
        let scoring = MatchMismatch::new(match_score, mismatch_score);
        let open_score = extension_score - open_extra;
        let gap = AffineGap::new(open_score, extension_score).expect("generated gap is valid");
        let aligner = XDrop::new(scoring, gap, threshold).expect("generated threshold is valid");
        let expected = oracle(&first, &second, &scoring, gap, threshold);
        let observed = aligner.score(&first, &second).expect("score succeeds");

        prop_assert_eq!(observed.score(), expected.score);
        prop_assert_eq!(observed.first_end(), expected.first_end);
        prop_assert_eq!(observed.second_end(), expected.second_end);
        prop_assert_eq!(observed.evaluated_cells(), expected.evaluated_cells);
        prop_assert_eq!(observed.retained_cells(), expected.retained_cells);
        prop_assert_eq!(observed.antidiagonals(), expected.antidiagonals);
        prop_assert_eq!(
            observed.maximum_frontier_width(),
            expected.maximum_frontier_width
        );
        prop_assert_eq!(
            observed.termination(),
            match expected.termination {
                OracleTermination::MatrixExhausted => XDropTermination::MatrixExhausted,
                OracleTermination::FrontiersExhausted => XDropTermination::FrontiersExhausted,
            }
        );

        let alignment = aligner.align(&first, &second).expect("traceback succeeds");
        prop_assert_eq!(alignment.summary(), observed);
        prop_assert_eq!(
            replay(&first, &second, alignment.alignment(), &scoring, gap),
            observed.score()
        );
        prop_assert_eq!(
            aligner.align(&first, &second).expect("repeat succeeds"),
            alignment
        );
    }

    #[test]
    fn symmetric_scoring_produces_symmetric_best_score(
        first_indices in prop::collection::vec(0_u8..4, 0..40),
        second_indices in prop::collection::vec(0_u8..4, 0..40),
        match_score in 0_i64..7,
        mismatch_score in -7_i64..2,
        extension_score in -3_i64..1,
        open_extra in 0_i64..7,
        threshold in 0_i64..30,
    ) {
        let first = dna(&first_indices);
        let second = dna(&second_indices);
        let scoring = MatchMismatch::new(match_score, mismatch_score);
        let gap = AffineGap::new(extension_score - open_extra, extension_score)
            .expect("generated gap is valid");
        let aligner = XDrop::new(scoring, gap, threshold).expect("valid configuration");
        let forward = aligner.score(&first, &second).expect("forward score succeeds");
        let transposed = aligner.score(&second, &first).expect("transposed score succeeds");
        prop_assert_eq!(forward.score(), transposed.score());
        prop_assert_eq!(forward.evaluated_cells(), transposed.evaluated_cells());
        prop_assert_eq!(forward.retained_cells(), transposed.retained_cells());
        prop_assert_eq!(forward.antidiagonals(), transposed.antidiagonals());
        prop_assert_eq!(
            forward.maximum_frontier_width(),
            transposed.maximum_frontier_width()
        );
    }

    #[test]
    fn forced_seed_score_is_the_sum_of_independent_sides_and_one_seed(
        first_indices in prop::collection::vec(0_u8..4, 1..24),
        second_indices in prop::collection::vec(0_u8..4, 1..24),
        first_seed_index in any::<usize>(),
        second_seed_index in any::<usize>(),
        direction_index in 0_u8..3,
        threshold in 0_i64..24,
    ) {
        let first = dna(&first_indices);
        let second = dna(&second_indices);
        let first_seed_index = first_seed_index % first.len();
        let second_seed_index = second_seed_index % second.len();
        let seed = XDropSeed::new(first_seed_index, second_seed_index);
        let direction = direction(direction_index);
        let scoring = MatchMismatch::new(3, -4);
        let gap = AffineGap::new(-6, -1).expect("valid gap");
        let aligner = XDrop::new(scoring, gap, threshold).expect("valid configuration");
        let summary = aligner
            .score_from_seed(&first, &second, seed, direction)
            .expect("seeded score succeeds");
        let seed_score = scoring.score(
            first.get(first_seed_index).expect("seed is in range"),
            second.get(second_seed_index).expect("seed is in range"),
        );
        let expected_score = seed_score
            + summary.upstream().map_or(0, align_fast::XDropScore::score)
            + summary.downstream().map_or(0, align_fast::XDropScore::score);
        prop_assert_eq!(summary.score(), expected_score);
        prop_assert!(summary.first_range().contains(&first_seed_index));
        prop_assert!(summary.second_range().contains(&second_seed_index));
        prop_assert!(summary.first_range().end <= first.len());
        prop_assert!(summary.second_range().end <= second.len());

        let alignment = aligner
            .align_from_seed(&first, &second, seed, direction)
            .expect("seeded traceback succeeds");
        prop_assert_eq!(alignment.summary(), summary);
        prop_assert_eq!(
            replay(&first, &second, alignment.alignment(), &scoring, gap),
            summary.score()
        );
        prop_assert_eq!(
            aligner
                .align_from_seed(&first, &second, seed, direction)
                .expect("repeat succeeds"),
            alignment
        );
    }

    #[test]
    fn zero_drop_identity_reaches_the_full_diagonal(
        indices in prop::collection::vec(0_u8..4, 0..100),
        match_score in 1_i64..8,
        mismatch_penalty in 1_i64..8,
        extension_penalty in 1_i64..5,
        open_extra in 0_i64..8,
    ) {
        let sequence = dna(&indices);
        let scoring = MatchMismatch::new(match_score, -mismatch_penalty);
        let extension_score = -extension_penalty;
        let gap = AffineGap::new(extension_score - open_extra, extension_score)
            .expect("generated gap is valid");
        let result = XDrop::new(scoring, gap, 0)
            .expect("valid configuration")
            .align(&sequence, &sequence)
            .expect("identity alignment succeeds");
        let expected_score = match_score * i64::try_from(sequence.len()).expect("length fits");
        prop_assert_eq!(result.score(), expected_score);
        prop_assert_eq!(result.first_end(), sequence.len());
        prop_assert_eq!(result.second_end(), sequence.len());
        prop_assert!(
            result
                .alignment()
                .operations()
                .iter()
                .all(|operation| matches!(operation, align_core::AlignmentOperation::Match))
        );
    }
}
