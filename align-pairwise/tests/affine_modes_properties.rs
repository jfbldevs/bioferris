//! Property-based tests for affine local and explicit four-end alignment.

use align_core::{AffineGap, DnaSequence, LinearGap, MatchMismatch};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, FreeEndGaps, Gotoh, SmithWaterman,
    verify_affine_end_gap_alignment, verify_local_affine_alignment,
};
use proptest::prelude::*;

fn dna_strategy(max_length: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')],
        0..=max_length,
    )
}

fn scoring() -> MatchMismatch {
    MatchMismatch::new(3, -4)
}

fn gap() -> AffineGap {
    AffineGap::new(-5, -1).expect("valid affine gap")
}

fn mirrored(free_ends: FreeEndGaps) -> FreeEndGaps {
    FreeEndGaps::new(
        free_ends.second_start(),
        free_ends.second_end(),
        free_ends.first_start(),
        free_ends.first_end(),
    )
}

proptest! {
    #![proptest_config(ProptestConfig::default())]

    #[test]
    fn every_four_end_configuration_replays_and_is_deterministic(
        first_raw in dna_strategy(12),
        second_raw in dna_strategy(12),
        first_start in any::<bool>(),
        first_end in any::<bool>(),
        second_start in any::<bool>(),
        second_end in any::<bool>(),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let ends = FreeEndGaps::new(first_start, first_end, second_start, second_end);
        let aligner = AffineEndGapAligner::new(scoring(), gap(), ends);
        let alignment = aligner.align(&first, &second).expect("affine alignment succeeds");
        let repeated = aligner.align(&first, &second).expect("repeated alignment succeeds");

        prop_assert_eq!(&alignment, &repeated);
        let report = verify_affine_end_gap_alignment(
            &first,
            &second,
            &alignment,
            aligner.substitution(),
            aligner.gap(),
            aligner.free_ends(),
        ).expect("traceback and score replay");
        prop_assert_eq!(report.recomputed_score(), alignment.score());
        prop_assert!(alignment.score() <= 3 * i64::try_from(first.len().min(second.len())).expect("small length"));
    }

    #[test]
    fn swapping_sequences_and_end_permissions_preserves_score(
        first_raw in dna_strategy(12),
        second_raw in dna_strategy(12),
        first_start in any::<bool>(),
        first_end in any::<bool>(),
        second_start in any::<bool>(),
        second_end in any::<bool>(),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let ends = FreeEndGaps::new(first_start, first_end, second_start, second_end);
        let forward = AffineEndGapAligner::new(scoring(), gap(), ends)
            .align(&first, &second).expect("forward alignment succeeds");
        let reverse = AffineEndGapAligner::new(scoring(), gap(), mirrored(ends))
            .align(&second, &first).expect("reverse alignment succeeds");

        prop_assert_eq!(forward.score(), reverse.score());
    }

    #[test]
    fn global_preset_remains_identical_to_gotoh(
        first_raw in dna_strategy(12),
        second_raw in dna_strategy(12),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let modes = AffineEndGapAligner::new(scoring(), gap(), FreeEndGaps::global())
            .align(&first, &second).expect("affine global mode succeeds");
        let gotoh = Gotoh::new(scoring(), gap())
            .align(&first, &second).expect("Gotoh succeeds");

        prop_assert_eq!(modes, gotoh);
    }

    #[test]
    fn local_results_replay_are_deterministic_and_are_symmetric_in_score(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = AffineLocalAligner::new(scoring(), gap());
        let forward = aligner.align(&first, &second).expect("local alignment succeeds");
        let repeated = aligner.align(&first, &second).expect("repeated alignment succeeds");
        let reverse = aligner.align(&second, &first).expect("reverse alignment succeeds");

        prop_assert_eq!(&forward, &repeated);
        prop_assert_eq!(forward.score(), reverse.score());
        prop_assert!(forward.score() >= 0);
        prop_assert!(forward.score() <= 3 * i64::try_from(first.len().min(second.len())).expect("small length"));
        let report = verify_local_affine_alignment(
            &first,
            &second,
            &forward,
            aligner.substitution(),
            aligner.gap(),
        ).expect("local traceback and score replay");
        prop_assert_eq!(report.recomputed_score(), forward.score());
    }

    #[test]
    fn equal_open_and_extension_local_score_reduces_to_smith_waterman(
        first_raw in dna_strategy(12),
        second_raw in dna_strategy(12),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let affine_gap = AffineGap::new(-2, -2).expect("linear affine specialization");
        let linear_gap = LinearGap::new(-2).expect("valid linear gap");
        let affine = AffineLocalAligner::new(scoring(), affine_gap)
            .align(&first, &second).expect("affine local succeeds");
        let linear = SmithWaterman::new(scoring(), linear_gap)
            .align(&first, &second).expect("linear local succeeds");

        prop_assert_eq!(affine.score(), linear.score());
        prop_assert_eq!(affine.first_range().end, linear.first_range().end);
        prop_assert_eq!(affine.second_range().end, linear.second_range().end);
    }
}
