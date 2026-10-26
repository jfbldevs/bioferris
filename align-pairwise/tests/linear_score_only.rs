//! Score-only linear kernels return the traceback kernels' scores and errors.

use align_core::{Dna, DnaSequence, LinearGap, MatchMismatch, Score, SubstitutionScorer, Symbol};
use align_pairwise::{NeedlemanWunsch, SmithWaterman};
use proptest::prelude::*;

/// A literal scorer without declared bounds, which selects the checked kernels.
#[derive(Clone, Copy)]
struct Unbounded(MatchMismatch);

impl SubstitutionScorer<Dna> for Unbounded {
    fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
        self.0.score(first, second)
    }
}

fn dna() -> impl Strategy<Value = DnaSequence> {
    proptest::collection::vec(prop::sample::select(b"ACGT".to_vec()), 0..40)
        .prop_map(|bytes| DnaSequence::from_ascii(&bytes).unwrap())
}

fn gap(score: Score) -> LinearGap {
    LinearGap::new(score).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn global_scores_equal_traceback_scores(
        first in dna(),
        second in dna(),
        matched in 0_i64..=6,
        mismatched in -6_i64..=0,
        gap_score in -6_i64..=0,
    ) {
        let scoring = MatchMismatch::new(matched, mismatched);
        let bounded = NeedlemanWunsch::new(scoring, gap(gap_score));
        prop_assert_eq!(
            bounded.score(&first, &second).unwrap(),
            bounded.align(&first, &second).unwrap().score()
        );
        let checked = NeedlemanWunsch::new(Unbounded(scoring), gap(gap_score));
        prop_assert_eq!(
            checked.score(&first, &second).unwrap(),
            checked.align(&first, &second).unwrap().score()
        );
    }

    #[test]
    fn local_scores_equal_traceback_scores(
        first in dna(),
        second in dna(),
        matched in 0_i64..=6,
        mismatched in -6_i64..=0,
        gap_score in -6_i64..=0,
    ) {
        let scoring = MatchMismatch::new(matched, mismatched);
        let bounded = SmithWaterman::new(scoring, gap(gap_score));
        prop_assert_eq!(
            bounded.score(&first, &second).unwrap(),
            bounded.align(&first, &second).unwrap().score()
        );
        let checked = SmithWaterman::new(Unbounded(scoring), gap(gap_score));
        prop_assert_eq!(
            checked.score(&first, &second).unwrap(),
            checked.align(&first, &second).unwrap().score()
        );
    }
}

#[test]
fn overflow_errors_name_the_same_cell() {
    let first = DnaSequence::from_ascii("ACGTACGTACGT").unwrap();
    let second = DnaSequence::from_ascii("ACGTACGT").unwrap();
    let positive = MatchMismatch::new(Score::MAX / 3, 0);
    let negative = MatchMismatch::new(0, Score::MIN / 3);

    for (scoring, gap_score) in [(positive, -1), (negative, Score::MIN / 3)] {
        let global = NeedlemanWunsch::new(scoring, gap(gap_score));
        let expected = format!("{:?}", global.align(&first, &second).unwrap_err());
        assert_eq!(
            format!("{:?}", global.score(&first, &second).unwrap_err()),
            expected
        );
        assert!(expected.starts_with("ScoreOverflow"), "{expected}");
    }

    let local = SmithWaterman::new(positive, gap(-1));
    let expected = format!("{:?}", local.align(&first, &second).unwrap_err());
    assert_eq!(
        format!("{:?}", local.score(&first, &second).unwrap_err()),
        expected
    );
    assert!(expected.starts_with("ScoreOverflow"), "{expected}");
    // The local zero floor keeps extreme negative scores representable.
    let floored = SmithWaterman::new(negative, gap(Score::MIN / 3));
    assert_eq!(
        floored.score(&first, &second).unwrap(),
        floored.align(&first, &second).unwrap().score()
    );
}

#[test]
fn longer_inputs_keep_equal_scores() {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut sequence = |length: usize| {
        let bytes = (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[usize::try_from(state % 4).unwrap()]
            })
            .collect::<Vec<_>>();
        DnaSequence::from_ascii(&bytes).unwrap()
    };
    let (first, second) = (sequence(2_000), sequence(2_500));
    let scoring = MatchMismatch::new(2, -1);
    let global = NeedlemanWunsch::new(scoring, gap(-2));
    assert_eq!(
        global.score(&first, &second).unwrap(),
        global.align(&first, &second).unwrap().score()
    );
    let local = SmithWaterman::new(scoring, gap(-2));
    assert_eq!(
        local.score(&first, &second).unwrap(),
        local.align(&first, &second).unwrap().score()
    );
}
