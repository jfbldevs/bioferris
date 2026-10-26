//! Exhaustive, non-DP score oracle for every short binary sequence pair.

use align_core::{Dna, DnaSequence, LinearGap, MatchMismatch, SubstitutionScorer, Symbol};
use align_pairwise::NeedlemanWunsch;

fn exhaustive_score(
    first: &[Symbol<Dna>],
    second: &[Symbol<Dna>],
    first_index: usize,
    second_index: usize,
    scorer: &MatchMismatch,
    gap: LinearGap,
) -> i64 {
    if first_index == first.len() && second_index == second.len() {
        return 0;
    }

    let mut best = i64::MIN;
    if first_index < first.len() && second_index < second.len() {
        best = best.max(
            scorer.score(first[first_index], second[second_index])
                + exhaustive_score(
                    first,
                    second,
                    first_index + 1,
                    second_index + 1,
                    scorer,
                    gap,
                ),
        );
    }
    if first_index < first.len() {
        best = best.max(
            gap.score()
                + exhaustive_score(first, second, first_index + 1, second_index, scorer, gap),
        );
    }
    if second_index < second.len() {
        best = best.max(
            gap.score()
                + exhaustive_score(first, second, first_index, second_index + 1, scorer, gap),
        );
    }
    best
}

fn all_binary_sequences(max_length: usize) -> Vec<DnaSequence> {
    let mut sequences = Vec::new();
    for length in 0..=max_length {
        for bits in 0..(1_usize << length) {
            let raw: Vec<u8> = (0..length)
                .map(|position| {
                    if bits & (1 << position) == 0 {
                        b'A'
                    } else {
                        b'C'
                    }
                })
                .collect();
            sequences.push(DnaSequence::from_ascii(raw).expect("generator emits valid DNA"));
        }
    }
    sequences
}

#[test]
fn dynamic_programming_matches_exhaustive_enumeration_through_length_four() {
    let sequences = all_binary_sequences(4);
    let scorer = MatchMismatch::new(2, -1);
    let gap = LinearGap::new(-2).expect("valid gap");
    let aligner = NeedlemanWunsch::new(scorer, gap);

    for first in &sequences {
        for second in &sequences {
            let expected =
                exhaustive_score(first.as_slice(), second.as_slice(), 0, 0, &scorer, gap);
            let actual = aligner
                .align(first, second)
                .expect("small alignment succeeds")
                .score();
            assert_eq!(
                actual, expected,
                "oracle disagreement for {first} versus {second}"
            );
        }
    }
}
