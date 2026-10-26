//! Exhaustive non-DP oracle over every alignment of every short substring pair.

use align_core::{Dna, DnaSequence, LinearGap, MatchMismatch, SubstitutionScorer, Symbol};
use align_pairwise::SmithWaterman;

fn exhaustive_global_score(
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
                + exhaustive_global_score(
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
                + exhaustive_global_score(
                    first,
                    second,
                    first_index + 1,
                    second_index,
                    scorer,
                    gap,
                ),
        );
    }
    if second_index < second.len() {
        best = best.max(
            gap.score()
                + exhaustive_global_score(
                    first,
                    second,
                    first_index,
                    second_index + 1,
                    scorer,
                    gap,
                ),
        );
    }
    best
}

fn exhaustive_local_score(
    first: &[Symbol<Dna>],
    second: &[Symbol<Dna>],
    scorer: &MatchMismatch,
    gap: LinearGap,
) -> i64 {
    let mut best = 0;
    for first_start in 0..=first.len() {
        for first_end in first_start..=first.len() {
            for second_start in 0..=second.len() {
                for second_end in second_start..=second.len() {
                    best = best.max(exhaustive_global_score(
                        &first[first_start..first_end],
                        &second[second_start..second_end],
                        0,
                        0,
                        scorer,
                        gap,
                    ));
                }
            }
        }
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
fn local_dynamic_programming_matches_exhaustive_enumeration_through_length_four() {
    let sequences = all_binary_sequences(4);
    let scorer = MatchMismatch::new(2, -1);
    let gap = LinearGap::new(-2).expect("valid gap");
    let aligner = SmithWaterman::new(scorer, gap);

    for first in &sequences {
        for second in &sequences {
            let expected =
                exhaustive_local_score(first.as_slice(), second.as_slice(), &scorer, gap);
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
