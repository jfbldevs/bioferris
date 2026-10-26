//! Exhaustive complete-path oracle for small affine-gap global alignments.

use align_core::{AffineGap, AlignmentOperation, DnaSequence, MatchMismatch};
use align_pairwise::{Gotoh, verify_global_affine_alignment};

const MATCH_SCORE: i64 = 2;
const MISMATCH_SCORE: i64 = -3;
const OPEN_SCORE: i64 = -4;
const EXTENSION_SCORE: i64 = -1;

fn sequences_of_length(length: usize) -> Vec<Vec<u8>> {
    let count = 1_usize << length;
    (0..count)
        .map(|bits| {
            (0..length)
                .map(|offset| {
                    if bits & (1_usize << offset) == 0 {
                        b'A'
                    } else {
                        b'C'
                    }
                })
                .collect()
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn enumerate_paths(
    first: &[u8],
    second: &[u8],
    first_index: usize,
    second_index: usize,
    previous: Option<AlignmentOperation>,
    score: i64,
    best: &mut Option<i64>,
) {
    if first_index == first.len() && second_index == second.len() {
        *best = Some(best.map_or(score, |current| current.max(score)));
        return;
    }

    if first_index < first.len() && second_index < second.len() {
        let substitution = if first[first_index] == second[second_index] {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        };
        enumerate_paths(
            first,
            second,
            first_index + 1,
            second_index + 1,
            Some(if substitution == MATCH_SCORE {
                AlignmentOperation::Match
            } else {
                AlignmentOperation::Substitution
            }),
            score + substitution,
            best,
        );
    }

    if first_index < first.len() {
        let gap_score = if previous == Some(AlignmentOperation::Deletion) {
            EXTENSION_SCORE
        } else {
            OPEN_SCORE
        };
        enumerate_paths(
            first,
            second,
            first_index + 1,
            second_index,
            Some(AlignmentOperation::Deletion),
            score + gap_score,
            best,
        );
    }

    if second_index < second.len() {
        let gap_score = if previous == Some(AlignmentOperation::Insertion) {
            EXTENSION_SCORE
        } else {
            OPEN_SCORE
        };
        enumerate_paths(
            first,
            second,
            first_index,
            second_index + 1,
            Some(AlignmentOperation::Insertion),
            score + gap_score,
            best,
        );
    }
}

fn exhaustive_score(first: &[u8], second: &[u8]) -> i64 {
    let mut best = None;
    enumerate_paths(first, second, 0, 0, None, 0, &mut best);
    best.expect("every sequence pair has at least one complete operation path")
}

#[test]
fn gotoh_matches_independent_complete_path_enumeration_through_length_four() {
    let substitution = MatchMismatch::new(MATCH_SCORE, MISMATCH_SCORE);
    let gap = AffineGap::new(OPEN_SCORE, EXTENSION_SCORE).expect("valid affine gap");
    let aligner = Gotoh::new(substitution, gap);
    let mut checked_pairs = 0_usize;

    for first_len in 0..=4 {
        for second_len in 0..=4 {
            for first_bytes in sequences_of_length(first_len) {
                for second_bytes in sequences_of_length(second_len) {
                    let expected = exhaustive_score(&first_bytes, &second_bytes);
                    let first = DnaSequence::from_ascii(&first_bytes).expect("binary DNA is valid");
                    let second =
                        DnaSequence::from_ascii(&second_bytes).expect("binary DNA is valid");
                    let observed = aligner
                        .align(&first, &second)
                        .expect("small alignment succeeds");
                    assert_eq!(
                        observed.score(),
                        expected,
                        "score mismatch for first={:?}, second={:?}",
                        String::from_utf8_lossy(&first_bytes),
                        String::from_utf8_lossy(&second_bytes)
                    );
                    verify_global_affine_alignment(&first, &second, &observed, &substitution, gap)
                        .expect("traceback independently replays");
                    checked_pairs += 1;
                }
            }
        }
    }

    assert_eq!(checked_pairs, 961);
}
