//! Independent complete-path oracles for affine local and four-end alignment.

use align_core::{AffineGap, AlignmentOperation, DnaSequence, LinearGap, MatchMismatch};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, FreeEndGaps, SmithWaterman,
    verify_affine_end_gap_alignment, verify_local_affine_alignment,
};

const MATCH_SCORE: i64 = 2;
const MISMATCH_SCORE: i64 = -3;
const OPEN_SCORE: i64 = -4;
const EXTENSION_SCORE: i64 = -1;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct OracleResult {
    score: i64,
    first_end: usize,
    second_end: usize,
}

impl OracleResult {
    fn observe(&mut self, score: i64, first_end: usize, second_end: usize) {
        let candidate_consumed = first_end + second_end;
        let selected_consumed = self.first_end + self.second_end;
        if score > self.score
            || (score == self.score && candidate_consumed > selected_consumed)
            || (score == self.score
                && candidate_consumed == selected_consumed
                && first_end > self.first_end)
        {
            *self = Self {
                score,
                first_end,
                second_end,
            };
        }
    }
}

fn sequences_of_length(length: usize) -> Vec<Vec<u8>> {
    (0..(1_usize << length))
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
fn enumerate_complete_paths(
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
        let pair_score = if first[first_index] == second[second_index] {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        };
        enumerate_complete_paths(
            first,
            second,
            first_index + 1,
            second_index + 1,
            Some(if pair_score == MATCH_SCORE {
                AlignmentOperation::Match
            } else {
                AlignmentOperation::Substitution
            }),
            score + pair_score,
            best,
        );
    }
    if first_index < first.len() {
        enumerate_complete_paths(
            first,
            second,
            first_index + 1,
            second_index,
            Some(AlignmentOperation::Deletion),
            score
                + if previous == Some(AlignmentOperation::Deletion) {
                    EXTENSION_SCORE
                } else {
                    OPEN_SCORE
                },
            best,
        );
    }
    if second_index < second.len() {
        enumerate_complete_paths(
            first,
            second,
            first_index,
            second_index + 1,
            Some(AlignmentOperation::Insertion),
            score
                + if previous == Some(AlignmentOperation::Insertion) {
                    EXTENSION_SCORE
                } else {
                    OPEN_SCORE
                },
            best,
        );
    }
}

fn complete_path_score(first: &[u8], second: &[u8]) -> i64 {
    let mut best = None;
    enumerate_complete_paths(first, second, 0, 0, None, 0, &mut best);
    best.expect("every pair of ranges has a complete path")
}

fn permitted_start(first_start: usize, second_start: usize, ends: FreeEndGaps) -> bool {
    if first_start == 0 && second_start == 0 {
        true
    } else if second_start == 0 {
        ends.first_start()
    } else if first_start == 0 {
        ends.second_start()
    } else {
        false
    }
}

fn permitted_end(
    first_end: usize,
    first_len: usize,
    second_end: usize,
    second_len: usize,
    ends: FreeEndGaps,
) -> bool {
    if first_end == first_len && second_end == second_len {
        true
    } else if second_end == second_len {
        ends.first_end()
    } else if first_end == first_len {
        ends.second_end()
    } else {
        false
    }
}

fn exhaustive_end_gap(first: &[u8], second: &[u8], ends: FreeEndGaps) -> OracleResult {
    let mut best: Option<OracleResult> = None;
    for first_start in 0..=first.len() {
        for second_start in 0..=second.len() {
            if !permitted_start(first_start, second_start, ends) {
                continue;
            }
            for first_end in first_start..=first.len() {
                for second_end in second_start..=second.len() {
                    if !permitted_end(first_end, first.len(), second_end, second.len(), ends) {
                        continue;
                    }
                    let score = complete_path_score(
                        &first[first_start..first_end],
                        &second[second_start..second_end],
                    );
                    match &mut best {
                        Some(current) => current.observe(score, first_end, second_end),
                        None => {
                            best = Some(OracleResult {
                                score,
                                first_end,
                                second_end,
                            });
                        }
                    }
                }
            }
        }
    }
    best.expect("the complete-range candidate is always present")
}

fn exhaustive_local(first: &[u8], second: &[u8]) -> i64 {
    let mut best = 0;
    for first_start in 0..=first.len() {
        for first_end in first_start..=first.len() {
            for second_start in 0..=second.len() {
                for second_end in second_start..=second.len() {
                    best = best.max(complete_path_score(
                        &first[first_start..first_end],
                        &second[second_start..second_end],
                    ));
                }
            }
        }
    }
    best
}

#[test]
fn all_sixteen_affine_end_modes_match_complete_paths_through_length_three() {
    let scorer = MatchMismatch::new(MATCH_SCORE, MISMATCH_SCORE);
    let gap = AffineGap::new(OPEN_SCORE, EXTENSION_SCORE).expect("valid gap");
    let mut checked = 0_usize;
    for bits in 0_u8..16 {
        let ends = FreeEndGaps::new(bits & 1 != 0, bits & 2 != 0, bits & 4 != 0, bits & 8 != 0);
        let aligner = AffineEndGapAligner::new(scorer, gap, ends);
        for first_len in 0..=3 {
            for second_len in 0..=3 {
                for first_bytes in sequences_of_length(first_len) {
                    for second_bytes in sequences_of_length(second_len) {
                        let expected = exhaustive_end_gap(&first_bytes, &second_bytes, ends);
                        let first = DnaSequence::from_ascii(&first_bytes).expect("first DNA");
                        let second = DnaSequence::from_ascii(&second_bytes).expect("second DNA");
                        let observed = aligner.align(&first, &second).expect("affine end-gap");
                        assert_eq!(
                            observed.score(),
                            expected.score,
                            "mode={bits:04b}, first={:?}, second={:?}",
                            String::from_utf8_lossy(&first_bytes),
                            String::from_utf8_lossy(&second_bytes)
                        );
                        assert_eq!(
                            (observed.first_range().end, observed.second_range().end),
                            (expected.first_end, expected.second_end),
                            "endpoint mode={bits:04b}, first={:?}, second={:?}",
                            String::from_utf8_lossy(&first_bytes),
                            String::from_utf8_lossy(&second_bytes)
                        );
                        verify_affine_end_gap_alignment(
                            &first, &second, &observed, &scorer, gap, ends,
                        )
                        .expect("independent affine end-gap replay");
                        checked += 1;
                    }
                }
            }
        }
    }
    assert_eq!(checked, 3_600);
}

#[test]
fn affine_local_matches_every_substring_path_through_length_four() {
    let scorer = MatchMismatch::new(MATCH_SCORE, MISMATCH_SCORE);
    let gap = AffineGap::new(OPEN_SCORE, EXTENSION_SCORE).expect("valid gap");
    let aligner = AffineLocalAligner::new(scorer, gap);
    let sequences = (0..=4).flat_map(sequences_of_length).collect::<Vec<_>>();
    let mut checked = 0_usize;
    for first_bytes in &sequences {
        for second_bytes in &sequences {
            let expected = exhaustive_local(first_bytes, second_bytes);
            let first = DnaSequence::from_ascii(first_bytes).expect("first DNA");
            let second = DnaSequence::from_ascii(second_bytes).expect("second DNA");
            let observed = aligner.align(&first, &second).expect("affine local");
            assert_eq!(observed.score(), expected, "first={first}, second={second}");
            verify_local_affine_alignment(&first, &second, &observed, &scorer, gap)
                .expect("independent affine local replay");
            checked += 1;
        }
    }
    assert_eq!(checked, 961);
}

#[test]
fn equal_open_and_extension_local_reduces_to_smith_waterman() {
    let scorer = MatchMismatch::new(2, -3);
    let affine_gap = AffineGap::new(-2, -2).expect("linear affine specialization");
    let linear_gap = LinearGap::new(-2).expect("linear gap");
    let affine = AffineLocalAligner::new(scorer, affine_gap);
    let linear = SmithWaterman::new(scorer, linear_gap);
    let sequences = (0..=4).flat_map(sequences_of_length).collect::<Vec<_>>();
    for first_bytes in &sequences {
        for second_bytes in &sequences {
            let first = DnaSequence::from_ascii(first_bytes).expect("first DNA");
            let second = DnaSequence::from_ascii(second_bytes).expect("second DNA");
            let affine_result = affine.align(&first, &second).expect("affine local");
            let linear_result = linear.align(&first, &second).expect("linear local");
            assert_eq!(
                affine_result.score(),
                linear_result.score(),
                "first={first}, second={second}"
            );
            assert_eq!(
                (
                    affine_result.first_range().end,
                    affine_result.second_range().end
                ),
                (
                    linear_result.first_range().end,
                    linear_result.second_range().end
                ),
                "endpoint first={first}, second={second}"
            );
        }
    }
}
