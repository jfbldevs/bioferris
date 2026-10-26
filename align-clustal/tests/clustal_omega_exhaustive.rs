//! Exhaustive small-alphabet checks with a path-enumeration distance oracle.

use align_clustal::{ClustalOmega, ClustalOmegaParameters, OmegaEngineMode};
use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;

const DISTANCE_SCALE: usize = 1_000_000_000;
const GAP_EVENT_PENALTY: isize = 3;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AlignmentStep {
    Both,
    FirstOnly,
    SecondOnly,
}

fn binary_strings(maximum_length: usize) -> Vec<Vec<u8>> {
    let mut strings = vec![Vec::new()];
    for length in 1..=maximum_length {
        for bits in 0..(1_usize << length) {
            strings.push(
                (0..length)
                    .map(|offset| {
                        if bits & (1 << offset) == 0 {
                            b'A'
                        } else {
                            b'C'
                        }
                    })
                    .collect(),
            );
        }
    }
    strings
}

fn score_path(first: &[u8], second: &[u8], path: &[AlignmentStep]) -> isize {
    let mut first_index = 0_usize;
    let mut second_index = 0_usize;
    let mut score = 0_isize;
    let mut previous = None;
    for step in path.iter().copied() {
        match step {
            AlignmentStep::Both => {
                score += isize::from(first[first_index] == second[second_index]);
                first_index += 1;
                second_index += 1;
            }
            AlignmentStep::FirstOnly => {
                if previous != Some(step) && second_index != 0 && second_index != second.len() {
                    score -= GAP_EVENT_PENALTY;
                }
                first_index += 1;
            }
            AlignmentStep::SecondOnly => {
                if previous != Some(step) && first_index != 0 && first_index != first.len() {
                    score -= GAP_EVENT_PENALTY;
                }
                second_index += 1;
            }
        }
        previous = Some(step);
    }
    assert_eq!(first_index, first.len());
    assert_eq!(second_index, second.len());
    score
}

fn enumerate_paths(
    first: &[u8],
    second: &[u8],
    first_index: usize,
    second_index: usize,
    path: &mut Vec<AlignmentStep>,
    best_score: &mut isize,
) {
    if first_index == first.len() && second_index == second.len() {
        *best_score = (*best_score).max(score_path(first, second, path));
        return;
    }
    if first_index < first.len() && second_index < second.len() {
        path.push(AlignmentStep::Both);
        enumerate_paths(
            first,
            second,
            first_index + 1,
            second_index + 1,
            path,
            best_score,
        );
        path.pop();
    }
    if first_index < first.len() {
        path.push(AlignmentStep::FirstOnly);
        enumerate_paths(
            first,
            second,
            first_index + 1,
            second_index,
            path,
            best_score,
        );
        path.pop();
    }
    if second_index < second.len() {
        path.push(AlignmentStep::SecondOnly);
        enumerate_paths(
            first,
            second,
            first_index,
            second_index + 1,
            path,
            best_score,
        );
        path.pop();
    }
}

#[allow(clippy::cast_precision_loss)]
fn exhaustive_distance(first: &[u8], second: &[u8]) -> f64 {
    let shorter = first.len().min(second.len());
    if shorter == 0 {
        return if first.is_empty() && second.is_empty() {
            0.0
        } else {
            1.0
        };
    }
    let mut best_score = isize::MIN;
    enumerate_paths(first, second, 0, 0, &mut Vec::new(), &mut best_score);
    let similarity = usize::try_from(best_score.max(0)).unwrap();
    let numerator = (shorter - similarity) * DISTANCE_SCALE;
    let fixed = (numerator + shorter / 2) / shorter;
    fixed as f64 / DISTANCE_SCALE as f64
}

#[test]
fn rolling_one_tuple_score_matches_every_global_path_through_length_four() {
    let strings = binary_strings(4);
    let aligner = ClustalOmega::<Dna>::published();
    let mut cases = 0_usize;
    for first in &strings {
        for second in &strings {
            let first_sequence = DnaSequence::from_ascii(first).unwrap();
            let second_sequence = DnaSequence::from_ascii(second).unwrap();
            let observed = aligner
                .pair_distance(&first_sequence, &second_sequence)
                .unwrap();
            let expected = exhaustive_distance(first, second);
            assert!(
                (observed - expected).abs() <= 5.0e-10,
                "one-tuple disagreement for {:?} versus {:?}: {observed} versus {expected}",
                String::from_utf8_lossy(first),
                String::from_utf8_lossy(second)
            );
            cases += 1;
        }
    }
    assert_eq!(cases, strings.len().pow(2));
}

#[test]
fn every_ordered_triple_through_length_two_has_a_valid_traceback() {
    let strings = binary_strings(2);
    let parameters =
        ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Viterbi);
    let aligner = ClustalOmega::<Dna>::with_parameters(parameters).unwrap();
    let mut cases = 0_usize;
    for first in &strings {
        for second in &strings {
            for third in &strings {
                let records = [
                    SequenceRecord::new("a", DnaSequence::from_ascii(first).unwrap()).unwrap(),
                    SequenceRecord::new("b", DnaSequence::from_ascii(second).unwrap()).unwrap(),
                    SequenceRecord::new("c", DnaSequence::from_ascii(third).unwrap()).unwrap(),
                ];
                let result = aligner.align(&records).unwrap();
                assert_eq!(
                    result.alignment().rows()[0].to_ungapped_ascii().as_bytes(),
                    first
                );
                assert_eq!(
                    result.alignment().rows()[1].to_ungapped_ascii().as_bytes(),
                    second
                );
                assert_eq!(
                    result.alignment().rows()[2].to_ungapped_ascii().as_bytes(),
                    third
                );
                assert_eq!(result.progressive_steps().len(), 2);
                cases += 1;
            }
        }
    }
    assert_eq!(cases, strings.len().pow(3));
}
