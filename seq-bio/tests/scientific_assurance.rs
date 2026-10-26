//! Independent path enumeration and cross-domain biological invariants.

use bioferris::align::advanced::core::{AlignmentOperation, PairwiseAlignment, Score};
use bioferris::prelude::*;

#[derive(Clone, Copy, Debug)]
struct Scoring {
    matching: Score,
    mismatch: Score,
    open: Score,
    extend: Score,
}

fn words(alphabet: &[u8], maximum: usize) -> Vec<String> {
    let mut words = vec![String::new()];
    let mut layer = vec![String::new()];
    for _ in 0..maximum {
        layer = layer
            .iter()
            .flat_map(|prefix| {
                alphabet
                    .iter()
                    .map(move |&byte| format!("{prefix}{}", char::from(byte)))
            })
            .collect();
        words.extend(layer.iter().cloned());
    }
    words
}

// Enumerate complete monotone paths, not the production DP recurrence. `local`
// permits termination at every vertex; callers separately enumerate all starts.
fn path_score(
    first: &[u8],
    second: &[u8],
    previous_gap: u8,
    scoring: Scoring,
    local: bool,
) -> Score {
    if first.is_empty() && second.is_empty() {
        return 0;
    }
    let mut best = if local { 0 } else { Score::MIN / 4 };
    if let (Some((&a, rest_a)), Some((&b, rest_b))) = (first.split_first(), second.split_first()) {
        let substitution = if a == b {
            scoring.matching
        } else {
            scoring.mismatch
        };
        best = best.max(substitution + path_score(rest_a, rest_b, 0, scoring, local));
    }
    if let Some((_, rest)) = first.split_first() {
        let gap = if previous_gap == 1 {
            scoring.extend
        } else {
            scoring.open
        };
        best = best.max(gap + path_score(rest, second, 1, scoring, local));
    }
    if let Some((_, rest)) = second.split_first() {
        let gap = if previous_gap == 2 {
            scoring.extend
        } else {
            scoring.open
        };
        best = best.max(gap + path_score(first, rest, 2, scoring, local));
    }
    best
}

fn optimum(first: &[u8], second: &[u8], scoring: Scoring, local: bool) -> Score {
    if !local {
        return path_score(first, second, 0, scoring, false);
    }
    (0..=first.len())
        .flat_map(|a| (0..=second.len()).map(move |b| (a, b)))
        .map(|(a, b)| path_score(&first[a..], &second[b..], 0, scoring, true))
        .max()
        .unwrap()
}

fn check_trace(
    result: &PairwiseAlignment,
    first: &[u8],
    second: &[u8],
    scoring: Scoring,
    local: bool,
) {
    let first_range = result.first_range();
    let second_range = result.second_range();
    assert!(first_range.start <= first_range.end && first_range.end <= first.len());
    assert!(second_range.start <= second_range.end && second_range.end <= second.len());
    if !local {
        assert_eq!(first_range, 0..first.len());
        assert_eq!(second_range, 0..second.len());
    }
    let (mut a, mut b, mut previous_gap, mut score) = (first_range.start, second_range.start, 0, 0);
    for &operation in result.operations() {
        match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                assert!(a < first_range.end && b < second_range.end);
                let equal = first[a] == second[b];
                assert_eq!(operation == AlignmentOperation::Match, equal);
                score += if equal {
                    scoring.matching
                } else {
                    scoring.mismatch
                };
                a += 1;
                b += 1;
                previous_gap = 0;
            }
            AlignmentOperation::Deletion => {
                assert!(a < first_range.end);
                score += if previous_gap == 1 {
                    scoring.extend
                } else {
                    scoring.open
                };
                a += 1;
                previous_gap = 1;
            }
            AlignmentOperation::Insertion => {
                assert!(b < second_range.end);
                score += if previous_gap == 2 {
                    scoring.extend
                } else {
                    scoring.open
                };
                b += 1;
                previous_gap = 2;
            }
        }
    }
    assert_eq!((a, b), (first_range.end, second_range.end));
    assert_eq!(score, result.score());
    assert_eq!(score, optimum(first, second, scoring, local));
}

#[test]
fn exhaustive_linear_paths_verify_scores_tracebacks_and_score_only_results() {
    let inputs = words(b"AC", 3);
    let mut comparisons = 0;
    for &(matching, mismatch, gap) in &[(2, -1, -2), (3, -2, -1), (1, 0, -1)] {
        let scoring = Scoring {
            matching,
            mismatch,
            open: gap,
            extend: gap,
        };
        for first in &inputs {
            for second in &inputs {
                let a = dna(first).unwrap();
                let b = dna(second).unwrap();
                for (aligner, local) in [
                    (Pairwise::needleman_wunsch(), false),
                    (Pairwise::smith_waterman(), true),
                ] {
                    let aligner = aligner.scores(matching, mismatch).gap(gap);
                    let result = aligner.align(&a, &b).unwrap();
                    check_trace(&result, first.as_bytes(), second.as_bytes(), scoring, local);
                    assert_eq!(aligner.score(&a, &b).unwrap(), result.score());
                    assert_eq!(
                        aligner.align(&a, &b).unwrap(),
                        result,
                        "deterministic traceback"
                    );
                    comparisons += 1;
                }
            }
        }
    }
    assert_eq!(comparisons, 1_350);
}

#[test]
fn exhaustive_affine_paths_verify_global_and_local_gap_run_costs() {
    let inputs = words(b"AC", 3);
    let mut comparisons = 0;
    for scoring in [
        Scoring {
            matching: 2,
            mismatch: -1,
            open: -3,
            extend: -1,
        },
        Scoring {
            matching: 3,
            mismatch: -2,
            open: -4,
            extend: -1,
        },
        Scoring {
            matching: 1,
            mismatch: 0,
            open: -1,
            extend: -1,
        },
    ] {
        for first in &inputs {
            for second in &inputs {
                let a = dna(first).unwrap();
                let b = dna(second).unwrap();
                for local in [false, true] {
                    let aligner = Pairwise::gotoh()
                        .scores(scoring.matching, scoring.mismatch)
                        .gaps(scoring.open, scoring.extend);
                    let aligner = if local { aligner.local() } else { aligner };
                    let result = aligner.align(&a, &b).unwrap();
                    check_trace(&result, first.as_bytes(), second.as_bytes(), scoring, local);
                    assert_eq!(aligner.score(&a, &b).unwrap(), result.score());
                    assert_eq!(aligner.align(&a, &b).unwrap(), result);
                    comparisons += 1;
                }
            }
        }
    }
    assert_eq!(comparisons, 1_350);
}

#[test]
fn exhaustive_edit_and_wavefront_objectives_match_independent_path_enumeration() {
    let inputs = words(b"AC", 3);
    let mut comparisons = 0;
    for first in &inputs {
        for second in &inputs {
            let a = dna(first).unwrap();
            let b = dna(second).unwrap();
            let edit_scoring = Scoring {
                matching: 0,
                mismatch: -1,
                open: -1,
                extend: -1,
            };
            let edit = Pairwise::myers().align(&a, &b).unwrap();
            check_trace(
                &edit,
                first.as_bytes(),
                second.as_bytes(),
                edit_scoring,
                false,
            );
            assert_eq!(
                Score::try_from(Pairwise::myers().distance(&a, &b).unwrap()).unwrap(),
                -edit.score()
            );
            for (mismatch, open, extend) in [(2, 3, 1), (3, 1, 2), (1, 0, 1)] {
                // WFA uses open + length*extension; the path oracle charges
                // open+extension for the first gap symbol, then extension.
                let scoring = Scoring {
                    matching: 0,
                    mismatch: -mismatch,
                    open: -(open + extend),
                    extend: -extend,
                };
                let wave = Pairwise::wfa().penalties(mismatch, open, extend);
                let result = wave.align(&a, &b).unwrap();
                check_trace(
                    result.alignment(),
                    first.as_bytes(),
                    second.as_bytes(),
                    scoring,
                    false,
                );
                assert_eq!(wave.penalty(&a, &b).unwrap(), -result.alignment().score());
                let adaptive = Pairwise::adaptive_affine().penalties(mismatch, open, extend);
                let result = adaptive.align(&a, &b).unwrap();
                check_trace(
                    result.alignment(),
                    first.as_bytes(),
                    second.as_bytes(),
                    scoring,
                    false,
                );
                assert_eq!(
                    adaptive.penalty(&a, &b).unwrap().penalty(),
                    -result.alignment().score()
                );
                comparisons += 2;
            }
            comparisons += 1;
        }
    }
    assert_eq!(comparisons, 1_575);
}

fn reverse_range(range: std::ops::Range<usize>, length: usize) -> (usize, usize) {
    (length - range.end, length - range.start)
}

#[test]
fn orf_coordinates_and_products_commute_with_whole_source_reverse_complement() {
    let codons = [
        ("GCT", 'A'),
        ("TGG", 'W'),
        ("AAA", 'K'),
        ("GGA", 'G'),
        ("CAT", 'H'),
    ];
    let mut comparisons = 0;
    for &(first, first_residue) in &codons {
        for &(second, second_residue) in &codons {
            for prefix in ["", "C", "CC"] {
                let input = dna(format!("{prefix}ATG{first}{second}TAAGG")).unwrap();
                let reverse = input.reverse_complement().unwrap();
                let original = input.find_orfs().unwrap();
                assert!(original.iter().any(|orf| orf.strand() == Strand::Forward
                    && orf.range() == (prefix.len()..prefix.len() + 9)
                    && orf.peptide() == format!("M{first_residue}{second_residue}")));
                let mut expected: Vec<_> = original
                    .iter()
                    .map(|orf| {
                        (
                            orf.strand() == Strand::Forward,
                            orf.reading_frame().offset(),
                            reverse_range(orf.range(), input.len()),
                            orf.stop_range()
                                .map(|range| reverse_range(range, input.len())),
                            orf.peptide().to_owned(),
                        )
                    })
                    .collect();
                let mut actual: Vec<_> = reverse
                    .find_orfs()
                    .unwrap()
                    .iter()
                    .map(|orf| {
                        (
                            orf.strand() == Strand::Reverse,
                            orf.reading_frame().offset(),
                            (orf.range().start, orf.range().end),
                            orf.stop_range().map(|range| (range.start, range.end)),
                            orf.peptide().to_owned(),
                        )
                    })
                    .collect();
                expected.sort_unstable();
                actual.sort_unstable();
                assert_eq!(actual, expected);
                assert_eq!(input.transcribe().unwrap().find_orfs().unwrap(), original);
                comparisons += 1;
            }
        }
    }
    assert_eq!(comparisons, 75);
}

#[test]
fn iupac_hits_obey_reverse_complement_bijection_without_collapsing_strands() {
    let source = dna("ARYNACGTBDHVWSKMNACGT").unwrap();
    let reverse = source.reverse_complement().unwrap();
    let mut comparisons = 0;
    for pattern in words(b"ACGTRYSWKMBDHVN", 2)
        .into_iter()
        .filter(|pattern| !pattern.is_empty())
    {
        let motif = dna_motif(pattern.as_bytes()).unwrap();
        for ambiguity in [MotifInputAmbiguity::Definite, MotifInputAmbiguity::Possible] {
            let options = MotifSearchOptions {
                ambiguity,
                ..Default::default()
            };
            let mut expected: Vec<_> = motif
                .find_with(&source, options)
                .unwrap()
                .iter()
                .map(|hit| {
                    (
                        reverse_range(hit.range(), source.len()),
                        hit.strand() == Strand::Forward,
                    )
                })
                .collect();
            let mut actual: Vec<_> = motif
                .find_with(&reverse, options)
                .unwrap()
                .iter()
                .map(|hit| {
                    (
                        (hit.range().start, hit.range().end),
                        hit.strand() == Strand::Reverse,
                    )
                })
                .collect();
            expected.sort_unstable();
            actual.sort_unstable();
            assert_eq!(actual, expected, "pattern={pattern}, policy={ambiguity:?}");
            comparisons += 1;
        }
    }
    assert_eq!(comparisons, 480);
}

#[test]
fn ordered_pairwise_batches_are_trace_exact_at_one_two_and_four_threads() {
    let inputs: Vec<_> = words(b"AC", 3)
        .iter()
        .map(|word| dna(word).unwrap())
        .collect();
    let pairs: Vec<_> = (0..32)
        .map(|index| {
            (
                &inputs[index % inputs.len()],
                &inputs[(index * 7 + 3) % inputs.len()],
            )
        })
        .collect();
    for (aligner, local) in [
        (Pairwise::needleman_wunsch(), false),
        (Pairwise::smith_waterman(), true),
    ] {
        let baseline = aligner
            .threads(1)
            .align_many(pairs.iter().copied())
            .unwrap();
        assert_eq!(baseline.len(), pairs.len());
        for (result, (first, second)) in baseline.iter().zip(&pairs) {
            check_trace(
                result,
                first.to_string().as_bytes(),
                second.to_string().as_bytes(),
                Scoring {
                    matching: 2,
                    mismatch: -1,
                    open: -2,
                    extend: -2,
                },
                local,
            );
        }
        for threads in [1, 2, 4] {
            assert_eq!(
                aligner
                    .threads(threads)
                    .align_many(pairs.iter().copied())
                    .unwrap(),
                baseline,
                "linear {local:?}, threads={threads}"
            );
        }
    }
    let baseline = Pairwise::gotoh()
        .threads(1)
        .align_many(pairs.iter().copied())
        .unwrap();
    assert_eq!(baseline.len(), pairs.len());
    for (result, (first, second)) in baseline.iter().zip(&pairs) {
        check_trace(
            result,
            first.to_string().as_bytes(),
            second.to_string().as_bytes(),
            Scoring {
                matching: 2,
                mismatch: -1,
                open: -5,
                extend: -1,
            },
            false,
        );
    }
    for threads in [1, 2, 4] {
        assert_eq!(
            Pairwise::gotoh()
                .threads(threads)
                .align_many(pairs.iter().copied())
                .unwrap(),
            baseline,
            "affine threads={threads}"
        );
    }
}

#[test]
fn published_scientific_example_output_matches_literal_inputs() {
    let first = dna("AC").unwrap();
    let second = dna("AG").unwrap();
    let alignment = Pairwise::needleman_wunsch().align(&first, &second).unwrap();
    let source = dna("TTATTTCAT").unwrap();
    let orfs = source.find_orfs().unwrap();
    let orf = &orfs[0];
    let observed = format!(
        "Global score: {}\nTrace columns: {}\nReverse ORF: {} {:?} stop={:?}\n",
        alignment.score(),
        alignment.columns(),
        orf.peptide(),
        orf.range(),
        orf.stop_range()
    );
    let guide = include_str!("../../docs/testing/contracts-and-science.md");
    let expected = guide
        .split_once("<!-- assurance-rust-output -->\n```text\n")
        .unwrap()
        .1
        .split_once("\n```")
        .unwrap()
        .0;
    assert_eq!(observed, format!("{expected}\n"));
}
