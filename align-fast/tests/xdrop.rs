//! Manual, exhaustive, boundary, seed, alphabet, and independent-oracle X-drop tests.

mod support;

use align_core::{
    AffineGap, AlignmentOperation, Dna, DnaSequence, MatchMismatch, ProteinSequence, RnaSequence,
    SubstitutionScorer, Symbol,
};
use align_fast::{XDrop, XDropDirection, XDropError, XDropSeed, XDropTermination};
use support::xdrop::{
    OracleSummary, OracleTermination, exhaustive_path_oracle, oracle, replay, unrestricted_oracle,
};

fn gap(open: i64, extension: i64) -> AffineGap {
    AffineGap::new(open, extension).expect("test gap scores are valid")
}

fn assert_matches_oracle(observed: align_fast::XDropScore, expected: OracleSummary) {
    assert_eq!(observed.score(), expected.score);
    assert_eq!(observed.first_end(), expected.first_end);
    assert_eq!(observed.second_end(), expected.second_end);
    assert_eq!(observed.evaluated_cells(), expected.evaluated_cells);
    assert_eq!(observed.retained_cells(), expected.retained_cells);
    assert_eq!(observed.antidiagonals(), expected.antidiagonals);
    assert_eq!(
        observed.maximum_frontier_width(),
        expected.maximum_frontier_width
    );
    assert_eq!(
        observed.termination(),
        match expected.termination {
            OracleTermination::MatrixExhausted => XDropTermination::MatrixExhausted,
            OracleTermination::FrontiersExhausted => XDropTermination::FrontiersExhausted,
        }
    );
}

fn binary_dna_strings(maximum_len: usize) -> Vec<DnaSequence> {
    let mut sequences = Vec::new();
    for length in 0..=maximum_len {
        for bits in 0..(1_usize << length) {
            let bytes: Vec<u8> = (0..length)
                .map(|position| {
                    if bits & (1 << position) == 0 {
                        b'A'
                    } else {
                        b'C'
                    }
                })
                .collect();
            sequences.push(DnaSequence::from_ascii(bytes).expect("generator emits valid DNA"));
        }
    }
    sequences
}

#[test]
fn empty_identity_difference_and_repetition_are_explicit() {
    let scoring = MatchMismatch::new(2, -3);
    let affine = gap(-5, -1);
    let aligner = XDrop::new(scoring, affine, 3).expect("valid configuration");
    let empty = DnaSequence::from_ascii(b"").expect("empty DNA is valid");
    let identity = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");

    let empty_result = aligner
        .align(&empty, &empty)
        .expect("empty alignment succeeds");
    assert_eq!(empty_result.score(), 0);
    assert_eq!(empty_result.first_end(), 0);
    assert_eq!(empty_result.second_end(), 0);
    assert!(empty_result.alignment().operations().is_empty());
    assert_eq!(empty_result.summary().evaluated_cells(), 1);
    assert_eq!(
        empty_result.termination(),
        XDropTermination::MatrixExhausted
    );

    let identity_result = aligner
        .align(&identity, &identity)
        .expect("identity alignment succeeds");
    assert_eq!(identity_result.score(), 8);
    assert_eq!(identity_result.first_end(), 4);
    assert_eq!(identity_result.second_end(), 4);
    assert_eq!(
        identity_result.alignment().operations(),
        &[AlignmentOperation::Match; 4]
    );
    assert_eq!(
        replay(
            &identity,
            &identity,
            identity_result.alignment(),
            &scoring,
            affine,
        ),
        8
    );

    let first = DnaSequence::from_ascii(b"AAAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TTTT").expect("valid DNA");
    let difference = XDrop::new(scoring, affine, 2)
        .expect("valid configuration")
        .align(&first, &second)
        .expect("different sequences are valid");
    assert_eq!(difference.score(), 0);
    assert!(difference.alignment().operations().is_empty());
    assert_eq!(
        difference.termination(),
        XDropTermination::FrontiersExhausted
    );

    let repeated_first = DnaSequence::from_ascii(b"AAAAAA").expect("valid DNA");
    let repeated_second = DnaSequence::from_ascii(b"AAAA").expect("valid DNA");
    let repeated = aligner
        .align(&repeated_first, &repeated_second)
        .expect("repetitive alignment succeeds");
    assert_eq!(repeated.score(), 8);
    assert_eq!(repeated.first_end(), 4);
    assert_eq!(repeated.second_end(), 4);
    replay(
        &repeated_first,
        &repeated_second,
        repeated.alignment(),
        &scoring,
        affine,
    );
}

#[test]
fn strict_threshold_equality_crosses_a_score_valley() {
    let first = DnaSequence::from_ascii(b"AAATAAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AAACAAA").expect("valid DNA");
    let scoring = MatchMismatch::new(3, -3);
    let affine = gap(-6, -1);

    let below = XDrop::new(scoring, affine, 2)
        .expect("valid configuration")
        .score(&first, &second)
        .expect("X-drop succeeds");
    assert_eq!(below.score(), 9);
    assert_eq!((below.first_end(), below.second_end()), (3, 3));

    let equality = XDrop::new(scoring, affine, 3)
        .expect("valid configuration")
        .align(&first, &second)
        .expect("threshold-equality path survives");
    assert_eq!(equality.score(), 15);
    assert_eq!((equality.first_end(), equality.second_end()), (7, 7));
    assert_eq!(
        equality.alignment().operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Substitution,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );
    replay(&first, &second, equality.alignment(), &scoring, affine);
}

#[test]
fn affine_internal_and_zero_cost_terminal_gaps_trace_correctly() {
    let first = DnaSequence::from_ascii(b"ACGGT").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
    let scoring = MatchMismatch::new(5, -4);
    let affine = gap(-2, -1);
    let result = XDrop::new(scoring, affine, 20)
        .expect("valid configuration")
        .align(&first, &second)
        .expect("affine X-drop succeeds");
    assert_eq!(result.score(), 18);
    let rendered = result
        .alignment()
        .render(&first, &second)
        .expect("alignment renders");
    assert_eq!(rendered.first(), "ACGGT");
    assert_eq!(rendered.second(), "ACG-T");
    replay(&first, &second, result.alignment(), &scoring, affine);

    let empty = DnaSequence::from_ascii(b"").expect("empty DNA is valid");
    let terminal = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let zero_gap = gap(0, 0);
    let terminal_result = XDrop::new(scoring, zero_gap, 0)
        .expect("valid configuration")
        .align(&empty, &terminal)
        .expect("zero-cost terminal gap succeeds");
    assert_eq!(terminal_result.score(), 0);
    assert_eq!(terminal_result.second_end(), 3);
    assert_eq!(
        terminal_result.alignment().operations(),
        &[AlignmentOperation::Insertion; 3]
    );
    replay(
        &empty,
        &terminal,
        terminal_result.alignment(),
        &scoring,
        zero_gap,
    );
}

#[test]
fn second_prior_diagonal_and_published_ties_are_deterministic() {
    let sequence = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let scoring = MatchMismatch::new(1, -1);
    let result = XDrop::new(scoring, gap(-2, -1), 0)
        .expect("valid configuration")
        .align(&sequence, &sequence)
        .expect("diagonal successors survive empty intervening frontier");
    assert_eq!(result.score(), 2);
    assert_eq!((result.first_end(), result.second_end()), (2, 2));

    let one = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let zero_scoring = MatchMismatch::new(0, 0);
    let tied = XDrop::new(zero_scoring, gap(0, 0), 0)
        .expect("valid configuration")
        .align(&one, &one)
        .expect("tied alignment succeeds");
    assert_eq!((tied.first_end(), tied.second_end()), (1, 1));
    assert_eq!(
        tied.alignment().operations(),
        &[AlignmentOperation::Deletion, AlignmentOperation::Insertion]
    );
    assert_eq!(
        XDrop::new(zero_scoring, gap(0, 0), 0)
            .expect("valid configuration")
            .align(&one, &one)
            .expect("repeat succeeds")
            .alignment(),
        tied.alignment()
    );
}

#[test]
fn dna_rna_and_protein_share_the_same_generic_api() {
    let scoring = MatchMismatch::new(4, -3);
    let affine = gap(-7, -1);
    let aligner = XDrop::new(scoring, affine, 8).expect("valid configuration");
    let dna = DnaSequence::from_ascii(b"ACGTN").expect("valid DNA");
    let rna = RnaSequence::from_ascii(b"ACGUN").expect("valid RNA");
    let protein = ProteinSequence::from_ascii(b"MKWVTF").expect("valid protein");

    let dna_result = aligner.align(&dna, &dna).expect("DNA alignment succeeds");
    let rna_result = aligner.align(&rna, &rna).expect("RNA alignment succeeds");
    let protein_result = aligner
        .align(&protein, &protein)
        .expect("protein alignment succeeds");
    assert_eq!(dna_result.score(), 20);
    assert_eq!(rna_result.score(), 20);
    assert_eq!(protein_result.score(), 24);
    replay(&dna, &dna, dna_result.alignment(), &scoring, affine);
    replay(&rna, &rna, rna_result.alignment(), &scoring, affine);
    replay(
        &protein,
        &protein,
        protein_result.alignment(),
        &scoring,
        affine,
    );
}

#[test]
fn forced_seed_supports_all_directions_and_counts_seed_once() {
    let first = DnaSequence::from_ascii(b"TTACGTAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"GGACGTCC").expect("valid DNA");
    let scoring = MatchMismatch::new(2, -3);
    let affine = gap(-6, -1);
    let aligner = XDrop::new(scoring, affine, 2).expect("valid configuration");
    let seed = XDropSeed::new(3, 3);

    let both = aligner
        .align_from_seed(&first, &second, seed, XDropDirection::Both)
        .expect("two-way seeded extension succeeds");
    assert_eq!(both.score(), 8);
    assert_eq!(both.first_range(), 2..6);
    assert_eq!(both.second_range(), 2..6);
    assert_eq!(
        both.alignment().operations(),
        &[AlignmentOperation::Match; 4]
    );
    assert!(both.summary().upstream().is_some());
    assert!(both.summary().downstream().is_some());
    replay(&first, &second, both.alignment(), &scoring, affine);

    let upstream = aligner
        .align_from_seed(&first, &second, seed, XDropDirection::Upstream)
        .expect("upstream extension succeeds");
    assert_eq!(upstream.score(), 4);
    assert_eq!(upstream.first_range(), 2..4);
    assert_eq!(upstream.second_range(), 2..4);
    assert!(upstream.summary().upstream().is_some());
    assert!(upstream.summary().downstream().is_none());

    let downstream = aligner
        .align_from_seed(&first, &second, seed, XDropDirection::Downstream)
        .expect("downstream extension succeeds");
    assert_eq!(downstream.score(), 6);
    assert_eq!(downstream.first_range(), 3..6);
    assert_eq!(downstream.second_range(), 3..6);
    assert!(downstream.summary().upstream().is_none());
    assert!(downstream.summary().downstream().is_some());

    let score_only = aligner
        .score_from_seed(&first, &second, seed, XDropDirection::Both)
        .expect("seeded score succeeds");
    assert_eq!(score_only, both.summary());
    assert_eq!(
        score_only.evaluated_cells(),
        score_only
            .upstream()
            .expect("upstream exists")
            .evaluated_cells()
            + score_only
                .downstream()
                .expect("downstream exists")
                .evaluated_cells()
            - 1
    );
}

#[test]
fn upstream_traceback_reverses_an_internal_gap_without_changing_its_side() {
    let first = DnaSequence::from_ascii(b"ACGGT").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
    let scoring = MatchMismatch::new(5, -4);
    let affine = gap(-2, -1);
    let result = XDrop::new(scoring, affine, 20)
        .expect("valid configuration")
        .align_from_seed(
            &first,
            &second,
            XDropSeed::new(4, 3),
            XDropDirection::Upstream,
        )
        .expect("upstream seeded alignment succeeds");
    assert_eq!(result.score(), 18);
    assert_eq!(result.first_range(), 0..5);
    assert_eq!(result.second_range(), 0..4);
    let rendered = result
        .alignment()
        .render(&first, &second)
        .expect("alignment renders");
    assert_eq!(rendered.first(), "ACGGT");
    assert_eq!(rendered.second(), "AC-GT");
    replay(&first, &second, result.alignment(), &scoring, affine);
}

#[test]
fn a_mismatching_seed_is_forced_even_when_its_score_is_negative() {
    let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"T").expect("valid DNA");
    let scoring = MatchMismatch::new(2, -3);
    let affine = gap(-5, -1);
    let result = XDrop::new(scoring, affine, 5)
        .expect("valid configuration")
        .align_from_seed(&first, &second, XDropSeed::new(0, 0), XDropDirection::Both)
        .expect("forced mismatch succeeds");
    assert_eq!(result.score(), -3);
    assert_eq!(result.first_range(), 0..1);
    assert_eq!(result.second_range(), 0..1);
    assert_eq!(
        result.alignment().operations(),
        &[AlignmentOperation::Substitution]
    );
    assert_eq!(result.summary().evaluated_cells(), 1);
    assert_eq!(result.summary().retained_cells(), 1);
    replay(&first, &second, result.alignment(), &scoring, affine);
}

#[test]
fn exhaustive_short_inputs_match_the_dense_oracle_and_traceback_replay() {
    let sequences = binary_dna_strings(4);
    let configurations = [
        (MatchMismatch::new(2, -3), gap(-5, -1), 0),
        (MatchMismatch::new(2, -3), gap(-5, -1), 3),
        (MatchMismatch::new(1, 0), gap(-2, 0), 2),
        (MatchMismatch::new(0, 0), gap(0, 0), 0),
    ];

    for first in &sequences {
        for second in &sequences {
            for &(scoring, affine, threshold) in &configurations {
                let aligner =
                    XDrop::new(scoring, affine, threshold).expect("configuration is valid");
                let expected = oracle(first, second, &scoring, affine, threshold);
                let observed = aligner.score(first, second).expect("score succeeds");
                assert_matches_oracle(observed, expected);
                let alignment = aligner.align(first, second).expect("traceback succeeds");
                assert_eq!(alignment.summary(), observed);
                replay(first, second, alignment.alignment(), &scoring, affine);
                assert_eq!(
                    aligner.align(first, second).expect("repeat succeeds"),
                    alignment
                );
            }
        }
    }
}

#[test]
fn large_threshold_matches_unrestricted_matrix_and_all_tiny_paths() {
    let sequences = binary_dna_strings(3);
    let configurations = [
        (MatchMismatch::new(3, -2), gap(-4, -1)),
        (MatchMismatch::new(0, 0), gap(0, 0)),
        (MatchMismatch::new(-1, -2), gap(-3, -1)),
    ];
    for first in &sequences {
        for second in &sequences {
            for &(scoring, affine) in &configurations {
                let matrix = unrestricted_oracle(first, second, &scoring, affine);
                let paths = exhaustive_path_oracle(first, second, &scoring, affine);
                assert_eq!(matrix, paths);
                let observed = XDrop::new(scoring, affine, 1_000)
                    .expect("valid configuration")
                    .score(first, second)
                    .expect("score succeeds");
                assert_eq!(observed.score(), matrix.score);
                assert_eq!(observed.first_end(), matrix.first_end);
                assert_eq!(observed.second_end(), matrix.second_end);
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct MaximumScorer;

impl SubstitutionScorer<Dna> for MaximumScorer {
    fn score(&self, _first: Symbol<Dna>, _second: Symbol<Dna>) -> i64 {
        i64::MAX
    }
}

#[test]
fn configuration_seed_dimension_and_score_failures_are_typed() {
    assert!(matches!(
        XDrop::new(MatchMismatch::new(1, -1), gap(-2, -1), -1),
        Err(XDropError::NegativeThreshold { threshold: -1 })
    ));
    let one = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let aligner =
        XDrop::new(MatchMismatch::new(1, -1), gap(-2, -1), 3).expect("valid configuration");
    assert!(matches!(
        aligner.score_from_seed(&one, &one, XDropSeed::new(1, 0), XDropDirection::Both),
        Err(XDropError::SeedOutOfBounds {
            first_index: 1,
            second_index: 0,
            first_len: 1,
            second_len: 1,
        })
    ));
    assert!(matches!(
        XDrop::<MatchMismatch>::memory_requirement(usize::MAX, 1),
        Err(XDropError::DimensionOverflow { .. })
    ));

    let two = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    assert!(matches!(
        XDrop::new(MaximumScorer, gap(-2, -1), i64::MAX)
            .expect("valid configuration")
            .score(&two, &two),
        Err(XDropError::ScoreOverflow {
            first_prefix: 2,
            second_prefix: 2,
        })
    ));
}
