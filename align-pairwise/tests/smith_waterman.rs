//! Manual edge cases and property-based tests for local linear-gap alignment.

use align_core::{
    AlignmentOperation, DenseSubstitutionMatrix, Dna, DnaSequence, LinearGap, MatchMismatch,
    ProteinSequence, RnaSequence,
};
use align_pairwise::{AlignmentError, SmithWaterman, verify_local_alignment};
use proptest::prelude::*;

fn scoring() -> MatchMismatch {
    MatchMismatch::new(2, -1)
}

fn gap() -> LinearGap {
    LinearGap::new(-2).expect("test gap is valid")
}

#[test]
fn empty_inputs_produce_the_canonical_empty_alignment() {
    let empty = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    let nonempty = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let aligner = SmithWaterman::new(scoring(), gap());

    for (first, second) in [(&empty, &empty), (&empty, &nonempty), (&nonempty, &empty)] {
        let alignment = aligner.align(first, second).expect("alignment succeeds");
        assert_eq!(alignment.score(), 0);
        assert!(alignment.operations().is_empty());
        assert_eq!(alignment.first_range(), 0..0);
        assert_eq!(alignment.second_range(), 0..0);
        let report = verify_local_alignment(first, second, &alignment, &scoring(), gap())
            .expect("alignment verifies");
        assert_eq!(report.recomputed_score(), 0);
        assert_eq!(report.columns(), 0);
    }
}

#[test]
fn identical_dna_rna_and_protein_sequences_are_supported() {
    let aligner = SmithWaterman::new(scoring(), gap());

    let dna = DnaSequence::from_ascii(b"ACGTRYN").expect("valid ambiguous DNA");
    let dna_alignment = aligner.align(&dna, &dna).expect("DNA alignment succeeds");
    assert_eq!(dna_alignment.score(), 14);
    assert_eq!(dna_alignment.first_range(), 0..dna.len());
    assert!(
        dna_alignment
            .operations()
            .iter()
            .all(|operation| *operation == AlignmentOperation::Match)
    );

    let rna = RnaSequence::from_ascii(b"ACGURYN").expect("valid ambiguous RNA");
    let rna_alignment = aligner.align(&rna, &rna).expect("RNA alignment succeeds");
    assert_eq!(rna_alignment.score(), 14);

    let protein =
        ProteinSequence::from_ascii(b"ARNDCQEGHILKMFPSTWYVBZXJUO").expect("valid protein");
    let protein_alignment = aligner
        .align(&protein, &protein)
        .expect("protein alignment succeeds");
    assert_eq!(protein_alignment.score(), 52);
    verify_local_alignment(&protein, &protein, &protein_alignment, &scoring(), gap())
        .expect("protein alignment verifies");
}

#[test]
fn completely_different_sequences_return_the_canonical_empty_alignment() {
    let first = DnaSequence::from_ascii(b"AAAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TTTT").expect("valid DNA");
    let aligner = SmithWaterman::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 0);
    assert!(alignment.operations().is_empty());
    assert_eq!(alignment.first_range(), 0..0);
    assert_eq!(alignment.second_range(), 0..0);
}

#[test]
fn zero_score_candidate_terminates_instead_of_extending_a_zero_path() {
    let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"T").expect("valid DNA");
    let scorer = MatchMismatch::new(1, 0);
    let gap = LinearGap::new(0).expect("zero gap is valid");
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 0);
    assert!(alignment.operations().is_empty());
    assert_eq!(alignment.first_range(), 0..0);
    assert_eq!(alignment.second_range(), 0..0);
}

#[test]
fn nonhomologous_prefixes_and_suffixes_are_excluded() {
    let first = DnaSequence::from_ascii(b"GGACAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TTACCC").expect("valid DNA");
    let scorer = MatchMismatch::new(2, -3);
    let gap = LinearGap::new(-2).expect("valid gap");
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");
    let rendered = alignment
        .render(&first, &second)
        .expect("alignment renders");

    assert_eq!(alignment.score(), 4);
    assert_eq!(alignment.first_range(), 2..4);
    assert_eq!(alignment.second_range(), 2..4);
    assert_eq!(rendered.first(), "AC");
    assert_eq!(rendered.second(), "AC");
    verify_local_alignment(&first, &second, &alignment, &scorer, gap).expect("alignment verifies");
}

#[test]
fn an_internal_gap_can_be_part_of_the_optimal_local_alignment() {
    let first = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AGT").expect("valid DNA");
    let scorer = MatchMismatch::new(2, -3);
    let gap = LinearGap::new(-1).expect("valid gap");
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");
    let rendered = alignment
        .render(&first, &second)
        .expect("alignment renders");

    assert_eq!(alignment.score(), 5);
    assert_eq!(alignment.first_range(), 0..4);
    assert_eq!(alignment.second_range(), 0..3);
    assert_eq!(rendered.first(), "ACGT");
    assert_eq!(rendered.second(), "A-GT");
    assert_eq!(
        alignment.operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );
}

#[test]
fn first_row_major_maximum_is_retained_for_repetitive_input() {
    let first = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let aligner = SmithWaterman::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 2);
    assert_eq!(alignment.first_range(), 0..1);
    assert_eq!(alignment.second_range(), 0..1);
    assert_eq!(alignment.operations(), &[AlignmentOperation::Match]);
}

#[test]
fn positive_predecessor_ties_prefer_diagonal() {
    let first = DnaSequence::from_ascii(b"AAB").expect("valid ambiguous DNA");
    let second = DnaSequence::from_ascii(b"AB").expect("valid ambiguous DNA");
    let scorer = MatchMismatch::new(1, 0);
    let gap = LinearGap::new(0).expect("zero gap is valid");
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");
    let rendered = alignment
        .render(&first, &second)
        .expect("alignment renders");

    assert_eq!(alignment.score(), 2);
    assert_eq!(alignment.first_range(), 1..3);
    assert_eq!(alignment.second_range(), 0..2);
    assert_eq!(rendered.first(), "AB");
    assert_eq!(rendered.second(), "AB");
}

#[test]
fn custom_dense_matrix_is_used_for_ambiguous_symbols() {
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        if first == second {
            3
        } else if first == b'N' || second == b'N' {
            1
        } else {
            -4
        }
    });
    let first = DnaSequence::from_ascii(b"TAN").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"GAA").expect("valid DNA");
    let gap = LinearGap::new(-5).expect("valid gap");
    let aligner = SmithWaterman::new(matrix.clone(), gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 4);
    assert_eq!(alignment.first_range(), 1..3);
    assert_eq!(alignment.second_range(), 1..3);
    verify_local_alignment(&first, &second, &alignment, &matrix, gap).expect("alignment verifies");
}

#[test]
fn score_overflow_is_reported() {
    let first = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let aligner = SmithWaterman::new(
        MatchMismatch::new(i64::MAX, 0),
        LinearGap::new(0).expect("valid gap"),
    );
    let error = aligner
        .align(&first, &second)
        .expect_err("unrepresentable score must fail");
    assert!(matches!(error, AlignmentError::ScoreOverflow { .. }));
}

#[test]
fn extreme_negative_candidate_resets_to_zero_without_overflow() {
    let first = DnaSequence::from_ascii(b"TA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"CA").expect("valid DNA");
    let scorer = MatchMismatch::new(2, i64::MIN);
    let gap = LinearGap::new(-1).expect("valid gap");
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner
        .align(&first, &second)
        .expect("negative candidate remains representable and nonoptimal");

    assert_eq!(alignment.score(), 2);
    assert_eq!(alignment.first_range(), 1..2);
    assert_eq!(alignment.second_range(), 1..2);
    verify_local_alignment(&first, &second, &alignment, &scorer, gap)
        .expect("selected local path verifies");
}

fn dna_strategy(max_length: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T'), Just(b'N')],
        0..=max_length,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(160))]

    #[test]
    fn every_traceback_replays_and_reconstructs_its_selected_ranges(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = SmithWaterman::new(scoring(), gap());
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        let report = verify_local_alignment(&first, &second, &alignment, &scoring(), gap())
            .expect("traceback and score verify");
        prop_assert_eq!(report.recomputed_score(), alignment.score());
        prop_assert_eq!(report.columns(), alignment.columns());
        prop_assert!(alignment.score() >= 0);

        let first_range = alignment.first_range();
        let second_range = alignment.second_range();
        prop_assert!(first_range.end <= first.len());
        prop_assert!(second_range.end <= second.len());
        let rendered = alignment.render(&first, &second).expect("alignment renders");
        let first_ascii = first.to_ascii_string();
        let second_ascii = second.to_ascii_string();
        prop_assert_eq!(rendered.first().replace('-', ""), &first_ascii[first_range]);
        prop_assert_eq!(rendered.second().replace('-', ""), &second_ascii[second_range]);
    }

    #[test]
    fn symmetric_scoring_produces_symmetric_optimal_scores(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = SmithWaterman::new(scoring(), gap());
        let forward = aligner.align(&first, &second).expect("alignment succeeds");
        let reverse = aligner.align(&second, &first).expect("alignment succeeds");
        prop_assert_eq!(forward.score(), reverse.score());
    }

    #[test]
    fn repeated_calls_are_bit_for_bit_deterministic(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = SmithWaterman::new(scoring(), gap());
        let first_run = aligner.align(&first, &second).expect("alignment succeeds");
        let second_run = aligner.align(&first, &second).expect("alignment succeeds");
        prop_assert_eq!(first_run, second_run);
    }
}
