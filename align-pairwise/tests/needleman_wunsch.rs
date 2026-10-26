//! Manual edge cases and property-based tests for global linear-gap alignment.

use align_core::{
    AlignmentOperation, DenseSubstitutionMatrix, Dna, DnaSequence, LinearGap, MatchMismatch,
    ProteinSequence, RnaSequence,
};
use align_pairwise::{AlignmentError, NeedlemanWunsch, verify_global_alignment};
use proptest::prelude::*;

fn scoring() -> MatchMismatch {
    MatchMismatch::new(2, -1)
}

fn gap() -> LinearGap {
    LinearGap::new(-2).expect("test gap is valid")
}

#[test]
fn empty_sequences_produce_the_empty_alignment() {
    let first = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    let second = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    let aligner = NeedlemanWunsch::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 0);
    assert!(alignment.operations().is_empty());
    assert_eq!(alignment.first_range(), 0..0);
    assert_eq!(alignment.second_range(), 0..0);
    let report = verify_global_alignment(&first, &second, &alignment, &scoring(), gap())
        .expect("alignment verifies");
    assert_eq!(report.recomputed_score(), 0);
    assert_eq!(report.columns(), 0);
    assert_eq!(report.matches(), 0);
    assert_eq!(report.substitutions(), 0);
    assert_eq!(report.deletions(), 0);
    assert_eq!(report.insertions(), 0);
}

#[test]
fn empty_first_sequence_is_all_insertions() {
    let first = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    let second = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let aligner = NeedlemanWunsch::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), -6);
    assert_eq!(alignment.operations(), &[AlignmentOperation::Insertion; 3]);
    let rendered = alignment
        .render(&first, &second)
        .expect("alignment renders");
    assert_eq!(rendered.first(), "---");
    assert_eq!(rendered.second(), "ACG");
}

#[test]
fn empty_second_sequence_is_all_deletions() {
    let first = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    let aligner = NeedlemanWunsch::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), -6);
    assert_eq!(alignment.operations(), &[AlignmentOperation::Deletion; 3]);
}

#[test]
fn identical_dna_rna_and_protein_sequences_are_supported() {
    let aligner = NeedlemanWunsch::new(scoring(), gap());

    let dna = DnaSequence::from_ascii(b"ACGTRYN").expect("valid ambiguous DNA");
    let dna_alignment = aligner.align(&dna, &dna).expect("DNA alignment succeeds");
    assert_eq!(dna_alignment.score(), 14);
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
}

#[test]
fn completely_different_sequences_prefer_substitutions_when_cheaper() {
    let first = DnaSequence::from_ascii(b"AAAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TTTT").expect("valid DNA");
    let aligner = NeedlemanWunsch::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), -4);
    assert_eq!(
        alignment.operations(),
        &[AlignmentOperation::Substitution; 4]
    );
}

#[test]
fn leading_internal_and_terminal_gaps_are_reconstructed() {
    let scorer = MatchMismatch::new(2, -3);
    let gap = LinearGap::new(-2).expect("valid gap");
    let aligner = NeedlemanWunsch::new(scorer, gap);

    let cases = [
        ("AC", "TAC", "-AC", "TAC"),
        ("AC", "AGC", "A-C", "AGC"),
        ("AC", "ACT", "AC-", "ACT"),
    ];
    for (first_raw, second_raw, expected_first, expected_second) in cases {
        let first = DnaSequence::from_ascii(first_raw).expect("valid DNA");
        let second = DnaSequence::from_ascii(second_raw).expect("valid DNA");
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        let rendered = alignment
            .render(&first, &second)
            .expect("alignment renders");
        assert_eq!(alignment.score(), 2);
        assert_eq!(rendered.first(), expected_first);
        assert_eq!(rendered.second(), expected_second);
        verify_global_alignment(&first, &second, &alignment, &scorer, gap)
            .expect("alignment verifies");
    }
}

#[test]
fn classic_gattaca_example_has_hand_calculated_score() {
    let first = ProteinSequence::from_ascii(b"GATTACA").expect("valid protein symbols");
    let second = ProteinSequence::from_ascii(b"GCATGCU").expect("valid protein symbols");
    let scorer = MatchMismatch::new(1, -1);
    let gap = LinearGap::new(-1).expect("valid gap");
    let aligner = NeedlemanWunsch::new(scorer, gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 0);
    verify_global_alignment(&first, &second, &alignment, &scorer, gap).expect("alignment verifies");
}

#[test]
fn repetitive_input_uses_documented_tie_breaking() {
    let first = DnaSequence::from_ascii(b"AAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let aligner = NeedlemanWunsch::new(scoring(), gap());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");
    let rendered = alignment
        .render(&first, &second)
        .expect("alignment renders");

    assert_eq!(alignment.score(), 2);
    assert_eq!(rendered.first(), "AAA");
    assert_eq!(rendered.second(), "-AA");
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
    let first = DnaSequence::from_ascii(b"AN").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let gap = LinearGap::new(-5).expect("valid gap");
    let aligner = NeedlemanWunsch::new(matrix.clone(), gap);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 4);
    verify_global_alignment(&first, &second, &alignment, &matrix, gap).expect("alignment verifies");
}

#[test]
fn score_overflow_is_reported() {
    let first = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let aligner = NeedlemanWunsch::new(
        MatchMismatch::new(i64::MAX, 0),
        LinearGap::new(0).expect("valid gap"),
    );
    let error = aligner
        .align(&first, &second)
        .expect_err("unrepresentable score must fail");
    assert!(matches!(error, AlignmentError::ScoreOverflow { .. }));
}

#[test]
fn unrepresentable_nonoptimal_candidate_does_not_abort_alignment() {
    let first = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TT").expect("valid DNA");
    let scorer = MatchMismatch::new(0, i64::MIN);
    let gap = LinearGap::new(-1).expect("valid gap");
    let aligner = NeedlemanWunsch::new(scorer, gap);
    let alignment = aligner
        .align(&first, &second)
        .expect("only a nonoptimal candidate is out of range");

    assert_eq!(alignment.score(), -4);
    verify_global_alignment(&first, &second, &alignment, &scorer, gap)
        .expect("selected optimal path remains representable");
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
    fn every_traceback_replays_and_rescores_independently(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = NeedlemanWunsch::new(scoring(), gap());
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        let report = verify_global_alignment(&first, &second, &alignment, &scoring(), gap())
            .expect("traceback and score verify");
        prop_assert_eq!(report.recomputed_score(), alignment.score());
        prop_assert_eq!(report.columns(), alignment.columns());

        let rendered = alignment.render(&first, &second).expect("alignment renders");
        prop_assert_eq!(rendered.first().replace('-', ""), first.to_ascii_string());
        prop_assert_eq!(rendered.second().replace('-', ""), second.to_ascii_string());
    }

    #[test]
    fn symmetric_scoring_produces_symmetric_optimal_scores(
        first_raw in dna_strategy(14),
        second_raw in dna_strategy(14),
    ) {
        let first = DnaSequence::from_ascii(&first_raw).expect("strategy emits valid DNA");
        let second = DnaSequence::from_ascii(&second_raw).expect("strategy emits valid DNA");
        let aligner = NeedlemanWunsch::new(scoring(), gap());
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
        let aligner = NeedlemanWunsch::new(scoring(), gap());
        let first_run = aligner.align(&first, &second).expect("alignment succeeds");
        let second_run = aligner.align(&first, &second).expect("alignment succeeds");
        prop_assert_eq!(first_run, second_run);
    }
}
