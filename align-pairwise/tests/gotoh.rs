//! Manual and property tests for scalar global Gotoh alignment.

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, DenseSubstitutionMatrix, Dna, DnaSequence,
    MatchMismatch, ProteinSequence, RnaSequence, Score, Sequence, SubstitutionScorer, Symbol,
};
use align_pairwise::{AlignmentError, Gotoh, NeedlemanWunsch, verify_global_affine_alignment};
use proptest::prelude::*;

fn dna(bytes: &[u8]) -> DnaSequence {
    DnaSequence::from_ascii(bytes).expect("test DNA must be valid")
}

#[derive(Clone)]
struct RowCachedMatrix<A: Alphabet>(DenseSubstitutionMatrix<A>);

impl<A: Alphabet> SubstitutionScorer<A> for RowCachedMatrix<A> {
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score {
        self.0.score(first, second)
    }

    fn fill_score_row(&self, first: Symbol<A>, output: &mut [Score]) -> bool {
        self.0.copy_score_row(first, output)
    }

    fn score_bounds(&self) -> Option<(Score, Score)> {
        self.0.score_bounds()
    }
}

fn standard_aligner() -> Gotoh<MatchMismatch> {
    Gotoh::new(
        MatchMismatch::new(2, -3),
        AffineGap::new(-5, -1).expect("valid affine gap"),
    )
}

#[test]
fn empty_inputs_have_the_unique_empty_global_alignment() {
    let result = standard_aligner()
        .align(&dna(b""), &dna(b""))
        .expect("empty alignment succeeds");
    assert_eq!(result.score(), 0);
    assert_eq!(result.first_range(), 0..0);
    assert_eq!(result.second_range(), 0..0);
    assert!(result.operations().is_empty());
}

#[test]
fn empty_boundaries_charge_one_open_and_only_extensions_afterward() {
    let empty = dna(b"");
    let sequence = dna(b"ACGT");
    let aligner = standard_aligner();

    let insertion = aligner
        .align(&empty, &sequence)
        .expect("empty first sequence succeeds");
    assert_eq!(insertion.score(), -8);
    assert_eq!(insertion.operations(), &[AlignmentOperation::Insertion; 4]);

    let deletion = aligner
        .align(&sequence, &empty)
        .expect("empty second sequence succeeds");
    assert_eq!(deletion.score(), -8);
    assert_eq!(deletion.operations(), &[AlignmentOperation::Deletion; 4]);

    verify_global_affine_alignment(
        &empty,
        &sequence,
        &insertion,
        aligner.substitution(),
        aligner.gap(),
    )
    .expect("insertion boundary score replays");
    verify_global_affine_alignment(
        &sequence,
        &empty,
        &deletion,
        aligner.substitution(),
        aligner.gap(),
    )
    .expect("deletion boundary score replays");
}

#[test]
fn identical_dna_uses_only_matches() {
    let sequence = dna(b"ACGTN");
    let result = standard_aligner()
        .align(&sequence, &sequence)
        .expect("alignment succeeds");
    assert_eq!(result.score(), 10);
    assert_eq!(result.operations(), &[AlignmentOperation::Match; 5]);
}

#[test]
fn completely_different_dna_prefers_substitutions_for_this_model() {
    let first = dna(b"AAA");
    let second = dna(b"CCC");
    let result = standard_aligner()
        .align(&first, &second)
        .expect("alignment succeeds");
    assert_eq!(result.score(), -9);
    assert_eq!(result.operations(), &[AlignmentOperation::Substitution; 3]);
}

#[test]
fn long_internal_gap_is_opened_once() {
    let first = dna(b"AAACCCCBBB");
    let second = dna(b"AAABBB");
    let aligner = Gotoh::new(
        MatchMismatch::new(2, -4),
        AffineGap::new(-3, -1).expect("valid affine gap"),
    );
    let result = aligner.align(&first, &second).expect("alignment succeeds");
    assert_eq!(result.score(), 6);
    assert_eq!(
        result.operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
            AlignmentOperation::Deletion,
            AlignmentOperation::Deletion,
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );
    let rendered = result.render(&first, &second).expect("render succeeds");
    assert_eq!(rendered.first(), "AAACCCCBBB");
    assert_eq!(rendered.second(), "AAA----BBB");
}

#[test]
fn initial_and_terminal_gap_runs_are_penalized() {
    let aligner = Gotoh::new(
        MatchMismatch::new(3, -4),
        AffineGap::new(-4, -1).expect("valid affine gap"),
    );
    let initial = aligner
        .align(&dna(b"GAC"), &dna(b"AC"))
        .expect("initial gap succeeds");
    assert_eq!(initial.score(), 2);
    assert_eq!(
        initial.operations(),
        &[
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );

    let terminal = aligner
        .align(&dna(b"ACG"), &dna(b"AC"))
        .expect("terminal gap succeeds");
    assert_eq!(terminal.score(), 2);
    assert_eq!(
        terminal.operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
        ]
    );
}

#[test]
fn rna_and_protein_use_the_same_affine_architecture() {
    let rna_first = RnaSequence::from_ascii(b"AUGCU").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"AUCU").expect("valid RNA");
    let rna_aligner = Gotoh::new(
        MatchMismatch::new(2, -3),
        AffineGap::new(-3, -1).expect("valid gap"),
    );
    let rna_result = rna_aligner
        .align(&rna_first, &rna_second)
        .expect("RNA alignment succeeds");
    assert_eq!(rna_result.score(), 5);
    assert_eq!(
        rna_result.operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );

    let protein_first = ProteinSequence::from_ascii(b"MKTLL").expect("valid protein");
    let protein_second = ProteinSequence::from_ascii(b"MKLL").expect("valid protein");
    let protein_aligner = Gotoh::new(
        MatchMismatch::new(3, -4),
        AffineGap::new(-5, -1).expect("valid gap"),
    );
    let protein_result = protein_aligner
        .align(&protein_first, &protein_second)
        .expect("protein alignment succeeds");
    assert_eq!(protein_result.score(), 7);
    assert_eq!(protein_result.operations()[2], AlignmentOperation::Deletion);
}

#[test]
fn dense_matrix_controls_ambiguity_symbol_scoring() {
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        if first == second {
            4
        } else if first == b'N' || second == b'N' {
            1
        } else {
            -5
        }
    });
    let aligner = Gotoh::new(matrix, AffineGap::new(-6, -1).expect("valid affine gap"));
    let result = aligner
        .align(&dna(b"AN"), &dna(b"AG"))
        .expect("matrix alignment succeeds");
    assert_eq!(result.score(), 5);
    assert_eq!(
        result.operations(),
        &[AlignmentOperation::Match, AlignmentOperation::Substitution]
    );
}

#[test]
fn dense_score_row_fast_path_matches_scalar_matrix_lookups() {
    let first = dna(b"ACGTRYSWKMBDHVNACGT");
    let second = dna(b"NVDHBMKWSRYTGCAN");
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        i64::from(first) * 3 - i64::from(second) * 2
    });
    let gap = AffineGap::new(-9, -2).expect("valid gap");
    let cached = Gotoh::new(RowCachedMatrix(matrix.clone()), gap)
        .align(&first, &second)
        .expect("cached dense alignment succeeds");
    let scalar = Gotoh::new(matrix, gap)
        .align(&first, &second)
        .expect("scalar dense alignment succeeds");
    assert_eq!(cached, scalar);
}

#[test]
fn diagonal_wins_a_three_way_score_tie() {
    let aligner = Gotoh::new(
        MatchMismatch::new(2, -2),
        AffineGap::new(-1, -1).expect("valid affine gap"),
    );
    let result = aligner
        .align(&dna(b"A"), &dna(b"C"))
        .expect("tie alignment succeeds");
    assert_eq!(result.score(), -2);
    assert_eq!(result.operations(), &[AlignmentOperation::Substitution]);
}

#[test]
fn repetitive_input_has_a_deterministic_canonical_optimum() {
    let aligner = Gotoh::new(
        MatchMismatch::new(2, -3),
        AffineGap::new(-2, -1).expect("valid affine gap"),
    );
    let first = dna(b"AAA");
    let second = dna(b"AA");
    let expected = aligner.align(&first, &second).expect("alignment succeeds");
    assert_eq!(expected.score(), 2);
    assert_eq!(
        expected.operations(),
        &[
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ]
    );
    for _ in 0..32 {
        assert_eq!(
            aligner.align(&first, &second).expect("repeat succeeds"),
            expected
        );
    }
}

#[test]
fn published_recurrence_allows_adjacent_opposite_gap_runs() {
    let aligner = Gotoh::new(
        MatchMismatch::new(1, -100),
        AffineGap::new(-2, -1).expect("valid affine gap"),
    );
    let first = dna(b"A");
    let second = dna(b"CCC");
    let result = aligner.align(&first, &second).expect("alignment succeeds");
    assert_eq!(result.score(), -6);
    assert_eq!(
        result.operations(),
        &[
            AlignmentOperation::Insertion,
            AlignmentOperation::Insertion,
            AlignmentOperation::Insertion,
            AlignmentOperation::Deletion,
        ]
    );
    verify_global_affine_alignment(
        &first,
        &second,
        &result,
        aligner.substitution(),
        aligner.gap(),
    )
    .expect("adjacent runs replay independently");
}

#[test]
fn exact_score_minimum_is_valid_but_a_later_underflow_is_reported() {
    let one = dna(b"A");
    let two = dna(b"AA");
    let exact_minimum = Gotoh::new(
        MatchMismatch::new(1, -1),
        AffineGap::new(Score::MIN, -1).expect("valid extreme gap"),
    );
    assert_eq!(
        exact_minimum
            .align(&dna(b""), &one)
            .expect("minimum score is representable")
            .score(),
        Score::MIN
    );

    let underflow = Gotoh::new(
        MatchMismatch::new(1, -1),
        AffineGap::new(Score::MIN, -1).expect("valid extreme gap"),
    )
    .align(&dna(b""), &two);
    assert!(matches!(
        underflow,
        Err(AlignmentError::ScoreOverflow { row: 0, column: 2 })
    ));
}

#[test]
fn checked_path_ignores_underflowed_auxiliary_states_when_best_is_representable() {
    let aligner = Gotoh::new(
        MatchMismatch::new(1, -1),
        AffineGap::new(Score::MIN, 0).expect("valid extreme gap"),
    );
    let sequence = dna(b"AA");
    let result = aligner
        .align(&sequence, &sequence)
        .expect("representable diagonal optimum succeeds");
    assert_eq!(result.score(), 2);
    assert_eq!(result.operations(), &[AlignmentOperation::Match; 2]);
}

#[test]
fn positive_optimal_score_overflow_is_reported() {
    let aligner = Gotoh::new(
        MatchMismatch::new(Score::MAX, -1),
        AffineGap::new(-5, -1).expect("valid affine gap"),
    );
    let sequence = dna(b"AA");
    assert!(matches!(
        aligner.align(&sequence, &sequence),
        Err(AlignmentError::ScoreOverflow { row: 2, column: 2 })
    ));
}

#[test]
fn negatively_underflowing_losing_diagonal_does_not_abort_alignment() {
    let aligner = Gotoh::new(
        MatchMismatch::new(1, Score::MIN),
        AffineGap::new(-5, -1).expect("valid affine gap"),
    );
    let result = aligner
        .align(&dna(b"A"), &dna(b"C"))
        .expect("representable two-gap optimum succeeds");
    assert_eq!(result.score(), -10);
    assert_eq!(
        result.operations(),
        &[AlignmentOperation::Insertion, AlignmentOperation::Deletion,]
    );
}

#[derive(Clone, Copy)]
struct UnboundedMatchMismatch;

impl SubstitutionScorer<Dna> for UnboundedMatchMismatch {
    fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
        if first == second { 3 } else { -4 }
    }
}

#[test]
fn scorer_without_declared_bounds_uses_checked_arithmetic() {
    let aligner = Gotoh::new(
        UnboundedMatchMismatch,
        AffineGap::new(-5, -1).expect("valid gap"),
    );
    let first = dna(b"ACGT");
    let second = dna(b"AGT");
    let result = aligner.align(&first, &second).expect("alignment succeeds");
    assert_eq!(result.score(), 4);
    verify_global_affine_alignment(
        &first,
        &second,
        &result,
        &UnboundedMatchMismatch,
        aligner.gap(),
    )
    .expect("checked-path result replays");
}

fn dna_from_indices(indices: &[u8]) -> Sequence<Dna> {
    let bytes: Vec<u8> = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index)])
        .collect();
    Sequence::<Dna>::from_ascii(&bytes).expect("generated DNA is valid")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(160))]

    #[test]
    fn affine_results_replay_are_symmetric_and_deterministic(
        first_indices in prop::collection::vec(0_u8..4, 0..9),
        second_indices in prop::collection::vec(0_u8..4, 0..9),
    ) {
        let first = dna_from_indices(&first_indices);
        let second = dna_from_indices(&second_indices);
        let aligner = standard_aligner();
        let forward = aligner.align(&first, &second).expect("generated alignment succeeds");
        let repeated = aligner.align(&first, &second).expect("repeated alignment succeeds");
        let reverse = aligner.align(&second, &first).expect("reverse alignment succeeds");

        prop_assert_eq!(&forward, &repeated);
        prop_assert_eq!(forward.first_range(), 0..first.len());
        prop_assert_eq!(forward.second_range(), 0..second.len());
        prop_assert_eq!(forward.score(), reverse.score());
        let report = verify_global_affine_alignment(
            &first,
            &second,
            &forward,
            aligner.substitution(),
            aligner.gap(),
        ).expect("generated result must replay");
        prop_assert_eq!(report.recomputed_score(), forward.score());
    }

    #[test]
    fn equal_open_and_extension_scores_reduce_to_linear_global_scores(
        first_indices in prop::collection::vec(0_u8..4, 0..9),
        second_indices in prop::collection::vec(0_u8..4, 0..9),
    ) {
        let first = dna_from_indices(&first_indices);
        let second = dna_from_indices(&second_indices);
        let substitution = MatchMismatch::new(2, -1);
        let affine = Gotoh::new(
            substitution,
            AffineGap::new(-2, -2).expect("linear affine specialization is valid"),
        );
        let linear = NeedlemanWunsch::new(
            substitution,
            align_core::LinearGap::new(-2).expect("valid linear gap"),
        );
        let affine_result = affine.align(&first, &second).expect("affine alignment succeeds");
        let linear_result = linear.align(&first, &second).expect("linear alignment succeeds");
        prop_assert_eq!(affine_result.score(), linear_result.score());
    }
}
