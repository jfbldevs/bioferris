//! Manual acceptance tests for affine local and explicit four-end alignment.

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, DenseSubstitutionMatrix, Dna, DnaSequence,
    MatchMismatch, PairwiseAlignment, ProteinSequence, RnaSequence, Score, SubstitutionScorer,
    Symbol,
};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, AlignmentError, FreeEndGaps, Gotoh, VerificationError,
    verify_affine_end_gap_alignment, verify_local_affine_alignment,
};

fn dna(raw: &[u8]) -> DnaSequence {
    DnaSequence::from_ascii(raw).expect("test DNA is valid")
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

fn scoring() -> MatchMismatch {
    MatchMismatch::new(3, -4)
}

fn gap() -> AffineGap {
    AffineGap::new(-5, -1).expect("test affine gap is valid")
}

#[test]
fn global_mode_is_identical_to_gotoh_for_manual_edge_cases() {
    let cases = [
        (b"".as_slice(), b"".as_slice()),
        (b"".as_slice(), b"ACGT".as_slice()),
        (b"ACGT".as_slice(), b"".as_slice()),
        (b"ACGT".as_slice(), b"ACGT".as_slice()),
        (b"AAAA".as_slice(), b"TTTT".as_slice()),
        (b"ACGGGT".as_slice(), b"ACT".as_slice()),
        (b"AAAAAA".as_slice(), b"AAA".as_slice()),
    ];
    let affine_modes = AffineEndGapAligner::new(scoring(), gap(), FreeEndGaps::global());
    let gotoh = Gotoh::new(scoring(), gap());

    for (first_raw, second_raw) in cases {
        let first = dna(first_raw);
        let second = dna(second_raw);
        let observed = affine_modes
            .align(&first, &second)
            .expect("global affine mode succeeds");
        let expected = gotoh
            .align(&first, &second)
            .expect("Gotoh alignment succeeds");
        assert_eq!(observed, expected, "first={first}, second={second}");
    }
}

#[test]
fn affine_boundaries_charge_one_open_then_extensions() {
    let empty = dna(b"");
    let sequence = dna(b"ACGT");
    let aligner = AffineEndGapAligner::new(scoring(), gap(), FreeEndGaps::global());

    let insertions = aligner
        .align(&empty, &sequence)
        .expect("insertion boundary succeeds");
    assert_eq!(insertions.score(), -8);
    assert_eq!(insertions.operations(), &[AlignmentOperation::Insertion; 4]);

    let deletions = aligner
        .align(&sequence, &empty)
        .expect("deletion boundary succeeds");
    assert_eq!(deletions.score(), -8);
    assert_eq!(deletions.operations(), &[AlignmentOperation::Deletion; 4]);
}

#[test]
fn internal_leading_and_trailing_gap_runs_have_manual_scores() {
    let strong_mismatch = MatchMismatch::new(3, -10);
    let aligner = AffineEndGapAligner::new(strong_mismatch, gap(), FreeEndGaps::global());
    let cases = [
        (
            b"ACGGGT".as_slice(),
            b"ACT".as_slice(),
            2,
            vec![
                AlignmentOperation::Match,
                AlignmentOperation::Match,
                AlignmentOperation::Deletion,
                AlignmentOperation::Deletion,
                AlignmentOperation::Deletion,
                AlignmentOperation::Match,
            ],
        ),
        (
            b"GGAC".as_slice(),
            b"AC".as_slice(),
            0,
            vec![
                AlignmentOperation::Deletion,
                AlignmentOperation::Deletion,
                AlignmentOperation::Match,
                AlignmentOperation::Match,
            ],
        ),
        (
            b"AC".as_slice(),
            b"ACGG".as_slice(),
            0,
            vec![
                AlignmentOperation::Match,
                AlignmentOperation::Match,
                AlignmentOperation::Insertion,
                AlignmentOperation::Insertion,
            ],
        ),
    ];

    for (first_raw, second_raw, expected_score, expected_operations) in cases {
        let first = dna(first_raw);
        let second = dna(second_raw);
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        assert_eq!(alignment.score(), expected_score);
        assert_eq!(alignment.operations(), expected_operations);
        verify_affine_end_gap_alignment(
            &first,
            &second,
            &alignment,
            &strong_mismatch,
            gap(),
            FreeEndGaps::global(),
        )
        .expect("manual affine path replays");
    }
}

#[test]
fn fitting_overlap_and_semiglobal_presets_select_expected_ranges() {
    let cases = [
        (
            b"ACGT".as_slice(),
            b"TTACGTGG".as_slice(),
            FreeEndGaps::fit_first_in_second(),
            0..4,
            2..6,
        ),
        (
            b"GGACGTAA".as_slice(),
            b"ACGT".as_slice(),
            FreeEndGaps::fit_second_in_first(),
            2..6,
            0..4,
        ),
        (
            b"TTACGT".as_slice(),
            b"ACGTGG".as_slice(),
            FreeEndGaps::first_suffix_second_prefix(),
            2..6,
            0..4,
        ),
        (
            b"ACGTGG".as_slice(),
            b"TTACGT".as_slice(),
            FreeEndGaps::second_suffix_first_prefix(),
            0..4,
            2..6,
        ),
        (
            b"TTACGT".as_slice(),
            b"ACGTGG".as_slice(),
            FreeEndGaps::semiglobal(),
            2..6,
            0..4,
        ),
    ];

    for (first_raw, second_raw, free_ends, first_range, second_range) in cases {
        let first = dna(first_raw);
        let second = dna(second_raw);
        let aligner = AffineEndGapAligner::new(scoring(), gap(), free_ends);
        let alignment = aligner.align(&first, &second).expect("alignment succeeds");
        assert_eq!(alignment.score(), 12);
        assert_eq!(alignment.first_range(), first_range);
        assert_eq!(alignment.second_range(), second_range);
        assert_eq!(alignment.operations(), &[AlignmentOperation::Match; 4]);
        verify_affine_end_gap_alignment(
            &first,
            &second,
            &alignment,
            aligner.substitution(),
            aligner.gap(),
            aligner.free_ends(),
        )
        .expect("free-end result replays");
    }
}

#[test]
fn all_empty_input_and_four_end_combinations_are_valid() {
    let empty = dna(b"");
    let sequence = dna(b"ACG");

    for bits in 0_u8..16 {
        let free_ends =
            FreeEndGaps::new(bits & 1 != 0, bits & 2 != 0, bits & 4 != 0, bits & 8 != 0);
        let aligner = AffineEndGapAligner::new(scoring(), gap(), free_ends);
        for (first, second) in [(&empty, &empty), (&empty, &sequence), (&sequence, &empty)] {
            let alignment = aligner
                .align(first, second)
                .expect("empty-edge alignment succeeds");
            verify_affine_end_gap_alignment(
                first,
                second,
                &alignment,
                aligner.substitution(),
                aligner.gap(),
                free_ends,
            )
            .expect("empty-edge result verifies");
        }
    }
}

#[test]
fn affine_local_excludes_nonhomologous_ends_and_keeps_an_internal_gap() {
    let scorer = MatchMismatch::new(6, -10);
    let affine_gap = AffineGap::new(-3, -1).expect("valid affine gap");
    let first = dna(b"TTACGGGTAA");
    let second = dna(b"CCACTCC");
    let aligner = AffineLocalAligner::new(scorer, affine_gap);
    let alignment = aligner
        .align(&first, &second)
        .expect("local alignment succeeds");

    assert_eq!(alignment.score(), 13);
    assert_eq!(alignment.first_range(), 2..8);
    assert_eq!(alignment.second_range(), 2..5);
    assert_eq!(
        alignment.operations(),
        &[
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
            AlignmentOperation::Deletion,
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
        ]
    );
    verify_local_affine_alignment(
        &first,
        &second,
        &alignment,
        aligner.substitution(),
        aligner.gap(),
    )
    .expect("local affine result replays");
}

#[test]
fn local_empty_different_identical_and_repetitive_cases_are_deterministic() {
    let aligner = AffineLocalAligner::new(scoring(), gap());
    let empty = dna(b"");
    let different_first = dna(b"AAAA");
    let different_second = dna(b"TTTT");

    for (first, second) in [
        (&empty, &empty),
        (&empty, &different_first),
        (&different_first, &empty),
        (&different_first, &different_second),
    ] {
        let alignment = aligner
            .align(first, second)
            .expect("local alignment succeeds");
        assert_eq!(alignment.score(), 0);
        assert_eq!(alignment.first_range(), 0..0);
        assert_eq!(alignment.second_range(), 0..0);
        assert!(alignment.operations().is_empty());
    }

    let identical = dna(b"ACGTN");
    let alignment = aligner
        .align(&identical, &identical)
        .expect("identical alignment succeeds");
    assert_eq!(alignment.score(), 15);
    assert_eq!(alignment.operations(), &[AlignmentOperation::Match; 5]);

    let repetitive_first = dna(b"AAAAAA");
    let repetitive_second = dna(b"AAAA");
    let first_run = aligner
        .align(&repetitive_first, &repetitive_second)
        .expect("repetitive alignment succeeds");
    for _ in 0..10 {
        assert_eq!(
            aligner
                .align(&repetitive_first, &repetitive_second)
                .expect("repeated alignment succeeds"),
            first_run
        );
    }
}

#[test]
fn adjacent_opposite_gaps_are_valid_when_better_than_a_substitution() {
    let first = dna(b"A");
    let second = dna(b"C");
    let scorer = MatchMismatch::new(2, -20);
    let affine_gap = AffineGap::new(-3, -1).expect("valid gap");
    let aligner = AffineEndGapAligner::new(scorer, affine_gap, FreeEndGaps::global());
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), -6);
    assert_eq!(alignment.columns(), 2);
    assert!(
        alignment
            .operations()
            .contains(&AlignmentOperation::Deletion)
    );
    assert!(
        alignment
            .operations()
            .contains(&AlignmentOperation::Insertion)
    );
    verify_affine_end_gap_alignment(
        &first,
        &second,
        &alignment,
        &scorer,
        affine_gap,
        FreeEndGaps::global(),
    )
    .expect("opposite-gap result replays");
}

#[test]
fn dna_rna_and_protein_use_the_same_public_api() {
    let ends = FreeEndGaps::fit_first_in_second();
    let aligner = AffineEndGapAligner::new(scoring(), gap(), ends);

    let dna_first = DnaSequence::from_ascii(b"ACGTRYN").expect("valid DNA");
    let dna_second = DnaSequence::from_ascii(b"TTACGTRYNGG").expect("valid DNA");
    assert_eq!(
        aligner
            .align(&dna_first, &dna_second)
            .expect("DNA alignment succeeds")
            .score(),
        21
    );

    let rna_first = RnaSequence::from_ascii(b"ACGURYN").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"UUACGURYNGG").expect("valid RNA");
    assert_eq!(
        aligner
            .align(&rna_first, &rna_second)
            .expect("RNA alignment succeeds")
            .score(),
        21
    );

    let protein_first = ProteinSequence::from_ascii(b"MTEYK").expect("valid protein");
    let protein_second = ProteinSequence::from_ascii(b"GGMTEYKPP").expect("valid protein");
    assert_eq!(
        aligner
            .align(&protein_first, &protein_second)
            .expect("protein alignment succeeds")
            .score(),
        15
    );

    let local = AffineLocalAligner::new(scoring(), gap());
    assert_eq!(
        local
            .align(&protein_first, &protein_second)
            .expect("local protein alignment succeeds")
            .score(),
        15
    );
}

#[test]
fn custom_substitution_matrix_is_honored() {
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        if first == second {
            4
        } else if first == b'N' || second == b'N' {
            1
        } else {
            -8
        }
    });
    let first = dna(b"AN");
    let second = dna(b"TTAACT");
    let ends = FreeEndGaps::fit_first_in_second();
    let aligner = AffineEndGapAligner::new(matrix.clone(), gap(), ends);
    let alignment = aligner.align(&first, &second).expect("alignment succeeds");

    assert_eq!(alignment.score(), 5);
    assert_eq!(alignment.second_range(), 3..5);
    verify_affine_end_gap_alignment(&first, &second, &alignment, &matrix, gap(), ends)
        .expect("matrix-scored result verifies");
}

#[test]
fn dense_score_row_fast_path_matches_scalar_matrix_lookups_in_both_kernels() {
    let first = dna(b"ACGTRYSWKMBDHVNACGT");
    let second = dna(b"NVDHBMKWSRYTGCAN");
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        i64::from(first) * 2 - i64::from(second) * 3
    });
    let gap = AffineGap::new(-9, -2).expect("valid gap");
    let free_ends = FreeEndGaps::new(true, false, false, true);

    let cached_end_gap = AffineEndGapAligner::new(RowCachedMatrix(matrix.clone()), gap, free_ends)
        .align(&first, &second)
        .expect("cached end-gap alignment succeeds");
    let scalar_end_gap = AffineEndGapAligner::new(matrix.clone(), gap, free_ends)
        .align(&first, &second)
        .expect("scalar end-gap alignment succeeds");
    assert_eq!(cached_end_gap, scalar_end_gap);

    let cached_local = AffineLocalAligner::new(RowCachedMatrix(matrix.clone()), gap)
        .align(&first, &second)
        .expect("cached local alignment succeeds");
    let scalar_local = AffineLocalAligner::new(matrix, gap)
        .align(&first, &second)
        .expect("scalar local alignment succeeds");
    assert_eq!(cached_local, scalar_local);
}

#[derive(Clone, Copy)]
struct UnboundedMatchMismatch;

impl SubstitutionScorer<Dna> for UnboundedMatchMismatch {
    fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
        if first == second { 3 } else { -4 }
    }
}

#[test]
fn unbounded_scorer_uses_checked_arithmetic() {
    let first = dna(b"ACGT");
    let second = dna(b"TTACGTGG");
    let ends = FreeEndGaps::fit_first_in_second();
    let alignment = AffineEndGapAligner::new(UnboundedMatchMismatch, gap(), ends)
        .align(&first, &second)
        .expect("checked affine alignment succeeds");
    assert_eq!(alignment.score(), 12);

    let local = AffineLocalAligner::new(UnboundedMatchMismatch, gap())
        .align(&first, &second)
        .expect("checked local alignment succeeds");
    assert_eq!(local.score(), 12);
}

#[test]
fn positive_score_overflow_is_reported() {
    let sequence = dna(b"AA");
    let affine_gap = AffineGap::new(0, 0).expect("zero affine gap is valid");

    for result in [
        AffineEndGapAligner::new(
            MatchMismatch::new(Score::MAX, 0),
            affine_gap,
            FreeEndGaps::semiglobal(),
        )
        .align(&sequence, &sequence),
        AffineLocalAligner::new(MatchMismatch::new(Score::MAX, 0), affine_gap)
            .align(&sequence, &sequence),
    ] {
        assert!(matches!(result, Err(AlignmentError::ScoreOverflow { .. })));
    }
}

#[test]
fn exact_minimum_single_gap_boundary_does_not_evaluate_a_later_extension() {
    let result = AffineEndGapAligner::new(
        MatchMismatch::new(1, -1),
        AffineGap::new(Score::MIN, -1).expect("valid extreme gap"),
        FreeEndGaps::global(),
    )
    .align(&dna(b""), &dna(b"A"))
    .expect("the represented length-one boundary is valid");
    assert_eq!(result.score(), Score::MIN);
}

#[test]
fn negatively_overflowing_losing_candidates_are_discarded() {
    let first = dna(b"AA");
    let second = dna(b"TT");
    let scorer = MatchMismatch::new(0, Score::MIN);
    let affine_gap = AffineGap::new(-1, -1).expect("valid affine gap");

    let global = AffineEndGapAligner::new(scorer, affine_gap, FreeEndGaps::global())
        .align(&first, &second)
        .expect("representable gap path succeeds");
    assert_eq!(global.score(), -4);

    let local = AffineLocalAligner::new(scorer, affine_gap)
        .align(&first, &second)
        .expect("zero local path succeeds");
    assert_eq!(local.score(), 0);
}

#[test]
fn independent_verifiers_reject_tampered_scores_and_operations() {
    let first = dna(b"ACGT");
    let second = dna(b"TTACGTGG");
    let ends = FreeEndGaps::fit_first_in_second();
    let alignment = AffineEndGapAligner::new(scoring(), gap(), ends)
        .align(&first, &second)
        .expect("alignment succeeds");

    let bad_score = PairwiseAlignment::try_new(
        alignment.score() + 1,
        alignment.first_range(),
        alignment.second_range(),
        alignment.operations().to_vec(),
    )
    .expect("structurally valid tampered result");
    assert_eq!(
        verify_affine_end_gap_alignment(&first, &second, &bad_score, &scoring(), gap(), ends),
        Err(VerificationError::ScoreMismatch {
            reported: 13,
            recomputed: 12,
        })
    );

    let mut bad_operations = alignment.operations().to_vec();
    bad_operations[0] = AlignmentOperation::Substitution;
    let bad_operation = PairwiseAlignment::try_new(
        alignment.score(),
        alignment.first_range(),
        alignment.second_range(),
        bad_operations,
    )
    .expect("structurally valid tampered operation");
    assert!(matches!(
        verify_affine_end_gap_alignment(&first, &second, &bad_operation, &scoring(), gap(), ends,),
        Err(VerificationError::OperationKindMismatch { column: 0, .. })
    ));
}

#[test]
fn local_verifier_rejects_noncanonical_zero_and_nonpositive_prefix() {
    let first = dna(b"TA");
    let second = dna(b"CA");
    let noncanonical =
        PairwiseAlignment::try_new(0, 0..1, 0..0, vec![AlignmentOperation::Deletion])
            .expect("structurally valid zero result");
    assert_eq!(
        verify_local_affine_alignment(&first, &second, &noncanonical, &scoring(), gap()),
        Err(VerificationError::NonCanonicalEmptyLocalAlignment)
    );

    let crosses_restart = PairwiseAlignment::try_new(
        1,
        0..2,
        0..2,
        vec![AlignmentOperation::Substitution, AlignmentOperation::Match],
    )
    .expect("structurally valid positive result");
    let prefix_scorer = MatchMismatch::new(2, -1);
    assert_eq!(
        verify_local_affine_alignment(&first, &second, &crosses_restart, &prefix_scorer, gap()),
        Err(VerificationError::NonPositiveLocalPrefix {
            column: 0,
            score: -1,
        })
    );
}
