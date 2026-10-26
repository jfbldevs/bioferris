//! Equivalence and reuse tests for caller-owned affine scratch storage.

use align_core::{
    AffineGap, Dna, DnaSequence, MatchMismatch, ProteinSequence, RnaSequence, Score, Sequence,
    SubstitutionScorer, Symbol,
};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, AffineWorkspace, AlignmentError, FreeEndGaps, Gotoh,
    verify_affine_end_gap_alignment, verify_global_affine_alignment, verify_local_affine_alignment,
};
use proptest::prelude::*;

fn scoring() -> MatchMismatch {
    MatchMismatch::new(3, -4)
}

fn gap() -> AffineGap {
    AffineGap::new(-5, -1).expect("valid affine gap")
}

fn dna_from_indices(indices: &[u8]) -> Sequence<Dna> {
    let symbols = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index)])
        .collect::<Vec<_>>();
    Sequence::from_ascii(symbols).expect("generated DNA is valid")
}

#[test]
fn one_workspace_is_reusable_across_shapes_modes_and_alphabets() {
    let mut workspace = AffineWorkspace::new();
    workspace
        .reserve(20, 30)
        .expect("explicit workspace reservation succeeds");
    assert!(workspace.score_row_capacity() >= 31);
    assert!(workspace.state_validity_capacity() >= 31);
    assert!(workspace.traceback_capacity() >= 21 * 31);

    let dna_first = DnaSequence::from_ascii(b"ACGTRYSWKMBDHVN").expect("valid DNA");
    let dna_second = DnaSequence::from_ascii(b"NVDHBMKWSRYTGCA").expect("valid DNA");
    let global = Gotoh::new(scoring(), gap());
    let retained_global = global
        .align_with_workspace(&dna_first, &dna_second, &mut workspace)
        .expect("workspace-backed global alignment succeeds");
    assert_eq!(
        retained_global,
        global
            .align(&dna_first, &dna_second)
            .expect("fresh global alignment succeeds")
    );

    let protein_first =
        ProteinSequence::from_ascii(b"QQMTEYKLVVVGAGGVGKSALTIQLIQNHF").expect("valid protein");
    let protein_second =
        ProteinSequence::from_ascii(b"GGMTEYKLVVVGAGGVGKSALTIQLIQPP").expect("valid protein");
    let local = AffineLocalAligner::new(scoring(), gap());
    let local_with_workspace = local
        .align_with_workspace(&protein_first, &protein_second, &mut workspace)
        .expect("workspace-backed local alignment succeeds");
    assert_eq!(
        local_with_workspace,
        local
            .align(&protein_first, &protein_second)
            .expect("fresh local alignment succeeds")
    );

    let rna_first = RnaSequence::from_ascii(b"AUGCUARYSWKMBDHVN").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"CCAUGCUARYSWKMBDHVNGG").expect("valid RNA");
    let free_ends = FreeEndGaps::fit_first_in_second();
    let end_gap = AffineEndGapAligner::new(scoring(), gap(), free_ends);
    let end_gap_with_workspace = end_gap
        .align_with_workspace(&rna_first, &rna_second, &mut workspace)
        .expect("workspace-backed end-gap alignment succeeds");
    assert_eq!(
        end_gap_with_workspace,
        end_gap
            .align(&rna_first, &rna_second)
            .expect("fresh end-gap alignment succeeds")
    );

    // A later workspace use cannot mutate an already returned owned alignment.
    verify_global_affine_alignment(
        &dna_first,
        &dna_second,
        &retained_global,
        global.substitution(),
        global.gap(),
    )
    .expect("retained result remains independently valid");
    verify_local_affine_alignment(
        &protein_first,
        &protein_second,
        &local_with_workspace,
        local.substitution(),
        local.gap(),
    )
    .expect("local workspace result verifies");
    verify_affine_end_gap_alignment(
        &rna_first,
        &rna_second,
        &end_gap_with_workspace,
        end_gap.substitution(),
        end_gap.gap(),
        free_ends,
    )
    .expect("end-gap workspace result verifies");
}

#[derive(Clone, Copy)]
struct UnboundedMaximum;

impl SubstitutionScorer<Dna> for UnboundedMaximum {
    fn score(&self, _first: Symbol<Dna>, _second: Symbol<Dna>) -> Score {
        Score::MAX
    }
}

#[test]
fn workspace_remains_valid_after_checked_arithmetic_failure() {
    let sequence = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let zero_gap = AffineGap::new(0, 0).expect("valid zero gap");
    let mut workspace = AffineWorkspace::new();
    let overflow = Gotoh::new(UnboundedMaximum, zero_gap).align_with_workspace(
        &sequence,
        &sequence,
        &mut workspace,
    );
    assert!(matches!(
        overflow,
        Err(AlignmentError::ScoreOverflow { .. })
    ));
    assert!(workspace.state_validity_capacity() > sequence.len());

    let recovered = Gotoh::new(scoring(), gap())
        .align_with_workspace(&sequence, &sequence, &mut workspace)
        .expect("bounded alignment after an error succeeds");
    assert_eq!(recovered.score(), 6);
}

#[test]
fn local_reuse_resets_every_top_boundary_trace() {
    let first = DnaSequence::from_ascii(b"C").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"GC").expect("valid DNA");
    let global = Gotoh::new(scoring(), gap());
    let local = AffineLocalAligner::new(scoring(), gap());
    let mut workspace = AffineWorkspace::new();

    global
        .align_with_workspace(&second, &first, &mut workspace)
        .expect("global call primes non-start boundary traces");
    let observed = local
        .align_with_workspace(&first, &second, &mut workspace)
        .expect("local reuse succeeds");
    let expected = local.align(&first, &second).expect("fresh local succeeds");
    assert_eq!(observed, expected);
    assert_eq!(observed.first_range(), 0..1);
    assert_eq!(observed.second_range(), 1..2);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn reused_and_fresh_affine_calls_are_exactly_identical(
        first_indices in prop::collection::vec(0_u8..4, 0..10),
        second_indices in prop::collection::vec(0_u8..4, 0..10),
        free_end_bits in 0_u8..16,
    ) {
        let first = dna_from_indices(&first_indices);
        let second = dna_from_indices(&second_indices);
        let free_ends = FreeEndGaps::new(
            free_end_bits & 1 != 0,
            free_end_bits & 2 != 0,
            free_end_bits & 4 != 0,
            free_end_bits & 8 != 0,
        );
        let global = Gotoh::new(scoring(), gap());
        let end_gap = AffineEndGapAligner::new(scoring(), gap(), free_ends);
        let local = AffineLocalAligner::new(scoring(), gap());
        let mut workspace = AffineWorkspace::new();

        let expected_global = global.align(&first, &second).expect("fresh global succeeds");
        let observed_global = global
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("reused global succeeds");
        prop_assert_eq!(&observed_global, &expected_global);

        // Change orientation and boundary policy before returning to the original
        // shape, exercising both shrink/grow and stale-boundary reset behavior.
        let expected_end_gap = end_gap.align(&second, &first).expect("fresh end-gap succeeds");
        let observed_end_gap = end_gap
            .align_with_workspace(&second, &first, &mut workspace)
            .expect("reused end-gap succeeds");
        prop_assert_eq!(observed_end_gap, expected_end_gap);

        let expected_local = local.align(&first, &second).expect("fresh local succeeds");
        let observed_local = local
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("reused local succeeds");
        prop_assert_eq!(observed_local, expected_local);

        let repeated_global = global
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("repeated reused global succeeds");
        prop_assert_eq!(&repeated_global, &expected_global);
    }
}
