//! Differential acceptance tests for the runtime-dispatched affine SIMD backend.

use core::mem::size_of;
use std::sync::OnceLock;

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, DenseSubstitutionMatrix, Dna, MatchMismatch, Protein,
    ProteinSequence, RnaSequence, Score, Sequence, SubstitutionScorer, Symbol,
};
use align_pairwise::{
    AffineBackend, AffineEndGapAligner, AffineKernelPolicy, AffineLocalAligner, AffineWorkspace,
    AlignmentError, FreeEndGaps, Gotoh, verify_affine_end_gap_alignment,
    verify_global_affine_alignment, verify_local_affine_alignment,
};
use proptest::prelude::*;

fn dna_from_indices(indices: &[u8]) -> Sequence<Dna> {
    let symbols = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index)])
        .collect::<Vec<_>>();
    Sequence::from_ascii(symbols).expect("generated DNA is valid")
}

fn workspace(policy: AffineKernelPolicy) -> AffineWorkspace {
    let mut workspace = AffineWorkspace::new();
    workspace.set_kernel_policy(policy);
    workspace
}

fn is_simd_backend(backend: Option<AffineBackend>) -> bool {
    matches!(
        backend,
        Some(AffineBackend::SimdI16 { lanes } | AffineBackend::SimdI32 { lanes }) if lanes > 1
    )
}

fn simd_is_available() -> bool {
    static AVAILABLE: OnceLock<bool> = OnceLock::new();
    *AVAILABLE.get_or_init(|| {
        let sequence = dna_from_indices(&[0]);
        let scoring = MatchMismatch::new(2, -3);
        let gap = AffineGap::new(-5, -1).expect("valid affine gap");
        let mut workspace = workspace(AffineKernelPolicy::Simd);
        match Gotoh::new(scoring, gap).align_with_workspace(&sequence, &sequence, &mut workspace) {
            Ok(_) => is_simd_backend(workspace.last_backend()),
            Err(AlignmentError::SimdUnavailable { .. }) => false,
            Err(error) => panic!("SIMD capability probe failed unexpectedly: {error}"),
        }
    })
}

#[derive(Clone, Debug)]
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

#[derive(Clone, Copy, Debug)]
struct UnknownBounds;

impl SubstitutionScorer<Dna> for UnknownBounds {
    fn score(&self, first: Symbol<Dna>, second: Symbol<Dna>) -> Score {
        if first == second { 3 } else { -4 }
    }
}

fn assert_all_modes_match<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    scoring: &S,
    gap: AffineGap,
    free_ends: FreeEndGaps,
) where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut scalar = workspace(AffineKernelPolicy::Scalar);
    let mut simd = workspace(AffineKernelPolicy::Simd);

    let global = Gotoh::new(scoring, gap);
    assert_eq!(
        global
            .align_with_workspace(first, second, &mut simd)
            .expect("SIMD global alignment succeeds"),
        global
            .align_with_workspace(first, second, &mut scalar)
            .expect("scalar global alignment succeeds")
    );

    let end_gap = AffineEndGapAligner::new(scoring, gap, free_ends);
    assert_eq!(
        end_gap
            .align_with_workspace(first, second, &mut simd)
            .expect("SIMD end-gap alignment succeeds"),
        end_gap
            .align_with_workspace(first, second, &mut scalar)
            .expect("scalar end-gap alignment succeeds")
    );

    let local = AffineLocalAligner::new(scoring, gap);
    assert_eq!(
        local
            .align_with_workspace(first, second, &mut simd)
            .expect("SIMD local alignment succeeds"),
        local
            .align_with_workspace(first, second, &mut scalar)
            .expect("scalar local alignment succeeds")
    );
    assert!(is_simd_backend(simd.last_backend()));
}

#[test]
fn strict_policy_reports_empty_and_wide_score_rejections() {
    let empty = dna_from_indices(&[]);
    let one = dna_from_indices(&[0]);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let mut strict = workspace(AffineKernelPolicy::Simd);
    let empty_result =
        Gotoh::new(MatchMismatch::new(2, -3), gap).align_with_workspace(&empty, &one, &mut strict);
    assert!(matches!(
        empty_result,
        Err(AlignmentError::SimdUnavailable { .. })
    ));
    assert_eq!(strict.last_backend(), None);

    let two = dna_from_indices(&[0, 0]);
    let wide = MatchMismatch::new(i64::from(i32::MAX), -1);
    let wide_result = Gotoh::new(wide, gap).align_with_workspace(&two, &two, &mut strict);
    assert!(matches!(
        wide_result,
        Err(AlignmentError::SimdUnavailable { .. })
    ));
    assert_eq!(strict.last_backend(), None);

    let mut automatic = workspace(AffineKernelPolicy::Auto);
    let fallback = Gotoh::new(wide, gap)
        .align_with_workspace(&two, &two, &mut automatic)
        .expect("automatic policy falls back to exact signed-64-bit scalar arithmetic");
    assert_eq!(fallback.score(), 2 * i64::from(i32::MAX));
    assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));

    let prefix_unsafe_gap =
        AffineGap::new(-1_000_000_000, -1_000_000_000).expect("valid large affine gap");
    let prefix_unsafe = Gotoh::new(MatchMismatch::new(1, -1), prefix_unsafe_gap);
    assert!(matches!(
        prefix_unsafe.align_with_workspace(&one, &one, &mut strict),
        Err(AlignmentError::SimdUnavailable { .. })
    ));
    let fallback = prefix_unsafe
        .align_with_workspace(&one, &one, &mut automatic)
        .expect("automatic policy avoids unsafe signed-32-bit prefix arithmetic");
    assert_eq!(fallback.score(), 1);
    assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));
}

#[test]
fn strict_simd_selects_the_narrowest_proven_exact_score_width() {
    if !simd_is_available() {
        return;
    }
    let sequence = dna_from_indices(&[0, 0]);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let mut simd = workspace(AffineKernelPolicy::Simd);
    let mut scalar = workspace(AffineKernelPolicy::Scalar);

    let narrow = Gotoh::new(MatchMismatch::new(8_000, -2), gap);
    let expected = narrow
        .align_with_workspace(&sequence, &sequence, &mut scalar)
        .expect("scalar narrow alignment succeeds");
    let observed = narrow
        .align_with_workspace(&sequence, &sequence, &mut simd)
        .expect("signed-16-bit SIMD alignment succeeds");
    assert_eq!(observed, expected);
    assert_eq!(
        simd.last_backend(),
        Some(AffineBackend::SimdI16 { lanes: 16 })
    );

    let wide = Gotoh::new(MatchMismatch::new(8_192, -2), gap);
    let expected = wide
        .align_with_workspace(&sequence, &sequence, &mut scalar)
        .expect("scalar wide alignment succeeds");
    let observed = wide
        .align_with_workspace(&sequence, &sequence, &mut simd)
        .expect("signed-32-bit SIMD fallback succeeds");
    assert_eq!(observed, expected);
    assert_eq!(
        simd.last_backend(),
        Some(AffineBackend::SimdI32 { lanes: 8 })
    );

    let large_extension = AffineGap::new(-5_000, -5_000).expect("valid large gap");
    Gotoh::new(MatchMismatch::new(1, -1), large_extension)
        .align_with_workspace(&sequence, &sequence, &mut simd)
        .expect("wide prefix arithmetic falls back to signed-32-bit SIMD");
    assert_eq!(
        simd.last_backend(),
        Some(AffineBackend::SimdI32 { lanes: 8 })
    );
}

#[test]
fn unknown_score_bounds_are_rejected_strictly_and_fall_back_automatically() {
    let first = dna_from_indices(&[0, 1, 2, 3]);
    let second = dna_from_indices(&[0, 2, 3]);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let aligner = Gotoh::new(UnknownBounds, gap);

    let mut strict = workspace(AffineKernelPolicy::Simd);
    assert!(matches!(
        aligner.align_with_workspace(&first, &second, &mut strict),
        Err(AlignmentError::SimdUnavailable { .. })
    ));
    assert_eq!(strict.last_backend(), None);

    let mut automatic = workspace(AffineKernelPolicy::Auto);
    let observed = aligner
        .align_with_workspace(&first, &second, &mut automatic)
        .expect("automatic checked scalar fallback succeeds");
    let mut scalar = workspace(AffineKernelPolicy::Scalar);
    let expected = aligner
        .align_with_workspace(&first, &second, &mut scalar)
        .expect("forced checked scalar alignment succeeds");
    assert_eq!(observed, expected);
    assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));
}

#[test]
fn simd_memory_model_and_reservation_match_active_buffers() {
    let first_len = 3;
    let second_len = 4;
    let estimate = AffineWorkspace::simd_memory_requirement::<Dna>(first_len, second_len)
        .expect("small SIMD memory estimate succeeds");
    let cells = (first_len + 1) * (second_len + 1);
    let final_score_row_bytes = (second_len + 1) * size_of::<Score>();
    let padded_first_len = 8;
    let score_column_bytes = 2 * (padded_first_len + 1) * size_of::<i32>();
    let substitution_profile_bytes = padded_first_len * Dna::SYMBOLS.len() * size_of::<i32>();
    let second_index_bytes = second_len * size_of::<usize>();
    let traceback_matrix_bytes = cells * size_of::<u8>();
    let maximum_operation_bytes = (first_len + second_len) * size_of::<AlignmentOperation>();

    assert_eq!(estimate.cells(), cells);
    assert_eq!(estimate.final_score_row_bytes(), final_score_row_bytes);
    assert_eq!(estimate.score_column_bytes(), score_column_bytes);
    assert_eq!(
        estimate.substitution_profile_bytes(),
        substitution_profile_bytes
    );
    assert_eq!(estimate.second_index_bytes(), second_index_bytes);
    assert_eq!(estimate.traceback_matrix_bytes(), traceback_matrix_bytes);
    assert_eq!(estimate.maximum_operation_bytes(), maximum_operation_bytes);
    assert_eq!(
        estimate.peak_buffer_bytes(),
        final_score_row_bytes
            + score_column_bytes
            + substitution_profile_bytes
            + second_index_bytes
            + traceback_matrix_bytes
            + maximum_operation_bytes
    );

    let mut workspace = AffineWorkspace::new();
    workspace
        .reserve_simd::<Dna>(first_len, second_len)
        .expect("SIMD reservation succeeds");
    assert!(workspace.simd_final_score_row_capacity() > second_len);
    assert!(workspace.simd_score_row_capacity() > padded_first_len);
    assert!(workspace.simd_profile_capacity() >= padded_first_len * Dna::SYMBOLS.len());
    assert!(workspace.simd_second_index_capacity() >= second_len);
    assert!(workspace.simd_traceback_capacity() >= cells);
    assert!(AffineWorkspace::simd_memory_requirement::<Dna>(usize::MAX, 0).is_err());
}

#[test]
fn automatic_policy_obeys_the_calibrated_cell_gate() {
    let scoring = MatchMismatch::new(2, -3);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let short = dna_from_indices(&[0; 63]);
    let threshold = dna_from_indices(&[0; 64]);
    let mut workspace = workspace(AffineKernelPolicy::Auto);

    Gotoh::new(scoring, gap)
        .align_with_workspace(&short, &short, &mut workspace)
        .expect("below-threshold alignment succeeds");
    assert_eq!(workspace.last_backend(), Some(AffineBackend::Scalar));

    Gotoh::new(scoring, gap)
        .align_with_workspace(&threshold, &threshold, &mut workspace)
        .expect("threshold alignment succeeds");
    if simd_is_available() {
        assert!(is_simd_backend(workspace.last_backend()));
    } else {
        assert_eq!(workspace.last_backend(), Some(AffineBackend::Scalar));
    }
}

#[test]
fn automatic_policy_obeys_the_calibrated_i16_shape_gate() {
    let scoring = MatchMismatch::new(2, -3);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let aligner = Gotoh::new(scoring, gap);
    let mut workspace = workspace(AffineKernelPolicy::Auto);

    for (first_len, second_len, simd_eligible) in [
        (5, 820, false),
        (6, 683, false),
        (7, 586, true),
        (8, 512, true),
        (9, 456, true),
        (11, 373, true),
        (12, 342, true),
    ] {
        let first = dna_from_indices(&vec![0; first_len]);
        let second = dna_from_indices(&vec![1; second_len]);
        aligner
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("automatic alignment succeeds");
        let expected = if simd_eligible && simd_is_available() {
            is_simd_backend(workspace.last_backend())
        } else {
            workspace.last_backend() == Some(AffineBackend::Scalar)
        };
        assert!(
            expected,
            "unexpected backend for {first_len}x{second_len} automatic alignment"
        );
    }
}

#[test]
fn automatic_policy_obeys_the_calibrated_i32_shape_gate() {
    let scoring = MatchMismatch::new(10_000, -4);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let aligner = Gotoh::new(scoring, gap);
    let mut workspace = workspace(AffineKernelPolicy::Auto);

    for (first_len, second_len, simd_eligible) in [
        (5, 820, false),
        (6, 683, true),
        (8, 512, true),
        (9, 456, false),
        (10, 410, false),
        (11, 373, true),
        (12, 342, true),
    ] {
        let first = dna_from_indices(&vec![0; first_len]);
        let second = dna_from_indices(&vec![1; second_len]);
        aligner
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("automatic signed-32-bit alignment succeeds");
        let expected = if simd_eligible && simd_is_available() {
            matches!(
                workspace.last_backend(),
                Some(AffineBackend::SimdI32 { lanes: 8 })
            )
        } else {
            workspace.last_backend() == Some(AffineBackend::Scalar)
        };
        assert!(
            expected,
            "unexpected signed-32-bit backend for {first_len}x{second_len} automatic alignment"
        );
    }
}

#[test]
fn simd_profile_supports_cached_matrices_rna_and_proteins() {
    if !simd_is_available() {
        return;
    }
    let gap = AffineGap::new(-7, -1).expect("valid affine gap");
    let free_ends = FreeEndGaps::new(true, false, false, true);

    let dna_first = dna_from_indices(&[0, 1, 2, 3, 3, 2, 1, 0, 0, 2, 1]);
    let dna_second = dna_from_indices(&[3, 1, 2, 0, 3, 2, 1, 0, 1]);
    let matrix = DenseSubstitutionMatrix::<Dna>::from_score_fn(|first, second| {
        if first == second {
            5
        } else {
            -i64::from(1 + (first ^ second) % 7)
        }
    });
    assert_all_modes_match(
        &dna_first,
        &dna_second,
        &RowCachedMatrix(matrix),
        gap,
        free_ends,
    );

    let rna_first = RnaSequence::from_ascii(b"AUGCURYSWKMBDHVNAUGCU").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"GCAUNVDHBMKWSRYAUG").expect("valid RNA");
    assert_all_modes_match(
        &rna_first,
        &rna_second,
        &MatchMismatch::new(4, -3),
        gap,
        free_ends,
    );

    let protein_first = ProteinSequence::from_ascii(b"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY")
        .expect("valid protein");
    let protein_second =
        ProteinSequence::from_ascii(b"GGMTEYKLVVVGAGGVGKSALTIQLIQPPFVDKYD").expect("valid protein");
    assert_all_modes_match(
        &protein_first,
        &protein_second,
        &MatchMismatch::new(6, -5),
        gap,
        free_ends,
    );
}

fn binary_dna_sequences(maximum_len: usize) -> Vec<Sequence<Dna>> {
    let mut sequences = Vec::new();
    for len in 1..=maximum_len {
        for encoded in 0..(1_usize << len) {
            let indices = (0..len)
                .map(|position| u8::from((encoded >> position) & 1 != 0))
                .collect::<Vec<_>>();
            sequences.push(dna_from_indices(&indices));
        }
    }
    sequences
}

#[test]
fn exhaustive_short_sequences_match_scalar_in_every_affine_mode() {
    if !simd_is_available() {
        return;
    }
    let sequences = binary_dna_sequences(3);
    let configurations = [
        (
            MatchMismatch::new(1, -1),
            AffineGap::new(-1, -1).expect("valid tie-heavy gap"),
        ),
        (
            MatchMismatch::new(3, -4),
            AffineGap::new(-5, -1).expect("valid affine gap"),
        ),
    ];
    let mut scalar = workspace(AffineKernelPolicy::Scalar);
    let mut simd = workspace(AffineKernelPolicy::Simd);

    for (scoring, gap) in configurations {
        let global = Gotoh::new(scoring, gap);
        let local = AffineLocalAligner::new(scoring, gap);
        for first in &sequences {
            for second in &sequences {
                let expected_global = global
                    .align_with_workspace(first, second, &mut scalar)
                    .expect("short scalar global alignment succeeds");
                let observed_global = global
                    .align_with_workspace(first, second, &mut simd)
                    .expect("short SIMD global alignment succeeds");
                assert_eq!(observed_global, expected_global);

                let expected_local = local
                    .align_with_workspace(first, second, &mut scalar)
                    .expect("short scalar local alignment succeeds");
                let observed_local = local
                    .align_with_workspace(first, second, &mut simd)
                    .expect("short SIMD local alignment succeeds");
                assert_eq!(observed_local, expected_local);

                for free_end_bits in 0_u8..16 {
                    let free_ends = FreeEndGaps::new(
                        free_end_bits & 1 != 0,
                        free_end_bits & 2 != 0,
                        free_end_bits & 4 != 0,
                        free_end_bits & 8 != 0,
                    );
                    let aligner = AffineEndGapAligner::new(scoring, gap, free_ends);
                    let expected = aligner
                        .align_with_workspace(first, second, &mut scalar)
                        .expect("short scalar end-gap alignment succeeds");
                    let observed = aligner
                        .align_with_workspace(first, second, &mut simd)
                        .expect("short SIMD end-gap alignment succeeds");
                    assert_eq!(observed, expected);
                }
            }
        }
    }
}

#[test]
fn one_simd_workspace_reuses_capacity_across_shapes_modes_and_alphabets() {
    if !simd_is_available() {
        return;
    }
    let scoring = MatchMismatch::new(3, -4);
    let gap = AffineGap::new(-5, -1).expect("valid affine gap");
    let mut workspace = workspace(AffineKernelPolicy::Simd);
    workspace
        .reserve_simd::<Protein>(64, 80)
        .expect("largest SIMD shape reserves successfully");
    let retained_capacities = (
        workspace.simd_final_score_row_capacity(),
        workspace.simd_score_row_capacity(),
        workspace.simd_profile_capacity(),
        workspace.simd_second_index_capacity(),
        workspace.simd_traceback_capacity(),
    );

    let dna_first = dna_from_indices(&[0, 1, 2, 3, 0, 1, 1, 2, 3, 0, 2, 1, 3]);
    let dna_second = dna_from_indices(&[3, 1, 2, 0, 1, 3, 2, 0, 1]);
    let global = Gotoh::new(scoring, gap);
    let expected_global = global
        .align(&dna_first, &dna_second)
        .expect("fresh global alignment succeeds");
    let retained_global = global
        .align_with_workspace(&dna_first, &dna_second, &mut workspace)
        .expect("workspace global alignment succeeds");
    assert_eq!(retained_global, expected_global);

    let protein_first =
        ProteinSequence::from_ascii(b"MTEYKLVVVGAGGVGKSALTIQLIQ").expect("valid protein");
    let protein_second =
        ProteinSequence::from_ascii(b"GGMTEYKLVVVGAGGVGKSA").expect("valid protein");
    let local = AffineLocalAligner::new(scoring, gap);
    assert_eq!(
        local
            .align_with_workspace(&protein_first, &protein_second, &mut workspace)
            .expect("workspace protein alignment succeeds"),
        local
            .align(&protein_first, &protein_second)
            .expect("fresh protein alignment succeeds")
    );

    let rna_first = RnaSequence::from_ascii(b"AUGCURYSWKMBDHVN").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"CCAUGCUARYSWKMBDHVNGG").expect("valid RNA");
    let end_gap = AffineEndGapAligner::new(scoring, gap, FreeEndGaps::fit_first_in_second());
    assert_eq!(
        end_gap
            .align_with_workspace(&rna_first, &rna_second, &mut workspace)
            .expect("workspace RNA alignment succeeds"),
        end_gap
            .align(&rna_first, &rna_second)
            .expect("fresh RNA alignment succeeds")
    );

    assert_eq!(
        retained_capacities,
        (
            workspace.simd_final_score_row_capacity(),
            workspace.simd_score_row_capacity(),
            workspace.simd_profile_capacity(),
            workspace.simd_second_index_capacity(),
            workspace.simd_traceback_capacity(),
        )
    );
    verify_global_affine_alignment(
        &dna_first,
        &dna_second,
        &retained_global,
        global.substitution(),
        global.gap(),
    )
    .expect("earlier owned result remains valid after workspace reuse");
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn simd_and_scalar_are_identical_across_scores_modes_and_boundaries(
        first_indices in prop::collection::vec(0_u8..4, 1..41),
        second_indices in prop::collection::vec(0_u8..4, 1..41),
        match_score in -2_i64..9,
        mismatch_score in -12_i64..4,
        raw_open in -14_i64..1,
        extension in -5_i64..1,
        free_end_bits in 0_u8..16,
    ) {
        if simd_is_available() {
            let first = dna_from_indices(&first_indices);
            let second = dna_from_indices(&second_indices);
            let scoring = MatchMismatch::new(match_score, mismatch_score);
            let gap = AffineGap::new(raw_open.min(extension), extension)
                .expect("generated affine gap is valid");
            let free_ends = FreeEndGaps::new(
                free_end_bits & 1 != 0,
                free_end_bits & 2 != 0,
                free_end_bits & 4 != 0,
                free_end_bits & 8 != 0,
            );

            let global = Gotoh::new(scoring, gap);
            let mut scalar_workspace = workspace(AffineKernelPolicy::Scalar);
            let expected_global = global
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar global alignment succeeds");
            let mut simd_workspace = workspace(AffineKernelPolicy::Simd);
            let observed_global = global
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD global alignment succeeds");
            prop_assert_eq!(&observed_global, &expected_global);
            let used_simd = is_simd_backend(simd_workspace.last_backend());
            prop_assert!(used_simd);
            verify_global_affine_alignment(
                &first,
                &second,
                &observed_global,
                &scoring,
                gap,
            )
            .expect("SIMD global result independently replays");

            let end_gap = AffineEndGapAligner::new(scoring, gap, free_ends);
            let expected_end_gap = end_gap
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar end-gap alignment succeeds");
            let observed_end_gap = end_gap
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD end-gap alignment succeeds");
            prop_assert_eq!(&observed_end_gap, &expected_end_gap);
            verify_affine_end_gap_alignment(
                &first,
                &second,
                &observed_end_gap,
                &scoring,
                gap,
                free_ends,
            )
            .expect("SIMD end-gap result independently replays");

            let local = AffineLocalAligner::new(scoring, gap);
            let expected_local = local
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar local alignment succeeds");
            let observed_local = local
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD local alignment succeeds");
            prop_assert_eq!(&observed_local, &expected_local);
            verify_local_affine_alignment(
                &first,
                &second,
                &observed_local,
                &scoring,
                gap,
            )
            .expect("SIMD local result independently replays");
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn signed_32_bit_simd_matches_scalar_across_vector_boundaries(
        first_indices in prop::collection::vec(0_u8..4, 1..35),
        second_indices in prop::collection::vec(0_u8..4, 1..35),
        match_score in 20_000_i64..30_001,
        mismatch_score in -20_000_i64..-99,
        raw_open in -15_000_i64..-4_999,
        extension in -5_000_i64..-99,
        free_end_bits in 0_u8..16,
    ) {
        if simd_is_available() {
            let first = dna_from_indices(&first_indices);
            let second = dna_from_indices(&second_indices);
            let scoring = MatchMismatch::new(match_score, mismatch_score);
            let gap = AffineGap::new(raw_open.min(extension), extension)
                .expect("generated affine gap is valid");
            let free_ends = FreeEndGaps::new(
                free_end_bits & 1 != 0,
                free_end_bits & 2 != 0,
                free_end_bits & 4 != 0,
                free_end_bits & 8 != 0,
            );
            let mut scalar_workspace = workspace(AffineKernelPolicy::Scalar);
            let mut simd_workspace = workspace(AffineKernelPolicy::Simd);

            let global = Gotoh::new(scoring, gap);
            let expected_global = global
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar signed-32-bit global alignment succeeds");
            let observed_global = global
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD signed-32-bit global alignment succeeds");
            prop_assert_eq!(&observed_global, &expected_global);
            prop_assert_eq!(
                simd_workspace.last_backend(),
                Some(AffineBackend::SimdI32 { lanes: 8 })
            );
            verify_global_affine_alignment(
                &first,
                &second,
                &observed_global,
                &scoring,
                gap,
            )
            .expect("signed-32-bit global result independently replays");

            let end_gap = AffineEndGapAligner::new(scoring, gap, free_ends);
            let expected_end_gap = end_gap
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar signed-32-bit end-gap alignment succeeds");
            let observed_end_gap = end_gap
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD signed-32-bit end-gap alignment succeeds");
            prop_assert_eq!(&observed_end_gap, &expected_end_gap);
            prop_assert_eq!(
                simd_workspace.last_backend(),
                Some(AffineBackend::SimdI32 { lanes: 8 })
            );
            verify_affine_end_gap_alignment(
                &first,
                &second,
                &observed_end_gap,
                &scoring,
                gap,
                free_ends,
            )
            .expect("signed-32-bit end-gap result independently replays");

            let local = AffineLocalAligner::new(scoring, gap);
            let expected_local = local
                .align_with_workspace(&first, &second, &mut scalar_workspace)
                .expect("scalar signed-32-bit local alignment succeeds");
            let observed_local = local
                .align_with_workspace(&first, &second, &mut simd_workspace)
                .expect("SIMD signed-32-bit local alignment succeeds");
            prop_assert_eq!(&observed_local, &expected_local);
            prop_assert_eq!(
                simd_workspace.last_backend(),
                Some(AffineBackend::SimdI32 { lanes: 8 })
            );
            verify_local_affine_alignment(
                &first,
                &second,
                &observed_local,
                &scoring,
                gap,
            )
            .expect("signed-32-bit local result independently replays");
        }
    }
}
