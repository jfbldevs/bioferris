//! Acceptance tests for exact affine score-only execution.

use core::mem::size_of;
use std::sync::OnceLock;

use align_core::{
    AffineGap, Alphabet, DenseSubstitutionMatrix, Dna, DnaSequence, MatchMismatch, Protein,
    ProteinSequence, RnaSequence, Score, Sequence, SubstitutionScorer, Symbol,
};
use align_pairwise::{
    AffineBackend, AffineEndGapAligner, AffineKernelPolicy, AffineLocalAligner, AffineWorkspace,
    AlignmentError, FreeEndGaps, Gotoh,
};
use proptest::prelude::*;

fn gap() -> AffineGap {
    AffineGap::new(-5, -1).expect("valid affine gap")
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
        let sequence = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let mut workspace = workspace(AffineKernelPolicy::Simd);
        match Gotoh::new(MatchMismatch::new(2, -3), gap()).score_with_workspace(
            &sequence,
            &sequence,
            &mut workspace,
        ) {
            Ok(_) => is_simd_backend(workspace.last_backend()),
            Err(AlignmentError::SimdUnavailable { .. }) => false,
            Err(error) => panic!("SIMD capability probe failed unexpectedly: {error}"),
        }
    })
}

fn binary_sequences(max_len: usize) -> Vec<DnaSequence> {
    let mut sequences = Vec::new();
    for length in 0..=max_len {
        for bits in 0..(1_usize << length) {
            let ascii = (0..length)
                .map(|offset| {
                    if bits & (1 << offset) == 0 {
                        b'A'
                    } else {
                        b'C'
                    }
                })
                .collect::<Vec<_>>();
            sequences.push(Sequence::from_ascii(ascii).expect("binary DNA is valid"));
        }
    }
    sequences
}

fn repeated_dna(symbol: u8, length: usize) -> DnaSequence {
    Sequence::from_ascii(vec![symbol; length]).expect("repeated canonical DNA is valid")
}

#[derive(Clone, Debug)]
struct CachedMatrix<A: Alphabet>(DenseSubstitutionMatrix<A>);

impl<A: Alphabet> SubstitutionScorer<A> for CachedMatrix<A> {
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

#[derive(Clone, Copy, Debug)]
struct MaximumMatch;

impl SubstitutionScorer<Dna> for MaximumMatch {
    fn score(&self, _first: Symbol<Dna>, _second: Symbol<Dna>) -> Score {
        Score::MAX
    }
}

fn assert_score_profiles_match<A, S>(first: &Sequence<A>, second: &Sequence<A>, substitution: &S)
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let mut alignment_workspace = workspace(AffineKernelPolicy::Scalar);
    let mut score_workspace = workspace(AffineKernelPolicy::Scalar);
    let global = Gotoh::new(substitution, gap());
    assert_eq!(
        global
            .score_with_workspace(first, second, &mut score_workspace)
            .expect("global score succeeds"),
        global
            .align_with_workspace(first, second, &mut alignment_workspace)
            .expect("global alignment succeeds")
            .score()
    );

    let local = AffineLocalAligner::new(substitution, gap());
    assert_eq!(
        local
            .score_with_workspace(first, second, &mut score_workspace)
            .expect("local score succeeds"),
        local
            .align_with_workspace(first, second, &mut alignment_workspace)
            .expect("local alignment succeeds")
            .score()
    );

    for flags in 0_u8..16 {
        let free_ends = FreeEndGaps::new(
            flags & 1 != 0,
            flags & 2 != 0,
            flags & 4 != 0,
            flags & 8 != 0,
        );
        let aligner = AffineEndGapAligner::new(substitution, gap(), free_ends);
        assert_eq!(
            aligner
                .score_with_workspace(first, second, &mut score_workspace)
                .expect("free-end score succeeds"),
            aligner
                .align_with_workspace(first, second, &mut alignment_workspace)
                .expect("free-end alignment succeeds")
                .score(),
            "free-end flags {flags:04b}"
        );
    }
    assert_eq!(score_workspace.last_backend(), Some(AffineBackend::Scalar));
}

#[test]
fn manual_global_local_free_end_and_empty_scores_are_exact() {
    let scoring = MatchMismatch::new(2, -3);
    let first = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"ACCT").expect("valid DNA");
    assert_eq!(
        Gotoh::new(scoring, gap()).score(&first, &second).unwrap(),
        3
    );

    let first = DnaSequence::from_ascii(b"TTACGTAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"GGACGTCC").expect("valid DNA");
    assert_eq!(
        AffineLocalAligner::new(scoring, gap())
            .score(&first, &second)
            .unwrap(),
        8
    );

    let first = DnaSequence::from_ascii(b"AAA").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"TTAAAGG").expect("valid DNA");
    assert_eq!(
        AffineEndGapAligner::new(scoring, gap(), FreeEndGaps::fit_first_in_second())
            .score(&first, &second)
            .unwrap(),
        6
    );

    let empty = DnaSequence::from_ascii(b"").expect("valid empty DNA");
    assert_eq!(
        Gotoh::new(scoring, gap()).score(&empty, &first).unwrap(),
        -7
    );
    assert_eq!(
        AffineLocalAligner::new(scoring, gap())
            .score(&empty, &first)
            .unwrap(),
        0
    );
    assert_eq!(
        AffineEndGapAligner::new(scoring, gap(), FreeEndGaps::new(false, false, false, true),)
            .score(&empty, &first)
            .unwrap(),
        0
    );
}

#[test]
fn typed_alphabets_and_cached_matrices_match_complete_alignment() {
    let dna_first = DnaSequence::from_ascii(b"ACGTNRYSACGT").expect("valid DNA");
    let dna_second = DnaSequence::from_ascii(b"AGGTNRYTACCT").expect("valid DNA");
    assert_score_profiles_match(&dna_first, &dna_second, &MatchMismatch::new(3, -4));

    let rna_first = RnaSequence::from_ascii(b"AUGCUNRYAUGC").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"AUGAUNRYUUCC").expect("valid RNA");
    assert_score_profiles_match(&rna_first, &rna_second, &MatchMismatch::new(4, -5));

    let protein_first = ProteinSequence::from_ascii(b"MTEYKLVVVGAGGVGKSAL").expect("protein");
    let protein_second = ProteinSequence::from_ascii(b"MTEFQLVVVGAAGVGKTAL").expect("protein");
    let matrix = CachedMatrix(DenseSubstitutionMatrix::<Protein>::from_score_fn(
        |first, second| if first == second { 7 } else { -5 },
    ));
    assert_score_profiles_match(&protein_first, &protein_second, &matrix);
}

#[test]
fn exhaustive_short_inputs_match_complete_scalar_and_strict_simd() {
    let sequences = binary_sequences(3);
    let scoring_profiles = [MatchMismatch::new(2, -3), MatchMismatch::new(0, 0)];
    let has_simd = simd_is_available();
    let mut alignment_workspace = workspace(AffineKernelPolicy::Scalar);
    let mut scalar_workspace = workspace(AffineKernelPolicy::Scalar);
    let mut simd_workspace = workspace(AffineKernelPolicy::Simd);

    for scoring in scoring_profiles {
        for first in &sequences {
            for second in &sequences {
                let global = Gotoh::new(scoring, gap());
                let expected = global
                    .align_with_workspace(first, second, &mut alignment_workspace)
                    .expect("scalar global alignment")
                    .score();
                assert_eq!(
                    global
                        .score_with_workspace(first, second, &mut scalar_workspace)
                        .expect("scalar global score"),
                    expected
                );
                if has_simd && !first.is_empty() && !second.is_empty() {
                    assert_eq!(
                        global
                            .score_with_workspace(first, second, &mut simd_workspace)
                            .expect("SIMD global score"),
                        expected
                    );
                }

                let local = AffineLocalAligner::new(scoring, gap());
                let expected = local
                    .align_with_workspace(first, second, &mut alignment_workspace)
                    .expect("scalar local alignment")
                    .score();
                assert_eq!(
                    local
                        .score_with_workspace(first, second, &mut scalar_workspace)
                        .expect("scalar local score"),
                    expected
                );
                if has_simd && !first.is_empty() && !second.is_empty() {
                    assert_eq!(
                        local
                            .score_with_workspace(first, second, &mut simd_workspace)
                            .expect("SIMD local score"),
                        expected
                    );
                }

                for flags in 0_u8..16 {
                    let free_ends = FreeEndGaps::new(
                        flags & 1 != 0,
                        flags & 2 != 0,
                        flags & 4 != 0,
                        flags & 8 != 0,
                    );
                    let aligner = AffineEndGapAligner::new(scoring, gap(), free_ends);
                    let expected = aligner
                        .align_with_workspace(first, second, &mut alignment_workspace)
                        .expect("scalar free-end alignment")
                        .score();
                    assert_eq!(
                        aligner
                            .score_with_workspace(first, second, &mut scalar_workspace)
                            .expect("scalar free-end score"),
                        expected,
                        "flags {flags:04b}, first {:?}, second {:?}",
                        first.as_slice(),
                        second.as_slice()
                    );
                    if has_simd && !first.is_empty() && !second.is_empty() {
                        assert_eq!(
                            aligner
                                .score_with_workspace(first, second, &mut simd_workspace)
                                .expect("SIMD free-end score"),
                            expected,
                            "SIMD flags {flags:04b}, first {:?}, second {:?}",
                            first.as_slice(),
                            second.as_slice()
                        );
                    }
                }
            }
        }
    }
    assert_eq!(scalar_workspace.traceback_capacity(), 0);
    assert_eq!(scalar_workspace.simd_traceback_capacity(), 0);
    if has_simd {
        assert_eq!(simd_workspace.traceback_capacity(), 0);
        assert_eq!(simd_workspace.simd_traceback_capacity(), 0);
    }
}

#[test]
fn strict_simd_uses_both_exact_lane_widths_without_traceback() {
    if !simd_is_available() {
        return;
    }
    let sequence = DnaSequence::from_ascii(b"ACGTACGTACGTACGT").expect("valid DNA");
    let mut workspace = workspace(AffineKernelPolicy::Simd);
    let narrow = Gotoh::new(MatchMismatch::new(100, -100), gap());
    assert_eq!(
        narrow
            .score_with_workspace(&sequence, &sequence, &mut workspace)
            .unwrap(),
        1_600
    );
    assert_eq!(
        workspace.last_backend(),
        Some(AffineBackend::SimdI16 { lanes: 16 })
    );

    let wide = Gotoh::new(MatchMismatch::new(50_000, -50_000), gap());
    assert_eq!(
        wide.score_with_workspace(&sequence, &sequence, &mut workspace)
            .unwrap(),
        800_000
    );
    assert_eq!(
        workspace.last_backend(),
        Some(AffineBackend::SimdI32 { lanes: 8 })
    );
    assert_eq!(workspace.traceback_capacity(), 0);
    assert_eq!(workspace.simd_traceback_capacity(), 0);
}

#[test]
fn automatic_score_policy_obeys_its_calibrated_cell_gate() {
    let aligner = Gotoh::new(MatchMismatch::new(2, -3), gap());
    let mut automatic = workspace(AffineKernelPolicy::Auto);

    for (first_len, second_len, simd_eligible) in [(63, 65, false), (64, 64, true)] {
        aligner
            .score_with_workspace(
                &repeated_dna(b'A', first_len),
                &repeated_dna(b'C', second_len),
                &mut automatic,
            )
            .expect("automatic score succeeds");
        if simd_eligible && simd_is_available() {
            assert!(is_simd_backend(automatic.last_backend()));
        } else {
            assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));
        }
    }
}

#[test]
fn automatic_score_policy_obeys_its_calibrated_i16_shape_gate() {
    let aligner = Gotoh::new(MatchMismatch::new(2, -3), gap());
    let mut automatic = workspace(AffineKernelPolicy::Auto);

    for (first_len, second_len, simd_eligible) in [
        (8, 512, false),
        (9, 456, false),
        (10, 410, true),
        (11, 373, true),
        (16, 256, true),
    ] {
        aligner
            .score_with_workspace(
                &repeated_dna(b'A', first_len),
                &repeated_dna(b'C', second_len),
                &mut automatic,
            )
            .expect("automatic signed-16-bit score succeeds");
        if simd_eligible && simd_is_available() {
            assert_eq!(
                automatic.last_backend(),
                Some(AffineBackend::SimdI16 { lanes: 16 })
            );
        } else {
            assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));
        }
    }
}

#[test]
fn automatic_score_policy_obeys_its_calibrated_i32_shape_gate() {
    let aligner = Gotoh::new(MatchMismatch::new(10_000, -4), gap());
    let mut automatic = workspace(AffineKernelPolicy::Auto);

    for (first_len, second_len, simd_eligible) in [
        (7, 586, false),
        (8, 512, true),
        (9, 456, false),
        (10, 410, false),
        (11, 373, true),
        (16, 256, true),
    ] {
        aligner
            .score_with_workspace(
                &repeated_dna(b'A', first_len),
                &repeated_dna(b'C', second_len),
                &mut automatic,
            )
            .expect("automatic signed-32-bit score succeeds");
        if simd_eligible && simd_is_available() {
            assert_eq!(
                automatic.last_backend(),
                Some(AffineBackend::SimdI32 { lanes: 8 })
            );
        } else {
            assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));
        }
    }
}

#[test]
fn unknown_bounds_fall_back_and_overflow_matches_complete_alignment() {
    let first = DnaSequence::from_ascii(b"ACGTACGT").expect("valid DNA");
    let second = DnaSequence::from_ascii(b"AGGTTCGT").expect("valid DNA");
    let mut automatic = workspace(AffineKernelPolicy::Auto);
    let aligner = Gotoh::new(UnknownBounds, gap());
    let score = aligner
        .score_with_workspace(&first, &second, &mut automatic)
        .expect("checked fallback succeeds");
    assert_eq!(score, aligner.align(&first, &second).unwrap().score());
    assert_eq!(automatic.last_backend(), Some(AffineBackend::Scalar));

    let two = DnaSequence::from_ascii(b"AA").expect("valid DNA");
    let overflowing = Gotoh::new(MaximumMatch, AffineGap::new(0, 0).expect("valid zero gap"));
    assert!(matches!(
        overflowing.score(&two, &two),
        Err(AlignmentError::ScoreOverflow { row: 2, column: 2 })
    ));
    assert!(matches!(
        overflowing.align(&two, &two),
        Err(AlignmentError::ScoreOverflow { row: 2, column: 2 })
    ));
}

#[test]
fn score_only_memory_and_reservation_exclude_traceback() {
    let scalar =
        Gotoh::<MatchMismatch>::score_memory_requirement(257, 509).expect("scalar dimensions fit");
    assert_eq!(scalar.columns(), 510);
    assert_eq!(scalar.score_buffer_bytes(), 510 * 2 * size_of::<Score>());
    assert_eq!(scalar.state_validity_buffer_bytes(), 510);
    assert_eq!(
        scalar.checked_peak_buffer_bytes(),
        scalar.score_buffer_bytes() + scalar.state_validity_buffer_bytes()
    );

    let simd = AffineWorkspace::simd_score_memory_requirement::<Dna>(257, 509)
        .expect("SIMD dimensions fit");
    assert_eq!(
        simd.peak_buffer_bytes(),
        simd.final_score_row_bytes()
            + simd.score_column_bytes()
            + simd.substitution_profile_bytes()
            + simd.second_index_bytes()
    );

    let mut workspace = AffineWorkspace::new();
    workspace.reserve_score(257, 509).unwrap();
    workspace.reserve_simd_score::<Dna>(257, 509).unwrap();
    assert!(workspace.score_row_capacity() >= 510);
    assert!(workspace.simd_score_row_capacity() >= 265);
    assert_eq!(workspace.traceback_capacity(), 0);
    assert_eq!(workspace.simd_traceback_capacity(), 0);
}

fn dna_strategy(max_length: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(0_u8..4, 0..=max_length)
}

fn dna_from_indices(indices: &[u8]) -> DnaSequence {
    let ascii = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index)])
        .collect::<Vec<_>>();
    Sequence::from_ascii(ascii).expect("generated DNA is valid")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn random_scalar_and_simd_scores_equal_complete_alignment(
        first_raw in dna_strategy(48),
        second_raw in dna_strategy(48),
        flags in 0_u8..16,
    ) {
        let first = dna_from_indices(&first_raw);
        let second = dna_from_indices(&second_raw);
        let scoring = MatchMismatch::new(4, -3);
        let free_ends = FreeEndGaps::new(
            flags & 1 != 0,
            flags & 2 != 0,
            flags & 4 != 0,
            flags & 8 != 0,
        );
        let mut alignment_workspace = workspace(AffineKernelPolicy::Scalar);
        let mut scalar_workspace = workspace(AffineKernelPolicy::Scalar);
        let mut simd_workspace = workspace(AffineKernelPolicy::Simd);

        let global = Gotoh::new(scoring, gap());
        let global_expected = global
            .align_with_workspace(&first, &second, &mut alignment_workspace)
            .expect("global alignment")
            .score();
        prop_assert_eq!(
            global.score_with_workspace(&first, &second, &mut scalar_workspace).unwrap(),
            global_expected
        );

        let local = AffineLocalAligner::new(scoring, gap());
        let local_expected = local
            .align_with_workspace(&first, &second, &mut alignment_workspace)
            .expect("local alignment")
            .score();
        prop_assert_eq!(
            local.score_with_workspace(&first, &second, &mut scalar_workspace).unwrap(),
            local_expected
        );

        let end_gap = AffineEndGapAligner::new(scoring, gap(), free_ends);
        let end_gap_expected = end_gap
            .align_with_workspace(&first, &second, &mut alignment_workspace)
            .expect("free-end alignment")
            .score();
        prop_assert_eq!(
            end_gap.score_with_workspace(&first, &second, &mut scalar_workspace).unwrap(),
            end_gap_expected
        );

        if simd_is_available() && !first.is_empty() && !second.is_empty() {
            prop_assert_eq!(
                global.score_with_workspace(&first, &second, &mut simd_workspace).unwrap(),
                global_expected
            );
            prop_assert_eq!(
                local.score_with_workspace(&first, &second, &mut simd_workspace).unwrap(),
                local_expected
            );
            prop_assert_eq!(
                end_gap.score_with_workspace(&first, &second, &mut simd_workspace).unwrap(),
                end_gap_expected
            );
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn signed_32_bit_score_kernel_matches_scalar(
        first_raw in prop::collection::vec(0_u8..4, 1..=36),
        second_raw in prop::collection::vec(0_u8..4, 1..=36),
        local in any::<bool>(),
    ) {
        if !simd_is_available() {
            return Ok(());
        }
        let first = dna_from_indices(&first_raw);
        let second = dna_from_indices(&second_raw);
        let scoring = MatchMismatch::new(50_000, -40_000);
        let wide_gap = AffineGap::new(-60_000, -20_000).expect("valid wide gap");
        let mut scalar = workspace(AffineKernelPolicy::Scalar);
        let mut simd = workspace(AffineKernelPolicy::Simd);
        let (expected, observed) = if local {
            let aligner = AffineLocalAligner::new(scoring, wide_gap);
            (
                aligner.score_with_workspace(&first, &second, &mut scalar).unwrap(),
                aligner.score_with_workspace(&first, &second, &mut simd).unwrap(),
            )
        } else {
            let aligner = Gotoh::new(scoring, wide_gap);
            (
                aligner.score_with_workspace(&first, &second, &mut scalar).unwrap(),
                aligner.score_with_workspace(&first, &second, &mut simd).unwrap(),
            )
        };
        prop_assert_eq!(observed, expected);
        prop_assert_eq!(
            simd.last_backend(),
            Some(AffineBackend::SimdI32 { lanes: 8 })
        );
    }
}
