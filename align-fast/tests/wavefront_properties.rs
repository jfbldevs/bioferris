//! Property tests for exact gap-affine wavefront alignment.

mod support;

use align_core::{AlignmentOperation, DnaSequence};
use align_fast::{
    WavefrontAffine, WavefrontAffinePenalties, WavefrontScore, WavefrontTraceRun,
    WavefrontWorkspace,
};
use proptest::prelude::*;
use support::wavefront::{dense_affine, replay_affine};

fn dna(indices: &[u8]) -> (Vec<u8>, DnaSequence) {
    let bytes: Vec<u8> = indices
        .iter()
        .map(|index| b"ACGT"[usize::from(*index) % 4])
        .collect();
    let sequence = DnaSequence::from_ascii(&bytes).expect("generator emits canonical DNA");
    (bytes, sequence)
}

fn same_work(first: WavefrontScore, second: WavefrontScore) -> bool {
    first.penalty() == second.penalty()
        && first.wavefronts() == second.wavefronts()
        && first.computed_offsets() == second.computed_offsets()
        && first.match_extensions() == second.match_extensions()
        && first.maximum_wavefront_width() == second.maximum_wavefront_width()
}

fn expand_runs(runs: impl Iterator<Item = WavefrontTraceRun>) -> Vec<AlignmentOperation> {
    runs.flat_map(|run| std::iter::repeat_n(run.operation(), run.count()))
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn affine_wavefront_matches_dense_oracle_and_replay(
        first_indices in prop::collection::vec(0_u8..4, 0..24),
        second_indices in prop::collection::vec(0_u8..4, 0..24),
        mismatch in 1_i64..9,
        opening in 0_i64..9,
        extension in 1_i64..6,
    ) {
        let (first_bytes, first) = dna(&first_indices);
        let (second_bytes, second) = dna(&second_indices);
        let expected = dense_affine(
            &first_bytes,
            &second_bytes,
            mismatch,
            opening,
            extension,
        );
        let penalties = WavefrontAffinePenalties::new(mismatch, opening, extension)
            .expect("generated penalties are valid");
        let aligner = WavefrontAffine::new(penalties);
        let penalty = aligner
            .penalty(&first, &second)
            .expect("lean affine penalty succeeds");
        let score = aligner.score(&first, &second).expect("affine score succeeds");
        let alignment = aligner
            .align(&first, &second)
            .expect("affine traceback succeeds");
        let mut workspace = WavefrontWorkspace::new();
        let compact = aligner
            .traceback_in_workspace(&first, &second, &mut workspace)
            .expect("compact affine traceback succeeds");
        let compact_penalty = compact.penalty();
        let compact_columns = compact.columns();
        let compact_operations = expand_runs(compact.runs());

        prop_assert_eq!(penalty, expected);
        prop_assert_eq!(score.penalty(), expected);
        prop_assert_eq!(alignment.penalty(), expected);
        prop_assert_eq!(compact_penalty, expected);
        prop_assert!(same_work(score, alignment.summary()));
        prop_assert_eq!(compact_columns, compact_operations.len());
        prop_assert_eq!(compact_operations, alignment.alignment().operations());
        prop_assert_eq!(
            replay_affine(
                &first_bytes,
                &second_bytes,
                alignment.alignment().operations(),
                mismatch,
                opening,
                extension,
            ),
            expected
        );
        prop_assert_eq!(
            aligner.align(&first, &second).expect("repeat succeeds"),
            alignment
        );
    }

    #[test]
    fn symmetric_affine_penalties_are_symmetric(
        first_indices in prop::collection::vec(0_u8..4, 0..40),
        second_indices in prop::collection::vec(0_u8..4, 0..40),
        mismatch in 1_i64..9,
        opening in 0_i64..9,
        extension in 1_i64..6,
    ) {
        let (_, first) = dna(&first_indices);
        let (_, second) = dna(&second_indices);
        let affine = WavefrontAffine::new(
            WavefrontAffinePenalties::new(mismatch, opening, extension)
                .expect("generated penalties are valid"),
        );
        prop_assert_eq!(
            affine
                .score(&first, &second)
                .expect("forward affine succeeds")
                .penalty(),
            affine
                .score(&second, &first)
                .expect("reverse affine succeeds")
                .penalty(),
        );
    }

}
