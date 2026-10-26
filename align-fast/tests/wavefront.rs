//! Exact gap-affine wavefront integration tests.

mod support;

use align_core::{AlignmentOperation, DnaSequence, ProteinSequence, RnaSequence};
use align_fast::{
    WavefrontAffine, WavefrontAffinePenalties, WavefrontError, WavefrontTraceRun,
    WavefrontWorkspace,
};
use support::wavefront::{binary_strings, dense_affine, exhaustive_affine, replay_affine};

fn dna(symbols: &[u8]) -> DnaSequence {
    DnaSequence::from_ascii(symbols).expect("test DNA is valid")
}

fn expand_runs(runs: impl Iterator<Item = WavefrontTraceRun>) -> Vec<AlignmentOperation> {
    runs.flat_map(|run| std::iter::repeat_n(run.operation(), run.count()))
        .collect()
}

fn deterministic_dna(length: usize) -> Vec<u8> {
    let bases = b"ACGT";
    (0..length)
        .map(|index| bases[(index * 7 + index / 19) % bases.len()])
        .collect()
}

fn insert_symbols(source: &[u8], position: usize, length: usize) -> Vec<u8> {
    let bases = b"TGCA";
    let mut result = Vec::with_capacity(source.len() + length);
    result.extend_from_slice(&source[..position]);
    result.extend((0..length).map(|index| bases[(index * 5 + index / 11) % bases.len()]));
    result.extend_from_slice(&source[position..]);
    result
}

fn assert_stress_case(first: &[u8], second: &[u8], expected: i64) {
    const MISMATCH: i64 = 4;
    const OPENING: i64 = 6;
    const EXTENSION: i64 = 2;

    assert_eq!(
        dense_affine(first, second, MISMATCH, OPENING, EXTENSION),
        expected,
        "dense oracle"
    );
    let first_sequence = dna(first);
    let second_sequence = dna(second);
    let aligner = WavefrontAffine::new(
        WavefrontAffinePenalties::new(MISMATCH, OPENING, EXTENSION)
            .expect("stress penalties are valid"),
    );
    let score = aligner
        .score(&first_sequence, &second_sequence)
        .expect("stress score succeeds");
    let owned = aligner
        .align(&first_sequence, &second_sequence)
        .expect("stress alignment succeeds");
    let mut workspace = WavefrontWorkspace::new();
    let compact = aligner
        .traceback_in_workspace(&first_sequence, &second_sequence, &mut workspace)
        .expect("stress compact traceback succeeds");
    let compact_penalty = compact.penalty();
    let compact_columns = compact.columns();
    let compact_run_count = compact.run_count();
    let compact_operations = expand_runs(compact.runs());

    assert_eq!(score.penalty(), expected);
    assert_eq!(owned.penalty(), expected);
    assert_eq!(compact_penalty, expected);
    assert_eq!(compact_columns, compact_operations.len());
    assert!(compact_run_count <= compact_columns);
    assert_eq!(compact_operations, owned.alignment().operations());
    assert_eq!(
        replay_affine(
            first,
            second,
            &compact_operations,
            MISMATCH,
            OPENING,
            EXTENSION,
        ),
        expected
    );
    let materialized = aligner
        .align_in_workspace(&first_sequence, &second_sequence, &mut workspace)
        .expect("stress materialized traceback succeeds");
    assert_eq!(materialized.penalty(), expected);
    assert_eq!(materialized.operations(), compact_operations);
    let materialized_summary = materialized.summary();
    let owned_summary = owned.summary();
    assert_eq!(
        materialized_summary.wavefronts(),
        owned_summary.wavefronts()
    );
    assert_eq!(
        materialized_summary.computed_offsets(),
        owned_summary.computed_offsets()
    );
    assert_eq!(
        materialized_summary.source_wavefront_lookups(),
        owned_summary.source_wavefront_lookups()
    );
    assert_eq!(
        materialized_summary.source_offset_lookups(),
        owned_summary.source_offset_lookups()
    );
    assert_eq!(
        materialized_summary.match_extensions(),
        owned_summary.match_extensions()
    );
    assert_eq!(
        materialized_summary.maximum_wavefront_width(),
        owned_summary.maximum_wavefront_width()
    );
}

#[test]
fn affine_penalties_follow_opening_plus_every_extension_symbol() {
    let penalties = WavefrontAffinePenalties::new(20, 3, 2).expect("valid penalties");
    assert_eq!(penalties.gap_penalty(0).expect("zero gap fits"), 0);
    assert_eq!(penalties.gap_penalty(1).expect("unit gap fits"), 5);
    assert_eq!(penalties.gap_penalty(3).expect("long gap fits"), 9);
    let aligner = WavefrontAffine::new(penalties);

    for (first, second, expected) in [
        (b"".as_slice(), b"".as_slice(), 0_i64),
        (b"".as_slice(), b"AAA".as_slice(), 9),
        (b"AAA".as_slice(), b"".as_slice(), 9),
        (b"ACGT".as_slice(), b"ACGT".as_slice(), 0),
        (b"ACGT".as_slice(), b"AGT".as_slice(), 5),
        (b"AC".as_slice(), b"TAC".as_slice(), 5),
        (b"AC".as_slice(), b"ACT".as_slice(), 5),
        (b"ACGT".as_slice(), b"AT".as_slice(), 7),
    ] {
        let first_sequence = dna(first);
        let second_sequence = dna(second);
        let score = aligner
            .score(&first_sequence, &second_sequence)
            .expect("affine score succeeds");
        let alignment = aligner
            .align(&first_sequence, &second_sequence)
            .expect("affine alignment succeeds");
        assert_eq!(score.penalty(), expected);
        assert_eq!(alignment.penalty(), expected);
        assert_eq!(alignment.alignment().score(), -expected);
        assert_eq!(
            replay_affine(
                first,
                second,
                alignment.alignment().operations(),
                penalties.mismatch(),
                penalties.gap_opening(),
                penalties.gap_extension(),
            ),
            expected
        );
    }
}

#[test]
fn published_gataca_example_reaches_penalty_eight() {
    let penalties = WavefrontAffinePenalties::new(4, 6, 2).expect("paper penalties are valid");
    let aligner = WavefrontAffine::new(penalties);
    let first = dna(b"GATACA");
    let second = dna(b"GAGATA");
    let score = aligner
        .score(&first, &second)
        .expect("paper score succeeds");
    let alignment = aligner
        .align(&first, &second)
        .expect("paper traceback succeeds");
    assert_eq!(score.penalty(), 8);
    assert_eq!(alignment.penalty(), 8);
    assert_eq!(alignment.alignment().score(), -8);
    assert_eq!(
        replay_affine(
            b"GATACA",
            b"GAGATA",
            alignment.alignment().operations(),
            4,
            6,
            2,
        ),
        8
    );
}

#[test]
fn documented_tie_prefers_substitution_over_two_unit_gaps() {
    let penalties = WavefrontAffinePenalties::new(2, 0, 1).expect("valid tie penalties");
    let alignment = WavefrontAffine::new(penalties)
        .align(&dna(b"A"), &dna(b"C"))
        .expect("tie alignment succeeds");
    assert_eq!(alignment.penalty(), 2);
    assert_eq!(
        alignment.alignment().operations(),
        &[AlignmentOperation::Substitution]
    );
}

#[test]
fn dna_rna_and_protein_share_literal_symbol_equality() {
    let affine =
        WavefrontAffine::new(WavefrontAffinePenalties::new(7, 5, 2).expect("valid penalties"));

    let dna_first = DnaSequence::from_ascii(b"ACGTN").expect("valid DNA");
    let dna_second = DnaSequence::from_ascii(b"ACCTN").expect("valid DNA");
    let rna_first = RnaSequence::from_ascii(b"ACGUN").expect("valid RNA");
    let rna_second = RnaSequence::from_ascii(b"ACCUN").expect("valid RNA");
    let protein_first = ProteinSequence::from_ascii(b"ACDEX").expect("valid protein");
    let protein_second = ProteinSequence::from_ascii(b"ACNEX").expect("valid protein");

    assert_eq!(affine.score(&dna_first, &dna_second).unwrap().penalty(), 7);
    assert_eq!(affine.score(&rna_first, &rna_second).unwrap().penalty(), 7);
    assert_eq!(
        affine
            .score(&protein_first, &protein_second)
            .unwrap()
            .penalty(),
        7
    );
}

#[test]
fn every_short_binary_affine_input_matches_independent_oracles() {
    let strings = binary_strings(3);
    for (mismatch, opening, extension) in [(1, 0, 1), (2, 0, 1), (4, 3, 1), (5, 2, 2)] {
        let penalties = WavefrontAffinePenalties::new(mismatch, opening, extension)
            .expect("oracle penalties are valid");
        let aligner = WavefrontAffine::new(penalties);
        let mut compact_workspace = WavefrontWorkspace::new();
        for first in &strings {
            for second in &strings {
                let expected = dense_affine(first, second, mismatch, opening, extension);
                let first_sequence = dna(first);
                let second_sequence = dna(second);
                let score = aligner
                    .score(&first_sequence, &second_sequence)
                    .expect("affine score succeeds");
                let penalty = aligner
                    .penalty(&first_sequence, &second_sequence)
                    .expect("lean affine penalty succeeds");
                let alignment = aligner
                    .align(&first_sequence, &second_sequence)
                    .expect("affine traceback succeeds");
                let compact = aligner
                    .traceback_in_workspace(
                        &first_sequence,
                        &second_sequence,
                        &mut compact_workspace,
                    )
                    .expect("compact affine traceback succeeds");
                let compact_columns = compact.columns();
                let compact_operations = expand_runs(compact.runs());
                assert_eq!(score.penalty(), expected, "{first:?} {second:?}");
                assert_eq!(penalty, expected, "{first:?} {second:?}");
                assert_eq!(alignment.penalty(), expected);
                assert_eq!(compact.penalty(), expected);
                assert_eq!(compact_columns, compact_operations.len());
                assert_eq!(
                    compact_operations,
                    alignment.alignment().operations(),
                    "{first:?} {second:?}"
                );
                assert_eq!(
                    replay_affine(
                        first,
                        second,
                        alignment.alignment().operations(),
                        mismatch,
                        opening,
                        extension,
                    ),
                    expected
                );
            }
        }
    }

    let tiny = binary_strings(2);
    for (mismatch, opening, extension) in [(2, 0, 1), (4, 3, 1), (5, 2, 2)] {
        for first in &tiny {
            for second in &tiny {
                assert_eq!(
                    dense_affine(first, second, mismatch, opening, extension),
                    exhaustive_affine(first, second, mismatch, opening, extension),
                    "dense oracle must independently agree with every tiny path"
                );
            }
        }
    }
}

#[test]
fn sparse_large_penalties_do_not_require_numeric_score_iteration() {
    let penalties = WavefrontAffinePenalties::new(1_000_003, 2_000_000, 999_983)
        .expect("large sparse penalties are valid");
    let aligner = WavefrontAffine::new(penalties);
    let mismatch = aligner
        .align(&dna(b"A"), &dna(b"C"))
        .expect("sparse mismatch succeeds");
    assert_eq!(mismatch.penalty(), 1_000_003);
    assert_eq!(mismatch.summary().wavefronts(), 2);

    let gap = aligner
        .align(&dna(b""), &dna(b"AAA"))
        .expect("sparse affine gap succeeds");
    assert_eq!(gap.penalty(), penalties.gap_penalty(3).unwrap());
    assert!(gap.summary().wavefronts() <= 4);
}

#[test]
fn scheduler_dispatch_boundaries_preserve_score_and_traceback() {
    let first = dna(b"A");
    let second = dna(b"C");

    for maximum_increment in [63_i64, 64, 4_096, 4_097] {
        let penalties = WavefrontAffinePenalties::new(maximum_increment, maximum_increment - 1, 1)
            .expect("scheduler-boundary penalties are valid");
        let aligner = WavefrontAffine::new(penalties);
        let penalty = aligner
            .penalty(&first, &second)
            .expect("scheduler-boundary lean penalty succeeds");
        let score = aligner
            .score(&first, &second)
            .expect("scheduler-boundary score succeeds");
        let alignment = aligner
            .align(&first, &second)
            .expect("scheduler-boundary traceback succeeds");

        assert_eq!(penalty, maximum_increment);
        assert_eq!(score.penalty(), maximum_increment);
        assert_eq!(alignment.penalty(), maximum_increment);
        assert_eq!(
            alignment.alignment().operations(),
            &[AlignmentOperation::Substitution]
        );
    }
}

#[test]
fn small_scheduler_history_index_wraps_and_reuses_safely() {
    let penalties = WavefrontAffinePenalties::new(1, 0, 1).expect("unit penalties are valid");
    let aligner = WavefrontAffine::new(penalties);
    let first_symbols = [b'A'; 96];
    let second_symbols = [b'C'; 96];
    let first = dna(&first_symbols);
    let second = dna(&second_symbols);
    let mut workspace = WavefrontWorkspace::new();

    let score = aligner
        .score_with_workspace(&first, &second, &mut workspace)
        .expect("score survives more than one history-index revolution");
    assert_eq!(score.penalty(), 96);
    assert!(score.wavefronts() > 64);

    {
        let alignment = aligner
            .align_in_workspace(&first, &second, &mut workspace)
            .expect("traceback survives more than one history-index revolution");
        assert_eq!(alignment.penalty(), 96);
        assert_eq!(
            replay_affine(
                &first_symbols,
                &second_symbols,
                alignment.operations(),
                penalties.mismatch(),
                penalties.gap_opening(),
                penalties.gap_extension(),
            ),
            96
        );
    }

    let short_first = dna(b"ACGT");
    let short_second = dna(b"AGGT");
    assert_eq!(
        aligner
            .score_with_workspace(&short_first, &short_second, &mut workspace)
            .expect("stale history-index slots are ignored")
            .penalty(),
        1
    );
}

#[test]
fn configuration_and_execution_statistics_are_explicit() {
    assert!(matches!(
        WavefrontAffinePenalties::new(0, 1, 1),
        Err(WavefrontError::NonPositiveMismatchPenalty { penalty: 0 })
    ));
    assert!(matches!(
        WavefrontAffinePenalties::new(-1, 1, 1),
        Err(WavefrontError::NonPositiveMismatchPenalty { penalty: -1 })
    ));
    assert!(matches!(
        WavefrontAffinePenalties::new(1, -1, 1),
        Err(WavefrontError::NegativeGapOpeningPenalty { penalty: -1 })
    ));
    assert!(matches!(
        WavefrontAffinePenalties::new(1, 0, 0),
        Err(WavefrontError::NonPositiveGapExtensionPenalty { penalty: 0 })
    ));
    assert!(matches!(
        WavefrontAffinePenalties::new(1, i64::MAX, 1),
        Err(WavefrontError::PenaltyOverflow { .. })
    ));
    let affine = WavefrontAffine::new(
        WavefrontAffinePenalties::new(4, 6, 2).expect("valid paper penalties"),
    );
    let affine_memory = affine
        .memory_requirement(7, 11, 8)
        .expect("small affine memory model fits");
    assert_eq!(affine_memory.maximum_wavefront_width(), 9);
    assert_eq!(affine_memory.maximum_score_wavefronts(), 6);
    assert_eq!(affine_memory.maximum_common_offsets(), 81);
    assert_eq!(affine_memory.maximum_state_offsets(), 243);
    assert!(matches!(
        affine.memory_requirement(7, 11, -1),
        Err(WavefrontError::NegativeTerminalPenalty { penalty: -1 })
    ));

    let identity = affine
        .score(&dna(b"ACGT"), &dna(b"ACGT"))
        .expect("identity succeeds");
    assert_eq!(identity.penalty(), 0);
    assert_eq!(identity.alignment_score(), 0);
    assert_eq!(identity.wavefronts(), 1);
    assert_eq!(identity.computed_offsets(), 1);
    assert_eq!(identity.source_wavefront_lookups(), 0);
    assert_eq!(identity.source_offset_lookups(), 0);
    assert_eq!(identity.match_extensions(), 4);
    assert_eq!(identity.maximum_wavefront_width(), 1);
    assert!(identity.modeled_peak_bytes() > 0);
}

#[test]
fn penalty_extremes_are_checked_without_saturation() {
    let penalties =
        WavefrontAffinePenalties::new(i64::MAX, 0, i64::MAX).expect("maximum deltas fit");
    assert!(penalties.gap_penalty(usize::MAX).is_err());
    let aligner = WavefrontAffine::new(penalties);
    let one_gap = aligner
        .align(&dna(b""), &dna(b"A"))
        .expect("one maximum-cost gap is representable");
    assert_eq!(one_gap.penalty(), i64::MAX);
    assert_eq!(one_gap.alignment().score(), -i64::MAX);
    assert!(aligner.score(&dna(b""), &dna(b"AA")).is_err());
}

#[test]
fn repeated_calls_are_deterministic_and_symmetric_in_penalty() {
    let first = dna(b"ACACACGT");
    let second = dna(b"CACACTGT");
    let affine =
        WavefrontAffine::new(WavefrontAffinePenalties::new(4, 6, 2).expect("valid penalties"));

    let affine_forward = affine.align(&first, &second).expect("affine succeeds");
    let affine_repeat = affine.align(&first, &second).expect("affine repeats");
    let affine_reverse = affine
        .align(&second, &first)
        .expect("reverse affine succeeds");
    assert_eq!(affine_forward, affine_repeat);
    assert_eq!(affine_forward.penalty(), affine_reverse.penalty());
}

#[test]
fn workspace_reuses_wavefronts_without_changing_results() {
    let affine =
        WavefrontAffine::new(WavefrontAffinePenalties::new(4, 6, 2).expect("valid penalties"));
    let first = dna(b"ACGTACGTACGTACGT");
    let second = dna(b"ACCTACGTTGGTACGT");
    let identity = dna(b"GATTACA");
    let mut workspace = WavefrontWorkspace::new();

    let expected_score = affine.score(&first, &second).expect("score succeeds");
    let expected_alignment = affine.align(&first, &second).expect("alignment succeeds");
    assert_eq!(
        affine
            .penalty_with_workspace(&first, &second, &mut workspace)
            .expect("workspace penalty succeeds"),
        expected_score.penalty()
    );
    assert_eq!(
        affine
            .score_with_workspace(&first, &second, &mut workspace)
            .expect("workspace score succeeds"),
        expected_score
    );
    assert_eq!(
        affine
            .align_with_workspace(&first, &second, &mut workspace)
            .expect("workspace alignment succeeds"),
        expected_alignment
    );
    let compact = affine
        .traceback_in_workspace(&first, &second, &mut workspace)
        .expect("workspace compact traceback succeeds");
    assert_eq!(compact.penalty(), expected_alignment.penalty());
    assert_eq!(compact.first_len(), first.len());
    assert_eq!(compact.second_len(), second.len());
    assert_eq!(compact.columns(), expected_alignment.alignment().columns());
    assert_eq!(
        expand_runs(compact.runs()),
        expected_alignment.alignment().operations()
    );
    let borrowed = affine
        .align_in_workspace(&first, &second, &mut workspace)
        .expect("borrowed workspace alignment succeeds");
    assert_eq!(borrowed.penalty(), expected_alignment.penalty());
    assert_eq!(
        borrowed.operations(),
        expected_alignment.alignment().operations()
    );
    assert_eq!(borrowed.first_len(), first.len());
    assert_eq!(borrowed.second_len(), second.len());
    let retained = workspace.retained_bytes();
    assert!(retained > 0);

    assert_eq!(
        affine
            .penalty_with_workspace(&identity, &identity, &mut workspace)
            .expect("identity workspace penalty succeeds"),
        0
    );
    assert_eq!(
        affine
            .score_with_workspace(&identity, &identity, &mut workspace)
            .expect("identity workspace score succeeds")
            .penalty(),
        0
    );
    assert_eq!(
        expand_runs(
            affine
                .traceback_in_workspace(&identity, &identity, &mut workspace)
                .expect("identity compact traceback succeeds")
                .runs()
        ),
        [AlignmentOperation::Match; 7]
    );
    assert_eq!(
        affine
            .align_in_workspace(&identity, &identity, &mut workspace)
            .expect("identity workspace alignment succeeds")
            .operations(),
        &[AlignmentOperation::Match; 7]
    );
    assert!(workspace.retained_bytes() >= retained);

    workspace.release();
    assert_eq!(workspace.retained_bytes(), 0);
}

#[test]
fn compact_traceback_covers_long_gaps_asymmetry_and_divergence() {
    let base = deterministic_dna(128);
    let penalties = WavefrontAffinePenalties::new(4, 6, 2).expect("stress penalties are valid");
    for gap_length in [1_usize, 8, 64, 512] {
        let inserted = insert_symbols(&base, 64, gap_length);
        let expected = penalties
            .gap_penalty(gap_length)
            .expect("stress gap penalty is representable");
        assert_stress_case(&base, &inserted, expected);
        assert_stress_case(&inserted, &base, expected);
    }

    let mut mixed = insert_symbols(&base, 32, 16);
    mixed.drain(96 + 16..96 + 16 + 8);
    for edit in 0..8 {
        let position = (edit + 1) * mixed.len() / 9;
        let symbol = mixed
            .get_mut(position)
            .expect("mixed stress substitution is in bounds");
        *symbol = match *symbol {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            _ => unreachable!("generator emits canonical DNA"),
        };
    }
    let mixed_expected = dense_affine(&base, &mixed, 4, 6, 2);
    assert_stress_case(&base, &mixed, mixed_expected);

    for percent in [1_usize, 5, 10, 20] {
        let first = deterministic_dna(256);
        let mut second = first.clone();
        let edits = (first.len() * percent / 100).max(1);
        for edit in 0..edits {
            let position = (edit + 1) * second.len() / (edits + 1);
            let symbol = second
                .get_mut(position)
                .expect("divergence stress substitution is in bounds");
            *symbol = match *symbol {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                _ => unreachable!("generator emits canonical DNA"),
            };
        }
        let expected = dense_affine(&first, &second, 4, 6, 2);
        assert_stress_case(&first, &second, expected);
    }
}
