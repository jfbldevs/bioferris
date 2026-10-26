//! Manual, exhaustive, block-boundary, and independent-oracle Myers tests.

use align_core::{
    AlignmentOperation, Dna, DnaSequence, LinearGap, MatchMismatch, PairwiseAlignment,
    ProteinSequence, RnaSequence,
};
use align_fast::{
    MyersPattern, MyersWorkspace, UnitEditVerificationError, verify_unit_edit_alignment,
};
use align_pairwise::NeedlemanWunsch;

fn reference_global(pattern: &[u8], text: &[u8]) -> usize {
    let mut previous: Vec<usize> = (0..=pattern.len()).collect();
    let mut current = vec![0; pattern.len() + 1];
    for (column, text_symbol) in text.iter().copied().enumerate() {
        current[0] = column + 1;
        for row in 1..=pattern.len() {
            current[row] = (previous[row - 1] + usize::from(pattern[row - 1] != text_symbol))
                .min(current[row - 1] + 1)
                .min(previous[row] + 1);
        }
        core::mem::swap(&mut previous, &mut current);
    }
    previous[pattern.len()]
}

fn reference_search(pattern: &[u8], text: &[u8]) -> Vec<usize> {
    let mut previous: Vec<usize> = (0..=pattern.len()).collect();
    let mut current = vec![0; pattern.len() + 1];
    let mut scores = Vec::with_capacity(text.len());
    for text_symbol in text.iter().copied() {
        current[0] = 0;
        for row in 1..=pattern.len() {
            current[row] = (previous[row - 1] + usize::from(pattern[row - 1] != text_symbol))
                .min(current[row - 1] + 1)
                .min(previous[row] + 1);
        }
        core::mem::swap(&mut previous, &mut current);
        scores.push(previous[pattern.len()]);
    }
    scores
}

fn binary_dna(max_len: usize) -> Vec<DnaSequence> {
    let mut sequences = Vec::new();
    for length in 0..=max_len {
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
            sequences.push(DnaSequence::from_ascii(bytes).expect("binary DNA is valid"));
        }
    }
    sequences
}

fn assert_alignment(pattern: &DnaSequence, text: &DnaSequence, expected: usize) {
    let myers = MyersPattern::new(pattern).expect("preprocessing succeeds");
    assert_eq!(
        myers.global_distance(text).expect("distance succeeds"),
        expected
    );
    let alignment = myers.align_global(text).expect("traceback succeeds");
    assert_eq!(
        alignment.score(),
        -i64::try_from(expected).expect("small distance")
    );
    let report =
        verify_unit_edit_alignment(pattern, text, &alignment).expect("independent replay succeeds");
    assert_eq!(report.distance(), expected);
}

#[test]
fn hand_calculated_empty_and_single_operation_cases() {
    let empty = DnaSequence::from_ascii(b"").expect("valid DNA");
    let acg = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    assert_alignment(&empty, &empty, 0);
    assert_alignment(&empty, &acg, 3);
    assert_alignment(&acg, &empty, 3);

    let insertions = MyersPattern::new(&empty)
        .expect("preprocessing succeeds")
        .align_global(&acg)
        .expect("traceback succeeds");
    assert_eq!(insertions.operations(), &[AlignmentOperation::Insertion; 3]);

    let deletions = MyersPattern::new(&acg)
        .expect("preprocessing succeeds")
        .align_global(&empty)
        .expect("traceback succeeds");
    assert_eq!(deletions.operations(), &[AlignmentOperation::Deletion; 3]);
}

#[test]
fn hand_calculated_internal_terminal_and_repetitive_cases() {
    for (pattern, text, distance) in [
        ("ACGT", "ACGT", 0),
        ("ACGT", "AGGT", 1),
        ("ACGT", "ACGGT", 1),
        ("ACGGT", "ACGT", 1),
        ("ACGT", "TACGT", 1),
        ("ACGT", "ACGTA", 1),
        ("AAAA", "TTTT", 4),
        ("AAAAAA", "AAA", 3),
        ("ACACACAC", "CACACACA", 2),
    ] {
        let pattern = DnaSequence::from_ascii(pattern).expect("valid DNA");
        let text = DnaSequence::from_ascii(text).expect("valid DNA");
        assert_alignment(&pattern, &text, distance);
    }
}

#[test]
fn classic_levenshtein_example_has_four_edits() {
    let pattern = ProteinSequence::from_ascii(b"GATTACA").expect("valid protein symbols");
    let text = ProteinSequence::from_ascii(b"GCATGCU").expect("valid protein symbols");
    let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
    assert_eq!(myers.global_distance(&text).expect("distance succeeds"), 4);
    let alignment = myers.align_global(&text).expect("traceback succeeds");
    assert_eq!(
        verify_unit_edit_alignment(&pattern, &text, &alignment)
            .expect("alignment verifies")
            .distance(),
        4
    );
}

#[test]
fn approximate_search_matches_hand_calculated_endpoints() {
    let pattern = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let text = DnaSequence::from_ascii(b"TTACGTACT").expect("valid DNA");
    let hits = MyersPattern::new(&pattern)
        .expect("preprocessing succeeds")
        .find_all(&text, 1)
        .expect("search succeeds");
    let observed: Vec<(usize, usize)> = hits
        .into_iter()
        .map(|hit| (hit.end(), hit.distance()))
        .collect();
    assert_eq!(observed, vec![(4, 1), (5, 0), (6, 1), (8, 1), (9, 1)]);
}

#[test]
fn empty_pattern_reports_only_scanned_endpoints() {
    let pattern = DnaSequence::from_ascii(b"").expect("valid DNA");
    let text = DnaSequence::from_ascii(b"ACG").expect("valid DNA");
    let hits = MyersPattern::new(&pattern)
        .expect("preprocessing succeeds")
        .find_all(&text, 0)
        .expect("search succeeds");
    assert_eq!(
        hits.iter()
            .map(|hit| (hit.end(), hit.distance()))
            .collect::<Vec<_>>(),
        vec![(1, 0), (2, 0), (3, 0)]
    );

    let empty_text = DnaSequence::from_ascii(b"").expect("valid DNA");
    assert!(
        MyersPattern::new(&pattern)
            .expect("preprocessing succeeds")
            .find_all(&empty_text, 0)
            .expect("search succeeds")
            .is_empty()
    );
}

#[test]
fn dna_rna_and_protein_use_literal_symbol_equality() {
    let dna_n = DnaSequence::from_ascii(b"N").expect("valid DNA");
    let dna_a = DnaSequence::from_ascii(b"A").expect("valid DNA");
    assert_eq!(
        MyersPattern::new(&dna_n)
            .expect("preprocessing succeeds")
            .global_distance(&dna_a)
            .expect("distance succeeds"),
        1
    );
    assert_eq!(
        MyersPattern::new(&dna_n)
            .expect("preprocessing succeeds")
            .global_distance(&dna_n)
            .expect("distance succeeds"),
        0
    );

    let rna = RnaSequence::from_ascii(b"ACGUUN").expect("valid RNA");
    let rna_other = RnaSequence::from_ascii(b"ACGUAN").expect("valid RNA");
    assert_eq!(
        MyersPattern::new(&rna)
            .expect("preprocessing succeeds")
            .global_distance(&rna_other)
            .expect("distance succeeds"),
        1
    );

    let protein = ProteinSequence::from_ascii(b"BZXJUO").expect("valid protein");
    let protein_other = ProteinSequence::from_ascii(b"BZAJUO").expect("valid protein");
    assert_eq!(
        MyersPattern::new(&protein)
            .expect("preprocessing succeeds")
            .global_distance(&protein_other)
            .expect("distance succeeds"),
        1
    );
}

#[test]
fn every_short_binary_input_matches_independent_matrices() {
    let sequences = binary_dna(4);
    for pattern in &sequences {
        let myers = MyersPattern::new(pattern).expect("preprocessing succeeds");
        for text in &sequences {
            let expected = reference_global(
                pattern.to_ascii_string().as_bytes(),
                text.to_ascii_string().as_bytes(),
            );
            assert_eq!(
                myers.global_distance(text).expect("distance succeeds"),
                expected
            );
            let alignment = myers.align_global(text).expect("traceback succeeds");
            assert_eq!(
                verify_unit_edit_alignment(pattern, text, &alignment)
                    .expect("traceback verifies")
                    .distance(),
                expected
            );

            let search_scores = reference_search(
                pattern.to_ascii_string().as_bytes(),
                text.to_ascii_string().as_bytes(),
            );
            let observed = myers.find_all(text, usize::MAX).expect("search succeeds");
            assert_eq!(
                observed
                    .iter()
                    .copied()
                    .map(align_fast::MyersMatch::distance)
                    .collect::<Vec<_>>(),
                search_scores
            );
        }
    }
}

#[test]
fn all_word_boundaries_match_independent_dynamic_programming() {
    for length in [63, 64, 65, 127, 128, 129] {
        let pattern_bytes: Vec<u8> = (0..length).map(|index| b"ACGT"[index % 4]).collect();
        let mut text_bytes = pattern_bytes.clone();
        text_bytes[length / 3] = b'T';
        text_bytes.insert(length / 2, b'G');
        text_bytes.remove((length * 3 / 4).min(text_bytes.len() - 1));

        let pattern = DnaSequence::from_ascii(&pattern_bytes).expect("valid DNA");
        let text = DnaSequence::from_ascii(&text_bytes).expect("valid DNA");
        let expected = reference_global(&pattern_bytes, &text_bytes);
        assert_alignment(&pattern, &text, expected);

        let expected_search = reference_search(&pattern_bytes, &text_bytes);
        let observed = MyersPattern::new(&pattern)
            .expect("preprocessing succeeds")
            .find_all(&text, usize::MAX)
            .expect("search succeeds");
        assert_eq!(
            observed
                .iter()
                .copied()
                .map(align_fast::MyersMatch::distance)
                .collect::<Vec<_>>(),
            expected_search
        );
    }
}

#[test]
fn long_banded_traceback_matches_independent_dynamic_programming() {
    for length in [512_usize, 513, 1_024] {
        let pattern_bytes: Vec<u8> = (0..length)
            .map(|index| b"ACGT"[(index * 11 + index / 17) % 4])
            .collect();
        let mut text_bytes = pattern_bytes.clone();
        for position in (13..text_bytes.len()).step_by(37) {
            let symbol = text_bytes
                .get_mut(position)
                .expect("generated mutation position is in bounds");
            *symbol = match *symbol {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        text_bytes.insert(length / 3, b'G');
        text_bytes.remove(length * 2 / 3);

        let pattern = DnaSequence::from_ascii(&pattern_bytes).expect("valid DNA");
        let text = DnaSequence::from_ascii(&text_bytes).expect("valid DNA");
        let expected = reference_global(&pattern_bytes, &text_bytes);
        assert_alignment(&pattern, &text, expected);
        let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
        assert_eq!(
            myers
                .global_distance_bounded(&text, expected)
                .expect("bounded distance succeeds"),
            Some(expected)
        );
        assert_eq!(
            myers
                .global_distance_bounded(&text, expected - 1)
                .expect("bounded rejection succeeds"),
            None
        );

        let identity = myers
            .align_global(&pattern)
            .expect("identity traceback succeeds");
        assert_eq!(identity.score(), 0);
        assert!(
            identity
                .operations()
                .iter()
                .all(|operation| *operation == AlignmentOperation::Match)
        );
    }
}

#[test]
fn scores_equal_independent_needleman_wunsch() {
    let scorer = MatchMismatch::new(0, -1);
    let gap = LinearGap::new(-1).expect("unit gap score is valid");
    let needleman_wunsch = NeedlemanWunsch::new(scorer, gap);
    for (pattern, text) in [
        ("", "ACGT"),
        ("ACGT", ""),
        ("ACGT", "AGT"),
        ("ACGTTGCA", "ACTTGGA"),
        ("NNNNACGT", "NANAACGA"),
    ] {
        let pattern = DnaSequence::from_ascii(pattern).expect("valid DNA");
        let text = DnaSequence::from_ascii(text).expect("valid DNA");
        let expected = needleman_wunsch
            .align(&pattern, &text)
            .expect("Needleman-Wunsch succeeds")
            .score();
        let observed = MyersPattern::new(&pattern)
            .expect("preprocessing succeeds")
            .align_global(&text)
            .expect("Myers succeeds")
            .score();
        assert_eq!(observed, expected);
    }
}

#[test]
fn tie_breaking_is_deterministic_on_repetitive_sequences() {
    let pattern = DnaSequence::from_ascii(b"AAAAAA").expect("valid DNA");
    let text = DnaSequence::from_ascii(b"AAAAA").expect("valid DNA");
    let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
    let first = myers.align_global(&text).expect("traceback succeeds");
    for _ in 0..20 {
        assert_eq!(
            myers.align_global(&text).expect("traceback succeeds"),
            first
        );
    }
}

#[test]
fn independent_verifier_rejects_an_incorrect_score() {
    let pattern = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let text = DnaSequence::from_ascii(b"A").expect("valid DNA");
    let invalid = PairwiseAlignment::try_new(-1, 0..1, 0..1, vec![AlignmentOperation::Match])
        .expect("structure is valid");
    assert_eq!(
        verify_unit_edit_alignment(&pattern, &text, &invalid),
        Err(UnitEditVerificationError::ScoreMismatch {
            reported: -1,
            recomputed: 0,
        })
    );
}

#[test]
fn memory_model_counts_blocks_and_rejects_overflow() {
    let estimate = MyersPattern::<Dna>::memory_requirement(65, 10).expect("dimensions fit");
    assert_eq!(estimate.pattern_symbols(), 65);
    assert_eq!(estimate.peq_words(), 30);
    assert_eq!(estimate.state_words(), 4);
    assert_eq!(estimate.band_score_records(), 0);
    assert_eq!(estimate.history_records(), 20);
    assert_eq!(estimate.history_column_records(), 0);
    assert_eq!(estimate.history_scratch_records(), 0);
    assert_eq!(estimate.operation_capacity(), 75);
    assert_eq!(estimate.search_match_capacity(), 10);
    assert_eq!(estimate.persistent_bytes(), 305);
    assert_eq!(estimate.distance_peak_bytes(), 337);
    assert_eq!(estimate.search_peak_bytes(), 497);
    assert_eq!(estimate.alignment_peak_bytes(), 732);

    let banded = MyersPattern::<Dna>::memory_requirement(512, 512).expect("dimensions fit");
    assert_eq!(banded.history_records(), 4_096);
    assert_eq!(banded.band_score_records(), 8);
    assert_eq!(banded.history_column_records(), 512);
    assert_eq!(banded.history_scratch_records(), 8);
    assert_eq!(banded.alignment_peak_bytes(), 76_544);

    assert!(MyersPattern::<Dna>::memory_requirement(usize::MAX, 0).is_err());
    assert!(MyersPattern::<Dna>::memory_requirement(usize::MAX, 1).is_err());
    assert!(MyersPattern::<Dna>::memory_requirement(1, usize::MAX).is_err());
}

#[test]
fn workspace_reuses_frontier_and_history_across_dimensions() {
    let long_pattern = DnaSequence::from_ascii("ACGT".repeat(128)).expect("valid DNA");
    let long_text = DnaSequence::from_ascii("ACGA".repeat(126)).expect("valid DNA");
    let short_pattern = DnaSequence::from_ascii(b"GATTACA").expect("valid DNA");
    let short_text = DnaSequence::from_ascii(
        b"GCATGCU".map(|symbol| if symbol == b'U' { b'T' } else { symbol }),
    )
    .expect("valid DNA");

    let long = MyersPattern::new(&long_pattern).expect("preprocessing succeeds");
    let short = MyersPattern::new(&short_pattern).expect("preprocessing succeeds");
    let mut workspace = MyersWorkspace::new();

    let expected_distance = long.global_distance(&long_text).expect("distance succeeds");
    assert_eq!(
        long.global_distance_with_workspace(&long_text, &mut workspace)
            .expect("workspace distance succeeds"),
        expected_distance
    );
    let expected_alignment = long.align_global(&long_text).expect("alignment succeeds");
    assert_eq!(
        long.align_global_with_workspace(&long_text, &mut workspace)
            .expect("workspace alignment succeeds"),
        expected_alignment
    );
    let retained_after_long = workspace.retained_bytes();
    assert!(retained_after_long > 0);

    assert_eq!(
        short
            .global_distance_with_workspace(&short_text, &mut workspace)
            .expect("short workspace distance succeeds"),
        short
            .global_distance(&short_text)
            .expect("short distance succeeds")
    );
    assert_eq!(
        short
            .find_all_with_workspace(&short_text, 3, &mut workspace)
            .expect("workspace search succeeds"),
        short.find_all(&short_text, 3).expect("search succeeds")
    );
    assert_eq!(workspace.retained_bytes(), retained_after_long);

    workspace.release();
    assert_eq!(workspace.retained_bytes(), 0);
}
