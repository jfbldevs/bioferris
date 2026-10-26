//! Hand-calculated tests for reference SP and TC metrics.

use align_core::{Dna, Symbol};
use align_msa::{
    AlignedRow, MsaError, MultipleAlignment, compare_alignments,
    compare_alignments_on_reference_columns,
};

fn row(identifier: &str, aligned: &str) -> AlignedRow<Dna> {
    AlignedRow::new(
        identifier,
        aligned
            .bytes()
            .map(|byte| {
                (byte != b'-').then(|| Symbol::<Dna>::from_ascii(byte).expect("valid test symbol"))
            })
            .collect(),
    )
    .unwrap()
}

fn alignment(rows: &[(&str, &str)]) -> MultipleAlignment<Dna> {
    MultipleAlignment::try_new(
        rows.iter()
            .map(|(identifier, symbols)| row(identifier, symbols))
            .collect(),
    )
    .unwrap()
}

#[test]
fn identical_alignment_scores_one_even_when_rows_are_reordered() {
    let reference = alignment(&[("a", "A-C"), ("b", "ATC"), ("c", "A-C")]);
    let candidate = alignment(&[("c", "A-C"), ("a", "A-C"), ("b", "ATC")]);
    let metrics = compare_alignments(&candidate, &reference).unwrap();
    assert!((metrics.sp_score() - 1.0).abs() < f64::EPSILON);
    assert!((metrics.tc_score() - 1.0).abs() < f64::EPSILON);
    assert_eq!(metrics.correct_pairs(), metrics.reference_pairs());
    assert_eq!(metrics.correct_columns(), metrics.reference_columns());
}

#[test]
fn shifted_gap_has_hand_calculated_sp_and_tc_scores() {
    let reference = alignment(&[("a", "A-C"), ("b", "ATC")]);
    let candidate = alignment(&[("a", "AC-"), ("b", "ATC")]);
    let metrics = compare_alignments(&candidate, &reference).unwrap();
    assert_eq!(metrics.correct_pairs(), 1);
    assert_eq!(metrics.reference_pairs(), 2);
    assert!((metrics.sp_score() - 0.5).abs() < f64::EPSILON);
    assert_eq!(metrics.correct_columns(), 1);
    assert_eq!(metrics.reference_columns(), 3);
    assert!((metrics.tc_score() - 1.0 / 3.0).abs() < f64::EPSILON);
}

#[test]
fn alignments_without_residue_pairs_have_defined_perfect_empty_score() {
    let reference = alignment(&[("empty", "")]);
    let candidate = alignment(&[("empty", "")]);
    let metrics = compare_alignments(&candidate, &reference).unwrap();
    assert!((metrics.sp_score() - 1.0).abs() < f64::EPSILON);
    assert!((metrics.tc_score() - 1.0).abs() < f64::EPSILON);
}

#[test]
fn metrics_reject_different_sequences_and_identifiers() {
    let reference = alignment(&[("a", "AC"), ("b", "AT")]);
    let different_symbols = alignment(&[("a", "AG"), ("b", "AT")]);
    assert_eq!(
        compare_alignments(&different_symbols, &reference),
        Err(MsaError::UngappedSequenceMismatch {
            identifier: "a".to_owned(),
        })
    );
    let different_identifier = alignment(&[("a", "AC"), ("c", "AT")]);
    assert_eq!(
        compare_alignments(&different_identifier, &reference),
        Err(MsaError::SequenceSetMismatch)
    );
}

#[test]
fn reference_column_mask_restricts_both_sp_and_tc_denominators() {
    let reference = alignment(&[("a", "A-C"), ("b", "ATC")]);
    let candidate = alignment(&[("a", "AC-"), ("b", "ATC")]);
    let metrics =
        compare_alignments_on_reference_columns(&candidate, &reference, &[true, false, true])
            .unwrap();
    assert_eq!(metrics.correct_pairs(), 1);
    assert_eq!(metrics.reference_pairs(), 2);
    assert!((metrics.sp_score() - 0.5).abs() < f64::EPSILON);
    assert_eq!(metrics.correct_columns(), 1);
    assert_eq!(metrics.reference_columns(), 2);
    assert!((metrics.tc_score() - 0.5).abs() < f64::EPSILON);
}

#[test]
fn reference_column_mask_validates_width_and_defines_empty_score() {
    let reference = alignment(&[("a", "AC"), ("b", "AT")]);
    let candidate = reference.clone();
    assert_eq!(
        compare_alignments_on_reference_columns(&candidate, &reference, &[true]),
        Err(MsaError::ReferenceColumnMaskLengthMismatch {
            expected: 2,
            actual: 1,
        })
    );
    let metrics =
        compare_alignments_on_reference_columns(&candidate, &reference, &[false, false]).unwrap();
    assert!((metrics.sp_score() - 1.0).abs() < f64::EPSILON);
    assert!((metrics.tc_score() - 1.0).abs() < f64::EPSILON);
}
