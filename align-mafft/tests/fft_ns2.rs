//! Integration checks for the public two-pass FFT-NS-2 API.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{MafftAlphabet, MafftFftNs1, MafftFftNs2};
use align_msa::{MultipleAlignment, SequenceRecord};

fn assert_valid<A: MafftAlphabet>(records: &[SequenceRecord<A>], alignment: &MultipleAlignment<A>) {
    assert_eq!(alignment.sequence_count(), records.len());
    assert!(
        alignment
            .rows()
            .iter()
            .all(|row| row.len() == alignment.width())
    );
    for (record, row) in records.iter().zip(alignment.rows()) {
        assert_eq!(record.identifier(), row.identifier());
        assert_eq!(record.sequence().to_ascii_string(), row.to_ungapped_ascii());
    }
    for column in 0..alignment.width() {
        assert!(
            alignment
                .column(column)
                .expect("column exists")
                .any(|cell| cell.is_some())
        );
    }
}

fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
}

fn protein(identifier: &str, sequence: &str) -> SequenceRecord<Protein> {
    SequenceRecord::new(identifier, ProteinSequence::from_ascii(sequence).unwrap()).unwrap()
}

#[test]
fn hand_calculated_internal_gap_is_recovered_after_both_passes() {
    let records = [dna("long", "ACGT"), dna("short", "AGT")];
    let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    assert_eq!(
        result.first_pass().alignment().rows()[0].to_ascii_with_gaps(),
        "ACGT"
    );
    assert_eq!(
        result.first_pass().alignment().rows()[1].to_ascii_with_gaps(),
        "A-GT"
    );
    assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
    assert_eq!(result.alignment().rows()[1].to_ascii_with_gaps(), "A-GT");
    assert_eq!(result.first_pass().merges().len(), 1);
    assert_eq!(result.second_merges().len(), 1);
    assert_valid(&records, result.alignment());
}

#[test]
fn empty_leading_internal_and_terminal_cases_are_valid() {
    let records = [
        dna("empty", ""),
        dna("leading", "CGT"),
        dna("full", "ACGTA"),
        dna("internal", "ACTA"),
        dna("terminal", "ACGT"),
    ];
    let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.first_pass().alignment());
    assert_valid(&records, result.alignment());
    assert!(
        result.alignment().rows()[0]
            .to_ascii_with_gaps()
            .bytes()
            .all(|byte| byte == b'-')
    );
}

#[test]
fn completely_different_sequences_have_finite_saturated_distances() {
    let records = [
        dna("a", "AAAAAAAA"),
        dna("c", "CCCCCCCC"),
        dna("g", "GGGGGGGG"),
    ];
    let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.alignment());
    for row in 0..records.len() {
        for column in 0..records.len() {
            let distance = result.second_distance_matrix().get(row, column).unwrap();
            assert!(distance.is_finite() && distance >= 0.0);
        }
    }
    assert!(result.second_root_merge_score().is_some_and(f64::is_finite));
}

#[test]
fn first_pass_is_exactly_the_standalone_fft_ns1_profile() {
    let records = [
        dna("a", "TACGTACGTACGT"),
        dna("b", "ACGTCGTACGT"),
        dna("c", "ACGTACGTACG"),
        dna("d", "TACGTTACGTACGT"),
    ];
    let standalone = MafftFftNs1::<Dna>::new().align(&records).unwrap();
    let nested = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    let first = nested.first_pass();
    assert_eq!(first.alignment(), standalone.alignment());
    assert_eq!(first.distance_matrix(), standalone.distance_matrix());
    assert_eq!(first.guide_tree(), standalone.guide_tree());
    assert_eq!(first.sequence_weights(), standalone.sequence_weights());
    assert_eq!(
        first.effective_tuple_length(),
        standalone.effective_tuple_length()
    );
    assert_eq!(first.merges(), standalone.merges());
    assert_eq!(first.root_merge_score(), standalone.root_merge_score());
}

#[test]
fn second_tree_corrects_a_dayhoff_tuple_collision() {
    let records = [
        protein("alanine", "AAAAAAAAAAAA"),
        protein("glycine", "GGGGGGGGGGGG"),
        protein("near_alanine", "AAAAAAAAAAAD"),
    ];
    let result = MafftFftNs2::<Protein>::new().align(&records).unwrap();
    assert_eq!(
        result.first_pass().guide_tree().postorder_leaves(),
        vec![0, 1, 2]
    );
    assert_eq!(result.second_guide_tree().postorder_leaves(), vec![0, 2, 1]);
    assert!(result.guide_tree_changed());
    assert!(
        result.second_distance_matrix().get(0, 2).unwrap()
            < result.second_distance_matrix().get(0, 1).unwrap()
    );
    assert_valid(&records, result.alignment());
}

#[test]
fn unchanged_tree_still_has_a_complete_second_pass() {
    let records = [
        dna("a", "ACGTACGTACGT"),
        dna("b", "ACGTACGTACGT"),
        dna("c", "ACGTACGTACGT"),
    ];
    let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    assert!(!result.guide_tree_changed());
    assert_eq!(result.first_pass().merges().len(), records.len() - 1);
    assert_eq!(result.second_merges().len(), records.len() - 1);
    assert_eq!(result.alignment(), result.first_pass().alignment());
}

#[test]
fn repetitive_input_is_fully_deterministic() {
    let records = [
        dna("a", "ACACACACACACACAC"),
        dna("b", "ACACACACACAC"),
        dna("c", "CACACACACACACA"),
    ];
    let aligner = MafftFftNs2::<Dna>::new();
    let first = aligner.align(&records).unwrap();
    let second = aligner.align(&records).unwrap();
    assert_eq!(first.alignment(), second.alignment());
    assert_eq!(
        first.first_pass().alignment(),
        second.first_pass().alignment()
    );
    assert_eq!(
        first.second_distance_matrix(),
        second.second_distance_matrix()
    );
    assert_eq!(first.second_guide_tree(), second.second_guide_tree());
    assert_eq!(
        first.second_sequence_weights(),
        second.second_sequence_weights()
    );
    assert_eq!(first.second_merges(), second.second_merges());
    assert_eq!(
        first.second_root_merge_score(),
        second.second_root_merge_score()
    );
}

#[test]
fn all_supported_ambiguity_symbols_are_consumed() {
    let proteins = [
        protein("p1", "ARNDCQEGHILKMFPSTWYVBZXJUO"),
        protein("p2", "ARNDCQEGHILKMFPSTWYV"),
    ];
    let dna_records = [dna("d1", "ACGTRYSWKMBDHVN"), dna("d2", "ACGTNYSWKMBDHVN")];
    let rna_records = [
        SequenceRecord::<Rna>::new("r1", RnaSequence::from_ascii("ACGURYSWKMBDHVN").unwrap())
            .unwrap(),
        SequenceRecord::<Rna>::new("r2", RnaSequence::from_ascii("ACGUNYSWKMBDHVN").unwrap())
            .unwrap(),
    ];
    assert_valid(
        &proteins,
        MafftFftNs2::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftFftNs2::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftFftNs2::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn second_weights_are_positive_finite_and_normalized() {
    let records = [
        dna("a", "ACGTACGT"),
        dna("b", "ACGTCGT"),
        dna("c", "ACGT"),
        dna("d", "TACGTACG"),
    ];
    let result = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    assert!(
        result
            .second_sequence_weights()
            .iter()
            .all(|weight| weight.is_finite() && *weight > 0.0)
    );
    assert!((result.second_sequence_weights().iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
}
