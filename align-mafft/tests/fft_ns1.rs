//! Integration checks for the public FFT-NS-1 API and important input classes.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{MafftAlphabet, MafftFftNs1};
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

#[test]
fn hand_calculated_internal_gap_is_recovered() {
    let records = [dna("long", "ACGT"), dna("short", "AGT")];
    let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
    assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "ACGT");
    assert_eq!(result.alignment().rows()[1].to_ascii_with_gaps(), "A-GT");
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
    let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.alignment());
    assert!(
        result.alignment().rows()[0]
            .to_ascii_with_gaps()
            .bytes()
            .all(|byte| byte == b'-')
    );
}

#[test]
fn completely_different_sequences_remain_consumption_valid() {
    let records = [
        dna("a", "AAAAAAAA"),
        dna("c", "CCCCCCCC"),
        dna("g", "GGGGGGGG"),
    ];
    let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.alignment());
    assert!(result.root_merge_score().is_some_and(f64::is_finite));
}

#[test]
fn highly_repetitive_input_is_deterministic() {
    let records = [
        dna("a", "ACACACACACACACAC"),
        dna("b", "ACACACACACAC"),
        dna("c", "CACACACACACACA"),
    ];
    let aligner = MafftFftNs1::<Dna>::new();
    let first = aligner.align(&records).unwrap();
    let second = aligner.align(&records).unwrap();
    assert_eq!(
        first
            .alignment()
            .rows()
            .iter()
            .map(align_msa::AlignedRow::to_ascii_with_gaps)
            .collect::<Vec<_>>(),
        second
            .alignment()
            .rows()
            .iter()
            .map(align_msa::AlignedRow::to_ascii_with_gaps)
            .collect::<Vec<_>>()
    );
}

#[test]
fn all_supported_ambiguity_symbols_are_consumed() {
    let proteins = [
        SequenceRecord::<Protein>::new(
            "p1",
            ProteinSequence::from_ascii("ARNDCQEGHILKMFPSTWYVBZXJUO").unwrap(),
        )
        .unwrap(),
        SequenceRecord::<Protein>::new(
            "p2",
            ProteinSequence::from_ascii("ARNDCQEGHILKMFPSTWYV").unwrap(),
        )
        .unwrap(),
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
        MafftFftNs1::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftFftNs1::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftFftNs1::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn memory_accounting_is_monotone_and_nonzero_for_nontrivial_input() {
    let records = [dna("a", "ACGTACGT"), dna("b", "ACGTCGT"), dna("c", "ACGT")];
    let result = MafftFftNs1::<Dna>::new().align(&records).unwrap();
    let memory = result.memory_estimate();
    assert!(memory.tuple_storage_bytes() > 0);
    assert!(memory.distance_matrix_bytes() >= 9 * size_of::<f64>());
    assert!(memory.returned_alignment_bytes() > 0);
    assert!(memory.estimated_peak_bytes() >= memory.returned_alignment_bytes());
}

#[test]
fn input_symbols_remain_typed() {
    fn alphabet_name<A: Alphabet>() -> &'static str {
        A::NAME
    }
    assert_eq!(alphabet_name::<Dna>(), "DNA");
}
