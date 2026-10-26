//! Integration checks for the public ancestral G-INS-i API.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{GinsiTermination, MafftAlphabet, MafftGinsi};
use align_msa::{MultipleAlignment, SequenceRecord};

fn assert_valid<A: MafftAlphabet>(records: &[SequenceRecord<A>], alignment: &MultipleAlignment<A>) {
    assert_eq!(alignment.sequence_count(), records.len());
    for (record, row) in records.iter().zip(alignment.rows()) {
        assert_eq!(row.identifier(), record.identifier());
        assert_eq!(row.len(), alignment.width());
        assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
    }
    for column in 0..alignment.width() {
        assert!(alignment.column(column).unwrap().any(|cell| cell.is_some()));
    }
}

fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
}

fn protein(identifier: &str, sequence: &str) -> SequenceRecord<Protein> {
    SequenceRecord::new(identifier, ProteinSequence::from_ascii(sequence).unwrap()).unwrap()
}

#[test]
fn all_pairs_tree_library_and_refinement_evidence_are_complete() {
    let records = [
        dna("a", "TACGTACGT"),
        dna("b", "ACGTCGT"),
        dna("c", "ACGTACG"),
        dna("d", "TACGTTACGT"),
    ];
    let result = MafftGinsi::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.initial_alignment());
    assert_valid(&records, result.alignment());
    assert_eq!(result.pair_evidence().len(), 6);
    assert_eq!(result.library_statistics().pair_count(), 6);
    assert_eq!(result.progressive_merges().len(), 3);
    assert_eq!(result.guide_tree().leaf_count(), 4);
    assert!(result.library_statistics().entry_count() > 0);
    assert!(result.library_statistics().segment_count() > 0);
    for pair in result.pair_evidence() {
        assert!(pair.first_sequence() < pair.second_sequence());
        assert_eq!(pair.generalized_skip_opening_count(), 0);
        assert_eq!(pair.generalized_skipped_pair_count(), 0);
        assert_eq!(pair.generalized_unpaired_residue_count(), 0);
        assert!((0.0..=1.0).contains(&pair.distance()));
        assert_eq!(
            result
                .distance_matrix()
                .get(pair.first_sequence(), pair.second_sequence()),
            Some(pair.distance())
        );
    }
    assert_eq!(
        result
            .pair_evidence()
            .iter()
            .map(align_mafft::GinsiPairEvidence::library_entry_count)
            .sum::<usize>(),
        result.library_statistics().entry_count()
    );
    for cycle in result.cycles() {
        assert_eq!(cycle.partition_count(), 2 * records.len() - 3);
        if cycle.accepted() {
            assert!(cycle.score_after() > cycle.score_before());
        } else {
            assert_eq!(
                cycle.score_after().to_bits(),
                cycle.score_before().to_bits()
            );
        }
    }
    if result.termination() == GinsiTermination::Converged {
        assert!(!result.cycles().last().unwrap().accepted());
    }
    assert!(result.final_score().is_finite());
}

#[test]
fn empty_leading_internal_terminal_identical_and_repetitive_inputs_are_valid() {
    for records in [
        vec![
            dna("empty", ""),
            dna("leading", "CGT"),
            dna("full", "ACGTA"),
            dna("internal", "ACTA"),
            dna("terminal", "ACGT"),
        ],
        vec![
            dna("a", "ACGTACGT"),
            dna("b", "ACGTACGT"),
            dna("c", "ACGTACGT"),
        ],
        vec![
            dna("a", "ACACACACACAC"),
            dna("b", "CACACACACA"),
            dna("c", "ACACACAC"),
        ],
    ] {
        let result = MafftGinsi::<Dna>::new().align(&records).unwrap();
        assert_valid(&records, result.alignment());
        assert!(result.final_score().is_finite());
    }
}

#[test]
fn results_and_all_audit_values_are_deterministic() {
    let records = [
        dna("a", "TTGCCCTTTGTC"),
        dna("b", "ATCATACA"),
        dna("c", "TAGTGATA"),
        dna("d", "TCGCCCCGGCAT"),
        dna("e", "ATCGGAGCTC"),
    ];
    let first = MafftGinsi::<Dna>::new().align(&records).unwrap();
    let second = MafftGinsi::<Dna>::new().align(&records).unwrap();
    assert_eq!(first.alignment(), second.alignment());
    assert_eq!(first.initial_alignment(), second.initial_alignment());
    assert_eq!(first.distance_matrix(), second.distance_matrix());
    assert_eq!(first.guide_tree(), second.guide_tree());
    assert_eq!(
        first.progressive_sequence_weights(),
        second.progressive_sequence_weights()
    );
    assert_eq!(first.pair_evidence(), second.pair_evidence());
    assert_eq!(first.progressive_merges(), second.progressive_merges());
    assert_eq!(first.pair_weights(), second.pair_weights());
    assert_eq!(
        first.final_score().to_bits(),
        second.final_score().to_bits()
    );
    assert_eq!(first.memory_estimate(), second.memory_estimate());
    assert_eq!(first.cycles().len(), second.cycles().len());
    for (left, right) in first.cycles().iter().zip(second.cycles()) {
        assert_eq!(
            left.score_before().to_bits(),
            right.score_before().to_bits()
        );
        assert_eq!(
            left.best_candidate_score().map(f64::to_bits),
            right.best_candidate_score().map(f64::to_bits)
        );
        assert_eq!(left.score_after().to_bits(), right.score_after().to_bits());
        assert_eq!(left.best_partition(), right.best_partition());
        assert_eq!(left.best_candidate_merge(), right.best_candidate_merge());
    }
}

#[test]
fn all_supported_alphabets_and_ambiguities_are_consumed() {
    let proteins = [
        protein("p1", "ARNDCQEGHILKMFPSTWYVBZXJUO"),
        protein("p2", "ARNDCQEGHILKMFPSTWYV"),
        protein("p3", "MFPSTWYVBZXJ"),
    ];
    let dna_records = [
        dna("d1", "ACGTRYSWKMBDHVN"),
        dna("d2", "ACGTNYSWKMBDHVN"),
        dna("d3", "RYSWKMBDHVN"),
    ];
    let rna_records = [
        SequenceRecord::<Rna>::new("r1", RnaSequence::from_ascii("ACGURYSWKMBDHVN").unwrap())
            .unwrap(),
        SequenceRecord::<Rna>::new("r2", RnaSequence::from_ascii("ACGUNYSWKMBDHVN").unwrap())
            .unwrap(),
        SequenceRecord::<Rna>::new("r3", RnaSequence::from_ascii("RYSWKMBDHVN").unwrap()).unwrap(),
    ];
    assert_valid(
        &proteins,
        MafftGinsi::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftGinsi::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftGinsi::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn gotoh_pair_weights_are_positive_symmetric_and_unit_mean() {
    let records = [
        dna("a", "ACGTACGT"),
        dna("b", "ACGTCGT"),
        dna("c", "ACGT"),
        dna("d", "TACGTACG"),
        dna("e", "TACGACG"),
    ];
    let result = MafftGinsi::<Dna>::new().align(&records).unwrap();
    let weights = result.pair_weights();
    for first in 0..records.len() {
        assert_eq!(weights.get(first, first), Some(0.0));
        for second in 0..first {
            let value = weights.get(first, second).unwrap();
            assert!(value.is_finite() && value > 0.0);
            assert_eq!(
                value.to_bits(),
                weights.get(second, first).unwrap().to_bits()
            );
        }
    }
    assert!((weights.sum_of_pairs() - 10.0).abs() < 1.0e-12);
    let memory = result.memory_estimate();
    assert!(memory.estimated_peak_bytes() >= memory.returned_final_alignment_bytes());
    assert!(memory.primary_library_bytes() > 0);
}
