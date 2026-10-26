//! Integration checks for the public ancestral E-INS-i API.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{EinsiTermination, MafftAlphabet, MafftEinsi};
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
fn generalized_pairs_tree_library_and_refinement_evidence_are_complete() {
    let records = [
        dna("a", "ACGTACGTACGTAAAAAAAAAAAAACGTACGTACGT"),
        dna("b", "ACGTACGTACGTCCCCCCCCCCCCCACGTACGTACGT"),
        dna("c", "ACGTACGTACGTGGGGGGGGGGGACGTACGTACGT"),
        dna("d", "ACGTACGTACGTTTTTTTTTTTTTTTACGTACGTACGT"),
    ];
    let result = MafftEinsi::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.initial_alignment());
    assert_valid(&records, result.alignment());
    assert_eq!(result.pair_evidence().len(), 6);
    assert_eq!(result.library_statistics().pair_count(), 6);
    assert_eq!(result.progressive_merges().len(), 3);
    assert_eq!(result.guide_tree().leaf_count(), 4);
    assert!(result.library_statistics().entry_count() > 0);
    assert!(result.library_statistics().segment_count() >= 12);
    let mut generalized_pair_count = 0usize;
    for pair in result.pair_evidence() {
        assert!(pair.first_sequence() < pair.second_sequence());
        assert_eq!(pair.fft_anchor_count(), 0);
        assert!((0.0..=1.0).contains(&pair.distance()));
        assert!(pair.alignment_score().is_finite() && pair.alignment_score() >= 0.0);
        assert_eq!(
            result
                .distance_matrix()
                .get(pair.first_sequence(), pair.second_sequence()),
            Some(pair.distance())
        );
        let first_range = pair.first_aligned_range();
        let second_range = pair.second_aligned_range();
        assert!(first_range.end <= records[pair.first_sequence()].sequence().len());
        assert!(second_range.end <= records[pair.second_sequence()].sequence().len());
        if pair.generalized_skip_opening_count() > 0 {
            generalized_pair_count += 1;
            assert!(
                pair.generalized_skipped_pair_count() >= 1,
                "generalized path did not consume a paired skip: {pair:?}"
            );
        }
        assert!(pair.segment_count() >= 2);
        assert!(
            pair.library_entry_count() + pair.generalized_skipped_pair_count()
                <= first_range.len().min(second_range.len())
        );
        assert!(pair.alignment_column_count() >= first_range.len().max(second_range.len()));
        assert!(pair.alignment_column_count() <= first_range.len() + second_range.len());
    }
    assert!(
        generalized_pair_count >= 1,
        "the multiple-motif case must exercise a generalized skip"
    );
    assert_eq!(
        result
            .pair_evidence()
            .iter()
            .map(align_mafft::EinsiPairEvidence::library_entry_count)
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
    if result.termination() == EinsiTermination::Converged {
        assert!(!result.cycles().last().unwrap().accepted());
    }
    assert!(result.final_score().is_finite());
}

#[test]
fn empty_no_hit_identical_repetitive_and_indel_cases_are_valid() {
    for records in [
        vec![dna("empty-a", ""), dna("empty-b", "")],
        vec![dna("empty", ""), dna("full", "ACGT")],
        vec![dna("a", "AAAA"), dna("b", "CCCC")],
        vec![dna("a", "ACGTACGT"), dna("b", "ACGTACGT")],
        vec![
            dna("a", "ACACACACAAAAAAAAACACACAC"),
            dna("b", "ACACACACCCCCCCCCCACACACAC"),
            dna("c", "ACACACACGGGGGGGGGACACACAC"),
        ],
        vec![
            dna("a", "TTACGTACGTAAAAAAAAAACGTACGTAA"),
            dna("b", "GGACGTCGTCCCCCCCCCCCACGTCGTCC"),
            dna("c", "CCACGTACGGGGGGGGGGACGTACGGG"),
        ],
    ] {
        let result = MafftEinsi::<Dna>::new().align(&records).unwrap();
        assert_valid(&records, result.initial_alignment());
        assert_valid(&records, result.alignment());
        assert!(result.final_score().is_finite());
        assert!(
            result
                .pair_evidence()
                .iter()
                .all(|pair| pair.alignment_score().is_finite())
        );
        for pair in result.pair_evidence() {
            if pair.alignment_score() == 0.0 {
                assert_eq!(pair.first_aligned_range(), 0..0);
                assert_eq!(pair.second_aligned_range(), 0..0);
                assert_eq!(pair.alignment_column_count(), 0);
                assert_eq!(pair.generalized_skip_opening_count(), 0);
                assert_eq!(pair.generalized_skipped_pair_count(), 0);
                assert_eq!(pair.generalized_unpaired_residue_count(), 0);
            }
        }
    }
}

#[test]
fn protein_dna_rna_and_ambiguities_are_consumed() {
    let proteins = [
        protein(
            "p1",
            "ARNDCQEGHILKMFPSTWYVAAAAAAAAAAAARNDCQEGHILKMFPSTWYVBZXJUO",
        ),
        protein(
            "p2",
            "ARNDCQEGHILKMFPSTWYVCCCCCCCCCCCCARNDCQEGHILKMFPSTWYVBZXJUO",
        ),
        protein(
            "p3",
            "ARNDCQEGHILKMFPSTWYVGGGGGGGGGGARNDCQEGHILKMFPSTWYVBZXJUO",
        ),
    ];
    let dna_records = [
        dna("d1", "ACGTRYSWKMBDHVNAAAAAAAAAACGTRYSWKMBDHVN"),
        dna("d2", "ACGTRYSWKMBDHVNCCCCCCCCCCACGTRYSWKMBDHVN"),
        dna("d3", "ACGTRYSWKMBDHVNGGGGGGGGGACGTRYSWKMBDHVN"),
    ];
    let rna_records = [
        SequenceRecord::<Rna>::new(
            "r1",
            RnaSequence::from_ascii("ACGURYSWKMBDHVNAAAAAAAAAACGURYSWKMBDHVN").unwrap(),
        )
        .unwrap(),
        SequenceRecord::<Rna>::new(
            "r2",
            RnaSequence::from_ascii("ACGURYSWKMBDHVNCCCCCCCCCCACGURYSWKMBDHVN").unwrap(),
        )
        .unwrap(),
        SequenceRecord::<Rna>::new(
            "r3",
            RnaSequence::from_ascii("ACGURYSWKMBDHVNGGGGGGGGGACGURYSWKMBDHVN").unwrap(),
        )
        .unwrap(),
    ];
    assert_valid(
        &proteins,
        MafftEinsi::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftEinsi::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftEinsi::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn all_outputs_audit_values_and_memory_are_deterministic() {
    let records = [
        dna("a", "ACGTACGTAAAAAAAAACGTACGT"),
        dna("b", "ACGTACGTCCCCCCCCCACGTACGT"),
        dna("c", "ACGTACGTGGGGGGGGGACGTACGT"),
        dna("d", "ACGTACGTTTTTTTTTTTACGTACGT"),
    ];
    let first = MafftEinsi::<Dna>::new().align(&records).unwrap();
    let second = MafftEinsi::<Dna>::new().align(&records).unwrap();
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

    let memory = first.memory_estimate();
    assert!(memory.pair_alignment_peak_bytes() > 0);
    assert!(memory.pair_score_rows_peak_bytes() > 0);
    assert!(memory.pair_traceback_peak_bytes() > 0);
    assert!(memory.pair_operation_peak_bytes() > 0);
    assert!(memory.pair_profile_peak_bytes() > 0);
    assert!(memory.primary_library_bytes() > 0);
    assert!(memory.estimated_peak_bytes() >= memory.returned_final_alignment_bytes());
}
