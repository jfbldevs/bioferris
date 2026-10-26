//! Integration checks for the public ancestral L-INS-i API.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{LinsiTermination, MafftAlphabet, MafftLinsi};
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
fn local_pairs_tree_library_and_refinement_evidence_are_complete() {
    let records = [
        dna("a", "TTACGTACGTAA"),
        dna("b", "GGACGTCGTCC"),
        dna("c", "CCACGTACGGG"),
        dna("d", "AAACGTTACGTTT"),
    ];
    let result = MafftLinsi::<Dna>::new().align(&records).unwrap();
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
        assert_eq!(pair.fft_anchor_count(), 0);
        assert_eq!(pair.generalized_skip_opening_count(), 0);
        assert_eq!(pair.generalized_skipped_pair_count(), 0);
        assert_eq!(pair.generalized_unpaired_residue_count(), 0);
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
        if pair.alignment_score() == 0.0 {
            assert_eq!(first_range, 0..0);
            assert_eq!(second_range, 0..0);
            assert_eq!(pair.alignment_column_count(), 0);
        } else {
            assert!(!first_range.is_empty());
            assert!(!second_range.is_empty());
            assert!(pair.alignment_column_count() >= first_range.len().max(second_range.len()));
            assert!(pair.alignment_column_count() <= first_range.len() + second_range.len());
        }
    }
    assert_eq!(
        result
            .pair_evidence()
            .iter()
            .map(align_mafft::LinsiPairEvidence::library_entry_count)
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
    if result.termination() == LinsiTermination::Converged {
        assert!(!result.cycles().last().unwrap().accepted());
    }
    assert!(result.final_score().is_finite());
}

#[test]
fn unrelated_terminal_flanks_are_excluded_from_pair_evidence_but_not_the_msa() {
    let records = [
        protein("a", "WWWWMKTAYIWWWW"),
        protein("b", "DDDDMKTAYIDDDD"),
        protein("c", "CCCCMKTAYICCCC"),
    ];
    let result = MafftLinsi::<Protein>::new().align(&records).unwrap();
    assert_valid(&records, result.initial_alignment());
    assert_valid(&records, result.alignment());
    for pair in result.pair_evidence() {
        assert_eq!(pair.first_aligned_range(), 4..10);
        assert_eq!(pair.second_aligned_range(), 4..10);
        assert_eq!(pair.alignment_column_count(), 6);
        assert_eq!(pair.segment_count(), 1);
        assert_eq!(pair.library_entry_count(), 6);
    }
    assert_eq!(result.library_statistics().entry_count(), 18);
    assert_eq!(result.library_statistics().segment_count(), 3);
}

#[test]
fn empty_no_hit_identical_repetitive_and_gap_cases_are_valid() {
    for records in [
        vec![dna("empty-a", ""), dna("empty-b", "")],
        vec![dna("empty", ""), dna("full", "ACGT")],
        vec![dna("a", "AAAA"), dna("b", "CCCC")],
        vec![
            dna("a", "TTACGTAA"),
            dna("b", "GGACGTCC"),
            dna("c", "CCACGTGG"),
        ],
        vec![
            dna("a", "ACACACACACAC"),
            dna("b", "CACACACACA"),
            dna("c", "ACACACAC"),
        ],
        vec![
            dna("a", "TTACGTTACGTAA"),
            dna("b", "GGACGTACGTCC"),
            dna("c", "CCACGTTACGGG"),
        ],
    ] {
        let result = MafftLinsi::<Dna>::new().align(&records).unwrap();
        assert_valid(&records, result.initial_alignment());
        assert_valid(&records, result.alignment());
        assert!(result.final_score().is_finite());
        assert!(
            result
                .pair_evidence()
                .iter()
                .all(|pair| pair.alignment_score().is_finite())
        );
    }
}

#[test]
fn all_outputs_and_audit_values_are_deterministic() {
    let records = [
        dna("a", "TTGCCCTTTGTC"),
        dna("b", "GGGCCCTTGG"),
        dna("c", "AAGCCCTTCC"),
        dna("d", "CCGCCCTTAA"),
        dna("e", "TTGCCCTTGG"),
    ];
    let first = MafftLinsi::<Dna>::new().align(&records).unwrap();
    let second = MafftLinsi::<Dna>::new().align(&records).unwrap();
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
fn parallel_refinement_is_bit_exact_and_reduces_in_partition_order() {
    let records = [
        dna("s0", "TTGCCCTTTGTC"),
        dna("s1", "ATCATACA"),
        dna("s2", "TAGTGATA"),
        dna("s3", "TCGCCCCGGCAT"),
        dna("s4", "ATCGGAGCTC"),
    ];
    let sequential = MafftLinsi::<Dna>::new()
        .with_thread_count(1)
        .unwrap()
        .align(&records)
        .unwrap();
    let parallel_aligner = MafftLinsi::<Dna>::new().with_thread_count(4).unwrap();
    assert_eq!(parallel_aligner.thread_count(), 4);
    let parallel = parallel_aligner.align(&records).unwrap();

    assert_eq!(sequential.alignment(), parallel.alignment());
    assert_eq!(sequential.initial_alignment(), parallel.initial_alignment());
    assert_eq!(sequential.distance_matrix(), parallel.distance_matrix());
    assert_eq!(sequential.guide_tree(), parallel.guide_tree());
    assert_eq!(sequential.pair_evidence(), parallel.pair_evidence());
    assert_eq!(
        sequential.progressive_merges(),
        parallel.progressive_merges()
    );
    assert_eq!(sequential.pair_weights(), parallel.pair_weights());
    assert_eq!(
        sequential.final_score().to_bits(),
        parallel.final_score().to_bits()
    );
    assert_eq!(sequential.cycles().len(), parallel.cycles().len());
    for (sequential, parallel) in sequential.cycles().iter().zip(parallel.cycles()) {
        assert_eq!(sequential.cycle_index(), parallel.cycle_index());
        assert_eq!(sequential.partition_count(), parallel.partition_count());
        assert_eq!(
            sequential.score_before().to_bits(),
            parallel.score_before().to_bits()
        );
        assert_eq!(
            sequential.best_candidate_score().map(f64::to_bits),
            parallel.best_candidate_score().map(f64::to_bits)
        );
        assert_eq!(
            sequential.score_after().to_bits(),
            parallel.score_after().to_bits()
        );
        assert_eq!(sequential.best_partition(), parallel.best_partition());
        assert_eq!(
            sequential.best_candidate_merge(),
            parallel.best_candidate_merge()
        );
        assert_eq!(sequential.accepted(), parallel.accepted());
    }
    assert!(
        parallel
            .memory_estimate()
            .candidate_parallel_workspace_peak_bytes()
            >= parallel.memory_estimate().candidate_alignment_peak_bytes()
    );
    assert!(matches!(
        MafftLinsi::<Dna>::new().with_thread_count(0),
        Err(align_mafft::MafftError::InvalidThreadCount { value: 0 })
    ));
}

#[test]
fn protein_dna_rna_and_ambiguities_are_consumed() {
    let proteins = [
        protein("p1", "XXARNDCQEGHILKMFPSTWYVBZXJUOXX"),
        protein("p2", "GGARNDCQEGHILKMFPSTWYVGG"),
        protein("p3", "AAMFPSTWYVBZXJAA"),
    ];
    let dna_records = [
        dna("d1", "TTACGTRYSWKMBDHVNAA"),
        dna("d2", "GGACGTNYSWKMBDHVNCC"),
        dna("d3", "CCRYSWKMBDHVNGG"),
    ];
    let rna_records = [
        SequenceRecord::<Rna>::new(
            "r1",
            RnaSequence::from_ascii("UUACGURYSWKMBDHVNAA").unwrap(),
        )
        .unwrap(),
        SequenceRecord::<Rna>::new(
            "r2",
            RnaSequence::from_ascii("GGACGUNYSWKMBDHVNCC").unwrap(),
        )
        .unwrap(),
        SequenceRecord::<Rna>::new("r3", RnaSequence::from_ascii("CCRYSWKMBDHVNGG").unwrap())
            .unwrap(),
    ];
    assert_valid(
        &proteins,
        MafftLinsi::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftLinsi::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftLinsi::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn pair_weights_and_local_memory_evidence_are_finite_and_complete() {
    let records = [
        dna("a", "TTACGTACGTAA"),
        dna("b", "GGACGTCGTCC"),
        dna("c", "CCACGTGG"),
        dna("d", "AAACGTTACGTTT"),
        dna("e", "GGACGACGCC"),
    ];
    let result = MafftLinsi::<Dna>::new().align(&records).unwrap();
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
    assert!(memory.pair_alignment_peak_bytes() > 0);
    assert!(memory.pair_score_rows_peak_bytes() > 0);
    assert!(memory.pair_traceback_peak_bytes() > 0);
    assert!(memory.pair_operation_peak_bytes() > 0);
    assert!(memory.pair_profile_peak_bytes() > 0);
    assert!(memory.primary_library_bytes() > 0);
    assert!(memory.estimated_peak_bytes() >= memory.returned_final_alignment_bytes());
}
