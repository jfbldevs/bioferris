//! Integration checks for the public iterative FFT-NS-i API.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence, Symbol};
use align_mafft::{FftNsiPairWeights, FftNsiTermination, MafftAlphabet, MafftFftNs2, MafftFftNsi};
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

type PairSignature<A> = Vec<(Option<Symbol<A>>, Option<Symbol<A>>)>;

fn pair_signature<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    first: usize,
    second: usize,
) -> PairSignature<A> {
    alignment.rows()[first]
        .symbols()
        .iter()
        .copied()
        .zip(alignment.rows()[second].symbols().iter().copied())
        .filter(|(left, right)| left.is_some() || right.is_some())
        .collect()
}

fn canonical_dna_index(symbol: Symbol<Dna>) -> usize {
    match symbol.as_ascii() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        byte => panic!("oracle requires a canonical DNA symbol, got {byte}"),
    }
}

fn independent_k80_matrix() -> [[f64; 4]; 4] {
    let exp_four_beta = (-2.0_f64).exp();
    let exp_alpha_beta = (-3.0_f64).exp();
    let same = 0.25 + 0.25 * exp_four_beta + 0.5 * exp_alpha_beta;
    let transition = 0.25 + 0.25 * exp_four_beta - 0.5 * exp_alpha_beta;
    let transversion = 0.25 - 0.25 * exp_four_beta;
    let mut raw = [[0.0; 4]; 4];
    for (first, row) in raw.iter_mut().enumerate() {
        for (second, value) in row.iter_mut().enumerate() {
            let probability = if first == second {
                same
            } else if matches!((first, second), (0, 2) | (2, 0) | (1, 3) | (3, 1)) {
                transition
            } else {
                transversion
            };
            *value = (probability / 0.25).ln();
        }
    }
    let average_identical = (0..4).map(|index| 0.25 * raw[index][index]).sum::<f64>();
    let average_independent = raw.iter().flatten().map(|score| score / 16.0).sum::<f64>();
    let scale = average_identical - average_independent;
    let mut normalized = [[0.0; 4]; 4];
    for (first, row) in normalized.iter_mut().enumerate() {
        for (second, value) in row.iter_mut().enumerate() {
            *value = (raw[first][second] - average_independent) / scale + 0.06;
        }
    }
    normalized
}

fn independent_dna_wsp(alignment: &MultipleAlignment<Dna>, weights: &FftNsiPairWeights) -> f64 {
    let matrix = independent_k80_matrix();
    let mut total = 0.0;
    for first in 0..alignment.sequence_count() {
        for second in 0..first {
            let mut pair_score = 0.0;
            let mut gap_orientation = 0i8;
            for (left, right) in alignment.rows()[first]
                .symbols()
                .iter()
                .copied()
                .zip(alignment.rows()[second].symbols().iter().copied())
            {
                match (left, right) {
                    (Some(left_symbol), Some(right_symbol)) => {
                        pair_score += matrix[canonical_dna_index(left_symbol)]
                            [canonical_dna_index(right_symbol)];
                        gap_orientation = 0;
                    }
                    (Some(_), None) => {
                        if gap_orientation != 1 {
                            pair_score -= 2.4;
                            gap_orientation = 1;
                        }
                    }
                    (None, Some(_)) => {
                        if gap_orientation != -1 {
                            pair_score -= 2.4;
                            gap_orientation = -1;
                        }
                    }
                    (None, None) => {}
                }
            }
            total += weights.get(first, second).unwrap() * pair_score;
        }
    }
    total
}

#[test]
fn fixed_best_first_case_improves_and_preserves_both_subalignments() {
    let records = [
        dna("s0", "TTGCCCTTTGTC"),
        dna("s1", "ATCATACA"),
        dna("s2", "TAGTGATA"),
        dna("s3", "TCGCCCCGGCAT"),
        dna("s4", "ATCGGAGCTC"),
    ];
    let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
    assert_eq!(result.accepted_cycle_count(), 1);
    assert_eq!(result.cycles().len(), 2);
    assert_eq!(result.termination(), FftNsiTermination::Converged);
    let accepted = &result.cycles()[0];
    assert!(accepted.accepted());
    assert!(accepted.score_after() > accepted.score_before());
    assert_eq!(accepted.partition_count(), 2 * records.len() - 3);
    let partition = accepted.accepted_partition().unwrap();
    for members in [partition.left_members(), partition.right_members()] {
        for first_slot in 0..members.len() {
            for second_slot in 0..first_slot {
                assert_eq!(
                    pair_signature(
                        result.initial().alignment(),
                        members[first_slot],
                        members[second_slot]
                    ),
                    pair_signature(
                        result.alignment(),
                        members[first_slot],
                        members[second_slot]
                    )
                );
            }
        }
    }
    assert!(result.changed_initial_alignment());
    assert_valid(&records, result.alignment());
}

#[test]
fn final_wsp_matches_an_independent_published_equation_oracle() {
    let records = [
        dna("a", "ACGTACGT"),
        dna("b", "ACGTCGT"),
        dna("c", "ACTTACGT"),
        dna("d", "ACGTAC"),
    ];
    let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
    let independently_replayed =
        independent_dna_wsp(result.alignment(), result.final_pair_weights());
    let tolerance = 1.0e-12 * (1.0 + result.final_score().abs());
    assert!((result.final_score() - independently_replayed).abs() <= tolerance);
}

#[test]
fn nested_initial_stage_is_exactly_standalone_fft_ns2() {
    let records = [
        dna("a", "TACGTACGTACGT"),
        dna("b", "ACGTCGTACGT"),
        dna("c", "ACGTACGTACG"),
        dna("d", "TACGTTACGTACGT"),
    ];
    let standalone = MafftFftNs2::<Dna>::new().align(&records).unwrap();
    let iterative = MafftFftNsi::<Dna>::new().align(&records).unwrap();
    let initial = iterative.initial();
    assert_eq!(initial.alignment(), standalone.alignment());
    assert_eq!(
        initial.first_pass().alignment(),
        standalone.first_pass().alignment()
    );
    assert_eq!(
        initial.second_distance_matrix(),
        standalone.second_distance_matrix()
    );
    assert_eq!(initial.second_guide_tree(), standalone.second_guide_tree());
    assert_eq!(
        initial.second_sequence_weights(),
        standalone.second_sequence_weights()
    );
    assert_eq!(initial.second_merges(), standalone.second_merges());
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
    let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
    assert_valid(&records, result.alignment());
    assert!(
        result.alignment().rows()[0]
            .to_ascii_with_gaps()
            .bytes()
            .all(|byte| byte == b'-')
    );
    assert!(result.final_score().is_finite());
}

#[test]
fn identical_divergent_and_repetitive_inputs_are_deterministic() {
    for records in [
        vec![
            dna("a", "ACGTACGT"),
            dna("b", "ACGTACGT"),
            dna("c", "ACGTACGT"),
        ],
        vec![
            dna("a", "AAAAAAAA"),
            dna("c", "CCCCCCCC"),
            dna("g", "GGGGGGGG"),
        ],
        vec![
            dna("a", "ACACACACACAC"),
            dna("b", "CACACACACA"),
            dna("c", "ACACACAC"),
        ],
    ] {
        let aligner = MafftFftNsi::<Dna>::new();
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_valid(&records, first.alignment());
        assert_eq!(first.alignment(), second.alignment());
        assert_eq!(first.final_pair_weights(), second.final_pair_weights());
        assert_eq!(
            first.final_score().to_bits(),
            second.final_score().to_bits()
        );
        assert_eq!(first.cycles().len(), second.cycles().len());
    }
}

#[test]
fn all_supported_alphabets_and_ambiguity_symbols_are_consumed() {
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
        MafftFftNsi::<Protein>::new()
            .align(&proteins)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &dna_records,
        MafftFftNsi::<Dna>::new()
            .align(&dna_records)
            .unwrap()
            .alignment(),
    );
    assert_valid(
        &rna_records,
        MafftFftNsi::<Rna>::new()
            .align(&rna_records)
            .unwrap()
            .alignment(),
    );
}

#[test]
fn final_pair_weights_are_positive_symmetric_and_normalized() {
    let records = [
        dna("a", "ACGTACGT"),
        dna("b", "ACGTCGT"),
        dna("c", "ACGT"),
        dna("d", "TACGTACG"),
        dna("e", "TACGACG"),
    ];
    let result = MafftFftNsi::<Dna>::new().align(&records).unwrap();
    let weights = result.final_pair_weights();
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
    #[allow(clippy::cast_precision_loss)]
    let expected = (records.len() * (records.len() - 1) / 2) as f64;
    assert!((weights.sum_of_pairs() - expected).abs() < 1.0e-12);
}
