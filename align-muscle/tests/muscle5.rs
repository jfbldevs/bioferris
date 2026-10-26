//! Hand-checked public API, edge-case, and invariant tests for MUSCLE 5.

use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_msa::SequenceRecord;
use align_muscle::{GuideTreePermutation, Muscle5, Muscle5Mode, Muscle5Parameters};

fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
}

fn fast_align() -> Muscle5<Dna> {
    Muscle5::with_parameters(Muscle5Parameters::published_align().with_iterations(1, 8)).unwrap()
}

fn fast_super5() -> Muscle5<Dna> {
    let parameters = Muscle5Parameters::published_super5()
        .with_iterations(1, 4)
        .with_super5_limits(0.0, 0.30, 0.10, 3, 16, 3)
        .unwrap();
    Muscle5::with_parameters(parameters).unwrap()
}

#[test]
fn empty_singleton_is_a_valid_zero_width_alignment_in_both_modes() {
    let records = [dna("empty", "")];
    for aligner in [
        Muscle5::<Dna>::published_align(),
        Muscle5::published_super5(),
    ] {
        let result = aligner.align(&records).unwrap();
        assert_eq!(result.alignment().sequence_count(), 1);
        assert_eq!(result.alignment().width(), 0);
        assert_eq!(result.alignment().row(0).unwrap().to_ungapped_ascii(), "");
        assert_eq!(result.evaluated_objective_pair_count(), 0);
        assert_eq!(result.possible_objective_pair_count(), 0);
    }
}

#[test]
fn terminal_and_internal_gaps_preserve_pairwise_residue_order() {
    let terminal = [dna("long", "AC"), dna("short", "C")];
    let terminal_result = fast_align().align(&terminal).unwrap();
    assert_eq!(
        terminal_result
            .alignment()
            .row(0)
            .unwrap()
            .to_ascii_with_gaps(),
        "AC"
    );
    assert_eq!(
        terminal_result
            .alignment()
            .row(1)
            .unwrap()
            .to_ascii_with_gaps(),
        "-C"
    );

    let internal = [dna("long", "ACGT"), dna("short", "ACT")];
    let internal_result = fast_align().align(&internal).unwrap();
    assert_eq!(
        internal_result
            .alignment()
            .row(0)
            .unwrap()
            .to_ascii_with_gaps(),
        "ACGT"
    );
    assert_eq!(
        internal_result
            .alignment()
            .row(1)
            .unwrap()
            .to_ascii_with_gaps(),
        "AC-T"
    );
}

#[test]
fn identical_disjoint_and_repetitive_inputs_remain_valid() {
    let identical = [dna("a", "ACGTAC"), dna("b", "ACGTAC"), dna("c", "ACGTAC")];
    let result = fast_align().align(&identical).unwrap();
    assert!(
        result
            .alignment()
            .rows()
            .iter()
            .all(|row| row.to_ascii_with_gaps() == "ACGTAC")
    );

    for records in [
        vec![dna("a", "AAAAAA"), dna("c", "CCCCCC"), dna("g", "GGGGGG")],
        vec![
            dna("r0", "ACACACACAC"),
            dna("r1", "CACACACACA"),
            dna("r2", "ACACGCACAC"),
        ],
    ] {
        let result = fast_align().align(&records).unwrap();
        for (row, record) in result.alignment().rows().iter().zip(&records) {
            assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
        }
        for column in 0..result.alignment().width() {
            assert!(
                result
                    .alignment()
                    .column(column)
                    .unwrap()
                    .any(|cell| cell.is_some())
            );
        }
    }
}

#[test]
fn ambiguity_symbols_are_preserved_for_every_supported_alphabet() {
    let dna_records = [
        SequenceRecord::new("dna1", DnaSequence::from_ascii("ACGTRYSWKMBDHVN").unwrap()).unwrap(),
        SequenceRecord::new("dna2", DnaSequence::from_ascii("ACGTNNNN").unwrap()).unwrap(),
    ];
    let rna_records = [
        SequenceRecord::new("rna1", RnaSequence::from_ascii("ACGURYSWKMBDHVN").unwrap()).unwrap(),
        SequenceRecord::new("rna2", RnaSequence::from_ascii("ACGUNNNN").unwrap()).unwrap(),
    ];
    let protein_records = [
        SequenceRecord::new(
            "protein1",
            ProteinSequence::from_ascii("ACDEFGHIKLMNPQRSTVWYBZXJUO").unwrap(),
        )
        .unwrap(),
        SequenceRecord::new("protein2", ProteinSequence::from_ascii("ACDXJO").unwrap()).unwrap(),
    ];
    let dna =
        Muscle5::<Dna>::with_parameters(Muscle5Parameters::published_align().with_iterations(0, 0))
            .unwrap()
            .align(&dna_records)
            .unwrap();
    let rna =
        Muscle5::<Rna>::with_parameters(Muscle5Parameters::published_align().with_iterations(0, 0))
            .unwrap()
            .align(&rna_records)
            .unwrap();
    let protein = Muscle5::<Protein>::with_parameters(
        Muscle5Parameters::published_align().with_iterations(0, 0),
    )
    .unwrap()
    .align(&protein_records)
    .unwrap();
    for (row, record) in dna.alignment().rows().iter().zip(&dna_records) {
        assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
    }
    for (row, record) in rna.alignment().rows().iter().zip(&rna_records) {
        assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
    }
    for (row, record) in protein.alignment().rows().iter().zip(&protein_records) {
        assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
    }
}

#[test]
fn posterior_distances_checksums_and_refinement_are_auditable() {
    let records = [
        dna("a", "ACGTAC"),
        dna("b", "ACGAC"),
        dna("c", "TCGTAC"),
        dna("d", "ACGTTC"),
    ];
    let result = fast_align().align(&records).unwrap();
    assert_ne!(result.raw_posterior_checksum(), 0);
    assert_ne!(result.consistent_posterior_checksum(), 0);
    assert!(result.raw_posterior_entry_count() > 0);
    assert!(result.consistent_posterior_entry_count() > 0);
    assert_eq!(result.expected_error_distances().size(), records.len());
    for first in 0..records.len() {
        for second in 0..records.len() {
            let forward = result
                .expected_error_distances()
                .get(first, second)
                .unwrap();
            let reverse = result
                .expected_error_distances()
                .get(second, first)
                .unwrap();
            assert!((0.0..=1.0).contains(&forward));
            assert!((forward - reverse).abs() <= 1.0e-12);
        }
    }
    assert_eq!(result.progressive_steps().len(), records.len() - 1);
    assert_eq!(result.refinement_attempts().len(), 8);
    for attempt in result.refinement_attempts() {
        assert!(!attempt.first_partition().is_empty());
        assert!(attempt.first_partition().len() < records.len());
        assert!(attempt.score_after() + 1.0e-9 >= attempt.score_before());
        assert_eq!(
            attempt.accepted(),
            attempt.score_after() > attempt.score_before()
        );
    }
    assert_eq!(result.evaluated_objective_pair_count(), 6);
    assert_eq!(result.possible_objective_pair_count(), 6);
}

#[test]
fn perturbations_and_tree_permutations_are_explicit_and_deterministic() {
    let records = [
        dna("a", "ACGTAC"),
        dna("b", "ACGAC"),
        dna("c", "TCGTAC"),
        dna("d", "ACGTTC"),
        dna("e", "CGTAC"),
        dna("f", "ACCTAC"),
    ];
    let unperturbed =
        Muscle5::<Dna>::with_parameters(Muscle5Parameters::published_align().with_iterations(0, 0))
            .unwrap()
            .align(&records)
            .unwrap();
    let perturbed_parameters = Muscle5Parameters::published_align()
        .with_iterations(0, 0)
        .with_perturbation_seed(7)
        .with_guide_tree_permutation(GuideTreePermutation::Bca);
    let aligner = Muscle5::<Dna>::with_parameters(perturbed_parameters).unwrap();
    let first = aligner.align(&records).unwrap();
    let second = aligner.align(&records).unwrap();
    assert_ne!(
        unperturbed.raw_posterior_checksum(),
        first.raw_posterior_checksum()
    );
    assert_eq!(first.alignment(), second.alignment());
    assert_eq!(first.guide_tree_permutation(), GuideTreePermutation::Bca);
    assert_eq!(first.perturbation_seed(), 7);
}

#[test]
fn super5_exposes_characteristic_pipeline_evidence() {
    let records = [
        dna("a", "AAAAAA"),
        dna("b", "CCCCCC"),
        dna("c", "GGGGGG"),
        dna("d", "TTTTTT"),
        dna("e", "ACGTAC"),
        dna("f", "TGCATG"),
    ];
    let result = fast_super5().align(&records).unwrap();
    assert_eq!(result.mode(), Muscle5Mode::Super5);
    let evidence = result.super5_evidence().unwrap();
    assert_eq!(evidence.decreasing_length_order().len(), records.len());
    assert!(!evidence.representatives().is_empty());
    assert!(
        evidence
            .clusters()
            .iter()
            .all(|cluster| cluster.representative_members().len() <= 3)
    );
    assert_eq!(result.guide_tree().leaf_count(), evidence.clusters().len());
    assert!(evidence.uclust_ee_comparison_count() > 0);
    assert!(evidence.direct_pair_posterior_count() > 0);
    assert!(result.evaluated_objective_pair_count() <= result.possible_objective_pair_count());
    assert!(result.memory_estimate().super5_auxiliary_bytes() > 0);
}

#[test]
fn ensembles_follow_stratified_and_diversified_schedules() {
    let records = [dna("a", "ACGT"), dna("b", "AGT"), dna("c", "ACT")];
    let aligner =
        Muscle5::<Dna>::with_parameters(Muscle5Parameters::published_align().with_iterations(0, 0))
            .unwrap();
    let stratified = aligner.stratified_ensemble(&records, 2).unwrap();
    assert_eq!(stratified.len(), 8);
    assert_eq!(
        stratified
            .members()
            .iter()
            .map(align_muscle::Muscle5EnsembleMember::name)
            .collect::<Vec<_>>(),
        [
            "none.0", "abc.0", "acb.0", "bca.0", "none.1", "abc.1", "acb.1", "bca.1"
        ]
    );
    let diversified = aligner.diversified_ensemble(&records, 5).unwrap();
    assert_eq!(diversified.len(), 5);
    assert_eq!(diversified.members()[4].name(), "none.4");
}

#[test]
fn duplicate_identifiers_and_invalid_ensemble_sizes_are_rejected() {
    let duplicate = [dna("same", "A"), dna("same", "C")];
    assert!(Muscle5::<Dna>::published_align().align(&duplicate).is_err());
    assert!(
        Muscle5::<Dna>::published_super5()
            .align(&duplicate)
            .is_err()
    );
    let records = [dna("a", "A")];
    assert!(
        Muscle5::<Dna>::published_align()
            .stratified_ensemble(&records, 0)
            .is_err()
    );
    assert!(
        Muscle5::<Dna>::published_align()
            .diversified_ensemble(&records, 0)
            .is_err()
    );
}

/// `SplitMix64` draws for the synthetic family below.
struct SplitMix64(u64);

impl SplitMix64 {
    const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    fn residue(&mut self) -> u8 {
        Self::RESIDUES[usize::try_from(self.next() % 20).unwrap()]
    }
}

/// Deterministic synthetic protein family: a random ancestor and `count`
/// descendants with 30% substitutions and 3% single-residue insertions and
/// deletions.
fn synthetic_family(seed: u64, count: usize, length: usize) -> Vec<SequenceRecord<Protein>> {
    let mut rng = SplitMix64(seed);
    let ancestor = (0..length).map(|_| rng.residue()).collect::<Vec<_>>();
    (0..count)
        .map(|index| {
            let mut sequence = Vec::with_capacity(length + length / 10);
            for &source in &ancestor {
                match rng.next() % 100 {
                    0..=2 => {}
                    3..=5 => {
                        sequence.push(source);
                        sequence.push(rng.residue());
                    }
                    6..=35 => sequence.push(rng.residue()),
                    _ => sequence.push(source),
                }
            }
            SequenceRecord::new(
                format!("s{index}"),
                ProteinSequence::from_ascii(&sequence).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

#[test]
fn refinement_survives_rounding_differences_between_summation_orders() {
    // Version 1 stopped on this family with "PPP refinement candidate decreased
    // its partition objective": a partition realigned into the same alignment,
    // but the current objective, summed row pair by row pair, rounded onto the
    // 1e-9 step above the DP's column-by-column sum.
    let records = synthetic_family(8, 6, 300);
    let result = Muscle5::<Protein>::published_align()
        .align(&records)
        .unwrap();
    assert_eq!(result.refinement_attempts().len(), 100);
    for attempt in result.refinement_attempts() {
        assert!(attempt.score_after() >= attempt.score_before());
        assert_eq!(
            attempt.accepted(),
            attempt.score_after() > attempt.score_before()
        );
    }
    for (row, record) in result.alignment().rows().iter().zip(&records) {
        assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
    }
}
