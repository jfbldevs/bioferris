//! Public API and manually derived acceptance cases for Clustal Omega.

use align_clustal::{
    ClustalOmega, ClustalOmegaError, ClustalOmegaParameters, OmegaEngine, OmegaEngineMode,
    OmegaGuideTreeMethod, OmegaGuideTreeMode,
};
use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence, Sequence};
use align_msa::{MsaError, SequenceRecord};

fn dna_records(inputs: &[(&str, &str)]) -> Vec<SequenceRecord<Dna>> {
    inputs
        .iter()
        .map(|(identifier, sequence)| {
            SequenceRecord::new(
                *identifier,
                DnaSequence::from_ascii(sequence.as_bytes()).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

#[test]
fn published_profiles_are_versioned_by_molecule_class() {
    let protein = ClustalOmega::<Protein>::published();
    let dna = ClustalOmega::<Dna>::published();
    let rna = ClustalOmega::<Rna>::published();

    assert_eq!(
        protein.parameters().molecule_profile().identifier(),
        "CLUSTALO-1.0.2-default-protein-v1"
    );
    assert_eq!(
        dna.parameters().molecule_profile().identifier(),
        "CLUSTALO-1.2.4-default-nucleotide-v1"
    );
    assert_eq!(
        rna.parameters().molecule_profile().identifier(),
        "CLUSTALO-1.2.4-default-nucleotide-v1"
    );
    assert_eq!(protein.parameters().one_tuple_gap_penalty(), 3);
    assert_eq!(dna.parameters().one_tuple_gap_penalty(), 3);
    assert_eq!(protein.parameters().thread_count(), 1);
    assert_eq!(
        protein.parameters().emission_pseudocounts(),
        (1.0, 1.5, 1.0)
    );
    assert_eq!(dna.parameters().emission_pseudocounts(), (0.55, 1.5, 1.0));
}

#[test]
fn one_tuple_distance_matches_hand_calculated_gap_event_cases() {
    let aligner = ClustalOmega::<Dna>::published();
    let sequence = |text: &str| DnaSequence::from_ascii(text.as_bytes()).unwrap();
    let assert_distance = |first: &str, second: &str, expected: f64| {
        let observed = aligner
            .pair_distance(&sequence(first), &sequence(second))
            .unwrap();
        assert!(
            (observed - expected).abs() < 1.0e-12,
            "distance {first:?}/{second:?}: {observed} versus {expected}"
        );
    };

    assert_distance("ACGT", "ACGT", 0.0);
    assert_distance("ACGT", "AGGT", 0.25);
    assert_distance("ACGTACGT", "ACGACGT", 0.428_571_429);
    assert_distance("ACGTACGT", "CGTACG", 0.0);
    assert_distance("", "", 0.0);
    assert_distance("", "A", 1.0);
}

#[test]
fn public_pipeline_accepts_protein_dna_and_rna_without_changing_sequences() {
    let protein_records = [
        SequenceRecord::new(
            "human",
            ProteinSequence::from_ascii(b"MKTAYIAKQRQISFVKSHFSRQ").unwrap(),
        )
        .unwrap(),
        SequenceRecord::new(
            "fish",
            ProteinSequence::from_ascii(b"MKTAYIAKQRQISFVKSFSRQ").unwrap(),
        )
        .unwrap(),
    ];
    let dna_records = dna_records(&[("long", "ACGTACGT"), ("short", "ACGACGT")]);
    let rna_records = [
        SequenceRecord::new("long", RnaSequence::from_ascii(b"ACGUACGU").unwrap()).unwrap(),
        SequenceRecord::new("short", RnaSequence::from_ascii(b"ACGACGU").unwrap()).unwrap(),
    ];

    let protein = ClustalOmega::<Protein>::published()
        .align(&protein_records)
        .unwrap();
    let dna = ClustalOmega::<Dna>::published()
        .align(&dna_records)
        .unwrap();
    let rna = ClustalOmega::<Rna>::published()
        .align(&rna_records)
        .unwrap();

    for (result, expected) in [
        (
            protein.alignment().rows()[0].to_ungapped_ascii(),
            "MKTAYIAKQRQISFVKSHFSRQ",
        ),
        (
            protein.alignment().rows()[1].to_ungapped_ascii(),
            "MKTAYIAKQRQISFVKSFSRQ",
        ),
        (dna.alignment().rows()[0].to_ungapped_ascii(), "ACGTACGT"),
        (dna.alignment().rows()[1].to_ungapped_ascii(), "ACGACGT"),
        (rna.alignment().rows()[0].to_ungapped_ascii(), "ACGUACGU"),
        (rna.alignment().rows()[1].to_ungapped_ascii(), "ACGACGU"),
    ] {
        assert_eq!(result, expected);
    }
}

#[test]
fn full_guide_and_viterbi_pipeline_expose_auditable_evidence() {
    let records = dna_records(&[
        ("reference", "ACGTACGT"),
        ("deletion", "ACGACGT"),
        ("middle", "CGTACG"),
        ("ambiguity", "ACGNACGT"),
    ]);
    let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
        .with_guide_tree_mode(OmegaGuideTreeMode::Full)
        .with_engine_mode(OmegaEngineMode::Viterbi);
    let result = ClustalOmega::<Dna>::with_parameters(parameters)
        .unwrap()
        .align(&records)
        .unwrap();

    assert_eq!(result.guide_evidence().method(), OmegaGuideTreeMethod::Full);
    let distances = result.guide_evidence().distances().unwrap();
    assert_eq!(distances.size(), 4);
    assert!((distances.get(0, 1).unwrap() - 0.428_571_429).abs() < 1.0e-12);
    assert_eq!(distances.get(0, 2), Some(0.0));
    assert_eq!(distances.get(0, 3), Some(0.125));
    assert_eq!(result.guide_tree().leaf_count(), records.len());
    assert_eq!(result.progressive_steps().len(), records.len() - 1);
    assert!(result.progressive_steps().iter().all(|step| {
        step.engine() == OmegaEngine::Viterbi
            && step.score().is_finite()
            && step.posterior_bounds().is_none()
            && step.result_width() > 0
    }));
    let memory = result.memory_estimate();
    assert!(memory.guide_peak_bytes() > 0);
    assert!(memory.hmm_peak_bytes() > 0);
    assert!(memory.returned_alignment_bytes() > 0);
    assert!(memory.estimated_peak_bytes() >= memory.returned_alignment_bytes());
}

#[test]
fn mac_and_auto_memory_fallback_select_the_documented_engines() {
    let records = dna_records(&[("a", "ACGTAC"), ("b", "ACGAC"), ("c", "ACGTTC")]);
    let mac_parameters =
        ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Mac);
    let mac = ClustalOmega::<Dna>::with_parameters(mac_parameters)
        .unwrap()
        .align(&records)
        .unwrap();
    assert!(mac.progressive_steps().iter().all(|step| {
        step.engine() == OmegaEngine::Mac
            && step.posterior_bounds().is_some_and(|(minimum, maximum)| {
                (0.0..=1.0).contains(&minimum)
                    && (0.0..=1.0).contains(&maximum)
                    && minimum <= maximum
            })
    }));

    let fallback_parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
        .with_engine_mode(OmegaEngineMode::Auto)
        .with_mac_memory_limit_bytes(0);
    let fallback = ClustalOmega::<Dna>::with_parameters(fallback_parameters)
        .unwrap()
        .align(&records)
        .unwrap();
    assert!(
        fallback
            .progressive_steps()
            .iter()
            .all(|step| step.engine() == OmegaEngine::Viterbi)
    );
}

#[test]
fn mbed_path_caps_leaf_clusters_and_covers_every_record_once() {
    let records = dna_records(&[
        ("a", "AAAAAAAA"),
        ("c", "CCCCCCCC"),
        ("g", "GGGGGGGG"),
        ("t", "TTTTTTTT"),
        ("acgt", "ACGTACGT"),
        ("tgca", "TGCATGCA"),
        ("split", "AAAACCCC"),
    ]);
    let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
        .with_maximum_cluster_size(2)
        .unwrap()
        .with_engine_mode(OmegaEngineMode::Viterbi);
    let result = ClustalOmega::<Dna>::with_parameters(parameters)
        .unwrap()
        .align(&records)
        .unwrap();

    let evidence = result.guide_evidence();
    assert_eq!(evidence.method(), OmegaGuideTreeMethod::Mbed);
    assert!(evidence.distances().is_none());
    assert!(!evidence.seed_indices().is_empty());
    let clusters = evidence.clusters().collect::<Vec<_>>();
    assert!(
        clusters
            .iter()
            .all(|cluster| (1..=2).contains(&cluster.len()))
    );
    let mut members = clusters
        .iter()
        .flat_map(|cluster| cluster.iter().copied())
        .collect::<Vec<_>>();
    members.sort_unstable();
    assert_eq!(members, (0..records.len()).collect::<Vec<_>>());
}

#[test]
fn invalid_collections_profiles_and_parameters_fail_explicitly() {
    let empty: Vec<SequenceRecord<Dna>> = Vec::new();
    assert!(matches!(
        ClustalOmega::<Dna>::published().align(&empty),
        Err(ClustalOmegaError::Msa(MsaError::EmptyCollection))
    ));

    let protein_parameters = ClustalOmegaParameters::protein_v1_0_2();
    assert!(matches!(
        ClustalOmega::<Dna>::with_parameters(protein_parameters),
        Err(ClustalOmegaError::MoleculeProfileMismatch { .. })
    ));
    assert!(matches!(
        ClustalOmegaParameters::nucleotide_v1_2_4().with_maximum_cluster_size(1),
        Err(ClustalOmegaError::InvalidClusterSize { value: 1 })
    ));
    assert!(matches!(
        ClustalOmegaParameters::nucleotide_v1_2_4().with_thread_count(0),
        Err(ClustalOmegaError::InvalidThreadCount { value: 0 })
    ));
    assert!(matches!(
        ClustalOmegaParameters::nucleotide_v1_2_4().with_mac_threshold(f64::NAN),
        Err(ClustalOmegaError::InvalidParameter { .. })
    ));

    let empty_sequence = Sequence::<Dna>::from_ascii(b"").unwrap();
    let singleton = [SequenceRecord::new("empty", empty_sequence).unwrap()];
    let result = ClustalOmega::<Dna>::published().align(&singleton).unwrap();
    assert_eq!(result.alignment().sequence_count(), 1);
    assert_eq!(result.alignment().width(), 0);
    assert!(result.progressive_steps().is_empty());
}

#[test]
fn deterministic_parallel_pipeline_is_identical_to_one_thread() {
    let records = dna_records(&[
        (
            "s0",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ),
        (
            "s1",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
        ),
        (
            "s2",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGT",
        ),
        (
            "s3",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGT",
        ),
        (
            "s4",
            "TCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ),
        (
            "s5",
            "ACGTACGTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ),
        (
            "s6",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGTACGTACGTACGTACGTACGTACGT",
        ),
        (
            "s7",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGTACGTACGT",
        ),
    ]);
    let base = ClustalOmegaParameters::nucleotide_v1_2_4();
    let serial = ClustalOmega::<Dna>::with_parameters(base.clone())
        .unwrap()
        .align(&records)
        .unwrap();
    let parallel = ClustalOmega::<Dna>::with_parameters(base.with_thread_count(8).unwrap())
        .unwrap()
        .align(&records)
        .unwrap();

    assert_eq!(serial.alignment(), parallel.alignment());
    assert_eq!(serial.guide_tree(), parallel.guide_tree());
    assert_eq!(serial.progressive_steps(), parallel.progressive_steps());
    assert_eq!(
        serial.guide_evidence().distances(),
        parallel.guide_evidence().distances()
    );
}
