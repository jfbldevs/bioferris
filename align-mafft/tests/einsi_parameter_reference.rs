//! Opaque-reference evidence for documented ancestral E-INS-i parameters.

use align_mafft::{
    EINSI_LOCAL_OFFSET_PENALTY, EINSI_PAIRED_SKIP_EXTENSION, EINSI_SKIP_GAP_OPEN, EinsiParameters,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/e-ins-i-parameter-sweep-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    suite: String,
    command_template: String,
    provenance: Provenance,
    cases: Vec<Case>,
    inference: Inference,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    version: String,
    source_inspected: bool,
    package_sha256: String,
    generator_sha256: String,
    raw_output_sha256: String,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    local: Observation,
    generalized: Observation,
    #[serde(rename = "LOP")]
    generalized_open_sweep: Vec<SweepObservation>,
    #[serde(rename = "LEXP")]
    paired_extension_sweep: InvariantSweep,
    lep: Vec<SweepObservation>,
}

#[derive(Clone, Copy, Debug, Deserialize)]
struct Observation {
    alignment_width: usize,
    homologies: usize,
    tree_length: f64,
}

#[derive(Clone, Copy, Debug, Deserialize)]
struct SweepObservation {
    value: f64,
    alignment_width: usize,
    homologies: usize,
    tree_length: f64,
}

#[derive(Debug, Deserialize)]
struct InvariantSweep {
    values: Vec<f64>,
    alignment_width: usize,
    homologies: usize,
    tree_length: f64,
}

#[derive(Debug, Deserialize)]
struct Inference {
    #[serde(rename = "ancestral_LOP")]
    ancestral_generalized_open: f64,
    #[serde(rename = "ancestral_LEXP")]
    ancestral_paired_extension: f64,
    public_lep: f64,
    effective_signed_matrix_offset: f64,
}

fn same_observation(first: Observation, second: Observation) -> bool {
    first.alignment_width == second.alignment_width
        && first.homologies == second.homologies
        && first.tree_length.to_bits() == second.tree_length.to_bits()
}

#[test]
fn opaque_sweeps_support_the_versioned_generalized_parameter_mapping() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(
        fixture.suite,
        "MAFFT-7.526-E-INS-i-ancestral-parameter-black-box-v1"
    );
    assert_eq!(
        fixture.command_template,
        "mafft --genafpair --oldgenafpair --maxiterate 0 --thread 1 --inputorder --treeout --quiet --ep 0 [OPTION VALUE] INPUT"
    );
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.package_sha256,
        "b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de"
    );
    assert_eq!(
        fixture.provenance.generator_sha256,
        "49e58ce9b5afd1aaf3da8209a69c387ac1e71775d39fa05997a788723da113ae"
    );
    assert_eq!(
        fixture.provenance.raw_output_sha256,
        "93bc2c55de6bf1a418f4adbf6308b70850b065bf320206f024113c79ff3fc700"
    );
    assert_eq!(fixture.cases.len(), 4);

    for case in &fixture.cases {
        assert!(!case.id.is_empty());
        assert!(case.generalized.tree_length < case.local.tree_length);
        assert_eq!(case.generalized_open_sweep.len(), 5);
        for pair in case.generalized_open_sweep.windows(2) {
            assert!(pair[0].value > pair[1].value);
            assert!(pair[0].tree_length <= pair[1].tree_length);
        }
        assert!(
            case.generalized_open_sweep.first().unwrap().tree_length
                < case.generalized_open_sweep.last().unwrap().tree_length
        );
        let ancestral = case
            .generalized_open_sweep
            .iter()
            .find(|observation| observation.value.to_bits() == (-6.0_f64).to_bits())
            .unwrap();
        assert!(same_observation(
            case.generalized,
            Observation {
                alignment_width: ancestral.alignment_width,
                homologies: ancestral.homologies,
                tree_length: ancestral.tree_length,
            }
        ));

        assert_eq!(case.paired_extension_sweep.values.len(), 8);
        assert!(case.paired_extension_sweep.values.contains(&0.0));
        assert!(same_observation(
            case.generalized,
            Observation {
                alignment_width: case.paired_extension_sweep.alignment_width,
                homologies: case.paired_extension_sweep.homologies,
                tree_length: case.paired_extension_sweep.tree_length,
            }
        ));

        assert_eq!(case.lep.len(), 5);
        for pair in case.lep.windows(2) {
            assert!(pair[0].value < pair[1].value);
            assert!(pair[0].tree_length < pair[1].tree_length);
            assert!(pair[0].homologies >= pair[1].homologies);
        }
    }

    let parameters = EinsiParameters::paper_2008();
    assert_eq!(
        fixture.inference.ancestral_generalized_open.to_bits(),
        EINSI_SKIP_GAP_OPEN.to_bits()
    );
    assert_eq!(
        parameters.skip_gap_open().to_bits(),
        EINSI_SKIP_GAP_OPEN.to_bits()
    );
    assert_eq!(
        fixture.inference.ancestral_paired_extension.to_bits(),
        EINSI_PAIRED_SKIP_EXTENSION.to_bits()
    );
    assert_eq!(
        parameters.paired_skip_extension().to_bits(),
        EINSI_PAIRED_SKIP_EXTENSION.to_bits()
    );
    assert_eq!(
        fixture.inference.public_lep.to_bits(),
        EINSI_LOCAL_OFFSET_PENALTY.to_bits()
    );
    assert_eq!(
        fixture.inference.effective_signed_matrix_offset.to_bits(),
        (-EINSI_LOCAL_OFFSET_PENALTY).to_bits()
    );
}
