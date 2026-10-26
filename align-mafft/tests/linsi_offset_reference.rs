//! Opaque-reference evidence for the sign of L-INS-i's public `lep` value.

use align_mafft::{LINSI_LOCAL_OFFSET_PENALTY, LinsiParameters};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/l-ins-i-offset-sign-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    suite: String,
    command_template: String,
    provenance: Provenance,
    case: Case,
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
    observations: Vec<Observation>,
}

#[derive(Clone, Copy, Debug, Deserialize)]
struct Observation {
    lep: f64,
    alignment_width: usize,
    total_pair_homologies: usize,
    tree_branch_length_sum: f64,
}

#[derive(Debug, Deserialize)]
struct Inference {
    public_lep: f64,
    effective_signed_matrix_offset: f64,
}

#[test]
fn opaque_mafft_sweep_supports_subtracting_the_positive_lep_penalty() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture.suite, "MAFFT-7.526-L-INS-i-lep-sign-black-box-v1");
    assert_eq!(
        fixture.command_template,
        "mafft --localpair --maxiterate 0 --thread 1 --inputorder --treeout --quiet --lep VALUE INPUT"
    );
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.package_sha256,
        "b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de"
    );
    assert_eq!(
        fixture.provenance.generator_sha256,
        "6b912ba8dc5f419ca4873e7e3a9e8eb529c477cefe981e0484a0763f8eb72838"
    );
    assert_eq!(
        fixture.provenance.raw_output_sha256,
        "9bccaa57c0e623f929f1675fed4520020ea3f22af79ba8952eac473182884f8c"
    );
    assert_eq!(fixture.case.id, "short_core_adversarial_flanks");
    assert_eq!(fixture.case.observations.len(), 7);

    for pair in fixture.case.observations.windows(2) {
        assert!(pair[0].lep < pair[1].lep);
        assert!(
            pair[0].tree_branch_length_sum < pair[1].tree_branch_length_sum,
            "tree length must increase as the public offset penalty increases"
        );
        assert!(pair[0].total_pair_homologies >= pair[1].total_pair_homologies);
    }
    let negative = fixture.case.observations.first().unwrap();
    let public_default = fixture
        .case
        .observations
        .iter()
        .find(|observation| observation.lep.to_bits() == 0.1_f64.to_bits())
        .unwrap();
    assert!(negative.alignment_width < public_default.alignment_width);
    assert!(negative.total_pair_homologies > public_default.total_pair_homologies);

    assert_eq!(
        fixture.inference.public_lep.to_bits(),
        LINSI_LOCAL_OFFSET_PENALTY.to_bits()
    );
    assert_eq!(
        LinsiParameters::paper_2005()
            .local_offset_penalty()
            .to_bits(),
        LINSI_LOCAL_OFFSET_PENALTY.to_bits()
    );
    assert_eq!(
        fixture.inference.effective_signed_matrix_offset.to_bits(),
        (-LINSI_LOCAL_OFFSET_PENALTY).to_bits()
    );
}
