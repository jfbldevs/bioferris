//! Exact replay of a synthetic modern MAFFT 7.526 E-INS-i case.

use align_core::ProteinSequence;
use align_mafft_compat::{MAFFT_7526_EINSI_PROFILE_IDENTIFIER, Mafft7526Einsi};
use align_msa::SequenceRecord;
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/einsi-modern-protein-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    profile: String,
    command: String,
    provenance: Provenance,
    case: Case,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    version: String,
    package_sha256: String,
    synthetic_seed: u64,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    sequences: Vec<SequenceValue>,
    reference_alignment: Vec<AlignedValue>,
}

#[derive(Debug, Deserialize)]
struct SequenceValue {
    id: String,
    sequence: String,
}

#[derive(Debug, Deserialize)]
struct AlignedValue {
    id: String,
    aligned: String,
}

#[test]
fn modern_generalized_affine_protein_case_is_exact() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture.profile, MAFFT_7526_EINSI_PROFILE_IDENTIFIER);
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert_eq!(fixture.provenance.synthetic_seed, 7_526);
    assert_eq!(
        fixture.provenance.package_sha256,
        "b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de"
    );
    assert_eq!(
        fixture.command,
        "mafft --genafpair --ep 0 --maxiterate 1000 --inputorder --quiet INPUT"
    );
    let records = fixture
        .case
        .sequences
        .iter()
        .map(|value| {
            SequenceRecord::new(
                value.id.clone(),
                ProteinSequence::from_ascii(value.sequence.as_bytes()).unwrap(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let actual = Mafft7526Einsi::new()
        .align(&records)
        .unwrap()
        .alignment()
        .rows()
        .iter()
        .map(|row| (row.identifier().to_owned(), row.to_ascii_with_gaps()))
        .collect::<Vec<_>>();
    let expected = fixture
        .case
        .reference_alignment
        .iter()
        .map(|row| (row.id.clone(), row.aligned.clone()))
        .collect::<Vec<_>>();
    assert_eq!(actual, expected, "{} diverged", fixture.case.id);
}
