//! Independent model checks and retained public-reference numerical observations.

use bioferris::prelude::*;
use proptest::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};

const OBSERVATIONS: &str = include_str!("fixtures/protein-properties-biopython-1.88.json");

fn assert_same_bits(actual: f64, expected: f64) {
    assert_eq!(actual.to_bits(), expected.to_bits());
}

fn assert_observation(actual: f64, observation: &Value, context: &str) {
    let expected_bits = u64::from_str_radix(observation["bits"].as_str().unwrap(), 16).unwrap();
    assert_eq!(actual.to_bits(), expected_bits, "{context}: {actual:?}");
    // The public decimal and retained bits must describe the same number.
    // serde_json's float_roundtrip feature is necessary for this guarantee.
    assert_eq!(
        observation["value"].as_f64().unwrap().to_bits(),
        expected_bits,
        "reference decimal parser: {context}"
    );
}

#[test]
fn retained_reference_digest_and_metadata_are_pinned() {
    assert_eq!(
        format!("{:x}", Sha256::digest(OBSERVATIONS.as_bytes())),
        "d7d9e17cc81dc2bbe1955684799c41d6b003450c81a8f6efecb2ec5e8805d476"
    );
    let reference: Value = serde_json::from_str(OBSERVATIONS).unwrap();
    assert_eq!(reference["reference"], "Biopython");
    assert_eq!(reference["version"], "1.88");
    assert_eq!(reference["interpreter"], "3.12.3");
    assert_eq!(reference["mass_model"], "average linear unmodified peptide");
    assert_eq!(reference["gravy_scale"], "KyteDoolitle");
    assert_eq!(reference["cases"].as_array().unwrap().len(), 432);
}

#[test]
fn all_7344_retained_reference_values_match_exact_binary64_bits() {
    let reference: Value = serde_json::from_str(OBSERVATIONS).unwrap();
    let mut comparisons = 0;
    for case in reference["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let input = case["sequence"].as_str().unwrap();
        let sequence = protein(input).unwrap();
        let analysis = sequence.analyze_protein().unwrap();
        for (name, actual) in [
            ("molecular_weight", analysis.molecular_weight()),
            ("isoelectric_point", analysis.isoelectric_point()),
            ("aromaticity", analysis.aromaticity()),
            ("gravy", analysis.gravy()),
        ] {
            assert_observation(actual, &case[name], &format!("{id}/{name}"));
            comparisons += 1;
        }
        let charges = case["charges"].as_array().unwrap();
        assert_eq!(charges.len(), 13);
        for entry in charges {
            let ph = entry["ph"].as_f64().unwrap();
            let actual = analysis.charge_at_ph(ph).unwrap();
            assert_observation(actual, &entry["charge"], &format!("{id}/charge/{ph}"));
            assert_eq!(
                sequence
                    .protein_properties_with(ProteinPropertiesOptions { ph })
                    .unwrap()
                    .charge()
                    .to_bits(),
                actual.to_bits()
            );
            comparisons += 1;
        }
        assert_eq!(sequence.to_string(), input, "immutable source: {id}");
    }
    assert_eq!(comparisons, 7_344);
}

#[test]
fn aggregate_and_borrowed_common_paths_preserve_models_and_options() {
    let sequence = protein("INGAR").unwrap();
    let analysis = ProteinAnalysis::new(&sequence).unwrap();
    let result = sequence.protein_properties().unwrap();
    assert_eq!(result, analysis.properties());
    assert_eq!(result.len(), 5);
    assert!(!result.is_empty());
    assert_eq!(analysis.len(), 5);
    assert!(!analysis.is_empty());
    assert_same_bits(result.ph(), 7.0);
    assert_same_bits(result.molecular_weight(), analysis.molecular_weight());
    assert_same_bits(result.charge(), analysis.charge_at_ph(7.0).unwrap());
    assert_same_bits(result.aromaticity(), analysis.aromaticity());
    assert_same_bits(result.gravy(), analysis.gravy());
    assert_eq!(result.pi_estimate(), analysis.isoelectric_point_estimate());
    assert_same_bits(result.isoelectric_point(), result.pi_estimate().value());
    assert_same_bits(ProteinPropertiesOptions::default().ph, 7.0);
}

#[test]
fn average_mass_subtracts_one_water_per_peptide_bond() {
    let sequence = protein("AG").unwrap();
    let result = sequence.protein_properties().unwrap();
    // Two independent free-residue masses minus one water molecule, in Da.
    assert!((result.molecular_weight() - 146.1445).abs() < 1e-12);
    let alanine = protein("A").unwrap().protein_properties().unwrap();
    assert_same_bits(alanine.molecular_weight(), 89.0932);
    assert_same_bits(alanine.gravy(), 1.8);
    assert_same_bits(alanine.aromaticity(), 0.0);
}

#[test]
fn aromaticity_and_gravy_have_independent_literal_oracles() {
    let result = protein("AFWY").unwrap().protein_properties().unwrap();
    assert_same_bits(result.aromaticity(), 0.75);
    assert!((result.gravy() - 0.6).abs() < 1e-15);
    let hydrophilic = protein("DDDD").unwrap().protein_properties().unwrap();
    assert_same_bits(hydrophilic.gravy(), -3.5);
}

#[test]
fn bounded_pi_exposes_nonbracketed_acidic_and_basic_roots() {
    for (input, boundary, bound) in [
        ("D".repeat(100), PiBoundary::Lower, 4.05),
        ("R".repeat(100), PiBoundary::Upper, 12.0),
    ] {
        let result = protein(input).unwrap().protein_properties().unwrap();
        let estimate = result.pi_estimate();
        assert_eq!(estimate.boundary(), boundary);
        assert!(!estimate.bracketed());
        assert!((estimate.value() - bound).abs() < 0.0001);
        assert_same_bits(result.isoelectric_point(), estimate.value());
    }
    let result = protein("INGAR").unwrap().protein_properties().unwrap();
    assert_eq!(result.pi_estimate().boundary(), PiBoundary::None);
    assert!(result.pi_estimate().bracketed());
}

#[test]
fn each_general_protein_ambiguity_or_noncanonical_residue_is_explicitly_rejected() {
    for residue in ['B', 'Z', 'J', 'X', 'U', 'O'] {
        let sequence = protein(format!("AG{residue}F")).unwrap();
        let error = sequence.protein_properties().unwrap_err();
        assert!(
            matches!(error, ProteinPropertyError::UnsupportedResidue { offset: 2, byte } if byte == u8::try_from(residue).unwrap())
        );
        assert_eq!(error.code(), "unsupported_protein_residue");
        assert!(!error.to_string().is_empty());
    }
}

#[test]
fn empty_input_is_not_a_zero_mass_or_division_by_zero_result() {
    let error = protein("").unwrap().protein_properties().unwrap_err();
    assert_eq!(error, ProteinPropertyError::EmptySequence);
    assert_eq!(error.code(), "empty_protein");
    assert!(std::error::Error::source(&error).is_none());
}

#[test]
fn invalid_ph_is_rejected_without_clamping_or_nan_results() {
    let sequence = protein("PETER").unwrap();
    let analysis = sequence.analyze_protein().unwrap();
    for ph in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.001, 14.001] {
        let error = analysis.charge_at_ph(ph).unwrap_err();
        assert!(matches!(error, ProteinPropertyError::InvalidPh { .. }));
        assert_eq!(error.code(), "invalid_ph");
        assert!(
            analysis
                .properties_with(ProteinPropertiesOptions { ph })
                .is_err()
        );
    }
    for ph in [0.0, -0.0, 14.0] {
        let result = analysis
            .properties_with(ProteinPropertiesOptions { ph })
            .unwrap();
        assert_eq!(result.ph().to_bits(), ph.to_bits());
        assert!(result.charge().is_finite());
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn charge_is_monotonic_and_results_are_finite(
        residues in prop::collection::vec(0usize..20, 1..300),
        first in 0.0f64..14.0,
        second in 0.0f64..14.0,
    ) {
        let canonical = b"ACDEFGHIKLMNPQRSTVWY";
        let input: String = residues.iter().map(|&i| char::from(canonical[i])).collect();
        let sequence = protein(&input).unwrap();
        let analysis = sequence.analyze_protein().unwrap();
        let (lower, upper) = (first.min(second), first.max(second));
        prop_assert!(analysis.charge_at_ph(lower).unwrap() >= analysis.charge_at_ph(upper).unwrap());
        let result = analysis.properties();
        prop_assert!(result.molecular_weight().is_finite() && result.molecular_weight() > 0.0);
        prop_assert!((0.0..=1.0).contains(&result.aromaticity()));
        prop_assert!(result.gravy().is_finite());
        prop_assert!((4.05..=12.0).contains(&result.isoelectric_point()));
        prop_assert_eq!(sequence.to_string(), input);
    }
}
