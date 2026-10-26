//! Public motif workflows and integrity of retained reference observations.

use bioferris::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};

#[test]
fn retained_matrix_reference_digest_and_runtime_metadata_are_pinned() {
    let observations = include_str!("../tests/fixtures/motifs-biopython-1.88.json");
    assert_eq!(
        format!("{:x}", Sha256::digest(observations.as_bytes())),
        "75e37c85cd6370b2b17a9773e5226a459071f2156d6b2a3285d7321f25f5c824"
    );
    let reference: Value = serde_json::from_str(observations).unwrap();
    assert_eq!(reference["schema_version"], 1);
    assert_eq!(reference["profile"], "biopython-1.88-motif-matrices-v1");
    assert_eq!(reference["reference"]["biopython"], "1.88");
    assert_eq!(reference["reference"]["numpy"], "2.3.3");
    assert_eq!(reference["reference"]["python"], "3.12.3");
    assert_eq!(reference["cases"].as_array().unwrap().len(), 360);
}

#[test]
fn prelude_iupac_conveniences_preserve_overlap_strand_and_input() {
    let sequence = dna("AGRYN").unwrap();
    let motif: DnaMotif = dna_motif("r").unwrap();
    assert_eq!(motif.pattern(), "R");
    assert_eq!(motif.width(), 1);
    let hits = motif.find_in(&sequence).unwrap();
    let observed: Vec<_> = hits.iter().map(|hit| (hit.range(), hit.strand())).collect();
    assert_eq!(
        observed,
        vec![
            (0..1, Strand::Forward),
            (1..2, Strand::Forward),
            (2..3, Strand::Forward),
            (3..4, Strand::Reverse),
        ]
    );
    assert_eq!(sequence.to_string(), "AGRYN");
    let overlaps = dna_motif("AA")
        .unwrap()
        .find_with(
            &dna("AAAA").unwrap(),
            MotifSearchOptions {
                strands: StrandSelection::Forward,
                ..MotifSearchOptions::default()
            },
        )
        .unwrap();
    assert_eq!(
        overlaps.iter().map(MotifMatch::range).collect::<Vec<_>>(),
        vec![0..2, 1..3, 2..4]
    );
}

#[test]
fn public_count_probability_score_pipeline_has_literal_independent_results() {
    let instances = ["AC", "AC", "AC", "CA"].map(|input| dna(input).unwrap());
    let counts = MotifCounts::from_instances(&instances).unwrap();
    assert_eq!(
        counts.columns(),
        &[[3.0, 1.0, 0.0, 0.0], [1.0, 3.0, 0.0, 0.0]]
    );
    let weights = counts.normalize([1.0; 4]).unwrap();
    assert_eq!(
        weights.columns(),
        &[[0.5, 0.25, 0.125, 0.125], [0.25, 0.5, 0.125, 0.125]]
    );
    let scores = weights.log_odds_uniform().unwrap();
    assert_eq!(scores.width(), 2);
    assert_eq!(
        scores.background().map(f64::to_bits),
        [0.25_f64.to_bits(); 4]
    );
    let sequence = dna("ACGTAC").unwrap();
    let scan = scores.calculate(&sequence, 5).unwrap();
    assert_eq!(scan, vec![2.0, -1.0, -2.0, -1.0, 2.0]);
    let hits = scores
        .find_with(
            &sequence,
            MatrixSearchOptions {
                threshold: 2.0,
                ..MatrixSearchOptions::default()
            },
        )
        .unwrap();
    assert_eq!(
        hits.iter()
            .map(|hit| (hit.range(), hit.strand(), hit.score().to_bits()))
            .collect::<Vec<_>>(),
        vec![
            (0..2, Strand::Forward, 2.0_f32.to_bits()),
            (2..4, Strand::Reverse, 2.0_f32.to_bits()),
            (4..6, Strand::Forward, 2.0_f32.to_bits()),
        ]
    );
    assert_eq!(sequence.to_string(), "ACGTAC");
}

#[test]
fn rna_uses_uracil_without_changing_count_or_scan_semantics() {
    let dna_instances = ["AC", "AT", "GC"].map(|input| dna(input).unwrap());
    let rna_instances = ["AC", "AU", "GC"].map(|input| rna(input).unwrap());
    let dna_counts = MotifCounts::from_instances(&dna_instances).unwrap();
    let rna_counts = MotifCounts::from_instances(&rna_instances).unwrap();
    assert_eq!(dna_counts.columns(), rna_counts.columns());
    let dna_scores = dna_counts
        .normalize([0.5; 4])
        .unwrap()
        .log_odds([0.4, 0.3, 0.2, 0.1])
        .unwrap();
    let rna_scores = rna_counts
        .normalize([0.5; 4])
        .unwrap()
        .log_odds([0.4, 0.3, 0.2, 0.1])
        .unwrap();
    let options = MatrixSearchOptions {
        threshold: -100.0,
        ..MatrixSearchOptions::default()
    };
    let dna_hits = dna_scores
        .find_with(&dna("ACGTAC").unwrap(), options)
        .unwrap();
    let rna_hits = rna_scores
        .find_with(&rna("ACGUAC").unwrap(), options)
        .unwrap();
    assert_eq!(dna_hits.len(), rna_hits.len());
    for (dna_hit, rna_hit) in dna_hits.iter().zip(&rna_hits) {
        assert_eq!(dna_hit.range(), rna_hit.range());
        assert_eq!(dna_hit.strand(), rna_hit.strand());
        assert_eq!(dna_hit.score().to_bits(), rna_hit.score().to_bits());
    }
    let pattern: RnaMotif = rna_motif("AUG").unwrap();
    assert!(rna_motif("ATG").is_err());
    let hits = pattern.find_in(&rna("AUGCAU").unwrap()).unwrap();
    assert_eq!(hits[0].range(), 0..3);
    assert_eq!(hits[0].strand(), Strand::Forward);
    assert_eq!(hits[1].range(), 3..6);
    assert_eq!(hits[1].strand(), Strand::Reverse);
}

#[test]
fn asymmetrical_background_reverse_scores_use_the_built_matrix() {
    let scores = MotifCounts::<Dna>::from_columns(vec![[1.0; 4]])
        .unwrap()
        .normalize_default()
        .unwrap()
        .log_odds([0.5, 0.25, 0.125, 0.125])
        .unwrap();
    let hits = scores
        .find_with(
            &dna("A").unwrap(),
            MatrixSearchOptions {
                strands: StrandSelection::Reverse,
                threshold: 1.0,
                ..MatrixSearchOptions::default()
            },
        )
        .unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].range(), 0..1);
    assert_eq!(hits[0].score().to_bits(), 1.0_f32.to_bits());
}

#[test]
fn public_negative_infinity_and_budget_failures_are_explicit() {
    let scores = MotifCounts::<Dna>::from_columns(vec![[1.0, 0.0, 0.0, 0.0]])
        .unwrap()
        .normalize_default()
        .unwrap()
        .log_odds_uniform()
        .unwrap();
    assert_eq!(
        scores.score_window(&dna("C").unwrap()).unwrap().to_bits(),
        f32::NEG_INFINITY.to_bits()
    );
    let error = scores
        .find_with(
            &dna("AA").unwrap(),
            MatrixSearchOptions {
                strands: StrandSelection::Forward,
                max_matches: 1,
                ..MatrixSearchOptions::default()
            },
        )
        .unwrap_err();
    assert!(matches!(
        error,
        MotifError::LimitExceeded {
            resource: "matches",
            limit: 1,
            observed: 2,
        }
    ));
    assert!(matches!(
        scores.calculate(&dna("AN").unwrap(), 0),
        Err(MotifError::AmbiguousInput {
            offset: 1,
            byte: b'N'
        })
    ));
}
