//! Static G-INS-i comparisons against the official opaque MAFFT 7.526 executable.

use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};
use align_mafft::{MafftAlphabet, MafftGinsi};
use align_msa::{
    AlignedRow, AlignmentMetrics, MultipleAlignment, SequenceRecord, compare_alignments,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/g-ins-i-globalpair-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    profile: String,
    provenance: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    version: String,
    source_inspected: bool,
    package_sha256: String,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    alphabet: String,
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

fn run_case<A: MafftAlphabet>(case: &Case) -> (AlignmentMetrics, bool) {
    let records = case
        .sequences
        .iter()
        .map(|value| {
            SequenceRecord::new(
                value.id.clone(),
                Sequence::<A>::from_ascii(value.sequence.as_bytes()).unwrap(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let candidate = MafftGinsi::<A>::new().align(&records).unwrap();
    let reference_rows = case
        .reference_alignment
        .iter()
        .map(|value| {
            let cells = value
                .aligned
                .bytes()
                .map(|byte| {
                    if byte == b'-' {
                        None
                    } else {
                        Some(Symbol::<A>::from_ascii(byte).expect("fixture symbol is valid"))
                    }
                })
                .collect();
            AlignedRow::new(value.id.clone(), cells).unwrap()
        })
        .collect();
    let reference = MultipleAlignment::try_new(reference_rows).unwrap();
    let metrics = compare_alignments(candidate.alignment(), &reference).unwrap();
    let exact = render(candidate.alignment()) == render(&reference);
    (metrics, exact)
}

fn render<A: Alphabet>(alignment: &MultipleAlignment<A>) -> Vec<(String, String)> {
    alignment
        .rows()
        .iter()
        .map(|row| (row.identifier().to_owned(), row.to_ascii_with_gaps()))
        .collect()
}

#[test]
#[allow(clippy::print_stderr)]
fn official_mafft_7526_ginsi_metrics_are_recorded_and_gated() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture.profile, "MAFFT-7.526-G-INS-i-globalpair-v1");
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.package_sha256,
        "7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7"
    );

    let mut sp_sum = 0.0;
    let mut tc_sum = 0.0;
    let mut exact_count = 0usize;
    for case in &fixture.cases {
        let (metrics, exact) = match case.alphabet.as_str() {
            "dna" => run_case::<Dna>(case),
            "rna" => run_case::<Rna>(case),
            "protein" => run_case::<Protein>(case),
            other => panic!("unexpected fixture alphabet {other}"),
        };
        eprintln!(
            "{}: SP={:.9}, TC={:.9}, exact={exact}",
            case.id,
            metrics.sp_score(),
            metrics.tc_score()
        );
        sp_sum += metrics.sp_score();
        tc_sum += metrics.tc_score();
        exact_count += usize::from(exact);
        if case.id.ends_with("identical") {
            assert!(exact, "{} should be version-independent", case.id);
        } else {
            assert!(
                metrics.sp_score() >= 0.90,
                "{} SP regression: {}",
                case.id,
                metrics.sp_score()
            );
            assert!(
                metrics.tc_score() >= 0.80,
                "{} TC regression: {}",
                case.id,
                metrics.tc_score()
            );
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let denominator = fixture.cases.len() as f64;
    let mean_sp = sp_sum / denominator;
    let mean_tc = tc_sum / denominator;
    eprintln!("mean: SP={mean_sp:.9}, TC={mean_tc:.9}, exact={exact_count}");
    assert!(mean_sp >= 0.95, "mean SP regression: {mean_sp}");
    assert!(mean_tc >= 0.85, "mean TC regression: {mean_tc}");
    assert!(exact_count >= 2);
}
