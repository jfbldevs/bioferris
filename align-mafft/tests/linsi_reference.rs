//! Static L-INS-i comparisons against the official opaque MAFFT 7.526 executable.

use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};
use align_mafft::{MafftAlphabet, MafftLinsi};
use align_msa::{
    AlignedRow, AlignmentMetrics, MultipleAlignment, SequenceRecord, compare_alignments,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/l-ins-i-localpair-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    suite: String,
    profile: String,
    command: String,
    provenance: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    version: String,
    source_inspected: bool,
    package_sha256: String,
    generator_sha256: String,
    case_specification_sha256: String,
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
    let candidate = MafftLinsi::<A>::new().align(&records).unwrap();
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
fn official_mafft_7526_linsi_metrics_are_recorded_and_gated() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture.suite, "MAFFT-7.526-L-INS-i-localpair-v1-black-box");
    assert_eq!(fixture.profile, "MAFFT-7.526-L-INS-i-localpair-v1");
    assert_eq!(
        fixture.command,
        "mafft --localpair --maxiterate 1000 --thread 1 --inputorder --quiet INPUT"
    );
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.package_sha256,
        "b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de"
    );
    assert_eq!(
        fixture.provenance.generator_sha256,
        "ab1a670183ce7ef70789da301f76221f5ed22fab1144134a37229717e87021cc"
    );
    assert_eq!(
        fixture.provenance.case_specification_sha256,
        "ad51ebccfd9214343d09a16c436180ae061c1bbf0bf0dc0439163d51c118e27e"
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
            assert!(metrics.sp_score() >= 0.96, "{} SP regression", case.id);
            assert!(metrics.tc_score() >= 0.92, "{} TC regression", case.id);
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let denominator = fixture.cases.len() as f64;
    let mean_sp = sp_sum / denominator;
    let mean_tc = tc_sum / denominator;
    eprintln!("mean: SP={mean_sp:.9}, TC={mean_tc:.9}, exact={exact_count}");
    assert!(mean_sp >= 0.99, "mean SP regression: {mean_sp}");
    assert!(mean_tc >= 0.98, "mean TC regression: {mean_tc}");
    assert!(exact_count >= 6);
}
