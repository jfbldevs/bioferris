//! Static E-INS-i comparisons against the official opaque MAFFT 7.526 executable.

use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};
use align_mafft::{MafftAlphabet, MafftEinsi};
use align_msa::{
    AlignedRow, AlignmentMetrics, MultipleAlignment, SequenceRecord, compare_alignments,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/mafft-7.526/e-ins-i-oldgenafpair-v1.json");

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
    let candidate = MafftEinsi::<A>::new().align(&records).unwrap();
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
fn official_mafft_7526_oldgenafpair_metrics_are_recorded_and_gated() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(
        fixture.suite,
        "MAFFT-7.526-E-INS-i-oldgenafpair-v1-black-box"
    );
    assert_eq!(fixture.profile, "MAFFT-7.526-E-INS-i-oldgenafpair-v1");
    assert_eq!(
        fixture.command,
        "mafft --genafpair --oldgenafpair --ep 0 --maxiterate 1000 --thread 1 --inputorder --quiet INPUT"
    );
    assert_eq!(fixture.provenance.version, "v7.526 (2024/Apr/26)");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.package_sha256,
        "b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de"
    );
    assert_eq!(
        fixture.provenance.generator_sha256,
        "7df575e8b9ed5d17fd823aced10b80a7935cc387a0ab44cf0479156a192e0ff9"
    );
    assert_eq!(
        fixture.provenance.case_specification_sha256,
        "cb91b85398dfad734cb0f552c4107dd8067f94b8130f9b231b16cb183aa9afe4"
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
        assert!(metrics.sp_score().is_finite());
        assert!(metrics.tc_score().is_finite());
        let (minimum_sp, minimum_tc, exact_required) = match case.id.as_str() {
            "dna_identical"
            | "rna_identical"
            | "protein_identical"
            | "dna_two_motifs_equal_internal" => (1.0, 1.0, true),
            // The modern executable's ancestral nucleotide route separates
            // entire repeats in this adversarial unequal-length case, whereas
            // the versioned Rust recurrence retains both conserved motifs.
            "dna_two_motifs_unequal_internal" => (0.52, 0.40, false),
            "rna_two_motifs_unequal_internal" => (0.97, 0.94, false),
            "protein_two_domains_equal_internal" => (0.99, 0.60, false),
            "protein_two_domains_unequal_internal" => (0.94, 0.65, false),
            other => panic!("unversioned E-INS-i reference case {other}"),
        };
        assert!(
            metrics.sp_score() >= minimum_sp,
            "{} SP regression: {} < {minimum_sp}",
            case.id,
            metrics.sp_score()
        );
        assert!(
            metrics.tc_score() >= minimum_tc,
            "{} TC regression: {} < {minimum_tc}",
            case.id,
            metrics.tc_score()
        );
        if exact_required {
            assert!(exact, "{} should be exactly reproducible", case.id);
        }
        sp_sum += metrics.sp_score();
        tc_sum += metrics.tc_score();
        exact_count += usize::from(exact);
    }
    #[allow(clippy::cast_precision_loss)]
    let denominator = fixture.cases.len() as f64;
    let mean_sp = sp_sum / denominator;
    let mean_tc = tc_sum / denominator;
    eprintln!("mean: SP={mean_sp:.9}, TC={mean_tc:.9}, exact={exact_count}");
    assert!(mean_sp >= 0.93, "mean SP regression: {mean_sp}");
    assert!(mean_tc >= 0.82, "mean TC regression: {mean_tc}");
    assert!(exact_count >= 4);
}
