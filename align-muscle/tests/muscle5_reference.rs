//! Inert black-box comparison against the pinned MUSCLE 5.1 executable.
//!
//! The executable and generator remain in the isolated sibling laboratory. This
//! test only decodes reviewed JSON; it never starts Python, WSL, or MUSCLE.

use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord, compare_alignments};
use align_muscle::{
    GuideTreePermutation, Muscle5, Muscle5Alphabet, Muscle5Mode, Muscle5Parameters,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/muscle-5.1/nature-default-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    schema_version: u32,
    suite: String,
    profile: String,
    provenance: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    name: String,
    version: String,
    source_inspected: bool,
    binary_sha256: String,
    generator_sha256: String,
    case_specification_sha256: String,
    platform: String,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    mode: String,
    alphabet: String,
    permutation: String,
    seed: u64,
    sequences: Vec<InputRecord>,
    reference_alignment: Vec<ReferenceRow>,
}

#[derive(Debug, Deserialize)]
struct InputRecord {
    id: String,
    sequence: String,
}

#[derive(Debug, Deserialize)]
struct ReferenceRow {
    id: String,
    aligned: String,
}

#[derive(Clone, Copy, Debug)]
struct ObservedMetrics {
    sp: f64,
    tc: f64,
    exact: bool,
}

fn acceptance(case_id: &str) -> (f64, f64, bool) {
    match case_id {
        "align_dna_internal_terminal" | "align_rna_internal_terminal" => (0.68, 0.40, false),
        "align_dna_repetitive" => (0.73, 0.46, false),
        "super5_dna" | "super5_rna" => (0.91, 0.71, false),
        "super5_divergent_dna" => (0.40, 0.0, false),
        "align_protein_globins"
        | "perm_none_protein"
        | "perm_abc_protein"
        | "perm_acb_protein"
        | "perm_bca_protein"
        | "perturb_1_protein"
        | "perturb_7_bca_protein"
        | "super5_protein"
        | "super5_redundant_protein" => (1.0, 1.0, true),
        _ => panic!("case {case_id}: acceptance profile is absent"),
    }
}

fn reference_alignment<A: Alphabet>(case: &Case) -> MultipleAlignment<A> {
    let rows = case
        .sequences
        .iter()
        .map(|source| {
            let row = case
                .reference_alignment
                .iter()
                .find(|row| row.id == source.id)
                .unwrap_or_else(|| {
                    panic!("case {}: reference row {} is absent", case.id, source.id)
                });
            let symbols = row
                .aligned
                .bytes()
                .map(|byte| {
                    if byte == b'-' {
                        None
                    } else {
                        Some(Symbol::<A>::from_ascii(byte).unwrap_or_else(|| {
                            panic!(
                                "case {}: invalid {} symbol {byte:?} in reference row {}",
                                case.id,
                                A::NAME,
                                row.id
                            )
                        }))
                    }
                })
                .collect();
            AlignedRow::new(&row.id, symbols).unwrap_or_else(|error| {
                panic!(
                    "case {}: invalid reference row {}: {error}",
                    case.id, row.id
                )
            })
        })
        .collect();
    MultipleAlignment::try_new(rows)
        .unwrap_or_else(|error| panic!("case {}: invalid reference alignment: {error}", case.id))
}

fn records<A: Alphabet>(case: &Case) -> Vec<SequenceRecord<A>> {
    case.sequences
        .iter()
        .map(|record| {
            SequenceRecord::new(
                &record.id,
                Sequence::<A>::from_ascii(record.sequence.as_bytes())
                    .unwrap_or_else(|error| panic!("case {}: invalid input: {error}", case.id)),
            )
            .unwrap_or_else(|error| panic!("case {}: invalid identifier: {error}", case.id))
        })
        .collect()
}

fn permutation(label: &str) -> GuideTreePermutation {
    match label {
        "none" => GuideTreePermutation::None,
        "abc" => GuideTreePermutation::Abc,
        "acb" => GuideTreePermutation::Acb,
        "bca" => GuideTreePermutation::Bca,
        _ => panic!("unknown guide-tree permutation {label}"),
    }
}

fn run_case<A: Muscle5Alphabet>(case: &Case) -> ObservedMetrics {
    let records = records::<A>(case);
    let reference = reference_alignment::<A>(case);
    for source in &case.sequences {
        assert_eq!(
            reference
                .row_by_identifier(&source.id)
                .expect("reviewed reference row exists")
                .to_ungapped_ascii(),
            source.sequence
        );
    }
    let mode = match case.mode.as_str() {
        "align" => Muscle5Mode::Align,
        "super5" => Muscle5Mode::Super5,
        other => panic!("case {}: unknown mode {other}", case.id),
    };
    let parameters = if mode == Muscle5Mode::Align {
        Muscle5Parameters::published_align()
    } else {
        Muscle5Parameters::published_super5()
    }
    .with_perturbation_seed(case.seed)
    .with_guide_tree_permutation(permutation(&case.permutation));
    let result = Muscle5::<A>::with_parameters(parameters)
        .unwrap_or_else(|error| panic!("case {}: invalid profile: {error}", case.id))
        .align(&records)
        .unwrap_or_else(|error| panic!("case {}: Rust alignment failed: {error}", case.id));
    let metrics = compare_alignments(result.alignment(), &reference)
        .unwrap_or_else(|error| panic!("case {}: metric comparison failed: {error}", case.id));
    ObservedMetrics {
        sp: metrics.sp_score(),
        tc: metrics.tc_score(),
        exact: result.alignment() == &reference,
    }
}

#[test]
fn publication_profile_agrees_with_isolated_muscle5_fixture() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("fixture JSON must be valid");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(
        fixture.suite,
        "MUSCLE-5.1.0-Nature-2022-default-v1-black-box"
    );
    // The fixture was recorded for v1. Version 2 changed only how refinement
    // scores the current alignment, so it is held to the same acceptance.
    assert_eq!(fixture.profile, "MUSCLE-5.1.0-Nature-2022-default-v1");
    assert_eq!(fixture.provenance.name, "MUSCLE");
    assert_eq!(fixture.provenance.version, "muscle 5.1.linux64 [12f0e2]");
    assert!(!fixture.provenance.source_inspected);
    assert_eq!(
        fixture.provenance.binary_sha256,
        "2bd82928a5cc72b5086b60928a3a820a0ce066ca9fbc0fecd1543eebca981e76"
    );
    assert_eq!(
        fixture.provenance.generator_sha256,
        "1124b017be24cb63405a3a088abea46fd7254bf904a14528227f6c4b04fd8277"
    );
    assert_eq!(
        fixture.provenance.case_specification_sha256,
        "b2e8cda7286208191c11c30021d2c544a75639bb0520e479a25d45e90849f400"
    );
    assert!(fixture.provenance.platform.contains("WSL2"));
    assert_eq!(fixture.cases.len(), 15);

    let mut sp_total = 0.0_f64;
    let mut tc_total = 0.0_f64;
    let mut exact_count = 0_usize;
    let mut failures = Vec::new();
    for case in &fixture.cases {
        let observed = match case.alphabet.as_str() {
            "dna" => run_case::<Dna>(case),
            "rna" => run_case::<Rna>(case),
            "protein" => run_case::<Protein>(case),
            alphabet => panic!("case {}: unknown alphabet {alphabet}", case.id),
        };
        println!(
            "{} mode={} perm={} seed={} SP={:.9} TC={:.9} exact={}",
            case.id,
            case.mode,
            case.permutation,
            case.seed,
            observed.sp,
            observed.tc,
            observed.exact
        );
        sp_total += observed.sp;
        tc_total += observed.tc;
        exact_count += usize::from(observed.exact);
        let (minimum_sp, minimum_tc, require_exact) = acceptance(&case.id);
        if observed.sp + f64::EPSILON < minimum_sp {
            failures.push(format!(
                "{}: SP {} is below {}",
                case.id, observed.sp, minimum_sp
            ));
        }
        if observed.tc + f64::EPSILON < minimum_tc {
            failures.push(format!(
                "{}: TC {} is below {}",
                case.id, observed.tc, minimum_tc
            ));
        }
        if require_exact && !observed.exact {
            failures.push(format!("{}: exact agreement is required", case.id));
        }
    }
    let count = f64::from(u32::try_from(fixture.cases.len()).expect("case count fits u32"));
    let mean_sp = sp_total / count;
    let mean_tc = tc_total / count;
    println!(
        "MUSCLE 5.1 black-box summary: cases={}, exact={exact_count}, mean_SP={mean_sp:.9}, mean_TC={mean_tc:.9}",
        fixture.cases.len()
    );

    if mean_sp < 0.88 {
        failures.push(format!("mean SP {mean_sp} is below 0.88"));
    }
    if mean_tc < 0.77 {
        failures.push(format!("mean TC {mean_tc} is below 0.77"));
    }
    if exact_count < 9 {
        failures.push(format!(
            "only {exact_count} cases agree exactly; expected at least 9"
        ));
    }
    assert!(
        failures.is_empty(),
        "MUSCLE 5.1 black-box acceptance failures:\n{}",
        failures.join("\n")
    );
}
