//! Offline black-box X-drop conformance against pinned Biotite 1.7.1 fixtures.

mod support;

use align_core::{AffineGap, Alphabet, MatchMismatch, Sequence};
use align_fast::{XDrop, XDropDirection, XDropSeed};
use serde::Deserialize;
use support::xdrop::replay;

const FIXTURE: &str = include_str!("fixtures/biotite-1.7.1/xdrop.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    schema_version: u32,
    profile: String,
    reference: Reference,
    inputs: Inputs,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Reference {
    name: String,
    package: String,
    version: String,
    license: String,
    requires_python: String,
    public_api_documentation: String,
    scientific_method_doi: String,
    source_code_inspected: bool,
    wheel: Wheel,
}

#[derive(Debug, Deserialize)]
struct Wheel {
    filename: String,
    sha256: String,
}

#[derive(Debug, Deserialize)]
struct Inputs {
    cases_file: String,
    cases_sha256: String,
    generator_sha256: String,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    alphabet: String,
    first: String,
    second: String,
    seed: [usize; 2],
    direction: String,
    threshold: i64,
    match_score: i64,
    mismatch_score: i64,
    gap_open: i64,
    gap_extend: i64,
    expected: Expected,
}

#[derive(Debug, Deserialize)]
struct Expected {
    score: i64,
    first_range: [usize; 2],
    second_range: [usize; 2],
    aligned_first: String,
    aligned_second: String,
    trace: Vec<[isize; 2]>,
}

fn direction(value: &str) -> XDropDirection {
    match value {
        "upstream" => XDropDirection::Upstream,
        "downstream" => XDropDirection::Downstream,
        "both" => XDropDirection::Both,
        value => panic!("unknown fixture X-drop direction {value}"),
    }
}

fn replay_reference(case: &Case) {
    let first_seed = isize::try_from(case.seed[0]).expect("fixture seed fits isize");
    let second_seed = isize::try_from(case.seed[1]).expect("fixture seed fits isize");
    assert_eq!(
        case.expected.aligned_first.len(),
        case.expected.aligned_second.len(),
        "case {}",
        case.id
    );
    let first_ungapped = case.expected.aligned_first.replace('-', "");
    let second_ungapped = case.expected.aligned_second.replace('-', "");
    assert_eq!(
        first_ungapped,
        case.first[case.expected.first_range[0]..case.expected.first_range[1]],
        "case {} first range",
        case.id
    );
    assert_eq!(
        second_ungapped,
        case.second[case.expected.second_range[0]..case.expected.second_range[1]],
        "case {} second range",
        case.id
    );
    assert_eq!(
        case.expected
            .trace
            .iter()
            .filter(|coordinate| { coordinate[0] == first_seed && coordinate[1] == second_seed })
            .count(),
        1,
        "case {} forced seed count",
        case.id
    );

    let mut score = 0_i64;
    let mut prior_gap = None;
    for (first, second) in case
        .expected
        .aligned_first
        .bytes()
        .zip(case.expected.aligned_second.bytes())
    {
        let increment = match (first, second) {
            (b'-', b'-') => panic!("case {} contains a double-gap column", case.id),
            (b'-', _) => {
                let increment = if prior_gap == Some(b'I') {
                    case.gap_extend
                } else {
                    case.gap_open
                };
                prior_gap = Some(b'I');
                increment
            }
            (_, b'-') => {
                let increment = if prior_gap == Some(b'D') {
                    case.gap_extend
                } else {
                    case.gap_open
                };
                prior_gap = Some(b'D');
                increment
            }
            _ => {
                prior_gap = None;
                if first == second {
                    case.match_score
                } else {
                    case.mismatch_score
                }
            }
        };
        score += increment;
    }
    assert_eq!(
        score, case.expected.score,
        "case {} fixture replay",
        case.id
    );
}

fn check_case<A: Alphabet>(case: &Case) -> Option<String> {
    let first = Sequence::<A>::from_ascii(case.first.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid first sequence: {error}", case.id));
    let second = Sequence::<A>::from_ascii(case.second.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid second sequence: {error}", case.id));
    let scoring = MatchMismatch::new(case.match_score, case.mismatch_score);
    let gap = AffineGap::new(case.gap_open, case.gap_extend)
        .unwrap_or_else(|error| panic!("case {} has invalid gap scores: {error}", case.id));
    let aligner = XDrop::new(scoring, gap, case.threshold)
        .unwrap_or_else(|error| panic!("case {} has invalid X-drop settings: {error}", case.id));
    let seed = XDropSeed::new(case.seed[0], case.seed[1]);
    let direction = direction(&case.direction);
    let score = aligner
        .score_from_seed(&first, &second, seed, direction)
        .unwrap_or_else(|error| panic!("case {} score failed: {error}", case.id));
    assert_eq!(score.score(), case.expected.score, "case {} score", case.id);

    let result = aligner
        .align_from_seed(&first, &second, seed, direction)
        .unwrap_or_else(|error| panic!("case {} traceback failed: {error}", case.id));
    assert_eq!(result.summary(), score, "case {} summary", case.id);
    assert_eq!(
        replay(&first, &second, result.alignment(), &scoring, gap),
        case.expected.score,
        "case {} Rust replay",
        case.id
    );

    let observed_first = result.first_range();
    let observed_second = result.second_range();
    let expected_first = case.expected.first_range[0]..case.expected.first_range[1];
    let expected_second = case.expected.second_range[0]..case.expected.second_range[1];
    (observed_first != expected_first || observed_second != expected_second).then(|| {
        format!(
            "{}: Rust {:?}/{:?}, Biotite {:?}/{:?}",
            case.id, observed_first, observed_second, expected_first, expected_second
        )
    })
}

#[test]
fn pinned_biotite_scores_and_unambiguous_ranges_are_reproduced() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("fixture JSON is valid");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.profile, "XDROP-forward-affine-v1");
    assert_eq!(fixture.reference.name, "Biotite align_local_gapped");
    assert_eq!(fixture.reference.package, "biotite");
    assert_eq!(fixture.reference.version, "1.7.1");
    assert_eq!(fixture.reference.license, "BSD-3-Clause");
    assert_eq!(fixture.reference.requires_python, ">=3.12");
    assert!(
        fixture
            .reference
            .public_api_documentation
            .ends_with("align.align_local_gapped.html")
    );
    assert_eq!(
        fixture.reference.scientific_method_doi,
        "10.1089/10665270050081478"
    );
    assert!(!fixture.reference.source_code_inspected);
    assert!(
        fixture
            .reference
            .wheel
            .filename
            .starts_with("biotite-1.7.1-")
    );
    assert_eq!(
        fixture.reference.wheel.sha256,
        "965824257f716a77de3df35616cdaf302483a9ab8df98dc1dd0096c987342c49"
    );
    assert_eq!(fixture.inputs.cases_file, "xdrop.json");
    assert_eq!(fixture.inputs.cases_sha256.len(), 64);
    assert_eq!(fixture.inputs.generator_sha256.len(), 64);
    assert_eq!(fixture.cases.len(), 24);

    let mut range_differences = Vec::new();
    for case in &fixture.cases {
        replay_reference(case);
        let difference = match case.alphabet.as_str() {
            "dna" => check_case::<align_core::Dna>(case),
            "protein" => check_case::<align_core::Protein>(case),
            alphabet => panic!("case {} has unknown alphabet {alphabet}", case.id),
        };
        if let Some(difference) = difference {
            range_differences.push(difference);
        }
    }
    assert_eq!(
        range_differences,
        ["dna_long_affine_gap: Rust 0..7/0..4, Biotite 0..3/0..3".to_owned()],
        "only the documented equal-score endpoint tie may differ"
    );
}
