//! Black-box local-alignment conformance against Biopython 1.87 static fixtures.

use align_core::{Alphabet, Dna, LinearGap, MatchMismatch, Protein, Rna, Sequence};
use align_pairwise::{SmithWaterman, verify_local_alignment};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/biopython-1.87/smith-waterman.json");

#[derive(Debug, Deserialize)]
struct ReferenceFixture {
    schema_version: u32,
    suite: String,
    provenance: Provenance,
    cases: Vec<ReferenceCase>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    generator_sha256: String,
    case_specification_sha256: String,
    python_version: String,
    platform: String,
    package: String,
    package_version: String,
    numpy_version: String,
    api: String,
    required_algorithm: String,
    mode: String,
    epsilon: f64,
}

#[derive(Debug, Deserialize)]
struct ReferenceCase {
    id: String,
    alphabet: String,
    first: String,
    second: String,
    scoring: ReferenceScoring,
    comparison: String,
    expected: ReferenceExpected,
}

#[derive(Clone, Copy, Debug, Deserialize)]
struct ReferenceScoring {
    #[serde(rename = "match")]
    match_score: i64,
    #[serde(rename = "mismatch")]
    mismatch_score: i64,
    gap: i64,
}

#[derive(Debug, Deserialize)]
struct ReferenceExpected {
    algorithm: String,
    score: i64,
    reported_optimal_alignment_count: Option<usize>,
    unique_alignment_count: Option<usize>,
    alignments_complete: bool,
    alignments: Vec<ReferenceAlignment>,
}

#[derive(Debug, Deserialize, Eq, PartialEq)]
struct ReferenceAlignment {
    first_start: usize,
    first_end: usize,
    second_start: usize,
    second_end: usize,
    first_aligned: String,
    second_aligned: String,
    columns: usize,
    matches: usize,
    substitutions: usize,
    deletions: usize,
    insertions: usize,
    score: i64,
}

fn compare_case<A: Alphabet>(case: &ReferenceCase) {
    let first = Sequence::<A>::from_ascii(case.first.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid first sequence: {error}", case.id));
    let second = Sequence::<A>::from_ascii(case.second.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid second sequence: {error}", case.id));
    let scorer = MatchMismatch::new(case.scoring.match_score, case.scoring.mismatch_score);
    let gap = LinearGap::new(case.scoring.gap)
        .unwrap_or_else(|error| panic!("case {} has invalid gap score: {error}", case.id));
    let aligner = SmithWaterman::new(scorer, gap);
    let alignment = aligner
        .align(&first, &second)
        .unwrap_or_else(|error| panic!("case {} failed to align: {error}", case.id));
    let report = verify_local_alignment(&first, &second, &alignment, &scorer, gap)
        .unwrap_or_else(|error| panic!("case {} failed independent replay: {error}", case.id));

    assert_eq!(
        alignment.score(),
        case.expected.score,
        "case {} disagrees with the Biopython optimal score",
        case.id
    );

    if case.expected.alignments_complete {
        assert_eq!(case.comparison, "all_optimal", "case {}", case.id);
        let reported_count = case
            .expected
            .reported_optimal_alignment_count
            .expect("complete fixture records the reference alignment count");
        let unique_count = case
            .expected
            .unique_alignment_count
            .expect("complete fixture records the unique alignment count");
        assert!(reported_count >= unique_count, "case {}", case.id);
        assert_eq!(
            unique_count,
            case.expected.alignments.len(),
            "case {} has inconsistent fixture counts",
            case.id
        );

        let rendered = alignment
            .render(&first, &second)
            .unwrap_or_else(|error| panic!("case {} failed to render: {error}", case.id));
        let first_range = alignment.first_range();
        let second_range = alignment.second_range();
        let observed = ReferenceAlignment {
            first_start: first_range.start,
            first_end: first_range.end,
            second_start: second_range.start,
            second_end: second_range.end,
            first_aligned: rendered.first().to_owned(),
            second_aligned: rendered.second().to_owned(),
            columns: report.columns(),
            matches: report.matches(),
            substitutions: report.substitutions(),
            deletions: report.deletions(),
            insertions: report.insertions(),
            score: report.recomputed_score(),
        };
        assert!(
            case.expected.alignments.contains(&observed),
            "case {} produced an optimal score but its ranges, traceback, and metrics were not \
             among Biopython's complete optimal set: {observed:?}",
            case.id
        );
    } else {
        assert_eq!(case.comparison, "score_only", "case {}", case.id);
        assert!(
            case.expected.reported_optimal_alignment_count.is_none(),
            "case {} must not claim an unenumerated alignment count",
            case.id
        );
        assert!(
            case.expected.unique_alignment_count.is_none(),
            "case {} must not claim an unenumerated unique count",
            case.id
        );
        assert!(case.expected.alignments.is_empty(), "case {}", case.id);
    }
}

#[test]
fn smith_waterman_matches_biopython_1_87_black_box_fixtures() {
    let fixture: ReferenceFixture =
        serde_json::from_str(FIXTURE).expect("the committed reference fixture is valid JSON");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.suite, "SW-local-linear-v1 Biopython conformance");
    assert_eq!(fixture.provenance.package, "biopython");
    assert_eq!(fixture.provenance.package_version, "1.87");
    assert_eq!(fixture.provenance.numpy_version, "2.5.1");
    assert_eq!(fixture.provenance.api, "Bio.Align.PairwiseAligner");
    assert_eq!(fixture.provenance.required_algorithm, "Smith-Waterman");
    assert_eq!(fixture.provenance.mode, "local");
    assert!((fixture.provenance.epsilon - 1e-9).abs() <= f64::EPSILON);
    assert_eq!(fixture.provenance.generator_sha256.len(), 64);
    assert_eq!(fixture.provenance.case_specification_sha256.len(), 64);
    assert_eq!(fixture.provenance.python_version, "3.13.14");
    assert!(fixture.provenance.platform.starts_with("Windows-11-"));
    assert_eq!(fixture.cases.len(), 11);
    assert_eq!(
        fixture
            .cases
            .iter()
            .filter(|case| case.expected.alignments_complete)
            .count(),
        8
    );

    for case in &fixture.cases {
        assert_eq!(
            case.expected.algorithm,
            fixture.provenance.required_algorithm
        );
        match case.alphabet.as_str() {
            "dna" => compare_case::<Dna>(case),
            "rna" => compare_case::<Rna>(case),
            "protein" => compare_case::<Protein>(case),
            alphabet => panic!("case {} uses unknown alphabet {alphabet:?}", case.id),
        }
    }
}
