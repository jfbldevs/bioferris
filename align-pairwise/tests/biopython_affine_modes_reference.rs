//! Static black-box conformance for affine local and four-end modes against Biopython 1.87.

use std::collections::BTreeSet;

use align_core::{AffineGap, Alphabet, Dna, MatchMismatch, Protein, Rna, Sequence};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, FreeEndGaps, verify_affine_end_gap_alignment,
    verify_local_affine_alignment,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/biopython-1.87/affine-modes.json");

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
    global_algorithm: String,
    local_algorithm: String,
    global_mode: String,
    local_mode: String,
    gap_convention: String,
    boundary_mapping: BoundaryMapping,
    epsilon: f64,
}

#[derive(Debug, Deserialize)]
struct BoundaryMapping {
    first_start: String,
    first_end: String,
    second_start: String,
    second_end: String,
}

#[derive(Debug, Deserialize)]
struct ReferenceCase {
    id: String,
    profile: String,
    alphabet: String,
    first: String,
    second: String,
    scoring: ReferenceScoring,
    free_ends: Option<ReferenceFreeEnds>,
    comparison: String,
    expected: ReferenceExpected,
}

#[derive(Clone, Copy, Debug, Deserialize)]
struct ReferenceScoring {
    #[serde(rename = "match")]
    identical: i64,
    #[serde(rename = "mismatch")]
    different: i64,
    open: i64,
    extend: i64,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd)]
#[allow(clippy::struct_excessive_bools)]
struct ReferenceFreeEnds {
    first_start: bool,
    first_end: bool,
    second_start: bool,
    second_end: bool,
}

impl ReferenceFreeEnds {
    const fn configuration(self) -> FreeEndGaps {
        FreeEndGaps::new(
            self.first_start,
            self.first_end,
            self.second_start,
            self.second_end,
        )
    }
}

#[derive(Debug, Deserialize)]
struct ReferenceExpected {
    algorithm: String,
    score: i64,
    reported_optimal_alignment_count: Option<usize>,
    profile_compatible_reference_alignment_count: Option<usize>,
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

fn assert_expected_set(case: &ReferenceCase, observed: &ReferenceAlignment) {
    if case.expected.alignments_complete {
        assert_eq!(case.comparison, "all_optimal", "case {}", case.id);
        let reported = case
            .expected
            .reported_optimal_alignment_count
            .expect("complete fixture records raw reference count");
        let compatible = case
            .expected
            .profile_compatible_reference_alignment_count
            .expect("complete fixture records compatible count");
        let unique = case
            .expected
            .unique_alignment_count
            .expect("complete fixture records unique count");
        assert!(reported > 0, "case {}", case.id);
        assert!(compatible >= unique, "case {}", case.id);
        assert_eq!(unique, case.expected.alignments.len(), "case {}", case.id);
        assert!(
            case.expected.alignments.contains(observed),
            "case {} produced the reference score, but its ranges, traceback, or metrics were not \
             in Biopython's complete compatible optimum set: {observed:?}",
            case.id,
        );
    } else {
        assert_eq!(case.comparison, "score_only", "case {}", case.id);
        assert!(
            case.expected.reported_optimal_alignment_count.is_none(),
            "case {}",
            case.id
        );
        assert!(
            case.expected
                .profile_compatible_reference_alignment_count
                .is_none(),
            "case {}",
            case.id
        );
        assert!(
            case.expected.unique_alignment_count.is_none(),
            "case {}",
            case.id
        );
        assert!(case.expected.alignments.is_empty(), "case {}", case.id);
    }
}

fn compare_case<A: Alphabet>(case: &ReferenceCase, provenance: &Provenance) {
    let first = Sequence::<A>::from_ascii(case.first.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid first sequence: {error}", case.id));
    let second = Sequence::<A>::from_ascii(case.second.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid second sequence: {error}", case.id));
    let scorer = MatchMismatch::new(case.scoring.identical, case.scoring.different);
    let gap = AffineGap::new(case.scoring.open, case.scoring.extend)
        .unwrap_or_else(|error| panic!("case {} has invalid affine scores: {error}", case.id));

    let (alignment, report) = match case.profile.as_str() {
        "end_gap" => {
            assert_eq!(case.expected.algorithm, provenance.global_algorithm);
            let free_ends = case
                .free_ends
                .expect("end-gap fixture must define four end permissions")
                .configuration();
            let alignment = AffineEndGapAligner::new(scorer, gap, free_ends)
                .align(&first, &second)
                .unwrap_or_else(|error| panic!("case {} failed to align: {error}", case.id));
            let report = verify_affine_end_gap_alignment(
                &first, &second, &alignment, &scorer, gap, free_ends,
            )
            .unwrap_or_else(|error| panic!("case {} failed independent replay: {error}", case.id));
            (alignment, report)
        }
        "local" => {
            assert_eq!(case.expected.algorithm, provenance.local_algorithm);
            assert!(case.free_ends.is_none(), "local case {}", case.id);
            let alignment = AffineLocalAligner::new(scorer, gap)
                .align(&first, &second)
                .unwrap_or_else(|error| panic!("case {} failed to align: {error}", case.id));
            let report = verify_local_affine_alignment(&first, &second, &alignment, &scorer, gap)
                .unwrap_or_else(|error| {
                    panic!("case {} failed independent replay: {error}", case.id)
                });
            (alignment, report)
        }
        profile => panic!("case {} has unknown profile {profile:?}", case.id),
    };

    assert_eq!(
        alignment.score(),
        case.expected.score,
        "case {} disagrees with the Biopython optimal score",
        case.id,
    );
    let rendered = alignment
        .render(&first, &second)
        .unwrap_or_else(|error| panic!("case {} failed to render: {error}", case.id));
    let first_range = alignment.first_range();
    let second_range = alignment.second_range();
    assert_expected_set(
        case,
        &ReferenceAlignment {
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
        },
    );
}

#[test]
fn affine_modes_match_biopython_1_87_black_box_fixtures() {
    let fixture: ReferenceFixture =
        serde_json::from_str(FIXTURE).expect("the committed reference fixture is valid JSON");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.suite, "AFFINE-MODES-v1 Biopython conformance");
    assert_eq!(fixture.provenance.package, "biopython");
    assert_eq!(fixture.provenance.package_version, "1.87");
    assert_eq!(fixture.provenance.numpy_version, "2.5.1");
    assert_eq!(fixture.provenance.api, "Bio.Align.PairwiseAligner");
    assert_eq!(
        fixture.provenance.global_algorithm,
        "Gotoh global alignment algorithm"
    );
    assert_eq!(
        fixture.provenance.local_algorithm,
        "Gotoh local alignment algorithm"
    );
    assert_eq!(fixture.provenance.global_mode, "global");
    assert_eq!(fixture.provenance.local_mode, "local");
    assert_eq!(
        fixture.provenance.gap_convention,
        "open + (length - 1) * extend"
    );
    assert_eq!(
        fixture.provenance.boundary_mapping.first_start,
        "left deletion"
    );
    assert_eq!(
        fixture.provenance.boundary_mapping.first_end,
        "right deletion"
    );
    assert_eq!(
        fixture.provenance.boundary_mapping.second_start,
        "left insertion"
    );
    assert_eq!(
        fixture.provenance.boundary_mapping.second_end,
        "right insertion"
    );
    assert!((fixture.provenance.epsilon - 1e-9).abs() <= f64::EPSILON);
    assert_eq!(
        fixture.provenance.generator_sha256,
        "dfe257d5e201a9f629e4de91e19d7dd161e1857b9f3d84594b9e628ad8ec8a96"
    );
    assert_eq!(
        fixture.provenance.case_specification_sha256,
        "cbda5eec88fa9364b941957814d933a7181f14ee626ea381438d4be8475cb9c8"
    );
    assert_eq!(fixture.provenance.python_version, "3.13.14");
    assert!(fixture.provenance.platform.starts_with("Windows-11-"));
    assert_eq!(fixture.cases.len(), 27);
    assert_eq!(
        fixture
            .cases
            .iter()
            .filter(|case| case.expected.alignments_complete)
            .count(),
        25
    );

    let configurations = fixture
        .cases
        .iter()
        .filter(|case| case.id.starts_with("endgap_mode_"))
        .map(|case| case.free_ends.expect("mode fixture defines ends"))
        .collect::<BTreeSet<_>>();
    assert_eq!(configurations.len(), 16);

    for case in &fixture.cases {
        match case.alphabet.as_str() {
            "dna" => compare_case::<Dna>(case, &fixture.provenance),
            "rna" => compare_case::<Rna>(case, &fixture.provenance),
            "protein" => compare_case::<Protein>(case, &fixture.provenance),
            alphabet => panic!("case {} uses unknown alphabet {alphabet:?}", case.id),
        }
    }
}
