//! Offline black-box conformance against Edlib 1.3.9.post1 fixtures.

use align_core::{Alphabet, DnaSequence, ProteinSequence, RnaSequence, Sequence};
use align_fast::{MyersPattern, verify_unit_edit_alignment};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/edlib-1.3.9.post1/myers.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    schema_version: u32,
    suite: String,
    generated_at_utc: String,
    provenance: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    generator_sha256: String,
    case_specification_sha256: String,
    python_version: String,
    platform: String,
    package: String,
    package_version: String,
    api: String,
    global_mode: String,
    global_task: String,
    search_mode: String,
    search_task: String,
    literal_equality: bool,
    distance_limit: i32,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    alphabet: String,
    pattern: String,
    text: String,
    search_threshold: usize,
    expected: Expected,
}

#[derive(Debug, Deserialize)]
struct Expected {
    global: GlobalExpected,
    search_scores: Vec<usize>,
    search_hits: Vec<SearchHit>,
}

#[derive(Debug, Deserialize)]
struct GlobalExpected {
    edit_distance: usize,
    cigar: Option<String>,
    pattern_aligned: String,
    text_aligned: String,
    columns: usize,
    matches: usize,
    substitutions: usize,
    deletions: usize,
    insertions: usize,
    edits: usize,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
struct SearchHit {
    end: usize,
    distance: usize,
}

fn check_case<A: Alphabet>(case: &Case) {
    let pattern = Sequence::<A>::from_ascii(case.pattern.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid pattern: {error}", case.id));
    let text = Sequence::<A>::from_ascii(case.text.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid text: {error}", case.id));
    let myers = MyersPattern::new(&pattern)
        .unwrap_or_else(|error| panic!("case {} preprocessing failed: {error}", case.id));

    let distance = myers
        .global_distance(&text)
        .unwrap_or_else(|error| panic!("case {} distance failed: {error}", case.id));
    assert_eq!(
        distance, case.expected.global.edit_distance,
        "case {}",
        case.id
    );

    let alignment = myers
        .align_global(&text)
        .unwrap_or_else(|error| panic!("case {} traceback failed: {error}", case.id));
    let report = verify_unit_edit_alignment(&pattern, &text, &alignment)
        .unwrap_or_else(|error| panic!("case {} replay failed: {error}", case.id));
    assert_eq!(
        report.distance(),
        case.expected.global.edit_distance,
        "case {}",
        case.id
    );

    let all_scores: Vec<usize> = myers
        .find_all(&text, usize::MAX)
        .unwrap_or_else(|error| panic!("case {} complete search failed: {error}", case.id))
        .into_iter()
        .map(align_fast::MyersMatch::distance)
        .collect();
    assert_eq!(all_scores, case.expected.search_scores, "case {}", case.id);

    let hits: Vec<SearchHit> = myers
        .find_all(&text, case.search_threshold)
        .unwrap_or_else(|error| panic!("case {} threshold search failed: {error}", case.id))
        .into_iter()
        .map(|hit| SearchHit {
            end: hit.end(),
            distance: hit.distance(),
        })
        .collect();
    assert_eq!(hits, case.expected.search_hits, "case {}", case.id);
}

#[test]
fn pinned_edlib_distances_paths_and_search_scores_are_reproduced() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("fixture JSON is valid");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.suite, "MYERS-bit-vector-v1 black-box conformance");
    assert!(fixture.generated_at_utc.ends_with("+00:00"));
    assert_eq!(fixture.provenance.generator_sha256.len(), 64);
    assert_eq!(fixture.provenance.case_specification_sha256.len(), 64);
    assert_eq!(fixture.provenance.python_version, "3.13.14");
    assert!(fixture.provenance.platform.starts_with("Windows-11"));
    assert_eq!(fixture.provenance.package, "edlib");
    assert_eq!(fixture.provenance.package_version, "1.3.9.post1");
    assert_eq!(fixture.provenance.api, "edlib.align");
    assert_eq!(fixture.provenance.global_mode, "NW");
    assert_eq!(fixture.provenance.global_task, "path");
    assert_eq!(fixture.provenance.search_mode, "reversed-prefix SHW");
    assert_eq!(fixture.provenance.search_task, "distance");
    assert!(fixture.provenance.literal_equality);
    assert_eq!(fixture.provenance.distance_limit, -1);
    assert_eq!(fixture.cases.len(), 20);

    for case in &fixture.cases {
        let global = &case.expected.global;
        assert_eq!(global.pattern_aligned.len(), global.text_aligned.len());
        assert_eq!(global.columns, global.pattern_aligned.len());
        assert_eq!(global.pattern_aligned.replace('-', ""), case.pattern);
        assert_eq!(global.text_aligned.replace('-', ""), case.text);
        assert_eq!(
            global.matches + global.substitutions + global.deletions + global.insertions,
            global.columns,
            "case {} fixture column metrics",
            case.id
        );
        assert_eq!(
            global.substitutions + global.deletions + global.insertions,
            global.edits,
            "case {} fixture edit metrics",
            case.id
        );
        assert_eq!(global.edits, global.edit_distance, "case {}", case.id);
        if let Some(cigar) = &global.cigar {
            assert!(!cigar.is_empty(), "case {} has an empty CIGAR", case.id);
        }
        assert_eq!(case.expected.search_scores.len(), case.text.len());

        match case.alphabet.as_str() {
            "dna" => {
                DnaSequence::from_ascii(case.pattern.as_bytes()).expect("fixture DNA is valid");
                check_case::<align_core::Dna>(case);
            }
            "rna" => {
                RnaSequence::from_ascii(case.pattern.as_bytes()).expect("fixture RNA is valid");
                check_case::<align_core::Rna>(case);
            }
            "protein" => {
                ProteinSequence::from_ascii(case.pattern.as_bytes())
                    .expect("fixture protein is valid");
                check_case::<align_core::Protein>(case);
            }
            alphabet => panic!("case {} has unknown alphabet {alphabet}", case.id),
        }
    }
}
