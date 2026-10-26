//! Offline black-box conformance against pinned WFA2-lib 2.3.6 observations.

mod support;

use std::fmt::Write as _;

use align_core::{AlignmentOperation, Alphabet, Sequence};
use align_fast::{
    WavefrontAffine, WavefrontAffinePenalties, WavefrontScore, WavefrontTraceRun,
    WavefrontWorkspace,
};
use serde::Deserialize;
use support::wavefront::replay_affine;

const FIXTURE: &str = include_str!("fixtures/wfa2-2.3.6/wavefront.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    schema_version: u32,
    profile: String,
    generated_at_utc: String,
    reference: Reference,
    environment: Environment,
    inputs: Inputs,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Reference {
    project: String,
    repository: String,
    tag: String,
    tag_object: String,
    commit: String,
    license: String,
    binary: String,
    binary_sha256: String,
    public_interface: String,
    source_code_inspected: bool,
    algorithms: Vec<String>,
    span: String,
    memory_mode: String,
    heuristic: Option<String>,
    threads: usize,
}

#[derive(Debug, Deserialize)]
struct Environment {
    host_platform: String,
    wsl_distribution: String,
    linux: String,
    compiler: String,
    build_tool: String,
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
    mode: String,
    alphabet: String,
    first: String,
    second: String,
    mismatch: Option<i64>,
    gap_opening: Option<i64>,
    gap_extension: Option<i64>,
    compare_cigar: bool,
    expected: Expected,
}

#[derive(Debug, Deserialize)]
struct Expected {
    raw_reference_score: i64,
    penalty: i64,
    cigar: String,
}

fn parse_reference_cigar(cigar: &str) -> Vec<AlignmentOperation> {
    let bytes = cigar.as_bytes();
    let mut operations = Vec::new();
    let mut index = 0_usize;
    while index < bytes.len() {
        let count_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        assert!(index > count_start, "CIGAR operation must have a count");
        let count = cigar[count_start..index]
            .parse::<usize>()
            .expect("CIGAR count is representable");
        assert!(count > 0, "CIGAR count must be positive");
        let operation = match bytes.get(index).copied() {
            Some(b'M') => AlignmentOperation::Match,
            Some(b'X') => AlignmentOperation::Substitution,
            Some(b'D') => AlignmentOperation::Deletion,
            Some(b'I') => AlignmentOperation::Insertion,
            value => panic!("unsupported reference CIGAR operation {value:?}"),
        };
        operations.extend(std::iter::repeat_n(operation, count));
        index += 1;
    }
    operations
}

fn encode_cigar(operations: &[AlignmentOperation]) -> String {
    let mut cigar = String::new();
    let mut start = 0_usize;
    while start < operations.len() {
        let operation = operations[start];
        let mut end = start + 1;
        while end < operations.len() && operations[end] == operation {
            end += 1;
        }
        let code = match operation {
            AlignmentOperation::Match => 'M',
            AlignmentOperation::Substitution => 'X',
            AlignmentOperation::Deletion => 'D',
            AlignmentOperation::Insertion => 'I',
        };
        write!(&mut cigar, "{}{}", end - start, code).expect("writing to String succeeds");
        start = end;
    }
    cigar
}

fn encode_trace_runs(runs: impl Iterator<Item = WavefrontTraceRun>) -> String {
    let mut cigar = String::new();
    for run in runs {
        let code = match run.operation() {
            AlignmentOperation::Match => 'M',
            AlignmentOperation::Substitution => 'X',
            AlignmentOperation::Deletion => 'D',
            AlignmentOperation::Insertion => 'I',
        };
        write!(&mut cigar, "{}{code}", run.count()).expect("writing to String succeeds");
    }
    cigar
}

fn assert_same_work(score: WavefrontScore, alignment: WavefrontScore, id: &str) {
    assert_eq!(score.penalty(), alignment.penalty(), "case {id} penalty");
    assert_eq!(
        score.wavefronts(),
        alignment.wavefronts(),
        "case {id} wavefronts"
    );
    assert_eq!(
        score.computed_offsets(),
        alignment.computed_offsets(),
        "case {id} offsets"
    );
    assert_eq!(
        score.source_wavefront_lookups(),
        alignment.source_wavefront_lookups(),
        "case {id} source wavefront lookups"
    );
    assert_eq!(
        score.source_offset_lookups(),
        alignment.source_offset_lookups(),
        "case {id} source offset lookups"
    );
    assert_eq!(
        score.match_extensions(),
        alignment.match_extensions(),
        "case {id} match extensions"
    );
    assert_eq!(
        score.maximum_wavefront_width(),
        alignment.maximum_wavefront_width(),
        "case {id} width"
    );
}

fn check_case<A: Alphabet>(case: &Case) {
    let first = Sequence::<A>::from_ascii(case.first.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid first sequence: {error}", case.id));
    let second = Sequence::<A>::from_ascii(case.second.as_bytes())
        .unwrap_or_else(|error| panic!("case {} has invalid second sequence: {error}", case.id));
    let reference_operations = parse_reference_cigar(&case.expected.cigar);

    let rust_operations = match case.mode.as_str() {
        "edit" => unreachable!("edit-WFA cases are filtered from the retained affine profile"),
        "affine" => {
            assert_eq!(-case.expected.raw_reference_score, case.expected.penalty);
            let mismatch = case.mismatch.expect("affine mismatch is present");
            let opening = case.gap_opening.expect("affine opening is present");
            let extension = case.gap_extension.expect("affine extension is present");
            assert_eq!(
                replay_affine(
                    case.first.as_bytes(),
                    case.second.as_bytes(),
                    &reference_operations,
                    mismatch,
                    opening,
                    extension,
                ),
                case.expected.penalty,
                "case {} external CIGAR replay",
                case.id
            );
            let penalties = WavefrontAffinePenalties::new(mismatch, opening, extension)
                .expect("fixture affine penalties are valid");
            let aligner = WavefrontAffine::new(penalties);
            let score = aligner
                .score(&first, &second)
                .unwrap_or_else(|error| panic!("case {} score failed: {error}", case.id));
            let result = aligner
                .align(&first, &second)
                .unwrap_or_else(|error| panic!("case {} alignment failed: {error}", case.id));
            let mut workspace = WavefrontWorkspace::new();
            let compact = aligner
                .traceback_in_workspace(&first, &second, &mut workspace)
                .unwrap_or_else(|error| {
                    panic!("case {} compact traceback failed: {error}", case.id)
                });
            let compact_penalty = compact.penalty();
            let compact_columns = compact.columns();
            let compact_cigar = encode_trace_runs(compact.runs());
            assert_eq!(score.penalty(), case.expected.penalty, "case {}", case.id);
            assert_eq!(result.penalty(), case.expected.penalty, "case {}", case.id);
            assert_eq!(compact_penalty, case.expected.penalty, "case {}", case.id);
            assert_same_work(score, result.summary(), &case.id);
            assert_eq!(
                compact_columns,
                result.alignment().columns(),
                "case {} compact columns",
                case.id
            );
            assert_eq!(
                compact_cigar,
                encode_cigar(result.alignment().operations()),
                "case {} compact CIGAR",
                case.id
            );
            assert_eq!(
                replay_affine(
                    case.first.as_bytes(),
                    case.second.as_bytes(),
                    result.alignment().operations(),
                    mismatch,
                    opening,
                    extension,
                ),
                case.expected.penalty,
                "case {} Rust CIGAR replay",
                case.id
            );
            result.alignment().operations().to_vec()
        }
        mode => panic!("case {} has unsupported mode {mode}", case.id),
    };

    if case.compare_cigar {
        assert_eq!(
            encode_cigar(&rust_operations),
            case.expected.cigar,
            "case {} unambiguous CIGAR",
            case.id
        );
    }
}

#[test]
fn pinned_wfa2_affine_scores_and_unambiguous_cigars_are_reproduced() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("fixture JSON is valid");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(fixture.profile, "WFA-global-affine-v1");
    assert!(fixture.generated_at_utc.ends_with("+00:00"));
    assert_eq!(fixture.reference.project, "WFA2-lib");
    assert_eq!(
        fixture.reference.repository,
        "https://github.com/smarco/WFA2-lib"
    );
    assert_eq!(fixture.reference.tag, "v2.3.6");
    assert_eq!(
        fixture.reference.tag_object,
        "0db345a8fe862fd7873d3354c499da385583a65a"
    );
    assert_eq!(
        fixture.reference.commit,
        "bcf473a6561fa297934a80eaa7c04b4ee412360c"
    );
    assert_eq!(fixture.reference.license, "MIT");
    assert_eq!(fixture.reference.binary, "bin/align_benchmark");
    assert_eq!(fixture.reference.binary_sha256.len(), 64);
    assert_eq!(
        fixture.reference.public_interface,
        "tools/README.md align-benchmark CLI"
    );
    assert!(!fixture.reference.source_code_inspected);
    assert_eq!(fixture.reference.algorithms, ["gap-affine-wfa"]);
    assert_eq!(fixture.reference.span, "global");
    assert_eq!(fixture.reference.memory_mode, "high");
    assert!(fixture.reference.heuristic.is_none());
    assert_eq!(fixture.reference.threads, 1);
    assert!(fixture.environment.host_platform.starts_with("Windows-11"));
    assert_eq!(fixture.environment.wsl_distribution, "Ubuntu-24.04");
    assert!(fixture.environment.linux.contains("Linux"));
    assert!(fixture.environment.compiler.starts_with("gcc (Ubuntu"));
    assert_eq!(fixture.environment.build_tool, "GNU Make 4.3");
    assert_eq!(fixture.inputs.cases_file, "wfa.json");
    assert_eq!(fixture.inputs.cases_sha256.len(), 64);
    assert_eq!(fixture.inputs.generator_sha256.len(), 64);
    assert_eq!(fixture.cases.len(), 19);

    for case in &fixture.cases {
        assert_eq!(case.mode, "affine");
        match case.alphabet.as_str() {
            "dna" => check_case::<align_core::Dna>(case),
            "rna" => check_case::<align_core::Rna>(case),
            "protein" => check_case::<align_core::Protein>(case),
            alphabet => panic!("case {} has unsupported alphabet {alphabet}", case.id),
        }
    }
}
