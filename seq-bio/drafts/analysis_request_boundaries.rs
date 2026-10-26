//! Actual-binary regressions for isolation between dedicated request namespaces.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

use serde_json::{Value, json};

struct FeatureCase {
    field: &'static str,
    request: Value,
    populated: Value,
    duplicate_key: &'static str,
}

fn missing_input() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/analysis-request-boundary-must-not-be-opened.missing");
    assert!(
        !path.exists(),
        "the rejection probe requires an absent input"
    );
    path.to_str()
        .expect("UTF-8 workspace fixture path")
        .to_owned()
}

fn feature_cases() -> [FeatureCase; 4] {
    [
        FeatureCase {
            field: "motif_parameters",
            request: json!({"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"A"}),
            populated: json!({"strands":"forward"}),
            duplicate_key: "max_matches",
        },
        FeatureCase {
            field: "analysis_parameters",
            request: json!({"schema_version":1,"operation":"alignment-stats","molecule":"dna","format":"aligned-fasta","input":missing_input()}),
            populated: json!({"max_rows":10}),
            duplicate_key: "max_rows",
        },
        FeatureCase {
            field: "blast_parameters",
            request: json!({"schema_version":1,"operation":"inspect-blast","format":"blast-tab6","input":missing_input()}),
            populated: json!({"max_queries":10}),
            duplicate_key: "max_queries",
        },
        FeatureCase {
            field: "compression_parameters",
            request: json!({"schema_version":1,"operation":"compress-file","compression":"gzip","input":missing_input(),"output":format!("{}.gz",missing_input())}),
            populated: json!({"level":1}),
            duplicate_key: "level",
        },
    ]
}

fn sequence_requests() -> [Value; 2] {
    [
        json!({"schema_version":1,"operation":"translate","molecule":"dna","sequence":"ATGGCTTAA","mode":"cds"}),
        json!({"schema_version":1,"operation":"protein-properties","molecule":"protein","sequence":"INGAR"}),
    ]
}

fn run(raw: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bioferris"))
        .args(["request", "--input", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("production binary starts");
    child.stdin.take().unwrap().write_all(raw).unwrap();
    child.wait_with_output().unwrap()
}

fn envelope(output: &Output, status: i32) -> Value {
    assert_eq!(output.status.code(), Some(status), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    assert_eq!(String::from_utf8_lossy(&output.stdout).lines().count(), 1);
    let value: Value = serde_json::from_slice(&output.stdout).expect("one complete JSON envelope");
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], status == 0);
    assert_eq!(value.get("result").is_some(), status == 0);
    assert_eq!(value.get("error").is_some(), status != 0);
    value
}

fn rejected(raw: &[u8]) -> Value {
    let first = run(raw);
    let value = envelope(&first, 2);
    let second = run(raw);
    envelope(&second, 2);
    assert_eq!(
        first.stdout, second.stdout,
        "error output must be deterministic"
    );
    assert_eq!(value["error"]["code"], "invalid_schema");
    value
}

fn rejected_field(request: &Value, field: &str) {
    let error = rejected(&serde_json::to_vec(request).unwrap());
    assert!(
        error["error"]["message"].as_str().unwrap().contains(field),
        "the error must identify the irrelevant field, not fail on input I/O: {error}"
    );
}

fn with_raw_field(request: &Value, field: &str, raw_value: &str) -> Vec<u8> {
    let serialized = request.to_string();
    let prefix = serialized.strip_suffix('}').unwrap();
    format!("{prefix},\"{field}\":{raw_value}}}").into_bytes()
}

#[test]
fn dedicated_objects_are_rejected_on_unrelated_existing_sequence_operations() {
    for base in sequence_requests() {
        for case in feature_cases() {
            for payload in [json!({}), case.populated] {
                let mut request = base.clone();
                request[case.field] = payload;
                rejected_field(&request, case.field);
            }
        }
        for payload in [json!({}), json!({"match_score":2})] {
            let mut request = base.clone();
            request["parameters"] = payload;
            rejected_field(&request, "parameters");
        }
    }
}

#[test]
fn dedicated_objects_do_not_cross_feature_boundaries_before_file_access() {
    for owner in feature_cases() {
        for foreign in feature_cases() {
            if owner.field == foreign.field {
                continue;
            }
            for payload in [json!({}), foreign.populated] {
                let mut request = owner.request.clone();
                request[foreign.field] = payload;
                rejected_field(&request, foreign.field);
            }
        }
        for (field, value) in [("ph", json!(7)), ("parameters", json!({}))] {
            let mut request = owner.request.clone();
            request[field] = value;
            rejected_field(&request, field);
        }
    }
}

#[test]
fn dedicated_nested_fields_require_strict_objects_and_preserve_duplicate_detection() {
    for case in feature_cases() {
        for wrong_type in ["[]", "[1]", "7", "false", "\"not an object\""] {
            let error = rejected(&with_raw_field(&case.request, case.field, wrong_type));
            assert!(
                error["error"]["message"]
                    .as_str()
                    .unwrap()
                    .contains("invalid type")
            );
        }
        let error = rejected(&with_raw_field(
            &case.request,
            case.field,
            r#"{"unexpected_nested_option":1}"#,
        ));
        assert!(
            error["error"]["message"]
                .as_str()
                .unwrap()
                .contains("unexpected_nested_option")
        );
        let duplicate = format!(
            "{{\"{}\":1,\"{}\":2}}",
            case.duplicate_key, case.duplicate_key
        );
        let error = rejected(&with_raw_field(&case.request, case.field, &duplicate));
        let message = error["error"]["message"].as_str().unwrap();
        assert!(message.contains("duplicate field") && message.contains(case.duplicate_key));
    }
}

#[test]
fn optional_null_fields_remain_absent_and_preserve_existing_result_bytes() {
    for base in sequence_requests() {
        let original = run(&serde_json::to_vec(&base).unwrap());
        envelope(&original, 0);
        let mut request = base;
        for field in [
            "motif_parameters",
            "analysis_parameters",
            "blast_parameters",
            "compression_parameters",
            "compression",
            "ph",
            "parameters",
        ] {
            request[field] = Value::Null;
        }
        let explicit_null = run(&serde_json::to_vec(&request).unwrap());
        envelope(&explicit_null, 0);
        assert_eq!(original.stdout, explicit_null.stdout);
    }
}

#[test]
fn unknown_duplicate_top_fields_and_positional_arrays_stay_rejected() {
    let [base, _] = sequence_requests();
    for raw in [
        b"[]".as_slice(),
        br#"[1,"translate","dna","ATGGCTTAA"]"#,
        br#"[{"schema_version":1,"operation":"translate","molecule":"dna","sequence":"ATGGCTTAA"}]"#,
    ] {
        rejected(raw);
    }
    for value in ["null", "{}", "1"] {
        let error = rejected(&with_raw_field(&base, "unknown_analysis_field", value));
        assert!(
            error["error"]["message"]
                .as_str()
                .unwrap()
                .contains("unknown_analysis_field")
        );
    }
    for field in [
        "motif_parameters",
        "analysis_parameters",
        "blast_parameters",
        "compression_parameters",
        "compression",
        "ph",
        "parameters",
    ] {
        let first = with_raw_field(&base, field, "null");
        let first: Value = serde_json::from_slice(&first).unwrap();
        let error = rejected(&with_raw_field(&first, field, "null"));
        let message = error["error"]["message"].as_str().unwrap();
        assert!(message.contains("duplicate field") && message.contains(field));
    }
}

#[test]
fn compression_selector_and_ph_are_not_general_request_options() {
    let mut unrelated = sequence_requests().to_vec();
    unrelated.extend(
        feature_cases()
            .into_iter()
            .filter(|case| case.request["operation"] != "compress-file")
            .map(|case| case.request),
    );
    unrelated.extend([
        json!({"schema_version":1,"operation":"inspect-file","molecule":"dna","format":"fasta","input":missing_input()}),
        json!({"schema_version":1,"operation":"index-bgzf","molecule":"dna","format":"fasta","input":missing_input(),"output":format!("{}.index.json",missing_input())}),
        json!({"schema_version":1,"operation":"get-bgzf-record","molecule":"dna","id":"one","input":missing_input(),"index":format!("{}.index.json",missing_input())}),
    ]);
    for mut request in unrelated {
        let mut engine_options = request.clone();
        engine_options["parameters"] = json!({});
        rejected_field(&engine_options, "parameters");
        request["compression"] = json!("gzip");
        rejected_field(&request, "compression");
    }
    let [mut translation, mut protein] = sequence_requests();
    translation["ph"] = json!(7);
    rejected_field(&translation, "ph");
    protein["ph"] = json!(7.4);
    let output = run(&serde_json::to_vec(&protein).unwrap());
    assert_eq!(envelope(&output, 0)["result"]["ph"], 7.4);
}
