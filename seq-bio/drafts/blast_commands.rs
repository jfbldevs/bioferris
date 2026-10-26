//! Actual-binary contracts for explicit BLAST result import profiles.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

use serde_json::{Value, json};

fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bioferris"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn envelope(output: &Output, status: i32) -> Value {
    assert_eq!(output.status.code(), Some(status), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], status == 0);
    assert_eq!(value.get("result").is_some(), status == 0);
    assert_eq!(value.get("error").is_some(), status != 0);
    value
}

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
        .to_str()
        .unwrap()
        .to_owned()
}

fn machine(format: &str, file: &str, extra: &[&str], status: i32) -> Value {
    let mut args = vec![
        "inspect-blast",
        "--format",
        format,
        "--input",
        file,
        "--json",
    ];
    args.extend_from_slice(extra);
    envelope(&run(&args, b""), status)
}

fn request(format: &str, file: &str, parameters: &Value, status: i32) -> Value {
    envelope(
        &run(
            &["request", "--input", "-"],
            &serde_json::to_vec(&&json!({
                "schema_version": 1, "operation": "inspect-blast", "input": file, "format": format,
                "blast_parameters": parameters,
            }))
            .unwrap(),
        ),
        status,
    )
}

struct Temporary(PathBuf);
impl Temporary {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-blast-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        file.write_all(bytes).unwrap();
        Self(path)
    }
    fn path(&self) -> &str {
        self.0.to_str().unwrap()
    }
}
impl Drop for Temporary {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).unwrap();
    }
}

#[test]
fn standard_tabular_json_preserves_all_source_fields_and_coordinates() {
    let result = machine("blast-tab6", &fixture("blast-ncbi-cookbook.tsv"), &[], 0);
    let result = &result["result"];
    assert_eq!(result["operation"], "inspect-blast");
    assert_eq!(
        result["coordinate_convention"],
        "zero-based-half-open-source-units"
    );
    assert_eq!(result["queries"][0]["id"], "AE000111.1");
    assert_eq!(result["queries"][0]["hits"].as_array().unwrap().len(), 3);
    let hsp = &result["queries"][0]["hits"][1]["hsps"][0];
    assert_eq!(
        hsp["hit_coordinates"],
        json!({"start": 6820, "end": 6928, "reported_start": 6928, "reported_end": 6821, "direction": "decreasing"})
    );
    assert_eq!(
        hsp["source_fields"][10],
        json!({"name": "evalue", "value": "7.69e-30"})
    );
    assert_eq!(hsp["hit_strand"], Value::Null);
    assert_eq!(result["report_metadata"], json!([]));
}

#[test]
fn commented_queries_keep_no_hits_headers_comments_and_processed_footer() {
    let result = machine("blast-tab7", &fixture("blast-example.tab7"), &[], 0);
    let result = &result["result"];
    assert_eq!(result["processed_queries"], 2);
    assert_eq!(result["queries"].as_array().unwrap().len(), 2);
    assert_eq!(result["queries"][0]["program"], "BLASTN");
    assert_eq!(result["queries"][0]["target"], "synthetic");
    assert_eq!(result["queries"][1]["id"], "absent");
    assert_eq!(result["queries"][1]["hits"], json!([]));
    assert_eq!(result["queries"][1]["identity_field"], "query-header");
    assert_eq!(result["queries"][1]["declared_rows"], 0);
    assert_eq!(result["comments"][0]["text"], "# BLAST processed 2 queries");
}

#[test]
fn human_summary_is_complete_and_marks_absent_values_with_a_dot() {
    let file = fixture("blast-example.tab7");
    let output = run(
        &["inspect-blast", "--format", "blast-tab7", "--input", &file],
        b"",
    );
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "query_id\thit_id\thsp_index\tquery_start\tquery_end\thit_start\thit_end\tevalue\tbit_score\nq\ts\t1\t0\t4\t4\t8\t0.00000000000000000001\t42\nabsent\t.\t.\t.\t.\t.\t.\t.\t.\n"
    );
}

#[test]
fn custom_fields_and_delimiter_are_literal_and_do_not_invent_optional_values() {
    let file = Temporary::new(b"q@s@custom value@1e-400\n");
    let direct = machine(
        "blast-tab6",
        file.path(),
        &["--fields", "qseqid,sseqid,note,evalue", "--delimiter", "@"],
        0,
    );
    let typed = request(
        "blast-tab6",
        file.path(),
        &json!({"fields": ["qseqid", "sseqid", "note", "evalue"], "delimiter": "@"}),
        0,
    );
    assert_eq!(direct, typed);
    let hsp = &direct["result"]["queries"][0]["hits"][0]["hsps"][0];
    assert_eq!(hsp["evalue"], 0.0);
    assert_eq!(hsp["source_fields"][3]["value"], "1e-400");
    for name in [
        "bit_score",
        "query_coordinates",
        "hit_coordinates",
        "query_strand",
        "identical",
    ] {
        assert_eq!(hsp[name], Value::Null);
    }
}

#[test]
fn all_four_formats_have_equal_direct_and_typed_request_results() {
    for (format, name) in [
        ("blast-tab6", "blast-ncbi-cookbook.tsv"),
        ("blast-tab7", "blast-example.tab7"),
        ("blast-xml", "blast-example.xml"),
        ("blast-xml2", "blast-example-xml2.xml"),
    ] {
        let file = fixture(name);
        assert_eq!(
            machine(format, &file, &[], 0),
            request(format, &file, &json!({}), 0)
        );
    }
}

#[test]
fn xml_report_metadata_uses_one_deterministic_pool_entry_per_shared_report() {
    let file = fixture("blast-example.xml");
    let one = machine("blast-xml", &file, &[], 0);
    let two = machine("blast-xml", &file, &[], 0);
    assert_eq!(one, two);
    let result = &one["result"];
    assert_eq!(result["report_metadata"].as_array().unwrap().len(), 1);
    assert_eq!(result["queries"][0]["report_metadata_index"], 0);
    assert_eq!(result["queries"][1]["report_metadata_index"], 0);
    assert!(
        !result["source_fields"]
            .as_array()
            .unwrap()
            .iter()
            .any(|field| field["name"] == "BlastOutput_reference")
    );
    let result = machine("blast-xml2", &fixture("blast-example-xml2.xml"), &[], 0);
    assert_eq!(
        result["result"]["report_metadata"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(result["result"]["queries"][0]["report_metadata_index"], 0);
    assert_eq!(result["result"]["queries"][1]["report_metadata_index"], 1);
}

#[test]
fn public_xml_profiles_preserve_translated_units_and_multiple_descriptions() {
    for (format, file) in [
        ("blast-xml", "blast-ncbi-2.9.0-translated.xml"),
        ("blast-xml2", "blast-ncbi-2.9.0-translated-xml2.xml"),
    ] {
        let value = machine(format, &fixture(file), &[], 0);
        assert_eq!(value["result"]["queries"][0]["program"], "blastx");
        let hsp = &value["result"]["queries"][0]["hits"][0]["hsps"][0];
        assert!(hsp["query_frame"].is_number());
        assert!(!hsp["query_alignment"].as_str().unwrap().is_empty());
        assert_eq!(hsp["identity_percent"], Value::Null);
    }
}

#[test]
fn gzip_is_detected_by_content_and_corrupt_trailers_are_not_successes() {
    let bytes = std::fs::read(fixture("blast-example.tab7")).unwrap();
    let mut compressor = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    compressor.write_all(&bytes).unwrap();
    let mut compressed = compressor.finish().unwrap();
    let file = Temporary::new(&compressed);
    assert_eq!(
        machine("blast-tab7", file.path(), &[], 0),
        machine("blast-tab7", &fixture("blast-example.tab7"), &[], 0)
    );
    let end = compressed.len();
    compressed[end - 8] ^= 1;
    let corrupt = Temporary::new(&compressed);
    let value = machine("blast-tab7", corrupt.path(), &[], 2);
    assert_eq!(value["error"]["code"], "io_error");
}

#[test]
fn all_common_resource_limits_and_xml_limits_fail_without_partial_results() {
    for flag in [
        "--max-line-bytes",
        "--max-total-bytes",
        "--max-queries",
        "--max-hits-per-query",
        "--max-hsps-per-query",
        "--max-query-bytes",
        "--max-query-id-bytes",
    ] {
        let value = machine(
            "blast-tab6",
            &fixture("blast-ncbi-cookbook.tsv"),
            &[flag, "0"],
            2,
        );
        assert_eq!(value["error"]["code"], "limit_exceeded", "{flag}: {value}");
    }
    for flag in ["--max-depth", "--max-event-bytes", "--max-metadata-bytes"] {
        let value = machine("blast-xml", &fixture("blast-example.xml"), &[flag, "0"], 2);
        assert_eq!(value["error"]["code"], "limit_exceeded", "{flag}: {value}");
    }
}

#[test]
fn resource_and_field_options_are_scoped_to_the_selected_profile() {
    for (format, file, options) in [
        (
            "blast-tab6",
            "blast-ncbi-cookbook.tsv",
            vec!["--max-depth", "64"],
        ),
        (
            "blast-tab7",
            "blast-example.tab7",
            vec!["--max-event-bytes", "8"],
        ),
        (
            "blast-tab6",
            "blast-ncbi-cookbook.tsv",
            vec!["--max-metadata-bytes", "8"],
        ),
        ("blast-xml", "blast-example.xml", vec!["--fields", "std"]),
        (
            "blast-xml2",
            "blast-example-xml2.xml",
            vec!["--delimiter", "tab"],
        ),
    ] {
        assert_eq!(
            machine(format, &fixture(file), &options, 2)["error"]["code"],
            "invalid_arguments"
        );
    }
}

#[test]
fn invalid_delimiters_repeated_flags_and_wrong_formats_are_explicit() {
    let path = fixture("blast-ncbi-cookbook.tsv");
    for args in [
        vec!["--delimiter", "##"],
        vec!["--delimiter", "β"],
        vec!["--max-queries", "-1"],
        vec!["--input", &path],
        vec!["--json"],
        vec!["--unknown", "x"],
    ] {
        assert_eq!(
            machine("blast-tab6", &path, &args, 2)["error"]["code"],
            "invalid_arguments"
        );
    }
    assert_eq!(
        machine("blast-tabular", &path, &[], 2)["error"]["code"],
        "invalid_arguments"
    );
    assert_eq!(
        machine("blast-tab6", "-", &[], 2)["error"]["code"],
        "invalid_arguments"
    );
    assert!(
        machine("blast-xml2", &fixture("blast-example.xml"), &[], 2)
            .get("error")
            .is_some()
    );
}

#[test]
fn malformed_input_and_missing_paths_have_stable_contextual_errors() {
    let malformed = Temporary::new(b"query\tsubject\n");
    let value = machine("blast-tab6", malformed.path(), &[], 2);
    assert_eq!(value["error"]["code"], "length_mismatch");
    assert_eq!(value["error"]["details"]["line"], 1);
    assert_eq!(value["error"]["details"]["record"], 1);
    let missing = machine("blast-tab6", &fixture("blast-file-does-not-exist"), &[], 2);
    assert_eq!(missing["error"]["code"], "io_error");
}

#[test]
fn nested_parameter_unknown_duplicate_and_wrong_types_fail_as_requests() {
    let path = serde_json::to_string(&fixture("blast-example.tab7")).unwrap();
    for parameters in [
        r#"{"unknown":1}"#,
        r#"{"max_queries":1,"max_queries":2}"#,
        r#"{"max_queries":"1"}"#,
        r#"{"fields":"std"}"#,
        r#"{"delimiter":9}"#,
    ] {
        let input = format!(
            r#"{{"schema_version":1,"operation":"inspect-blast","format":"blast-tab7","input":{path},"blast_parameters":{parameters}}}"#
        );
        let value = envelope(&run(&["request", "--input", "-"], input.as_bytes()), 2);
        assert_eq!(value["error"]["code"], "invalid_schema");
    }
}

fn documented_block<'a>(document: &'a str, marker: &str, language: &str) -> &'a str {
    let tail = document.split_once(marker).unwrap().1;
    let open = format!("```{language}\n");
    tail.split_once(&open)
        .unwrap()
        .1
        .split_once("\n```")
        .unwrap()
        .0
}

#[test]
fn contract_recipes_check_complete_literal_inputs_outputs_and_failure_envelopes() {
    let document = include_str!("../../docs/algorithms/blast-input-v1.md");
    let input = documented_block(
        document,
        "### Commented result with a hit and a no-hit query",
        "text",
    );
    assert_eq!(
        std::fs::read_to_string(fixture("blast-example.tab7")).unwrap(),
        format!("{input}\n")
    );
    let expected = documented_block(document, "Observed stdout; empty stderr, exit 0:", "text");
    let path = fixture("blast-example.tab7");
    let output = run(
        &["inspect-blast", "--format", "blast-tab7", "--input", &path],
        b"",
    );
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("{expected}\n")
    );
    let expected: Value = serde_json::from_str(documented_block(
        document,
        "The following complete stdout value was observed",
        "json",
    ))
    .unwrap();
    let actual = machine(
        "blast-tab6",
        &fixture("blast-minimal.tsv"),
        &["--fields", "qseqid,sseqid"],
        0,
    );
    assert_eq!(actual, expected);
    let expected: Value = serde_json::from_str(documented_block(
        document,
        "Complete stdout; empty stderr, exit 2:",
        "json",
    ))
    .unwrap();
    assert_eq!(
        machine(
            "blast-tab6",
            &fixture("blast-ncbi-cookbook.tsv"),
            &["--max-queries", "0"],
            2
        ),
        expected
    );
}
