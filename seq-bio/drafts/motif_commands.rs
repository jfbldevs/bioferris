//! Actual-binary contracts for IUPAC matching and nucleotide matrix workflows.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

use serde_json::{Value, json};

fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bioferris"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("actual CLI must start");
    child
        .stdin
        .take()
        .expect("stdin pipe")
        .write_all(input)
        .expect("test input");
    child.wait_with_output().expect("CLI must exit")
}

fn envelope(output: &Output, status: i32) -> Value {
    assert_eq!(output.status.code(), Some(status), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    let value: Value = serde_json::from_slice(&output.stdout).expect("one complete JSON envelope");
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], status == 0);
    assert_eq!(value.get("result").is_some(), status == 0);
    assert_eq!(value.get("error").is_some(), status != 0);
    value
}

fn machine(args: &[&str], status: i32) -> Value {
    let mut args = args.to_vec();
    args.push("--json");
    envelope(&run(&args, b""), status)
}

fn request(value: &Value, status: i32) -> Value {
    envelope(
        &run(
            &["request", "--input", "-"],
            &serde_json::to_vec(value).expect("request JSON"),
        ),
        status,
    )
}

fn fixture() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/motif-instances-dna.fa")
        .to_str()
        .expect("fixture path is UTF-8")
        .to_owned()
}

fn model_request(operation: &str, parameters: Value) -> Value {
    let mut request = json!({"schema_version": 1, "operation": operation, "molecule": "dna"});
    request["motif_parameters"] = parameters;
    if operation == "motif-scan" {
        request["sequence"] = json!("ACGTAC");
    }
    request
}

fn hit(strand: &str, start: usize, end: usize) -> Value {
    json!({"strand": strand, "range": {"start": start, "end": end}})
}

#[test]
fn iupac_definite_and_possible_modes_have_original_coordinates_and_sorted_strands() {
    let definite = machine(&["motif-find", "dna", "AGRYN", "r"], 0);
    assert_eq!(definite["result"]["pattern"], "R");
    assert_eq!(definite["result"]["ambiguity"], "definite");
    assert_eq!(
        definite["result"]["matches"],
        json!([
            hit("forward", 0, 1),
            hit("forward", 1, 2),
            hit("forward", 2, 3),
            hit("reverse", 3, 4)
        ])
    );
    let possible = machine(
        &["motif-find", "dna", "AGRYN", "R", "--ambiguity", "possible"],
        0,
    );
    assert_eq!(
        possible["result"]["matches"].as_array().expect("matches")[4..],
        [hit("forward", 4, 5), hit("reverse", 4, 5)]
    );
    let rejected = machine(
        &["motif-find", "dna", "AGRYN", "R", "--ambiguity", "reject"],
        2,
    );
    assert_eq!(rejected["error"]["code"], "ambiguous_motif_input");
    assert_eq!(rejected["error"]["details"]["offset"], 2);
}

#[test]
fn overlaps_palindromes_short_inputs_and_rna_are_explicit() {
    let overlap = machine(
        &["motif-find", "dna", "AAAA", "AA", "--strands", "forward"],
        0,
    );
    assert_eq!(
        overlap["result"]["matches"],
        json!([
            hit("forward", 0, 2),
            hit("forward", 1, 3),
            hit("forward", 2, 4)
        ])
    );
    for molecule in ["dna", "rna"] {
        let text = if molecule == "dna" { "ATAT" } else { "AUAU" };
        let pattern = if molecule == "dna" { "AT" } else { "AU" };
        let value = machine(&["motif-find", molecule, text, pattern], 0);
        assert_eq!(
            value["result"]["matches"],
            json!([
                hit("forward", 0, 2),
                hit("reverse", 0, 2),
                hit("forward", 2, 4),
                hit("reverse", 2, 4)
            ])
        );
    }
    for source in ["", "A"] {
        assert_eq!(
            machine(&["motif-find", "dna", source, "AAA"], 0)["result"]["matches"],
            json!([])
        );
    }
    assert_eq!(
        machine(&["motif-find", "dna", "AC", ""], 2)["error"]["code"],
        "empty_motif_pattern"
    );
    assert_eq!(
        machine(&["motif-find", "rna", "AC", "AT"], 2)["error"]["code"],
        "invalid_motif_pattern"
    );
}

#[test]
fn human_pattern_matrix_and_scan_reports_are_deterministic() {
    let output = run(
        &["motif-find", "dna", "AAAA", "AA", "--strands", "forward"],
        b"",
    );
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(
        output.stdout,
        b"strand\tstart\tend\nforward\t0\t2\nforward\t1\t3\nforward\t2\t4\n"
    );
    let output = run(
        &[
            "motif-scan",
            "dna",
            "ACGTAC",
            "--format",
            "fasta",
            "--input",
            &fixture(),
            "--pseudocounts",
            "1,1,1,1",
            "--threshold",
            "2",
        ],
        b"",
    );
    assert!(output.status.success(), "{output:?}");
    assert_eq!(
        output.stdout,
        b"strand\tstart\tend\tscore\nforward\t0\t2\t2\nreverse\t2\t4\t2\nforward\t4\t6\t2\n"
    );
    let output = run(
        &[
            "motif-matrix",
            "dna",
            "--format",
            "fasta",
            "--input",
            "-",
            "--pseudocounts",
            "1,1,1,1",
        ],
        b">one\nA\n",
    );
    assert!(output.status.success(), "{output:?}");
    let text = String::from_utf8(output.stdout).expect("TSV");
    assert!(text.starts_with("matrix\tposition\tbase\tvalue\ncounts\t0\tA\t1\n"));
    assert_eq!(text.lines().count(), 13);
}

#[test]
fn file_instances_build_all_three_matrices_and_normalize_background() {
    let value = machine(
        &[
            "motif-matrix",
            "dna",
            "--format",
            "fasta",
            "--input",
            &fixture(),
            "--pseudocounts",
            "1,1,1,1",
            "--background",
            "1,1,1,1",
        ],
        0,
    );
    let model = &value["result"]["model"];
    assert_eq!(model["alphabet"], json!(["A", "C", "G", "T"]));
    assert_eq!(model["width"], 2);
    assert_eq!(model["source"]["instance_count"], 4);
    assert_eq!(model["source"]["total_symbols"], 8);
    assert_eq!(
        model["counts"],
        json!([[3.0, 1.0, 0.0, 0.0], [1.0, 3.0, 0.0, 0.0]])
    );
    assert_eq!(
        model["weights"],
        json!([[0.5, 0.25, 0.125, 0.125], [0.25, 0.5, 0.125, 0.125]])
    );
    assert_eq!(
        model["scores"],
        json!([[1.0, 0.0, -1.0, -1.0], [0.0, 1.0, -1.0, -1.0]])
    );
    assert_eq!(model["background"], json!([0.25, 0.25, 0.25, 0.25]));
    assert_eq!(model["precision"]["normalization_sum"], "neumaier");
}

#[test]
fn matrix_default_zero_pseudocounts_preserve_negative_infinity_as_strings() {
    let value = envelope(
        &run(
            &[
                "motif-matrix",
                "dna",
                "--format",
                "fasta",
                "--input",
                "-",
                "--json",
            ],
            b">one\nA\n",
        ),
        0,
    );
    assert_eq!(
        value["result"]["model"]["pseudocounts"],
        json!([0.0, 0.0, 0.0, 0.0])
    );
    assert_eq!(
        value["result"]["model"]["scores"],
        json!([[2.0, "-inf", "-inf", "-inf"]])
    );
    assert!(
        !String::from_utf8_lossy(
            &serde_json::to_vec(&value["result"]["model"]["scores"]).expect("score JSON")
        )
        .contains("null")
    );
}

#[test]
fn instance_and_count_column_json_sources_match_the_file_model() {
    let file = machine(
        &[
            "motif-matrix",
            "dna",
            "--format",
            "fasta",
            "--input",
            &fixture(),
            "--pseudocounts",
            "1,1,1,1",
        ],
        0,
    );
    for source in [
        json!({"instances": ["AC", "AC", "AC", "CA"], "pseudocounts": [1,1,1,1]}),
        json!({"count_columns": [[3,1,0,0], [1,3,0,0]], "pseudocounts": [1,1,1,1]}),
    ] {
        let value = request(&model_request("motif-matrix", source), 0);
        for field in [
            "counts",
            "weights",
            "scores",
            "background",
            "alphabet",
            "precision",
        ] {
            assert_eq!(
                value["result"]["model"][field],
                file["result"]["model"][field]
            );
        }
    }
}

#[test]
fn worked_scan_matches_inline_request_and_checked_in_recipe() {
    let direct = machine(
        &[
            "motif-scan",
            "dna",
            "ACGTAC",
            "--format",
            "fasta",
            "--input",
            &fixture(),
            "--pseudocounts",
            "1,1,1,1",
            "--threshold",
            "2",
        ],
        0,
    );
    let expected = json!([
        {"strand":"forward", "range":{"start":0,"end":2}, "score":2.0},
        {"strand":"reverse", "range":{"start":2,"end":4}, "score":2.0},
        {"strand":"forward", "range":{"start":4,"end":6}, "score":2.0}
    ]);
    assert_eq!(direct["result"]["matches"], expected);
    let recipe: Value = serde_json::from_str(include_str!("../tests/fixtures/motif-request.json"))
        .expect("recipe JSON");
    let inline = request(&recipe, 0);
    assert_eq!(inline["result"]["matches"], expected);
    assert_eq!(
        inline["result"]["model"]["scores"],
        direct["result"]["model"]["scores"]
    );
}

#[test]
fn asymmetric_background_and_one_window_scans_are_supported() {
    let mut input = model_request(
        "motif-scan",
        json!({"count_columns":[[1,1,1,1]],"background":[4,2,1,1]}),
    );
    input["sequence"] = json!("A");
    let value = request(&input, 0);
    assert_eq!(
        value["result"]["matches"],
        json!([{"strand":"reverse", "range":{"start":0,"end":1}, "score":1.0}])
    );
    input["molecule"] = json!("rna");
    input["sequence"] = json!("U");
    let value = request(&input, 0);
    assert_eq!(
        value["result"]["model"]["alphabet"],
        json!(["A", "C", "G", "U"])
    );
    assert_eq!(value["result"]["matches"][0]["strand"], "forward");
}

#[test]
fn ambiguous_matrix_scan_or_instance_is_rejected_without_partial_success() {
    let mut input = model_request("motif-scan", json!({"instances":["AC","AC"]}));
    input["sequence"] = json!("ACN");
    let value = request(&input, 2);
    assert_eq!(value["error"]["code"], "ambiguous_motif_input");
    assert_eq!(value["error"]["details"]["offset"], 2);
    input = model_request("motif-matrix", json!({"instances":["AC","AN"]}));
    assert_eq!(request(&input, 2)["error"]["code"], "ambiguous_motif_input");
}

#[test]
fn model_dimensions_values_and_zero_total_columns_fail_explicitly() {
    for (parameters, code) in [
        (json!({"instances":[]}), "invalid_motif_matrix"),
        (json!({"instances":[""]}), "invalid_motif_matrix"),
        (json!({"instances":["AC","A"]}), "motif_dimension_mismatch"),
        (json!({"count_columns":[]}), "invalid_motif_matrix"),
        (json!({"count_columns":[[0,0,0,0]]}), "invalid_motif_matrix"),
        (
            json!({"count_columns":[[-1,0,0,0]]}),
            "invalid_motif_matrix",
        ),
        (
            json!({"instances":["A"],"pseudocounts":[-1,0,0,0]}),
            "invalid_motif_matrix",
        ),
        (
            json!({"instances":["A"],"background":[0,1,1,1]}),
            "invalid_motif_matrix",
        ),
    ] {
        assert_eq!(
            request(&model_request("motif-matrix", parameters), 2)["error"]["code"],
            code
        );
    }
    let value = request(
        &model_request(
            "motif-matrix",
            json!({"count_columns":[[0,0,0,0]],"pseudocounts":[1,1,1,1]}),
        ),
        0,
    );
    assert_eq!(
        value["result"]["model"]["weights"],
        json!([[0.25, 0.25, 0.25, 0.25]])
    );
}

#[test]
fn every_search_and_model_budget_is_enforced_globally() {
    for (flag, value, resource) in [
        ("--max-windows", "5", "windows"),
        ("--max-comparisons", "11", "comparisons"),
        ("--max-matches", "0", "matches"),
    ] {
        let result = machine(&["motif-find", "dna", "AAAA", "AA", flag, value], 2);
        assert_eq!(result["error"]["code"], "limit_exceeded");
        assert_eq!(result["error"]["details"]["resource"], resource);
    }
    for (flag, value) in [
        ("--max-instances", "3"),
        ("--max-symbols", "7"),
        ("--max-width", "1"),
    ] {
        assert_eq!(
            machine(
                &[
                    "motif-matrix",
                    "dna",
                    "--format",
                    "fasta",
                    "--input",
                    &fixture(),
                    flag,
                    value
                ],
                2
            )["error"]["code"],
            "limit_exceeded"
        );
    }
    for parameters in [
        json!({"instances":["AC"],"max_windows":9}),
        json!({"instances":["AC"],"max_comparisons":19}),
        json!({"instances":["AC"],"max_matches":0}),
    ] {
        assert_eq!(
            request(&model_request("motif-scan", parameters), 2)["error"]["code"],
            "limit_exceeded"
        );
    }
    assert_eq!(
        request(
            &model_request(
                "motif-matrix",
                json!({"count_columns":[[1,1,1,1],[1,1,1,1]],"max_width":1})
            ),
            2
        )["error"]["details"]["resource"],
        "width"
    );
}

#[test]
fn malformed_fasta_gaps_and_fastq_are_not_coerced_into_instances() {
    for text in [b"AC\n".as_slice(), b">one\nA-\n", b"@one\nAC\n+\nII\n"] {
        let value = envelope(
            &run(
                &[
                    "motif-matrix",
                    "dna",
                    "--format",
                    "fasta",
                    "--input",
                    "-",
                    "--json",
                ],
                text,
            ),
            2,
        );
        assert!(value.get("result").is_none());
    }
    let unequal = envelope(
        &run(
            &[
                "motif-matrix",
                "dna",
                "--format",
                "fasta",
                "--input",
                "-",
                "--json",
            ],
            b">one\nAC\n>two\nA\n",
        ),
        2,
    );
    assert_eq!(unequal["error"]["code"], "motif_dimension_mismatch");
    assert_eq!(unequal["error"]["details"]["instance_index"], 1);
    assert_eq!(unequal["error"]["details"]["model_source"]["kind"], "fasta");
}

#[test]
fn unknown_duplicate_and_inapplicable_direct_flags_fail() {
    for args in [
        vec!["motif-find", "dna", "A", "A", "--unknown", "1"],
        vec![
            "motif-find",
            "dna",
            "A",
            "A",
            "--strands",
            "forward",
            "--strands",
            "reverse",
        ],
        vec!["motif-find", "dna", "A", "A", "--threshold", "0"],
        vec!["motif-find", "dna", "A", "A", "--json", "--json"],
        vec!["motif-find", "dna", "A", "A", "--max-matches", "-1"],
        vec!["motif-find", "protein", "A", "A"],
        vec![
            "motif-matrix",
            "dna",
            "--format",
            "fasta",
            "--input",
            "-",
            "--strands",
            "forward",
        ],
        vec![
            "motif-matrix",
            "dna",
            "--format",
            "fasta",
            "--input",
            "-",
            "--pseudocounts",
            "1,2,3",
        ],
        vec![
            "motif-scan",
            "dna",
            "A",
            "--format",
            "fasta",
            "--input",
            "-",
            "--ambiguity",
            "possible",
        ],
        vec![
            "motif-scan",
            "dna",
            "A",
            "--format",
            "fasta",
            "--input",
            "-",
            "--threshold",
            "NaN",
        ],
    ] {
        assert_eq!(
            machine(&args, 2)["error"]["code"],
            "invalid_arguments",
            "{args:?}"
        );
    }
}

#[test]
fn strict_nested_request_fields_reject_unknown_duplicate_and_wrong_types() {
    for raw in [
        r#"{"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"A","motif_parameters":{"unknown":1}}"#,
        r#"{"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"A","motif_parameters":{"strands":"both","strands":"forward"}}"#,
        r#"{"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"A","motif_parameters":[]}"#,
        r#"{"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"A","motif_parameters":{"max_matches":"1"}}"#,
        r#"{"schema_version":1,"operation":"motif-matrix","molecule":"dna","motif_parameters":{"count_columns":[[1,2,3]]}}"#,
        r#"{"schema_version":1,"operation":"motif-matrix","molecule":"dna","motif_parameters":{"instances":[1]}}"#,
        r#"{"schema_version":1,"operation":"motif-matrix","molecule":"dna","motif_parameters":{"pseudocounts":[1,2,3,4,5]}}"#,
    ] {
        let value = envelope(&run(&["request", "--input", "-"], raw.as_bytes()), 2);
        assert_eq!(value["error"]["code"], "invalid_schema", "{raw}");
    }
}

#[test]
fn typed_request_strings_are_data_and_model_sources_are_exclusive() {
    let mut input = json!({"schema_version":1,"operation":"motif-find","molecule":"dna","sequence":"A","pattern":"--json"});
    assert_eq!(request(&input, 2)["error"]["code"], "invalid_motif_pattern");
    input["pattern"] = json!("A");
    input["motif_parameters"] = json!({"threshold":0});
    assert_eq!(request(&input, 2)["error"]["code"], "invalid_arguments");
    for parameters in [
        json!({}),
        json!({"instances":["A"],"count_columns":[[1,0,0,0]]}),
        json!({"count_columns":[[1,0,0,0]],"max_instances":1}),
        json!({"instances":["--json"]}),
    ] {
        let value = request(&model_request("motif-matrix", parameters), 2);
        assert!(matches!(
            value["error"]["code"].as_str(),
            Some("invalid_arguments" | "invalid_sequence")
        ));
    }
    let mut nested = model_request("motif-matrix", json!({}));
    nested["input"] = json!("-");
    nested["format"] = json!("fasta");
    assert_eq!(request(&nested, 2)["error"]["code"], "invalid_schema");
}

#[test]
fn gzip_instance_streams_and_integrity_errors_use_the_existing_reader() {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder
        .write_all(include_bytes!("../tests/fixtures/motif-instances-dna.fa"))
        .expect("synthetic gzip payload");
    let compressed = encoder.finish().expect("complete gzip member");
    let args = [
        "motif-matrix",
        "dna",
        "--format",
        "fasta",
        "--input",
        "-",
        "--json",
    ];
    let value = envelope(&run(&args, &compressed), 0);
    assert_eq!(value["result"]["model"]["source"]["instance_count"], 4);
    let broken = envelope(&run(&args, &compressed[..compressed.len() - 3]), 2);
    assert!(broken.get("result").is_none());
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn scan_json_retains_float32_score_bits_and_threshold_rounding() {
    let parameters =
        json!({"count_columns":[[1,2,3,5]], "background":[0.3,0.2,0.2,0.3], "threshold":-10});
    let mut input = model_request("motif-scan", parameters);
    input["sequence"] = json!("AA");
    let value = request(&input, 0);
    let score = value["result"]["matches"][0]["score"]
        .as_f64()
        .expect("finite score");
    let expected = f64::from_bits(0xbffb_8f38_88b9_310a) as f32;
    assert_eq!((score as f32).to_bits(), expected.to_bits());
    input["motif_parameters"]["threshold"] = json!(score + 1e-12);
    let rounded = request(&input, 0);
    assert_eq!(
        rounded["result"]["matches"].as_array().expect("hits").len(),
        4
    );
    let effective = rounded["result"]["effective_threshold"]
        .as_f64()
        .expect("finite threshold");
    assert_eq!((effective as f32).to_bits(), expected.to_bits());
    // Negative finite f32 values increase numerically when their bit pattern decreases.
    assert!(expected < 0.0);
    let next = f32::from_bits(expected.to_bits() - 1);
    input["motif_parameters"]["threshold"] = json!(f64::from(next));
    let higher = request(&input, 0);
    let hits = higher["result"]["matches"]
        .as_array()
        .expect("higher threshold hits");
    assert_eq!(hits.len(), 2);
    assert!(hits.iter().all(|hit| hit["strand"] == "reverse"));
}

#[test]
fn named_json_file_source_and_optional_null_policies_are_supported() {
    let mut input = model_request("motif-matrix", json!({"pseudocounts":[1,1,1,1]}));
    input["format"] = json!("fasta");
    input["input"] = json!(fixture());
    let value = request(&input, 0);
    assert_eq!(value["result"]["model"]["source"]["kind"], "fasta");
    let pattern = json!({"schema_version":1,"operation":"motif-find","molecule":"dna", "sequence":"A", "pattern":"A", "motif_parameters":{"strands":null,"threshold":null}});
    assert_eq!(
        request(&pattern, 0)["result"]["matches"],
        json!([hit("forward", 0, 1)])
    );
}
