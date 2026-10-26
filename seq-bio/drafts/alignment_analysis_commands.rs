//! Actual-binary contracts awaiting the parent stage 4 target promotion.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

use serde_json::{Value, json};

fn run(args: &[&str], input: &[u8]) -> Output {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut child = Command::new(env!("CARGO_BIN_EXE_bioferris"))
        .args(args)
        .current_dir(manifest.parent().expect("workspace root"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("CLI starts");
    child
        .stdin
        .take()
        .expect("stdin pipe")
        .write_all(input)
        .expect("test input");
    child.wait_with_output().expect("CLI exits")
}

fn envelope(output: &Output, code: i32) -> Value {
    assert_eq!(output.status.code(), Some(code), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    let value: Value = serde_json::from_slice(&output.stdout).expect("complete JSON envelope");
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], code == 0);
    assert_eq!(value.get("result").is_some(), code == 0);
    assert_eq!(value.get("error").is_some(), code != 0);
    value
}

fn machine(args: &[&str], status: i32) -> Value {
    let mut args = args.to_vec();
    args.push("--json");
    envelope(&run(&args, b""), status)
}

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
        .to_str()
        .expect("UTF-8 fixture path")
        .to_owned()
}

fn request(operation: &str, parameters: Value, status: i32) -> Value {
    let mut input = json!({"schema_version":1,"operation":operation,"molecule":"dna","format":"aligned-fasta","input":fixture("alignment-analysis.fa")});
    input["analysis_parameters"] = parameters;
    envelope(
        &run(
            &["request", "--input", "-"],
            &serde_json::to_vec(&input).expect("request JSON"),
        ),
        status,
    )
}

fn stdin(
    operation: &str,
    molecule: &str,
    format: &str,
    text: &[u8],
    policies: &[&str],
    status: i32,
) -> Value {
    let mut args = vec![
        operation, molecule, "--format", format, "--input", "-", "--json",
    ];
    args.extend_from_slice(policies);
    envelope(&run(&args, text), status)
}

#[test]
fn statistics_match_hand_counted_cells_pairs_gap_runs_and_identity() {
    let value = machine(
        &[
            "alignment-stats",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            &fixture("alignment-analysis.fa"),
        ],
        0,
    );
    let statistics = &value["result"]["statistics"];
    assert_eq!(statistics["rows"], 3);
    assert_eq!(statistics["columns"], 5);
    assert_eq!(
        statistics["cells"],
        json!({"residues":10,"gaps":5,"missing":0,"ambiguous":1})
    );
    assert_eq!(
        statistics["gap_runs"],
        json!({"leading":0,"internal":3,"trailing":1})
    );
    assert_eq!(
        statistics["pairs"],
        json!({"total":15,"missing":0,"excluded_ambiguity":0,"double_gap":3,"single_gap":4,"residue_pairs":8,"identities":7,"mismatches":1})
    );
    assert_eq!(statistics["denominator"], 8);
    assert_eq!(statistics["identity"], 0.875);
    assert_eq!(
        statistics["row_statistics"].as_array().expect("rows").len(),
        3
    );
    assert_eq!(
        statistics["column_statistics"]
            .as_array()
            .expect("columns")
            .len(),
        5
    );
    assert_eq!(statistics["column_statistics"][2]["pairs"]["double_gap"], 3);
    assert_eq!(statistics["row_statistics"][2]["identifier"], "ambiguous");
}

#[test]
fn statistics_policy_denominators_preserve_all_disjoint_pair_categories() {
    for (denominator, expected) in [
        ("residue-pairs", 7),
        ("residue-and-single-gap", 10),
        ("all-pairs", 15),
    ] {
        let value = request(
            "alignment-stats",
            json!({"ambiguity":"ignore","denominator":denominator}),
            0,
        );
        let statistics = &value["result"]["statistics"];
        assert_eq!(statistics["denominator"], expected);
        assert_eq!(statistics["pairs"]["identities"], 7);
        assert_eq!(statistics["pairs"]["excluded_ambiguity"], 2);
        assert_eq!(statistics["pairs"]["single_gap"], 3);
        assert_eq!(statistics["pairs"]["mismatches"], 0);
    }
    let rejected = request("alignment-stats", json!({"ambiguity":"reject"}), 2);
    assert_eq!(rejected["error"]["code"], "rejected_cell");
    assert_eq!(rejected["error"]["details"]["row"], 2);
    assert_eq!(rejected["error"]["details"]["column"], 1);
}

#[test]
fn every_explicit_alignment_profile_and_molecule_uses_the_existing_reader() {
    for molecule in ["dna", "rna", "protein"] {
        let datatype = match molecule {
            "dna" => "DNA",
            "rna" => "RNA",
            _ => "PROTEIN",
        };
        let nexus = format!(
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE={datatype} GAP=- MISSING=?;\nMATRIX\na AC-G\nb ACCG\n;\nEND;\n"
        );
        let cases = [
            ("aligned-fasta", ">a\nAC-G\n>b\nACCG\n"),
            ("clustal", "CLUSTAL W\n\na AC-G\nb ACCG\n"),
            (
                "phylip-strict-sequential",
                "2 4\na         AC\n-G\nb         AC\nCG\n",
            ),
            (
                "phylip-strict-interleaved",
                "2 4\na         AC\nb         AC\n\n-G\nCG\n",
            ),
            ("phylip-relaxed-sequential", "2 4\na AC\n-G\nb AC\nCG\n"),
            ("phylip-relaxed-interleaved", "2 4\na AC\nb AC\n\n-G\nCG\n"),
            ("stockholm", "# STOCKHOLM 1.0\na AC-G\nb ACCG\n//\n"),
            ("nexus", nexus.as_str()),
        ];
        for (format, text) in cases {
            let value = stdin("alignment-stats", molecule, format, text.as_bytes(), &[], 0);
            assert_eq!(
                value["result"]["statistics"]["rows"], 2,
                "{molecule}/{format}"
            );
            assert_eq!(
                value["result"]["statistics"]["columns"], 4,
                "{molecule}/{format}"
            );
            assert_eq!(
                value["result"]["statistics"]["pairs"]["identities"], 3,
                "{molecule}/{format}"
            );
        }
    }
}

#[test]
fn no_comparison_evidence_is_null_and_all_gap_runs_are_leading_only() {
    let value = stdin(
        "alignment-stats",
        "dna",
        "aligned-fasta",
        b">a\n---\n>b\n---\n",
        &[],
        0,
    );
    assert_eq!(value["result"]["statistics"]["identity"], Value::Null);
    assert_eq!(value["result"]["statistics"]["denominator"], 0);
    assert_eq!(
        value["result"]["statistics"]["gap_runs"],
        json!({"leading":2,"internal":0,"trailing":0})
    );
    assert_eq!(value["result"]["statistics"]["pairs"]["double_gap"], 3);
}

#[test]
fn consensus_reports_every_column_decision_and_empty_evidence() {
    let value = request("alignment-consensus", json!({}), 0);
    let consensus = &value["result"]["consensus"];
    assert_eq!(consensus["sequence"], "ANNGT");
    assert_eq!(
        consensus["columns"][1],
        json!({"column":1,"symbol":"N","status":"tie","winning_count":1,"denominator":2,"support":0.5,"winners":["C","N"]})
    );
    assert_eq!(
        consensus["columns"][2],
        json!({"column":2,"symbol":"N","status":"no-evidence","winning_count":0,"denominator":0,"support":null,"winners":[]})
    );
    assert_eq!(
        request("alignment-consensus", json!({"ambiguity":"ignore"}), 0)["result"]["consensus"]["sequence"],
        "ACNGT"
    );
    assert_eq!(
        request("alignment-consensus", json!({"gaps":"vote"}), 0)["result"]["consensus"]["sequence"],
        "AN-GT"
    );
    assert_eq!(
        request("alignment-consensus", json!({"gaps":"vote"}), 0)["result"]["consensus"]["columns"]
            [2]["status"],
        "gap"
    );
}

#[test]
fn tie_policy_precedes_threshold_and_ascii_resolution_must_meet_support() {
    let unknown = request("alignment-consensus", json!({"min_support":0.9}), 0);
    assert_eq!(
        unknown["result"]["consensus"]["columns"][1]["status"],
        "tie"
    );
    let rejected = request(
        "alignment-consensus",
        json!({"min_support":0.9,"ties":"error"}),
        2,
    );
    assert_eq!(rejected["error"]["code"], "consensus_tie");
    assert_eq!(rejected["error"]["details"]["winners"], json!(["C", "N"]));
    let low = request(
        "alignment-consensus",
        json!({"ties":"ascii-order","min_support":0.6}),
        0,
    );
    assert_eq!(
        low["result"]["consensus"]["columns"][1]["status"],
        "below-threshold"
    );
    let accepted = request(
        "alignment-consensus",
        json!({"ties":"ascii-order","min_support":0.5}),
        0,
    );
    assert_eq!(accepted["result"]["consensus"]["columns"][1]["symbol"], "C");
    assert_eq!(
        accepted["result"]["consensus"]["columns"][1]["status"],
        "tie-resolved"
    );
}

#[test]
fn invalid_support_errors_preserve_nonfinite_values_as_explicit_strings() {
    for (text, expected) in [("NaN", "NaN"), ("inf", "+inf"), ("-inf", "-inf")] {
        let value = stdin(
            "alignment-consensus",
            "dna",
            "aligned-fasta",
            b"not an alignment",
            &["--min-support", text],
            2,
        );
        assert_eq!(value["error"]["code"], "invalid_minimum_support");
        assert_eq!(value["error"]["details"]["value"], expected);
    }
    for value in [-0.1, 1.1] {
        assert_eq!(
            request("alignment-consensus", json!({"min_support":value}), 2)["error"]["code"],
            "invalid_minimum_support"
        );
    }
}

#[test]
fn missing_and_custom_gap_glyphs_remain_distinct_and_cannot_guess_coordinates() {
    let path = fixture("alignment-analysis.nex");
    let value = machine(
        &[
            "alignment-stats",
            "dna",
            "--format",
            "nexus",
            "--input",
            &path,
        ],
        0,
    );
    assert_eq!(value["result"]["source"]["gap_symbols"], json!(["~"]));
    assert_eq!(value["result"]["source"]["missing_symbol"], "!");
    assert_eq!(
        value["result"]["statistics"]["cells"],
        json!({"residues":5,"gaps":1,"missing":2,"ambiguous":0})
    );
    assert_eq!(value["result"]["statistics"]["pairs"]["missing"], 1);
    for column in ["0", "1", "2"] {
        let value = machine(
            &[
                "alignment-map",
                "dna",
                "--format",
                "nexus",
                "--input",
                &path,
                "--row",
                "0",
                "--column",
                column,
            ],
            2,
        );
        assert_eq!(value["error"]["code"], "missing_coordinates");
    }
    let sliced = machine(
        &[
            "alignment-slice",
            "dna",
            "--format",
            "nexus",
            "--input",
            &path,
            "--column-end",
            "1",
        ],
        0,
    );
    assert_eq!(sliced["result"]["alignment"]["rows"][0]["aligned"], "A");
    assert_eq!(
        sliced["result"]["alignment"]["rows"][0]["coordinates"],
        Value::Null
    );
    for (flag, kind) in [("--gaps", "gap"), ("--missing", "missing")] {
        let value = machine(
            &[
                "alignment-consensus",
                "dna",
                "--format",
                "nexus",
                "--input",
                &path,
                flag,
                "reject",
            ],
            2,
        );
        assert_eq!(value["error"]["code"], "rejected_cell");
        assert_eq!(value["error"]["details"]["kind"], kind);
    }
}

#[test]
fn coordinate_maps_distinguish_gaps_unaligned_positions_and_invalid_positions() {
    let value = request("alignment-map", json!({"row":0,"column":2}), 0);
    assert_eq!(value["result"]["mapping"]["position"], Value::Null);
    assert_eq!(value["result"]["mapping"]["column"], 2);
    let value = request("alignment-map", json!({"row":0,"position":2}), 0);
    assert_eq!(value["result"]["mapping"]["column"], 3);
    for (parameters, code) in [
        (json!({"row":3,"column":0}), "invalid_row"),
        (json!({"row":0,"column":5}), "invalid_column"),
        (json!({"row":0,"position":4}), "invalid_source_position"),
    ] {
        assert_eq!(
            request("alignment-map", parameters, 2)["error"]["code"],
            code
        );
    }
}

#[test]
fn explicit_reverse_origins_are_coordinate_assertions_without_reversing_cells() {
    let origin = json!({"row":0,"source_length":20,"start":10,"end":14,"strand":"reverse"});
    for (column, position) in [(0, 13), (1, 12), (3, 11), (4, 10)] {
        let value = request(
            "alignment-map",
            json!({"row":0,"column":column,"origins":[origin]}),
            0,
        );
        assert_eq!(value["result"]["mapping"]["position"], position);
        assert_eq!(
            value["result"]["mapping"]["coordinates"]["strand"],
            "reverse"
        );
    }
    let unaligned = request(
        "alignment-map",
        json!({"row":0,"position":0,"origins":[origin]}),
        0,
    );
    assert_eq!(unaligned["result"]["mapping"]["column"], Value::Null);
    let slice = request(
        "alignment-slice",
        json!({"row_end":1,"column_start":1,"column_end":4,"origins":[origin]}),
        0,
    );
    let row = &slice["result"]["alignment"]["rows"][0];
    assert_eq!(row["aligned"], "C-G");
    assert_eq!(row["coordinates"]["range"], json!({"start":11,"end":13}));
    assert_eq!(row["coordinates"]["source_length"], 20);
}

#[test]
fn origin_arrays_reject_duplicates_unknown_or_invalid_ranges_and_allow_distinct_cli_rows() {
    let origin = json!({"row":0,"source_length":20,"start":10,"end":14,"strand":"reverse"});
    let duplicate = request(
        "alignment-map",
        json!({"row":0,"column":0,"origins":[origin,origin]}),
        2,
    );
    assert_eq!(duplicate["error"]["code"], "invalid_arguments");
    for bad in [
        json!({"row":0,"source_length":4,"start":0,"end":5,"strand":"forward"}),
        json!({"row":0,"source_length":20,"start":10,"end":13,"strand":"reverse"}),
    ] {
        assert_eq!(
            request(
                "alignment-map",
                json!({"row":0,"column":0,"origins":[bad]}),
                2
            )["error"]["code"],
            "invalid_coordinates"
        );
    }
    let value = machine(
        &[
            "alignment-map",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            &fixture("alignment-analysis.fa"),
            "--row",
            "1",
            "--column",
            "0",
            "--origin",
            "0,20,10,14,reverse",
            "--origin",
            "1,200,100,103,forward",
        ],
        0,
    );
    assert_eq!(value["result"]["mapping"]["position"], 100);
}

#[test]
fn slices_preserve_source_positions_and_empty_dimension_shapes() {
    let value = request(
        "alignment-slice",
        json!({"row_end":2,"column_start":1,"column_end":4}),
        0,
    );
    let rectangle = &value["result"]["alignment"];
    assert_eq!(rectangle["row_count"], 2);
    assert_eq!(rectangle["width"], 3);
    assert_eq!(rectangle["rows"][0]["aligned"], "C-G");
    assert_eq!(rectangle["rows"][1]["aligned"], "--G");
    assert_eq!(
        rectangle["rows"][0]["coordinates"]["range"],
        json!({"start":1,"end":3})
    );
    assert_eq!(
        rectangle["rows"][1]["coordinates"]["range"],
        json!({"start":1,"end":2})
    );
    let empty_rows = request(
        "alignment-slice",
        json!({"row_start":1,"row_end":1,"column_end":3}),
        0,
    );
    assert_eq!(empty_rows["result"]["alignment"]["width"], 3);
    assert_eq!(empty_rows["result"]["alignment"]["rows"], json!([]));
    let empty_columns = request(
        "alignment-slice",
        json!({"column_start":2,"column_end":2}),
        0,
    );
    assert_eq!(empty_columns["result"]["alignment"]["row_count"], 3);
    assert_eq!(empty_columns["result"]["alignment"]["width"], 0);
    assert_eq!(
        empty_columns["result"]["alignment"]["rows"][0]["aligned"],
        ""
    );
    assert_eq!(
        empty_columns["result"]["alignment"]["rows"][0]["coordinates"]["range"],
        json!({"start":2,"end":2})
    );
}

#[test]
fn annotated_slices_default_to_rejection_and_discard_returns_full_ordered_evidence() {
    let path = fixture("alignment-analysis.sto");
    let identity = machine(
        &[
            "alignment-slice",
            "dna",
            "--format",
            "stockholm",
            "--input",
            &path,
            "--annotation-policy",
            "discard",
        ],
        0,
    );
    let annotations = identity["result"]["alignment"]["annotations"]
        .as_array()
        .expect("ordered annotations");
    assert_eq!(annotations.len(), 5);
    assert_eq!(identity["result"]["discarded_annotations"], json!([]));
    assert_eq!(
        annotations[0],
        json!({"scope":"file","sequence":null,"tag":"ID","line":2,"value":"original"})
    );
    assert_eq!(annotations[1]["scope"], "comment");
    assert_eq!(annotations[2]["scope"], "sequence");
    assert_eq!(annotations[3]["scope"], "column");
    assert_eq!(annotations[4]["scope"], "residue");
    let rejected = machine(
        &[
            "alignment-slice",
            "dna",
            "--format",
            "stockholm",
            "--input",
            &path,
            "--row-end",
            "1",
        ],
        2,
    );
    assert_eq!(rejected["error"]["code"], "annotation_loss");
    assert_eq!(
        rejected["error"]["details"]["annotations"],
        json!(annotations)
    );
    let discarded = machine(
        &[
            "alignment-slice",
            "dna",
            "--format",
            "stockholm",
            "--input",
            &path,
            "--row-end",
            "1",
            "--column-start",
            "1",
            "--column-end",
            "3",
            "--annotation-policy",
            "discard",
        ],
        0,
    );
    assert_eq!(
        discarded["result"]["discarded_annotations"],
        json!(annotations)
    );
    assert_eq!(discarded["result"]["alignment"]["annotations"], json!([]));
    assert_eq!(discarded["result"]["alignment"]["rows"][0]["aligned"], "C-");
}

#[test]
fn row_column_and_cell_limits_reject_without_partial_results() {
    for (field, value, resource) in [
        ("max_rows", 2, "rows"),
        ("max_columns", 4, "columns"),
        ("max_cells", 14, "cells"),
    ] {
        let mut parameters = json!({});
        parameters[field] = json!(value);
        let result = request("alignment-stats", parameters, 2);
        assert_eq!(result["error"]["code"], "limit_exceeded");
        assert_eq!(result["error"]["details"]["resource"], resource);
        assert_eq!(result["error"]["details"]["rows"], 3);
        assert_eq!(result["error"]["details"]["columns"], 5);
    }
    assert_eq!(
        request("alignment-slice", json!({"column_end":6}), 2)["error"]["code"],
        "invalid_range"
    );
    assert_eq!(
        request("alignment-slice", json!({"row_start":2,"row_end":1}), 2)["error"]["code"],
        "invalid_arguments"
    );
}

#[test]
fn strict_nested_origins_and_policies_reject_unknown_duplicate_and_wrong_json_types() {
    let prefix = format!(
        "{{\"schema_version\":1,\"operation\":\"alignment-map\",\"molecule\":\"dna\",\"format\":\"aligned-fasta\",\"input\":{},\"analysis_parameters\":",
        json!(fixture("alignment-analysis.fa"))
    );
    for parameters in [
        "[]",
        r#"{"row":0,"column":0,"unknown":1}"#,
        r#"{"row":0,"row":1,"column":0}"#,
        r#"{"row":"0","column":0}"#,
        r#"{"row":0,"column":0,"origins":[[0,20,10,14,"reverse"]]}"#,
        r#"{"row":0,"column":0,"origins":[{"row":0,"source_length":20,"start":10,"end":14,"strand":"reverse","unknown":1}]}"#,
        r#"{"row":0,"column":0,"origins":[{"row":0,"row":1,"source_length":20,"start":10,"end":14,"strand":"reverse"}]}"#,
        r#"{"row":0,"column":0,"origins":[{"row":0,"source_length":20,"start":10,"end":14}]}"#,
    ] {
        let text = format!("{prefix}{parameters}}}");
        assert_eq!(
            envelope(&run(&["request", "--input", "-"], text.as_bytes()), 2)["error"]["code"],
            "invalid_schema",
            "{parameters}"
        );
    }
}

#[test]
fn inapplicable_options_missing_map_selectors_and_duplicate_flags_are_rejected() {
    for (operation, parameters) in [
        ("alignment-stats", json!({"gaps":"ignore"})),
        ("alignment-consensus", json!({"denominator":"all-pairs"})),
        ("alignment-slice", json!({"min_support":0.5})),
        ("alignment-map", json!({"row":0})),
        ("alignment-map", json!({"row":0,"column":0,"position":0})),
        ("alignment-map", json!({"column":0})),
    ] {
        assert_eq!(
            request(operation, parameters, 2)["error"]["code"],
            "invalid_arguments"
        );
    }
    for args in [
        vec![
            "alignment-stats",
            "dna",
            "--format",
            "aligned-fasta",
            "--format",
            "clustal",
            "--input",
            "-",
        ],
        vec![
            "alignment-stats",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            "-",
            "--output",
            "new.fa",
        ],
        vec![
            "alignment-stats",
            "dna",
            "--format",
            "fasta",
            "--input",
            "-",
        ],
        vec![
            "alignment-map",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            "-",
            "--row",
            "0",
            "--column",
            "0",
            "--origin",
            "0,4,0,4,forward",
            "--origin",
            "0,4,0,4,forward",
        ],
        vec![
            "alignment-stats",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            "-",
            "--json",
            "--json",
        ],
        vec![
            "alignment-stats",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            "-",
            "--max-rows",
            "-1",
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
fn gzip_streams_and_errors_preserve_complete_result_envelopes() {
    let source = include_bytes!("../tests/fixtures/alignment-analysis.fa");
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(source).expect("gzip fixture");
    let bytes = encoder.finish().expect("finished gzip");
    let value = stdin("alignment-stats", "dna", "aligned-fasta", &bytes, &[], 0);
    assert_eq!(value["result"]["statistics"]["cells"]["residues"], 10);
    let broken = stdin(
        "alignment-stats",
        "dna",
        "aligned-fasta",
        &bytes[..bytes.len() - 3],
        &[],
        2,
    );
    assert!(broken.get("result").is_none());
}

#[test]
fn human_reports_are_labeled_and_checked_in_request_keeps_reverse_origins() {
    let path = fixture("alignment-analysis.fa");
    let output = run(
        &[
            "alignment-map",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            &path,
            "--row",
            "0",
            "--column",
            "2",
        ],
        b"",
    );
    assert!(output.status.success(), "{output:?}");
    assert!(output.stderr.is_empty());
    assert_eq!(output.stdout,b"row\tidentifier\tdirection\tcolumn\tposition\n0\t\"reference\"\tcolumn-to-sequence\t2\t.\n");
    let output = run(
        &[
            "alignment-consensus",
            "dna",
            "--format",
            "aligned-fasta",
            "--input",
            &path,
        ],
        b"",
    );
    assert!(String::from_utf8(output.stdout).expect("human TSV").starts_with("consensus\tANNGT\ncolumn\tsymbol\tstatus\twinning_count\tdenominator\tsupport\twinners\n"));
    let recipe = include_bytes!("../tests/fixtures/alignment-analysis-request.json");
    let value = envelope(&run(&["request", "--input", "-"], recipe), 0);
    assert_eq!(value["result"]["alignment"]["rows"][0]["aligned"], "C-G");
    assert_eq!(
        value["result"]["alignment"]["rows"][0]["coordinates"]["range"],
        json!({"start":11,"end":13})
    );
}

#[test]
fn rna_and_protein_consensus_use_molecule_specific_unknown_symbols() {
    let rna = stdin(
        "alignment-consensus",
        "rna",
        "aligned-fasta",
        b">a\nAC\n>b\nUC\n",
        &[],
        0,
    );
    assert_eq!(rna["result"]["consensus"]["sequence"], "NC");
    let protein = stdin(
        "alignment-consensus",
        "protein",
        "aligned-fasta",
        b">a\nA-\n>b\nC-\n",
        &[],
        0,
    );
    assert_eq!(protein["result"]["consensus"]["sequence"], "XX");
    assert_eq!(
        protein["result"]["consensus"]["columns"][0]["status"],
        "tie"
    );
    assert_eq!(
        protein["result"]["consensus"]["columns"][1]["status"],
        "no-evidence"
    );
    let ambiguity = stdin(
        "alignment-stats",
        "protein",
        "aligned-fasta",
        b">a\nAX\n>b\nA-\n",
        &[],
        0,
    );
    assert_eq!(ambiguity["result"]["statistics"]["cells"]["ambiguous"], 1);
}

#[test]
fn malformed_documents_return_context_without_a_partial_analysis() {
    for text in [
        b">a\nAC-\n>b\nAC\n".as_slice(),
        b">a\nAC$\n",
        b"not an alignment\n",
    ] {
        let value = stdin("alignment-stats", "dna", "aligned-fasta", text, &[], 2);
        assert_eq!(value["error"]["details"]["input"], "-");
        assert!(value.get("result").is_none());
    }
}
