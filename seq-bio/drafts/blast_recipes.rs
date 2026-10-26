//! Execute complete documented BLAST inputs, outputs, requests and Cargo example.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

use serde_json::Value;

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/blast-input-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- blast-{name} -->");
    assert_eq!(guide.matches(&marker).count(), 1, "unique marker {marker}");
    guide
        .split_once(&marker)
        .unwrap()
        .1
        .split_once(&format!("```{language}\n"))
        .unwrap()
        .1
        .split_once("```")
        .unwrap()
        .0
}

#[test]
fn complete_literal_inputs_program_and_request_match_promoted_files() {
    let guide = guide();
    for (name, file, language) in [
        ("input", "blast-example.tab6", "text"),
        ("commented-input", "blast-example.tab7", "text"),
        ("minimal-input", "blast-minimal.tsv", "text"),
        ("request-input", "blast-request.json", "json"),
    ] {
        assert_eq!(
            fence(&guide, name, language),
            fs::read_to_string(root().join("seq-bio/tests/fixtures").join(file)).unwrap()
        );
    }
    let program = fs::read_to_string(root().join("seq-bio/examples/blast_analysis.rs")).unwrap();
    assert_eq!(fence(&guide, "rust-program", "rust"), program);
    let literal = program
        .split_once("const INPUT: &str = \"")
        .unwrap()
        .1
        .split_once("\";")
        .unwrap()
        .0;
    let expected = fence(&guide, "input", "text")
        .replace('\t', "\\t")
        .replace('\n', "\\n");
    assert_eq!(literal, expected);
}

#[test]
fn all_six_documented_cli_recipes_match_complete_streams_status_and_schema() {
    let guide = guide();
    for (command, output, language, status) in [
        ("human", "human", "text", 0),
        ("json", "json", "json", 0),
        ("request", "json", "json", 0),
        ("xml", "xml", "text", 0),
        ("xml2", "xml", "text", 0),
        ("limit", "limit", "json", 2),
    ] {
        let invocation = fence(&guide, &format!("{command}-command"), "bash").trim();
        assert_eq!(invocation.lines().count(), 1);
        let arguments = invocation
            .strip_prefix("./target/debug/bioferris ")
            .unwrap();
        let observed = Command::new(env!("CARGO_BIN_EXE_bioferris"))
            .args(arguments.split_whitespace())
            .current_dir(root())
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            observed.status.code(),
            Some(status),
            "{invocation}: {observed:?}"
        );
        assert!(observed.stderr.is_empty(), "{invocation}: {observed:?}");
        let expected = fence(&guide, &format!("{output}-output"), language);
        assert_eq!(observed.stdout, expected.as_bytes(), "{invocation}");
        if language == "json" {
            assert_eq!(String::from_utf8_lossy(&observed.stdout).lines().count(), 1);
            let value: Value = serde_json::from_slice(&observed.stdout).unwrap();
            assert_eq!(value["schema_version"], 1);
            assert_eq!(value["ok"], status == 0);
            assert_eq!(value.get("result").is_some(), status == 0);
            assert_eq!(value.get("error").is_some(), status != 0);
        }
    }
}

#[test]
fn real_cargo_alias_matches_the_complete_manual_program_output() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "rust-command", "bash").trim(),
        "cargo bio blast_analysis"
    );
    let observed = Command::new(env!("CARGO"))
        .args(["bio", "blast_analysis", "--offline", "--locked"])
        .current_dir(root())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(observed.status.success(), "{observed:?}");
    assert!(observed.stderr.is_empty(), "{observed:?}");
    assert_eq!(
        observed.stdout,
        fence(&guide, "rust-output", "text").as_bytes()
    );
}
