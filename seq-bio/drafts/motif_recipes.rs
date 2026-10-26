//! Verify complete documented programs, literal inputs and real command results.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/motifs-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- motifs-{name} -->");
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
fn literal_program_and_fixtures_match_their_complete_documented_inputs() {
    let guide = guide();
    let program = fs::read_to_string(root().join("seq-bio/examples/motif_analysis.rs")).unwrap();
    assert_eq!(fence(&guide, "rust-program", "rust"), program);
    for (marker, language, filename) in [
        ("instances-input", "text", "motif-instances-dna.fa"),
        ("request-input", "json", "motif-request.json"),
        ("columns-input", "json", "motif-recipe-columns-request.json"),
    ] {
        let input =
            fs::read_to_string(root().join("seq-bio/tests/fixtures").join(filename)).unwrap();
        assert_eq!(fence(&guide, marker, language), input, "{filename}");
    }
}

#[test]
fn all_eight_cli_recipes_match_byte_exact_output_status_and_streams() {
    let guide = guide();
    for (name, language, status) in [
        ("find-human", "text", 0),
        ("possible-json", "json", 0),
        ("reject-json", "json", 2),
        ("scan-human", "text", 0),
        ("matrix-human", "text", 0),
        ("scan-request", "json", 0),
        ("columns-request", "json", 0),
        ("budget-json", "json", 2),
    ] {
        let invocation = fence(&guide, &format!("{name}-command"), "bash").trim();
        assert_eq!(invocation.lines().count(), 1);
        let args = invocation
            .strip_prefix("./target/debug/bioferris ")
            .unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_bioferris"))
            .args(args.split_whitespace())
            .current_dir(root())
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(status),
            "{invocation}: {output:?}"
        );
        assert!(output.stderr.is_empty(), "{invocation}: {output:?}");
        assert_eq!(
            output.stdout,
            fence(&guide, &format!("{name}-output"), language).as_bytes(),
            "{invocation}"
        );
    }
}

#[test]
fn actual_cargo_alias_matches_all_documented_values_boundaries_and_errors() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "rust-command", "bash").trim(),
        "cargo bio motif_analysis"
    );
    let output = Command::new(env!("CARGO"))
        .args(["bio", "motif_analysis", "--offline", "--locked"])
        .current_dir(root())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert!(output.stderr.is_empty(), "{output:?}");
    assert_eq!(
        output.stdout,
        fence(&guide, "rust-output", "text").as_bytes()
    );
}
