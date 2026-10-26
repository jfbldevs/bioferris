//! Keep the complete protein-property program, its invocation and observed result together.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/protein-properties-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- protein-properties-{name} -->");
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
fn literal_program_matches_its_complete_documented_source() {
    assert_eq!(
        fence(&guide(), "rust-program", "rust"),
        include_str!("../examples/protein_analysis.rs")
    );
}

#[test]
fn actual_cargo_alias_matches_all_documented_values_boundaries_and_errors() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "rust-command", "bash").trim(),
        "cargo bio protein_analysis"
    );
    let output = Command::new(env!("CARGO"))
        .args(["bio", "protein_analysis", "--offline", "--locked"])
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
