//! Verify the complete documented motif program and its real `cargo bio` output.

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
fn literal_program_matches_its_complete_documented_source() {
    let program = fs::read_to_string(root().join("seq-bio/examples/motif_analysis.rs")).unwrap();
    assert_eq!(fence(&guide(), "rust-program", "rust"), program);
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
