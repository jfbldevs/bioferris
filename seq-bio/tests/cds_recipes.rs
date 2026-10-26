//! Check the documented annotated-CDS program and run it through the real Cargo alias.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/cds-translation-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- cds-{name} -->");
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
fn literal_inputs_and_complete_program_match_their_checked_in_files() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "input", "text"),
        include_str!("fixtures/annotated-example.gb")
    );
    let program = include_str!("../examples/translate_annotations.rs");
    assert_eq!(fence(&guide, "rust-program", "rust"), program);
    let input = program
        .split_once("const INPUT: &str = r#\"")
        .unwrap()
        .1
        .split_once("\"#;")
        .unwrap()
        .0;
    assert_eq!(input, fence(&guide, "input", "text"));
}

#[test]
fn actual_cargo_alias_matches_the_complete_documented_analysis() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "rust-command", "bash").trim(),
        "cargo bio translate_annotations"
    );
    let output = Command::new(env!("CARGO"))
        .args(["bio", "translate_annotations", "--offline", "--locked"])
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
