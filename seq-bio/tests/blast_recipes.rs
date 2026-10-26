//! Check the documented BLAST inputs and programs, and run the real Cargo examples.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

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
fn complete_literal_inputs_and_program_match_their_files() {
    let guide = guide();
    for (name, file, language) in [
        ("input", "blast-example.tab6", "text"),
        ("commented-input", "blast-example.tab7", "text"),
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

#[test]
fn streaming_program_and_its_real_cargo_run_match_the_guide() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "stream-program", "rust"),
        fs::read_to_string(root().join("seq-bio/examples/blast_stream.rs")).unwrap()
    );
    assert_eq!(
        fence(&guide, "stream-command", "bash").trim(),
        "cargo bio blast_stream"
    );
    let observed = Command::new(env!("CARGO"))
        .args(["bio", "blast_stream", "--offline", "--locked"])
        .current_dir(root())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(observed.status.success(), "{observed:?}");
    assert!(observed.stderr.is_empty(), "{observed:?}");
    assert_eq!(
        observed.stdout,
        fence(&guide, "stream-output", "text").as_bytes()
    );
}
