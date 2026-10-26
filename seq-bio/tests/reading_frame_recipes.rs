//! Keep the reading-frame guide's complete program and output executable.

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/reading-frames-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- reading-frame-{name} -->");
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
fn literal_program_matches_its_documentation() {
    assert_eq!(
        fence(&guide(), "rust-program", "rust"),
        include_str!("../examples/reading_frames.rs")
    );
}

#[test]
fn actual_cargo_alias_matches_the_complete_manual_output() {
    let guide = guide();
    let invocation = fence(&guide, "rust-command", "bash").trim();
    assert_eq!(invocation, "cargo bio reading_frames");
    let output = Command::new(env!("CARGO"))
        .args(["bio", "reading_frames", "--offline", "--locked"])
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
