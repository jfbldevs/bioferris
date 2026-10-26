//! Verify the actual Markdown program and its real `cargo bio` output.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn root() -> PathBuf {
    std::env::var_os("BIOFERRIS_WORKSPACE_ROOT").map_or_else(
        || {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .to_path_buf()
        },
        PathBuf::from,
    )
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/api/compressed-io-v1.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- compressed-io-{name} -->");
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
        fs::read_to_string(root().join("seq-bio/examples/compressed_io.rs")).unwrap()
    );
}

#[test]
fn actual_cargo_alias_matches_complete_documented_output() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "rust-command", "bash").trim(),
        "cargo bio compressed_io"
    );
    let output = Command::new(env!("CARGO"))
        .args(["bio", "compressed_io", "--offline", "--locked"])
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
