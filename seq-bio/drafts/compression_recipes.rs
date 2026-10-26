//! Verify the actual Markdown program, literal files and complete CLI outputs.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

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

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-compression-recipe-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn prepare(&self, guide: &str) {
        fs::write(self.0.join("input.fa"), fence(guide, "fasta-input", "text")).unwrap();
        fs::write(
            self.0.join("index-request.json"),
            fence(guide, "request-input", "json"),
        )
        .unwrap();
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn execute(
    guide: &str,
    directory: &Directory,
    command: &str,
    result: &str,
    language: &str,
    status: i32,
) {
    let invocation = fence(guide, &format!("{command}-command"), "bash").trim();
    assert_eq!(invocation.lines().count(), 1);
    let args = invocation.strip_prefix("bioferris ").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_bioferris"))
        .args(args.split_whitespace())
        .current_dir(&directory.0)
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
        fence(guide, &format!("{result}-output"), language).as_bytes(),
        "{invocation}"
    );
}

#[test]
fn literal_program_and_inputs_match_their_complete_documented_fences() {
    let guide = guide();
    let promoted = root().join("seq-bio/examples/compressed_io.rs");
    let source = if promoted.exists() {
        promoted
    } else {
        root().join("seq-bio/drafts/compressed_io.rs")
    };
    assert_eq!(
        fence(&guide, "rust-program", "rust"),
        fs::read_to_string(source).unwrap()
    );
    assert_eq!(
        fence(&guide, "fasta-input", "text"),
        fs::read_to_string(root().join("seq-bio/tests/fixtures/compression-example.fa")).unwrap()
    );
    assert_eq!(
        fence(&guide, "request-input", "json"),
        fs::read_to_string(root().join("seq-bio/tests/fixtures/compression-index-request.json"))
            .unwrap()
    );
}

#[test]
fn all_nine_direct_cli_recipes_match_byte_exact_output_status_and_streams() {
    let guide = guide();
    let directory = Directory::new();
    directory.prepare(&guide);
    for (name, language, status) in [
        ("gzip", "text", 0),
        ("bgzf", "json", 0),
        ("index", "json", 0),
        ("record-human", "text", 0),
        ("record-json", "json", 0),
        ("missing", "json", 2),
        ("occupied", "json", 2),
        ("limit", "json", 2),
        ("ordinary", "json", 2),
    ] {
        execute(&guide, &directory, name, name, language, status);
    }
    assert!(!directory.0.join("limit.gz").exists());
    assert!(!directory.0.join("invalid.index.json").exists());
}

#[test]
fn actual_request_matches_the_same_typed_index_result_without_reparsing_strings() {
    let guide = guide();
    let directory = Directory::new();
    directory.prepare(&guide);
    execute(&guide, &directory, "bgzf", "bgzf", "json", 0);
    execute(&guide, &directory, "request", "index", "json", 0);
    execute(&guide, &directory, "record-json", "record-json", "json", 0);
}

#[test]
fn nested_compression_parameters_reject_arrays_duplicates_unknown_and_untyped_numbers() {
    let directory = Directory::new();
    for parameters in [
        "[]",
        "[6]",
        r#"{"level":6,"level":6}"#,
        r#"{"unknown":1}"#,
        r#"{"level":-1}"#,
        r#"{"level":6.5}"#,
        r#"{"max_input_bytes":"1"}"#,
    ] {
        let request = format!(
            r#"{{"schema_version":1,"operation":"compress-file","compression":"gzip","input":"input.fa","output":"never.gz","compression_parameters":{parameters}}}"#
        );
        fs::write(directory.0.join("request.json"), request).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_bioferris"))
            .args(["request", "--input", "request.json", "--json"])
            .current_dir(&directory.0)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stderr.is_empty());
        let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(
            value["error"]["code"], "invalid_schema",
            "{parameters}: {value}"
        );
        assert!(!directory.0.join("never.gz").exists());
    }
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
