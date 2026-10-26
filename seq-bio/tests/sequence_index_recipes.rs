//! Run the documented indexing programs and verify their full inputs and outputs.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::io::index::{IndexOptions, SequenceIndex, SourceCheck};
use bioferris::prelude::Dna;

const FASTQ: &[u8] = include_bytes!("fixtures/index-example.fastq");
const FIXTURE_PATH: &str = "seq-bio/tests/fixtures/index-example.fastq";

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/examples/sequence-index-recipes.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- index-recipe-{name} -->");
    assert_eq!(guide.matches(&marker).count(), 1, "unique {marker}");
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
        for _ in 0..128 {
            let path = std::env::temp_dir().join(format!(
                "bioferris-index-recipes-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => {
                    let directory = Self(path);
                    let fixtures = directory.0.join("seq-bio/tests/fixtures");
                    fs::create_dir_all(&fixtures).unwrap();
                    fs::write(fixtures.join("index-example.fastq"), FASTQ).unwrap();
                    return directory;
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(error) => panic!("cannot create recipe directory: {error}"),
            }
        }
        panic!("cannot allocate a unique recipe directory");
    }

    fn verify_unchanged_inputs(&self) {
        assert_eq!(fs::read(self.0.join(FIXTURE_PATH)).unwrap(), FASTQ);
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn assert_output(output: &Output, expected: &str) {
    assert_eq!(
        output.status.code(),
        Some(0),
        "stdout: {}; stderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty(), "{:?}", output.stderr);
    assert_eq!(output.stdout, expected.as_bytes());
}

/// Runs a documented one-line `cargo bio EXAMPLE -- PATHS` command with its file
/// arguments rebased to an isolated copy of the fixture.
fn run_documented(guide: &str, marker: &str, example: &str, expected_paths: &[&str]) -> Directory {
    let directory = Directory::new();
    let command = fence(guide, marker, "bash").trim();
    assert_eq!(command.lines().count(), 1);
    let (invocation, paths) = command.split_once(" -- ").unwrap();
    let mut arguments = invocation.split_whitespace();
    assert_eq!(arguments.next(), Some("cargo"));
    assert_eq!(arguments.clone().collect::<Vec<_>>(), ["bio", example]);
    let paths: Vec<_> = paths.split_whitespace().collect();
    assert_eq!(paths, expected_paths);
    // Cargo resolves the checked-in alias from the real root, while the file
    // arguments point to this test's isolated copy rather than user artifacts.
    let output = Command::new(env!("CARGO"))
        .args(arguments)
        .args(["--offline", "--locked", "--"])
        .args(paths.iter().map(|path| directory.0.join(path)))
        .current_dir(root())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let expected = marker.replace("command", "output");
    assert_output(&output, fence(guide, &expected, "text"));
    directory
}

#[test]
fn documented_literal_inputs_and_complete_programs_match_the_source_files() {
    let guide = guide();
    assert_eq!(fence(&guide, "fastq-input", "fastq").as_bytes(), FASTQ);
    assert_eq!(
        fence(&guide, "rust-program", "rust,no_run"),
        include_str!("../examples/indexed_records.rs")
    );
    assert_eq!(
        fence(&guide, "binary-program", "rust,no_run"),
        include_str!("../examples/binary_index.rs")
    );
}

#[test]
fn manual_recipe_runs_the_actual_cargo_alias_with_isolated_file_arguments() {
    let guide = guide();
    let directory = run_documented(
        &guide,
        "rust-command",
        "indexed_records",
        &[FIXTURE_PATH, "rust-reads.index.json"],
    );
    let mut index = SequenceIndex::<Dna>::open(
        directory.0.join(FIXTURE_PATH),
        directory.0.join("rust-reads.index.json"),
        IndexOptions::default(),
    )
    .unwrap();
    assert_eq!(index.len(), 2);
    assert!(index.get("absent").unwrap().is_none());
    directory.verify_unchanged_inputs();
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 2);
}

#[test]
fn binary_recipe_runs_the_actual_cargo_alias_with_isolated_file_arguments() {
    let guide = guide();
    let directory = run_documented(
        &guide,
        "binary-command",
        "binary_index",
        &[FIXTURE_PATH, "rust-reads.index.json", "rust-reads.bfidx"],
    );
    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut index = SequenceIndex::<Dna>::open_with(
            directory.0.join(FIXTURE_PATH),
            directory.0.join("rust-reads.bfidx"),
            IndexOptions::default(),
            check,
        )
        .unwrap();
        assert_eq!(index.ids().collect::<Vec<_>>(), ["read1", "read2"]);
        assert!(index.get("absent").unwrap().is_none());
    }
    directory.verify_unchanged_inputs();
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 3);
}
