//! Execute complete annotated-output documentation through the real Cargo alias.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::Dna;
use bioferris::io::annotated::parse_genbank;

const GENBANK: &[u8] = include_bytes!("fixtures/annotated-example.gb");
const EMBL: &[u8] = include_bytes!("fixtures/annotated-example.embl");

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()
}

fn guide() -> String {
    fs::read_to_string(root().join("docs/examples/annotated-output-recipes.md")).unwrap()
}

fn fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    let marker = format!("<!-- annotated-output-recipe-{name} -->");
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
                "bioferris-annotated-recipes-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(error) => panic!("cannot create recipe directory: {error}"),
            }
        }
        panic!("cannot allocate a unique recipe directory");
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn run(guide: &str, name: &str, destination: Option<&Path>) -> Output {
    let text = fence(guide, &format!("{name}-command"), "bash").trim();
    assert_eq!(text.lines().count(), 1);
    let (invocation, file_argument) = text
        .split_once(" -- ")
        .map_or((text, None), |(left, right)| (left, Some(right)));
    let mut arguments = invocation.split_whitespace();
    assert_eq!(arguments.next(), Some("cargo"));
    assert_eq!(
        arguments.clone().collect::<Vec<_>>(),
        ["bio", "annotated_records"]
    );
    let mut command = Command::new(env!("CARGO"));
    command.args(arguments).args(["--offline", "--locked"]);
    match (file_argument, destination) {
        (Some("annotated-slice.gb"), Some(path)) => {
            command.arg("--").arg(path);
        }
        (None, None) => {}
        _ => panic!("unexpected documented output argument"),
    }
    command
        .current_dir(root())
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn assert_result(output: &Output, status: i32, stdout: &[u8], stderr: &[u8]) {
    assert_eq!(
        output.status.code(),
        Some(status),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, stdout);
    assert_eq!(output.stderr, stderr);
    assert_eq!(
        fs::read(root().join("seq-bio/tests/fixtures/annotated-example.gb")).unwrap(),
        GENBANK
    );
    assert_eq!(
        fs::read(root().join("seq-bio/tests/fixtures/annotated-example.embl")).unwrap(),
        EMBL
    );
}

#[test]
fn literal_inputs_and_entire_example_source_match_the_recipe() {
    let guide = guide();
    assert_eq!(
        fence(&guide, "genbank-input", "genbank").as_bytes(),
        GENBANK
    );
    assert_eq!(fence(&guide, "embl-input", "embl").as_bytes(), EMBL);
    assert_eq!(
        fence(&guide, "program", "rust,no_run"),
        include_str!("../examples/annotated_records.rs")
    );
}

#[test]
fn actual_manual_alias_matches_the_complete_transformation_output() {
    let guide = guide();
    assert_result(
        &run(&guide, "manual", None),
        0,
        fence(&guide, "manual-output", "text").as_bytes(),
        b"",
    );
}

#[test]
fn named_output_is_exact_and_existing_destination_is_preserved() {
    let guide = guide();
    let directory = Directory::new();
    let destination = directory.0.join("annotated-slice.gb");
    let expected_stdout = format!(
        "{}{}",
        fence(&guide, "manual-output", "text"),
        fence(&guide, "save-suffix", "text")
    );
    assert_result(
        &run(&guide, "save", Some(&destination)),
        0,
        expected_stdout.as_bytes(),
        b"",
    );
    let bytes = fs::read(&destination).unwrap();
    assert_eq!(bytes, fence(&guide, "saved-file", "genbank").as_bytes());
    let parsed = parse_genbank::<Dna>(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(parsed.len(), 1);
    assert_eq!(parsed[0].record().id(), "LOCAL_cds");
    assert_eq!(parsed[0].record().source_sequence(), "ATGGCTTAA");
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
    assert_eq!(
        fence(&guide, "save-command", "bash"),
        fence(&guide, "existing-command", "bash")
    );
    assert_result(
        &run(&guide, "existing", Some(&destination)),
        1,
        fence(&guide, "manual-output", "text").as_bytes(),
        fence(&guide, "existing-stderr", "text").as_bytes(),
    );
    assert_eq!(fs::read(&destination).unwrap(), bytes);
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
}
