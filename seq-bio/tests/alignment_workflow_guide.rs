//! The combined alignment program documented in the Rust alignment recipes.

use std::path::Path;
use std::process::Command;

fn guide_block<'a>(guide: &'a str, name: &str) -> &'a str {
    guide
        .split_once(&format!("<!-- {name}:start -->"))
        .unwrap()
        .1
        .split_once(&format!("<!-- {name}:end -->"))
        .unwrap()
        .0
        .split_once('\n')
        .unwrap()
        .1
        .split_once('\n')
        .unwrap()
        .1
        .rsplit_once("```")
        .unwrap()
        .0
}

#[test]
fn manual_alignment_program_and_observed_stdout_match_the_guide() {
    let guide = include_str!("../../docs/examples/alignment-recipes.md");
    assert_eq!(
        guide_block(guide, "alignment-rust-program"),
        include_str!("../examples/alignment_workflow.rs")
    );
    let command: Vec<_> = guide_block(guide, "alignment-rust-command")
        .split_whitespace()
        .collect();
    let output = Command::new(command[0])
        .args(&command[1..])
        .current_dir(Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap())
        .output()
        .unwrap();
    assert!(output.status.success(), "{:?}", output.stderr);
    assert!(output.stderr.is_empty());
    assert_eq!(
        output.stdout,
        guide_block(guide, "alignment-rust-output").as_bytes()
    );
}
