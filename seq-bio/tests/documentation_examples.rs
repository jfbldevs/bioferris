//! Keep the quick-start and identity programs, inputs, outputs and exit codes executable.

use std::error::Error;
use std::path::Path;
use std::process::Command;

#[test]
fn identity_guide_matches_the_verified_example() {
    let guide = include_str!("../../docs/api/bioferris-identity.md");
    let program = include_str!("../examples/my_analysis.rs");
    let (_, program) = program.split_once("\n\n").unwrap();
    let program_block = guide.split_once("<!-- identity-program -->").unwrap().1;
    assert_eq!(fence(program_block, "rust").trim(), program.trim());
    let output_block = guide.split_once("<!-- identity-output -->").unwrap().1;
    let quick_start = include_str!("../../docs/quick-start.md");
    assert_eq!(
        fence(output_block, "text"),
        fence(marked(quick_start, "quick-start-sequence-output"), "text")
    );
}

#[test]
fn quick_start_uses_the_actual_programs_and_fasta_input() {
    let guide = include_str!("../../docs/quick-start.md");
    for (name, program) in [
        (
            "quick-start-sequence-program",
            include_str!("../examples/my_analysis.rs"),
        ),
        (
            "quick-start-alignment-program",
            include_str!("../examples/quick_start.rs"),
        ),
        (
            "quick-start-invalid-program",
            include_str!("../examples/invalid_input.rs"),
        ),
    ] {
        let (_, program) = program
            .split_once("\n\n")
            .expect("example module documentation");
        assert_eq!(fence(marked(guide, name), "rust").trim(), program.trim());
    }

    let (_, input) = guide
        .split_once("<!-- quick-start-fasta-input -->")
        .expect("documented FASTA input");
    let (input, _) = input
        .split_once("<!-- /quick-start-fasta-input -->")
        .expect("end of documented FASTA input");
    assert_eq!(
        fence(input, "fasta"),
        include_str!("fixtures/phix174-e.fasta")
    );
}

fn marked<'a>(body: &'a str, name: &str) -> &'a str {
    let start = format!("<!-- {name} -->");
    let end = format!("<!-- /{name} -->");
    assert_eq!(body.matches(&start).count(), 1, "unique start: {name}");
    assert_eq!(body.matches(&end).count(), 1, "unique end: {name}");
    let (_, rest) = body.split_once(&start).expect("missing start marker");
    let (content, _) = rest.split_once(&end).expect("missing end marker");
    content
}

#[test]
fn cargo_alias_matches_documented_success_and_error() -> Result<(), Box<dyn Error>> {
    let guide = include_str!("../../docs/quick-start.md");
    assert_eq!(
        fence(marked(guide, "quick-start-cargo-alias"), "toml"),
        include_str!("../../.cargo/config.toml")
    );

    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root");
    for (name, documented, exit) in [
        ("quick-start-sequence", &["bio", "my_analysis"][..], 0),
        ("quick-start-alignment", &["bio", "quick_start"], 0),
        ("quick-start-invalid", &["bio", "invalid_input"], 0),
        (
            "quick-start-invalid-dna",
            &["bio", "invalid_input", "--", "ACGU"],
            1,
        ),
        (
            "quick-start-unknown-example",
            &["bio", "nonexistent_analysis"],
            101,
        ),
    ] {
        let command = fence(marked(guide, &format!("{name}-command")), "bash");
        assert_eq!(command.lines().count(), 1);
        let mut args = command.split_whitespace();
        assert_eq!(args.next(), Some("cargo"));
        assert_eq!(args.collect::<Vec<_>>(), documented);
        // Cargo options must precede the `--` that starts program arguments.
        let (cargo_args, program_args) = documented
            .iter()
            .position(|arg| *arg == "--")
            .map_or((documented, &[][..]), |split| documented.split_at(split));

        // Run the real checked-in alias, not a stale example binary. This invokes
        // `run`, not `test`, so it does not recursively execute this test suite.
        let output = Command::new(env!("CARGO"))
            .args(cargo_args)
            .args(["--offline", "--locked"])
            .args(program_args)
            .current_dir(root)
            .output()?;
        assert_eq!(
            output.status.code(),
            Some(exit),
            "{name}: unexpected exit status: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let (actual, unused) = if exit == 0 {
            (&output.stdout, &output.stderr)
        } else {
            (&output.stderr, &output.stdout)
        };
        assert!(unused.is_empty(), "{name}: unexpected data in other stream");
        assert_eq!(
            actual,
            fence(marked(guide, &format!("{name}-output")), "text").as_bytes()
        );
    }
    Ok(())
}

fn fence<'a>(body: &'a str, language: &str) -> &'a str {
    let marker = format!("```{language}\n");
    let (_, rest) = body.split_once(&marker).expect("missing recipe code fence");
    let (content, _) = rest.split_once("```").expect("unclosed recipe code fence");
    content
}
