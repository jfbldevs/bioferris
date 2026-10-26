# Contract and scientific assurance tests

## Contents

- [Capability and scope index](#capability-and-scope-index)
- [Run a complete Rust analysis](#run-a-complete-rust-analysis)
- [Independent scientific inventory](#independent-scientific-inventory)
- [Run the checks](#run-the-checks)
- [Failures, limits and compatibility](#failures-limits-and-compatibility)
- [Removed CLI conformance suite](#removed-cli-conformance-suite)

[API index](../README.md#api-index) · Testing policy ·
Existing integration evidence

## Capability and scope index

| Assurance | Rust access |
|---|---|
| Independent small-input alignment objectives | Public `Pairwise` facade, checked against exhaustive monotone paths |
| Biological orientation invariants | `find_orfs`, transcription and IUPAC motif searches |
| Ordered batch determinism | Pairwise batches at 1, 2 and 4 requested threads |

These tests add assurance, not biological functionality. The
[scientific assurance](../../seq-bio/tests/scientific_assurance.rs) source calls
the Rust library directly, with existing dependencies and the pinned Rust 1.85.0
toolchain. The suite that checked the removed `bioferris` CLI is described under
[Removed CLI conformance suite](#removed-cli-conformance-suite).

## Run a complete Rust analysis

Literal inputs are DNA `AC` and `AG`, with match 2, mismatch -1 and linear gap -2,
plus DNA `TTATTTCAT`. Its reverse complement is the complete CDS `ATGAAATAA`.
This complete program calls the library:

<!-- assurance-rust-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let first = dna("AC")?;
    let second = dna("AG")?;
    let alignment = Pairwise::needleman_wunsch().align(&first, &second)?;
    assert_eq!(alignment.score(), 1);
    assert_eq!(alignment.columns(), 2);
    println!("Global score: {}", alignment.score());
    println!("Trace columns: {}", alignment.columns());

    let source = dna("TTATTTCAT")?;
    let orfs = source.find_orfs()?;
    let orf = orfs.first().ok_or("missing expected reverse ORF")?;
    assert_eq!(orf.strand(), Strand::Reverse);
    assert_eq!(orf.peptide(), "MK");
    assert_eq!(orf.range(), 3..9);
    assert_eq!(orf.stop_range(), Some(0..3));
    println!("Reverse ORF: {} {:?} stop={:?}", orf.peptide(), orf.range(), orf.stop_range());
    Ok(())
}
```

From the repository root, compile and execute this Markdown program as a Rust doctest:

```bash
cargo build -p bioferris --lib --offline --locked --jobs 2
rustdoc --test docs/testing/contracts-and-science.md --edition 2024 --extern bioferris=target/debug/libbioferris.rlib -L dependency=target/debug/deps -D warnings --test-args=--nocapture
```

The program's complete stdout is below; its exit status is 0 and stderr is empty.
The Rustdoc test runner additionally prints its own test summary.

<!-- assurance-rust-output -->
```text
Global score: 1
Trace columns: 2
Reverse ORF: MK 3..9 stop=Some(0..3)
```

One match and one substitution give `2 + (-1) = 1`. ORF coordinates refer to
the original source, are zero-based and half-open, and exclude the separately
reported terminal stop. This biological example is not evidence of gene prediction
accuracy. The scientific test reads the output fence and recomputes these values.

## Independent scientific inventory

The alignment oracle enumerates monotone paths recursively. It does not call a
production verifier or reproduce a dynamic-programming table. Local cases
enumerate all possible starting vertices and allow termination at each vertex.
Every returned traceback is independently checked for bounds, symbol equality,
complete coordinate consumption, gap-run cost and optimal objective.

| Test family | Finite scope and independent checks |
|---|---|
| Linear global/local | All 225 ordered pairs of A/C words of lengths 0..3, three scoring schemes, two modes: 1,350 trace/score cases; score-only equality and deterministic traceback |
| Affine global/local | Same 225 pairs, three affine schemes, two modes: 1,350 cases; first gap symbol uses `open`, subsequent symbols `extension`; score-only equality and deterministic traceback |
| Edit/WFA/adaptive | 225 unit-edit cases plus 675 WFA and 675 adaptive-affine cases: 1,575 cases; WFA's `open + length*extension` is explicitly converted into the oracle's first-symbol cost |
| ORF orientation | 25 independently known two-codon products at three source offsets: 75 cases; reverse-complement coordinate/strand bijection, retained terminal-stop coordinates and DNA/RNA parity |
| IUPAC orientation | All 240 nonempty IUPAC patterns of lengths 1..2 against a fixed ambiguous source, definite/possible policies: 480 cases; range/strand bijection preserves overlap and palindromic labels |
| Ordered batch determinism | 32 ordered A/C input pairs; Needleman-Wunsch, Smith-Waterman and global Gotoh produce identical complete results at 1, 2 and 4 requested threads; 96 baseline traces are also independently rescored and checked for optimality |

The 4,275 alignment cases use exact integer objectives. They establish optimality
and valid deterministic traces in this finite scope, not equality to an external
implementation's choice among tied optimal traces. Biological metamorphic
relations complement the existing exhaustive codon tables and retained exact
numerical observations; they do not replace those tests or introduce tolerances.

## Run the checks

From the repository root, after the pinned toolchain and locked dependencies are
available locally:

```bash
cargo test -p bioferris --test scientific_assurance --offline --locked --jobs 2
cargo clippy -p bioferris --test scientific_assurance --offline --locked --jobs 2 -- -D warnings
rustfmt --edition 2024 --check seq-bio/tests/scientific_assurance.rs
```

The first command runs only the scientific target, not the workspace/native-MSA
campaign. A successful test process returns exit 0; assertion/build failures
return nonzero. Observed test-binary stdout on 2026-10-21, exit 0 (Cargo's build
progress goes to stderr); the Clippy and rustfmt commands also returned exit 0:

```text
running 7 tests
test published_scientific_example_output_matches_literal_inputs ... ok
test ordered_pairwise_batches_are_trace_exact_at_one_two_and_four_threads ... ok
test orf_coordinates_and_products_commute_with_whole_source_reverse_complement ... ok
test iupac_hits_obey_reverse_complement_bijection_without_collapsing_strands ... ok
test exhaustive_linear_paths_verify_scores_tracebacks_and_score_only_results ... ok
test exhaustive_affine_paths_verify_global_and_local_gap_run_costs ... ok
test exhaustive_edit_and_wavefront_objectives_match_independent_path_enumeration ... ok

test result: ok. 7 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s
```

Test order and timing can vary between runs.

## Failures, limits and compatibility

The scientific tests run in memory; they write no files. On failure, retain the
test name, scoring scheme and input pair. Investigate the implementation or contract discrepancy; do not refresh
expected output automatically, weaken exact comparisons or silently omit cases.
These tests do not qualify untested platforms, all numeric values, broader format
dialects, filesystem crash durability, gene prediction or universal Biopython/native
equivalence. Historical reference corpora and closed scientific evidence remain
unchanged.

The separate QA-runner adversaries
also executed the actual `bioferris-qa` binary in isolated Git fixtures:
`cargo test -p bioferris-qa --test runner_adversarial --offline --locked --jobs 2 -- --nocapture`
returned exit 0 with 7 passed and no failed, ignored, measured or filtered tests
(six adversarial scenarios plus the ordinarily inert Rust self-worker).
A tracked `sealed_module.rs` changed from `const INITIAL: u8 = 1;` to
`const CHANGED: u8 = 2;`: the child really exited 0, but source-hash mismatch
prevented overall PASS. Other cases observed failed-and-cleaned inherited
descendants, INCOMPLETE receipts after QA-parent interruption, actual SIGTERM
and timeout SIGKILL statuses, and OUTPUT_LIMIT for a strictly bounded writer of
64 MiB plus 8,192 bytes. Parent interruption does not imply automatic child
cleanup: the test retains and validates the owned worker's PID/process group,
executable and start identity before terminating it. These results were observed
on the recorded Linux x86_64 host; non-Linux Unix fallback and Windows execution
remain unexecuted locally, and process groups are not a secure sandbox.

## Removed CLI conformance suite

Until the `bioferris` CLI was removed (decision 029),
`seq-bio/tests/discovery_conformance.rs` derived a 38-operation whitelist from the
CLI's `--help --json` and `formats --json` discovery and required every operation
to belong to exactly one row of the requirements matrix. It
checked 1,565 foreign operation/field pairs, required-field omissions, raw
duplicates, 259 wrong-type cases and 73 optional-null cases in the strict JSON
requests, ran this page's former `request` examples under five-second deadlines
and one-MiB output caps, and compared direct CLI, strict request and public Rust
results for ten routes. Those adapter checks guarded an interface that no longer
exists, so no replacement suite runs; the Rust cores they reached keep their own
tests. The suite and its examples remain in product commit `e9b3e644`, and their
results remain in the validation record.

The requirements matrix is now schema version 2 without CLI operations.
`requirements_traceability.rs`
checks that its 23 rows have unique IDs, the expected fields, nonempty
capability/acceptance/scope text and existing test files:

```bash
cargo test -p bioferris --test requirements_traceability --offline --locked
```

A row naming a missing test file fails with `<id>: missing <path>`.
