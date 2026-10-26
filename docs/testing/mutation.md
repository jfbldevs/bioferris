# Curated semantic mutation audit

## Contents

- [Scope and access](#scope-and-access)
- [Literal faults and targeted tests](#literal-faults-and-targeted-tests)
- [Run the audit](#run-the-audit)
- [Evidence and classification](#evidence-and-classification)
- [Source and process safety](#source-and-process-safety)
- [Limits and verification](#limits-and-verification)

[Documentation index](../README.md) · Audit implementation

## Scope and access

This repository-local Rust audit deliberately breaks five independently tested
scientific or safety rules **only in new sandbox copies**. It checks whether the
existing tests actually reject those faults. It is not a new biological feature,
a biological API operation or a comprehensive mutation-testing tool.
No external mutation package, runtime installation or new dependency is required.

| Access | Interface |
|---|---|
| Maintainer command | `bioferris-qa mutation --output NEW_DIRECTORY` |
| Internal Rust helper | `mutation::run(&Path) -> Result<bool, Box<dyn Error>>` |
| Agent evidence | Outer durable receipts plus `audit/summary.json`, source hashes, exact patches and logs |

The public QA wrapper captures a real helper exit status and enforces its outer
deadline. The internal helper is not an API replacement for that wrapper.
A missing final summary is **incomplete**, never a passing audit.

## Literal faults and targeted tests

Each replacement must have exactly one occurrence in the hashed source copy.
Zero/multiple matches or an unchanged replacement abort instead of silently
testing another fault. Every mutant starts from the same unmodified snapshot;
faults are never accumulated.

| ID | File under workspace | Exact original text | Exact replacement |
|---|---|---|---|
| `bgzf-crc-bypass` | `seq-bio/src/io/bgzf.rs` | `if actual_crc != expected_crc {` | `if false && actual_crc != expected_crc {` |
| `index-span-digest-bypass` | `seq-bio/src/io/index.rs` | `\|\| hashing.digest.finalize().as_slice() != location.sha256` | `\|\| false && hashing.digest.finalize().as_slice() != location.sha256` |
| `motif-exclusive-threshold` | `seq-bio/src/motifs.rs` | `if score >= threshold {` | `if score > threshold {` |
| `blast-coordinate-off-by-one` | `seq-bio/src/io/blast.rs` | `self.reported_start.min(self.reported_end) - 1..self.reported_start.max(self.reported_end)` | `self.reported_start.min(self.reported_end)..self.reported_start.max(self.reported_end)` |
| `alignment-missing-counted-as-gap` | `seq-bio/src/alignment_analysis.rs` | `AlignmentCell::Missing(_) => self.missing = add(self.missing, 1)?,` | `AlignmentCell::Missing(_) => self.gaps = add(self.gaps, 1)?,` |

The index row's text changed on 2026-10-22, when record digests became raw
bytes for [binary sidecars](../algorithms/binary-sequence-index-v1.md); the
mutation and its targeted test are the same. The
recorded audit
used the earlier text,
`|| format!("{:x}", hashing.digest.finalize()) != entry.sha256`.

The corresponding input scenarios and assertions are:

- BGZF: an independently constructed `abc` member has its declared CRC altered;
  the test expects typed CRC failure with the original actual/declared values.
- Index: stable plain source `>a\nACGT\n` is first read with its correct selected
  span digest and must return exactly record `a`, sequence `ACGT`. Only the
  expected SHA-256 is then replaced by 64 zeroes. The real `read_entry` path must
  return typed `StaleSource` with its span-integrity message while source bytes
  and metadata still match. This isolates the digest layer: the existing
  restored-mtime source-edit integration test can also be rejected by Unix ctime
  guards, so it cannot by itself prove that the digest check remains effective.
- Motif: one-column scores `[0.1; 4]` scan DNA `AA`, with a double threshold
  infinitesimally above the rounded float32 score but rounding to that same value.
  The inclusive boundary must yield two hits, not zero.
- BLAST: the retained public NCBI cookbook tabular fixture's reported coordinates
  must convert to the literal zero-based half-open ranges `5564..5671` and
  `6820..6928`, preserving original reverse endpoints.
- Alignment: rows `A-N?`, `A--A`, `C--A`, `N.-A` must count eight residues,
  seven gaps and one missing cell, with the independently counted denominators.

All five unmodified targeted commands must pass first, in the snapshot. These
are the actual argument forms; `--exact` prevents a broad substring selection:

```sh
cargo test --offline --locked -p bioferris --lib io::bgzf::tests::crc_and_isize_errors_retain_actual_and_declared_values --jobs 2 -- --exact --nocapture
cargo test --offline --locked -p bioferris --lib io::index::tests::selected_span_digest_is_checked_independently_of_source_metadata --jobs 2 -- --exact --nocapture
cargo test --offline --locked -p bioferris --lib motifs::tests::thresholds_round_to_scan_precision_before_inclusive_comparison --jobs 2 -- --exact --nocapture
cargo test --offline --locked -p bioferris --lib io::blast::tests::ncbi_cookbook_rows_preserve_scores_ids_and_reverse_endpoints --jobs 2 -- --exact --nocapture
cargo test --offline --locked -p bioferris --lib alignment_analysis::tests::exact_statistics_cover_missing_ambiguity_gaps_and_all_denominators --jobs 2 -- --exact --nocapture
```

Each mutant reruns only its corresponding command in its independent new source
tree. These are ordinary Rust unit/integration tests, not a documentation
navigation check or a scientific qualification/campaign rerun.

## Run the audit

Prerequisites: repository root, pinned Rust 1.85.0, already cached locked
dependencies and a trusted writable destination parent. The audit is offline;
a missing dependency is an error, not permission to install or change packages.

```sh
qa_parent="$(mktemp -d /tmp/bioferris-mutation-example.XXXXXX)"
cargo run --locked --quiet -p bioferris-qa --jobs 2 -- mutation --output "$qa_parent/run"
```

`qa_parent/run` must not already exist, even as an empty directory. The outer
runner owns that new directory; the helper creates its new `audit` child.
All artifacts remain for inspection. There is no implicit deletion of evidence,
user sources, original files or previous audit runs.

Verification results for this increment are recorded below after actual audit
execution; no expected mutant outcome is substituted for a real test failure.

## Evidence and classification

The outer receipt directory preserves its command plan, Git/toolchain/source
provenance, actual exit/signal, stdout/stderr and terminal status. Its `audit`
child contains:

- `started.json`: the five literal faults, scope and start timestamp.
- `source-sha256.json`: every explicitly selected source/manifest/fixture hash.
- `baseline-source/`: the unmodified reproducible source snapshot used by baseline
  commands; `mutants/ID/`: five fresh independently patched copies.
- `build/`: one shared isolated Cargo target directory, reused sequentially.
- `logs/NAME.started.json`: exact argv, working directory, pinned toolchain,
  target-directory environment and 600-second per-command deadline.
- `logs/NAME.stdout.log`, `.stderr.log`, `.result.json`: complete subprocess
  streams, actual exit/signal, elapsed time, hashes and executed test names.
- `logs/ID.patch.json`: exact before/after text, one-match proof and both file
  hashes. Paths are relative to the copied workspace.
- `summary.json`: baseline executions, all mutant states, killed/expected counts,
  elapsed time and the original-source unchanged check.

| State | Required evidence / meaning |
|---|---|
| Baseline `PASS` | Actual exit `0`, selected test reported `ok`, one passing test, no failed/ignored/measured tests |
| `KILLED` | Actual exit `101`, selected test `FAILED`, one failing runtime test, no passed/ignored/measured tests, and a panic diagnostic from an assertion or an unexpected `Ok` at `unwrap_err()` |
| `SURVIVED` | The mutated selected test genuinely passed; the audit fails |
| `INVALID_MUTANT` | Compiler diagnostic / compilation failure; not a kill |
| `TIMEOUT` | Command exceeded 600 seconds; terminated and reaped; not a kill |
| `EXECUTION_FAILED` | Spawn/runtime/signal/log evidence does not establish the selected expectation failure or pass; not a kill |
| Missing summary | Interrupted, aborted or incomplete; not a pass |

The `unwrap_err()` case is a runtime test expectation, not a compilation error.
Arbitrary panics are deliberately not enough to count as a killed mutant. A
zero-test successful command also cannot pass. Filtered unrelated test counts
are expected because the selectors intentionally execute exactly one test.

An audit returns success only if the baseline passed all five selectors, all five
mutants were killed by those runtime expectations and original selected source
hashes still match. Nested mutant Cargo exit `101` is expected evidence of a kill;
it is retained, not rewritten to zero. Only the aggregate QA wrapper succeeds
when those reviewed expectations are satisfied.

## Source and process safety

The copy allowlist is explicit: workspace member `Cargo.toml` and `src` trees,
root Cargo manifest/lock/toolchain/local alias configuration and thirteen compile-required public
biological fixtures. The explicitly declared pairwise example and benchmark
target source directories are also copied so Cargo can resolve their manifests;
none of those examples or benchmarks is executed. Source-tree static calibration data used by the existing
crates is copied unchanged. Symlinks and special source files are refused.
Every copied file is checked against its SHA-256 before use; original selected
files are re-enumerated and rehashed at completion.

No root data tree, scripts, documentation evidence archive, sealed/reserved MSA
manifest or frozen qualification collection is traversed. No scientific kernel
or default is changed in the real worktree. Snapshot sources include workspace
members needed by Cargo resolution; this does not authorize their campaigns or
execute all their tests.

Commands run with `RUSTUP_TOOLCHAIN=1.85.0`, `--offline --locked --jobs 2`,
`RUST_TEST_THREADS=1`, no interactive stdin and file-backed output logs. Each
Cargo command has its own deadline. Unix children get a dedicated process group
through safe Rust's `CommandExt`; the OS `kill` utility terminates that owned
group on timeout. Windows uses `taskkill /T`; unsupported platforms can terminate
the direct child only and do not gain a portable process-tree guarantee. There
is no unsafe Rust or direct FFI.

Concurrent edits to selected real source during an audit are correctly treated
as source drift, even if unrelated to the injected fault. Run after source edits
have stopped; documentation-only work outside the copied allowlist can continue.
An interrupted copy or missing dependency leaves an incomplete/failing audit,
not permission to infer a result from an old log.

## Limits and verification

Five transparent cases provide specific evidence that important CRC, selected
span integrity, numerical threshold, coordinate and missing-data assertions are
effective. They do **not** measure a representative automatic mutation score,
exhaust all branches or prove universal scientific correctness. Equivalent
mutants, broad operator generation, higher-order faults and systematic mutation
of all algorithms remain outside this increment.

The harness unit tests independently reject empty/duplicate/no-op patch matches,
zero-test success, unrelated runtime failures, compile failure and timeout as
false kill evidence. Whole-workspace QA and historical scientific acceptance
remain separate; this audit neither modifies nor expands sealed MSA qualification.
