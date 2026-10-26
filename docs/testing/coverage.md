# Scoped Rust source coverage

## Contents

- [Access and status](#access-and-status)
- [Inputs and invocation](#inputs-and-invocation)
- [Measurement pipeline](#measurement-pipeline)
- [Artifacts and result interpretation](#artifacts-and-result-interpretation)
- [Errors and limits](#errors-and-limits)
- [Verification and boundaries](#verification-and-boundaries)

## Access and status

| Access | Availability |
|---|---|
| Repository QA command | `bioferris-qa coverage --output NEW_DIRECTORY`; integrated run verified on Linux x86_64 on 2026-10-20 |
| Biological Rust API | No new biological method; ordinary public API tests are measured |

[API index](../README.md) · [Robustness corpus](robustness.md) ·
Source

This diagnostic measures a declared ordinary regression scope with the pinned
Rust 1.85.0 LLVM tools. It does not imply whole-workspace coverage, native
reference fidelity or renewed scientific qualification. No `cargo-llvm-cov`
plugin, new crate, compiler upgrade or external-language runtime is required.

## Inputs and invocation

The input is the current repository source, lockfile and test inventory, with
these six explicit Cargo packages:

```text
bioferris
align-core
align-pairwise
align-fast
align-fast-kernels
bioferris-qa
```

Run from the repository root. `target/qa-coverage-v5` must not already exist:

```bash
LANG=C.UTF-8 LC_ALL=C.UTF-8 CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 cargo +1.85.0 run -p bioferris-qa --offline --locked --jobs 2 -- coverage --output target/qa-coverage-v5
```

The environment prefix matches the CI resource policy (no incremental state or
dev/test debug symbols); the receipts record it.

The `llvm-tools-preview` component must already be installed for Rust 1.85.0.
The workflow locates `llvm-cov` and `llvm-profdata` inside that exact toolchain;
it never installs components, downloads plugins or chooses system LLVM tools.
Unset ambient `RUSTFLAGS`, `CARGO_ENCODED_RUSTFLAGS` and `CARGO_BUILD_TARGET`
before this workflow; conflicting instrumentation/target settings fail rather
than silently changing the measurement.

**Observed result** (2026-10-20, Linux x86_64, Rust 1.85.0, `target/qa-coverage-v5`):
the invocation above returned exit 0 after 531.5 seconds of runner time. This
run predates the CLI removal (decision 029): its counts include the deleted `bioferris`
binary and the integration tests that ran it, and it has not been repeated since.
Complete runner stdout:

```text
QA coverage: PASS
```

Complete stdout of the report step:

```text
coverage report: 130 objects (14 nested examples), 5987 profiles; scoped JSON, LCOV and summary ready
```

All 1,090 listed tests ran and passed in 114 harnesses, with none failed,
ignored, measured or filtered, and source hashes were unchanged. The objects are
116 Cargo inventory executables plus 14 examples that documentation tests built
through `cargo bio`; all 130 profile binary IDs were matched. The report maps
266 source files, all inside the workspace. `TOTAL` row of `coverage-summary.txt`:

| Measure | Count | Missed | Covered |
|---|---:|---:|---:|
| Regions | 54,005 | 13,288 | 75.39% |
| Functions | 8,032 | 1,678 | 79.11% |
| Lines | 104,635 | 17,487 | 83.29% |
| Branches | 0 | 0 | not measured |

These are structural measurements of the declared six-package scope on this
host, not a fidelity or correctness score. Per-file rows are in the retained
`coverage-summary.txt`, and the
results archive
lists the run's receipts and reports. Earlier attempts `qa-coverage-v1` to `v4`
are described in the execution record;
none is counted as this result.

## Measurement pipeline

1. The durable QA runner creates a new output directory and starts a receipt.
   A separate `coverage-prepare` child checks the compiler/LLVM tools, then creates
   fresh `raw/` and `build/` children. Existing directories are never cleaned or
   reused.
2. Inventory compiles/runs the selected packages' all-target test listing with
   `--message-format=json`, `--offline --locked --jobs 2` and `-- --list --color never`,
   for the explicit host target. Build scripts and proc macros are therefore compiled
   without instrumentation; they run during compilation, not as measured tests.
3. The same selectors execute the actual tests with `-- --color never`.
   The runner compares listed/executed counts and refuses unexplained failed,
   ignored, measured or filtered tests. Expensive MSA compatibility unit suites
   are not among the selected packages; ordinary facade integration tests can
   still exercise transitive algorithms.
4. `coverage-report` reads executable paths from Cargo's actual
   `compiler-artifact` JSON records in `01-coverage-inventory.stdout.log`.
   Test harnesses, binaries and examples are retained and classified;
   binary filenames/hashes are never guessed. At least one selected test artifact
   is required, and every object must exist under this run's isolated build tree.
5. Raw profiles are merged with the pinned `llvm-profdata`. Its binary-ID listing
   is compared with each object's GNU build ID, read with the pinned `llvm-readobj`.
   An ID that no inventory object provides must belong to an example built in this
   run's isolated tree by a documentation test's nested `cargo bio` call; such
   examples are added as objects and listed in the manifest. Any other ID fails.
6. `llvm-cov` exports JSON, LCOV and a text report using every object, with
   `--check-binary-ids` so no profile is left without its executable.

Only the inventory/test steps set these explicit measurement variables:

```text
CARGO_TARGET_DIR=<absolute new output directory>/build
CARGO_INCREMENTAL=0
CARGO_BUILD_TARGET=<rustc host triple>
RUSTFLAGS=-C instrument-coverage
LLVM_PROFILE_FILE=<absolute new output directory>/raw/%p-%m.profraw
```

The absolute profile filename includes both process ID and binary signature, so
test subprocesses contribute independent profiles. Inventory execution can
also emit profiles; listing activity is not confused with having executed the
test bodies. A fresh isolated build prevents stale, uninstrumented artifacts
from being treated as this run's executable objects. Without the explicit target,
Cargo also passes `RUSTFLAGS` to build scripts and proc macros, whose compile-time
profiles have no test object; the first integrated report failed for that reason
(see the execution record).

## Artifacts and result interpretation

The new output directory retains runner logs/receipts plus:

| Artifact | Meaning |
|---|---|
| `coverage-environment.json` | Declared package scope, compiler/LLVM versions, instrumentation and explicit exclusions |
| `raw/*.profraw` | Original per-process/module profiles, retained rather than deleted |
| `raw-profiles.txt` | Sorted exact input list used for merging |
| `coverage.profdata` | Indexed merged profile |
| `coverage-binary-ids.txt` | Binary IDs recorded in the merged profile, matched against object build IDs |
| `coverage.json` | Complete LLVM JSON source coverage export |
| `coverage.lcov` | Complete LCOV export |
| `coverage-summary.txt` | LLVM text summary, including its denominator |
| `llvm-*.json` and `llvm-*.stderr.log` | Actual LLVM command arguments, statuses, output hashes and diagnostics |
| `coverage-manifest.json` | Cargo-derived executable classifications, nested example objects with build IDs, object/profile/report hashes, profile and binary-ID counts, mapped-file count and scope |

Vendor sources, registry sources and Rust sysroot-library mappings are explicitly
excluded from the source report. Standard-library paths are excluded in both forms:
remapped `/rustc/<commit>/library/` and an installed `rust-src` component's
`lib/rustlib/src/rust/library/`. After export, a Rust check mirroring the pattern
fails the report if any mapping still falls under an excluded location. The
complete exclusion regex is retained in both coverage manifests. All other
emitted source mappings remain included,
including tests and transitive workspace algorithm code. No BioFerris reader,
error path, parser or production core is silently excluded to improve a score.
The report fails if any of the six selected package directories has no source
mapping after export. This does not mean every source file was linked or every
reachable function executed.

There is deliberately no arbitrary percentage threshold in this initial
measurement. Review uncovered critical paths and per-file data first; a later
ratchet needs a reproducible baseline, explicit source/test denominators and an
approved policy. Region/line/function percentages are structural measurements,
not scientific correctness probabilities. Branch/MC/DC coverage is not promised
by the pinned stable instrumentation or these exports.

## Errors and limits

- Output must be a new child directory under workspace `target`; traversal,
  newline-containing/non-UTF-8 measurement paths and overwrite/reuse are refused.
- A missing pinned LLVM tool, wrong compiler, unknown/malformed artifact record,
  no selected test artifacts, missing/foreign objects, no profiles, empty profile
  or a non-profile raw-directory entry fails closed.
- A profile binary ID that matches neither an inventory object nor an example
  built in this run's isolated tree fails and names the ID.
- An exported mapping under an excluded vendor, registry or standard-library
  location fails and names the file.
- At most 20,000 raw profile files are accepted. The workflow does not claim a
  strict process-memory/disk budget; LLVM reports, isolated build artifacts and
  profiles can be large.
- Inventory and tests each have a 3,600-second runner timeout; preparation/report
  helpers each have a 1,200-second timeout. A helper failure or timeout is not
  turned into coverage success merely because a partial export exists.
- Every LLVM subprocess has separate captured output and a retained actual exit
  status. Missing/empty reports or missing selected-package mappings fail.
- The hidden preparation/report subcommands are implementation helpers, not a
  stable public biological API or permission to reuse another measurement's
  output directory.

Never regenerate scientific golden outputs or reopen a sealed campaign to improve
coverage. Investigate a failed test or missing source mapping instead. A missing
final runner summary means the run is incomplete, not PASS.

## Verification and boundaries

Observed with Rust 1.85.0:

```bash
cargo +1.85.0 clippy -p bioferris-qa --all-targets --offline --locked --jobs 2 -- -D warnings
cargo +1.85.0 test -p bioferris-qa coverage::tests:: --offline --locked --jobs 2
```

Both commands returned status 0. The focused coverage checks passed seven tests:
Cargo artifact parsing/classification, malformed/empty inventory refusal,
output-path rejection, the explicit host build target, parsing of the observed
binary-ID listing, fail-closed binary-ID matching and recognition of excluded
sources in both standard-library path forms. Other QA/receipt tests were
intentionally filtered by this focused invocation, not counted as covered
acceptance evidence. The integrated measurement under
[Inputs and invocation](#inputs-and-invocation) is the observed result; these
module checks do not substitute for it.

No numerical scientific assertion, algorithmic tie-breaking rule, biological
result, production compiler or package version is changed by coverage tooling.
Coverage is scoped ordinary test evidence. Cross-platform LLVM execution,
instrumented CPU-backend comparisons, full-workspace coverage and baseline
thresholds require separate recorded work. Existing frozen campaign evidence
remains historical and untouched.
