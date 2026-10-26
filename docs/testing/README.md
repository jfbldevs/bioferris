# Professional testing strategy

## Access index

| Purpose | Rust access | Development CLI |
|---|---|---|
| Ordinary regression | Native Rust tests across the workspace | `bioferris-qa run` |
| Durable evidence | Rust process runner and integration tests | `bioferris-qa exec` |
| Robustness and scientific contracts | Synthetic generators, independent oracles and library integration tests | `bioferris-qa run --profile robustness` |
| Diagnostic LLVM coverage | Rust orchestration, pinned compiler instrumentation | `bioferris-qa coverage` |
| Curated mutation audit | Rust source-copy/fault harness | `bioferris-qa mutation` |
| Synthetic scaling baseline | Release-mode Rust workload and bounded sampler | `bioferris-qa performance` |

The development CLI is the repository's QA tool, not a biological interface; the
biological API is the Rust library only ({D}).

## Contents

- [Acceptance and scope](#acceptance-and-scope)
- [Execution tiers](#execution-tiers)
- [Durable result contract](#durable-result-contract)
- [Manual and agent use](#manual-and-agent-use)
- [Requirements and specialist guides](#requirements-and-specialist-guides)
- [Tooling and platform policy](#tooling-and-platform-policy)
- [Implementation state](#implementation-state)

[API documentation index](../README.md) · Testing entry point

## Acceptance and scope

The user approved risk-based assurance across reproducible execution,
requirements traceability, scientific fidelity, parsers, publication/resources,
Rust/CLI/agent contracts (the CLI and agent contracts ended with decision 029),
coverage/mutation effectiveness and portability/performance.
All remain part of the objective; completing a small fast suite does not close
the plan. Each release requires current-source evidence, not copied historical
counts. Preserve the frozen once-only MSA campaign and its reserved inputs.

Requirements are mapped in requirements.json;
`requirements_traceability.rs`
requires every row to name existing test files. The map links tests and explicit
acceptance/scope, not a percentage of Biopython parity.
New bugs require minimized regression cases. Do not relax exact numerical
assertions, automatically refresh expected results, hide failures behind retries,
or silently skip difficult suites. Coverage is diagnostic, not proof of fidelity.

## Execution tiers

| Tier | Execution | Intended cadence |
|---|---|---|
| Fast | Formatting, strict toolkit/QA Clippy, named all-target test inventory and tests, separate toolkit doctests | Each change / pull request |
| Full | Strict all-feature workspace Clippy, ordinary workspace tests and inventory, all-feature target compilation/doctests/release build/strict docs | Nightly / integration |
| Release | Full gates with all features enabled for the workspace runtime tests | Before release |
| Robustness | New parser/contract/science/resource targets and publication fault unit tests | Focused development; also contained in ordinary regression |
| Coverage | Instrumented six-package ordinary regression and pinned LLVM reports | Scheduled diagnostic baseline |
| Mutation | Baseline plus five fresh isolated semantic faults with exact targeted tests | Scheduled assertion-strength audit |
| Optimized | Toolkit/QA/kernel inventories and all-target tests in release profile | Before release / optimization changes |
| Dependencies | Pinned auditor, all-feature policy checks and dated database provenance | Each change / scheduled / release |
| Performance | Three fixed synthetic sizes, warmup and five measured samples per size, descriptive phase/RSS aggregates | Scheduled baseline; no timing thresholds |

The canonical Cargo gate definitions are in
gates.rs. Separate doctests are intentional:
`--all-targets` does not replace them. Separate doctest steps also list their
test names first and enforce the same zero-omission and name/count checks as
ordinary tests. Optional `scientific-json` workers are
compiled by the all-feature target gate. Existing documentation-driven recipe
tests run real Cargo aliases and compare complete published results.
Dependency policy scanning is included in full/release and separately runnable;
its [setup and observed security results](security.md) are never inferred from a
configured `deny.toml` alone. Set `BIOFERRIS_QA_DENY_BIN` or install the pinned
auditor in PATH before invoking full/release. No installation occurs implicitly.

## Durable result contract

QA output must be a new directory in a trusted local filesystem. Existing files
and directories are refused, including empty directories; no cleanup or overwrite
is performed. Schema version 1 receipts retain plan, commit/dirty status, source
hashes, toolchain/platform, selected environment, exact argument vectors, real
child exit/signal, elapsed time and separate stdout/stderr checksums.

`initializing.json` and `started.json` are explicitly `INCOMPLETE`. A missing `summary.json` never means
success. Each command has append-only started/child/result files. The final
summary is `PASS` only after every planned step succeeds, source hashes remain
unchanged, and test steps match the listed test-name multiset with no failed,
ignored, measured or filtered tests. A successful log audit cannot replace the
actual process exit status. Unexpected source changes invalidate the run.
Source comparison is a start/end check, not a concurrent-filesystem snapshot;
freeze sources and maintained guides during qualification. Metadata commands
have separate 60-second deadlines and receipts. Source hashing checks its
120-second budget between reads; it cannot interrupt a blocked kernel file read.
Each captured log is limited to 64 MiB, checked while polling and after exit.

Timeouts and signals are recorded separately from assertion failures. Logs are
regular files rather than pipes, avoiding output-buffer deadlock. Unix timeout
cleanup uses a new process group and the OS `kill` utility; Windows uses
`taskkill /T`. Unix also rejects a successful leader that leaves an inherited
process group alive; normal-exit tree quiescence is not established by the
Windows standard-library implementation, and the receipt states that limitation.
This is trusted development tooling, not a sandbox for malicious
programs that detach from their process group. Abrupt termination of the runner
can leave an incomplete run and live children: inspect retained PIDs/processes
before deciding whether to wait or restart. Never restart solely because an
observation handle expires. Filesystem responsiveness and OS utility availability
are prerequisites; this is not a power-loss durability or hostile-directory claim.

## Manual and agent use

From the repository root, build the pinned development tool:

```bash
cargo build -p bioferris-qa --offline --locked --jobs 2
```

Choose a fresh output path for every invocation:

```bash
./target/debug/bioferris-qa run --profile fast --output target/qa-fast
```

Full and release tiers use the same interface with their corresponding profile.
The command returns 0 only for a passing plan; 1 means a failed plan, and 2 means
an invalid invocation or an infrastructure error. Agents must inspect receipts,
not parse progress messages or infer success from Cargo stdout. Progress goes to
stderr; the final stdout line names the profile and state.

### Complete observed runner example

Literal input: the tool's diagnostic `probe streams` case, which writes one line
to each child stream and exits 0. This is a runner fixture, not a biological result.
After the build above, the following invocation ran from the repository root:

<!-- qa-pass-command -->
```bash
./target/debug/bioferris-qa exec --output target/qa-example-pass-v1 --timeout-seconds 5 -- ./target/debug/bioferris-qa probe streams
```

Observed complete runner stdout, exit 0:

<!-- qa-pass-output -->
```text
QA custom: PASS
```

Observed complete runner stderr (elapsed time is a measurement, not a constant):

```text
QA command: ./target/debug/bioferris-qa probe streams
QA command: PASS (5 ms)
```

Saved child stdout in `00-command.stdout.log`:

<!-- qa-child-stdout -->
```text
probe: stdout
```

Saved child stderr in `00-command.stderr.log`:

<!-- qa-child-stderr -->
```text
probe: stderr
```

The summary records the child's actual exit code 0; log text alone does not decide
PASS. To exercise a relevant failure, use the literal `probe fail` case:

<!-- qa-fail-command -->
```bash
./target/debug/bioferris-qa exec --output target/qa-example-fail-v1 --timeout-seconds 5 -- ./target/debug/bioferris-qa probe fail
```

Observed complete runner stdout, exit 1:

<!-- qa-fail-output -->
```text
QA custom: FAILED
```

Observed complete runner stderr:

```text
QA command: ./target/debug/bioferris-qa probe fail
QA command: FAILED (3 ms)
```

The child stdout log is empty; complete child stderr is:

<!-- qa-child-failure -->
```text
probe: deliberate failure
```

These create-only paths already contain the observed runs in the development
checkout. A repeat must choose a new path; reuse is refused with exit 2 and does
not modify the earlier run. The integration tests compare the literal Markdown
results against fresh isolated invocations, without copying timing constants.

Coverage and mutation have separate
[coverage](coverage.md) and [mutation](mutation.md) recipes. The workload example
in [faults and resources](faults-and-resources.md) shows literal inputs, complete
JSON output and the distinction between host measurements and enforced limits.

## Requirements and specialist guides

- [Parser robustness and the interrupted-read regression](robustness.md).
- [Independent scientific oracles and the removed CLI conformance suite](contracts-and-science.md).
- [Create-only publication faults and resource measurements](faults-and-resources.md).
- [Scoped LLVM coverage, inclusions and exclusions](coverage.md).
- [Curated mutation patches and runtime-failure classification](mutation.md).
- [Actual portable/AVX2 execution and unqualified platform lanes](portability.md).
- [Dated dependency checks and XML security remediation](security.md).
- [Controlled synthetic scaling and sampling policy](performance.md).

## Tooling and platform policy

The product and new harnesses remain Rust edition 2024 on Rust 1.85.0. No Python,
JavaScript, FFI, new unsafe application code or implicit compiler upgrades.
Coverage uses already-installed matching `llvm-tools`; no coverage plugin is
required. Stable Rust generators provide bounded adversarial testing. Optional
nightly/libFuzzer tooling is not installed or represented as completed work.

The CI workflow configures Linux x86_64,
Linux ARM64, macOS ARM64 and Windows x86_64 lanes with pinned action revisions,
read-only repository permissions, deadlines and retained incomplete artifacts.
It has not been pushed or executed remotely by this increment. Runner availability,
repository Actions settings and branch protections require repository-owner setup.
Configured platform support is not measured numerical portability qualification.

Maintained source/guides use explicit LF checkout/formatting; byte-exact fixtures
and frozen evidence remain exempt from normalization. This prevents host-native
formatting from redefining documentation examples without rewriting old evidence.

## Implementation state

The execution record tracks the full
requirement-by-requirement completion audit and checkpoint distinctions.

Final integrated counts and the release, optimized, coverage, mutation,
performance and dependency outcomes are recorded only in that execution record
and its results archive. Both are outside the hashed source selection, so
recording an outcome does not change the qualified inputs. The audit found and
reproduced transient `Interrupted` read handling and identified two upstream XML
denial-of-service advisories; both were fixed in a separate commit with focused
regression tests. Unavailable platform evidence remains an explicit limitation,
not PASS.
