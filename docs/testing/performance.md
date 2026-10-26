# Initial synthetic scaling baseline

## Contents

- [Access and status](#access-and-status)
- [Literal workload definition](#literal-workload-definition)
- [Invocation and results](#invocation-and-results)
- [Observed baseline](#observed-baseline)
- [Sampling and interpretation](#sampling-and-interpretation)
- [Artifacts, errors and limits](#artifacts-errors-and-limits)
- [Verification and scope](#verification-and-scope)

## Access and status

| Access | Availability |
|---|---|
| Repository QA command | `bioferris-qa performance --output NEW_DIRECTORY`, controlled baseline verified for the recorded host and inputs |
| Rust workload driver | Existing `qa_workload` example, invoking real public library APIs |
| Biological API | No new biological operation, performance promise or result schema |

[API index](../README.md) · [QA runner](README.md) ·
[Driver](../../seq-bio/examples/qa_workload.rs) ·
Sampling controller

This is an initial, repeatable, descriptive baseline for a declared synthetic
workload, not a reference-implementation comparison. Correct content/count
verification is mandatory before retaining a sample as valid. Timing thresholds,
speedups, dedicated-hardware claims, affinity claims and statistical significance
are deliberately absent. Historical algorithm benchmarks remain separate.

## Literal workload definition

The complete generated inputs are fixed by these three pairs, in this order:

```text
records=128  symbols_per_record=64
records=2048 symbols_per_record=256
records=4096 symbols_per_record=4096
```

For each pair `(N, S)`, write exactly `N` FASTA records, in ascending order
`i = 0..N`. Each record has these complete bytes:

```text
>record_{i formatted as six decimal digits} synthetic\n
{first S bytes of the infinite repetition aCgTaCgTaCgT...}\n
```

The `\n` above denotes one LF byte, not a literal backslash or an extra empty
line. No wrapping, randomness, clock values or unknown biological data are used.
For example, the complete first record of the smallest case is:

```fasta
>record_000000 synthetic
aCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgTaCgT
```

Every record contains exactly `S + 26` source bytes, including header/newlines.
The case-preserving source is DNA. The driver builds, saves, fully verifies on
open and retrieves all records from both plain and BGZF indexes; every returned
header/sequence must match. An absent-ID lookup must remain absent. BGZF output
uses the existing default writer options and explicit completion.

The BLAST input contains exactly these outfmt-6 custom two-field rows, one per
ascending `i = 0..N`:

```text
q{i formatted as six decimal digits}\ts{i formatted as six decimal digits}\n
```

Here `\t` is one tab byte and `\n` is one LF. Fields are explicitly `qseqid`
and `sseqid`; every parsed query must have the expected identifier, one hit and
one HSP. No score, coordinate, alignment or search execution is fabricated.

## Invocation and results

The actual controlled invocation, from the repository root with the already built
QA executable, was:

```bash
./target/debug/bioferris-qa performance --output target/qa-performance-v1
```

The reproducible Cargo alternative below additionally builds/launches that QA
executable using pinned Rust 1.85.0 and cached locked dependencies. It is an
alternate invocation, not the command recorded for this measurement. Select a
different new directory if `target/qa-performance-v1` already exists:

```bash
cargo +1.85.0 run -p bioferris-qa --offline --locked --jobs 2 -- performance --output target/qa-performance-v1
```

No benchmark, compile or test job should be running concurrently during the
measurement. The controller runs its workload children sequentially; it cannot
establish exclusive use of the machine or the absence of unrelated OS activity.

The controlled run completed with **PASS**, process status **0**, and
`source_unchanged: true` in the durable runner summary. All 18 child receipts
also report PASS, exit 0 and no signal: three warmups and 15 measured samples.
No competing build/test job was run during sampling; this does not establish a
dedicated machine or eliminate unrelated OS activity.

Complete runner stdout:

```text
QA performance: PASS
```

The controller's complete stdout, retained separately in
`target/qa-performance-v1/01-performance-run.stdout.log`, is:

```text
performance baseline: 3 sizes, 5 measured samples each; descriptive summary ready
```

Controller stderr was empty. The release-build step took 11,443 ms and the
controller step 5,872 ms; both returned 0. Those orchestration durations are not
individual library-phase timings. The final retained results tree was
122,286,157 bytes, below the 512-MiB bound. The longest complete child receipt was
822,868 µs, below its 60-second deadline.

### Observed baseline

The following values are copied from the complete measured summary; every
phase/case has five valid measured samples, with the warmup excluded.

| Phase (µs, min / median / max; n = 5) | 128 × 64 | 2,048 × 256 | 4,096 × 4,096 |
|---|---:|---:|---:|
| `source_write` | 27 / 29 / 32 | 232 / 236 / 403 | 4992 / 5076 / 5635 |
| `plain_build` | 208 / 213 / 226 | 5605 / 5726 / 5990 | 102893 / 104260 / 106625 |
| `plain_save` | 1824 / 1854 / 2033 | 20674 / 21141 / 23339 | 37800 / 39506 / 41662 |
| `plain_open` | 281 / 283 / 289 | 7903 / 8105 / 9255 | 139176 / 140445 / 148125 |
| `plain_lookup` | 353 / 357 / 363 | 7513 / 7749 / 8295 | 75209 / 76688 / 78527 |
| `bgzf_write` | 832 / 841 / 925 | 3313 / 3658 / 4005 | 22801 / 23631 / 25108 |
| `bgzf_build` | 273 / 281 / 323 | 4801 / 5020 / 6469 | 83396 / 84763 / 88491 |
| `bgzf_save` | 2002 / 2217 / 2555 | 24624 / 25272 / 30318 | 47576 / 50398 / 50627 |
| `bgzf_open` | 320 / 326 / 333 | 6081 / 6530 / 6770 | 86804 / 88956 / 90062 |
| `bgzf_lookup` | 1130 / 1153 / 1262 | 27047 / 28251 / 28852 | 170772 / 174276 / 174999 |
| `blast_write` | 23 / 23 / 25 | 123 / 148 / 152 | 223 / 232 / 272 |
| `blast_read` | 171 / 215 / 251 | 2369 / 2524 / 2661 | 4752 / 4792 / 4908 |

| Process high-water RSS | 128 × 64 | 2,048 × 256 | 4,096 × 4,096 |
|---|---:|---:|---:|
| Maximum observed KiB across five measured processes | 3068 | 4332 | 6420 |

RSS was available for all five samples at every phase. The raw JSON retains each
phase's maximum high-water RSS as well as the overall process maximum above.
The host metadata records Linux/x86_64, Intel Core i9-14900HX, 32 available
logical CPUs and Rust 1.85.0 / LLVM 19.1.7. CPU0 governor metadata was unavailable
(null); affinity and dedicated-machine status were not established.

Both record count and record length change between these cases. Three points
and five descriptive samples do not establish an asymptotic complexity model,
statistical significance, superiority to another implementation or a resource
guarantee for different data.

The first runner step builds `qa_workload` with `--release`, `--message-format=json`,
`--offline --locked --jobs 2`. The next step invokes the controller against the
new output directory. It derives the exact executable from the first step's
`compiler-artifact` JSON, checking the example source, name, non-test profile and
optimization level. It never guesses a target filename.

## Sampling and interpretation

Each size runs one warmup child followed by five measured children, all with
fresh create-only data directories. Warmup data/results remain inspectable but
are excluded from aggregates. The controller never overwrites or silently
retries a failed sample. There are 18 total child executions and 15 measured
samples. The largest case contains exactly 16 MiB of residue bytes per child.

For each size, output records `n = 5` and min/median/max `elapsed_us` for each
ordered driver phase:

```text
source_write
plain_build
plain_save
plain_open
plain_lookup
bgzf_write
bgzf_build
bgzf_save
bgzf_open
bgzf_lookup
blast_write
blast_read
```

Phase durations come from the driver's own monotonic timer. Parent process
elapsed time, process startup, orchestration and JSON parsing are retained
separately and are not substituted for phase timings. No time subtraction or
normalization silently changes the measured work.

Linux `process_peak_rss_kib` is the process-wide cumulative high-water RSS from
`/proc/self/status` when available, not per-phase allocations or a phase-specific
increment. Aggregation reports the maximum available observation and the number
of available samples. Unavailable values remain null, not zero. Other operating
systems or unavailable procfs can legitimately have no RSS measurement.

The environment artifact records the pinned compiler, OS, architecture,
available parallelism and available Linux CPU-model/CPU0-governor/online-CPU
metadata. Missing metadata remains null. Reading CPU0's governor is not proof
that the child ran there, that all cores share it, that turbo is disabled or
that frequency/temperature remained fixed. There is no affinity or dedicated
machine guarantee, confidence interval, p-value, or native-reference speedup.

## Artifacts, errors and limits

The new output directory contains the durable runner's logs/receipts plus:

| Artifact | Contents |
|---|---|
| `00-performance-build.stdout.log` | Original Cargo JSON used to locate the actual release example |
| `performance-environment.json` | Executable/build-log SHA-256, compiler/host metadata, sizes and complete sampling/resource policy |
| `results/NxS/warmup-00/` | Retained warmup input/output and actual status; excluded from statistics |
| `results/NxS/sample-00/` through `sample-04/` | Each complete raw stdout JSON, stderr log, started record, child receipt and private data files |
| `performance-summary.json` | Success/failure state and descriptive per-size/per-phase aggregates |

The actual raw evidence is retained locally under `target/qa-performance-v1/`.
This ignored build-output directory is not guaranteed to exist in a fresh
checkout. Preserve or export its receipts before cleaning build artifacts.
The exact summary/environment/runner hashes observed for this run are:

```text
69920f8f4bfd24d1541b64123876ad29ce9d9490ae36729a880507bc58a6fef4  target/qa-performance-v1/performance-summary.json
2907dce5a1dd3a907842c7c9777ed671c3e1748b627be90a609ac84c0dc42075  target/qa-performance-v1/performance-environment.json
dcd13a9840ea0c87123ea039a229781db4f66e9f59cbb3b14b2c168d7b7cfbf9  target/qa-performance-v1/summary.json
```

The environment record additionally pins the actual release example SHA-256 to
`ab6b452eb38b7496b05f05f50004bbf4dfc71bea4d5f783b9a7d5b67ce22d85a`
and the Cargo build-log SHA-256 to
`cc2830e0b774bc660ec860a2489dd3215894ad49e34c2399c01eb117747f4c99`.
Each of the 18 `results/NxS/{warmup,sample}-NN/receipt.json` files retains the
actual child status, elapsed time and stdout/stderr hashes; complete individual
results and data are not replaced by the compact table above.

Each sample's `data/` holds the FASTA source, plain sidecar, BGZF source, BGZF
sidecar and BLAST TSV. A successful child must exit 0, emit no stderr, satisfy
the strict report schema, have the exact ordered 12 phases and match all declared
dimensions, source-byte counts and verified record/query counts. Unknown,
duplicate, missing, mistyped or inconsistent fields are not silently accepted.

Each child has a 60-second deadline, a one-MiB cap per captured output stream
and a private directory. The aggregate results tree is checked against a
512-MiB dataset/log budget and a 512-entry limit while children run and after
they exit. Log/disk polling is cooperative and may observe a short write burst
beyond the bound before termination; it is not an OS filesystem quota. The
reviewed finite driver additionally caps `records * symbols` at 16 MiB.
Release-build artifacts are outside this retained dataset/log budget.

Timeout, nonzero status, unexpected stderr, malformed report or resource failure
is retained as failure, not included in the summary as a successful sample.
Partial files remain available; no cleanup, overwrite, threshold-based retry or
selection of the fastest five is performed. Missing final summary means
incomplete evidence, not success.

Ambient `RUSTFLAGS`, `CARGO_ENCODED_RUSTFLAGS`, `CARGO_BUILD_TARGET` and
`CARGO_TARGET_DIR` are refused before constructing the build plan. This keeps
coverage instrumentation, CPU-specific overrides and foreign target trees out
of the default release baseline. The workflow installs nothing and requires no
network service, external reference runtime, unsafe/FFI, sealed data or campaign.

## Verification and scope

Observed on Rust 1.85.0:

```bash
cargo +1.85.0 clippy -p bioferris-qa --all-targets --offline --locked --jobs 2 -- -D warnings
cargo +1.85.0 test -p bioferris-qa performance::tests:: --offline --locked --jobs 2
```

Both commands returned status 0. The focused run passed three performance-module
tests: Cargo example selection/ambiguity refusal, five-sample aggregation with
unavailable RSS retained as null, and malformed/inconsistent phase-report
rejection. Other QA tests were deliberately filtered in that focused invocation;
they are not claimed as execution evidence here. The separate controlled
experiment above now records the measured baseline; these unit/lint results
are not substituted for it.

These synthetic sizes establish an initial local scaling observation only.
They do not model all real sequence distributions, adversarial inputs, cold-cache
storage, parallel throughput, ordinary-gzip indexing, concurrent writers or
publication-grade benchmarking. Five descriptive samples are not a statistical
significance claim. API correctness, scientific-reference comparisons, coverage,
mutation testing and robustness remain separate evidence layers.
