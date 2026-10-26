# Publication faults and bounded resource assurance

## Contents

- [Scope and access](#scope-and-access)
- [Publication failure matrix](#publication-failure-matrix)
- [Independent selected-span digest check](#independent-selected-span-digest-check)
- [Run the library workload](#run-the-library-workload)
- [Bounds and error example](#bounds-and-error-example)
- [Isolated adversarial and deadline checks](#isolated-adversarial-and-deadline-checks)
- [Measurement provenance and interpretation](#measurement-provenance-and-interpretation)
- [Verification and limitations](#verification-and-limitations)

[Documentation index](../README.md) · [Testing source](../../seq-bio/tests/resource_assurance.rs)

## Scope and access

| Capability | Rust access | Command access |
|---|---|---|
| Create-only publication fault injection | Internal `cfg(test)` seam in [publication.rs](../../seq-bio/src/publication.rs); no public/runtime hook | Targeted Rust unit tests |
| Plain/BGZF sequence-index build, open and lookup | Existing public `SequenceIndex` / `BgzfSequenceIndex` | `cargo bio qa_workload` development example |
| Streaming BLAST import and bounded adversarial input | Existing public readers and resource options | Rust integration-test subprocesses |
| Time and Linux peak RSS observations | `Instant` and `/proc` reads in the development workload | JSON QA report, not a biological result schema |

This increment adds test infrastructure, not biological algorithms, publication
semantics or a new supported `bioferris` command. Manual use remains a Rust
program: the complete [qa_workload example](../../seq-bio/examples/qa_workload.rs)
imports the real library. The integration harness imports that same example;
it does not simulate index construction, compression, validation or retrieval.
See the existing [plain indexing](../api/sequence-index-v1.md),
[compressed I/O](../api/compressed-io-v1.md) and
[BLAST input](../api/blast-input-v1.md) contracts for supported biological profiles.

## Publication failure matrix

The shared publisher writes a newly created temporary file in the trusted
destination directory, runs the callback, synchronizes the file, then publishes
it using a create-only hard link. Its temporary name is removed on scope exit.
There is no overwrite or copy fallback. Tests use real temporary files and real
publication code, with fault routing compiled only under `cfg(test)`.
Non-test builds retain the native filesystem calls and production API.

| Failure or boundary | Verified outcome |
|---|---|
| Destination inspection / temporary-file creation fails | Original `Io` error; callback not called; no published file |
| Every byte-prefix of `>dna\nATGCatgcNN\n` fails during writing | Callback `Write` error; no sync/link; temporary data removed when removal succeeds |
| Zero-length accepted write | `write_all` returns `WriteZero`; no publication |
| Interrupted and one-byte short writes | Complete unchanged payload is published |
| Callback explicitly flushes and flush fails | Original callback `Write` error; no publication |
| Callback returns its own typed error or panics | Error type / panic preserved; unwinding removes temporary data when possible |
| File synchronization or hard link fails | Original `Io` error; no named partial result |
| Hard link reports `AlreadyExists` | `DestinationExists`, not a generic I/O error |
| All 128 temporary candidates already exist | `TemporaryUnavailable`; all existing sentinel files unchanged |
| Existing file, directory or Unix dangling/live symlink | Rejected before callback; existing object unchanged |
| Destination created during callback | Competitor's bytes retained; our publication rejected |
| Four actual concurrent publishers | Exactly one complete winner; three `DestinationExists`; no temporary residue when cleanup succeeds |
| Temporary removal itself fails | Existing best-effort behavior: original result retained and an unpublished temporary alias may remain |

The last row is an explicit limitation, not a newly repaired guarantee. After a
successful hard link, failed removal may leave a second name for the completed
file; after an earlier error, it may leave an unpublished private temporary file.
The publisher cannot guarantee cleanup when the filesystem refuses it. It does
not promise hostile-directory safety, snapshots or crash/power-loss durability.
It synchronizes its `File`; it does not add a new implicit callback flush or
finalize caller-owned compression/buffering. Compression completion remains
explicit under the [compressed-output contract](../algorithms/compressed-output-v1.md).

Run from the repository root:

```sh
cargo test --locked -p bioferris --lib publication::tests --jobs 2
```

Observed result, summary excerpt from stdout; exit status `0`:

```text
test result: ok. 14 passed; 0 failed; 0 ignored; 0 measured; 173 filtered out; finished in 0.01s
```

Eleven tests were added to the three existing publication tests. The filtered
count and elapsed time describe this checkpoint, not a permanent test-count API.
Fault plans are thread-local and reset by guards; parallel tests do not share an
injected runtime filesystem policy. The test-only writer forwards accepted bytes
to the actual temporary file, rather than an in-memory substitute.

## Independent selected-span digest check

The [index unit test](../../seq-bio/src/io/index.rs) calls the real `read_entry`
path on this stable literal source:

```fasta
>a
ACGT
```

First, its correct SHA-256 must return record ID `a`, source sequence `ACGT`.
Then only the expected entry digest changes to 64 zeroes. Reading must return
`IndexError::StaleSource` with message `record span bytes differ from the verified
SHA-256`. File contents, length, identity and timestamps are checked unchanged
before and after. No production index guard or API is altered.

This isolates an important layered-validation boundary: a same-size source edit
with restored modification time can be rejected by Unix ctime checks before the
digest path runs. That integration scenario remains valuable, but cannot alone
demonstrate that the independent raw-span digest comparison is effective.
The [mutation audit](mutation.md) targets this isolated assertion.

```sh
cargo test --locked -p bioferris --lib io::index::tests::selected_span_digest_is_checked_independently_of_source_metadata --jobs 2 -- --exact
```

Observed test stdout below; exit status `0` (Cargo build messages are on stderr):

```text
running 1 test
test io::index::tests::selected_span_digest_is_checked_independently_of_source_metadata ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 187 filtered out; finished in 0.00s
```

## Run the library workload

Prerequisites: the pinned Rust 1.85.0 toolchain, the repository root as the working
directory and a trusted writable temporary directory. This Bash invocation
creates a unique empty directory; it never overwrites existing artifacts:

<!-- qa-workload-command -->
```sh
qa_directory="$(mktemp -d /tmp/bioferris-qa-example.XXXXXX)"
cargo bio qa_workload --jobs 2 -- --directory "$qa_directory" --records 2 --symbols 4
```

The literal DNA source generated as `sequences.fa` is:

<!-- qa-workload-fasta -->
```fasta
>record_000000 synthetic
aCgT
>record_000001 synthetic
aCgT
```

The literal two-column BLAST tabular source generated as `blast.tsv` is below.
The separators are actual tab bytes; columns are explicitly `qseqid,sseqid`:

<!-- qa-workload-blast -->
```text
q000000	s000000
q000001	s000001
```

Each plain and BGZF index is built, saved create-only, dropped and reopened with
full verification. Every record is retrieved and compared for source case and
complete header; a missing key must return `None`. BLAST is read query-by-query;
every query ID, hit ID, HSP count and final query count is checked. This imports
reported results; it does not execute BLAST. Separate phases distinguish build,
save, verified open and all-record lookup rather than conflating them.

Complete observed stdout below; stderr was empty and exit status was `0`.
The literal timings/RSS are one real debug-build observation, not golden limits:

<!-- qa-workload-report -->
```json
{
  "bgzf_bytes": 103,
  "kind": "bioferris-qa-workload",
  "phases": [
    {
      "elapsed_us": 39,
      "phase": "source_write",
      "process_peak_rss_kib": 3972
    },
    {
      "elapsed_us": 96,
      "phase": "plain_build",
      "process_peak_rss_kib": 4404
    },
    {
      "elapsed_us": 6255,
      "phase": "plain_save",
      "process_peak_rss_kib": 4468
    },
    {
      "elapsed_us": 130,
      "phase": "plain_open",
      "process_peak_rss_kib": 4476
    },
    {
      "elapsed_us": 44,
      "phase": "plain_lookup",
      "process_peak_rss_kib": 4480
    },
    {
      "elapsed_us": 1434,
      "phase": "bgzf_write",
      "process_peak_rss_kib": 4960
    },
    {
      "elapsed_us": 389,
      "phase": "bgzf_build",
      "process_peak_rss_kib": 4960
    },
    {
      "elapsed_us": 1276,
      "phase": "bgzf_save",
      "process_peak_rss_kib": 4960
    },
    {
      "elapsed_us": 375,
      "phase": "bgzf_open",
      "process_peak_rss_kib": 5004
    },
    {
      "elapsed_us": 106,
      "phase": "bgzf_lookup",
      "process_peak_rss_kib": 5004
    },
    {
      "elapsed_us": 19,
      "phase": "blast_write",
      "process_peak_rss_kib": 5004
    },
    {
      "elapsed_us": 43,
      "phase": "blast_read",
      "process_peak_rss_kib": 5004
    }
  ],
  "process_peak_rss_kib": 5004,
  "records": 2,
  "schema_version": 1,
  "source_bytes": 60,
  "symbols_per_record": 4,
  "verified_bgzf_records": 2,
  "verified_blast_queries": 2,
  "verified_plain_records": 2
}
```

The example retains its five create-only artifacts (`sequences.fa`,
`sequences.index.json`, `sequences.fa.bgz`, `sequences.bgzf.index.json`,
`blast.tsv`) for inspection. It performs no implicit cleanup of caller files and
is not a multi-file transaction: a later error may leave earlier outputs. Use a
new directory for another run. A failed source write can leave its own incomplete
diagnostic source; named BGZF/index publication uses the shared atomic create-only
publisher. The integration tests clean only their own newly created directories.

## Bounds and error example

The workload accepts `--directory PATH`, `--records N` (default 2048) and
`--symbols N` (default 256). The directory must already exist. Duplicate/unknown
options, missing values and invalid integers are rejected. Inclusive bounds:

| Resource | Bound |
|---|---|
| Records | `1..=100000` |
| Symbols per record | `1..=4096` |
| Records × symbols | At most `16777216` |
| Reader physical line / record symbols | `8192` bytes / `4096` symbols |
| Decoded FASTA source | `32 MiB` |
| Aggregate index ID bytes / sidecar bytes | `2 MiB` / `64 MiB` |
| Compressed BGZF bytes / members | `32 MiB` / `4096` |
| BLAST hits / HSPs per query | `1` / `1` |
| BLAST bytes per query | `128` |

These input/count budgets are not allocator or operating-system memory caps.
BLAST's other limits retain their documented defaults. Input generation remains
bounded before file creation; it does not accept arbitrary unbounded user data.

Literal invalid input: zero records, four symbols, the same existing directory:

<!-- qa-workload-error-command -->
```sh
cargo bio qa_workload --jobs 2 -- --directory "$qa_directory" --records 0 --symbols 4
```

Observed stdout: empty. Complete stderr; exit status `2`:

<!-- qa-workload-error -->
```text
qa_workload: workload bounds: 1..=100000 records, 1..=4096 symbols, product <=16777216
```

Invalid size arguments are checked before files are created. Existing output
paths are refused without modification; native filesystem messages are
platform-dependent. This development example uses a human stderr error, not a
typed library error. Agents should call the Rust library for biology and check
Cargo exit status when running QA workflows.

## Isolated adversarial and deadline checks

Run all new resource checks:

```sh
cargo test --locked -p bioferris --test resource_assurance --jobs 2 -- --nocapture
```

The complete documented small-input/count report is checked against the example;
timings/RSS are type-checked observations, not reproducible numerical expectations.
The suite also runs these real library operations in its own test executable as
a subprocess, with file-backed logs rather than potentially blocking pipes:

- Plain FASTA and BGZF index build/save/open/all-record lookup, plus streaming
  BLAST at `32 × 64` and `2048 × 256` record/symbol sizes; 30-second deadlines.
- Two-record FASTA `>first\naCgT\n>second\nNN\n` with `max_records=1`;
  a BGZF encoding of exactly that source with `max_compressed_bytes=1`.
- 4096 tabular rows `q\ts0\n` through `q\ts4095\n`, all for one query, with
  at most four hits/HSPs; the reader fails and becomes exhausted without emitting
  an incomplete query.
- On-demand `>id\n` followed by 1048576 `A` bytes without a newline;
  a 32-byte line budget, 64-symbol record budget, 128-byte total budget and
  64-byte transport buffer. The parser must reject without draining the source.
- A child acknowledges readiness, then waits; its parent imposes a 50 ms deadline,
  terminates only that owned child and waits for it. Worker startup has a separate
  10-second bound. The test verifies that the deadline path is actually exercised.

The adversarial worker has a 10-second deadline. Errors are typed and checked as
`limit_exceeded`. Observed diagnostic excerpts from stderr in one complete run:

```text
adversarial_limits: {"bgzf_index_error":"limit_exceeded","blast_error":"limit_exceeded","blast_successful_queries":0,"fasta_error":"limit_exceeded","fasta_transport_bytes":68,"plain_index_error":"limit_exceeded","process_peak_rss_kib":5628,"sidecars_created":0}
resource_workload: records=32 symbols=64 elapsed_ms=19 sampled_peak_rss_kib=Some(5800) worker_peak_rss_kib=5800
resource_workload: records=2048 symbols=256 elapsed_ms=932 sampled_peak_rss_kib=Some(7320) worker_peak_rss_kib=7052
```

The 68 observed transport bytes include bounded `BufReader` read-ahead; they are
not a universal physical-I/O bound for arbitrary upstream buffering. The suite
asserts at most 128 bytes for this precise transport. It verifies no index
sidecars were created by rejected builds.

## Measurement provenance and interpretation

Observed on 2026-10-18, working tree based on `24b67664` with this test-only
increment and concurrent parser-assurance fixes, using:

- `rustc 1.85.0 (4d91de4e4 2025-02-17)` and
  `cargo 1.85.0 (d73d2caf9 2024-12-31)`; debug builds, no new dependencies.
- Linux `6.18.33.2-microsoft-standard-WSL2`, `x86_64`;
  Intel Core i9-14900HX; kernel-reported memory `16230096 KiB`, swap `4194304 KiB`.
- `/proc/self/status` `VmHWM` for the worker's cumulative process high-water RSS;
  parent samples of `/proc/<owned-child>/status` `VmRSS` and `VmHWM` every roughly
  2 ms while waiting. No privileged operations or external measurement runtime.
- `Instant` wall time for each library phase, excluding the measurement read
  after that phase. Parent elapsed time includes subprocess startup and polling,
  but not Cargo compilation. The complete example JSON is a separate run.

Per-phase RSS is **cumulative process high water**, not isolated phase allocation
or a delta. Parent sampling may observe later peaks during report serialization
and test teardown that the worker's last self-report misses. Kernel granularity,
page cache, allocator state, scheduling, debug code and other host workloads all
affect the values. Missing/inaccessible procfs or non-Linux platforms report
`null`/`None`; logical boundedness/deadline tests do not need RSS measurements.
No hard portable RSS cap, benchmark ranking, scaling law or production throughput
claim is established by these small synthetic observations.

## Verification and limitations

Targeted publication tests, resource tests, the small Cargo example, its zero-size
error and strict targeted Clippy were executed. Root integration owns the final
all-target/workspace gates. The tests use no runtime installation, unsafe code,
FFI or external scientific subprocess, and never read or rerun sealed MSA data.

This is systematic coverage of selected failure boundaries, not an exhaustive
filesystem fault model. It does not emulate every OS error, arbitrary callback
misbehavior, disk corruption, power loss, all supported parser formats or every
possible resource-exhaustion attack. Deadlines contain test subprocess lifetime;
they are not public library cancellation or OS resource limits. Existing
scientific compatibility qualification remains unchanged.
