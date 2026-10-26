# Bounded reader and compressed-stream robustness

## Contents

- [Access and scope](#access-and-scope)
- [Run the suite](#run-the-suite)
- [Worked input and actual result](#worked-input-and-actual-result)
- [Test inventory](#test-inventory)
- [Corpus provenance and limits](#corpus-provenance-and-limits)
- [Interrupted-read regression](#interrupted-read-regression)
- [Verification and limitations](#verification-and-limitations)

## Access and scope

| Access | Availability |
|---|---|
| Rust public readers and writers | Exercised directly by this suite |
| Manual or CI execution | Cargo integration-test target `robust_streaming` |

[API index](../README.md) · [File profiles](../api/file-io-v1.md) ·
[BLAST profiles](../api/blast-input-v1.md) ·
[Compressed streams](../api/compressed-io-v1.md)

These finite, deterministic tests extend the existing format, reference-fixture
and property suites. They check transport invariance, preservation, termination,
typed errors and integrity boundaries. They are not a new biological algorithm,
coverage percentage, performance benchmark, universal compatibility claim or
scientific qualification campaign.

## Run the suite

From the repository root, with the existing pinned Rust 1.85.0 toolchain and cached
lockfile dependencies:

```bash
cargo +1.85.0 test -p bioferris --test robust_streaming --offline --locked --jobs 2
```

Success has Cargo status 0; a failed assertion has status 101. Cargo's build
diagnostics use stderr, and its test harness prints the result inventory to
stdout. No network access, reference interpreter, new dependency, reserved MSA
input or external alignment engine is required.

Source: [tests](../../seq-bio/tests/robust_streaming.rs) and
[bounded transport/generator support](../../seq-bio/tests/robust_support/mod.rs).

## Worked input and actual result

The source is Sanger FASTQ, DNA, with four Phred+33 qualities. Input bytes use LF
and include the final newline:

<!-- robustness-input -->
```fastq
@sample
AcGT
+
!I~+
```

The complete runnable invocation below executes the checked recipe in the source
test, using the actual public `FastqReader` and seven short-read schedules:

<!-- robustness-command -->
```bash
cargo +1.85.0 test -p bioferris --test robust_streaming --offline --locked --jobs 2 documented_robustness_recipe -- --exact --nocapture
```

Complete recipe payload on stdout (Cargo compilation/harness progress is not
part of this payload):

<!-- robustness-output -->
```text
records: 1
id: sample
source: AcGT
normalized: ACGT
phred: [0, 40, 93, 10]
chunk schedules: 7
fused after invalid quality: true
invalid quality code: invalid_quality
```

The recipe also checks this complete malformed input, whose quality contains a
space. `FastqReader` returns `IoErrorKind::InvalidQuality`, yields no record for
that read, and stays exhausted without touching its source again:

<!-- robustness-invalid -->
```fastq
@sample
AcGT
+
!I ~
```

The regression test reads both literal inputs, this invocation and the expected
payload from this page. Successful quality decoding preserves source case and
separately normalizes the computational sequence. The expected error is handled
by the test, so this recipe command still exits 0.

## Test inventory

| Test area | Independent expectation and exercised variation |
|---|---|
| Generated FASTA/FASTQ | 16 three-record cases; mixed case/IUPAC, UTF-8 headers, CRLF, wrapping and generated qualities; compare complete `Record` values built independently |
| Sequence mutations | Seeded deletion/insertion/bit/control edits and arbitrary bytes; complete record/typed-error context invariant across chunks; accepted records round-trip without case/quality loss |
| Reader limits | Every limit from zero through fixture length + 1 for line, record and total bytes, across all chunk schedules |
| Sequence transport | One transient interruption or permanent failure at every source offset, including EOF; partial-line retry and error fusion |
| Alignment import | Six generated rectangles in eight profiles: FASTA, CLUSTAL, Stockholm, NEXUS and all four explicit PHYLIP combinations; independent row/width checks plus mutations |
| Graph/result import | Eight generated linear GFA/DOT graphs and single-deletion native distance documents; independent paths/topology/JSON, byte mutations and transport invariance |
| Eager document transport | Alignment/GFA/DOT/result JSON retry interruptions and reject permanent errors at every byte offset |
| Annotated input | 23 generated lengths for GenBank and EMBL; exact case, qualifiers and complement/join extraction against an independent mapping; declared-length mismatch and grammar edits |
| BLAST tabular 6/7 | Generated HSP groups, descending source coordinates, raw numerical lexemes and no-hit metadata; arbitrary/edited bytes and fused errors |
| BLAST legacy XML/XML2 | Independently generated schema-profile text, UTF-8/entity boundaries, hit/HSP coordinates, whole-model invariance, no-hit legacy query and exact decoded schemaLocation tab/newline preservation |
| BLAST transport | Every tabular offset and 17 evenly spread XML offsets including EOF; transient retry and permanent typed failure |
| XML security bounds | Both XML profiles: ordinary queries survive 256 total namespace declarations; 257 declarations, duplicates and 512 unsupported ordinary attributes are refused in a child with a 10-second watchdog |
| Independent gzip/BGZF | Hand-framed RFC stored-DEFLATE members and independent bitwise CRC32; concatenation, embedded empty members, all compressed-byte interruption offsets |
| BGZF integrity boundary | Every nonempty proper truncation of a member and high-bit edit at every member byte, behind a valid prefix; rejected member never releases payload or changes destination buffer |
| Compression writers | Seeded arbitrary bytes of lengths 0, 1, 7, 31, 255, 1,024 and 4,096; seven write schedules and bounded BGZF blocks; decoded bytes must match exactly |
| Gzip trailer failures | Every proper truncation from the recognized two-byte signature, and every CRC/ISIZE trailer byte corrupted; complete consumption must fail |
| Documentation | Literal FASTQ, invalid-quality input, invocation and full recipe payload checked by `documented_robustness_recipe` |

There are 18 Rust test functions; table rows describe areas, not additional
independent qualifications. Several functions exercise multiple corpus families.

## Corpus provenance and limits

All new corpus content is synthetic and generated by this test source. No public
reference fixture, reference implementation source, sealed manifest, held-out
campaign or network response is consumed. Format interpretation follows the
existing linked BioFerris contracts; the suite does not expand their profiles.

Xorshift32 uses explicit nonzero seeds, with no clock or ambient entropy. The
fixed seeds are `0x6d2b79f5`, `0x925a10f3`, `0x81e423a7`, `0x549bc621`,
`0x65b8a2d9`, `0x415f23cd`, `0x931c62ab`, `0x3f92a8e1`, `0xe1539c42`,
`0xa319827d`, `0xf713c4a9`, `0xc6217b3a` and `0x78ba301f`. Each mutation batch
contains the original input, four truncations and 32 single-site edits. The
arbitrary-byte batch contains 48 cases, each at most 96 bytes.

Transport source buffers are capped at 8,192 bytes; decoded gzip/BGZF output is
capped at 16 KiB. Read schedules cap each source read at 1, 2, 3, 5, 8, 13 or 64
bytes, behind a seven-byte `BufReader` where applicable. The one-byte schedule
crosses every physical byte seam, including UTF-8, CRLF, XML entities and framing.
BLAST tests additionally cap total input at 8 KiB, line/query/event/metadata
buffers at 4 KiB, depth at 32, queries/hits at eight and HSPs at 16 per query.
The namespace/attribute security probe raises its bounded physical-line, event
and metadata allowances to 8 KiB so the test reaches the parser's declaration
limit rather than being stopped by a smaller test-only transport allowance.
These are fixture/test bounds, not promises of exact process-memory limits.

Iteration and decoded-output checks have explicit finite budgets. Successful
stream exhaustion and terminal errors are followed by two further iterator
calls; both must return `None` without another underlying source read. A transient
interruption fires once and then permits progress. Hard failures remain errors;
they are not replaced with EOF or retried indefinitely.

## Interrupted-read regression

The initial new-target run returned status 101: 13 tests passed and three failed.
Two failures reproduced a production defect at source offset zero:

```text
interrupted_sequence_reads_retry_without_losing_a_partial_line:
  observed terminal Io/fasta at record 1, line 1 instead of two valid records
blast_transient_and_permanent_transport_failures_have_explicit_outcomes:
  observed terminal Io/BLAST tabular at record 1, line 1 instead of a valid query
```

The third failure was a test-preparation error: the synthetic legacy XML lacked
required root query fields. Those fields were added to make the input conform to
the existing contract; no product validation was weakened.

The shared line reader now retries `ErrorKind::Interrupted`, retaining any bytes
already accumulated for the current line. All other I/O errors remain terminal.
That retry was the only intended production behavior change at the first
robustness checkpoint: no reader API, grammar, biological semantics or resource
limit changes were needed for those transport regressions.
The original success expectations remain ordinary tests, not ignored tests or
expected failures. Every-offset tests also cover interruption inside a partial
line, not only the initial read. After the fix, all 16 then-existing tests passed
with status 0; the documentation recipe is the seventeenth check.

### Subsequent XML dependency security regression

The integration security audit then identified RustSec advisories
`RUSTSEC-2026-0194` (attribute handling) and `RUSTSEC-2026-0195` (namespace
declarations) for the previously pinned `quick-xml` 0.38.3. The integration owner
selected exact `quick-xml` 0.41.0, whose declared MSRV is 1.79, without changing
Rust 1.85.0 or disabling duplicate/namespace validation. Historical records of
the earlier dependency remain unchanged.

The eighteenth test checks the patched parser's 256-total-declaration boundary,
including XML2's existing default declaration, while preserving ordinary query
semantics. A duplicate declaration remains malformed; hundreds of unsupported
ordinary attributes remain explicitly unsupported, not silently dropped. A
self-subprocess watchdog limits this finite probe to 10 seconds. This bounded
regression does not establish an asymptotic complexity proof or replace a fresh
advisory scan. The targeted probe and all 18 tests pass after the minimal public
dependency API adaptation. A schemaLocation case separately pins literal tab/newline
and escaped tab/newline handling so the newer attribute API cannot silently add
whitespace normalization to this source-preserving profile.

## Verification and limitations

Observed on Rust 1.85.0 during this increment:

- Full `robust_streaming` invocation above, after the XML security update: status 0, **18 passed**, zero failed,
  ignored, measured or filtered tests.
- Exact documented recipe invocation above: status 0, one passed and 17 filtered;
  the complete payload matches the literal block on this page.
- `cargo +1.85.0 clippy -p bioferris --test robust_streaming --offline --locked --jobs 2 -- -D warnings`:
  status 0.
- Final `cargo +1.85.0 clippy -p bioferris --all-targets --offline --locked --jobs 2 -- -D warnings`:
  status 0, including the namespace watchdog and attribute-preservation checks.
- `rustfmt +1.85.0 --edition 2024 --check seq-bio/tests/robust_streaming.rs seq-bio/tests/robust_support/mod.rs`:
  status 0.

The first package-wide all-target Clippy attempt returned 101 for two concurrent
publication-test instrumentation findings outside this target (`unused_self` and
`map_unwrap_or`). Those changes belong to the publication increment; this page
does not relabel that earlier failed checkpoint. The final passing package-wide
check is listed separately above.
No whole-workspace runtime result or historical scientific count is implied by
this target. Final integration verification remains a separate gate.

`parse_genbank`/`parse_embl` expose string parsing, not public buffered-stream
constructors. Their generated grammar test therefore does not claim injected
read/invalid-UTF-8 transport coverage; separate file/profile tests retain their
existing scope. Gzip may release an earlier payload prefix before discovering a
bad trailer; only BGZF promises whole-member verification before release. EOF
must be consumed to validate unread transport suffixes.

This finite corpus is not exhaustive fuzzing or a proof of memory safety,
optimality, numerical fidelity, conformance to every format dialect, resistance
to every decompression bomb or concurrent mutation, or native compressed-byte
identity. It does not reopen closed scientific evidence. Larger fuzz campaigns,
coverage measurement, mutation testing, OS/CPU matrices and performance workloads
remain separately scoped testing layers.
