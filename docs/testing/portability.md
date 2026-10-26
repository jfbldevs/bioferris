# Numerical, SIMD and thread portability assurance

## Contents

- [Capability and scope index](#capability-and-scope-index)
- [Use the safe extension API](#use-the-safe-extension-api)
- [Explicit backend inventory](#explicit-backend-inventory)
- [Reproduce the local checks](#reproduce-the-local-checks)
- [Observed host and executed lanes](#observed-host-and-executed-lanes)
- [Numerical reference portability](#numerical-reference-portability)
- [Thread determinism and limits](#thread-determinism-and-limits)

[API index](../README.md#api-index) · [Contract/scientific assurance](contracts-and-science.md) ·
[Kernel source and tests](../../align-fast-kernels/src/lib.rs)

## Capability and scope index

| Capability | Rust access |
|---|---|
| Checked wavefront byte extension | `align_fast_kernels::WavefrontExtension` safe slice interface; alignment facades use this internal workspace crate |
| Explicit Portable/AVX2 assurance | Test-only construction of private backend kinds, after runtime AVX2 detection; no backend-selection production API |
| Exact retained numerical profiles | Existing public protein/motif APIs and offline reference-bit tests; no reference runtime is required |
| Ordered thread determinism | Existing `Pairwise` batch APIs with 1/2/4 requested threads |

This increment adds tests and documentation only. It changes no production
computation, dependency, toolchain, numeric tolerance or reference snapshot.
It neither reruns nor reads the closed MSA evaluation.

## Use the safe extension API

Literal inputs are eight bytes `ACGTACGT` and `ACGTTTGT`. At diagonal zero,
the first four bytes agree and the fifth differs. This complete Rust program
uses the normal runtime-selected safe API; backend forcing remains test-only:

<!-- portability-rust-program -->
```rust
use align_fast_kernels::WavefrontExtension;

fn main() {
    let kernel = WavefrontExtension::new();
    let mut offsets = [0];
    let result = kernel.extend_with_matches(b"ACGTACGT", b"ACGTTTGT", 0, &mut offsets);
    assert_eq!(offsets, [4]);
    assert_eq!(result.matches(), 4);
    println!("Offsets: {offsets:?}");
    println!("Diagonal range: {:?}", result.range().map(|range| (range.lower(), range.upper())));
    println!("Matched bytes: {}", result.matches());

    let mut invalid = [-1, 99];
    let result = kernel.extend_with_matches(b"ACGTACGT", b"ACGTTTGT", 0, &mut invalid);
    assert_eq!(invalid, [-1, 99]);
    assert_eq!(result.matches(), 0);
    println!("Invalid offsets retained: {invalid:?}");
}
```

Run from the repository root, with the pinned toolchain and dependencies cached:

```bash
cargo build -p align-fast-kernels --lib --offline --locked --jobs 2
rustdoc --test docs/testing/portability.md --edition 2024 --extern align_fast_kernels=target/debug/libalign_fast_kernels.rlib -L dependency=target/debug/deps -D warnings --test-args=--nocapture
```

The program's complete stdout is below; exit 0, empty stderr. Rustdoc additionally
prints its own harness summary. The kernel unit test also reads this output fence.

<!-- portability-rust-output -->
```text
Offsets: [4]
Diagonal range: Some((0, 0))
Matched bytes: 4
Invalid offsets retained: [-1, 99]
```

Offsets are byte positions, not biological sequence coordinates or scores. An
unreachable negative offset and an out-of-input offset remain unchanged; they do
not create matching bytes. This low-level interface does not return a biological
validation error. Its reported diagonal span includes representable nonnegative
entries even if their input coordinates are malformed, so that span alone is
not proof of a legal alignment. Ordinary callers should use the typed alignment
facade rather than invent wavefront state.

## Explicit backend inventory

The existing tests compared runtime dispatch with the portable implementation.
On an AVX2 machine that does not directly establish forced-Portable execution.
The new `explicit_portability_tests` module instead constructs Portable explicitly
and compares it against a separate byte-at-a-time scalar oracle. That oracle
uses wider coordinate arithmetic and checked byte lookup, not production word
loads, vector intrinsics or extension helpers.

| Family | Cases per backend | Boundaries exercised |
|---|---:|---|
| Dimensions/offsets | 1,904 | Input lengths 0..65 around 4/8/16/32/64-byte edges; 0/1/7/8/9/15/16/17 offsets; negative/unreachable/outside entries; `i32` extreme diagonals and overflowing intervals |
| Long unaligned/tails | 432 | Lengths 255/256/257 and 8,191/8,192/8,193/8,257/16,385; independently shifted slices; full match or mismatches at word/vector/tail boundaries; 17 offsets to include a scalar tail |
| Seeded arbitrary bytes | 2,048 | Fixed seed `0xdaf1_a55e_5eed_0123`, arbitrary bytes, unequal input lengths, mixed offsets and up to 40 lanes; reproducible without a new RNG dependency |
| Every word-byte difference | 2,040 | Each of eight mismatch positions and every nonzero byte XOR difference; checked against independent byte walking |

All 6,424 cases compare updated offsets, diagonal range and aggregate match count.
Both counting and non-counting public entry points are exercised. The long cases
straddle the existing 8,192-byte contiguous-AVX2 threshold; extreme-invalid
intervals also test the dispatcher's portable fallback guards.

On x86-64, the explicit AVX2 test first calls
`std::arch::is_x86_feature_detected!("avx2")`. Only then does it construct the
private AVX2 kind. Each AVX2 case is compared with both the independent scalar
expectation and explicitly forced Portable. No new `unsafe` code is introduced.
When runtime capability is absent, the test prints
`AVX2 lane: NOT EXECUTED (CPU/OS capability unavailable)` and returns; a passing
harness is not evidence that the unavailable lane ran. On other architectures,
the x86-64 test is not compiled and only Portable can be exercised locally.

## Reproduce the local checks

From the repository root:

```bash
cargo test -p align-fast-kernels --offline --locked --jobs 2 -- --nocapture
cargo test -p align-fast-kernels --release --offline --locked --jobs 2 -- --nocapture
cargo test -p bioferris --test protein_properties --test motifs --test scientific_assurance --offline --locked --jobs 2
cargo clippy -p align-fast-kernels --all-targets --offline --locked --jobs 2 -- -D warnings
```

The lane messages are explicit observed stdout records, not a native-compatibility
percentage. Success returns exit 0. A mismatched offset, range, match count,
reference bit or complete batch result fails an assertion and returns nonzero;
do not silently select a different backend, skip the case or update expectations.
Cargo build diagnostics use stderr and are separate from program stdout.

## Observed host and executed lanes

Local observation date: 2026-10-18.

| Property | Observed value |
|---|---|
| Compiler | `rustc 1.85.0 (4d91de4e4 2025-02-17)` |
| Compiler host / edition | `x86_64-unknown-linux-gnu` / 2024 |
| LLVM | 19.1.7 |
| OS/kernel | Linux `6.18.33.2-microsoft-standard-WSL2`, x86-64 |
| C library | Ubuntu glibc 2.39 |
| CPU | Intel Core i9-14900HX |
| Runtime AVX2 capability | Present; the guarded AVX2 lane actually executed |

Portable and guarded AVX2 execution on this same machine are distinct backend
checks, not two operating-system or architecture qualifications. Linux ARM,
macOS, Windows, non-AVX2 x86-64 and other compiler/libm combinations were **not
executed locally**. A configured future CI lane is not a passing remote result.
The parent verification record must retain actual commands, profiles, lane
messages, outcomes and any unexecuted gates.

## Numerical reference portability

Existing protein observations pin Biopython 1.88/CPython 3.12.3 on Linux x86-64,
including 7,344 retained binary64 values. Existing motif observations additionally
pin NumPy 2.3.3 and the recorded WSL2/glibc environment. See
[protein provenance](../../seq-bio/tests/fixtures/protein-properties-provenance.md)
and [motif provenance](../../seq-bio/tests/fixtures/motifs-provenance.md).

The existing tests compare retained bits exactly, not rounded displays or an
epsilon. Motif `ln` and protein `powf` can depend on platform math behavior;
exactness established on the recorded lane does not establish every platform.
The extension tests here operate on integer coordinates and raw byte equality;
they do not prove cross-platform floating-point identity. Preserve reference
snapshots and exact assertions. Any portability discrepancy needs a separately
reviewed platform/model decision, not a wider tolerance or a silent rebaseline.

## Thread determinism and limits

The new [scientific suite](../../seq-bio/tests/scientific_assurance.rs) submits
32 ordered small pairs through Needleman-Wunsch, Smith-Waterman and global Gotoh
with 1, 2 and 4 requested threads. It compares complete results, including ranges
and traceback operations, not just scores, and independently checks the 96
one-thread baseline traces. This is reproducibility through existing supported
APIs, not proof that a particular scheduler ran workers simultaneously, a speedup
measurement or a new thread-safety guarantee for all other algorithms.

All new probes are finite synthetic cases. No enormous allocation is required to
exercise integer-coordinate extremes. Safe test code does not establish absence
of every memory-safety bug in the pre-existing intrinsic implementation. Fuzzing,
sanitizers, broader CPUs, alternate operating systems and numerical environments
remain separate assurance lanes. No scientific kernel default, sealed campaign,
reference artifact, public backend selector or runtime dependency is changed.
