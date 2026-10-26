# MAFFT 7.526 compatibility profiles

Status: validated on the retained 2026-08-19 suites; expanded diagnostics on
2026-10-11 now pass after the M2 nucleotide pairwise gap-scale correction.
Exactness is corpus-scoped; the joint MSA fidelity gate remains closed.

## Purpose and boundary

`align-mafft-compat` provides modern MAFFT behavior for callers who need the
same alignment as the pinned native release rather than a reconstruction of an
ancestral publication. It is intentionally separate from `align-mafft`:

- `align-mafft` retains the six independently reconstructed and formally
  accepted historical profiles;
- `align-mafft-compat` delegates the numerical pipeline to the source-guided,
  pure-Rust `mafft = 0.1.2` port and adapts its result to this workspace's typed
  records, validation, errors, and reporting.

The compatibility path is pure Rust and runs in-process. It uses no C code,
FFI, subprocess, Python runtime, sidecar, or installed MAFFT executable. It is
nevertheless a delegated port, not a second independent rewrite from the
papers. Claims about the independent historical profiles must not be applied to
this crate.

## Versioned public profiles

| Profile | Rust entry point | Native invocation represented | Alphabet scope |
|---|---|---|---|
| `MAFFT-7.526-L-INS-i-localpair-v1` | `Mafft7526Linsi<A>` | `mafft --localpair --maxiterate 1000 --inputorder --quiet INPUT` | Protein, DNA, RNA |
| `MAFFT-7.526-E-INS-i-genafpair-protein-v1` | `Mafft7526Einsi` | `mafft --genafpair --ep 0 --maxiterate 1000 --inputorder --quiet INPUT` | Protein only |

The ergonomic equivalents are `Multiple::mafft_7526_linsi()` and
`Multiple::mafft_7526_einsi()`.

E-INS-i is deliberately protein-only. The pinned port matched native MAFFT on
the protein fixtures and BAliBASE corpus, but its generalized-affine DNA/RNA
route diverged on adversarial unequal-internal-region cases. Excluding those
types at compile time is more accurate than publishing an unqualified
compatibility claim.

## Adapter invariants

Before and after delegation, the adapter requires:

1. a nonempty record collection with unique identifiers;
2. input order retained exactly;
3. one rectangular output row for every input row;
4. only symbols valid for the selected Rust alphabet;
5. finite engine evidence; and
6. exact recovery of every ungapped source sequence.

The upstream port drops source content for some mixed empty/nonempty inputs.
The adapter therefore removes empty records before execution and reinserts them
as all-gap rows. This is a safe workspace extension, not a claim about native
MAFFT's empty-sequence behavior.

## Exact native validation

The comparator was the official Linux amd64 MAFFT v7.526 package, SHA-256
`b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de`.
The Rust engine is pinned to crates.io package `mafft = 0.1.2`, checksum
`1933406ecd145d4b61eaf5c6d4029b7212a0af879ad71a6dac0157d14ca722ed`.

One family from every BAliBASE 3 R1-R5 category was run through the native
commands above and the release Rust adapter. Ordered identifiers and complete
gapped rows were compared after only ASCII-case and dot-gap normalization.

| Profile | Exact native alignments | Exact at 1 vs 4 Rust workers | Macro SP | Macro TC |
|---|---:|---:|---:|---:|
| L-INS-i | 6 / 6 | 6 / 6 | 0.921134 | 0.701348 |
| protein E-INS-i | 6 / 6 | 6 / 6 | 0.912156 | 0.658429 |
| Combined | 12 / 12 | 12 / 12 | — | — |

The difficult RV30/BB30027 L-INS-i result now matches native MAFFT exactly and
scores SP `0.897608`, TC `0.818182`. The independent historical L-INS-i
reconstruction scores SP `0.606938`, TC `0.000000` on that case. This closes
the observed modern-output gap for the new compatibility profile without
changing or relabeling the historical profile.

The checked-in black-box tests add:

- seven exact L-INS-i cases across protein, DNA, and RNA; and
- one synthetic protein E-INS-i case selected because modern `--genafpair`
  differs from the ancestral `--oldgenafpair` route.

The complete case metrics, canonical alignment hashes, commands, artifact
identities, and worker comparisons are in
`../benchmarks/results/mafft-7526-rust-compat-balibase3-subset.json`.
The reproducible external-data runner is
[`../../align-mafft-compat/examples/balibase.rs`](../../align-mafft-compat/examples/balibase.rs).

## Determinism

For a fixed adapter version, input order, profile, and parameters, the Rust
result is deterministic in the tested scope. The adapter creates a local Rayon
pool with an explicit positive worker count, and all 12 BAliBASE outputs were
byte-identical between one and four workers. Unit tests also gate this property.

Native MAFFT's `--thread` setting is not the same API as the adapter's worker
count. Iterative native runs can select a different but valid alignment when an
explicit native thread count changes floating-point scheduling. The native
commands represented by these profile identifiers therefore omit `--thread`.
The adapter's `.threads(n)` controls its own deterministic execution and does
not promise to emulate every native threading configuration.

## Performance evidence

A matched Linux fresh-process campaign compares native MAFFT with one- and
four-worker Rust executions on all six pinned BAliBASE cases. All 12 profile/case
comparisons remain exact. Rust wins every paired case: geometric-mean speedups are
`3.86x`/`4.10x` for L-INS-i and `4.04x`/`4.40x` for protein E-INS-i at one/four
workers. The complete protocol, limitations, and portfolio context are in the
project-wide performance assessment,
with raw samples in
`mafft-7526-process-performance-linux.json`.

## What is and is not achieved

Within the pinned commands and validation domain, these are successful Rust
ports of the native behavior: they run without the original executable and
produce the same tested alignments. The scope does not establish equivalence
for every possible sequence collection, option, platform-specific tie, or
future MAFFT release. It also does not cover the rest of MAFFT's strategy and
option surface.

The initial separate MSA native-fidelity checkpoint
added a frozen diagnostic expansion: L-INS-i matched 17/19 cases, with divergent
DNA and RNA counterexamples; protein E-INS-i matches 5/5. All valid outputs
preserve their sources and repeat deterministically. A native-arithmetic flag
experiment changed none of the candidate logical outputs and was rejected.
These synthetic observations are not BAliBASE SP/TC or new performance results.
At that historical checkpoint the global native-fidelity gate remained closed;
those observations did not establish
that all protein E-INS-i inputs are exact.

Subsequent M2 fixes the delegated engine's nucleotide pairwise gap-penalty scale
in a versioned runtime fork. The same frozen corpus now yields L-INS-i **19/19**
and protein E-INS-i **5/5**, with no regressions, deterministic source-preserving
repeats, and all 24 outputs exact at one/four workers (96 fresh Rust runs).
Permanent regression fixtures retain the two previously failing native rows.
The source-guided checkpoint
records provenance, the first divergent stage and reproducible offline evidence.
Neither independent holdouts nor new performance measurements are used here.

The subsequent separately labeled native-compatibility research path passes the
complete frozen v1 qualification, including **16/16 once-only reserved MAFFT
inputs** (12 L-INS-i and four protein E-INS-i), with exact logical rows and cold
determinism. All nine actual final gates, 496 Python tests and 930 Rust test
functions pass. See the final evidence and declared scope.
This later result does not rewrite the historical M1/M2 checkpoints, promote
scientific defaults, qualify all possible inputs/options, or establish an
independent MAFFT kernel reconstruction. No new performance campaign is run.

Accordingly, the accurate project-level statement is:

> The workspace contains six independent historical MAFFT reconstructions and
> a separate pure-Rust compatibility adapter that exactly matches native MAFFT
> 7.526 L-INS-i and protein E-INS-i on its pinned reference suites.

It would still be inaccurate to call every implementation here a perfect or
universal replacement for all MAFFT versions and modes.

## Dependency and license

The adapter exactly pins `mafft = 0.1.2` and `rayon = 1.11.0`. The Rust port is
licensed `MIT AND BSD-3-Clause`; both licenses are permitted by `deny.toml`.
The adapter code itself remains `MIT OR Apache-2.0`. Distributors must preserve
the notices and terms of all transitive dependencies in the usual way.

The workspace now patches only `mafft-core 0.1.2` to `vendor/mafft-core`.
Its original runtime file hashes, exact author revision and local changes are
recorded in `PROVENANCE.json`; both upstream license notices are included.
The registry cache and independent historical MAFFT implementation are unchanged.
