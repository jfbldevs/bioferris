# Ergonomic Rust API v1

## Contents

- [Full API index](../README.md#api-index)
- [Runnable alignment recipes with inputs and outputs](../examples/alignment-recipes.md)
- [Scope](#scope)
- [Getting started](#getting-started)
- [Public selection surface](#public-selection-surface)
- [Choose a multiple-alignment profile](#choose-a-multiple-alignment-profile)
- [Convenience defaults](#convenience-defaults)
- [Results](#results)
- [Human-readable reports](#human-readable-reports)
- [FASTA](#fasta)
- [Threads and unsupported options](#threads-and-unsupported-options)
- [Advanced access](#advanced-access)
- [Validation contract](#validation-contract)

Status: implemented as the pre-1.0 `align` facade.

For complete input sequences, runnable programs and captured results across all
ten families, see the [Rust alignment recipes](../examples/alignment-recipes.md).

Acceptance evidence is recorded in the
Rust API v1 acceptance matrix.

## Scope

`align` is the small Rust-native entry point for the frozen ten-family algorithm
portfolio. It does not reimplement a recurrence, translate an objective, invoke an
external executable, or require Python. The detailed algorithm crates remain
available for callers that need workspaces, memory models, verification evidence,
or non-default scientific parameters.

A future Python package must call this facade or an equally narrow Rust adapter. It
must not duplicate algorithm or validation logic.

## Getting started

During workspace development:

```toml
[dependencies]
align = { path = "align" }
```

The common import is intentionally small:

```rust
use align::prelude::*;
```

Typed alphabets keep invalid molecules out of every kernel:

```rust
let target = DnaSequence::try_from("GATTACA")?;
let query = DnaSequence::try_from("GCATGCT")?;

let run = Pairwise::gotoh()
    .scores(2, -1)
    .gaps(-5, -1)
    .run(&target, &query)?;

println!("{}", run.report().labels("Target", "Query"));
# Ok::<(), align::Error>(())
```

Configured builders are immutable values and can be reused:

```rust
# use align::prelude::*;
# let first = DnaSequence::try_from("ACGT")?;
# let second = DnaSequence::try_from("AGT")?;
let aligner = Pairwise::needleman_wunsch()
    .scores(2, -1)
    .gap(-2);

let one = aligner.align(&first, &second)?;
let batch = aligner
    .threads(2)
    .align_many([(&first, &second), (&second, &first)])?;
assert_eq!(one, batch[0]);
# Ok::<(), align::Error>(())
```

`score()` returns only the optimum. For Needleman--Wunsch and Smith--Waterman it
runs the [linear score-only kernels](../algorithms/linear-score-only.md), whose
memory is one score row of `second.len() + 1` values instead of a full traceback
matrix; the value equals `align()`'s score.

## Public selection surface

| Family | Facade constructor |
|---|---|
| Needleman--Wunsch | `Pairwise::needleman_wunsch()` |
| Smith--Waterman | `Pairwise::smith_waterman()` |
| Gotoh | `Pairwise::gotoh()` |
| Myers | `Pairwise::myers()` |
| WFA | `Pairwise::wfa()` |
| Adaptive exact affine | `Pairwise::adaptive_affine()` |
| X-drop | `Pairwise::xdrop()` |
| Clustal Omega | `Multiple::clustal_omega()` |
| MUSCLE 5 Align | `Multiple::muscle5()` |
| MUSCLE 5 Super5 | `Multiple::muscle5_super5()` |
| MAFFT historical profiles | `Multiple::mafft(strategy)` |
| MAFFT 7.526 L-INS-i compatibility | `Multiple::mafft_7526_linsi()` |
| MAFFT 7.526 protein E-INS-i compatibility | `Multiple::mafft_7526_einsi()` |
| POASTA | `Graph::poasta()` |

The six MAFFT values are `FftNs1`, `FftNs2`, `FftNsi`, `Ginsi`, `Linsi`, and
`Einsi`. `MafftStrategy::profile_identifier()` returns the exact frozen profile
identifier.

The two `mafft_7526_*` constructors are a post-v1 compatibility extension. They
delegate to the separately pinned, source-guided pure-Rust port and report that
provenance explicitly. They do not alter the frozen six-profile selector or its
acceptance status. L-INS-i accepts protein, DNA, and RNA records; modern E-INS-i is
protein-only because DNA/RNA have not passed the exact native-output gate. See the
[compatibility specification](../algorithms/mafft-7526-compatibility.md).

Gotoh exposes global, local, semiglobal, fitting, suffix/prefix overlap, and exact
four-end configuration without counting them as separate algorithm families.

## Choose a multiple-alignment profile

The multiple-alignment constructors are of two kinds:

- **Compatibility profiles**, `Multiple::mafft_7526_linsi()` and
  `Multiple::mafft_7526_einsi()`, run a pinned pure-Rust port of MAFFT 7.526
  L-INS-i and E-INS-i. Choose them when results should agree with current native
  MAFFT.
- **Scientific profiles**, `Multiple::clustal_omega()`, `Multiple::muscle5()`,
  `Multiple::muscle5_super5()` and `Multiple::mafft(strategy)`, are frozen
  reconstructions of the published methods: Clustal Omega, MUSCLE 5.1 Align and
  Super5, and the MAFFT strategies as published in 2002–2005. They follow the
  papers, not every later change in the native tools, and Clustal Omega's
  unpublished HHalign values are reconstructed. Choose them to study or cite
  those methods.

The publication benchmark,
measured at `e9b3e644` on one laptop, gives these v1 numbers: mean BAliBASE 3 sum
of pairs (SP) and PREFAB 4 Q over the cases both tools completed, and the ratio
of native to BioFerris single-CPU time on BAliBASE, for a typical case
(geometric mean) and in total (dominated by the largest cases). A ratio above 1
means BioFerris was faster.

| Constructor | Native tool | BAliBASE SP, BioFerris / native | PREFAB Q, BioFerris / native | Time ratio, typical / total |
|---|---|---:|---:|---:|
| `mafft_7526_linsi()` | MAFFT 7.526 L-INS-i | 0.848 / 0.848 | 0.825 / 0.833 | 1.07 / 0.52 |
| `mafft_7526_einsi()` | MAFFT 7.526 E-INS-i | 0.837 / 0.836 | 0.825 / 0.826 | 1.06 / 0.50 |
| `muscle5()` | MUSCLE 5.1 `-align` | 0.854 / 0.856 | 0.650 / 0.671 | 0.54 / 0.50 |
| `muscle5_super5()` | MUSCLE 5.1 `-super5` | 0.844 / 0.845 | 0.732 / 0.727 | 0.43 / 0.30 |
| `clustal_omega()` | Clustal Omega 1.2.4 | 0.748 / 0.806 | 0.769 / 0.778 | 1.17 / 1.57 |
| `mafft(MafftStrategy::FftNs1)` | MAFFT 7.526 FFT-NS-1 | 0.767 / 0.781 | 0.757 / 0.779 | 5.43 / 3.61 |
| `mafft(MafftStrategy::FftNs2)` | MAFFT 7.526 FFT-NS-2 | 0.784 / 0.793 | 0.767 / 0.779 | 4.32 / 2.75 |
| `mafft(MafftStrategy::FftNsi)` | MAFFT 7.526 FFT-NS-i | 0.784 / 0.830 | 0.749 / 0.800 | 1.35 / 0.45 |
| `mafft(MafftStrategy::Ginsi)` | MAFFT 7.526 G-INS-i | 0.788 / 0.857 | 0.744 / 0.829 | 0.69 / 0.10 |
| `mafft(MafftStrategy::Linsi)` | MAFFT 7.526 L-INS-i | 0.792 / 0.848 | 0.766 / 0.833 | 0.73 / 0.11 |
| `mafft(MafftStrategy::Einsi)` | MAFFT 7.526 E-INS-i, `--oldgenafpair` | 0.792 / 0.848 | 0.758 / 0.832 | 0.75 / 0.12 |

In short:

- The 7.526 profiles score as native MAFFT does on average, within 0.008, but
  are not byte-exact reproductions: they reproduced native output on 39–54% of
  the cases they were not developed on
  (exactness).
  They are as fast as native on a typical case and take about twice as long on
  the largest.
- The MUSCLE 5 profiles score within 0.002 SP of native MUSCLE 5.1 on the cases
  both completed, at 2 to 3.3 times its total time. Profile v1, which the
  benchmark measured, failed on 61 cases (`muscle5`) and 40 (`muscle5-super5`);
  [profile v2](../algorithms/muscle5.md#version-history) completes them (see the
  follow-up below).
- The Clustal Omega reconstruction scores 0.057 SP below native on BAliBASE,
  0.464 against 0.608 on RV11, and 0.126 Q below on HomFam. Prefer native
  `clustalo` where accuracy matters.
- The historical MAFFT strategies FFT-NS-1 and FFT-NS-2 come within 0.015 SP of
  native and are 2.8 and 3.6 times faster in total; the iterative ones are
  0.046 to 0.069 SP below native and up to ten times slower. Use them for the
  published methods, not for the most accurate alignment.

These are v1 measurements on one host, with bootstrap intervals, coverage and
per-set means on the results page. Since then, profile v2 completes every case,
and Super5 reuses repeated pair posteriors within a run, with the same output and
less time. The
follow-up benchmark
measured both MUSCLE 5 profiles again at `4166ec1e`, on the same host. Its scores
cover every complete case, 231 BAliBASE and 26 PREFAB, so its means are not
comparable with v1's, which left out the cases profile v1 failed on:

| Constructor | Native tool | BAliBASE SP, BioFerris / native | PREFAB Q, BioFerris / native | Time ratio, typical / total |
|---|---|---:|---:|---:|
| `muscle5()` | MUSCLE 5.1 `-align` | 0.869 / 0.872 | 0.812 / 0.831 | 0.56 / 0.45 |
| `muscle5_super5()` | MUSCLE 5.1 `-super5` | 0.854 / 0.855 | 0.814 / 0.812 | 0.77 / 0.47 |

No run failed, every score difference has a 95% interval that includes zero, and
both profiles take 2.1 to 2.2 times native total time on one CPU. They match
native accuracy on average; choose native MUSCLE 5.1 where speed matters.

## Convenience defaults

Defaults belong only to the facade and are visible through each builder's `Debug`
representation and getters. They do not rename or modify a scientific profile.

| Builder | Default |
|---|---|
| Linear pairwise | match `2`, mismatch `-1`, gap `-2` |
| Gotoh | global; match `2`, mismatch `-1`, open `-5`, extension `-1` |
| WFA | mismatch cost `4`, gap opening `6`, gap extension `2` |
| X-drop | match `2`, mismatch `-1`, open `-5`, extension `-1`, threshold `20` |
| Clustal Omega | molecule-specific frozen published profile |
| MUSCLE 5 | frozen published Align or Super5 profile selected by constructor |
| MAFFT | exact strategy selected by `MafftStrategy` |
| POASTA | mismatch cost `4`, gap opening `6`, gap extension `2` |

Gotoh/X-drop use `open + (length - 1) * extension`. WFA/POASTA use their published
positive-cost convention `opening + length * extension`. Their builders therefore
use `gaps(...)` for signed score contributions and `penalties(...)`/`costs(...)` for
positive costs.

## Results

The facade returns the native algorithm result, retaining all scientific evidence.
It does not coerce distance, score, and penalty into one misleading number.

- score-maximizing alignments expose `score()`;
- Myers exposes `distance()`;
- WFA exposes `penalty()`;
- algorithm-specific MSA and POASTA evidence remains available on the native result.

`PairwiseAlignmentExt` adds extended CIGAR, exact-match count, and identity to every
owned pairwise result. `MultipleAlignmentExt` adds aligned FASTA output to MSA
results. Facade `Run` values additionally expose direct human, JSON, alignment,
phylogenetic-matrix, tabular, and graph file methods. The required types and
extension traits are imported by `prelude::*`.

## Human-readable reports

Every public family has a `.run(...)` counterpart that invokes the same frozen
kernel as `.align(...)` while retaining the inputs and configuration needed to
explain the result. The native result remains available through `.result()` and
`.into_result()`.

```rust
# use align::prelude::*;
# let target = DnaSequence::try_from("GATTACA")?;
# let query = DnaSequence::try_from("GCATGCT")?;
let run = Pairwise::gotoh().run(&target, &query)?;
let report = run
    .report()
    .labels("Target", "Query")
    .block_width(60)?
    .verbosity(Verbosity::Detailed);
println!("{report}");
# Ok::<(), align::Error>(())
```

Pairwise reports show configuration, ranges, the correct objective name, identity,
substitution/indel statistics, extended CIGAR, coordinates, and a match line. Myers
also exposes `run_distance` and `run_find_all`; WFA says `Penalty`, and seeded X-drop
retains its seed, direction, stopping reason, and traversal counters.

Callers that do not need a traceback can still obtain a complete, unambiguous
summary. Linear and Gotoh builders expose `run_score`, WFA exposes `run_penalty`,
X-drop exposes `run_score` and `run_score_from_seed`, and Myers exposes
`run_distance`. These runs show the native objective and execution evidence but do
not fabricate coordinates, identity, CIGAR, or aligned rows.

```rust
# use align::prelude::*;
# let target = DnaSequence::try_from("GATTACA")?;
# let query = DnaSequence::try_from("GCATGCT")?;
let run = Pairwise::wfa().run_penalty(&target, &query)?;
assert_eq!(run.result().penalty(), 16);
println!("{}", run.report().verbosity(Verbosity::Detailed));
# Ok::<(), align::Error>(())
```

MSA reports show the exact profile identifier, rows, columns, conserved columns,
gap characters, independent row coordinates, and block conservation. Detailed mode
adds retained Clustal Omega, MUSCLE 5, or MAFFT evidence. `Graph::poasta().run(...)`
adds graph dimensions to its MSA, while `run_query(...)` renders the selected graph
path against the query and retains stable node identifiers.

`ReportOptions` controls block width, verbosity, and human versus native Rust
coordinates. Deterministic `Display` output never contains ANSI color escapes. The
precise definitions are frozen in the
[human-readable output contract](output-contract-v1.md).

Interchange writers remain separate from human text:

```rust,no_run
# use align::prelude::*;
# let records = read_fasta::<Dna>("input.fasta")?;
let run = Graph::poasta().run(&records)?;
run.write_text("alignment.txt")?;
run.write_json("alignment.json")?;
run.write_clustal("alignment.aln")?;
run.write_phylip("alignment.phy")?;
run.write_nexus("alignment.nex")?;
run.write_stockholm("alignment.sto")?;
run.write_gfa("alignment.gfa")?;
run.write_dot("alignment.dot")?;
# Ok::<(), align::Error>(())
```

Additional methods cover aligned pairwise FASTA and Myers TSV. Lower-level format
functions have path and `_to` stream variants in `align-formats`. PHYLIP supports
relaxed/strict and sequential/interleaved options; NEXUS selects its datatype from
the Rust alphabet. The versioned JSON shapes and exact phylogenetic-format rules
are documented in the
[phylogenetic and structured output contract](phylogenetic-output-v1.md).

## FASTA

This section documents the unchanged legacy reader. For metadata-preserving
streaming FASTA/Sanger FASTQ and gzip, see the additive
[`bioferris` file-I/O API](file-io-v1.md). For aligned inputs use the separate
[alignment readers](alignment-input-v1.md).

```rust,no_run
use align::prelude::*;

let records = read_fasta::<Protein>("proteins.fasta")?;
let run = Multiple::mafft(MafftStrategy::Linsi).run(&records)?;
run.write_fasta("aligned.fasta")?;
# Ok::<(), align::Error>(())
```

The parser retains the complete trimmed header as the identifier, ignores ASCII
whitespace in sequence lines, normalizes lowercase symbols through the selected
alphabet, and rejects malformed structure or invalid symbols with `FastaError`.

The workspace's runnable tutorials combine these
inputs with pairwise, MSA, Myers, POASTA, JSON, and phylogenetic-output examples.
Each tutorial program is independently executable and comments the complete data
flow.

## Threads and unsupported options

Every pairwise builder accepts `.threads(count)` for ordered parallel `align_many`;
the default is one worker and zero is rejected. Results retain input order. Clustal
Omega accepts `.threads(count)`, and MAFFT FFT-NS-i and L-INS-i accept the same
setting. L-INS-i evaluates independent fixed-tree refinement candidates in
parallel, then reduces them in canonical partition order; scores, ties, and output
are bit-exact across worker counts. Other retained MAFFT profiles reject `.threads(...)` with
`Error::UnsupportedOption`; an option is never silently ignored.

Both MAFFT 7.526 compatibility builders accept `.threads(count)`. Their local Rayon
worker count is deterministic in the pinned tests and BAliBASE validation (1 versus
4 workers, 12/12 exact outputs). It does not emulate native MAFFT's `--thread`
option, whose iterative scheduling can select a different alignment.

## Advanced access

Detailed types are grouped without polluting the beginner path:

```rust
use align::advanced::{core, mafft_compat, pairwise};
```

This exposes the original crates and their reusable workspaces, score-only kernels,
memory estimates, verification functions, evidence records, and complete parameter
objects. `mafft_compat` is visibly separate from the independent `mafft` module.

## Validation contract

Facade tests invoke every facade entry point and its frozen kernel directly on the
same typed inputs. They require exact equality of the canonical pairwise traceback,
native WFA/X-drop result, rectangular MSA, MAFFT profile identifier, POASTA graph,
and POASTA query alignment. Format tests cover multiline input, normalization,
malformed structure, and aligned/unaligned output. Doctests compile the beginner
path as part of the workspace gate.
