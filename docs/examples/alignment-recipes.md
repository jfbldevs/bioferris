# Alignment recipes with verified input and output

## Contents

- [API and recipe index](#api-and-recipe-index)
- [Run the examples](#run-the-examples)
- [Global, local and affine pairwise alignment](#global-local-and-affine-pairwise-alignment)
- [Edit distance, approximate search and affine costs](#edit-distance-approximate-search-and-affine-costs)
- [Seeded X-drop extension](#seeded-x-drop-extension)
- [Clustal Omega, MUSCLE 5 and MAFFT](#clustal-omega-muscle-5-and-mafft)
- [POASTA graph construction and GFA](#poasta-graph-construction-and-gfa)
- [FASTA round-trip and alignment export](#fasta-round-trip-and-alignment-export)
- [Combine distance, MSA and a graph query](#combine-distance-msa-and-a-graph-query)
- [Errors and compatibility limits](#errors-and-compatibility-limits)
- [Verify the documented outputs](#verify-the-documented-outputs)

## API and recipe index

| Workflow | Rust entry points | Recipe |
|---|---|---|
| Global linear, local linear and global affine alignment | `Pairwise::needleman_wunsch()`, `smith_waterman()`, `gotoh()` | `pairwise` |
| Edit distance and approximate search | `Pairwise::myers()` | `edit-search` |
| Affine cost minimization and adaptive dispatch | `Pairwise::wfa()`, `adaptive_affine()` | `edit-search` |
| Seeded heuristic extension | `Pairwise::xdrop()` | `extension` |
| Multiple sequence alignment | `Multiple::clustal_omega()`, `muscle5()`, `mafft(strategy)` | `multiple` |
| Partial-order graph alignment | `Graph::poasta()` | `graph` |
| FASTA and alignment formats | `parse_fasta`, `bioferris::io::*_to`, run writers | `formats` |

This page exercises all ten alignment families. Adaptive affine dispatch reuses
existing exact kernels; it is not an eleventh family. See the
[full API capability index](../README.md#api-index) for biological operations.

## Run the examples

Run each command below from the repository root. They use the pinned
Rust toolchain and local dependency cache. If the cache has not been populated,
run `cargo fetch --locked` once before using `--offline`.

The source is
[`seq-bio/examples/alignment_recipes.rs`](../../seq-bio/examples/alignment_recipes.rs).
Its six recipe names choose the literal inputs documented here; they are example
program arguments. The program
computes every alignment and serialization through the existing API; it does not
print saved scientific answers. Each successful command exits with status 0,
writes the shown stdout with a final newline, and creates no files.

The output blocks below are complete captured stdout, including whitespace.
They are checked against freshly computed results by the tests at the end of
this page. Small examples demonstrate API behavior; they do not establish
universal reference equivalence, biological optimality or performance.

Pairwise `target` is the first input and `query` is the second. All displayed
ranges use Rust's zero-based, half-open convention: `2..9` includes indices 2
through 8. Extended CIGAR uses `=` for an exact match, `X` for a substitution,
`D` for a target symbol opposite a query gap, and `I` for a target gap opposite
a query symbol. Higher scores are better; smaller distances and penalties are
better. These are different objectives and must not be compared as one scale.

## Global, local and affine pairwise alignment

Input sequences:

```text
target: GATTACA
query:  GCATGCT

additional affine-gap input:
target: AACCCGG
query:  AAGG
```

The [`pairwise` source function](../../seq-bio/examples/alignment_recipes.rs)
calls `.run(&target, &query)` on each builder, then obtains the score, ranges,
CIGAR and aligned strings from the typed run. Needleman-Wunsch and
Smith-Waterman use `.scores(2, -1).gap(-2)`; Gotoh uses
`.scores(2, -1).gaps(-5, -1).global()`.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- pairwise
```

Expected stdout:

<!-- output: pairwise -->
```text
Needleman-Wunsch: score=2
target=0..7 query=0..7 cigar=1=2X1=1X1=1X
target GATTACA
query  GCATGCT
Smith-Waterman: score=4
target=1..3 query=2..4 cigar=2=
target AT
query  AT
Gotoh: score=2
target=0..7 query=0..7 cigar=1=2X1=1X1=1X
target GATTACA
query  GCATGCT
[affine gap example]
Gotoh: score=1
target=0..7 query=0..4 cigar=2=3D2=
target AACCCGG
query  AA---GG
```

The global alignments consume both complete sequences. Smith-Waterman selects
the local `AT` segment, so its ranges are shorter and start at different input
positions. A tie policy selects one best local result; this is not a list of
every optimal alignment. The final Gotoh example has four matches and one
length-three deletion: `4 * 2 + (-5 + (3 - 1) * -1) = 1`. Its affine convention
is `open + (length - 1) * extension`.

See the [pairwise modes](../algorithms/affine-alignment-modes.md) and
[output contract](../api/output-contract-v1.md) for additional modes and reports.

## Edit distance, approximate search and affine costs

Inputs:

```text
pattern:    GATTACA
comparison: GACTACA
search text: TTTGATTACATTTGACTACA
maximum search distance: 1
WFA/adaptive costs: mismatch=4, opening=6, extension=2
```

The [`edit_search` source function](../../seq-bio/examples/alignment_recipes.rs)
uses `Pairwise::myers().distance(...)`, `.run(...)` for an edit traceback,
and `.find_all(&pattern, &text, 1)` for endpoint search. WFA and adaptive affine
use `.penalties(4, 6, 2).run(...)` on the same pattern/comparison pair.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- edit-search
```

Expected stdout:

<!-- output: edit-search -->
```text
global_distance=1
Myers bit-vector: distance=1
target=0..7 query=0..7 cigar=2=1X4=
target GATTACA
query  GACTACA
search: maximum_distance=1
text_end=9 distance=1
text_end=10 distance=0
text_end=11 distance=1
text_end=20 distance=1
WFA: penalty=4
target=0..7 query=0..7 cigar=2=1X4=
target GATTACA
query  GACTACA
Adaptive exact affine: penalty=4
target=0..7 query=0..7 cigar=2=1X4=
target GATTACA
query  GACTACA
```

Myers counts one unit-cost substitution. WFA charges the configured mismatch
cost 4; a gap of length `k` would cost `opening + k * extension`, different
from the signed Gotoh scoring convention above. Adaptive affine chooses an
existing exact backend under the same cost contract; this example makes no
speed claim.

Search results are **exclusive text endpoints**, not start/end intervals or
nonoverlapping occurrences. Endpoint 10 is the exact match ending after index 9;
the same nearby region also has accepted endpoints 9 and 11 at distance 1.
The API returns the best edit distance to a suffix ending at each accepted
endpoint, permits overlap, and does not return a starting coordinate. Do not
infer a start as `text_end - pattern.len()` when indels are possible.

See the [Myers contract](../algorithms/myers-bit-vector.md) and
[WFA cost contract](../algorithms/wfa.md).

## Seeded X-drop extension

Inputs:

```text
target: TTGATTACAGG
query:  CCGACTACATT
seed: target index 5, query index 5 (the sixth symbol, T, in each)
direction: Both
match=2, mismatch=-3, gap open=-6, gap extension=-1, X-drop threshold=8
```

The [`extension` source function](../../seq-bio/examples/alignment_recipes.rs)
uses `Pairwise::xdrop().scores(2, -3).gaps(-6, -1).threshold(8)` and
`.run_from_seed(&target, &query, XDropSeed::new(5, 5), XDropDirection::Both)`.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- extension
```

Expected stdout:

<!-- output: extension -->
```text
X-drop: score=9
target=2..9 query=2..9 cigar=2=1X4=
target GATTACA
query  GACTACA
```

The forced seed is included in an extension in both directions. Six matches and
one mismatch score `6 * 2 - 3 = 9`; the unrelated terminal symbols remain
outside the returned ranges. X-drop prunes paths using the threshold and seed:
it is not a guarantee of unrestricted Smith-Waterman optimality and does not
discover seeds or produce BLAST statistics. Seed positions are symbol indices,
not boundaries between symbols. See the [X-drop contract](../algorithms/x-drop.md).

## Clustal Omega, MUSCLE 5 and MAFFT

Input is typed protein FASTA with these complete records, in this order:

```text
>human
MKTAYIAKQRQ
>mouse
MKTAYIAKQKQ
>fish
MKTAYIAQRQ
```

The [`multiple` source function](../../seq-bio/examples/alignment_recipes.rs)
parses those records with `parse_fasta::<Protein>`, then calls
`Multiple::clustal_omega().threads(1).run(&records)`,
`Multiple::muscle5().run(&records)`, and
`Multiple::mafft(MafftStrategy::FftNs1).run(&records)`.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- multiple
```

Expected stdout:

<!-- output: multiple -->
```text
Clustal Omega
profile=CLUSTALO-1.0.2-default-protein-v1
rows=3 columns=11
human MKTAYIAKQRQ
mouse MKTAYIAKQKQ
fish MKTAYIAQRQ-
MUSCLE 5
profile=MUSCLE-5.1.0-Nature-2022-default-v2
rows=3 columns=11
human MKTAYIAKQRQ
mouse MKTAYIAKQKQ
fish MKTAYIA-QRQ
MAFFT
profile=MAFFT-2002-FFT-NS-1-v1
rows=3 columns=11
human MKTAYIAKQRQ
mouse MKTAYIAKQKQ
fish MKTAYIA-QRQ
```

Three input sequences become a rectangular alignment with 11 columns. `-` is
an inserted gap; removing gaps recovers each original sequence. The algorithms
place the fish gap differently here. This is an observed difference between
the selected scientific profiles, not a claim that the methods are equivalent.
These illustrative names do not identify externally sourced biological records.

These constructors select the retained publication profiles shown in stdout.
The MUSCLE 5 profile is at version 2 since 2026-10-22, which changed only how
refinement compares scores; this alignment is the same as under version 1
([version history](../algorithms/muscle5.md#version-history)).
In particular, `FftNs1` is the historical 2002 MAFFT profile; it does not mean
native MAFFT 7.526 output. The six historical MAFFT strategies and separately
named `mafft_7526_*` compatibility constructors have different contracts. Their
evidence is described in the [selection reference](../api/rust-api-v1.md#public-selection-surface)
and [MAFFT compatibility specification](../algorithms/mafft-7526-compatibility.md).
The short example does not extend those claims to other versions or options.
Historical FFT-NS-1 does not accept `.threads(...)`, even for 1: the facade
returns `UnsupportedOption` instead of pretending to parallelize it.

## POASTA graph construction and GFA

Input DNA records, in progressive insertion order:

```text
>reference
acgt acgt
>allele
ACGTTCGT
>deletion
ACGTCGT
```

Parsing removes sequence-line whitespace and normalizes case. The
[`graph` source function](../../seq-bio/examples/alignment_recipes.rs) calls
`Graph::poasta().costs(4, 6, 2).run(&records)`, then passes the authoritative
`run.result().graph()` to `bioferris::io::write_poa_gfa_to`.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- graph
```

Expected stdout (GFA fields are separated by literal tab characters):

<!-- output: graph -->
```text
POASTA
profile=LEE-2002-POASTA-2025-global-affine-ASTAR-DFS-superbubble-iterative-v1
rows=3 columns=8
reference ACGTACGT
allele ACGTTCGT
deletion ACGT-CGT
[GFA]
H	VN:Z:1.0
S	n0	A	CL:i:0
S	n1	C	CL:i:1
S	n2	G	CL:i:2
S	n3	T	CL:i:3
S	n4	A	CL:i:4
S	n5	C	CL:i:5
S	n6	G	CL:i:6
S	n7	T	CL:i:7
S	n8	T	CL:i:4
L	n0	+	n1	+	0M	RC:i:3
L	n1	+	n2	+	0M	RC:i:3
L	n2	+	n3	+	0M	RC:i:3
L	n3	+	n4	+	0M	RC:i:1
L	n3	+	n5	+	0M	RC:i:1
L	n3	+	n8	+	0M	RC:i:1
L	n4	+	n5	+	0M	RC:i:1
L	n5	+	n6	+	0M	RC:i:3
L	n6	+	n7	+	0M	RC:i:3
L	n8	+	n5	+	0M	RC:i:1
P	reference	n0+,n1+,n2+,n3+,n4+,n5+,n6+,n7+	*
P	allele	n0+,n1+,n2+,n3+,n8+,n5+,n6+,n7+	*
P	deletion	n0+,n1+,n2+,n3+,n5+,n6+,n7+	*
```

The rectangular alignment summarizes the graph's columns. GFA retains the
alternative `A`/`T` nodes and the deletion path explicitly: `S` defines segments,
`L` defines links, and `P` records the original sequence paths. `n0` etc. are
zero-based node identifiers, not sequence positions in every path. `CL` retains
the alignment class and `RC` the link support count. This graph is not a
phylogenetic tree. POASTA minimizes nonnegative affine costs using
`opening + length * extension`; the recipe does not collapse progressive
insertions into an invented single MSA score. Record order is part of the
progressive workflow. See the [POASTA contract](../algorithms/poasta.md).

## FASTA round-trip and alignment export

The FASTA input is the same three DNA records shown above, including the
lowercase, space-containing `acgt acgt` line. The separate JSON distance input
is `GATTACA` versus `GACTACA`.

The [`formats` source function](../../seq-bio/examples/alignment_recipes.rs)
reads records with `parse_fasta::<Dna>`, writes them with
`bioferris::io::write_fasta_to`, reparses them and compares the typed records.
It runs historical MAFFT FFT-NS-1 once and gives that same alignment to
`write_phylip_to`, `write_nexus_to` and `write_stockholm_to`. Finally it calls
`Pairwise::myers().run_distance(...).write_json_to(...)`. All writers use
in-memory byte streams; the displayed labels are added by the example, not
part of the individual file formats.

```bash
cargo run --quiet --offline --locked -p bioferris --example alignment_recipes -- formats
```

Expected stdout:

<!-- output: formats -->
```text
records_equal_after_round_trip=true
[normalized FASTA]
>reference
ACGTACGT
>allele
ACGTTCGT
>deletion
ACGTCGT
[PHYLIP]
3 8
reference ACGTACGT
allele ACGTTCGT
deletion ACGT-CGT
[NEXUS]
#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=-;
  MATRIX
    'reference' ACGTACGT
    'allele' ACGTTCGT
    'deletion' ACGT-CGT
  ;
END;
[Stockholm]
# STOCKHOLM 1.0
reference ACGTACGT
allele    ACGTTCGT
deletion  ACGT-CGT
//
[Myers distance JSON]
{
  "schema": "align-output-v1",
  "kind": "myers_distance",
  "algorithm": "Myers bit-vector",
  "mode": "global unit edit distance",
  "alphabet": "DNA",
  "pattern": {"length": 7, "sequence": "GATTACA"},
  "text": {"length": 7, "sequence": "GACTACA"},
  "distance": 1
}
```

`records_equal_after_round_trip=true` means typed identifiers and sequences
survive this FASTA round-trip. It does not mean the original bytes are retained:
case, whitespace and wrapping are normalized. The current parser is eager and
keeps the complete trimmed FASTA header as the identifier. It has no separate
description, annotation or quality fields. `write_fasta_to` wraps at 80 symbols.

PHYLIP defaults to relaxed sequential output; its first line is `rows columns`.
NEXUS chooses `DATATYPE=DNA` from the Rust type. Stockholm carries the same
gapped rows, without annotations that the record model does not contain. These
are alignment writers, not inferred-tree outputs or demonstrated alignment-file
readers. A changed file extension does not perform a format conversion.

Native facade JSON uses `schema: "align-output-v1"`. That is distinct from the
`schema_version: 1` success/error envelope of the removed `bioferris` CLI.
This distance-only result has no fabricated traceback or coordinates. For file
output, the facade also provides methods such as `run.write_phylip(path)` and
`run.write_json(path)`; stream writers avoid choosing a filesystem destination
in these examples. See the [format and JSON contracts](../api/phylogenetic-output-v1.md).

## Combine distance, MSA and a graph query

[`alignment_workflow.rs`](../../seq-bio/examples/alignment_workflow.rs) composes
Myers distance, MAFFT and a POASTA query through the existing facades in one
ordinary Rust program. Literal inputs are DNA `ACGT` and `AGT`, named
`reference` and `deletion` for the collection workflows.

<!-- alignment-rust-program:start -->
```rust
//! Myers distance, a MAFFT alignment and a POASTA graph query in one Rust program.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let target = dna("ACGT")?;
    let query = dna("AGT")?;
    let distance = Pairwise::myers().distance(&target, &query)?;
    println!("Edit distance: {distance}");

    let records = vec![
        SequenceRecord::new("reference", target)?,
        SequenceRecord::new("deletion", query.clone())?,
    ];
    let multiple = Multiple::mafft(MafftStrategy::FftNs1).align(&records)?;
    println!("MSA records: {}", multiple.alignment().sequence_count());

    let graph = Graph::poasta().align(&records)?;
    let aligned_query = Graph::poasta().run_query(graph.graph(), &query)?;
    println!("Graph paths: {}", graph.graph().path_count());
    println!("{}", aligned_query.report());
    Ok(())
}
```
<!-- alignment-rust-program:end -->

From the repository root with its pinned compiler and cached dependencies:

<!-- alignment-rust-command:start -->
```bash
cargo bio alignment_workflow --offline --locked
```
<!-- alignment-rust-command:end -->

Observed stdout, exit 0 and empty stderr:

<!-- alignment-rust-output:start -->
```text
Edit distance: 1
MSA records: 2
Graph paths: 2
Algorithm            : POASTA
Profile              : POASTA-2025-global-affine-ASTAR-DFS-superbubble-v1
Mode                 : global affine sequence-to-graph
Scoring              : mismatch=4 gap_open=6 gap_extend=2 (positive costs)
Score (negated cost) : 0
Backend              : Sparse
Graph                : nodes=4 edges=4 classes=4 paths=2
Query length         : 3
Trace columns        : 3
Aligned pairs        : 3
Trace                : 3A

Graph  n0  AGT  n3
           |||
Query   1  AGT   3

Legend: | exact match, . substitution, blank graph/query gap; nX is a stable graph node
```
<!-- alignment-rust-output:end -->

One deletion gives edit distance 1. The graph has two source paths; querying an
existing path costs zero and does not add a third path. Report query coordinates
are human-oriented, while `n0`/`n3` are graph node labels, not sequence offsets.
Invalid DNA such as `AC-GT` returns a typed sequence error before alignment.
This program was previously documented with the removed alignment CLI
(decision 029).

## Errors and compatibility limits

All computational calls propagate typed errors with `?`; callers can inspect
the original Rust error instead of scraping a report. Invalid example names
or missing/extra arguments fail with a usage message and nonzero exit status.

| Input or configuration problem | Relevant behavior |
|---|---|
| Invalid typed symbols, such as `-` in a raw DNA sequence | `SequenceError`; gaps belong in alignment results, not unaligned inputs |
| Positive signed gap scores or invalid cost parameters | Structured scoring/algorithm errors; validate parameters rather than changing their sign implicitly |
| X-drop seed index outside either sequence | X-drop error; indices identify real symbols |
| Empty/duplicate record collections or invalid record identifiers | Record/algorithm validation errors under the selected contract |
| `.threads(...)` on historical MAFFT FFT-NS-1 | `UnsupportedOption`; the example deliberately leaves it unset |
| Whitespace in PHYLIP or Stockholm identifiers | Format error; names are not silently rewritten; NEXUS can quote ordinary spaces |
| More than ten identifier bytes in strict PHYLIP | Format error rather than truncation; the example uses relaxed PHYLIP |
| I/O failure | Writer error; no successful export is claimed |

Not every family accepts every parameter or molecule/profile combination. Use
the [Rust API reference](../api/rust-api-v1.md) and individual specifications for
supported options. None of these examples calls a native reference executable,
opens the reserved MSA evaluation corpus, or changes the existing scientific
defaults. Historical fidelity evidence retains its original scope.

## Verify the documented outputs

The example contains six tests. Each executes a recipe, locates its named
stdout fence in this Markdown file and compares the complete computed string,
including tabs, spaces and the final newline. Changing an algorithm result or
documented output independently makes the corresponding test fail.

```bash
cargo test --offline --locked -p bioferris --example alignment_recipes
cargo clippy --offline --locked -p bioferris --example alignment_recipes -- -D warnings
```

[`alignment_workflow_guide.rs`](../../seq-bio/tests/alignment_workflow_guide.rs)
also compares the combined program above with its source file and runs the real
`cargo bio alignment_workflow` alias against the documented output.

These checks validate the documented examples and their use of existing APIs.
They do not replace scientific differential tests or reopen a closed
qualification campaign.
