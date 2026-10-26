# Phylogenetic and structured output v1

## Contents

- [Full API index](../README.md#api-index)
- [Worked file and alignment export examples](../examples/alignment-recipes.md)
- [Recommended format by use](#recommended-format-by-use)
- [Facade example](#facade-example)
- [PHYLIP](#phylip)
- [NEXUS](#nexus)
- [Stockholm 1.0](#stockholm-10)
- [JSON schema](#json-schema-align-output-v1)
- [Why there is no Newick writer yet](#why-there-is-no-newick-writer-yet)
- [Compatibility rules](#compatibility-rules)

Status: implementation contract for the Rust-native `align` facade.

This layer serializes alignments for downstream phylogenetic programs. It does
not infer a phylogenetic tree, attach an evolutionary model, or reinterpret an
MSA guide tree as phylogenetic evidence.

Additive readers now live in `bioferris`:
[alignment import profiles](alignment-input-v1.md) and
[graph/result import profiles](graph-result-input-v1.md). Their documented
subsets and loss checks are separate from these unchanged export contracts.

The removed alignment CLI export adapter
reused these writers with create-only named publication and additional
representability guards: its CLUSTAL required non-whitespace ASCII identifiers,
and its CLUSTAL/Stockholm rejected zero-column alignments. Those guards never
changed the Rust writers or claimed support for new dialects.

## Recommended format by use

| Need | Format | Facade method | Conventional extension |
|---|---|---|---|
| Broad sequence-tool compatibility | aligned FASTA | `run.write_fasta(...)` | `.fasta`, `.fa` |
| Human review and alignment viewers | CLUSTAL | `run.write_clustal(...)` | `.aln` |
| Maximum compatibility with classic phylogeny tools | PHYLIP | `run.write_phylip(...)` | `.phy`, `.phylip` |
| Typed matrix plus room for future phylogenetic blocks | NEXUS | `run.write_nexus(...)` | `.nex`, `.nexus` |
| Alignment exchange in Stockholm-aware tools | Stockholm 1.0 | `run.write_stockholm(...)` | `.sto`, `.stockholm` |
| Stable application or future binding payload | JSON v1 | `run.write_json(...)` | `.json` |
| Authoritative POASTA graph topology | GFA 1.0 | `run.write_gfa(...)` | `.gfa` |
| POASTA graph visualization | Graphviz DOT | `run.write_dot(...)` | `.dot` |

The implementations follow the published [PHYLIP sequence-format
rules](https://phylipweb.github.io/phylip/doc/sequence.html), the [NEXUS
format](https://plewis.github.io/nexus/), and the [Stockholm 1.0
specification](https://sonnhammer.sbc.su.se/Stockholm.html).

## Facade example

```rust,no_run
use align::prelude::*;

let records = read_fasta::<Dna>("input.fasta")?;
let run = Multiple::mafft(MafftStrategy::Linsi).run(&records)?;

run.write_fasta("alignment.fasta")?;
run.write_clustal("alignment.aln")?;
run.write_phylip("alignment.phy")?;
run.write_nexus("alignment.nex")?;
run.write_stockholm("alignment.sto")?;
run.write_json("alignment.json")?;
run.write_text("alignment.txt")?;
# Ok::<(), align::Error>(())
```

All facade methods return `align::Result<()>`. The lower-level `align-formats`
crate also exposes path writers and `_to` stream variants.

## PHYLIP

The default is relaxed sequential PHYLIP. It writes the exact row count and MSA
width, retains stable row order, and requires portable ASCII identifiers without
whitespace or the PHYLIP-reserved punctuation `()[]:;,`.

`PhylipOptions` controls the two independent compatibility choices:

- `PhylipLayout::Sequential` or `PhylipLayout::Interleaved`;
- `PhylipNameStyle::Relaxed` or `PhylipNameStyle::Strict`.

Strict mode uses the original fixed ten-byte taxon field. Identifiers longer than
ten bytes are rejected; they are never truncated. Leading or trailing whitespace
is also rejected because padding would make the encoded identity ambiguous.
Interleaved mode uses a positive configurable block width, defaulting to 60.

```rust,no_run
# use align::prelude::*;
# let records = read_fasta::<Dna>("input.fasta")?;
# let run = Multiple::mafft(MafftStrategy::Linsi).run(&records)?;
let options = PhylipOptions::new()
    .with_layout(PhylipLayout::Interleaved)
    .with_name_style(PhylipNameStyle::Strict)
    .with_block_width(60)?;
run.write_phylip_with("alignment-strict.phy", options)?;
# Ok::<(), align::Error>(())
```

Relaxed PHYLIP is the safe default for meaningful identifiers. Strict PHYLIP is
appropriate only when every taxon already has a unique portable name of at most
ten bytes.

## NEXUS

The writer emits one self-contained `DATA` block:

```text
#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX=... NCHAR=...;
  FORMAT DATATYPE=DNA MISSING=? GAP=-;
  MATRIX
    'taxon' ACGT
  ;
END;
```

`DATATYPE` is selected from the Rust type: `DNA`, `RNA`, or `PROTEIN`. Every taxon
label is single-quoted; an embedded apostrophe is doubled. Spaces and punctuation
are therefore retained without lossy renaming. Control characters are rejected.
The current writer intentionally emits only the aligned character matrix. It does
not invent `TREES`, `ASSUMPTIONS`, or model blocks.

## Stockholm 1.0

The writer emits the required `# STOCKHOLM 1.0` header, one unwrapped aligned row
per input sequence, and the `//` terminator. The current typed MSA has no column or
sequence annotation model, so no synthetic `#=GF`, `#=GS`, `#=GR`, or `#=GC`
records are emitted. Sequence names containing whitespace are rejected because
Stockholm uses whitespace to separate the name from aligned symbols. Names that
start with `#` or equal `//` are also rejected because they would collide with
markup or the document terminator.

## JSON schema `align-output-v1`

Every JSON document contains:

- `schema: "align-output-v1"`;
- a `kind` discriminator;
- exact algorithm/profile and molecule metadata;
- native zero-based coordinates and typed objective values;
- the data required by that result kind, without parsing a human report.

| `kind` | Typed payload |
|---|---|
| `pairwise_alignment` | inputs, half-open ranges, objective, statistics, extended CIGAR, aligned rows, retained details |
| `pairwise_value` | inputs and score/distance/penalty only; `has_traceback` is false |
| `multiple_alignment` | profile, dimensions, statistics, ordered identifiers and aligned rows, retained evidence |
| `myers_distance` | pattern, text, and exact global unit edit distance |
| `myers_search` | pattern, text, threshold, and ordered exclusive text endpoints |
| `poasta_graph_query` | graph dimensions, query, costs, backend, search evidence, operations, and aligned pairs |

JSON coordinates preserve native Rust semantics: ranges are zero-based and
half-open, POASTA node/query indices are zero-based, and Myers `text_end` is an
exclusive endpoint. Pairwise identity is a ratio in `[0, 1]` or `null` for an
empty alignment. `score`, `distance`, and `penalty` remain distinct objectives,
and each objective declares `maximize` or `minimize`.

Output is deterministic UTF-8 with one final newline. Strings use complete JSON
escaping, including control characters and JavaScript line separators. Core
scientific fields are typed JSON values. Algorithm-specific `summary` and
`details` arrays are display evidence and must not be parsed to reconstruct native
results. A breaking schema change requires a different schema identifier.

The JSON writer is implemented directly in Rust without adding a production
serialization dependency. `serde_json` is used only by tests to parse and validate
every emitted result kind.

## Why there is no Newick writer yet

Newick represents a tree, not an alignment. The current portfolio computes
pairwise alignments, MSAs, and partial-order graphs; it does not expose a supported
phylogenetic inference result with branch lengths, rooting semantics, or an
evolutionary model. Internal MSA guide trees organize progressive alignment and
must not be presented as inferred phylogenies.

A future tree layer should first introduce a distinct validated phylogenetic tree
type and inference contract. Newick and a NEXUS `TREES` block can then serialize
that type without conflating alignment strategy with evolutionary inference.

## Compatibility rules

- Writers preserve row order, alignment width, symbols, and gaps exactly.
- No writer truncates, hashes, deduplicates, or silently rewrites identifiers.
- A format-specific identifier restriction returns `FormatError` before creating
  a path output file whenever validation is possible in advance.
- Human reports are not parsed to produce any interchange output.
- Changing a file extension never changes the selected algorithm or profile.
