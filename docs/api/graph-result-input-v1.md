# Graph and result imports, version 1

## Contents

- [Capability index](#capability-index)
- [Worked Rust program](#worked-rust-program)
- [Graph contracts](#graph-contracts)
- [Result JSON contract](#result-json-contract)
- [Errors and resource limits](#errors-and-resource-limits)
- [Compatibility and verification](#compatibility-and-verification)

[Documentation index](../README.md) · [Rust API](rust-api-v1.md) ·
Removed CLI contract

## Capability index

| Input | Rust entry points | Output |
|---|---|---|
| Project GFA 1.0 residue-graph subset | `io::graph::{parse_gfa, parse_gfa_reader, read_gfa}` | Immutable `GfaDocument<A>` with nodes, classes, edges and paths |
| Project DOT export dialect | `io::graph::{parse_dot, parse_dot_reader, read_dot}` | Immutable `DotTopology<A>`; **no source paths** |
| Native `align-output-v1` result | `io::result_json::{parse_result_json, parse_result_json_reader, read_result_json}` | `ImportedResult::Alignment` |
| Biological response written by the removed CLI, `schema_version: 1` | Same result JSON entry points | `ImportedResult::Cli`, including structured errors |

Importing never executes anything. The `bioferris` CLI was removed
(decision 029); its archived responses stay importable.

The `parse_*` graph functions take a UTF-8 `&str` and an explicit alphabet such
as `Dna`, `Rna` or `Protein`. The `*_reader` functions take `impl BufRead`.
The `read_*` functions take a filesystem path and detect gzip by signature, not
extension; concatenated gzip members are decoded. These are bounded, eager
document readers, not streaming record iterators. The JSON reader obtains the
alphabet from the validated document discriminator.

## Worked Rust program

This complete program imports a two-residue GFA path, the corresponding DOT
topology and an independently specified JSON distance result. The strings are
the **literal input**, including GFA tab delimiters and DOT escaped label lines.
It also demonstrates two expected failures without discarding information.

From the repository root, save the following manifest as
`analysis/graph-import/Cargo.toml`, creating those directories if necessary:

```toml
[package]
name = "graph-import-analysis"
version = "0.1.0"
edition = "2024"

[workspace]

[dependencies]
bioferris = { path = "../../seq-bio" }
```

The empty `[workspace]` keeps this user analysis package separate from the
toolkit's workspace members. Save the program below as
`analysis/graph-import/src/main.rs`; no external graph program, interpreter or
alignment executable is required.

<!-- graph-result-program -->
```rust
use bioferris::{Dna, io::graph::{parse_dot, parse_gfa}};
use bioferris::io::result_json::{ImportedResult, parse_result_json};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gfa = "H\tVN:Z:1.0\nS\tn0\tA\tCL:i:0\nS\tn1\tC\tCL:i:1\n\
               L\tn0\t+\tn1\t+\t0M\tRC:i:1\nP\tread1\tn0+,n1+\t*\n";
    let dot = "digraph POASTA {\n  rankdir=LR;\n\
               n0 [label=\"n0\\nA\\nclass 0\"];\n\
               n1 [label=\"n1\\nC\\nclass 1\"];\n\
               n0 -> n1 [label=\"1\"];\n}\n";
    let json = r#"{"schema":"align-output-v1","kind":"myers_distance",
        "algorithm":"Myers bit-vector","mode":"global unit edit distance",
        "alphabet":"DNA","pattern":{"length":4,"sequence":"ACGT"},
        "text":{"length":3,"sequence":"AGT"},"distance":1}"#;

    let graph = parse_gfa::<Dna>(gfa)?;
    println!("GFA: {} nodes, {} edge, {} classes",
        graph.nodes().len(), graph.edges().len(), graph.class_count());
    println!("path {}: {}", graph.paths()[0].name(),
        graph.path_sequence(0).ok_or("missing path")?);

    let topology = parse_dot::<Dna>(dot)?;
    println!("DOT: {} nodes, {} edge; paths unavailable",
        topology.nodes().len(), topology.edges().len());

    let ImportedResult::Alignment(result) = parse_result_json(json)? else {
        return Err("expected a native result".into());
    };
    println!("JSON: {:?}, distance={}", result.kind(), result.value()["distance"]);

    let bad_gfa = parse_gfa::<Dna>(&gfa.replace("0M", "1M")).unwrap_err();
    println!("bad GFA: {}", bad_gfa.code());
    let bad_json = parse_result_json(&json.replace("\"distance\":1", "\"distance\":0"))
        .unwrap_err();
    println!("bad JSON: {}", bad_json.code());
    Ok(())
}
```

Run it from the repository root:

```bash
cargo run --quiet --manifest-path analysis/graph-import/Cargo.toml
```

Observed stdout from the program's checked workflow; successful process exit `0`:

<!-- graph-result-output -->
```text
GFA: 2 nodes, 1 edge, 2 classes
path read1: AC
DOT: 2 nodes, 1 edge; paths unavailable
JSON: MyersDistance, distance=1
bad GFA: unsupported_format_feature
bad JSON: invalid_schema
```

The path is `A` followed by `C`. The DOT file describes the same two nodes and
one edge, but its contents cannot reveal the path name or which sequence used
that edge. The JSON number `1` passes dimensional/equality checks; the importer
does **not** establish that it is the optimal edit distance. The bad GFA requests
an unsupported overlap, and the bad JSON falsely claims zero distance between
different sequences. Library errors are returned as values; this example handles
them, so they do not become process failures.

To verify the program directly from this guide in the repository root:

```bash
cargo build -p bioferris --lib --offline --locked
rustdoc --test docs/api/graph-result-input-v1.md --edition 2024 --extern bioferris=target/debug/libbioferris.rlib -L dependency=target/debug/deps -D warnings
cargo test -p bioferris --test graph_result_import --offline --locked
```

## Graph contracts

### GFA subset

The supported profile is narrower than general GFA: exactly one
`H\tVN:Z:1.0` header, then `S` nodes, then optional `L` links, then optional `P`
paths. Blank lines and `#` comments are ignored. CRLF and LF work. Required
fields are tab-separated; arbitrary extra tags are rejected, not dropped.

- Segment names are contiguous `n0`, `n1`, … in that order. Each segment has one
  valid residue and exactly one `CL:i:<class>` tag. Residues are normalized to
  the selected alphabet's uppercase representation; original case is not kept.
- Class identifiers are contiguous from zero. Distinct nodes in a class must
  have distinct residues. Links within a class and cycles in the class quotient
  are invalid.
- Links are sorted and unique by `(source, target)`, use `+`/`+` orientation,
  `0M` overlap and exactly one positive `RC:i:<support>` tag. All targets exist.
- Paths retain their file order and exact file-level names. Steps are forward
  `n<integer>+` references, with `*` in the overlap field. Every consecutive pair
  must have a link, a path cannot revisit a class, every node must be covered,
  and edge support must equal the number of paths traversing that link.
- Names are unique across paths and segments, follow GFA ASCII naming rules,
  and are not percent-decoded. The existing writer encodes a source identifier
  `reference one` as `reference%20one`; the imported name stays encoded.
- As a project-specific extension, `P\tempty\t*\t*` represents an empty path.
  A header-only empty graph is also accepted. These do not imply acceptance by
  unrelated GFA implementations.

Multi-residue/unknown-sequence segments, reverse traversal, nonzero or unknown
link overlap, additional optional tags, GFA 1.1/1.2 records and arbitrary segment
names are unsupported. Integers must be canonical unsigned decimal (no leading
zeroes except `0`) and fit `usize`. File-level name collisions are rejected even
if an existing writer can emit them: for example, a source identifier `n0`
conflicts with the first segment name. This is not a universal writer-round-trip
claim for every possible source identifier.
The removed alignment CLI export adapter
rejected these collisions before publication; the Rust writer remains unchanged.

`GfaDocument::nodes/edges/paths` expose immutable slices; node identifiers are
zero-based, edge support counts paths, and `path_sequence(index)` returns the
normalized source sequence or `None` for an absent index. No conversion to
`PoaGraph` is provided: the computational graph's public API cannot restore all
these identifiers, classes and occurrences directly without changing identity.

### DOT subset

Only the line-oriented dialect emitted by `write_poa_dot_to` is accepted:
`digraph POASTA {`, `rankdir=LR;`, numbered nodes with the three-part label shown
above, numbered edges with a positive integer support label, and a closing `}`.
Leading/trailing line whitespace is accepted. Nodes precede edges. The same
node/class/order/acyclicity checks apply as for GFA, but support cannot be checked
against paths because the file has none.

General Graphviz syntax, subgraphs, extra attributes, arbitrary identifiers,
HTML labels and comments are unsupported. DOT imports cannot recover source
paths, original sequence identifiers or residue occurrences; no fake paths are
created. This reader is a topology inspector, not a general DOT parser or a
computational graph loader.

## Result JSON contract

Exactly one JSON object is required. The parser rejects trailing documents,
duplicate keys (including escaped aliases), unexpected fields in supported
objects, missing fields, mismatched JSON types and unknown schema versions.
`schema: "align-output-v1"` and `schema_version: 1` are distinct discriminators;
having both is an error. Uppercase canonical DNA/RNA/protein symbols are required,
and only aligned rows accept `-` gaps. Metadata arrays are preserved as declared
strings; they are not converted into fresh execution evidence.

| Native kind | Checks performed in addition to the schema |
|---|---|
| `pairwise_alignment` | Input lengths/alphabet, zero-based half-open ranges, rectangular rows, ungapped rows versus input ranges, no double-gap columns, exact extended CIGAR, gap/match statistics, identity (absolute tolerance `1e-12`) and score-field agreement |
| `pairwise_value` | Inputs, objective direction/sign and `has_traceback: false`; no fabricated rows |
| `multiple_alignment` | Nonempty row collection, unique identifiers, rectangular rows, no all-gap column, row/column/conserved/gap counts |
| `myers_distance` | Input lengths, declared algorithm/mode, distance bounds and zero-distance/equality agreement; not a recomputation of optimal distance |
| `myers_search` | Declared algorithm/mode, one-based ranks, ordered unique exclusive text endpoints, endpoint/distance bounds; not start coordinates or a completeness proof |
| `poasta_graph_query` | Existing typed affine-cost constraints, graph dimensions, query consumption/order, unique bounded node IDs, operation/trace/aligned-pair agreement and selected search-counter bounds |

For responses of the removed CLI, the importer accepts its six original
operations: `translate`, `complement`, `reverse-complement`, `transcribe`,
`back-transcribe` and `pairwise`, plus well-formed `ok: false` errors. Transforms
and table-1 translations are recomputed with the typed library to check output
content and lengths. Pairwise responses must agree with their nested native
document, supported algorithm and the CLI's fixed scoring/mode contract. Help,
format-discovery and file-operation responses are **not** supported.

This profile accepts translation only at code 1, frame 0; the CLI's configured
translations with code 2/11 or frame 1/2 are rejected as
`unsupported_format_feature`. Responses of the CLI's later
alignment operations, motif, alignment-analysis, BLAST
and compressed-I/O operations are likewise outside it. Alignment `native-json`
exports remain importable as native results; reading them validates structure,
not historical execution claims.

`ImportedResult` distinguishes native results from CLI responses; each variant
provides a typed `kind()` and an immutable borrowed `serde_json::Value` through
`value()`. Imported errors remain errors in the source envelope, not successful
biological computations. Unknown operation/kind errors are explicit. These
functions never read execution requests, run arbitrary commands, modify files,
contact a network or infer authorization from result metadata.

**Validation is not scientific certification.** Native objective numbers and
algorithm/profile labels are assertions in the file. An admissible but suboptimal
Myers distance can pass; graph JSON does not contain topology, so edge continuity
and graph-dependent costs cannot be verified. Search completeness, optimality,
performance counters and historical execution provenance cannot be established
from the result alone. No imported document becomes a verified computational
run, and no old qualification evidence is rewritten.

## Errors and resource limits

All readers return `BioIoError`; branch on `kind()` or stable `code()`, not message
text. Graph parse errors attach one-based document number `1` and physical line
where available; global invariants have no precise line. JSON syntax errors have
their parser-reported line, while schema/content errors generally do not.

| Code | Meaning |
|---|---|
| `io_error` | Opening, reading or gzip validation failed |
| `invalid_header` | Required GFA header absent |
| `invalid_sequence` | A graph residue is invalid for the selected alphabet |
| `invalid_structure` | Graph names, grammar, classes, edges, paths or support are inconsistent |
| `truncated_input` | DOT closing brace absent |
| `unsupported_format_feature` | A graph feature, alphabet, result kind or CLI response operation is outside this contract |
| `invalid_schema` | JSON grammar, version, shape or supported content is inconsistent |
| `limit_exceeded` | Decompressed graph/result text exceeds 64 MiB |

Each graph/result document has a fixed **64 MiB decompressed text** ceiling,
checked before parsing. Reader functions consume at most the ceiling plus one
byte for oversize detection. This bounds source size, **not total resident
memory**: parsed strings, JSON values, graph collections and validation indexes
require additional allocations. Input is eagerly buffered. There is no custom
limit override, record indexing, output writer for imported documents or
performance guarantee in this increment. The JSON parser's recursion limit also
rejects excessively nested input. For files that pass the size limit, readers
consume through EOF, so gzip trailer/checksum failures are surfaced.

## Compatibility and verification

Implemented Rust-only subset, not complete GFA/Graphviz support or Biopython
parity. Syntax references: [GFA 1 specification](https://gfa-spec.github.io/GFA-spec/GFA1.html),
[Graphviz DOT language](https://graphviz.org/doc/info/lang.html) and
[RFC 8259](https://datatracker.ietf.org/doc/html/rfc8259), consulted 2026-10-15.
The project-specific tags, writer dialect and native JSON fields are governed by
the existing public Rust writers and [Rust API contract](rust-api-v1.md).

`seq-bio/tests/graph_result_import.rs` checks independent fixtures, existing graph
writer readback, all six native JSON writer kinds, seven pairwise modes/engines,
archived CLI responses/errors and exact transform envelopes, gzip path input,
malformed inputs, size/stream failures and the published output above. The program also runs as a Rust documentation
test. These are import compatibility/regression checks, not a reopened scientific
campaign or a performance measurement.
