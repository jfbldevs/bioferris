# Human-readable output contract v1

## Contents

- [Full API index](../README.md#api-index)
- [Complete example reports](../examples/alignment-recipes.md#global-local-and-affine-pairwise-alignment)
- [Layers](#layers)
- [Coordinates](#coordinates)
- [Objective names](#objective-names)
- [Pairwise statistics and symbols](#pairwise-statistics-and-symbols)
- [Multiple alignments](#multiple-alignments)
- [Detail and deterministic text](#detail-and-deterministic-text)
- [File formats](#file-formats)
- [Compatibility](#compatibility)

Status: implementation contract for the Rust-native `align` facade.

This contract governs presentation only. Algorithm results, objective functions,
tracebacks, profile identifiers, graph evidence, and validation rules remain
unchanged and authoritative.

## Layers

Output is deliberately split into four layers:

```text
algorithm -> native result -> run/report view -> text or interchange format
```

An algorithm never writes to standard output. A `Run` retains its native result
and the input references or metadata required to explain it. A report borrows that
run and renders incrementally; creating a report does not clone a sequence,
alignment, or graph.

## Coordinates

- Native Rust ranges remain zero-based and half-open.
- Normal human reports use one-based inclusive coordinates.
- An empty consumed range is written as `-` rather than inventing a coordinate.
- Pairwise block endpoints count biological symbols, never gap characters.
- MSA row endpoints advance independently and ignore gaps.
- POASTA query coordinates count query residues. Graph nodes retain their stable
  zero-based `PoaNodeId` values because a branched graph has no single biological
  coordinate axis.

## Objective names

Reports must not collapse different objectives into an ambiguous `Result` value:

- Needleman--Wunsch, Smith--Waterman, Gotoh, and X-drop report `Score`;
- Myers reports `Distance`;
- WFA reports `Penalty` and may additionally expose `Alignment score` in detailed
  output;
- POASTA reports its maximizing `Score` and identifies it as the negated affine
  path cost;
- MSA algorithms report their profile identifier and only display an objective
  when that native result actually retains one.

Score-only reports explicitly include `score-only` in their mode. Because no
traceback exists, they do not display ranges, identity, CIGAR, or aligned rows.
They retain the untouched native summary through `.result()`.

## Pairwise statistics and symbols

`Columns` is the number of traceback operations. `Identity` is exact matches
divided by all alignment columns, including substitutions and gaps. The report
separately counts:

- exact matches;
- substitutions;
- inserted residues and insertion runs;
- deleted residues and deletion runs.

Extended CIGAR uses `=`, `X`, `I`, and `D`. Its direction follows the shared
`PairwiseAlignment` convention: `I` consumes only the second sequence and `D`
consumes only the first.

The middle display line uses:

- `|` for an exact match;
- `.` for a substitution;
- a space for a gap column.

Ambiguity symbols are compared literally, exactly as the kernels compare them.

## Multiple alignments

MSAs render in stable row order and fixed-width column blocks. The conservation
line contains `*` only when every row has the same literal biological symbol and
none has a gap. All other columns contain a space. No biochemical similarity
classes are implied.

The summary reports row count, column count, exact conserved-column count, and
total gap characters. Identifiers are never truncated by default.

## Detail and deterministic text

Three verbosity levels are defined:

- `Compact`: identity of the algorithm/profile and its primary dimensions/value;
- `Normal`: summary, statistics, and block alignment;
- `Detailed`: normal output plus retained execution evidence specific to the
  selected algorithm.

The default block width is 60 alignment columns. Width zero is rejected before
rendering. `Display` output is deterministic UTF-8, uses `\n`, contains no ANSI
escape sequences, and has no trailing newline. This makes `println!("{report}")`
behave consistently. The provided text-file helper adds exactly one final newline.
Control characters in caller-supplied display labels are escaped visibly rather
than injected into the output. Terminal color is an optional presentation concern
and must never change redirected or golden output.

## File formats

Human reports are not machine interchange formats. The format crate separately
provides aligned pairwise FASTA, aligned MSA FASTA, CLUSTAL, PHYLIP, NEXUS,
Stockholm 1.0, Myers TSV, POASTA GFA, and POASTA Graphviz DOT writers. Those
writers consume native typed results and do not parse human report text.

Facade runs expose the same output workflow directly through `.write_text(...)`,
`.write_json(...)`, and result-appropriate methods such as `.write_fasta(...)`,
`.write_phylip(...)`, `.write_nexus(...)`, `.write_stockholm(...)`,
`.write_tsv(...)`, `.write_gfa(...)`, and `.write_dot(...)`.

JSON now uses the separately frozen `align-output-v1` schema. It preserves typed
objective semantics and native coordinates for pairwise, MSA, Myers, and POASTA
query results. Consumers must still not treat human text as a JSON substitute or
parse it back into scientific results. The complete interchange contract is in
[`phylogenetic-output-v1.md`](phylogenetic-output-v1.md).

## Compatibility

The pre-existing `.align()`, `.score()`, `.distance()`, `.penalty()`, and
`.find_all()` methods remain supported. New `.run()`, `.run_score()`,
`.run_penalty()`, `.run_distance()`, and `.run_find_all()` methods are additive and
call the same frozen kernels. A future Python facade may consume the structured
runs or native results, but must not reproduce report statistics or algorithm
logic.
