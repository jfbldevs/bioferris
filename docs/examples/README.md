# Worked API examples

## Example index

| Workflow | Input and result shown | Runnable guide |
|---|---|---|
| First complete Rust program | DNA `ATGGCTTAA`, complement, reverse complement, RNA and protein; run with `cargo bio my_analysis` | [Quick start](../quick-start.md#use-the-rust-library) |
| Checked regions and literal patterns | Typed slices, concatenation, overlap counts, RNA/protein search and preserved FASTQ qualities | [Sequence operation examples](../api/sequence-operations-v1.md#run-a-rust-analysis) |
| Sequence and streaming statistics | Literal DNA/RNA/protein, GC policies, bounded k-mers and complete file summaries | [Statistics examples](../api/sequence-statistics-v1.md) |
| Combine library operations | Coding DNA, translated protein, two DNA sequences and a full alignment report | [Translation and alignment](../quick-start.md#combine-translation-and-alignment) |
| DNA/RNA transformations | Literal sequences, complement/reverse complement, transcription and back-transcription | [Quick start](../quick-start.md#use-the-rust-library) and [supported operations](../api/seq-bio-v1.md#supported-operations) |
| Translation policies | Full, first-stop and complete-CDS translation, alternative starts and ambiguous codons | [Translation modes and ambiguity](../api/genetic-codes-v1.md#translation-modes-and-ambiguity) |
| Genetic codes and forward frames | NCBI tables 1/2/11 in the examples (every current table supported), frame metadata, table-specific starts/stops and typed errors | [Configured translation examples](../api/genetic-codes-v1.md) |
| Six-frame translation and ORFs | Both strands, original-input coordinates, alternative starts, partial candidates and rejected ambiguity | [Complete Rust examples](../api/reading-frames-v1.md) |
| Annotated CDS translation | Literal GenBank CDS, annotation-selected code/frame, source protein matching, mismatch/ignore behavior and result budgets | [Complete Rust examples](../api/cds-translation-v1.md) |
| Canonical protein properties | Literal INGAR/PETER/DD/AX, mass/charge/pI/aromaticity/GRAVY, boundary status and errors | [Complete Rust examples](../api/protein-properties-v1.md) |
| IUPAC and matrix motifs | Literal ambiguous DNA, both strands, count/PWM/PSSM, pseudocount/background, float32 scores and -inf; `cargo bio motif_analysis` | [Complete motif workflows](../api/motifs-v1.md) |
| Alignment analysis | Literal FASTA/Stockholm/NEXUS, strict slices, reverse/local origins, gap/identity counts, consensus and full annotation-loss reports; `cargo bio alignment_analysis` | [Complete analysis workflows](../api/alignment-analysis-v1.md) |
| BLAST result import | Tabular 6/7 and legacy/XML2 examples, reported reverse endpoints, absent values and bounded errors; `cargo bio blast_analysis`; streaming large XML reports, `cargo bio blast_stream` | [Complete BLAST reader workflows](../api/blast-input-v1.md) |
| Compressed records and BGZF indexes | Preserved case/qualities, gzip/BGZF output, virtual seek, verified sidecars and create-only errors; `cargo bio compressed_io` | [Complete compressed-I/O workflows](../api/compressed-io-v1.md) |
| Machine-readable output | Alignment JSON (`align-output-v1`), Myers distance JSON and the other export formats written by Rust run writers | [Export recipes](alignment-recipes.md#fasta-round-trip-and-alignment-export) |
| Errors and recovery | Partial codon matched as a typed error, an explicit ignore policy, and the stderr and exit status of an unhandled invalid DNA input | [Handle invalid input](../quick-start.md#handle-invalid-input) |
| Pairwise alignment | Global affine, global linear and local reports with scores, ranges, identity and CIGAR | [Rust recipes](alignment-recipes.md#global-local-and-affine-pairwise-alignment) |
| Combined alignment workflow | Distance, a multiple alignment and a graph query in one Rust program | [Combined Rust program](alignment-recipes.md#combine-distance-msa-and-a-graph-query) |
| Edit distance and search | Myers distance/search endpoints, WFA costs and adaptive affine results | [Rust alignment recipes](alignment-recipes.md) |
| Seeded extension | X-drop seed, direction and resulting alignment | [Rust alignment recipes](alignment-recipes.md) |
| Multiple alignment | Protein FASTA inputs and actual Clustal Omega, MUSCLE and MAFFT outputs | [Rust alignment recipes](alignment-recipes.md) |
| Graph alignment | DNA records, POASTA graph alignment and GFA output | [Rust alignment recipes](alignment-recipes.md) |
| File conversion and export | FASTA round-trip and alignment matrix/JSON output | [Rust alignment recipes](alignment-recipes.md) |
| Streaming records | Literal FASTQ, decoded qualities and FASTA conversion | [Sequence file recipes](../api/file-io-v1.md) |
| Persistent record retrieval | Literal FASTQ, Rust build/save/reopen, canonical records, missing keys and create-only refusal; binary sidecars with full and quick source checks, `cargo bio binary_index` | [Index recipes](sequence-index-recipes.md) |
| Alignment import | Literal aligned FASTA, CLUSTAL, PHYLIP, Stockholm and NEXUS with checked results | [Alignment reader recipes](../api/alignment-input-v1.md) |
| Annotated records | Literal GenBank/EMBL, metadata and strand-aware extraction; streaming a file with a record-error policy, `cargo bio stream_annotations` | [Annotation recipes](../api/annotated-input-v1.md) and [record streams](../api/annotated-input-v1.md#stream-records-from-large-files) |
| Annotated writing and transforms | Literal input/output, complete CDS slicing, reverse complement, loss reports and create-only files | [Annotated output recipes](annotated-output-recipes.md) |
| Graph/result import | Literal GFA, DOT and existing JSON with validated outputs | [Graph/result recipes](../api/graph-result-input-v1.md) |
| Real sequence workflow | Complete 276-base RefSeq FASTA input and its 91-residue protein output | [Real FASTA example](../quick-start.md#work-with-a-real-fasta-file) |

## How to use these examples

Start from the repository root; see the
[repository setup](../quick-start.md#enter-the-repository). Each worked
case identifies its literal input, complete executable command, expected result
and the meaning of that result. Output excerpts, when used, are explicitly
labeled; complete captured outputs are not shortened with ellipses.

The quick start introduces writing and running your own Rust program, using the
repository's [`cargo bio` alias](../quick-start.md#short-commands-for-your-programs).
The Rust alignment page demonstrates alignment workflows. The
CLI recipes and alignment agent recipes
remain only as historical records of the removed CLI; their commands no longer
run. All commands run locally. The examples use synthetic data unless a linked fixture provenance says
otherwise; none is a new scientific-fidelity or performance benchmark.
Motif, alignment-analysis, BLAST and compressed-I/O examples are implemented and
verified for their documented scope; see the integration validation.
Protein-property verification is limited to its documented models and reference
corpus. Passing recipes does not establish universal reference equivalence.

## Result conventions

- Scores are maximized; edit distances and WFA/POASTA costs are minimized. Compare
  objectives only when their scoring model and parameters agree.
- Human alignment reports use one-based inclusive coordinates. Native Rust/JSON
  ranges use zero-based, end-exclusive coordinates.
- A program run with `cargo bio` exits with status `0` on success. When `main`
  returns an error, Rust prints `Error: ` and the error's debug form to stderr
  and exits with status `1`; Cargo's own errors exit with `101`.
- FASTA input and output are shown separately. Legacy `SequenceRecord` retains the
  complete trimmed header as its identifier. New `Record` splits ID/description
  and retains header/case/quality; annotated imports use a separate rich document.
- Motif PSSM scores favor the model when positive; binary64 matrices round each
  window once to binary32. Zero-probability matrix scores are negative infinity
  (printed `-inf`), not missing values.
- Alignment missing data is not a gap; consensus evidence/status and identity
  denominator policies accompany results. BLAST endpoints remain in the reported
  coordinate units; absent scores/strand are not synthesized.
- Compressed output preserves decoded bytes, not identical native compression.
  BGZF virtual positions are not linear offsets, and verified index reopen reads
  the complete source before later retained-handle lookups.

The [API index](../README.md#api-index) and feature matrix
describe the full supported scope.

## Keep examples verified

Rust integration tests compare each guide's documented program with its source
file and its documented output with a real run. The quick-start test also invokes
the checked-in Cargo alias offline with the lockfile and compares exit status,
both streams and the full output of `my_analysis`, `quick_start` and both
`invalid_input` runs, plus the documented Cargo error for a nonexistent example.
It needs the pinned Cargo toolchain and cached dependencies, but no prebuilt
example binary. Rust alignment recipes have their own computed-output checks.
Run:

```bash
cargo test -p bioferris --test documentation_examples --test alignment_workflow_guide --offline --locked
cargo test -p bioferris --example alignment_recipes --offline --locked
cargo test -p bioferris --test sequence_operations --test statistics --test genetic_codes --test reading_frame_recipes --offline --locked
cargo test -p bioferris --test file_io_guide --test alignment_import --test annotated_import --test annotated_streaming --test graph_result_import --offline --locked
cargo test -p bioferris --test sequence_index --test binary_sequence_index --test sequence_index_recipes --offline --locked
cargo test -p bioferris --test annotated_output --test annotated_output_recipes --offline --locked
cargo test -p bioferris --test cds_translation --test cds_recipes --offline --locked
cargo test -p bioferris --test protein_properties --test protein_properties_recipes --offline --locked
cargo test -p bioferris --test motifs --test motif_recipes --offline --locked
cargo test -p bioferris --test alignment_analysis_public --test alignment_analysis_recipes --offline --locked
cargo test -p bioferris --test blast_recipes --test compression_recipes --offline --locked
cargo test -p bioferris --lib --offline --locked
```

The checks fail when the documented input/output contract drifts. They do not
update expected results automatically; a changed result requires investigation
and an explained documentation update under the
documentation policy.
