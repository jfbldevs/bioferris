# BioFerris documentation

<!-- Public documentation index of the exported library; scripts/export_public_library.py copies it to docs/README.md. -->

Rust-native bioinformatics, from sequences to workflows.

## Contents

- [API index](#api-index)
- [Guides](#guides)
- [Scope and limits](#scope-and-limits)
- [Algorithm specifications](#algorithm-specifications)
- [Testing](#testing)
- [Type and method reference](#type-and-method-reference)

[Project overview](../README.md) · [Quick start](quick-start.md) · [Worked examples](examples/README.md)

## API index

| Capability | Rust entry points | Reference |
|---|---|---|
| Typed sequences | `dna`, `rna`, `protein` | [Sequence construction](api/seq-bio-v1.md#short-path-for-manual-rust-use) |
| Checked sequence regions and concatenation | `SequenceExt::slice`, `concat`; `Record::slice` preserves case/quality | [Sequence operations](api/sequence-operations-v1.md) |
| Literal search and counting | `find`, `contains`, `count`, `count_overlapping`, lazy `find_all` | [Search semantics](api/sequence-operations-v1.md#literal-search-and-overlap) |
| Composition, GC and bounded k-mers | `statistics`, `statistics_with`, `summarize_sequences`, streaming `summarize_records` | [Statistics contract and examples](api/sequence-statistics-v1.md) |
| Protein properties | `protein_properties`, `protein_properties_with`, borrowed `analyze_protein`; average mass, charge/pI, aromaticity and GRAVY | [Models, precision and complete results](api/protein-properties-v1.md) |
| IUPAC nucleotide motifs | `dna_motif`, `rna_motif`, `IupacMotif::find_with`; explicit ambiguity and strand policies | [Motif workflows](api/motifs-v1.md) |
| Count, probability and log-odds matrices | `MotifCounts` → `MotifWeights` → `MotifScores`; explicit pseudocount/background, float32 scan scores | [Models and complete results](api/motifs-v1.md#build-and-score-matrices) |
| Complement | `complement()` | [Sequence operations](api/seq-bio-v1.md#supported-operations) |
| Reverse complement | `reverse_complement()` | [Sequence operations](api/seq-bio-v1.md#supported-operations) |
| Transcription | `transcribe()`, `back_transcribe()` | [Sequence operations](api/seq-bio-v1.md#supported-operations) |
| Standard-code translation | `translate()`, `translate_to_stop()`, `translate_cds()`, `translate_standard(options)` | [Translation policies](api/seq-bio-v1.md#translation-policies) |
| Configured genetic code and forward frame | `translate_with(TranslationConfig)`; every current NCBI table, frames 0/1/2 | [Genetic codes and frames](api/genetic-codes-v1.md) |
| Six-frame translation | `translate_six_frames`, `translate_six_frames_with`; both strands with original-input coordinates | [Reading-frame workflows](api/reading-frames-v1.md) |
| ORF discovery | `find_orfs`, `find_orfs_with(OrfOptions)`; explicit starts/stops/strands and count budget | [ORF policies and examples](api/reading-frames-v1.md#policies-and-coordinates) |
| Annotation-aware CDS translation | `Feature::translate_cds`, `translate_cds_with`; qualifiers, extraction and exact source-protein checks | [Complete CDS workflows](api/cds-translation-v1.md) |
| Validated protein conversion | `Translation::to_protein()` | [Translation result contract](api/seq-bio-v1.md#supported-operations) |
| Linear pairwise alignment | `Pairwise::needleman_wunsch()`, `Pairwise::smith_waterman()` | [Pairwise recipes](examples/alignment-recipes.md#global-local-and-affine-pairwise-alignment) |
| Affine pairwise alignment | `Pairwise::gotoh()`, `Pairwise::adaptive_affine()` | [Modes and defaults](examples/alignment-recipes.md#global-local-and-affine-pairwise-alignment) |
| Edit distance and approximate search | `Pairwise::myers()` | [Distance and endpoints](examples/alignment-recipes.md#edit-distance-approximate-search-and-affine-costs) |
| Wavefront alignment | `Pairwise::wfa()` | [Objective conventions](examples/alignment-recipes.md#edit-distance-approximate-search-and-affine-costs) |
| Seeded extension | `Pairwise::xdrop()` | [Heuristic extension](examples/alignment-recipes.md#seeded-x-drop-extension) |
| Clustal Omega MSA | `Multiple::clustal_omega()` | [MSA profiles](examples/alignment-recipes.md#clustal-omega-muscle-5-and-mafft) and [which profile to choose](api/rust-api-v1.md#choose-a-multiple-alignment-profile) |
| MUSCLE 5 MSA | `Multiple::muscle5()`, `Multiple::muscle5_super5()` | [MSA profiles](examples/alignment-recipes.md#clustal-omega-muscle-5-and-mafft) and [which profile to choose](api/rust-api-v1.md#choose-a-multiple-alignment-profile) |
| MAFFT MSA | `Multiple::mafft(strategy)`, `Multiple::mafft_7526_linsi()`, `Multiple::mafft_7526_einsi()` | [Profiles](examples/alignment-recipes.md#clustal-omega-muscle-5-and-mafft), [which profile to choose](api/rust-api-v1.md#choose-a-multiple-alignment-profile) and [thread limits](api/rust-api-v1.md#threads-and-unsupported-options) |
| Graph alignment and construction | `Graph::poasta()` | [Graph workflows](examples/alignment-recipes.md#poasta-graph-construction-and-gfa) |
| Alignment analysis | `AlignmentAnalysis`: strict slices, source-coordinate maps, gap/identity evidence and policy-based consensus | [Gaps, missing data and annotation-loss policies](api/alignment-analysis-v1.md) |
| Legacy records and FASTA | Unchanged `SequenceRecord`, `read_fasta`, `parse_fasta`, `write_fasta` | [Legacy FASTA contract](api/rust-api-v1.md#fasta) |
| Rich records, FASTA and Sanger FASTQ | `Record`, `io::fasta`, `io::fastq`, streaming readers and checked writers | [File I/O and examples](api/file-io-v1.md) |
| Persistent record indexes | `io::index::SequenceIndex`: build, JSON or binary sidecars (`save_new`, `save_new_as`), full verified open, opt-in quick check (`open_with`, `SourceCheck::Quick`), get, source-order IDs | [Index contract](api/sequence-index-v1.md), [binary layout](algorithms/binary-sequence-index-v1.md) and [worked examples](examples/sequence-index-recipes.md) |
| Compressed byte output and BGZF streams | `io::gzip`, `io::bgzf`: explicit writer finish, checked block reads and virtual seeks | [Compressed-I/O workflows](api/compressed-io-v1.md) |
| Verified BGZF record indexes | `io::bgzf_index::BgzfSequenceIndex`: build/save/full verified open and checked retrieval | [Bounds and preservation](api/compressed-io-v1.md#limits-and-compatibility) |
| BLAST result import | `io::blast`: streaming/collecting tabular 6/7 and bounded legacy/XML2 profiles; `BlastXmlReader` for large reports | [Reported values and complete examples](api/blast-input-v1.md) and [streaming large reports](api/blast-input-v1.md#stream-large-xml-reports) |
| Alignment import | `io::alignment`: aligned FASTA, CLUSTAL, explicit PHYLIP dialects, Stockholm, NEXUS DATA subset | [Alignment input](api/alignment-input-v1.md) |
| Annotated records and feature extraction | `io::annotated`: GenBank/EMBL subsets, exact ranges/join/complement; `genbank`/`embl` and `AnnotatedReader` stream one record at a time | [Annotated input](api/annotated-input-v1.md) |
| Checked annotated writing and transformations | `write_genbank`, `write_embl`, create-only path helpers; annotated `slice`/`reverse_complement` and typed change reports | [Annotated output](api/annotated-output-v1.md) and [complete recipes](examples/annotated-output-recipes.md) |
| Graph and result imports | `io::graph`, `io::result_json`: GFA, DOT topology, result-JSON subsets | [Graph/result input](api/graph-result-input-v1.md) |
| Human reports and native results | `run.report()`, `run.result()` | [Output contract](api/output-contract-v1.md) |
| Alignment and graph export | Run writers for FASTA, CLUSTAL, PHYLIP, NEXUS, Stockholm, JSON, Myers TSV, GFA and DOT, depending on result type | [Export recipes](examples/alignment-recipes.md#fasta-round-trip-and-alignment-export) and [formats](api/phylogenetic-output-v1.md) |
| Typed errors | Typed Rust error values; stable `code()` strings where the error type provides them | [Invalid input](quick-start.md#handle-invalid-input) and each guide's errors section |

Use `bioferris::prelude::*` for the common sequence and alignment types; readers and
writers are explicit imports from `bioferris::io` and its documented submodules. Each
row describes an implemented, tested subset, documented in its guide; none establishes
universal equivalence with Biopython or with the native programs.

## Guides

| Topic | Document | Contents |
|---|---|---|
| First use | [Quick start](quick-start.md) | Write/run a Rust analysis, `cargo bio` shortcut and real FASTA input/output |
| Worked examples | [Input/output recipe index](examples/README.md) | Exact inputs, runnable commands, complete results, error streams and interpretation |
| Sequences | [`bioferris` biological API](api/seq-bio-v1.md) | Complements, transcription, translation, ambiguity and limits |
| Genetic codes and frames | [Configured translation](api/genetic-codes-v1.md) | Every current NCBI table, explicit forward offsets, CDS checks, complete Rust examples |
| Six frames and ORFs | [Reading-frame workflows](api/reading-frames-v1.md) | Both strands, original-input ranges, stop coordinates, policies, budgets and verified examples |
| Annotated CDS | [Annotation-aware translation](api/cds-translation-v1.md) | GenBank/EMBL qualifiers, complete-CDS checks, exact reference comparison and full Rust results |
| Regions and literal search | [Checked sequence operations](api/sequence-operations-v1.md) | Strict ranges, concatenation, overlap and preserved quality slices |
| Alignment workflows | [Alignment recipes](examples/alignment-recipes.md) | Pairwise, search, explicit MSA profiles, POASTA, FASTA input and create-only exports in Rust |
| Sequence files | [File I/O v1](api/file-io-v1.md) | Streaming records, FASTQ qualities, gzip and safe conversion |
| Record indexing | [Persistent indexes](api/sequence-index-v1.md) | Stable plain FASTA/FASTQ, strict JSON or binary sidecars, full or opt-in quick source checks, preserved retrieval and errors |
| Sequence statistics | [Statistics v1](api/sequence-statistics-v1.md) | Count-first composition/GC, ambiguity, bounded k-mers and streaming file summaries |
| Protein properties | [Property workflows](api/protein-properties-v1.md) | Canonical peptides, explicit models/pH, pI boundaries, exact reference scope and complete Rust results |
| Nucleotide motifs | [IUPAC and matrix workflows](api/motifs-v1.md) | Pattern ambiguity, count/PWM/PSSM, pseudocounts/background, reverse scoring, precision, budgets and complete results |
| Alignment analysis | [Analysis workflows](api/alignment-analysis-v1.md) | Strict slices, original-source maps, exact identity/gap evidence, consensus policies and complete annotation-loss reports |
| BLAST results | [Tabular and XML workflows](api/blast-input-v1.md) | Reported optional values/endpoints, outfmt 6/7, legacy XML and namespace-aware XML2 subsets, streaming and resource bounds |
| Compressed files | [Gzip/BGZF and indexed retrieval](api/compressed-io-v1.md) | Explicit finalization, byte preservation, virtual offsets, verified sequence sidecars and create-only output |
| Alignment input | [Alignment readers](api/alignment-input-v1.md) | Supported dialects, annotations, missing cells and checked adapters |
| Annotation input | [GenBank/EMBL readers](api/annotated-input-v1.md) | Metadata, qualifiers, exact extraction, preserved unsupported locations and record streams with an explicit error policy |
| Annotation output | [Checked writers and transforms](api/annotated-output-v1.md) | Same-profile output, new local identities, strict coordinates and explicit annotation loss |
| Graph/result input | [Graph and JSON imports](api/graph-result-input-v1.md) | Project profiles, topology limits and no provenance claims |
| Alignments | [`align` Rust API](api/rust-api-v1.md) | Pairwise, multiple and graph alignment, and which multiple-alignment profile to choose |
| Reports | [Output contract](api/output-contract-v1.md) | Human-readable reports, coordinates and native results |
| Alignment formats | [Export and structured JSON](api/phylogenetic-output-v1.md) | FASTA, CLUSTAL, PHYLIP, NEXUS, Stockholm and graphs |

## Scope and limits

- Typed sequences are normalized to uppercase and reject gaps. Configured translation
  supports every current NCBI table and forward offsets 0, 1 and 2, with explicit stop,
  CDS and partial-codon policies; the convenience methods use table 1, frame 0.
- Streaming FASTA and Sanger FASTQ records keep their header, source case and
  qualities; gzip input is detected. Alignment, annotated, graph and JSON readers
  support the explicit subsets their guides document.
- Named file output is created new, and writers that would drop information, such as
  qualities, require an explicit policy.
- Opening an index verifies its source in full by default; the opt-in quick check of a
  binary sidecar is faster but cannot detect every change. The BGZF index is separate.
- BLAST imports do not run searches or invent absent values; XML readers bound depth,
  event size and query counts and never resolve remote entities.
- Motif matrices fix their numerical precision; they do not infer binding significance.
- The multiple-alignment reconstructions and the MAFFT 7.526 compatibility profiles
  follow different contracts; see [which profile to choose](api/rust-api-v1.md#choose-a-multiple-alignment-profile).
- Broader formats and other Biopython domains are not covered yet.

## Algorithm specifications

Each specification states what the method computes, its sources, parameters,
complexity and limits.

| Method or contract | Specification |
|---|---|
| Needleman–Wunsch | [needleman-wunsch.md](algorithms/needleman-wunsch.md) |
| Smith–Waterman | [smith-waterman.md](algorithms/smith-waterman.md) |
| Gotoh | [gotoh.md](algorithms/gotoh.md) |
| Affine local and end-gap modes | [affine-alignment-modes.md](algorithms/affine-alignment-modes.md) |
| Affine score-only execution | [affine-score-only.md](algorithms/affine-score-only.md) |
| Affine SIMD backend | [affine-simd.md](algorithms/affine-simd.md) |
| Linear score-only execution | [linear-score-only.md](algorithms/linear-score-only.md) |
| Myers bit-vector edit distance | [myers-bit-vector.md](algorithms/myers-bit-vector.md) |
| Wavefront alignment | [wfa.md](algorithms/wfa.md) |
| X-drop extension | [x-drop.md](algorithms/x-drop.md) |
| Clustal Omega | [clustal-omega.md](algorithms/clustal-omega.md) |
| MUSCLE 5 | [muscle5.md](algorithms/muscle5.md) |
| MAFFT FFT-NS-1 | [mafft-fft-ns-1.md](algorithms/mafft-fft-ns-1.md) |
| MAFFT FFT-NS-2 | [mafft-fft-ns-2.md](algorithms/mafft-fft-ns-2.md) |
| MAFFT FFT-NS-i | [mafft-fft-ns-i.md](algorithms/mafft-fft-ns-i.md) |
| MAFFT G-INS-i | [mafft-g-ins-i.md](algorithms/mafft-g-ins-i.md) |
| MAFFT L-INS-i | [mafft-l-ins-i.md](algorithms/mafft-l-ins-i.md) |
| MAFFT E-INS-i | [mafft-e-ins-i.md](algorithms/mafft-e-ins-i.md) |
| MAFFT 7.526 compatibility | [mafft-7526-compatibility.md](algorithms/mafft-7526-compatibility.md) |
| POASTA partial-order alignment | [poasta.md](algorithms/poasta.md) |
| Alignment analysis and coordinates | [alignment-analysis-v1.md](algorithms/alignment-analysis-v1.md) |
| Six-frame translation and ORFs | [reading-frames-v1.md](algorithms/reading-frames-v1.md) |
| Annotation-aware CDS translation | [cds-translation-v1.md](algorithms/cds-translation-v1.md) |
| Protein properties | [protein-properties-v1.md](algorithms/protein-properties-v1.md) |
| Motifs and scoring matrices | [motifs-v1.md](algorithms/motifs-v1.md) |
| BLAST result input | [blast-input-v1.md](algorithms/blast-input-v1.md) |
| Binary sequence index sidecar | [binary-sequence-index-v1.md](algorithms/binary-sequence-index-v1.md) |
| BGZF streams and virtual offsets | [bgzf-v1.md](algorithms/bgzf-v1.md) |
| Verified BGZF indexes | [bgzf-index-v1.md](algorithms/bgzf-index-v1.md) |
| Deterministic compressed output | [compressed-output-v1.md](algorithms/compressed-output-v1.md) |

## Testing

The [testing guide](testing/README.md) describes how the library is tested: unit,
property and integration tests, reference fixtures from Biopython and other tools,
documented examples checked against real runs, and security and robustness checks,
with what each establishes and what it does not.

## Type and method reference

Generate the HTML reference from the repository root with
`cargo doc -p bioferris --no-deps --open`. The entry point is
`target/doc/bioferris/index.html`.
