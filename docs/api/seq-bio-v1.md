# Biological sequence workflows, P0

## Contents

- [Composition, GC and bounded k-mer statistics](sequence-statistics-v1.md)
- [Canonical protein properties and explicit numerical models](protein-properties-v1.md)
- [Additional analysis and interoperability](#additional-analysis-and-interoperability)
- [Full API index](../README.md#api-index)
- [Worked input/output examples](../examples/README.md)
- [Short path for manual Rust use](#short-path-for-manual-rust-use)
- [Explicit options and existing alignment API](#explicit-options-and-existing-alignment-api)
- [Supported operations](#supported-operations)
- [Translation policies](#translation-policies)
- [Offline example with a real sequence](#offline-example-with-a-real-sequence)
- [Deliberate limits and differences](#deliberate-limits-and-differences)
- [Public references](#public-references)

`bioferris` is an additive Rust-native facade on the existing alignment workspace.
It is a first functional subset inspired by Biopython **1.88**, not a complete
replacement, Python binding, or claim of identical Python API behavior.

The additive [checked sequence operations](sequence-operations-v1.md) provide
strict slicing, typed concatenation, literal find/count and preserved record
slices. They do not change the original transforms and translation below.

For literal inputs and complete results, start with the
[complete Rust workflow](../quick-start.md#use-the-rust-library).

## Short path for manual Rust use

```rust
use bioferris::prelude::*;

let product = dna("ATGGCTTAA")?.translate_cds()?;
assert_eq!(product.as_str(), "MA");
let protein = product.to_protein()?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`dna`, `rna`, and `protein` are validated constructors, not new sequence types.
`translate()`, `translate_to_stop()`, and `translate_cds()` are short forms of
the explicit standard-code options below. Partial codons still fail by default;
short syntax does not silently weaken validation. See `seq-bio/tests/ergonomics.rs`
for executable checks of this common path and generic error handling.
The Rust library is the only interface: the `bioferris` CLI and its JSON
agent interface were removed (decision 029).

## Explicit options and existing alignment API

```rust
use bioferris::prelude::*;

let dna = DnaSequence::try_from("ATGGCTTAA")?;
let reverse = dna.reverse_complement()?;
assert_eq!(reverse.to_string(), "TTAAGCCAT");
let rna = dna.transcribe()?;
let product = rna.translate_standard(TranslationOptions {
    mode: TranslationMode::Cds,
    ..TranslationOptions::default()
})?;
let protein = product.to_protein()?;
assert_eq!(protein.to_string(), "MA");
let run = Pairwise::gotoh().run(&protein, &protein)?;
println!("{}", run.report());
# Ok::<(), Box<dyn std::error::Error>>(())
```

This example is also an executable Rust doctest in `seq-bio/src/lib.rs`.
`bioferris::align` re-exports the complete existing alignment facade; the prelude
exposes its unchanged `Pairwise`, `Multiple`, `Graph`, records, and format helpers.
`bioferris::io` exposes the unchanged format APIs. Reuse is direct: no new alignment
algorithm implementation, runtime process, or copied scoring model.
Alignment-specific error/result aliases are named `AlignmentError` and
`AlignmentResult` in the toolkit prelude so they do not shadow the standard
generic `Result` in mixed I/O/translation workflows. The old `align` prelude is
unchanged.

The additive [file-I/O layer](file-io-v1.md) now shares this namespace: richer
`Record` values, streaming FASTA/Sanger FASTQ, gzip and explicit import
submodules. Legacy re-exports remain available. Annotated and alignment documents
are distinct from simple sequences and are not silently coerced into them.

## Supported operations

| Operation | Rust surface | Contract |
|---|---|---|
| Complement | `NucleotideSequenceExt::complement` | Complete typed DNA/RNA IUPAC alphabet, same molecule, immutable |
| Reverse complement | `NucleotideSequenceExt::reverse_complement` | Reverse then complement, no input mutation |
| Transcription | `DnaSequenceExt::transcribe` | Coding strand, T to U, returns typed RNA |
| Back-transcription | `RnaSequenceExt::back_transcribe` | U to T, returns typed DNA |
| Translation | `NucleotideSequenceExt::translate_standard` | NCBI table 1 only, reading frame zero |
| Configured translation | `ConfiguredTranslationExt::translate_with` | Every current NCBI table and forward offsets 0/1/2; [complete contract](genetic-codes-v1.md) |
| Six-frame translation and ORFs | `ReadingFramesExt` | Both strands, original-input coordinates, explicit ORF policies; [complete contract and examples](reading-frames-v1.md) |
| Annotated CDS translation | `Feature::translate_cds`, `translate_cds_with` | GenBank/EMBL extraction, annotation-selected code/frame and exact protein-reference checks; [complete workflows](cds-translation-v1.md) |
| Protein properties | `protein_properties`, `protein_properties_with`, `analyze_protein` | Verified canonical-20, nonempty linear-peptide models; average mass, charge/pI, aromaticity and GRAVY; [complete inputs/results and limits](protein-properties-v1.md) |
| Nucleotide motifs | `dna_motif`, `rna_motif`, `MotifCounts`, `MotifWeights`, `MotifScores` | Implemented and verified IUPAC matching and explicit count/PWM/PSSM models; [complete guide and scope](motifs-v1.md) |
| Alignment analysis | `AlignmentAnalysis` | Implemented and verified slices, original-source maps, gap/identity evidence and policy-based consensus; [complete guide and scope](alignment-analysis-v1.md) |
| Alignment input | `Translation::to_protein` | Validates stop-free protein; fails at first `*` |

The existing validator constructs transform results, so these APIs return its
structured `SequenceError` rather than introduce unchecked constructors or
production `unwrap`. No valid typed DNA/RNA transform is expected to fail.

## Additional analysis and interoperability

The subsequent motif, alignment-analysis, BLAST and compressed-I/O stages are
implemented and verified for their documented scope, with runnable manual
programs. The integration validation
records the checks and evidence limitations. They do not change the original P0
translation contract below or establish universal Biopython equivalence.

| Rust program | Inputs, outputs and explicit boundaries |
|---|---|
| [`motif_analysis`](../../seq-bio/examples/motif_analysis.rs) | [Complete DNA/RNA motif guide](motifs-v1.md): definite/possible/reject ambiguity, both strands, explicit counts/pseudocounts/background, binary64 models and final binary32 scan scores; no significance inference |
| [`alignment_analysis`](../../seq-bio/examples/alignment_analysis.rs) | [Complete alignment guide](alignment-analysis-v1.md): strict slices and source origins, separate gap/missing states, exact pair counts, policy-qualified consensus and refusal/full reports for annotation loss |
| [`blast_analysis`](../../seq-bio/examples/blast_analysis.rs) | [`io::blast` result readers](blast-input-v1.md): reported tabular 6/7 and bounded legacy/XML2 profiles, optional values retained without inference; not a search engine |
| [`compressed_io`](../../seq-bio/examples/compressed_io.rs) | [`io::gzip`, `io::bgzf`, `io::bgzf_index`](compressed-io-v1.md): byte-preserving finalized output, checked virtual seeks, fully verified BGZF FASTA/Sanger FASTQ sidecars and create-only publication |

Run each with `cargo bio NAME` from the repository root. The linked guides supply
literal inputs, complete programs, actual stdout, errors and budgets.

Reference comparisons qualify explicit models, versions, configurations and
retained corpora; they do not establish universal compatibility. XML processing
has no remote DTD/entity resolution and rejects unsupported semantic profiles.
BGZF blocks are bounded to 65,536 compressed and uncompressed bytes; virtual
positions are not linear offsets. Encoded bytes need not match native compressors,
and a BAM fixture used as compressed bytes does not provide a BAM parser.

### Translation policies

The standard-code methods below remain table 1/frame 0. The additive
[`TranslationConfig` API](genetic-codes-v1.md) selects other supported codes and
frames without changing existing options, defaults or standard-code errors.

`TranslationMode::Full` retains definite stops as `*` and continues translating.
`ToStop` emits no stop marker and ends before the first definite stop. `Cds`
requires complete codons, at least a distinct start and terminal stop, and no
internal definite stop. Its start must be ATG, TTG or CTG, or an IUPAC codon whose
every expansion is a valid start. The first residue becomes M; only the terminal
stop is omitted. GTG is **not** a valid start under table 1.

Synonymous ambiguity resolves to that amino acid. D/N, E/Q and I/L sets use B, Z
and J; only-stop ambiguity gives `*`. Mixed sense/stop and other unresolved sets
give X, which does not terminate `ToStop`. For example GCN -> A, RAY -> B,
SAR -> Z, WTA -> J, TAR/TRA -> `*`, TAN/NNN -> X. These rules derive from the
published code and ambiguity sets; ambiguous-CDS acceptance is a specified
all-expansions rule, not a claimed live differential result against Biopython.

`PartialCodonPolicy::Reject` is the default for incomplete trailing codons.
Selecting `Ignore` translates only complete codons and reports the original
remainder through `ignored_trailing_bases()`, even if `ToStop` stops earlier.
Complete codons after an early stop are not included in that count. CDS never
ignores incomplete codons. Error precedence is incomplete length, too-short CDS,
invalid start, missing terminal stop, then internal stop. Offsets are zero-based
nucleotide positions; no warning channel is emulated.

Empty DNA/RNA is valid for transformations and non-CDS translation. Invalid raw
symbols fail before translation through the existing sequence parser.

## Offline example with a real sequence

From the workspace root:

```bash
cargo run --offline --locked -p bioferris --example translate_fasta -- \
  seq-bio/tests/fixtures/phix174-e.fasta
```

This translates a 276-base public RefSeq CDS to its independently annotated
91-residue endolysin. Provenance and the table 11/table 1 distinction are in
`seq-bio/tests/fixtures/README.md`. No internet access or interpreter is needed
to run the example or the tests.

## Deliberate limits and differences

- Reused types normalize to uppercase and reject gaps, mixed T/U alphabets, and
  unknown characters; Biopython's more permissive string-like `Seq` is different.
- Standard methods remain table 1/frame 0; the separate configured API supports
  every current NCBI table and forward offsets 0/1/2. No custom codes, custom
  stop character, gaps, partial/undefined sequence content or mutable sequences.
- Partial codons error by default, rather than the reference's warning/truncation.
- Existing FASTA identifiers keep the complete trimmed header and output wraps
  at 80 symbols. No claim of `SeqIO` metadata/byte-output equivalence.
- P0 itself did not include annotated records or FASTQ; those now have explicit
  [file-I/O subsets](file-io-v1.md). BLAST result import has a separate verified
  [bounded input profile](blast-input-v1.md). BLAST search execution,
  network services, structures and new phylogenetic algorithms remain outside
  the implemented toolkit.
- No performance equivalence claim. Tests use documentation/specification values,
  not a runtime differential campaign or universal compatibility certificate.

## Public references

Consulted 2026-10-15; biological tables are data, not copied implementation code.

- [Biopython 1.88 tutorial](https://biopython.org/docs/1.88/Tutorial/index.html)
- [Bio.Seq public API](https://biopython.org/docs/1.88/api/Bio.Seq.html)
- [NCBI genetic code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1)
- [INSDC/DDBJ base and amino-acid codes](https://www.ddbj.nig.ac.jp/ddbj/code-e.html)
