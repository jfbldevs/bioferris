# BioFerris quick start

## Contents

- [API capability index](README.md#api-index)
- [Worked examples with complete inputs and outputs](examples/README.md)
- [Enter the repository](#enter-the-repository)
- [Use the Rust library](#use-the-rust-library)
- [Short commands for your programs](#short-commands-for-your-programs)
- [Combine translation and alignment](#combine-translation-and-alignment)
- [Work with a real FASTA file](#work-with-a-real-fasta-file)
- [Read more file formats](#read-more-file-formats)
- [Use analysis and interoperability workflows](#use-analysis-and-interoperability-workflows)
- [Use BioFerris from an AI agent](#use-bioferris-from-an-ai-agent)
- [Handle invalid input](#handle-invalid-input)
- [Browse the full reference](#browse-the-full-reference)

[Documentation index](README.md) · Current status

Start by writing a Rust program that imports `bioferris`, calls methods and combines
the results. This is the library workflow, analogous to writing a Biopython
analysis script, but in Rust: Cargo compiles and runs your `.rs` program.
The library is the only interface, for people and AI agents alike: there is no
command-line program (decision 029) and there are no Python bindings.
For complete worked cases, use the [worked example index](examples/README.md)
and the [Rust alignment recipes](examples/alignment-recipes.md). Each states its
input, command, exact output and interpretation; failure examples include exit
status and output stream.

## Enter the repository

Run the commands below from the root of this repository. You need Cargo and the toolchain pinned in
[`rust-toolchain.toml`](../rust-toolchain.toml): Rust 1.85.0, edition 2024.
The first build may download dependencies; once the toolchain and dependencies
are available locally, you can add `--offline` to Cargo commands.
Running a compiled program requires neither network access nor a Python interpreter.

## Use the Rust library

See [more file formats](#read-more-file-formats) for streaming records once you
have tried this basic program.

Input: coding DNA `ATGGCTTAA`. The complete program is already saved in
[`seq-bio/examples/my_analysis.rs`](../seq-bio/examples/my_analysis.rs):

<!-- quick-start-sequence-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATGGCTTAA")?;

    println!("DNA: {sequence}");
    println!("Complement: {}", sequence.complement()?);
    println!("Reverse complement: {}", sequence.reverse_complement()?);
    println!("RNA: {}", sequence.transcribe()?);
    println!("Protein: {}", sequence.translate_cds()?);

    Ok(())
}
```
<!-- /quick-start-sequence-program -->

Run it from the repository root:

<!-- quick-start-sequence-command -->
```bash
cargo bio my_analysis
```
<!-- /quick-start-sequence-command -->

Complete stdout (exit status `0`, empty stderr in a warning-free build):

<!-- quick-start-sequence-output -->
```text
DNA: ATGGCTTAA
Complement: TACCGAATT
Reverse complement: TTAAGCCAT
RNA: AUGGCUUAA
Protein: MA
```
<!-- /quick-start-sequence-output -->

`dna()` validates the input; the methods return new values without changing the
original sequence. Complement changes bases without reversing their order;
reverse complement also reverses that order. `transcribe()` replaces
coding-strand T with U.
`translate_cds()` validates a complete CDS and omits the terminal stop, producing
`MA`. Use `translate()` instead when you want the full result `MA*`.
These convenience methods use genetic code table 1, reading frame 0. Use the
[configured translation API](api/genetic-codes-v1.md) for tables 2/11 or explicit
forward offsets 1/2; the existing example and its defaults do not change.

For both strands and candidate coding regions, see the
[six-frame translation and ORF program](api/reading-frames-v1.md#run-a-rust-program),
including its complete input/output and error cases.

Edit the sequence or add ordinary Rust variables, loops and functions in this
file, then run the same command again. Cargo rebuilds when needed. The `?`
operator returns an error instead of continuing after invalid input: for example,
`dna("ACGU")?` rejects U in DNA, and translating `ATGA` rejects its incomplete
codon. Library errors are typed; [Handle invalid input](#handle-invalid-input)
shows how a program inspects them and what a failing `main` prints.

## Short commands for your programs

The repository includes this configuration in
[`.cargo/config.toml`](../.cargo/config.toml):

<!-- quick-start-cargo-alias -->
```toml
[alias]
bio = "run --quiet -p bioferris --example"
```
<!-- /quick-start-cargo-alias -->

Thus `cargo bio my_analysis` is a standard Cargo alias for
`cargo run --quiet -p bioferris --example my_analysis`. No plugin, shell wrapper,
new runtime or global installation is needed. Cargo's `--quiet` hides its normal
progress messages, not compiler warnings or errors.

For another program, save a `.rs` file in `seq-bio/examples/` and pass its name
without `.rs` to `cargo bio`. The file must be a valid Rust program with a `main`
function. A nonexistent example is a Cargo error, not a biological API error.
For instance, this command names a file that does not exist:

<!-- quick-start-unknown-example-command -->
```bash
cargo bio nonexistent_analysis
```
<!-- /quick-start-unknown-example-command -->

Complete stderr with the pinned Cargo version (exit status `101`, empty stdout):

<!-- quick-start-unknown-example-output -->
```text
error: no example target named `nonexistent_analysis`.
Available example targets:
    alignment_analysis
    alignment_recipes
    alignment_workflow
    annotated_records
    binary_index
    blast_analysis
    blast_stream
    compressed_io
    genetic_codes
    index_benchmark
    indexed_records
    invalid_input
    motif_analysis
    my_analysis
    protein_analysis
    qa_workload
    quick_start
    reading_frames
    sequence_operations
    sequence_statistics
    stream_annotations
    translate_annotations
    translate_fasta

```
<!-- /quick-start-unknown-example-output -->

Use an existing example name or add the corresponding Rust program. This Cargo
diagnostic depends on the toolchain and available targets; it is not a BioFerris
library error.

The alias selects named examples in this package; it does not execute arbitrary
Rust files by path. Run it inside this repository so Cargo finds the configuration.

Cargo options go before `--`; arguments for your program go after it. For example,
the [FASTA workflow](#work-with-a-real-fasta-file) uses
`cargo bio translate_fasta -- seq-bio/tests/fixtures/phix174-e.fasta`.
For a reproducible offline run after dependencies have been downloaded, use
`cargo bio my_analysis --locked --offline`.

In a separate analysis project already configured to depend on this local
library, your program can live in `src/main.rs` and run with plain `cargo run`
(or `cargo run -q`). That is a different project layout, not what `cargo run`
alone selects in this multi-package workspace. Both the package and Rust import
are named `bioferris`; the package is not published to crates.io yet. Its source
directory retains the legacy name `seq-bio/`. See the
[identity and migration guide](api/bioferris-identity.md).
`cargo run my_analysis` does not select an example: it passes `my_analysis` as a
program argument. See the official [Cargo run documentation](https://doc.rust-lang.org/cargo/commands/cargo-run.html)
and [alias configuration](https://doc.rust-lang.org/cargo/reference/config.html#alias).

## Combine translation and alignment

The workspace package is named `bioferris` and imported as `bioferris`.
This program demonstrates the complete workflow:

Input: coding DNA `ATGGCTTAA`, alignment target `GATTACA`, and query `GCATGCT`.
The exact program is saved as
[`seq-bio/examples/quick_start.rs`](../seq-bio/examples/quick_start.rs). Run it with:

<!-- quick-start-alignment-command -->
```bash
cargo bio quick_start
```
<!-- /quick-start-alignment-command -->

<!-- quick-start-alignment-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let coding = dna("ATGGCTTAA")?;
    assert_eq!(coding.reverse_complement()?.to_string(), "TTAAGCCAT");
    assert_eq!(coding.transcribe()?.to_string(), "AUGGCUUAA");

    let product = coding.translate_cds()?;
    assert_eq!(product.as_str(), "MA");
    let translated_protein = product.to_protein()?;
    assert_eq!(translated_protein, protein("MA")?);
    println!("{product}");

    let target = dna("GATTACA")?;
    let query = dna("GCATGCT")?;
    let run = Pairwise::gotoh()
        .scores(2, -1)
        .gaps(-5, -1)
        .run(&target, &query)?;
    println!("{}", run.report());
    Ok(())
}
```
<!-- /quick-start-alignment-program -->

Complete stdout (exit status `0`, empty stderr):

<!-- quick-start-alignment-output -->
```text
MA
Algorithm     : Gotoh
Mode          : global affine
Scoring       : match=2 mismatch=-1 gap_open=-5 gap_extend=-1
Sequences     : First=7 Second=7
Range         : First 1..7, Second 1..7
Score         : 2
Columns       : 7
Identity      : 3/7 (42.86%)
Matches       : 3
Substitutions : 4
Insertions    : 0 residues, 0 runs
Deletions     : 0 residues, 0 runs
CIGAR         : 1=2X1=1X1=1X

First   1  GATTACA  7
           |..|.|.
Second  1  GCATGCT  7

Legend: | match, . substitution, blank gap
```
<!-- /quick-start-alignment-output -->

The first line is the translated protein. The following report shows the global
alignment of the two DNA inputs, not an alignment of that protein. Three matches
and four substitutions give `3 × 2 + 4 × (−1) = 2`.

The `?` operator propagates typed errors. `Translation` may contain `*`, but
`ProteinSequence` rejects it: `to_protein()` validates the conversion and never
silently removes stops. No separate sequence model is needed to pass results to
`Pairwise`, `Multiple` or `Graph`.

`bioferris` is not published to crates.io yet (`publish = false`); these examples
use the local workspace. Do not assume `cargo add bioferris` installs this project.

The alignment API also supports custom scoring, Needleman–Wunsch, Smith–Waterman,
Myers, WFA, X-drop, multiple alignment and graphs. The
[alignment recipes](examples/alignment-recipes.md) give each with literal input
and complete output; the [alignment API](api/rust-api-v1.md) lists the options.

## Work with a real FASTA file

The existing example reads DNA from FASTA, translates each record as a complete
CDS and emits proteins as FASTA:

Input file [`phix174-e.fasta`](../seq-bio/tests/fixtures/phix174-e.fasta):

<!-- quick-start-fasta-input -->
```fasta
>NC_001422.1:568-843 Escherichia phage phiX174, complete genome
ATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGC
TTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAAC
ATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAG
GAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAGTGA
```
<!-- /quick-start-fasta-input -->

Command:

```bash
cargo bio translate_fasta -- \
  seq-bio/tests/fixtures/phix174-e.fasta
```

Expected stdout (exit status `0`):

```fasta
>NC_001422.1:568-843 Escherichia phage phiX174, complete genome
MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNLRKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE
```

The included sample contains 276 nucleotides and produces a 91-residue protein.
Its [provenance and annotation](../seq-bio/tests/fixtures/README.md) are documented.
You can replace the path with a FASTA file containing complete CDS sequences
valid under table 1. The example exits with an error on an invalid record; output
is not atomic if a later record fails after earlier records have been emitted.

This example uses the legacy eager reader and its complete trimmed-header
identifiers. It remains unchanged. The new streaming readers, FASTQ qualities and
annotated input are covered in the next section.

## Read more file formats

The [file-I/O guide](api/file-io-v1.md) gives complete Rust programs, literal
FASTA/FASTQ input, preserved fields/qualities and actual output. It also covers
gzip input and checked conversion between formats. Use
`bioferris::io::fasta::<Dna>(path)?` or `fastq::<Dna>(path)?` inside your program.

For other inputs, see the Rust library
[alignment readers](api/alignment-input-v1.md),
[GenBank/EMBL readers](api/annotated-input-v1.md),
[annotated writing and transform recipes](examples/annotated-output-recipes.md), and
[graph/result readers](api/graph-result-input-v1.md). These implement explicit
subsets and preserve or reject unsupported information; do not infer universal
format coverage from a filename extension. Alignment analysis reuses the
alignment reader profiles.

## Use analysis and interoperability workflows

The four new manual programs use the Rust library directly. Their guides place
literal inputs, complete source, actual results and errors together:

| Program from the repository root | Complete workflow |
|---|---|
| `cargo bio motif_analysis` | [DNA/RNA IUPAC and count/PWM/PSSM analysis](api/motifs-v1.md), including reverse scores, pseudocount/background and explicit precision |
| `cargo bio alignment_analysis` | [Alignment slices, original-source maps, gap/identity statistics and consensus](api/alignment-analysis-v1.md), including missing-data refusal and complete annotation-loss reports |
| `cargo bio blast_analysis` | [Reported BLAST tabular/XML results](api/blast-input-v1.md), preserved endpoints and absent fields; no BLAST search is run |
| `cargo bio compressed_io` | [Gzip/BGZF output and verified BGZF retrieval](api/compressed-io-v1.md), preserved case/quality, explicit finalization and create-only files |

These stages are implemented and verified for their documented scope, not
universal Biopython equivalents. The integration validation
records the checks and evidence limits.
Protein properties have a separate [verified model guide](api/protein-properties-v1.md).

## Use BioFerris from an AI agent

AI agents use the library the same way people do: write or edit a Rust program,
run it with Cargo and read what it prints. Invoke Cargo with an argument array
such as `["cargo", "bio", "my_analysis", "--locked", "--offline"]`, with the
repository root as the working directory; do not build shell commands by
interpolating unescaped sequences. Pass inputs as program arguments after `--`,
as in `["cargo", "bio", "invalid_input", "--locked", "--offline", "--", "ATGA"]`.

Check the exit status and both streams. Status 0 means the program finished;
status 1 with `Error: ` on stderr means `main` returned an error, as
[shown below](#handle-invalid-input); status 101 with `error:` on stderr is a
Cargo failure such as the [unknown example](#short-commands-for-your-programs)
above. Decide inside the program which errors to handle: match the typed error,
or use its stable `code()` where the error type provides one. Message text is
not a stable control interface. When another tool must parse results, print them
in a format you choose, such as tab-separated lines, or use the library's
writers: FASTA/FASTQ, alignment JSON `align-output-v1`, GFA and the other
[export formats](api/phylogenetic-output-v1.md).

For repeated record access, use the
[indexing Rust program](examples/sequence-index-recipes.md). `SequenceIndex`
builds and persists indexes for stable plain FASTA/Sanger FASTQ. Reopening checks
the entire source, so keep an index open when making many lookups. That API
remains plain-file-only. The separate [BGZF workflow](api/compressed-io-v1.md)
uses `BgzfSequenceIndex`; it fully verifies the compressed source and entries on
open. Ordinary gzip is not BGZF random access.

The former `bioferris` command-line program and its JSON requests were removed
on 2026-10-21; their historical contract records what
they did.

## Handle invalid input

Input: DNA `ATGA`, one base longer than a complete codon. The complete program is
saved in [`seq-bio/examples/invalid_input.rs`](../seq-bio/examples/invalid_input.rs):

<!-- quick-start-invalid-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // An optional argument replaces the default input; `?` ends the program on
    // input that is not DNA.
    let input = std::env::args().nth(1).unwrap_or_else(|| "ATGA".to_owned());
    let sequence = dna(&input)?;

    match sequence.translate_cds() {
        Ok(protein) => println!("CDS {sequence}: {protein}"),
        Err(error @ TranslationError::PartialCodon { .. }) => {
            println!("Rejected {sequence}: {error}");
        }
        Err(error) => return Err(error.into()),
    }

    let translation = sequence.translate_standard(TranslationOptions {
        partial_codons: PartialCodonPolicy::Ignore,
        ..TranslationOptions::default()
    })?;
    println!(
        "Remainder ignored: {translation} (ignored bases: {})",
        translation.ignored_trailing_bases()
    );
    Ok(())
}
```
<!-- /quick-start-invalid-program -->

Run it from the repository root:

<!-- quick-start-invalid-command -->
```bash
cargo bio invalid_input
```
<!-- /quick-start-invalid-command -->

Complete stdout (exit status `0`, empty stderr):

<!-- quick-start-invalid-output -->
```text
Rejected ATGA: sequence length 4 leaves 1 incomplete codon bases
Remainder ignored: M (ignored bases: 1)
```
<!-- /quick-start-invalid-output -->

The program matches the typed `TranslationError::PartialCodon` instead of parsing
the message, so it continues. Discarding the trailing base must be requested
explicitly with `PartialCodonPolicy::Ignore`; that policy translates the complete
codon `ATG` and never relaxes CDS validation. With `ATGGCTTAA` as the argument,
the same program prints `CDS ATGGCTTAA: MA` and `Remainder ignored: MA* (ignored
bases: 0)`.

An error the program does not handle ends it through `?`. RNA `ACGU` is not DNA:

<!-- quick-start-invalid-dna-command -->
```bash
cargo bio invalid_input -- ACGU
```
<!-- /quick-start-invalid-dna-command -->

Complete stderr (exit status `1`, empty stdout):

<!-- quick-start-invalid-dna-output -->
```text
Error: SequenceError { alphabet: "DNA", index: 3, byte: 85 }
```
<!-- /quick-start-invalid-dna-output -->

`dna()` rejected byte 85 (`U`) at zero-based index 3, before any translation ran.
Rust printed the error's `Debug` form because `main` returned it. Handle the
errors your analysis expects, as above, and let `?` end the program for the rest.

## Browse the full reference

```bash
cargo doc -p bioferris --no-deps --locked --open
```

This generates and opens the Rust API reference at
`target/doc/bioferris/index.html`. To generate the alignment API references too,
run `cargo doc --workspace --no-deps --locked`.

- [Biological contract: ambiguity, CDS, policies and differences](api/seq-bio-v1.md).
- Historical record of the removed CLI and JSON interface.
- [Alignment API, formats and advanced options](api/rust-api-v1.md).
- Implemented and planned capabilities.
- Tests and recorded results.

Every new feature must update the relevant documentation alongside its
implementation, following the documentation policy.
