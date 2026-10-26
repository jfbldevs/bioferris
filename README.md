# BioFerris

Rust-native bioinformatics, from sequences to workflows.

BioFerris is a Rust library for sequence analysis that follows the common
workflows of Biopython while making biological assumptions explicit. It reads
and writes sequence, annotation, alignment and search-result files, translates
coding sequences, analyses sequences and proteins, and aligns sequences with
pairwise, multiple and partial-order methods, all from ordinary Rust programs.

## Contents

- [What it offers](#what-it-offers)
- [Design](#design)
- [Requirements](#requirements)
- [Install](#install)
- [First program](#first-program)
- [Documentation](#documentation)
- [Multiple-alignment methods and their provenance](#multiple-alignment-methods-and-their-provenance)
- [Status and platforms](#status-and-platforms)
- [License](#license)

## What it offers

| Area | Capabilities |
|---|---|
| Sequences | Typed DNA, RNA and protein sequences; complement, reverse complement, transcription; checked slicing, concatenation and literal search |
| Translation | Every current NCBI genetic code, forward frames, six-frame translation, open reading frames, and annotation-aware CDS translation |
| Files | Streaming FASTA and Sanger FASTQ, plain or gzip; gzip and BGZF output; persistent record indexes with JSON or compact binary sidecars |
| Annotated records | GenBank and EMBL reading, record by record or whole, with features and qualifiers; checked writing and transformations |
| Search results | BLAST tabular (outfmt 6 and 7), legacy XML and XML2, collected or streamed one query at a time |
| Analyses | Composition, GC content and k-mers; protein mass, charge, isoelectric point, aromaticity and GRAVY; IUPAC motifs and position-specific scoring matrices |
| Pairwise alignment | Needleman–Wunsch, Smith–Waterman, Gotoh with every end-gap mode, Myers edit distance and search, wavefront alignment, X-drop extension, and linear-memory score-only kernels |
| Multiple alignment | MUSCLE 5 (Align and Super5), Clustal Omega, six historical MAFFT strategies and two MAFFT 7.526 compatibility profiles |
| Graph alignment | POASTA partial-order alignment and graph construction, with GFA output |
| Alignment files | Import of aligned FASTA, CLUSTAL, PHYLIP, Stockholm and NEXUS; export to the same formats, JSON and graphs |

The [API index](docs/README.md#api-index) lists every entry point.

## Design

- **Typed molecules.** A `Sequence<Dna>` cannot be passed where protein is
  expected, and every symbol is checked against its alphabet when it is read.
- **Explicit limits and errors.** Readers enforce declared size limits, and
  errors are typed values that name the record and line they come from.
- **Refusal instead of guessing.** An operation that cannot honour everything
  that changes its result refuses with a reason. For example, a coding sequence
  with fuzzy boundaries or an annotated exception is not translated, rather than
  translated into a peptide that may contradict the record.
- **Safe, self-contained Rust.** No `unsafe` code, no foreign-function
  interface, no external programs and no Python at run time. Named outputs are
  created new and never overwrite a file.

## Requirements

- Rust 1.85 or newer, installed with [rustup](https://rustup.rs); the library uses
  the 2024 edition. On Windows, rustup also asks for the Visual Studio C++ build
  tools.
- Internet access the first time you build, to download dependencies; after
  that, builds work offline.

## Install

BioFerris is not yet on crates.io. Add it to your project's `Cargo.toml` from
its repository, together with the patch described below:

```toml
[dependencies]
bioferris = { git = "https://github.com/jfbldevs/bioferris" }

[patch.crates-io]
mafft-core = { git = "https://github.com/jfbldevs/bioferris" }
```

or, from a local copy of this repository:

```toml
[dependencies]
bioferris = { path = "/path/to/bioferris/seq-bio" }

[patch.crates-io]
mafft-core = { path = "/path/to/bioferris/vendor/mafft-core" }
```

The patch matters for the two MAFFT 7.526 compatibility profiles. They run the
fork of `mafft-core` in `vendor/mafft-core`, which changes nucleotide gap
scaling and ambiguity-code scoring, and Cargo applies a patch only from the
project being built. Without it, your project gets the unmodified `mafft-core`
from crates.io, and those profiles may give different nucleotide results; every
other feature is unaffected.

The crate and its Rust import are both `bioferris`; the source directory is
`seq-bio/`.

## First program

Input: the coding DNA `ATGGCTTAA`. Create a project with `cargo new
my-analysis`, add the dependency above, and write `src/main.rs`:

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

Run it with `cargo run --release`. Complete standard output, exit status 0:

```text
DNA: ATGGCTTAA
Complement: TACCGAATT
Reverse complement: TTAAGCCAT
RNA: AUGGCUUAA
Protein: MA
```

`dna()` validates its input, and each method returns a new value.
`translate_cds()` checks a complete coding sequence and omits the final stop
codon, giving `MA`; `translate()` would give `MA*`. Invalid input becomes a typed
error: `dna("ACGU")` rejects U in DNA. The
[quick start](docs/quick-start.md#handle-invalid-input) shows how a program
inspects such errors and continues.

The same program is in this repository as
[`seq-bio/examples/my_analysis.rs`](seq-bio/examples/my_analysis.rs). From the
repository root, `cargo bio my_analysis` runs it, and `cargo bio <name>` runs
any of the other examples in `seq-bio/examples/`.

## Documentation

| Guide | Content |
|---|---|
| [Quick start](docs/quick-start.md) | First programs, real files and error handling |
| [Documentation index](docs/README.md) | API index and every guide |
| [Worked examples](docs/examples/README.md) | Literal inputs, complete programs, exact outputs and their interpretation |
| [Alignment guide](docs/api/rust-api-v1.md) | Pairwise, multiple and graph alignment, and which multiple-alignment profile to choose |
| [Algorithm specifications](docs/README.md#algorithm-specifications) | What each method computes, its parameters and its limits |
| [Testing](docs/testing/README.md) | How the library is tested, and what the tests do and do not establish |

The documented outputs are checked by the test suite against real runs. For the
type and method reference, run `cargo doc -p bioferris --no-deps --open` from the
repository root.

## Multiple-alignment methods and their provenance

The MUSCLE 5, Clustal Omega and historical MAFFT strategies are implemented in
Rust from their publications, supplementary material, manuals and public data;
they do not contain or derive from the native programs' source code. Where a
publication leaves a value unstated, the specification declares the choice
made. They are validated against the native programs' recorded outputs and
reference alignments, and they may differ from the native programs: see
[which profile to choose](docs/api/rust-api-v1.md#choose-a-multiple-alignment-profile)
for measured accuracy and speed.

The two MAFFT 7.526 compatibility profiles instead run a fork of the rust-MAFFT
port, vendored in `vendor/mafft-core` under its own licenses, MIT AND
BSD-3-Clause, with its provenance recorded.

## Status and platforms

BioFerris covers the capabilities above within the scopes their guides
document; it is not a complete replacement for Biopython. It is tested on Linux
x86_64 with Rust 1.85. Windows, macOS and ARM have not been qualified yet. A
manuscript describing the library and its evaluation is in preparation.

## License

Licensed under either the [Apache License, Version 2.0](LICENSE-APACHE) or the
[MIT license](LICENSE-MIT), at your option. The vendored `mafft-core` is licensed
MIT AND BSD-3-Clause; distributions must also preserve its terms.
