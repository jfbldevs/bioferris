# Six-frame translation and ORF discovery

## Capability index

| Capability | Rust API |
|---|---|
| All frames on both strands | `translate_six_frames`, `translate_six_frames_with` |
| Explicitly configured ORF candidates | `find_orfs`, `find_orfs_with(OrfOptions)` |
| Coordinates, stop positions and peptides | `FrameTranslation`, `OpenReadingFrame` getters |

The `translate-six-frames` and `find-orfs` CLI operations were removed with the
`bioferris` CLI (decision 029).

## Contents

- [Run a Rust program](#run-a-rust-program)
- [Policies and coordinates](#policies-and-coordinates)
- [Alternative starts and partial candidates](#alternative-starts-and-partial-candidates)
- [Errors and limits](#errors-and-limits)
- [Compatibility and verification](#compatibility-and-verification)

[API index](../README.md#api-index) · [Quick start](../quick-start.md) ·
[Algorithm specification](../algorithms/reading-frames-v1.md)

## Run a Rust program

Run commands from the repository root with its pinned Rust toolchain.
The prelude exports `ReadingFramesExt`, configuration enums and result types.
The complete program below is
[`reading_frames.rs`](../../seq-bio/examples/reading_frames.rs). It uses literal
DNA inputs `ATGAAATAA`, its reverse complement `TTATTTCAT`, the deliberately
ambiguous `ATGNNNTAA` and the unterminated `ATGAA`, plus RNA `GUGAAAUAA`.
Coordinates always refer to the original input.

<!-- reading-frame-rust-program -->
```rust
//! Translate all six frames and discover strand-aware, stop-terminated ORFs.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATGAAATAA")?;
    println!("DNA: {sequence}");
    for frame in sequence.translate_six_frames() {
        println!(
            "{:?} frame {} {:?}: {} (leading {}, trailing {})",
            frame.strand(),
            frame.reading_frame().offset(),
            frame.range(),
            frame.peptide(),
            frame.skipped_leading_bases(),
            frame.ignored_trailing_bases()
        );
    }
    for input in ["ATGAAATAA", "TTATTTCAT"] {
        for orf in dna(input)?.find_orfs()? {
            println!(
                "ORF {input}: {:?} frame {} {:?}, stop {:?}, protein {}",
                orf.strand(),
                orf.reading_frame().offset(),
                orf.range(),
                orf.stop_range(),
                orf.to_protein()?
            );
        }
    }
    match dna("ATGNNNTAA")?.find_orfs() {
        Err(error) => println!("Ambiguous input: {error}"),
        Ok(_) => return Err("expected ambiguous ORF input to be rejected".into()),
    }

    let bacterial = OrfOptions {
        genetic_code: GeneticCode::BacterialArchaealPlastid,
        strands: StrandSelection::Forward,
        start_policy: OrfStartPolicy::GeneticCode,
        ..OrfOptions::default()
    };
    for orf in rna("GUGAAAUAA")?.find_orfs_with(bacterial)? {
        println!(
            "GUG start, table 11: {:?}, stop {:?}, protein {}",
            orf.range(),
            orf.stop_range(),
            orf.to_protein()?
        );
    }
    let partial = OrfOptions {
        strands: StrandSelection::Forward,
        stop_policy: OrfStopPolicy::AllowPartial,
        ..OrfOptions::default()
    };
    for orf in dna("ATGAA")?.find_orfs_with(partial)? {
        println!(
            "Unterminated ATGAA: {:?}, stop {:?}, protein {}",
            orf.range(),
            orf.stop_range(),
            orf.to_protein()?
        );
    }
    let no_results = OrfOptions {
        max_orfs: 0,
        ..OrfOptions::default()
    };
    if let Err(error) = dna("ATGAAATAA")?.find_orfs_with(no_results) {
        println!("Zero result budget: {error}");
    }
    Ok(())
}
```

Complete invocation:

<!-- reading-frame-rust-command -->
```bash
cargo bio reading_frames
```

Complete stdout, exit 0 and no stderr:

<!-- reading-frame-rust-output -->
```text
DNA: ATGAAATAA
Forward frame 0 0..9: MK* (leading 0, trailing 0)
Forward frame 1 1..7: *N (leading 1, trailing 2)
Forward frame 2 2..8: EI (leading 2, trailing 1)
Reverse frame 0 0..9: LFH (leading 0, trailing 0)
Reverse frame 1 2..8: YF (leading 1, trailing 2)
Reverse frame 2 1..7: IS (leading 2, trailing 1)
ORF ATGAAATAA: Forward frame 0 0..6, stop Some(6..9), protein MK
ORF TTATTTCAT: Reverse frame 0 3..9, stop Some(0..3), protein MK
Ambiguous input: ORF scanning requires concrete bases; found N at nucleotide 3
GUG start, table 11: 0..6, stop Some(6..9), protein MK
Unterminated ATGAA: 0..3, stop None, protein M
Zero result budget: ORF result count exceeds limit 0
```

Six-frame translation retains stops as `*`. ORFs exclude their terminal stop from
both peptide and coding range, recording its coordinates separately. `to_protein()`
constructs the existing stop-free protein type for alignment. Finding no candidate
is a successful empty vector, not an error. An ORF candidate is not a predicted gene.

## Policies and coordinates

`translate_six_frames()` uses table 1; `translate_six_frames_with(GeneticCode)`
also supports tables 2 and 11. It translates complete codons without initiation
recoding, resolves IUPAC symbols through the existing translation engine, and
reports actual skipped-leading and incomplete-trailing bases in each frame.
Even empty inputs have six entries. Leading counts are capped at input length;
the selected frame is still offset 0/1/2. This fixed partial-codon policy does not
change the existing single-frame methods' default rejection of partial codons.

`find_orfs()` uses the following public `OrfOptions` fields. Configure individual
fields with `find_orfs_with(OrfOptions { ..., ..OrfOptions::default() })`.

| Field | Rust choices | Default |
|---|---|---|
| `genetic_code` | Any current NCBI `GeneticCode` except tables 27, 28 and 31 | `Standard` (1) |
| `strands` | `StrandSelection::{Both, Forward, Reverse}` | `Both` |
| `start_policy` | `OrfStartPolicy::{Atg, GeneticCode, AnySense}` | `Atg` |
| `stop_policy` | `OrfStopPolicy::{Require, AllowPartial}` | `Require` |
| `min_amino_acids` | Positive `usize` | 1 |
| `max_orfs` | `usize`; zero permits only an empty result | 100000 |

Start-based modes retain only the earliest eligible
start before each stop and recode its first residue to M. `AnySense` emits the
whole non-stop segment without recoding; it does not require an initiator.
`AllowPartial` also permits an unterminated final candidate, but never invents
a missing start or includes a partial codon. Minimum length excludes stop and
is inclusive. The scanner requires concrete A/C/G/T DNA or A/C/G/U RNA.

Ranges are zero-based, half-open and ascending in the original input.
Reverse-oriented `[s, e)` maps to `[n-e, n-s)` for input length `n`.
`reading_frame()` is offset 0/1/2 in the oriented sequence, not a signed frame
or a genomic modulo-three value. Order is forward frames 0/1/2, then reverse
frames 0/1/2; candidates within a frame follow translation direction. Reverse
candidates therefore appear in decreasing genomic position. No sorting by length
or silent nested-start enumeration occurs.

## Alternative starts and partial candidates

The program's RNA `GUGAAAUAA` uses table 11, the forward strand and
table-specific initiation (`OrfStartPolicy::GeneticCode`). The supported
bacterial initiation set accepts GUG and recodes that initial residue to M, so
the candidate `0..6` translates to `MK`. `AnySense` instead retains its ordinary
V assignment; default `Atg` does not use GUG as a start.

DNA `ATGAA` has one complete start codon and two trailing bases. Explicit
`OrfStopPolicy::AllowPartial` permits that one-residue candidate `0..3` with no
stop (`stop None`); the two trailing bases stay outside it. Default `Require`
would return no ORFs, not silently relax the policy.


## Errors and limits

ORF discovery returns `Err(OrfError)` and never a partial list:

| Variant | Cause | Recovery |
|---|---|---|
| `AmbiguousBase { offset, byte }` | The first ambiguous base, such as N in the program's `ATGNNNTAA` at nucleotide 3 | Resolve the input uncertainty, or use six-frame translation to inspect its ambiguous products; the scanner does not guess the base |
| `ResultLimitExceeded { limit }` | Another qualifying candidate would exceed `max_orfs`, as the program's zero budget shows | Increase the explicit budget if appropriate; no truncated success is returned |
| `ConditionalStops { genetic_code }` | NCBI table 27, 28 or 31 was selected: no codon in it always ends translation, so ORF boundaries are undefined | Translate with that table instead, or choose the table that defines the ORFs of interest |
| `InvalidMinimumLength` | `min_amino_acids` is zero | Use a positive minimum length |

Supported genetic codes are the `GeneticCode` variants; converting a number that
is not a current NCBI table with `GeneticCode::try_from` fails before any scan.
Six-frame translation accepts every table, including 27, 28 and 31. Invalid sequence bytes
also fail earlier, when `dna()` or `rna()` builds the sequence.

The scanner retains at most `max_orfs` products, but this is not an exact memory
budget. Input copying, reverse complements and retained peptides are linear in
sequence length. No circular scan, all-nested-start mode, ambiguous ORF scan,
annotation-aware CDS translation or recoding exception is implemented here.
Inputs and the retained products must fit in memory.

## Compatibility and verification

The [algorithm contract](../algorithms/reading-frames-v1.md) records the public
NCBI and Biopython 1.88 references, independent implementation and exact ordering.
These tests verify the stated BioFerris subset, not universal Biopython fidelity
or biological gene prediction. No live Biopython execution or closed scientific
campaign is part of this increment. Existing single-frame APIs, alignment kernels,
defaults and historical evidence are unchanged.

Rust tests check independent codon/ORF expectations, strand/frame coordinates,
boundaries and real-record input. Recipe tests compare this guide's complete
program and run the actual Cargo alias, comparing complete output rather than
silently updating snapshots. The validation record
states the commands executed for the original increment, including the removed
CLI, and their results.
