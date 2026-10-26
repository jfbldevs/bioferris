# Genetic codes and forward reading frames

## Contents

- [Capabilities and defaults](#capabilities-and-defaults)
- [Run a Rust analysis](#run-a-rust-analysis)
- [Supported genetic codes](#supported-genetic-codes)
- [Translation modes and ambiguity](#translation-modes-and-ambiguity)
- [Frames and result metadata](#frames-and-result-metadata)
- [Errors and input validation](#errors-and-input-validation)
- [Compatibility, limits and verification](#compatibility-limits-and-verification)

[API capability index](../README.md#api-index) · [Quick start](../quick-start.md)

## Capabilities and defaults

| Capability | Rust access |
|---|---|
| Standard table 1, frame 0 | Existing `translate`, `translate_to_stop`, `translate_cds`, `translate_standard` |
| Any current NCBI table: 1–6, 9–16 and 21–33 | `translate_with(TranslationConfig)` |
| Forward offset 0, 1 or 2 | `TranslationConfig::reading_frame` |
| Typed result and leading/trailing metadata | `Translation` getters |

For simultaneous translation of both strands and ORF discovery, use the separate
[reading-frame workflow API](reading-frames-v1.md). It shares these genetic codes
without changing the single-frame behavior documented here.

`ConfiguredTranslationExt`, `TranslationConfig`, `GeneticCode`, `ReadingFrame`
and the existing translation types are exported by `bioferris::prelude::*`.
This guide introduces manual use as an ordinary Rust program. The CLI
`--genetic-code`/`--frame` selectors and JSON requests were removed with the
`bioferris` CLI (decision 029). To translate every record of a
FASTA or FASTQ file, read it with the [file readers](file-io-v1.md) and call
`translate_with` on each record's sequence; the
[`translate_fasta`](../../seq-bio/examples/translate_fasta.rs) example shows the
reading loop.

`TranslationConfig::default()` means standard genetic code 1, forward frame 0,
`TranslationMode::Full` and `PartialCodonPolicy::Reject`. All four fields are
public. Configuration is explicit: no organism, strand or annotation qualifier
is inferred. DNA and RNA use equivalent codon assignments, with T and U
respectively. Typed sequences normalize lowercase and validate the entire input.

## Run a Rust analysis

Run from the repository root using its pinned Rust toolchain. The
complete program is checked in as
[`genetic_codes.rs`](../../seq-bio/examples/genetic_codes.rs). Its literal inputs
are DNA `ATATGAAGA`, `GTGAAATAA`, `CCATGAAAT`, `AATGTGATAA` and `AC`, plus RNA
`GUGAAAUAA`. It compares tables, validates a bacterial CDS, selects a frame,
handles an internal stop and translates an empty suffix.

<!-- genetic-codes-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATATGAAGA")?;
    println!("DNA: {sequence}");
    for genetic_code in [
        GeneticCode::Standard,
        GeneticCode::VertebrateMitochondrial,
        GeneticCode::BacterialArchaealPlastid,
    ] {
        let translation = sequence.translate_with(TranslationConfig {
            genetic_code,
            ..TranslationConfig::default()
        })?;
        println!(
            "Table {} ({}): {translation}",
            genetic_code.id(),
            genetic_code.name()
        );
    }

    let cds_config = TranslationConfig {
        genetic_code: GeneticCode::BacterialArchaealPlastid,
        mode: TranslationMode::Cds,
        ..TranslationConfig::default()
    };
    println!(
        "Table 11 CDS GTGAAATAA: {}",
        dna("GTGAAATAA")?.translate_with(cds_config)?
    );
    println!(
        "Table 11 RNA GUGAAAUAA: {}",
        rna("GUGAAAUAA")?.translate_with(cds_config)?
    );

    let framed = dna("CCATGAAAT")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::Two,
        partial_codons: PartialCodonPolicy::Ignore,
        ..TranslationConfig::default()
    })?;
    println!(
        "Frame {} CCATGAAAT: {framed}; skipped: {}; ignored: {}",
        framed.reading_frame().offset(),
        framed.skipped_leading_bases(),
        framed.ignored_trailing_bases()
    );
    println!("Stop-free protein: {}", framed.to_protein()?);

    let stopped = dna("AATGTGATAA")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::One,
        mode: TranslationMode::Cds,
        ..TranslationConfig::default()
    });
    match stopped {
        Err(TranslationError::InternalStop { offset }) => {
            println!("CDS error: internal stop at original nucleotide {offset}");
        }
        Err(error) => return Err(error.into()),
        Ok(_) => return Err("expected an internal-stop error".into()),
    }

    let empty = dna("AC")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::Two,
        ..TranslationConfig::default()
    })?;
    println!(
        "Frame 2 AC: length {}; skipped: {}",
        empty.len(),
        empty.skipped_leading_bases()
    );
    Ok(())
}
```

Complete invocation:

<!-- genetic-codes-command -->
```bash
cargo bio genetic_codes
```

Complete observed stdout, exit 0, no stderr:

<!-- genetic-codes-output -->
```text
DNA: ATATGAAGA
Table 1 (Standard): I*R
Table 2 (Vertebrate Mitochondrial): MW*
Table 11 (Bacterial, Archaeal and Plant Plastid): I*R
Table 11 CDS GTGAAATAA: MK
Table 11 RNA GUGAAAUAA: MK
Frame 2 CCATGAAAT: MK; skipped: 2; ignored: 1
Stop-free protein: MK
CDS error: internal stop at original nucleotide 4
Frame 2 AC: length 0; skipped: 2
```

Table 2 changes all three codons of `ATATGAAGA`. In CDS mode, table 11 recodes
the initial GTG/GUG to M and omits the final TAA/UAA stop. The framed sequence
skips `CC`, translates `ATGAAA` to `MK`, and explicitly ignores one trailing T.
The internal-stop error refers to position 4 in the original `AATGTGATAA`, not
position 3 in its selected suffix. That error is handled, so the example exits
successfully. A frame equal to input length produces an empty non-CDS product.

## Supported genetic codes

These are the published NCBI assignments and initiation sets, not an automatic
choice based on a species name. Codons below use DNA notation; replace T with U
for RNA. `GeneticCode::id()` returns the `u16` table number and `name()` its NCBI
name; `GeneticCode::ALL` lists every table in identifier order.

| NCBI ID and Rust variant | Assignments that differ from table 1 | Allowed CDS start codons | Codons that end translation |
|---|---|---|---|
| 1, `Standard` | None | TTG, CTG, ATG | TAA, TAG, TGA |
| 2, `VertebrateMitochondrial` | TGA → W; ATA → M; AGA/AGG → stop | ATT, ATC, ATA, ATG, GTG | TAA, TAG, AGA, AGG |
| 3, `YeastMitochondrial` | TGA → W; CTT/CTC/CTA/CTG → T; ATA → M | ATA, ATG, GTG | TAA, TAG |
| 4, `MoldProtozoanCoelenterateMitochondrial` | TGA → W | TTA, TTG, CTG, ATT, ATC, ATA, ATG, GTG | TAA, TAG |
| 5, `InvertebrateMitochondrial` | TGA → W; ATA → M; AGA/AGG → S | TTG, ATT, ATC, ATA, ATG, GTG | TAA, TAG |
| 6, `CiliateDasycladaceanHexamitaNuclear` | TAA/TAG → Q | ATG | TGA |
| 9, `EchinodermFlatwormMitochondrial` | TGA → W; AAA → N; AGA/AGG → S | ATG, GTG | TAA, TAG |
| 10, `EuplotidNuclear` | TGA → C | ATG | TAA, TAG |
| 11, `BacterialArchaealPlastid` | None | TTG, CTG, ATT, ATC, ATA, ATG, GTG | TAA, TAG, TGA |
| 12, `AlternativeYeastNuclear` | CTG → S | CTG, ATG | TAA, TAG, TGA |
| 13, `AscidianMitochondrial` | TGA → W; ATA → M; AGA/AGG → G | TTG, ATA, ATG, GTG | TAA, TAG |
| 14, `AlternativeFlatwormMitochondrial` | TAA → Y; TGA → W; AAA → N; AGA/AGG → S | ATG | TAG |
| 15, `BlepharismaMacronuclear` | TAG → Q | ATG | TAA, TGA |
| 16, `ChlorophyceanMitochondrial` | TAG → L | ATG | TAA, TGA |
| 21, `TrematodeMitochondrial` | TGA → W; ATA → M; AAA → N; AGA/AGG → S | ATG, GTG | TAA, TAG |
| 22, `ScenedesmusObliquusMitochondrial` | TCA → stop; TAG → L | ATG | TCA, TAA, TGA |
| 23, `ThraustochytriumMitochondrial` | TTA → stop | ATT, ATG, GTG | TTA, TAA, TAG, TGA |
| 24, `RhabdopleuridaeMitochondrial` | TGA → W; AGA → S; AGG → K | TTG, CTG, ATG, GTG | TAA, TAG |
| 25, `CandidateDivisionSr1Gracilibacteria` | TGA → G | TTG, ATG, GTG | TAA, TAG |
| 26, `PachysolenTannophilusNuclear` | CTG → A | CTG, ATG | TAA, TAG, TGA |
| 27, `KaryorelictNuclear` | TAA/TAG → Q; TGA → W | ATG | TGA (else W) |
| 28, `CondylostomaNuclear` | TAA/TAG → Q; TGA → W | ATG | TAA (else Q), TAG (else Q), TGA (else W) |
| 29, `MesodiniumNuclear` | TAA/TAG → Y | ATG | TGA |
| 30, `PeritrichNuclear` | TAA/TAG → E | ATG | TGA |
| 31, `BlastocrithidiaNuclear` | TAA/TAG → E; TGA → W | ATG | TAA (else E), TAG (else E) |
| 32, `BalanophoraceaePlastid` | TAG → W | TTG, CTG, ATT, ATC, ATA, ATG, GTG | TAA, TGA |
| 33, `CephalodiscidaeMitochondrial` | TAA → Y; TGA → W; AGA → S; AGG → K | TTG, CTG, ATG, GTG | TAG |

In tables 27, 28 and 31, every codon that can end translation also codes an
amino acid, shown in parentheses; none ends translation unconditionally.
`GeneticCode::has_conditional_stops()` is true for exactly these three. Such a
codon translates as its amino acid in every mode, so `ToStop` never stops there;
only a complete CDS may end at one, as its terminal codon. ORF discovery refuses
these tables, because its boundaries need codons that always end translation.

The data are transcribed from NCBI's
[`gc.prt`](https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt), version 4.6,
downloaded on October 22, 2026, which also matches the
[NCBI Genetic Codes page](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
The implementation and its tests use this public biological data, not reference
implementation source; a test parses the same file independently. Identifiers
that are not current NCBI tables, such as 0, 7, 8, 17 to 20 and 34 or more, are
rejected. `GeneticCode::try_from(id)` is the checked numeric entry point; there
is no fallback to table 1. `GeneticCode` is non-exhaustive, so a `match` on it
needs a wildcard arm: NCBI adds tables from time to time.

## Translation modes and ambiguity

| `TranslationMode` | Behavior |
|---|---|
| `Full` | Translate every complete codon; retain definite stops as `*` |
| `ToStop` | Stop before the first definite stop; do not emit its `*` |
| `Cds` | Require a complete CDS, recode the first codon to M and omit the terminal stop |

A complete CDS has at least two codons, a universally accepted start, a last
codon that can end translation in the selected table and no definite internal
stop. All concrete expansions
of an ambiguous first codon must belong to the selected table's start set.
For example, ATN is a valid CDS start in tables 2/11 but not table 1. Initiation
recoding applies only to the first codon; an internal GTG still produces V.

Ambiguous codons are expanded under the selected table. A unique residue is
retained; D/N becomes B, E/Q becomes Z, and I/L becomes J. Other mixed outcomes
become X, including mixed sense/stop codons such as TAN. A stop is definite only
when every concrete expansion is a stop. Therefore X neither ends `ToStop`
translation nor satisfies the terminal-stop requirement of a CDS. A mixed
sense/stop codon inside a CDS is retained as X, not rejected as a definite stop.

Partial codons are checked before scanning for an early stop. `Reject` returns
an error for any trailing one or two bases in the selected suffix, even if an
earlier stop would terminate `ToStop`. `Ignore` permits that suffix remainder
for `Full`/`ToStop` and records its length. CDS always rejects incomplete codons;
setting `Ignore` does not weaken its structural checks.

## Frames and result metadata

`ReadingFrame::{Zero, One, Two}` represent zero-based forward offsets 0, 1 and 2.
`offset()` returns `usize`. `ReadingFrame::try_from(offset)` rejects every other
value; it does not clip, wrap or apply modulo three. Frame selection is not a
strand selection. To translate a reverse complement, construct that sequence
explicitly; offsets then refer to that new sequence, not the original strand.

Given input length `n` and selected frame `f`, the translated suffix is `f..n`.
`f > n` is an error; `f == n` is an empty suffix. Empty suffixes succeed in Full
and ToStop but fail CDS's minimum-length check. Source sequences are immutable.

| `Translation` method | Meaning |
|---|---|
| `as_str()`, `Display` | Peptide text, including any retained `*` markers |
| `len()`, `is_empty()` | Emitted residue/stop count, not nucleotide length |
| `genetic_code()` | Selected `GeneticCode`, even when the product is empty |
| `reading_frame()` | Selected `ReadingFrame` |
| `skipped_leading_bases()` | Exactly `f`, separate from ignored trailing bases |
| `ignored_trailing_bases()` | `(n - f) % 3` for successful translation; does not count complete codons after an early stop |
| `to_protein()` | Validate the product as the existing stop-free protein type; return `SequenceError` on `*` |

Translation equality includes code/frame metadata, not only peptide text.
`to_protein()` preserves ambiguity residues but does not silently remove stops.
It returns a sequence for downstream alignment, not a metadata-bearing
translation result.

## Errors and input validation

`dna(...)`/`rna(...)` first validate the entire source, including the prefix
that will be skipped and the suffix that might be ignored. A `SequenceError`
exposes its original zero-based byte index, byte and alphabet. Gaps, whitespace,
FASTA headers and wrong-molecule T/U are not accepted as sequence symbols.
Parsing files is a separate [reader workflow](file-io-v1.md).

`translate_with` returns `Result<Translation, TranslationError>`. The error enum
is non-exhaustive, so callers matching it should include a fallback arm, as in
the runnable example. Numeric selectors also use this error type.

| Error variant | Meaning and coordinate convention |
|---|---|
| `UnsupportedGeneticCode { id }` | Numeric code is not a current NCBI table |
| `InvalidReadingFrame { frame }` | Numeric frame is not 0, 1 or 2 |
| `FrameOutOfBounds { frame, length }` | Frame exceeds the original input length |
| `PartialCodon { length, remainder }` | Selected suffix length is not divisible by three under the requested policy |
| `CdsTooShort { length }` | Selected suffix has fewer than six bases after partial-codon checks |
| `InvalidStartCodon { codon }` | Historical table-1 start error; preserved unchanged |
| `InvalidStartCodonForCode { codon, genetic_code }` | Invalid start for a table other than 1 |
| `MissingTerminalStop { codon }` | Last codon cannot end translation in the selected table |
| `InternalStop { offset }` | Definite internal CDS stop at a zero-based nucleotide offset in the original source |

Codon error fields retain normalized DNA or RNA letters, so RNA reports U, not
T. For a valid typed input, checks run in this order: frame bounds, partial
codons, CDS minimum length, start, terminal stop, then internal stops. No partial
successful peptide is returned on a translation error. Do not infer a fallback
frame from `FrameOutOfBounds`, such as frame 2 on the one-base input `A`; choose
a biologically intended offset or correct the input. Match on the typed
variants for control flow, not on message text.

## Compatibility, limits and verification

The separate `ConfiguredTranslationExt` trait keeps existing implementations of
`NucleotideSequenceExt` source-compatible. Existing two-field
`TranslationOptions { mode, partial_codons }` literals remain valid.
`TranslationConfig::from(options)` selects table 1/frame 0 while retaining both
legacy policies. Existing convenience methods retain their table-1/frame-0
behavior and error wording through the shared engine.

The supported behavior follows the public
[Biopython 1.88 translation contract](https://biopython.org/docs/1.88/api/Bio.Seq.html#Bio.Seq.translate)
within BioFerris's typed alphabets and the current NCBI tables. This is not the
entire Biopython API: partial codons are rejected by default instead of
truncating with a warning, input is fully validated, frames are explicit, and
stop markers remain `*` rather than being customizable. No gaps, undefined
sequences, arbitrary codon tables, automatic `/transl_table` interpretation,
ORF discovery, six-frame helper, RNA editing, context-dependent recoding,
programmed frameshifts or polyadenylation-based completion of partial stops are
implemented by this API. Selecting table 2 applies its static codon assignments,
not every biological exception discussed on the NCBI page.

Translation uses fixed-table codon expansion with at most 64 concrete outcomes
per codon, giving O(n) time for source length n. The current sequence adapter
copies n validated ASCII bytes and reserves peptide capacity for the selected
suffix, so additional memory is O(n), including in ToStop mode. There is no
streaming peptide iterator or library translation-size budget. File readers
have separate input limits; allocator exhaustion is not a recoverable
translation error. No performance benchmark is claimed here.

The independent Rust checks in
[`tests/genetic_codes.rs`](../../seq-bio/tests/genetic_codes.rs) use grouped
codons, explicit table-2 reassignments, independently listed initiation sets and
string-set ambiguity expansion. Their scope includes all 64 concrete codons
and all 3,375 IUPAC codons in each of three tables, in both DNA and RNA; CDS
start/terminal/internal-stop decisions; B/Z/J/X outcomes; frame boundaries and
metadata; legacy compatibility; generated-input properties; and the actual
Rust example's source and complete stdout. These are specification-derived
checks, not a live differential run or universal compatibility claim.

The existing [phiX174 fixture](../../seq-bio/tests/fixtures/phix174-e.fasta)
also validates its annotated table-11 CDS against the recorded 91-residue
RefSeq protein NP_040709.1. Its provenance is recorded in the
[fixture guide](../../seq-bio/tests/fixtures/README.md#real-sequence). That CDS
starts with ATG and cannot alone establish the broader table-11 initiation set;
the independent exhaustive tests cover those starts separately. Tests use
local data and require neither network access nor a Python runtime.

The result-JSON importer accepts translation results written by the removed CLI
only at table 1/frame 0; see the separately scoped
[result reader](graph-result-input-v1.md#result-json-contract).
