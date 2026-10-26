# Annotation-aware CDS translation

## Contents

- [Availability and policies](#availability-and-policies)
- [Run a Rust analysis](#run-a-rust-analysis)
- [Errors, limits and compatibility](#errors-limits-and-compatibility)

[API capability index](../README.md) · [Normative contract](../algorithms/cds-translation-v1.md)

## Availability and policies

Translate an exact local `CDS` feature from an existing GenBank/EMBL record into
a checked protein. Rust methods are `Feature::translate_cds(parent)` and
`translate_cds_with(parent, CdsTranslationOptions)`. The `translate-cds-features`
CLI/JSON operation was removed with the `bioferris` CLI
(decision 029).

`codon_start` is a one-based annotation value (1/2/3), applied after extracting
joins and reverse complements into coding orientation. `transl_table` selects
any current [NCBI table](genetic-codes-v1.md#supported-genetic-codes) and takes
precedence over the configured fallback in `CdsTranslationOptions`. A complete
CDS needs a supported initiator, complete codons, a last codon that can end
translation and no internal definite stop. The returned peptide begins with M and
excludes the terminal stop.

A CDS with the flag `/ribosomal_slippage` is translated from its annotated
location exactly as written. Records usually encode the frameshift there, as a
join whose parts repeat or skip bases, for example
`join(13469..13546,13546..21555)`; some flag a CDS whose own location has no
shift. The flag must have no value and appear once, otherwise the CDS is refused
as `invalid_cds_qualifier`. `CdsTranslation::ribosomal_slippage()` reports the
flag. The complete-CDS rules and the reference check apply as for any CDS, so a
product read without the shift it needs fails them.

| Reference policy | Behavior |
|---|---|
| `CdsReferencePolicy::CheckIfPresent` (default) | Check exact `/translation` equality if present; otherwise return `Absent` |
| `CdsReferencePolicy::RequireMatch` | Require a present and exactly matching reference |
| `CdsReferencePolicy::Ignore` | Do not interpret the reference; return `Ignored` |

Checking removes only ASCII whitespace from a wrapped reference. Case, residues,
stops and differences are not corrected. A successful check returns `Matched`.
Explicit `Ignore` disables only this comparison, not biological/extraction checks.

## Run a Rust analysis

Work in the repository root, using the [repository setup](../quick-start.md#enter-the-repository).
The checked-in [complete example](../../seq-bio/examples/translate_annotations.rs)
contains this literal GenBank input, also saved as
[`annotated-example.gb`](../../seq-bio/tests/fixtures/annotated-example.gb):

<!-- cds-input -->
```text
LOCUS       DEMO 15 bp DNA circular SYN 15-OCT-2026
DEFINITION  Original coding record.
ACCESSION   DEMO
VERSION     DEMO.2
FEATURES             Location/Qualifiers
     source          1..15
                     /organism="Synthetic construct"
     CDS             4..12
                     /gene="demo"
                     /codon_start=1
                     /transl_table=1
                     /translation="MA"
ORIGIN
        1 aaaATGGCTT AAccc
//
```

Complete program, with an intentionally mismatching cloned reference to demonstrate
error handling and explicit recovery:

<!-- cds-rust-program -->
```rust
//! Translate an annotated CDS and check its source protein with explicit policies.

use bioferris::{io::annotated::parse_genbank, prelude::*};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const INPUT: &str = r#"LOCUS       DEMO 15 bp DNA circular SYN 15-OCT-2026
DEFINITION  Original coding record.
ACCESSION   DEMO
VERSION     DEMO.2
FEATURES             Location/Qualifiers
     source          1..15
                     /organism="Synthetic construct"
     CDS             4..12
                     /gene="demo"
                     /codon_start=1
                     /transl_table=1
                     /translation="MA"
ORIGIN
        1 aaaATGGCTT AAccc
//
"#;
    let records = parse_genbank::<Dna>(INPUT)?;
    let record = records.first().ok_or("the example requires one record")?;
    let feature = record
        .features()
        .iter()
        .find(|feature| feature.key == "CDS")
        .ok_or("the example requires a CDS")?;
    let product = feature.translate_cds(record.record().sequence())?;
    println!("Record: {}", record.record().id());
    println!("Location: {}", feature.location.raw());
    println!(
        "Extracted DNA: {}",
        feature.extract(record.record().sequence())?
    );
    println!("Genetic code: {}", product.genetic_code().id());
    println!("Codon start: {}", product.codon_start());
    println!("Protein: {}", product.peptide());
    println!("Reference: {:?}", product.reference_status());

    let mut mismatch = feature.clone();
    let reference = mismatch
        .qualifiers
        .iter_mut()
        .find(|qualifier| qualifier.name == "translation")
        .ok_or("the example requires a reference translation")?;
    reference.value = Some("MV".to_owned());
    match mismatch.translate_cds(record.record().sequence()) {
        Err(error) => println!("Mismatch: {error}"),
        Ok(_) => return Err("a mismatching reference must fail".into()),
    }
    let recovered = mismatch.translate_cds_with(
        record.record().sequence(),
        CdsTranslationOptions {
            reference_policy: CdsReferencePolicy::Ignore,
            ..CdsTranslationOptions::default()
        },
    )?;
    println!(
        "Explicitly ignored reference: {:?}",
        recovered.reference_status()
    );
    println!("Computed protein: {}", recovered.peptide());
    Ok(())
}
```

Run the ordinary Rust program:

<!-- cds-rust-command -->
```bash
cargo bio translate_annotations
```

Actual complete stdout; exit 0, stderr empty:

<!-- cds-rust-output -->
```text
Record: DEMO.2
Location: 4..12
Extracted DNA: ATGGCTTAA
Genetic code: 1
Codon start: 1
Protein: MA
Reference: Matched
Mismatch: CDS /translation differs at residue 1: reference MV, computed MA
Explicitly ignored reference: Ignored
Computed protein: MA
```

File coordinates `4..12` are one-based inclusive. Extraction produces nine
bases in coding orientation; the two-amino-acid protein excludes TAA. The original
record stays unchanged when the cloned qualifier is modified. The mismatch is a
handled Rust error; its position is a zero-based protein offset (residue 1 is
the second residue, A computed against V in the reference). `Ignore` makes the
unchecked reference explicit; unsupported exception qualifiers and invalid
complete CDS products still fail under `Ignore`.
`product.translation().to_protein()?` provides the existing typed protein for
alignment and other sequence workflows.

To translate every CDS of a file, read it with `read_genbank` or `read_embl`
(plain or gzip input, detected by signature) and call `translate_cds` on each
feature whose key is `CDS`, in record and feature order. Any translation-error
nucleotide offset is relative to the extracted coding-oriented sequence, not
genomic coordinates.

## Errors, limits and compatibility

| Error code | Meaning |
|---|---|
| `invalid_cds_feature` | Rust method called on a non-`CDS` feature |
| `unsupported_cds_qualifier` | Pseudogene or biological translation exception requires unsupported semantics |
| `invalid_cds_qualifier` | Duplicate, missing-value or invalid numeric/reference qualifier |
| `unsupported_genetic_code` | Requested/annotated table is not a current NCBI table |
| `missing_cds_translation` | `require-match` requires a `/translation` qualifier |
| `cds_translation_mismatch` | Reference and computed protein differ in residues or length |
| `cds_translation_error` | Complete-CDS checks failed; the `Translation` variant holds the original translation error |
| `unsupported_format_feature`, `length_mismatch` | Extraction is unsupported or out of bounds |
| `limit_exceeded` | Existing parser/extraction limit exceeded |

`CdsTranslationError::code()` returns these identifiers; inspect the typed
variants for structured details. Parsing errors retain existing format/line
context.

The current conservative GenBank/EMBL input profiles and their eager 1 GiB
decompressed-document cap remain applicable. Input size is not a peak-memory
bound. Exact local points/ranges/joins/complements are supported; fuzzy/remote
locations remain unsupported. Explicitly joined origin-crossing features work.
Pseudogenes, `transl_except`, `exception` and `trans_splicing` are rejected; no
exception-based recoding or partial-CDS inference is performed.
Duplicate interpreted qualifiers are refused, even when values agree. When
reference comparison is ignored, `/translation` is not interpreted at all.

Biopython 1.88's documented feature translation informs extraction/table/frame
semantics. Automatic exact reference comparison and explicit conservative refusals
are BioFerris policies. Tests use independent biological expectations and all 11
published proteins in an unchanged phiX174 GenBank snapshot. The products also
match an independently captured Biopython 1.88/CPython 3.12.3 behavioral snapshot
for those 11 CDS features. Rust tests consume that snapshot offline; no interpreter
is needed. See [fixture provenance](../../seq-bio/tests/fixtures/cds-provenance.md)
and the validation record. This targeted
agreement does not establish universal annotation compatibility.
