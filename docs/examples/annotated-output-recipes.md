# Annotated output worked recipes

## Contents

- [Setup and literal inputs](#setup-and-literal-inputs)
- [Complete manual Rust program](#complete-manual-rust-program)
- [Observed transformation output](#observed-transformation-output)
- [Create a named GenBank file](#create-a-named-genbank-file)
- [Refuse an existing destination](#refuse-an-existing-destination)
- [Verification and limits](#verification-and-limits)

[API index](../README.md#api-index) · [Annotated output contract](../api/annotated-output-v1.md)
· [Annotated input contract](../api/annotated-input-v1.md) · [Recipe index](README.md)

These workflows are **Rust-only**: write a program that imports `bioferris`, then
run it through the repository's `cargo bio` alias. They do not add annotated
operations to the optional `bioferris` CLI or its JSON request contract. An agent
can use the Rust API and typed errors/report directly, without parsing this
example's human-readable stdout.

## Setup and literal inputs

Run from the repository root with the pinned Rust 1.85 toolchain on `PATH` and the
locked dependencies cached. `cargo bio annotated_records` builds and runs the
checked-in example. It reads these two synthetic DNA fixtures without changing
them. No network service or external-language runtime is involved.

Complete [`annotated-example.gb`](../../seq-bio/tests/fixtures/annotated-example.gb):

<!-- annotated-output-recipe-genbank-input -->
```genbank
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

Complete [`annotated-example.embl`](../../seq-bio/tests/fixtures/annotated-example.embl):

<!-- annotated-output-recipe-embl-input -->
```embl
ID   DEMO; SV 2; circular; genomic DNA; STD; SYN; 15 BP.
AC   DEMO;
DE   Original coding record.
DR   TEST; original-record.
FH   Key             Location/Qualifiers
FH
FT   source          1..15
FT                   /organism="Synthetic construct"
FT   CDS             4..12
FT                   /gene="demo"
FT                   /codon_start=1
FT                   /transl_table=1
FT                   /translation="MA"
FT   tRNA            4..6
FT                   /anticodon=(pos:4..6,aa:Met)
SQ   Sequence 15 BP; 6 A; 4 C; 2 G; 3 T; 0 other;
     aaaATGGCTT AAccc 15
//
```

Both sources contain `aaaATGGCTTAAccc`, with the same CDS at one-based inclusive
coordinates `4..12`. The EMBL input deliberately includes an unreviewed `DR`
field and a feature with an unreviewed coordinate-bearing `/anticodon`
qualifier. Default transformation policy must not silently retain or discard
them.

## Complete manual Rust program

This is the entire executable
[`seq-bio/examples/annotated_records.rs`](../../seq-bio/examples/annotated_records.rs):

<!-- annotated-output-recipe-program -->
```rust,no_run
//! Checked annotated writing and derivation with explicit identity and loss policy.

use bioferris::io::annotated::{
    AnnotationChange, AnnotationLossPolicy, AnnotationTransformOptions, parse_embl, parse_genbank,
    write_embl, write_genbank, write_genbank_new,
};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = std::fs::read_to_string("seq-bio/tests/fixtures/annotated-example.gb")?;
    let source = parse_genbank::<Dna>(&input)?;
    let source = &source[0];
    let options = AnnotationTransformOptions::new("LOCAL_cds")?.with_date("16-OCT-2026")?;
    let sliced = source.slice(3..12, options)?;
    let cds = sliced.record().features()[1].extract(sliced.record().record().sequence())?;
    println!("CDS extraction: {cds}");
    println!("Protein: {}", cds.translate_cds()?);
    println!(
        "Original unchanged: {} {} symbols",
        source.record().id(),
        source.record().sequence().len()
    );
    let mut bytes = Vec::new();
    write_genbank(&mut bytes, sliced.record())?;
    println!("GenBank slice:");
    print!("{}", String::from_utf8(bytes)?);

    let partial = source.slice(
        4..12,
        AnnotationTransformOptions::new("LOCAL_partial")?.with_date("16-OCT-2026")?,
    );
    match partial {
        Ok(_) => return Err("partial CDS unexpectedly accepted".into()),
        Err(error) => println!("Partial CDS refused: {}", error.code()),
    }

    let input = std::fs::read_to_string("seq-bio/tests/fixtures/annotated-example.embl")?;
    let entries = parse_embl::<Dna>(&input)?;
    let entry = &entries[0];
    match entry.reverse_complement(AnnotationTransformOptions::new("LOCAL_rc")?) {
        Ok(_) => return Err("unreviewed metadata unexpectedly accepted".into()),
        Err(error) => println!("Default metadata policy: {}", error.code()),
    }
    let reversed = entry.reverse_complement(
        AnnotationTransformOptions::new("LOCAL_rc")?
            .with_loss_policy(AnnotationLossPolicy::Discard),
    )?;
    for change in &reversed.report().changes {
        match change {
            AnnotationChange::MetadataChanged {
                key, after: None, ..
            } if key == "DR" => println!("Removed metadata: {key}"),
            AnnotationChange::FeatureRemoved { index, key, .. } => {
                println!("Removed feature: {index} {key}");
            }
            _ => {}
        }
    }
    let cds = reversed.record().features()[1].extract(reversed.record().record().sequence())?;
    println!("Reverse-complement CDS extraction: {cds}");
    let mut bytes = Vec::new();
    write_embl(&mut bytes, reversed.record())?;
    println!("EMBL reverse complement:");
    print!("{}", String::from_utf8(bytes)?);

    let mut arguments = std::env::args_os().skip(1);
    if let Some(path) = arguments.next() {
        if arguments.next().is_some() {
            return Err("usage: cargo bio annotated_records -- [NEW_GENBANK_PATH]".into());
        }
        match write_genbank_new(path, sliced.record()) {
            Ok(()) => println!("Saved new GenBank file."),
            Err(error) => {
                eprintln!("Save refused: {}", error.code());
                std::process::exit(1);
            }
        }
    }
    Ok(())
}
```

The program performs an exact `slice(3..12, ...)`, extracts its complete CDS,
translates that sequence explicitly, and then reverse-complements the EMBL
record under an explicit discard policy. Annotation writing itself does not
translate or infer biological qualifiers.

## Observed transformation output

Run the program with no output path; this invocation creates no files:

<!-- annotated-output-recipe-manual-command -->
```bash
cargo bio annotated_records
```

Complete stdout, exit 0; stderr is empty:

<!-- annotated-output-recipe-manual-output -->
```text
CDS extraction: ATGGCTTAA
Protein: MA
Original unchanged: DEMO.2 15 symbols
GenBank slice:
LOCUS       LOCAL_cds 9 bp DNA linear SYN 16-OCT-2026
DEFINITION  BioFerris slice [3,12) of DEMO.2.
FEATURES             Location/Qualifiers
     source          1..9
                     /organism="Synthetic construct"
     CDS             1..9
                     /gene="demo"
                     /codon_start=1
                     /transl_table=1
                     /translation="MA"
ORIGIN
        1 ATGGCTTAA
//
Partial CDS refused: data_loss
Default metadata policy: data_loss
Removed metadata: DR
Removed feature: 2 tRNA
Reverse-complement CDS extraction: ATGGCTTAA
EMBL reverse complement:
ID   LOCAL_rc; SV 1; circular; genomic DNA; STD; SYN; 15 BP.
DE   BioFerris reverse complement of DEMO.2.
FH   Key             Location/Qualifiers
FH
FT   source          1..15
FT                   /organism="Synthetic construct"
FT   CDS             complement(4..12)
FT                   /gene="demo"
FT                   /codon_start=1
FT                   /transl_table=1
FT                   /translation="MA"
SQ   Sequence 15 BP; 3 A; 2 C; 4 G; 6 T; 0 other;
     gggTTAAGCC ATttt                                                         15
//
```

Rust ranges are zero-based and half-open; the source CDS `4..12` corresponds to
`[3, 12)`, giving nine bases and protein `MA` under the standard complete-CDS
translation mode. The proper slice of a circular source becomes linear. The
whole-record `source` feature is updated to cover the derived sequence; the CDS
is retained in full, not clipped.

The attempted `[4, 12)` slice would lose one CDS base, so it fails with
`data_loss`, including under a discard policy. On reverse complement, the
retained CDS becomes `complement(4..12)` and still extracts `ATGGCTTAA`.
Source case travels with its nucleotide, producing `gggTTAAGCCATttt`.

The default policy rejects the EMBL transformation. Explicit `Discard` permits
the reported removal of `DR` and the entire tRNA feature, not a guessed rewrite
of `/anticodon`. Report feature index 2 is zero-based in the original model.
The program prints selected removal entries; `report().changes` contains the
full typed account of metadata, location and feature changes.

Derived records use caller-selected local identities rather than impersonating
the source database accession/version. The GenBank output ID is `LOCAL_cds`;
modern EMBL's local `SV 1` yields ID `LOCAL_rc.1`. The GenBank date is explicitly
provided, never inferred from the clock. Original records remain unchanged.

## Create a named GenBank file

Start with `annotated-slice.gb` absent in a trusted destination directory. Pass
that new path after `--` to the same Rust program:

<!-- annotated-output-recipe-save-command -->
```bash
cargo bio annotated_records -- annotated-slice.gb
```

Exit 0 and empty stderr. Complete stdout is the full
[transformation output above](#observed-transformation-output), followed by this
exact additional line:

<!-- annotated-output-recipe-save-suffix -->
```text
Saved new GenBank file.
```

The new `annotated-slice.gb` contains this complete document, not the surrounding
example report:

<!-- annotated-output-recipe-saved-file -->
```genbank
LOCUS       LOCAL_cds 9 bp DNA linear SYN 16-OCT-2026
DEFINITION  BioFerris slice [3,12) of DEMO.2.
FEATURES             Location/Qualifiers
     source          1..9
                     /organism="Synthetic construct"
     CDS             1..9
                     /gene="demo"
                     /codon_start=1
                     /transl_table=1
                     /translation="MA"
ORIGIN
        1 ATGGCTTAA
//
```

`write_genbank_new` validates the complete model and rendered document before
publishing a new file. It uses same-directory create-only publication; there is
no overwrite or copy fallback. `write_embl_new` provides the same path-level
guarantee for an EMBL record. Stream writers such as `write_genbank(&mut bytes,
record)` preflight semantic failures but cannot roll back a partial stream after
an I/O error; callers own flushing and output handling.

## Refuse an existing destination

Run the exact same command after `annotated-slice.gb` exists:

<!-- annotated-output-recipe-existing-command -->
```bash
cargo bio annotated_records -- annotated-slice.gb
```

The example exits 1. Stdout is the same full transformation output shown above,
with **no** `Saved new GenBank file.` suffix. Complete stderr is:

<!-- annotated-output-recipe-existing-stderr -->
```text
Save refused: io_error
```

The existing file and both input fixtures remain byte-for-byte unchanged. The
example deliberately returns its own exit status 1; this is not the CLI's
schema-version-1 JSON/exit-2 error contract. A Rust caller receives
`BioIoError` with category `Io` and code `io_error`. Choose a new destination if
another derived copy is needed; do not delete an unrelated file to run a recipe.

## Verification and limits

[`annotated_output_recipes.rs`](../../seq-bio/tests/annotated_output_recipes.rs)
reads the literal inputs, complete program, invocations, expected stdout/stderr
and saved GenBank document directly from this page. It runs the actual checked-in
Cargo alias, verifies the create-only refusal, parses the named output again,
and checks that fixtures remain unchanged. Named-output paths are rebased into
unique test-created temporary directories; tests never use a user's
`annotated-slice.gb`.

```bash
cargo test -p bioferris --test annotated_output_recipes --offline --locked
```

These examples verify a conservative same-format model, not universal INSDC
validation, exact source-byte reproduction or Biopython fidelity. Writers retain
supported metadata/feature/qualifier order and source case; regenerate sequence
rows, lengths and standard SQ counts; and reparse to check semantic preservation
before external output. Raw caches are not the writer's data source. Qualifier
physical newlines remain significant and are not arbitrarily reflowed.

Defaults bound a physical output line to 8 MiB including its newline, a record to
100 million symbols, and the complete document to 1 GiB. Rendering is eager, not
a constant-memory streaming encoder. A single rendered location must fit one
bounded physical line, even when the input location was wrapped.

Fuzzy/remote/other opaque locations may be preserved by an unchanged same-format
writer, but transformations reject their unsupported coordinate semantics.
Cross-format conversion is rejected. SQ may contain length alone or the standard
A/C/G/T/other counts; unknown SQ tails are rejected rather than discarded. RNA
storage remains strict: no implicit T/U substitution, and literal U contributes
to the standard SQ `other` count. Transformation restrictions and reviewed
metadata/qualifier policies are specified in the
[annotated output contract](../api/annotated-output-v1.md).
