# Alignment input, version 1

## Contents

- [Availability and scope](#availability-and-scope)
- [Read aligned FASTA](#read-aligned-fasta)
- [Select a format explicitly](#select-a-format-explicitly)
- [Preserve annotations and missing data](#preserve-annotations-and-missing-data)
- [Supported profiles](#supported-profiles)
- [Errors and resource limits](#errors-and-resource-limits)
- [Verification and references](#verification-and-references)

[API capability index](../README.md) · [Rust quick start](../quick-start.md)

## Availability and scope

`bioferris::io::alignment` reads aligned FASTA, CLUSTAL, four explicit PHYLIP
profiles, Stockholm 1.0, and a restricted NEXUS DATA block through the Rust
library. An AI agent calls the same typed Rust API without parsing human output.

An `AlignmentDocument<Dna>`, `AlignmentDocument<Rna>` or
`AlignmentDocument<Protein>` retains rectangular rows, input row order, named
annotations, explicit missing cells and all-gap columns. Alphabet selection is
mandatory. Valid residues become uppercase; source wrapping, whitespace, case,
and general format headers are not a byte-preservation contract. Comments and
supported annotation text are retained, not executed. Annotations have a scope,
optional sequence name, tag, value and one-based source line.

`parse_alignment(reader, format)` accepts a buffered reader. `read_alignment(path,
format)` opens plain or gzip input through the shared signature-detecting input
layer. Both read one complete document into memory; neither is a streaming
iterator. They do not infer molecule or format from a filename.

For statistics, consensus, checked slices and source-coordinate maps, the
separate [alignment-analysis workflow](alignment-analysis-v1.md) reuses these
explicit reader profiles; it is not a generic alignment-conversion tool.

## Read aligned FASTA

<!-- alignment-fasta-recipe -->
Literal DNA alignment input:

```fasta
>reference
AC-G
>sample
ACTG
```

Complete Rust program, using the library directly:

```rust
use bioferris::{Dna, io::alignment::{AlignmentFormat, parse_alignment}};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = b">reference\nAC-G\n>sample\nACTG\n";
    let document = parse_alignment::<Dna>(&input[..], AlignmentFormat::Fasta)?;
    println!("{} rows, {} columns", document.rows().len(), document.width());
    for row in document.rows() {
        println!("{}: {}", row.identifier(), row.aligned());
    }
    let alignment = document.to_multiple_alignment()?;
    assert_eq!(alignment.rows()[0].to_ascii_with_gaps(), "AC-G");
    Ok(())
}
```

Observed stdout, exit 0:

```text
2 rows, 4 columns
reference: AC-G
sample: ACTG
```

The third column is a gap in the reference, not a deletion of the column itself.
The checked adapter produces the existing `MultipleAlignment<Dna>` without
rerunning an alignment algorithm. FASTA row identifiers retain the complete
trimmed header, unlike the new sequence-record API's split identifier.
<!-- /alignment-fasta-recipe -->

To compile and execute all complete Rust blocks in this page from the toolkit
root (the pinned toolchain and offline dependencies must already be available):

```bash
cargo build -p bioferris --lib --offline --locked
rustdoc --test docs/api/alignment-input-v1.md --edition 2024 --extern bioferris=target/debug/libbioferris.rlib -L dependency=target/debug/deps -D warnings
```

To run the integration check that prints and verifies the first result:

```bash
cargo test -p bioferris --test alignment_import documented_fasta_input_output --offline --locked -- --exact --nocapture
```

Cargo adds its test-harness report around the biological stdout. For an ordinary
analysis program, use the [library-first setup](../quick-start.md#use-the-rust-library)
and call the same functions.

## Select a format explicitly

These complete, independent inputs each describe the same two-row, four-column
alignment. All PHYLIP dialect decisions are explicit: `Strict` reserves the first
ten ASCII bytes for the name; `Relaxed` uses a whitespace-delimited name.
`Sequential` finishes a row before the next name; `Interleaved` reads equal-width
fragments in repeated row order. Continuation lines have no repeated names.

```rust
use bioferris::{Dna, io::alignment::{AlignmentFormat, parse_alignment}};
use bioferris::prelude::{PhylipLayout, PhylipNameStyle};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inputs = [
        ("CLUSTAL", AlignmentFormat::Clustal,
         "CLUSTAL W\n\na AC-G\nb ACTG\n"),
        ("PHYLIP relaxed sequential", AlignmentFormat::Phylip {
            names: PhylipNameStyle::Relaxed,
            layout: PhylipLayout::Sequential,
         }, "2 4\na AC\n-G\nb AC\nTG\n"),
        ("PHYLIP strict interleaved", AlignmentFormat::Phylip {
            names: PhylipNameStyle::Strict,
            layout: PhylipLayout::Interleaved,
         }, "2 4\na         AC\nb         AC\n\n-G\nTG\n"),
        ("Stockholm", AlignmentFormat::Stockholm,
         "# STOCKHOLM 1.0\na AC-G\nb ACTG\n//\n"),
        ("NEXUS", AlignmentFormat::Nexus,
         "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA;\nMATRIX\na AC-G\nb ACTG\n;\nEND;\n"),
    ];
    for (label, format, input) in inputs {
        let document = parse_alignment::<Dna>(input.as_bytes(), format)?;
        let alignment = document.to_multiple_alignment()?;
        assert_eq!(alignment.rows()[0].to_ascii_with_gaps(), "AC-G");
        assert_eq!(alignment.rows()[1].to_ascii_with_gaps(), "ACTG");
        println!("{label}: {} rows, {} columns", document.rows().len(), document.width());
    }
    Ok(())
}
```

Observed stdout, exit 0:

<!-- alignment-formats-output -->
```text
CLUSTAL: 2 rows, 4 columns
PHYLIP relaxed sequential: 2 rows, 4 columns
PHYLIP strict interleaved: 2 rows, 4 columns
Stockholm: 2 rows, 4 columns
NEXUS: 2 rows, 4 columns
```
<!-- /alignment-formats-output -->

Run this block with the `rustdoc --test` command above. The integration suite
also verifies all four PHYLIP combinations with wrapped independent fixtures,
not just these two examples.

## Preserve annotations and missing data

The computational `MultipleAlignment` type has no annotation or explicit
missing-cell representation. Its default conversion therefore rejects potential
information loss. Permission to omit annotations is a separate, named method;
the borrowed source document stays unchanged. It does not permit changing missing
cells to gaps, replacing them with `N`/`X`, or deleting all-gap columns.

```rust
use bioferris::{Dna, io::alignment::{AlignmentFormat, parse_alignment}};
use bioferris::io::IoErrorKind;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = b"# STOCKHOLM 1.0\n#=GF ID example\na AC.G\nb ACTG\n//\n";
    let document = parse_alignment::<Dna>(&input[..], AlignmentFormat::Stockholm)?;
    let error = document.to_multiple_alignment().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::DataLoss);
    println!("Strict conversion: {}", error.code());
    println!("Retained annotation: {}={}",
             document.annotations()[0].tag(), document.annotations()[0].value());
    let alignment = document.to_multiple_alignment_discarding_annotations()?;
    println!("Explicit conversion: {}", alignment.rows()[0].to_ascii_with_gaps());
    println!("Unchanged source: {}", document.rows()[0].aligned());

    let missing = b">a\nAC?G\n>b\nACTG\n";
    let document = parse_alignment::<Dna>(&missing[..], AlignmentFormat::Fasta)?;
    let error = document.to_multiple_alignment_discarding_annotations().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::DataLoss);
    println!("Missing-cell conversion: {}", error.code());
    Ok(())
}
```

Observed stdout, exit 0 (errors are handled by the program):

<!-- alignment-loss-output -->
```text
Strict conversion: data_loss
Retained annotation: ID=example
Explicit conversion: AC-G
Unchanged source: AC.G
Missing-cell conversion: data_loss
```
<!-- /alignment-loss-output -->

`.` and `-` are both Stockholm gap glyphs; `?` is this profile's explicit
missing-data marker. Converting a dot gap to the kernel's gap cell does not change
its biological meaning. Unsupported semantics are not implicitly repaired.

## Supported profiles

| Input | Accepted subset | Explicit limits / differences |
|---|---|---|
| Aligned FASTA | Nonempty set of unique full-header names; wrapped sequences; blank lines and sequence ASCII whitespace; `-` gaps, `?` missing cells | Equal final widths required; no semicolon comments, dot-match notation, or protein `*`; distinct from ungapped sequence FASTA |
| CLUSTAL | `CLUSTAL` header; blank-separated equal-width blocks in stable name/order; optional cumulative nongap residue counts; `*:.` consensus with spaces | Not `MUSCLE` headers; no changing row order; consensus retained as a column annotation, so default conversion rejects it; missing `?` not supported |
| PHYLIP | Exactly two header dimensions; strict/relaxed names and sequential/interleaved layout selected explicitly; sequence whitespace; `-` gaps, `?` missing cells | No auto-detection, repeated continuation names, numbered residues, dot-match notation, options, trees or multiple datasets; narrower alphabet than historical PHYLIP |
| Stockholm 1.0 | Exact first-line version; final `//`; wrapped rows; `-`/`.` gaps; `?` missing cells; GF/GS/GC/GR and comments | One alignment only; unknown markup prefixes rejected; free-text GF/GS repeats retained separately; GC/GR fragments concatenate by target/tag and must match final width; no annotation-specific biological interpretation |
| NEXUS | Exactly `BEGIN DATA`, `DIMENSIONS`, `FORMAT`, `MATRIX`, `END`/`ENDBLOCK`; explicit DNA/RNA/PROTEIN datatype; one row per line; optional `INTERLEAVE=NO`; comments; quoted names with doubled apostrophes | No TAXA/CHARACTERS/TREES/SETS/other blocks, wrapped/interleaved matrices, `MATCHCHAR`, `EQUATE`, `TOKENS`, `TRANSPOSE`, or polymorphism braces; unknown commands/options rejected |

NEXUS names without quotes translate underscores into spaces; quoted underscores
remain literal. Default NEXUS glyphs are gap `-` and missing `?`; distinct custom
unquoted ASCII punctuation glyphs are accepted except syntax delimiters
`[];'"=`. NEXUS comments are retained as document annotations, including nested
comment text. All formats reject duplicate final names. Zero-width rows are
accepted if the document/profile expresses them; zero rows are not.

No annotation-preserving document writer is provided in this increment. Existing
alignment writers accept `MultipleAlignment`, so converting then writing is a
semantic round trip only for representable data; it cannot recreate annotations,
source layout, or omitted original whitespace. Historical writers and algorithms
are unchanged.

## Errors and resource limits

All functions return `BioIoError`. Branch on `kind()` or stable `code()`, not the
free-text explanation. `format()` identifies the selected input profile.
`record()` and `line()` provide one-based context where available; alignment
record context means a row ordinal, not a separate document. Some whole-document
or decompression errors do not have row/line context. `AlignmentRow::line()` is
the first source line for a row, not the line of every wrapped fragment.

| Code | Typical cause |
|---|---|
| `invalid_header` | Missing magic/header, empty identifier or invalid name field |
| `invalid_sequence` | Byte outside the selected alphabet; NEXUS datatype mismatch |
| `length_mismatch` | Unequal row widths, declared dimensions/counts, or annotation width |
| `invalid_structure` | Duplicate names, malformed block/order, unknown annotation target, or kernel-incompatible all-gap column |
| `truncated_input` | Missing PHYLIP continuation, Stockholm terminator, or NEXUS command terminator/quote/comment closure |
| `unsupported_format_feature` | A rejected dialect, command, markup prefix or trailing dataset |
| `data_loss` | Implicit annotation omission or unrepresentable missing cells |
| `io_error` | Filesystem, reader or gzip integrity failure |
| `limit_exceeded` | Document exceeds the common 1 GiB decompressed-input cap, or an integer count overflows |

These eager readers use the shared document byte cap, not the streaming readers'
per-line/per-record controls. They can retain multiple text/row copies, so the
byte cap is not a peak-memory guarantee. There is no random-access index or
whole-collection iterator. Gzip integrity is checked because these readers
consume the full input; the shared decoder accepts concatenated members.

## Verification and references

`seq-bio/tests/alignment_import.rs` checks independent inputs, four PHYLIP
profiles, case normalization, exact dimensions, missing cells versus gaps,
all-gap columns, annotation scope/content/width, malformed input, truncation,
source positions, injected reader errors, plain/gzip paths and gzip checksum
failure. A 128-case property test feeds
arbitrary UTF-8 both directly and inside each profile's header/framing to check
for panics. It also checks semantic readback
of the existing FASTA, CLUSTAL, PHYLIP, Stockholm and NEXUS writers. Documentation
tests verify literal inputs and the full printed outputs above. Run:

```bash
cargo test -p bioferris --test alignment_import --offline --locked
```

This is supported-subset verification, not full Biopython compatibility, a live
differential comparison, performance evidence, or renewed scientific algorithm
qualification. Design references consulted 2026-10-15:

- [PHYLIP 3.69 sequence format documentation](https://phylipweb.github.io/phylip/doc/sequence.html): sequential/interleaved structure and strict names. This implementation deliberately rejects legacy numbered symbols and unsupported alphabet extensions.
- [Stockholm 1.0 specification](https://sonnhammer.sbc.su.se/Stockholm.html): gap glyphs and four annotation scopes. The reader's typed alphabet and missing-cell policy are explicit restrictions.
- [Paul Lewis's NEXUS format guide](https://plewis.github.io/nexus/): DATA block, quoted names, comments and format options. Only the above restricted block is implemented.
