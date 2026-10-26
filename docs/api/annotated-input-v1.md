# GenBank and EMBL input, version 1

## Contents

- [Availability and scope](#availability-and-scope)
- [Read annotated records](#read-annotated-records)
- [Stream records from large files](#stream-records-from-large-files)
- [Record and metadata contract](#record-and-metadata-contract)
- [Locations and extraction](#locations-and-extraction)
- [Errors and limits](#errors-and-limits)
- [Verification and references](#verification-and-references)

[API capability index](../README.md) · [Rust API](rust-api-v1.md)

## Availability and scope

The Rust module `bioferris::io::annotated` reads bounded GenBank and modern
seven-field EMBL documents, preserving source metadata and feature qualifiers.
It is an explicit subset, not a general INSDC validator or complete Biopython
`SeqIO` replacement. General feature extraction is a Rust API. Separate
[checked writers/transforms](annotated-output-v1.md) and
[annotation-aware CDS translation](cds-translation-v1.md) reuse these models.

| Input | Rust entry points | Supported profile |
|---|---|---|
| GenBank text | `parse_genbank::<Dna>(text)` | `LOCUS`, optional metadata/features, literal `ORIGIN` sequence, `//` |
| GenBank path | `read_genbank::<Dna>(path)` | Same profile; plain/gzip detected by content |
| EMBL text | `parse_embl::<Dna>(text)` | Modern `ID`, optional metadata/`FT`, literal `SQ` sequence, `//` |
| EMBL path | `read_embl::<Dna>(path)` | Same profile; plain/gzip detected by content |
| Record stream | `genbank::<Dna>(path)`, `embl::<Dna>(path)` or `AnnotatedReader::new(reader, format)` | Same profiles, one record at a time; plain/gzip detected by content |
| Feature location | `Location::parse(text)` and `extract(&sequence)` | Exact points/ranges, ordered joins and reverse complements |

All document readers return `Result<Vec<AnnotatedRecord<A>>, BioIoError>` and
buffer the bounded document. [`AnnotatedReader`](#stream-records-from-large-files)
yields the same records one at a time; neither is a random-access index. The
caller explicitly selects `Dna`, `Rna` or, for GenBank `aa` records, `Protein`.
Empty/blank-only input returns an empty vector or no records.

## Read annotated records

The complete synthetic DNA inputs are:

<!-- annotated-genbank-input -->
```text
LOCUS       DEMO 9 bp DNA linear SYN 15-OCT-2026
DEFINITION  Example coding sequence.
ACCESSION   DEMO
VERSION     DEMO.1
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="demo"
ORIGIN
        1 atggcttaa
//
```

<!-- annotated-embl-input -->
```text
ID   DEMO; SV 1; linear; DNA; STD; SYN; 9 BP.
AC   DEMO;
DE   Example coding sequence.
FH   Key             Location/Qualifiers
FT   CDS             1..9
FT                   /gene="demo"
SQ   Sequence 9 BP;
     atggcttaa 9
//
```

Use these strings directly in an ordinary Rust analysis program. Save the
following complete program as `seq-bio/examples/read_annotations.rs` in your
checkout, then run `cargo bio read_annotations` from the repository root. This
example is supplied in the guide, not installed as an additional example target.
The existing [library-first quick start](../quick-start.md#use-the-rust-library)
explains the Cargo alias.

<!-- annotated-program -->
```rust
use bioferris::{Dna, io::annotated::{Location, parse_embl, parse_genbank}};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const GENBANK: &str = r#"LOCUS       DEMO 9 bp DNA linear SYN 15-OCT-2026
DEFINITION  Example coding sequence.
ACCESSION   DEMO
VERSION     DEMO.1
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="demo"
ORIGIN
        1 atggcttaa
//
"#;
    const EMBL: &str = r#"ID   DEMO; SV 1; linear; DNA; STD; SYN; 9 BP.
AC   DEMO;
DE   Example coding sequence.
FH   Key             Location/Qualifiers
FT   CDS             1..9
FT                   /gene="demo"
SQ   Sequence 9 BP;
     atggcttaa 9
//
"#;
    let genbank = parse_genbank::<Dna>(GENBANK)?;
    let embl = parse_embl::<Dna>(EMBL)?;
    let entry = &genbank[0];
    println!("GenBank: {} {}", entry.record().id(), entry.record().sequence());
    println!("EMBL: {} {}", embl[0].record().id(), embl[0].record().sequence());
    println!("CDS: {}", entry.features()[0].extract(entry.record().sequence())?);
    let reverse = Location::parse("complement(join(1..3,4..6))")?;
    println!("Reverse feature: {}", reverse.extract(entry.record().sequence())?);
    let fuzzy = Location::parse("<1..9")?;
    match fuzzy.extract(entry.record().sequence()) {
        Err(error) => println!("Boundary: {}", error.code()),
        Ok(_) => return Err("a fuzzy location must not be extracted exactly".into()),
    }
    Ok(())
}
```

Complete program stdout, successful exit:

<!-- annotated-output -->
```text
GenBank: DEMO.1 ATGGCTTAA
EMBL: DEMO.1 ATGGCTTAA
CDS: ATGGCTTAA
Reverse feature: AGCCAT
Boundary: unsupported_format_feature
```

Both formats identify the same nine-base sequence. This `extract` call returns
DNA and does not interpret translation qualifiers. Use the separate
[feature translation methods](cds-translation-v1.md) to apply those qualifiers
and obtain a checked protein. The reverse
feature joins `ATG` and `GCT`, then reverse-complements the whole result.
The final line is a handled Rust error.

To use files, save either literal input above unchanged as `demo.gb` or
`demo.embl` and replace the corresponding parse call with
`bioferris::io::annotated::read_genbank::<Dna>("demo.gb")?` or
`read_embl::<Dna>("demo.embl")?`. Results are unchanged. Filenames do not select
the format; callers do. Gzip streams are consumed through EOF so CRC/truncation
failures surface rather than returning a partially accepted document.

## Stream records from large files

`AnnotatedReader` reads one record at a time, so memory holds the current
record, not the whole file. It delimits records exactly as the document readers
do and parses each one with the same code: every record equals the
corresponding element of `parse_genbank`/`parse_embl`, raw source text included.
`genbank::<A>(path)` and `embl::<A>(path)` open a plain or gzip file;
`AnnotatedReader::new(reader, format)` reads any buffered reader. Each item is a
`Result<AnnotatedRecord<A>, BioIoError>`.

This input is the repository file `seq-bio/tests/fixtures/annotated-stream.gb`.
Its second record declares RNA:

<!-- annotated-stream-input -->
```text
LOCUS       DEMO1 9 bp DNA linear SYN 22-OCT-2026
DEFINITION  First record: a coding sequence.
ACCESSION   DEMO1
VERSION     DEMO1.1
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="demo"
ORIGIN
        1 atggcttaa
//
LOCUS       DEMO2 6 bp RNA linear SYN 22-OCT-2026
DEFINITION  Second record: declared as RNA.
ACCESSION   DEMO2
VERSION     DEMO2.1
ORIGIN
        1 augcuu
//
LOCUS       DEMO3 12 bp DNA linear SYN 22-OCT-2026
DEFINITION  Third record: no features.
ACCESSION   DEMO3
VERSION     DEMO3.1
ORIGIN
        1 ggatccgaat tc
//
```

The [`stream_annotations` example](../../seq-bio/examples/stream_annotations.rs)
reads it as DNA and asks the reader to continue past records that fail to parse:

<!-- annotated-stream-program -->
```rust
use bioferris::Dna;
use bioferris::io::annotated::{RecordErrorPolicy, genbank};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The documented three-record file; pass a path to read another one.
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "seq-bio/tests/fixtures/annotated-stream.gb".to_owned());
    let records = genbank::<Dna>(&path)?.with_record_errors(RecordErrorPolicy::Continue);
    for result in records {
        match result {
            Ok(entry) => println!(
                "{}: {} bp, {} feature(s)",
                entry.record().id(),
                entry.record().sequence().len(),
                entry.features().len()
            ),
            Err(error) => println!(
                "Skipped record {} at line {}: {}",
                error.record().unwrap_or_default(),
                error.line().unwrap_or_default(),
                error.code()
            ),
        }
    }
    Ok(())
}
```

Run it from the repository root:

<!-- annotated-stream-command -->
```bash
cargo bio stream_annotations
```

Complete stdout, exit status 0:

<!-- annotated-stream-output -->
```text
DEMO1.1: 9 bp, 1 feature(s)
Skipped record 2 at line 11: invalid_sequence
DEMO3.1: 12 bp, 0 feature(s)
```

Record and line numbers are one-based positions in the input. The DNA reader
rejects DEMO2 at its `LOCUS` line because it declares RNA; it never converts the
molecule. With the default `RecordErrorPolicy::Stop`, that error is the last
item and DEMO3 is not read. `Continue` moves on only after an error inside a
complete record, from its `LOCUS`/`ID` line through `//`. Errors that leave the
next record boundary unknown end iteration under either policy: text outside a
record, a missing `//`, an exceeded line or total byte limit, invalid UTF-8, and
I/O or gzip failures. A second pass with `genbank::<Rna>` reads DEMO2 and skips
the two DNA records the same way. Gzip integrity is checked when the stream
reaches its end, so a corrupt trailer arrives as a final `io_error` item after
the last record; read to the end before trusting a compressed file.

`AnnotatedReader::with_limits(io::open_input(path)?, AnnotatedFormat::GenBank,
limits)` accepts a custom `ReadLimits`. `max_line_bytes` and `max_total_bytes`
bound the input as it is read; the total counts every byte of the stream, so
raise it from its 1 GiB default to read larger files. `max_record_symbols`
bounds each record's declared and observed sequence length, and a record over it
is a record error that `Continue` can skip. The reader keeps the current
record's source text and parsed values, and reuses one text buffer sized by the
largest record so far.

## Record and metadata contract

`AnnotatedRecord<A>` exposes:

- `record()`: validated uppercase sequence, original sequence case, derived
  header/id/description; no quality values are invented.
- `locus()` and `declared_molecule()`: original identifier and molecule text.
- `accessions()`: ordered accessions, including repetitions.
- `metadata()`: ordered physical header fields with key, trimmed payload,
  continuation flag and original physical line. Unknown fields remain present.
- `features()`: ordered keys, locations, repeated qualifiers and original lines.
- `raw_record()`: original text from the record header through `//`, including
  its original LF/CRLF endings. Blank lines between records are excluded.

The derived identity prefers version, then first accession, then locus. Modern
EMBL derives its version from `ID`'s `SV` number. Conflicting explicit versions
are errors. GenBank `DEFINITION` continuations and EMBL `DE` lines form the
description, joined with spaces. The derived record header is not a literal
GenBank/EMBL header; exact source text is retained separately.

Qualifiers keep ordering and repeated keys. A flag such as `/pseudo` has
`value: None`; quoted strings decode doubled quotes. Wrapped payloads retain
newlines in `value`; `raw` retains slash, quotes and wrapped payloads, while
`Feature::raw_lines` retains original column indentation. The parser does not
flatten translation strings or apply controlled-vocabulary semantics.
Feature/qualifier names follow the INSDC character grammar and 15/20-character
limits, respectively; unknown names within that grammar remain available.

GenBank `bp` molecule tokens include DNA, RNA, mRNA, tRNA, rRNA, uRNA, snRNA,
snoRNA, scRNA, pre-RNA, ncRNA, tmRNA and cRNA, case-insensitively, optionally
prefixed by `ss-`, `ds-` or `ms-`. Modern EMBL also accepts `genomic DNA`,
`genomic RNA` and the exact case-insensitive token `other DNA`. GenBank `aa`
selects the protein alphabet; EMBL protein flat files are not supported. Unknown
molecule tokens are rejected.
`other DNA` denotes synthetic DNA under the INSDC vocabulary and uses the DNA
alphabet without changing stored residues or the original molecule text. This
narrow addition is covered by the public pUC19 M77789.2 fixture; it does not add
`other RNA`, `unassigned DNA` or arbitrary strings containing `DNA`.

Important RNA limitation: declared molecule and requested alphabet must agree,
and sequence symbols are validated without T/U substitution. INSDC RNA records
using T in their stored sequence are rejected by `Rna`; selecting `Dna` does
not bypass the RNA declaration. This is a known compatibility gap, not a claim
that those source records are invalid. No automatic transcription is performed.

## Locations and extraction

File coordinates are one-based inclusive. `1..9` becomes
`LocationExpression::Range(0..9)`, a Rust zero-based half-open range; point
`9` becomes `8..9`. Exact local ranges are checked against sequence length
during document parsing and again during extraction.

`join(a,b,...)` concatenates operands in their written order, including
overlaps/repetitions. `complement(x)` reverse-complements all of `x`, so a
complemented join reverses segment order as well as orientation. DNA and RNA
use their own complement alphabets; a protein complement is unsupported.
Circular origin-spanning features work only when explicitly represented as
in-bounds joins; topology metadata does not trigger implicit wrapping.

Fuzzy boundaries, remote references, between-base sites, `order`, uncertain
positions and unknown expressions are preserved with an unsupported reason.
They are not silently approximated or certified as valid INSDC syntax.
Extraction returns `unsupported_format_feature`; no remote/network resolver is
invoked. Exact operators with malformed arguments, invalid/overflowing
coordinates or unbalanced parentheses are rejected.

No CDS translation is implicit. `/codon_start`, `/transl_table`,
`/transl_except`, and `/translation` remain annotations. Callers must decide
how to handle them; the current standard-code translation API is not a general
annotation-aware translation engine.

## Errors and limits

Use `BioIoError::kind()` or `code()` for control flow, not message parsing.
Record parsing attaches one-based record/line context. A failure found while
finalizing a record may identify its first line; it is not always the exact
offending field. Standalone location errors use format `insdc-location`;
document parse errors use `genbank` or `embl`.

| Code | Typical cause |
|---|---|
| `invalid_header` | Malformed identifier, length/unit, conflicting versions |
| `invalid_sequence` | Wrong molecule/type, invalid symbol, missing position counter |
| `length_mismatch` | Sequence/declared length, counter or feature bounds disagree |
| `truncated_input` | Missing `//`, interrupted record, unterminated quoted qualifier |
| `invalid_structure` | Malformed supported feature/location/qualifier syntax |
| `unsupported_format_feature` | Legacy EMBL header, unknown molecule, CONTIG-only record, unsupported extraction |
| `limit_exceeded` | Document, line, declared/observed sequence or location nesting bound |
| `io_error` | File/stream or gzip integrity failure |

Defaults are 1 GiB of decompressed document bytes, 8 MiB per physical line,
100 million sequence symbols per record/extracted feature, and 64 levels of
location nesting. For the document readers, document bytes are bounded while
reading and line/record bounds are checked after document buffering; they do not
accept custom limits. `AnnotatedReader` checks line and total bytes as it reads,
and its [`with_limits`](#stream-records-from-large-files) constructor accepts
custom limits. Returned records duplicate source text and parsed values, so
these are input bounds, not a maximum process-memory guarantee.

Literal ORIGIN/SQ sequence and valid position counters are required. Metadata,
features and qualifier syntax follow the documented fixed-column profile;
qualifiers must start on their own continuation line. This is not an import
path for assembly-only CONTIG records, legacy EMBL variants, Swiss-Prot, DDBJ
dialects, arbitrary INSDC extensions or unknown-length sequences. Physical
sequence whitespace and counters are not part of the biological sequence.
EMBL SQ nucleotide-count summaries are retained but not independently validated.
Only its total length and sequence-row counters are checked.

Same-profile [annotated writers and checked transformations](annotated-output-v1.md)
are a separate bounded API; there is no implicit conversion to a metadata-free record.
Borrowing `record()` explicitly accesses only the sequence-level view; callers
must not claim that writing that view as FASTA preserves feature annotations.
Raw source retention alone is not editable, semantically validated round-trip writing.
The separate writer renders and validates the supported model; transformed
records regenerate their raw caches. Cross-format annotation conversion remains
unsupported.

## Verification and references

The Rust integration suite `seq-bio/tests/annotated_import.rs` checks independently
authored fixtures, metadata and repeated qualifiers, all exact subranges of an
IUPAC DNA fixture against an independent complement table, joined/strand-aware
extraction, malformed/truncated inputs, limits, strict RNA behavior, plain/gzip
paths and corrupt gzip CRCs. It reads the literal inputs/output above and compares
them to computed results. The complete Rust program is also executable as a
Markdown doctest; no reference runtime or alignment qualification is involved.

From the repository root on the pinned Rust 1.85 toolchain:

```bash
cargo test -p bioferris --test annotated_import --offline --locked
cargo test -p bioferris --test annotated_streaming --offline --locked
cargo test -p bioferris --doc --offline --locked
cargo build -p bioferris --lib --offline --locked
rustdoc --test docs/api/annotated-input-v1.md --edition 2024 --extern bioferris=target/debug/libbioferris.rlib -L dependency=target/debug/deps -D warnings
```

Scoped verification on 2026-10-15 passed: 17 annotated integration tests, the
module doctest and the complete Markdown program doctest. The program was also
compiled and executed separately; its stdout is reproduced above. Strict Clippy
for the annotated integration target passed. These checks do not replace the
whole-workspace integration checks recorded by the parent file-I/O increment.

The streaming reader was added on 2026-10-22. `seq-bio/tests/annotated_streaming.rs`
checks that streamed GenBank and EMBL records equal the document readers' records,
including CRLF, blank lines and a missing final newline; that framing errors
carry the same kind, record and line as the document readers' and end iteration
under both policies; the `Stop` and `Continue` policies, including a second pass
with the other molecule; line, total and record-symbol limits; plain and gzip
paths, a late gzip CRC failure and invalid UTF-8. It also checks that the input,
program and command above are the actual fixture, example and alias, and runs
`cargo bio stream_annotations` to compare its stdout with the output shown.

This is specification-based subset verification, not a live Biopython differential
campaign, performance measurement, or universal INSDC compatibility claim.
Public format references consulted on 2026-10-15:

- [INSDC Feature Table Definition 11.4, April 2026](https://www.insdc.org/submitting-standards/feature-table/)
- [NCBI GenBank sample record and field descriptions](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
- [ENA flat-file example](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html)
- [ENA user manual, modern ID line molecule field](https://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/doc/usrman.txt)

The worked inputs above and the initial boundary fixtures are synthetic.
Additional [versioned public-data fixtures and provenance](../../seq-bio/tests/fixtures/annotated-provenance.md)
cover the real DQ067560.1 genomic-DNA record from both NCBI and ENA: 534 bases,
source-declared nucleotide counts and exact exon 6 extraction. Its fuzzy CDS
is preserved and explicitly not extracted. The M77789.2 pUC19 fixture verifies
`other DNA` import, its 2,686 bases and source-declared base counts. X56734.1
records the unchanged RNA/T-storage limitation.
All subsequent fixture tests are offline; public downloads are not part of CI.
