# File I/O and records v1

## Capability index

| Input | Rust | Preservation / supported profile |
|---|---|---|
| Unaligned FASTA | Streaming read/write | Full header, split identity/description, source case |
| Sanger FASTQ | Streaming read/write | Above plus Phred+33 scores; wrapped input |
| Plain/gzip input | Sequential read, by signature | Concatenated members; errors on damaged/truncated streams |
| Aligned FASTA, CLUSTAL, PHYLIP, Stockholm, NEXUS DATA | [Readers](alignment-input-v1.md), [analysis utilities](alignment-analysis-v1.md) | Explicit dialects; gaps/missing cells retained; nontrivial annotated slices reject or report explicit loss |
| GenBank, EMBL | [Readers](annotated-input-v1.md), [same-profile writers/transforms](annotated-output-v1.md), [CDS translation](cds-translation-v1.md) | Conservative annotated records, exact extraction, checked CDS proteins and explicit transformation losses |
| GFA, DOT, existing result JSON | [Imports](graph-result-input-v1.md) | Validated documents, not new algorithm execution |
| Streaming statistics | `summarize_records` | FASTA/Sanger FASTQ; aggregate composition/GC and bounded k-mers |
| Persistent plain indexes | `io::index::SequenceIndex` | Stable plain FASTA/Sanger FASTQ; full verified reopen, then direct lookup |
| Motif-instance FASTA | `FastaReader` with `MotifCounts` | Nonempty equal-width concrete DNA/RNA sites; [matrix/scan guide](motifs-v1.md) |
| BLAST tabular 6/7, legacy XML, single-file XML2 | `io::blast` streaming/collected imports | Explicit profiles, raw values and no-hit queries retained; [full guide](blast-input-v1.md); no search engine |
| Gzip/BGZF output | `io::gzip::GzipWriter`, `io::bgzf::BgzfWriter` | Explicit finalization; create-only named publication; [compressed I/O guide](compressed-io-v1.md) |
| BGZF byte access | `io::bgzf::BgzfReader`, `BgzfVirtualOffset` | Checked members, CRC/ISIZE, explicit EOF policy and virtual positions; not arbitrary gzip seeking |
| Persistent BGZF indexes | `io::bgzf_index::BgzfSequenceIndex` | Stable BGZF FASTA/Sanger FASTQ; source/record/block verification on open; [full guide](compressed-io-v1.md) |

The motif, alignment-analysis, BLAST and compressed/BGZF interfaces are accepted
for their documented scopes in the
integration validation. The
file CLI (`formats`, `inspect-file`, `convert`, `translate-file` and the other
file operations) and its JSON requests were removed with the `bioferris` CLI
(decision 029).

## Contents

- [Rust-first record workflow](#rust-first-record-workflow)
- [FASTA streaming and singleton reads](#fasta-streaming-and-singleton-reads)
- [FASTQ and conversion rules](#fastq-and-conversion-rules)
- [Compression, streaming and limits](#compression-streaming-and-limits)
- [Persistent indexing](sequence-index-v1.md)
- [Analysis, BLAST and compressed workflows](#analysis-blast-and-compressed-workflows)
- [Errors and output safety](#errors-and-output-safety)
- [Verification and references](#verification-and-references)

[API index](../README.md) · [Quick start](../quick-start.md)

## Rust-first record workflow

Manual use means writing a Rust program that imports `bioferris`. From this
repository, save the complete program below in
`seq-bio/examples/read_records.rs` and run `cargo bio read_records`.
That is a user-created example, not a checked-in example target. Outside this
repository use the [path-dependency setup](../quick-start.md#use-the-rust-library).

Literal Sanger FASTQ input (also stored in
[io-example.fastq](../../seq-bio/tests/fixtures/io-example.fastq)):

<!-- io-fastq-input -->
```text
@read info
aCgT
+
!I+~
```

Complete program; using `Cursor` makes this input self-contained:

```rust
use std::io::Cursor;
use bioferris::prelude::*;
use bioferris::io::{FastqReader, one_record};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = "@read info\naCgT\n+\n!I+~\n";
    let record = one_record(FastqReader::<_, Dna>::new(Cursor::new(input)))?;
    println!("ID: {}", record.id());
    println!("Description: {}", record.description().unwrap_or(""));
    println!("Source: {}", record.source_sequence());
    println!("DNA: {}", record.sequence());
    println!("Phred: {:?}", record.qualities().unwrap_or(&[]));
    assert_eq!(record.sequence().to_string(), "ACGT");
    assert_eq!(record.qualities(), Some([0, 40, 10, 93].as_slice()));
    Ok(())
}
```

Actual stdout:

<!-- io-rust-record-output -->
```text
ID: read
Description: info
Source: aCgT
DNA: ACGT
Phred: [0, 40, 10, 93]
```

`Record<A>` is immutable. `header()` retains the complete header after `>` or
`@`, including spacing; `id()` is the first whitespace-delimited token;
`description()` is the trimmed remaining text. `source_sequence()` retains
letter case; `sequence()` exposes the existing uppercase computational type.
Qualities are decoded integers, not printable characters. New constructors
reject empty identities, header controls other than tab, invalid residues,
quality-length mismatches and scores outside 0–93.

This additive model does not change legacy `SequenceRecord`: its identifier is
still the complete trimmed FASTA header. Existing eager `read_fasta`,
`parse_fasta` and alignment algorithms are unchanged. The new record does not
silently become an alignment record or discard qualities.

## FASTA streaming and singleton reads

Use `bioferris::io::fasta::<Dna>(path)?` or `fastq::<Dna>(path)?` for path input;
iterate with `for record in reader { let record = record?; ... }`. Use `Rna`
or `Protein` explicitly for those molecules. No extension/sequence guessing.
For a supplied `BufRead`, use `FastaReader::<_, Dna>::new(reader)`.

Complete streaming example with literal input:

```rust
use std::io::Cursor;
use bioferris::prelude::*;
use bioferris::io::FastaReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = ">gene description\nATGGCTTAA\n>second\nATGTAA\n";
    for record in FastaReader::<_, Dna>::new(Cursor::new(input)) {
        let record = record?;
        println!("{}: {}", record.id(), record.sequence().translate_cds()?);
    }
    Ok(())
}
```

Actual stdout:

<!-- io-rust-fasta-output -->
```text
gene: MA
second: M
```

FASTA ignores blank lines and ASCII whitespace in sequence lines. Headers must
start with `>` in column one. Empty sequences are allowed; an empty stream yields
zero records (legacy eager parsing still rejects an empty document). Duplicate
IDs remain separate records: no implicit dictionary replacement.
`one_record(reader)` requires exactly one record and consumes its EOF; it fails
on zero/multiple records and propagates encountered parsing/decompression errors.

## FASTQ and conversion rules

Only Sanger/Phred+33 is implemented: ASCII 33–126 maps to Phred 0–93. No encoding
autodetection: a Phred+64 file can contain legal +33 characters but mean different
scores. The caller is responsible for selecting this documented encoding.
Wrapped sequence and quality lines and CRLF are supported. Quality termination
uses sequence length, so `@` and `+` within qualities are not delimiters.
A nonempty `+` header must exactly repeat the original `@` header.
Empty reads, blank lines, sequence whitespace and unequal lengths are rejected.

`write_record_fastq(writer, &record)` requires real quality scores: it never
fabricates scores, even with permission to lose data.
`write_record_fasta(writer, &record, LossPolicy::Reject)` rejects quality loss
before writing; choose `LossPolicy::Allow` deliberately to drop them. Both
writers preserve header and source case, wrap symbols at 80 columns, and use
LF. FASTQ emits a bare `+` separator. Wrapping, line endings and redundant
separator text are not byte-preserved; round trips preserve record content.
Arbitrary writers can fail after partial output.

The existing computational `Protein` type still rejects `*` and gaps. Full-mode
translation may write stop-containing peptide FASTA, which that strict protein
reader cannot reload without an explicit biological policy. CDS mode normally
avoids this issue; do not silently remove stops.

## Compression, streaming and limits

Buffered FASTA/FASTQ and BLAST-tabular line parsing retries transient
`std::io::ErrorKind::Interrupted` without losing an accumulated line prefix or
resetting limits. Other read errors retain their typed context and fuse streaming
iteration. A permanently blocking or endlessly interrupting custom reader still
requires a caller-managed deadline; byte budgets are not I/O timeouts. See the
[verified transport regression](../testing/security.md#interrupted-read-regression).

`open_input(path)` and `decode_input(reader)` inspect the gzip signature, not
the filename extension. They support concatenated gzip members using pinned
`flate2 1.1.9` with its Rust backend. No external decompressor or C backend.
Integrity is verified only as data are consumed; stopping early cannot certify
an unread gzip trailer. Consume the entire iterator, propagating every error,
when full-file integrity matters. Unsupported compression is not auto-converted.

Sequence readers perform no reading in their constructors. FASTA retains the
current record and one look-ahead header; FASTQ retains the current record.
Both produce `Result<Record<A>, BioIoError>` and become permanently exhausted
after their first error. This is bounded-record streaming, not constant memory
independent of the largest record. A record keeps both its source text and its
normalized sequence, one byte per symbol each, so reading peaks at about twice
the largest record; the reader memory record
has the GRCh38 measurement.

`with_limits(reader, ReadLimits { .. })` lets Rust callers override:

| Limit | Default | Applies to |
|---|---|---|
| `max_line_bytes` | 8 MiB | Physical line, including its line ending |
| `max_record_symbols` | 100 million | One sequence |
| `max_total_bytes` | 1 GiB | Decompressed sequence-source bytes |

Existing alignment/annotated helpers eagerly collect documents under a 1 GiB
source cap; graph/result JSON have a 64 MiB cap. BLAST additionally provides
bounded complete-query iterators and collecting helpers, with its own explicit
source/query/metadata budgets. Source limits are **not peak-memory guarantees**:
records, parsed structures, compression headers and decompression buffers have
additional costs. This is not a general sandbox for hostile compressed data.
The separate [index API](sequence-index-v1.md) supports stable, uncompressed
regular FASTA/Sanger FASTQ files with source-ordered ID lookup and fully verified
reopening. These generic streaming gzip readers do not provide seekable indexes.
Separate gzip/BGZF writers, checked BGZF virtual access and verified BGZF indexes
are implemented through the [compressed I/O API](compressed-io-v1.md). Ordinary
gzip cannot use BGZF virtual offsets or a BGZF sequence sidecar.

## Analysis, BLAST and compressed workflows

The guides below show literal inputs, complete Rust programs, their actual
results, errors and limits. No runtime Python interpreter or external
search/alignment process is introduced by these file-analysis workflows.

| Workflow | Rust access |
|---|---|
| IUPAC motif matching and count/PWM/PSSM models | [Motif API](motifs-v1.md), direct instances/counts or declared FASTA sites |
| Alignment coordinates, slicing, statistics and consensus | [Alignment analysis](alignment-analysis-v1.md) over existing checked adapters/documents |
| BLAST results | [BLAST input](blast-input-v1.md), `BlastTabularReader` / `BlastXmlReader` |
| Compressed byte output and BGZF indexing | [Compressed I/O](compressed-io-v1.md), composed `Write` sinks and `BgzfSequenceIndex` |

Declare format/molecule where the operation requires them. Compression is a
container, not a biological format: the gzip and BGZF writers preserve the bytes
written to them without parsing or guessing FASTA/FASTQ, whereas a BGZF sequence
index requires a supported biological profile and validates every block and
record. General sequential gzip decoding does not certify BGZF member structure
or create virtual seek metadata.

Named compressed outputs and sidecars remain create-only. Rust stream sinks
require explicit checked `finish`; dropping them is not successful finalization.
Do not interpret a sidecar as a trusted shortcut: opening both plain and BGZF
indexes revalidates the supplied stable source, and the two index formats are
not interchangeable.

## Errors and output safety

New Rust errors expose `kind()`, stable `code()`, `format()`, `message()`, and
optional one-based `record()`/`line()`. I/O failures preserve a source error.
Codes are `io_error`, `invalid_header`, `invalid_sequence`, `invalid_quality`,
`length_mismatch`, `truncated_input`, `unsupported_format_feature`,
`limit_exceeded`, `invalid_structure`, `data_loss` and `invalid_schema`.
Do not parse human messages for control flow.

Named output is create-only. The functions that publish a named file, such as
`SequenceIndex::save_new`, `write_gzip_new`, `write_bgzf_new` and the annotated
writers' create-only paths, write a unique temporary file in the destination
directory, synchronize it, and publish it via a no-overwrite hard link. Parsing/translation failure does not publish a destination; ordinary
cleanup removes the temporary file. Filesystems without hard-link support return
`output_error`; there is no destructive fallback. Existing files, including an
input used as the output path, are never overwritten. A crash or cleanup I/O
failure may leave a temporary file; this is not a directory-durability guarantee.
Use a trusted destination directory: this is not a defense against a concurrent
filesystem adversary who can rename/replace temporary paths in that directory.
A writer given an arbitrary `Write` sink streams and cannot roll back: a late
error can follow earlier records, so publish to a named file when you need an
all-or-nothing result.

## Verification and references

Rust integration/property tests cover preserved fields, legacy behavior, wrapped
qualities, CRLF, singleton cardinality, fused/contextual failures, byte/symbol
limits, gzip members/CRC/truncation, typed molecules, conversion loss and output
publication. Markdown input/output blocks are checked by tests; the complete
Rust programs are rustdoc-tested.
This validates the established sequence/file subsets, not universal Biopython
equivalence or unmeasured performance. See the historical
increment validation report. The added motif,
alignment-analysis, BLAST and compressed/BGZF guides separately identify their
retained reference/oracle evidence and current final-integration status; old
file-I/O validation is not relabeled as acceptance of these new interfaces.

Public behavior references:
[Biopython 1.88 SeqIO](https://biopython.org/docs/1.88/Tutorial/chapter_seqio.html),
[QualityIO](https://biopython.org/docs/1.88/api/Bio.SeqIO.QualityIO.html),
[NCBI SRA FASTQ guidance](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/),
[RFC 1952](https://datatracker.ietf.org/doc/html/rfc1952),
[flate2 1.1.9 MultiGzDecoder](https://docs.rs/flate2/1.1.9/flate2/bufread/struct.MultiGzDecoder.html).
These informed independently written Rust behavior; no Python runtime or
reference-source copying is involved.
