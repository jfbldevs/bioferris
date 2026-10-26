# Sequence indexing worked recipes

## Contents

- [Setup and literal input](#setup-and-literal-input)
- [Manual Rust program](#manual-rust-program)
- [Binary sidecar and quick check](#binary-sidecar-and-quick-check)
- [Verification and limits](#verification-and-limits)

[API index](../README.md#api-index) · [Indexing contract](../api/sequence-index-v1.md)
· [Recipe index](README.md)

Manual use means writing a Rust program that imports the library. The earlier
CLI recipes on this page were removed with the `bioferris` CLI
(decision 029). This synthetic recipe verifies a defined
FASTA/Sanger FASTQ subset, not universal Biopython equivalence or a performance
claim.

## Setup and literal input

Run from the repository root with the pinned Rust 1.85 toolchain on `PATH` and the
locked dependencies cached. The Rust example uses the repository-local
`cargo bio` alias.

Start without `rust-reads.index.json` and `rust-reads.bfidx`. The examples create
those sidecars and never overwrite them. For another run, use new destination
names; do not remove an unrelated existing file just to run an example.

The program reads this complete checked-in DNA fixture,
[`index-example.fastq`](../../seq-bio/tests/fixtures/index-example.fastq):

<!-- index-recipe-fastq-input -->
```fastq
@read1 coding
atgGCT
+read1 coding
IIIIII
@read2 ambiguous
Nn
+
@+
```

The source contains two records and 67 bytes, including LF line endings and the
final newline. Keys are the case-sensitive first header tokens, `read1` and
`read2`. This is Sanger Phred+33 FASTQ: the second read's `@+` quality line means
scores 31 and 10, not a new record or a separator. Its `Nn` sequence retains
source case alongside the normalized computational sequence `NN`.

The output below is complete, byte-checked stdout, with exit 0 and empty stderr.

## Manual Rust program

This is the complete source of
[`seq-bio/examples/indexed_records.rs`](../../seq-bio/examples/indexed_records.rs):

<!-- index-recipe-rust-program -->
```rust,no_run
//! Build, persist, reopen and query a stable Sanger FASTQ source from Rust.

use std::io;
use std::path::PathBuf;

use bioferris::io::index::{IndexOptions, IndexedFormat, SequenceIndex};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [source, destination] = paths.as_slice() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "usage: cargo bio indexed_records -- SOURCE.fastq NEW_INDEX.json",
        )
        .into());
    };

    let options = IndexOptions::default();
    let index = SequenceIndex::<Dna>::build(source, IndexedFormat::FastqSanger, options)?;
    index.save_new(destination)?;
    if let Err(error) = index.save_new(destination) {
        println!("Existing destination refused: {}", error.code());
    }
    drop(index);

    let mut index = SequenceIndex::<Dna>::open(source, destination, options)?;
    println!("Records: {}", index.len());
    println!("IDs: {}", index.ids().collect::<Vec<_>>().join(", "));
    if let Some(record) = index.get("read2")? {
        println!("Header: {}", record.header());
        println!("Source sequence: {}", record.source_sequence());
        println!("Normalized sequence: {}", record.sequence());
        println!("Phred scores: {:?}", record.qualities());
    }
    println!("Missing ID: {}", index.get("absent")?.is_none());
    Ok(())
}
```

The program expects two command-line file arguments. Run it against the literal
fixture above, creating a new `rust-reads.index.json`:

<!-- index-recipe-rust-command -->
```bash
cargo bio indexed_records -- seq-bio/tests/fixtures/index-example.fastq rust-reads.index.json
```

Complete stdout, exit 0:

<!-- index-recipe-rust-output -->
```text
Existing destination refused: output_error
Records: 2
IDs: read1, read2
Header: read2 ambiguous
Source sequence: Nn
Normalized sequence: NN
Phred scores: Some([31, 10])
Missing ID: true
```

The program builds, saves, drops the original handle, and reopens the persisted
index before lookup. The saved file is a version-1 BioFerris sequence-index
sidecar, not a biological export or a copy of the reads: it holds the source
fingerprint, source-ordered keys and verified record spans, and the original
FASTQ is unchanged. Retrieval needs both the source and its matching index.

`save_new` is create-only: the second save to the same destination fails with
`output_error`, leaving the first sidecar unchanged and no temporary file
behind. Choose a new path if a separate index is needed; there is no overwrite
switch. `ids()` preserves source order. A missing key returns `Ok(None)`, so
`is_none()` is true; absence is not an I/O error. Keys are case-sensitive and do
not include the description: `read2 ambiguous` and `Read2` are not the stored key
`read2`. A retrieved `Record` keeps the header, source case and qualities next
to the normalized sequence.

## Binary sidecar and quick check

A binary sidecar holds the same entries as the JSON one in less space, and an
index opened from it keeps only the identifiers in memory, reading each record's
location from the sidecar when it is fetched. It also records the source
modification time, so an opt-in quick check can compare it instead of rereading
the whole source. This is the complete source of
[`seq-bio/examples/binary_index.rs`](../../seq-bio/examples/binary_index.rs):

<!-- index-recipe-binary-program -->
```rust,no_run
//! Save JSON and binary sidecars, then reopen the binary one with full and quick
//! source checks.

use std::io;
use std::path::PathBuf;

use bioferris::io::index::{
    IndexOptions, IndexedFormat, SequenceIndex, SidecarFormat, SourceCheck,
};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [source, json, binary] = paths.as_slice() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "usage: cargo bio binary_index -- SOURCE.fastq NEW_INDEX.json NEW_INDEX.bfidx",
        )
        .into());
    };

    let options = IndexOptions::default();
    let index = SequenceIndex::<Dna>::build(source, IndexedFormat::FastqSanger, options)?;
    index.save_new(json)?;
    index.save_new_as(binary, SidecarFormat::Binary)?;
    drop(index);
    println!("JSON sidecar: {} bytes", std::fs::metadata(json)?.len());
    println!("Binary sidecar: {} bytes", std::fs::metadata(binary)?.len());

    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut index = SequenceIndex::<Dna>::open_with(source, binary, options, check)?;
        let record = index.get("read2")?.ok_or("missing read2")?;
        println!(
            "{check:?} check: {} records; read2 {} {:?}",
            index.len(),
            record.source_sequence(),
            record.qualities()
        );
    }
    if let Err(error) = SequenceIndex::<Dna>::open_with(source, json, options, SourceCheck::Quick) {
        println!("Quick check of the JSON sidecar: {}: {error}", error.code());
    }
    Ok(())
}
```

Run it against the same literal fixture, creating two new sidecars:

<!-- index-recipe-binary-command -->
```bash
cargo bio binary_index -- seq-bio/tests/fixtures/index-example.fastq rust-reads.index.json rust-reads.bfidx
```

Complete stdout, exit 0:

<!-- index-recipe-binary-output -->
```text
JSON sidecar: 591 bytes
Binary sidecar: 276 bytes
Full check: 2 records; read2 Nn Some([31, 10])
Quick check: 2 records; read2 Nn Some([31, 10])
Quick check of the JSON sidecar: invalid_index: invalid sequence index: quick source checks need a binary sidecar, which records the source modification time
```

Both sidecars come from one build. The binary one is 276 bytes: an 88-byte
header, two 64-byte record entries, the identifiers `read1` and `read2` each
followed by a line feed, two 8-byte lookup positions and a 32-byte SHA-256
checksum. The [binary layout](../algorithms/binary-sequence-index-v1.md) lists
every field.

`SourceCheck::Full`, which `open` also uses, hashes the whole 67-byte source and
reparses both records before any lookup. `SourceCheck::Quick` instead requires
the source length and modification time to equal the saved ones and verifies up
to 32 records, here both, exactly as a lookup does. Both return the same record as
the JSON sidecar.

A JSON sidecar does not record the modification time, so a quick check refuses it
with `invalid_index`; open it with the full check, or save a binary sidecar from
the opened index. A quick check also cannot see every change: an edit that keeps
the length and keeps or restores the modification time, outside the sampled
records, passes it. Lookups still refuse any record whose bytes changed. See
[what a quick check misses](../api/sequence-index-v1.md#what-a-quick-check-misses)
before using it.

## Verification and limits

[`sequence_index_recipes.rs`](../../seq-bio/tests/sequence_index_recipes.rs)
reads the literal fixture, commands, complete programs and expected outputs
directly from this Markdown, and runs the actual `cargo bio indexed_records` and
`cargo bio binary_index` aliases. The test copies the fixture to a unique
temporary directory. Each alias runs from the repository root so Cargo can
resolve it, with its file arguments rebased to that copy; no source fixture or
user output is written.

Run the recipe checks:

```bash
cargo test -p bioferris --test sequence_index_recipes --offline --locked
```

The supported profile is one stable, uncompressed regular FASTA or Sanger FASTQ
source with explicit DNA/RNA/protein typing. Gzip/BGZF, stdin, multiple sources,
custom key functions, other FASTQ encodings and Biopython `index_db` files are
outside this profile. Duplicate first-token IDs fail with `duplicate_record_id`.
An invalid sidecar fails with `invalid_index`; a changed source fails with
`stale_source` when detected. Neither is a missing-key result.

`IndexOptions` defaults are 100,000 records, 16 MiB of key bytes and a 64 MiB
sidecar; its `read_limits` default to 8 MiB physical lines, 100 million symbols
per record and 1 GiB of total input. Limits are inclusive, with zero permitted.
They bound counts/bytes, not exact peak memory.

Verified opening scans the entire source and rederives its entries, unless the
caller opts into the quick check of a binary sidecar. Only subsequent lookups on
the retained index handle seek to selected records. Keep the
source unchanged: metadata/digest checks are not a filesystem snapshot or a
guarantee against hostile concurrent modification. SHA-256 fingerprints check
observed bytes, not source authenticity. The recipes demonstrate the documented
record/persistence contract, not exact Biopython internals or a speed advantage.
