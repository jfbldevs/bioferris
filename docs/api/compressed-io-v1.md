# Compressed I/O: Rust workflows

## Contents

- [Scope and status](#scope-and-status)
- [Write a Rust analysis](#write-a-rust-analysis)
- [Errors and recovery](#errors-and-recovery)
- [Limits and compatibility](#limits-and-compatibility)
- [Verification](#verification)

[API index](../README.md) · [Gzip contract](../algorithms/compressed-output-v1.md) ·
[BGZF contract](../algorithms/bgzf-v1.md) · [BGZF index contract](../algorithms/bgzf-index-v1.md)

## Scope and status

The safe-Rust cores, public modules and Cargo example are implemented and
verified on Rust 1.85.0. Stage 6 is accepted for the documented stream/index
profiles; see the main validation record.
That record also covers the `compress-file`, `index-bgzf` and `get-bgzf-record`
CLI operations, strict JSON requests and discovery, which were removed with the
`bioferris` CLI (decision 029). The complete program and real Cargo
alias below pass the tests. Earlier isolated/staging runs remain separate
checkpoints, not a package-release or universal native-parity claim.

Rust access is library-first: ordinary gzip and BGZF `Write` sinks compose with
existing biological writers. BGZF additionally has checked byte readers,
virtual seeks and verified FASTA/Sanger FASTQ sequence indexes. Ordinary gzip is not seekable by the BGZF virtual-offset API and cannot use a
BGZF sequence sidecar.

## Write a Rust analysis

The literal biological input is:

```text
@read demo
aCgT
+
!5I~
```

The Sanger quality symbols decode to Phred scores 0, 20, 40 and 93. The separate
eight-byte BGZF stream demonstration uses exactly `ACGTACGT`.

The complete runnable program is
[`seq-bio/examples/compressed_io.rs`](../../seq-bio/examples/compressed_io.rs).
It uses a new unique temporary directory;
all named output is create-only. Its explicit final cleanup removes only the
three files this run created, followed by a nonrecursive empty-directory removal.
No input or preexisting user file is removed. On earlier failure, partial
temporary example files may remain for inspection; publication helpers remove
their own unpublished temporary output.

<!-- compressed-io-rust-program -->
```rust
//! Lossless compressed records, explicit finalization and verified BGZF retrieval.

use std::io::{Cursor, Read, Write};
use std::path::PathBuf;

use bioferris::io::bgzf::{
    BgzfOutputError, BgzfReader, BgzfVirtualOffset, BgzfWriteOptions, BgzfWriter, write_bgzf_new,
};
use bioferris::io::bgzf_index::{BgzfIndexOptions, BgzfSequenceIndex};
use bioferris::io::gzip::{GzipWriteOptions, write_gzip_new};
use bioferris::io::index::IndexedFormat;
use bioferris::io::{FastqReader, one_record, open_input, write_record_fastq};
use bioferris::prelude::*;

fn unique_directory() -> std::io::Result<PathBuf> {
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(std::io::Error::other)?
        .as_nanos();
    for attempt in 0..128 {
        let directory = std::env::temp_dir().join(format!(
            "bioferris-compressed-io-{}-{timestamp}-{attempt}",
            std::process::id()
        ));
        match std::fs::create_dir(&directory) {
            Ok(()) => return Ok(directory),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => return Err(error),
        }
    }
    Err(std::io::Error::other(
        "no unique example directory available",
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let directory = unique_directory()?;
    let gzip_path = directory.join("read.fastq.gz");
    let bgzf_path = directory.join("read.fastq.bgz");
    let index_path = directory.join("read.index.json");
    let record = Record::<Dna>::from_ascii("read demo", "aCgT", Some(vec![0, 20, 40, 93]))?;
    println!("Input case: {}", record.source_sequence());

    write_gzip_new(&gzip_path, GzipWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    })?;
    let decoded = one_record(FastqReader::<_, Dna>::new(open_input(&gzip_path)?))?;
    println!(
        "Gzip record: {} {} {:?}",
        decoded.id(),
        decoded.source_sequence(),
        decoded.qualities().ok_or("missing qualities")?
    );

    write_bgzf_new(&bgzf_path, BgzfWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    })?;
    let index = BgzfSequenceIndex::<Dna>::build(
        &bgzf_path,
        IndexedFormat::FastqSanger,
        BgzfIndexOptions::default(),
    )?;
    println!("BGZF IDs: {}", index.ids().collect::<Vec<_>>().join(","));
    index.save_new(&index_path)?;
    let mut reopened =
        BgzfSequenceIndex::<Dna>::open(&bgzf_path, &index_path, BgzfIndexOptions::default())?;
    let found = reopened.get("read")?.ok_or("missing expected record")?;
    println!(
        "BGZF record: {} {} {:?}",
        found.id(),
        found.source_sequence(),
        found.qualities().ok_or("missing qualities")?
    );
    println!("Normalized: {}", found.sequence());
    println!("Missing ID: {}", reopened.get("absent")?.is_none());

    let refused = write_bgzf_new(&bgzf_path, BgzfWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    });
    println!(
        "Overwrite refused: {}",
        matches!(refused, Err(BgzfOutputError::DestinationExists(_)))
    );
    if let Err(error) = BgzfSequenceIndex::<Dna>::build(
        &gzip_path,
        IndexedFormat::FastqSanger,
        BgzfIndexOptions::default(),
    ) {
        println!("Ordinary gzip as BGZF index source: {}", error.code());
    }

    let mut writer = BgzfWriter::with_options(
        Vec::new(),
        BgzfWriteOptions {
            block_size: 4,
            ..Default::default()
        },
    )?;
    writer.write_all(b"ACGTACGT")?;
    let mut reader = BgzfReader::new(Cursor::new(writer.finish()?));
    reader.seek_virtual(BgzfVirtualOffset::new(0, 2)?)?;
    let mut suffix = String::new();
    reader.read_to_string(&mut suffix)?;
    println!("BGZF suffix: {suffix}");
    println!("BGZF EOF: {:?}", reader.finish()?);
    if let Err(error) = BgzfVirtualOffset::new(0, 65_536) {
        println!("Unrepresentable offset: {}", error.code());
    }

    // Remove only the three files this run created in its own new directory.
    // Never recursively remove a workspace or clean unrelated/user input files.
    drop(reopened);
    drop(index);
    std::fs::remove_file(gzip_path)?;
    std::fs::remove_file(bgzf_path)?;
    std::fs::remove_file(index_path)?;
    std::fs::remove_dir(directory)?;
    Ok(())
}
```

Run from the workspace root:

<!-- compressed-io-rust-command -->
```bash
cargo bio compressed_io
```

Complete stdout, exit 0, empty stderr:

<!-- compressed-io-rust-output -->
```text
Input case: aCgT
Gzip record: read aCgT [0, 20, 40, 93]
BGZF IDs: read
BGZF record: read aCgT [0, 20, 40, 93]
Normalized: ACGT
Missing ID: true
Overwrite refused: true
Ordinary gzip as BGZF index source: invalid_bgzf_header
BGZF suffix: GTACGT
BGZF EOF: Present
Unrepresentable offset: invalid_bgzf_virtual_offset
```

Original sequence case, headers and qualities survive both compression paths.
Computational sequences remain normalized. `finish` is checked explicitly;
dropping a writer performs no hidden finalization. Other compressors need not
produce identical compressed bytes; only the decoded bytes are preserved.

The BGZF index validates all members, including the final empty member, parses
all records and fingerprints the complete compressed source and raw record spans.
Reopening it repeats full verification; a sidecar is never trusted merely because
its offsets point at decodable bytes. A retrieved record keeps its source case,
header and qualities, while `sequence()` is normalized; FASTA records have no
qualities rather than fabricated scores.

BGZF seeking uses a checked 48-bit compressed member address plus a 16-bit
within-block position. The 65,536-byte endpoint must use the next member at
offset zero.

## Errors and recovery

The program above shows each refusal without changing its input:

- A key absent from a valid index returns `Ok(None)` (`Missing ID: true`), which
  is distinct from an index or source that fails verification.
- Named output is create-only: `write_bgzf_new` and `write_gzip_new` refuse an
  existing destination with `DestinationExists`, keep the original file and
  remove their own unpublished temporary output. Choose a new destination.
- Ordinary gzip is valid compressed input for the sequential sequence readers
  but not a BGZF index source; building an index from it fails with
  `invalid_bgzf_header`.
- A within-block position of 65,536 cannot be represented and fails with
  `invalid_bgzf_virtual_offset`.

Errors expose stable `code()` identifiers. Compression and finalization I/O
failures preserve typed causes and abort publication. A corrupted member is
checked fully before any bytes of that member are returned; earlier valid stream
prefixes may already have been read.

## Limits and compatibility

`GzipWriteOptions` and `BgzfWriteOptions` default to compression level 6 and
accept 0..=9. `BgzfWriteOptions::block_size` defaults to 65,280 bytes and accepts
1..=65,536.

For indexes, `BgzfIndexOptions` defaults are compressed bytes 1 GiB, decompressed bytes 1 GiB,
1,000,000 members, 100,000 records, 16 MiB aggregate ID bytes, 64 MiB serialized
sidecar, 8 MiB physical lines and 100 million symbols per record. All are
inclusive, with explicit zero behavior. Limits bound distinct resources, not an
exact peak-memory or execution-time guarantee. See the linked contracts for
every error, field and supported format rule.

Sources must remain stable regular files. Handle/path metadata plus digests are
not a filesystem snapshot, mutation lock or authentication guarantee. No general
BAM parser, ordinary-gzip index, SQLite parity, automatic format detection,
quality loss or compressed-byte identity with another backend is claimed.
First-token IDs are case-sensitive.

## Verification

The checks cover all compression levels, malformed headers/trailers,
CRC/ISIZE, DEFLATE expansion bounds, read/write interruptions, failed sinks,
explicit finish/drop behavior, offsets, embedded empty members, full 65,536-byte
boundaries, exact/zero budgets, sidecar forgery, stale sources and create-only
cleanup. Three native Biopython 1.88 fixtures match seven records, their case,
headers, qualities, full raw-span digests and block tables.
[Provenance](../../seq-bio/tests/fixtures/bgzf-provenance.md) retains versions,
hashes and intentional differences.

The checks in `seq-bio/tests/compression_recipes.rs` read this Markdown's
program and result fences and run the real `cargo bio compressed_io` alias. The
main validation record records
accepted Stage 6 scope, toolkit/optimized checks, complete workspace harness
output and terminal audit, documentation and debug/release response/file pairs.
Earlier isolated/staging/reference results remain separate checkpoints. Final
source/evidence review and the local commit accompany this acceptance. No package
release, universal native fidelity, snapshot safety or performance claim follows
from scoped acceptance.
