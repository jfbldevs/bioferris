# BGZF byte streams and virtual offsets

## Contents

- [Scope and status](#scope-and-status)
- [References](#references)
- [Reader contract](#reader-contract)
- [Offsets and seeking](#offsets-and-seeking)
- [Writer contract](#writer-contract)
- [Worked Rust example](#worked-rust-example)
- [Errors and limits](#errors-and-limits)
- [Acceptance](#acceptance)

[API index](../README.md) · Implementation plan

## Scope and status

The Stage 6 BGZF core provides public `io::bgzf` exports for safe Rust byte
compression, decompression and virtual seeks. Gzip output, create-only
publication and indexed FASTA/FASTQ retrieval are also implemented and verified;
the CLI/strict-JSON adapters were removed with the `bioferris` CLI
(decision 029). Stage 6 is accepted for the
documented stream/index profiles in the
main validation record.
Prior isolated/staging results remain separate checkpoints.
A BGZF byte reader is not
a BAM parser, sequence parser, authenticated index or general gzip reader.

## References

- [SAM/BAM specification 1.6, revision b5341fb, 2025-08-12](https://samtools.github.io/hts-specs/SAMv1.pdf), sections 4.1–4.1.2; the block-layout table and EOF discussion on printed pages 13–14 were rendered and inspected on 2026-10-17.
- [RFC 1952](https://www.rfc-editor.org/rfc/rfc1952): gzip member framing and CRC/size trailer.
- [Biopython 1.88 Bio.bgzf public API](https://biopython.org/docs/1.88/api/Bio.bgzf.html): binary streams and virtual-position interoperability.
- [HTSlib bgzip manual](https://www.htslib.org/doc/bgzip.html), version 1.24: BGZF interoperability and the distinction between compressed offsets and index entries.

BGZF concatenates gzip members, each bounded to 65,536 compressed bytes including
framing and 65,536 uncompressed bytes. Each member carries one `BC` extra subfield
with a two-byte little-endian payload containing its compressed size minus one.
Virtual offsets combine a 48-bit compressed block address and a 16-bit offset
within its decompressed data. They are positions, not linear byte counts.
The recommended terminal marker is a particular 28-byte empty member. Embedded
empty members do not end a stream; later members must still be read.

## Reader contract

`BgzfReader<R: Read>` decodes incrementally with at most one validated data block
retained. `Read` and `BufRead` adapters expose the same computation as typed
`read_bytes`/`finish` methods. A block is checked completely before its bytes are
published. The caller can receive a valid earlier prefix before a later error.

This profile requires compression method 8 and gzip flags exactly 4. Other gzip
optional-header forms are refused. Additional well-formed extra subfields are
accepted; missing/duplicate `BC`, wrong `SLEN`, truncated subfields and inconsistent
block lengths fail. Raw DEFLATE must finish exactly at the compressed payload end,
remain within the uncompressed block bound and agree with both CRC32 and ISIZE.

The default EOF policy requires the final validated member to be empty. Any valid
terminal empty member qualifies; its compressed bytes need not match the canonical
marker. Explicit legacy permission accepts a physical EOF after a nonempty member
or an empty physical stream and reports `Missing`. It never excuses a partial
member, corrupt trailer, ordinary gzip input or bytes after an embedded marker.

Construction is lazy. EOF status is `Unchecked` until physical EOF is observed.
`finish` drains and validates the remaining input. After a random seek it validates
only the remaining suffix, not unread earlier members. Returning the underlying
reader without finishing does not establish complete validation.

## Offsets and seeking

`BgzfVirtualOffset` validates both components, exposes the raw `u64` encoding and
supports comparison without implementing arithmetic. An offset value alone does
not establish that its compressed component points to a real block boundary.

Positions default to compressed offset zero at the supplied stream start. Explicit
origins support a caller-positioned source. A seekable underlying stream must use
the same byte coordinate system; ordinary files should start at physical zero.

`seek_virtual` is available only for `Read + Seek`. It checks the target against
physical length, decodes its block and checks its uncompressed offset. Exact
ends of short blocks normalize to the next block at offset zero. The end of a
65,536-byte data block is representable only by that next-block position. Gaps
between members, middle-of-header offsets, beyond-EOF positions and within-block
offsets beyond decoded data fail. A decodable embedded member is not proof of a
source-wide authenticated boundary; indexing must verify the complete source.

## Writer contract

`BgzfWriter<W: Write>` preserves arbitrary input bytes. Defaults select compression
level 6 and 65,280-byte data blocks; explicit block sizes 1..=65,536 and compression
levels 0..=9 are validated before writing. If a candidate's compressed framing
would exceed the format bound, it is split and recompressed under bounded buffers.
Every accepted input byte appears once and in order.

`flush_blocks` emits pending data and flushes the sink without terminating it.
`virtual_position` also flushes pending data and returns a stable next-block
position. This explicit choice avoids publishing offsets into a candidate that
later requires splitting. It may alter block boundaries and compression ratio.
The writer does not promise identical compressed bytes or boundaries to zlib,
Biopython or HTSlib.

Consuming `finish` writes pending data, emits the canonical 28-byte EOF marker,
flushes the sink and returns it only on success. Dropping a writer does not
silently finalize it. Write/finalization errors are observable and can leave a
partial underlying output. `write_bgzf_new(path, options, callback)` now uses the
existing create-only transaction and publishes only after successful finalization;
see the [compressed-output contract](compressed-output-v1.md).

## Worked Rust example

Literal input is the eight bytes `ACGTACGT`. This complete program writes two
four-byte blocks plus the terminal empty member, then seeks within the first
block. It was executed against the isolated pinned Rust modules. Public exports
are wired in the main checkout, where the guide programs and executable recipes
also pass.

```rust
use bioferris::io::bgzf::{BgzfReader, BgzfVirtualOffset, BgzfWriteOptions, BgzfWriter};
use std::io::{Cursor, Read, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = BgzfWriter::with_options(Vec::new(), BgzfWriteOptions {
        block_size: 4,
        ..BgzfWriteOptions::default()
    })?;
    writer.write_all(b"ACGTACGT")?;
    let encoded = writer.finish()?;
    let mut reader = BgzfReader::new(Cursor::new(encoded));
    reader.seek_virtual(BgzfVirtualOffset::new(0, 2)?)?;
    let mut suffix = String::new();
    reader.read_to_string(&mut suffix)?;
    println!("Suffix: {suffix}");
    println!("EOF: {:?}", reader.finish()?);
    let invalid = BgzfVirtualOffset::new(0, 65_536).unwrap_err();
    println!("Invalid offset: {}", invalid.code());
    Ok(())
}
```

Observed stdout, exit 0, empty stderr:

```text
Suffix: GTACGT
EOF: Present
Invalid offset: invalid_bgzf_virtual_offset
```

The second offset component is a byte position within one block, not an absolute
source position. A 65,536-byte block's end uses the next member at offset zero,
never a wrapped 16-bit component.

## Errors and limits

Typed errors distinguish source/sink I/O, truncation, unsupported headers, malformed
extra fields, block sizes, raw DEFLATE, expansion limits, CRC, ISIZE, EOF policy,
virtual components, seek targets, invalid options and failed stream state.
They retain block offsets and observed/expected values where applicable.
Standard `io::Error` adapters retain the typed cause for downcasting. After a
failed operation the stream refuses further decoding/writing until a successful
explicit seek recovers a reader. No failed block's unverified bytes are exposed.

Bounds constrain individual compressed/uncompressed buffers, not total streamed
input/output volume or exact allocator overhead. Compressed address arithmetic
must remain within 48 bits, including canonical next-block positions. Empty
output is represented by the canonical EOF member; all-zero input bytes are data,
whereas an all-zero compressed file is invalid.

## Acceptance

Verification must cover independent stored-DEFLATE/CRC expectations, public native
BGZF data, interoperability through permitted Biopython 1.88 public API executions,
all block boundaries and virtual-offset extremes, short reads/writes, interruptions,
empty/embedded/final markers, malformed fields/trailers, expansion bombs, failed
I/O and explicit finalization. Public Rust interfaces, examples and
main-checkout verification support accepted documented scope in the
main validation record.
Earlier isolated/staged compression/index checks remain separate checkpoints.

Prior isolated stream evidence: 28 Rust tests passed, including 16,065 byte
mutations. The combined harness also checks create-only publication, 16 index
tests, 11 gzip tests and the shared CLI adapter. Strict standalone all-target
pedantic Clippy and formatting passed at that checkpoint. The retained
[fixture provenance](../../seq-bio/tests/fixtures/bgzf-provenance.md) identifies
public/native data and both directions of observed Biopython interoperability.
These historical counts are not final main-checkout results. The
[`compressed_io.rs`](../../seq-bio/examples/compressed_io.rs) example and its
[recipe checks](../../seq-bio/tests/compression_recipes.rs) are promoted, with
complete workflows in the [manual guide](../api/compressed-io-v1.md); the CLI
process checks were removed with the CLI.
Main-checkout execution, documentation checks and complete workspace harness
output/terminal audit now support the accepted scope. Final source/evidence
review and the local commit accompany this acceptance; isolated sequence-index
evidence remains separately retained.
