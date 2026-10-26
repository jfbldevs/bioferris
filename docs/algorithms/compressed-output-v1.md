# Deterministic compressed output

## Contents

- [Scope and status](#scope-and-status)
- [Gzip writer](#gzip-writer)
- [Create-only publication](#create-only-publication)
- [Worked Rust example](#worked-rust-example)
- [Errors](#errors)
- [Verification requirements](#verification-requirements)

[API index](../README.md) · [BGZF contract](bgzf-v1.md)

## Scope and status

The Stage 6 gzip core and public `io::gzip` exports are implemented and
verified. Stage 6 is accepted for the documented stream/index profiles; see the
main validation record. Its
`compress-file` CLI and strict-JSON adapter were removed with the `bioferris` CLI
(decision 029).
Prior isolated/staging checks remain separate checkpoints.
It preserves arbitrary bytes and supports existing biological
writers through `std::io::Write`. This does not introduce ordinary-gzip seeking
or indexes. BGZF uses its separate block and virtual-position contract.

## Gzip writer

`GzipWriter<W: Write>` uses the existing pure-Rust flate2 1.1.9 raw compressor,
CRC32 and a deterministic [RFC 1952](https://www.rfc-editor.org/rfc/rfc1952)
member header. Default compression level is 6; valid explicit levels are 0..=9.
MTIME is zero, OS is 255 and filename/comment/extra fields are absent. XFL is 2
for level 9, 4 for level 1 and zero otherwise. No clock or path enters output.

The writer uses a bounded reusable compressed-output buffer rather than retaining
the entire input. Typed writes and standard `Write` operations share the same
computation. `flush` emits pending compressed bytes and flushes the sink while
leaving the stream open. Consuming `finish` completes raw DEFLATE, writes CRC32
and ISIZE modulo 2^32, flushes the sink and returns it only on success.

Dropping the writer performs no hidden finalization or output I/O. No compressor
or owned sink is intentionally leaked. A failed write/flush leaves a failed state;
later writes refuse it. Standard I/O errors preserve the typed compression cause.
Existing sink bytes can contain a valid prefix before a later failure.

The output is a single gzip member. Input bytes, case, NULs and line endings are
preserved. Determinism is qualified by the pinned backend, options, write/flush
schedule and input; this does not promise compressed-byte identity with zlib or
other software. The public [flate2 compressor API](https://docs.rs/flate2/1.1.9/flate2/struct.Compress.html)
provides the compression engine; no compressor algorithm is recreated here.

## Create-only publication

`write_gzip_new(path, options, callback)` and `write_bgzf_new` reuse the existing
same-directory temporary file, content synchronization and create-only hard-link
publisher. They refuse existing files and dangling symlinks, retain callback
errors separately from compression/publication failures, finalize before
publishing and remove failed temporary output. The destination directory must be
trusted. No overwrite or copy fallback is permitted.

Existing biological writers compose directly. `write_record_fastq(output,
&record)` preserves the header, sequence case and qualities;
`write_record_fasta(output, &record, LossPolicy::Reject)` refuses to discard
qualities. Selecting `LossPolicy::Discard` is explicit biological loss, not an
effect of compression. Compression never invents quality scores.

## Worked Rust example

Literal input is the five bytes `aCgT\n`. This complete program was executed
against the isolated pinned Rust 1.85 modules. Those public modules are now wired
in the main checkout, where the guide programs and executable recipes also pass.

```rust
use bioferris::io::decode_input;
use bioferris::io::gzip::{GzipWriteOptions, GzipWriter};
use std::io::{Cursor, Read, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = GzipWriter::new(Vec::new());
    writer.write_all(b"aCgT\n")?;
    let encoded = writer.finish()?;
    println!("MTIME: {:?}", &encoded[4..8]);
    let mut decoded = String::new();
    decode_input(Cursor::new(encoded))?.read_to_string(&mut decoded)?;
    print!("Round-trip: {decoded}");
    let invalid = GzipWriteOptions { compression_level: 10 }.validate().unwrap_err();
    println!("Invalid level: {}", invalid.code());
    Ok(())
}
```

Observed stdout, exit 0, empty stderr:

```text
MTIME: [0, 0, 0, 0]
Round-trip: aCgT
Invalid level: invalid_gzip_options
```

Case and newline survive unchanged. `finish` is essential: successful writes or
dropping a writer alone do not produce a verified complete gzip stream.

## Errors

`GzipError` distinguishes `invalid_gzip_options`, `gzip_io_error`,
`gzip_deflate_error`, `gzip_count_overflow` and `gzip_failed_state`.
`GzipOutputError<E>` and `BgzfOutputError<E>` retain callback `E`, typed
compression failures, publication I/O, occupied destinations and unavailable
temporary names separately. `Error::source` retains original causes. Standard
`Write` adapters wrap, not erase, the typed error. Publication failures never
authorize overwrite.

## Verification requirements

Tests cover independent gzip decoding and trailer expectations, all compression
levels, empty/all-byte/large input, repeated flushes, failed/short/interrupted
output, explicit finish, drop without finish and atomic publication refusal/
cleanup. Existing FASTA/Sanger FASTQ writers compose directly with the
compressed writer and retain biological metadata. Prior isolated checks passed
11 gzip tests, 28 BGZF stream tests, 16 index tests, four adapter unit tests and
eight actual-process adapter tests. These historical counts are separate from
final main-checkout verification. The promoted
[`compressed_io.rs`](../../seq-bio/examples/compressed_io.rs) and its
[recipe checks](../../seq-bio/tests/compression_recipes.rs) remain; the CLI
process tests were removed with the CLI.
The [manual guide](../api/compressed-io-v1.md) shows complete workflows;
their main-checkout execution, documentation checks and complete workspace
harness output/terminal audit support accepted Stage 6 scope in the
main validation record.
Final source/evidence review and the local commit accompany this acceptance.
