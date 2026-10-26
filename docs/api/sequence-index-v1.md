# Persistent sequence indexes

## Capability index

| Capability | Rust |
|---|---|
| Build and validate a source index | `SequenceIndex::<Dna/Rna/Protein>::build` |
| Publish a new JSON sidecar | `index.save_new(path)` |
| Publish a new binary sidecar | `index.save_new_as(path, SidecarFormat::Binary)` |
| Fully verify and reopen either sidecar | `SequenceIndex::open(source, sidecar, options)` |
| Reopen a binary sidecar with a quick source check | `SequenceIndex::open_with(source, sidecar, options, SourceCheck::Quick)` |
| Seek and parse one record | `index.get(id)` |
| Source-order identifiers / metadata | `ids`, `len`, `is_empty`, `format`, `source_len`, `source_sha256` |

The `index` and `get-record` CLI operations and their JSON requests were removed
with the `bioferris` CLI (decision 029).

## Contents

- [Scope and identity](#scope-and-identity)
- [Small Rust example](#small-rust-example)
- [Persistence and verification](#persistence-and-verification)
- [Binary sidecars and quick checks](#binary-sidecars-and-quick-checks)
- [Limits](#limits)
- [Errors and recovery](#errors-and-recovery)
- [Compatibility and verification](#compatibility-and-verification)

[API index](../README.md#api-index) · [Complete worked recipes](../examples/sequence-index-recipes.md) ·
[Streaming file contract](file-io-v1.md)

## Scope and identity

Index **one stable, uncompressed regular file** with an explicit molecule and
`IndexedFormat::Fasta` or `IndexedFormat::FastqSanger`. There is no format,
molecule or quality-encoding guessing. Existing streaming reader grammar and
validation remain authoritative. FASTA allows empty records; Sanger FASTQ
requires nonempty sequences and exact Phred+33 quality lengths. Wrapped lines,
CRLF, missing final newline and Unicode headers retain their existing meanings.

Keys are case-sensitive `Record::id()` values: the first whitespace-delimited
header token, not the full description. Duplicate keys are errors, never replaced
or renamed. `ids()` borrows identifiers in source order and supports exact length
and reverse iteration. `get` takes `&mut self` because the retained file is seeked.
Missing keys return `Ok(None)` after source metadata checks.

Indexed retrieval equals sequential reading for header, ID, description,
original sequence case, normalized typed sequence and decoded quality scores.
Original wrapping, line endings and repeated FASTQ `+` text are not fields in
`Record` and are not reproduced by its writers. Source data is never modified.
FASTA output refuses quality loss; FASTQ output retains scores.

The public module is `bioferris::io::index`. `IndexAlphabet` is sealed to DNA,
RNA and protein; the existing alignment alphabet and legacy eager FASTA API are
unchanged. Paths can contain non-UTF-8 bytes.

For stable BGZF FASTA/Sanger FASTQ, use the separate
[`BgzfSequenceIndex` workflow](compressed-io-v1.md#write-a-rust-analysis).
This plain-file API continues to
reject compressed sources; ordinary gzip does not acquire BGZF virtual offsets
or verified sequence-index support.

## Small Rust example

Run from the repository root. The checked-in
[input file](../../seq-bio/tests/fixtures/index-example.fastq) contains exactly:

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

Save this complete program as `seq-bio/examples/index_lookup.rs`. It builds a
verified in-memory index without writing a sidecar. This is a user-created example
name; the checked-in persistence example is `indexed_records`.

```rust
use bioferris::io::index::{IndexOptions, IndexedFormat, SequenceIndex};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut index = SequenceIndex::<Dna>::build(
        "seq-bio/tests/fixtures/index-example.fastq",
        IndexedFormat::FastqSanger,
        IndexOptions::default(),
    )?;
    let record = index.get("read2")?.ok_or("missing read2")?;
    println!("{}: {}", record.id(), record.sequence());
    println!("Source: {}", record.source_sequence());
    println!("Phred: {:?}", record.qualities());
    assert_eq!(index.ids().collect::<Vec<_>>(), ["read1", "read2"]);
    assert_eq!(record.source_sequence(), "Nn");
    assert_eq!(record.qualities(), Some([31, 10].as_slice()));
    assert!(index.get("absent")?.is_none());
    Ok(())
}
```

Run from the repository root:

```bash
cargo bio index_lookup --offline --locked
```

Observed stdout (the executable guide check also verifies the assertions):

```text
read2: NN
Source: Nn
Phred: Some([31, 10])
```

`Nn` is preserved while the computational sequence is `NN`. Quality characters
`@` and `+` represent scores 31 and 10; they are not interpreted as record
delimiters. An absent ID is not a fabricated empty record.

For the complete build/save/reopen program, its `cargo bio indexed_records`
invocation, observed output and a create-only failure, use the
[Rust persistence recipe](../examples/sequence-index-recipes.md#manual-rust-program).
This example does not invoke Python or a reference executable.

## Persistence and verification

`save_new` writes compact deterministic JSON without a trailing newline;
`save_new_as(path, SidecarFormat::Binary)` writes the
[binary layout](#binary-sidecars-and-quick-checks) instead. Both use a
same-directory temporary file, synchronization of its contents, and a create-only
hard link. It refuses existing destinations, including dangling symlinks. Failed
publication cleans up its temporary file and does not replace a prior destination.
The destination directory must be trusted; no hostile-directory or power-loss
durability guarantee is made. No overwrite or copy fallback exists.
The temporary name is also checked against the destination name; a caller may
use a filename resembling the internal temporary convention without exposing a
partially written final file.

The JSON sidecar is a BioFerris JSON object, not a SQLite database:

| Field | Contract |
|---|---|
| `format_id` / `schema_version` | `bioferris-sequence-index` / integer `1` |
| `molecule` / `source_format` | `dna`, `rna` or `protein`; `fasta` or `fastq-sanger` |
| `key_policy` | `record-id-first-token-v1` |
| `source_len` / `source_sha256` | Raw byte length and 64-character lowercase hexadecimal SHA-256 |
| `entries` | Source-ordered objects with `id`, `offset`, `byte_length`, `sequence_length`, `record`, `header_line`, `sha256` |

Offsets are zero-based raw bytes; lengths count bytes or biological symbols as
named. Record ordinals and physical header lines are one-based. Each entry digest
covers its raw span, including line endings. Leading FASTA blanks precede its
first span; trailing blanks before the next header belong to the previous span.
All bytes, including such gaps, are covered by the full-source digest.
Unknown/duplicate fields, positional arrays, invalid types/versions/digests,
duplicate IDs, overlapping/out-of-file spans and inconsistent metadata fail.
The source path is **not stored or followed from the sidecar**: callers supply
both source and sidecar explicitly.

Building parses every record and fingerprints the source and spans. Opening
first checks the complete source digest, then reparses and rederives all entries
for exact comparison. A correct whole-file digest alone cannot validate a forged
offset or a truncated-but-parseable FASTA record. Consequently, build and open
scan the source, and reopening may cost more than building. This full check is
the default for both encodings. The only alternative is the opt-in
[quick check](#binary-sidecars-and-quick-checks) of a binary sidecar; there is no
unchecked load.

After opening, lookups perform ID search, a direct seek, bounded record parsing,
span hashing and source metadata checks. Keep one index open for repeated
lookups; every full reopening repeats the whole-source verification.
A built index, or one opened from a JSON sidecar, holds every entry in memory: the
identifiers, their sorted order and 64 bytes per record. One opened from a binary
sidecar holds only the identifiers and their order, and reads the selected
record's 64-byte entry from the sidecar at each lookup. Neither holds the source
sequences; construction still holds one bounded record at a time.

Metadata checks cover the retained handle and the original absolute path before
and after scans/lookups. Unix checks additionally compare device, inode and
change timestamps. Selected-span bytes, parsed ID and sequence length are checked
at retrieval. Detected changes are errors, including path replacement/unlink.
Metadata-only changes can conservatively invalidate an existing open handle;
reopening against identical bytes can succeed.

**Sources must remain stable.** These checks are not a filesystem snapshot,
mutation lock, authenticity check or defense against every hostile race. A
sidecar and source modified together can describe different valid data; SHA-256
does not identify who produced it. Ordinary special-file types are rejected
before opening, but hostile path replacement between checks is unsupported.

## Binary sidecars and quick checks

A binary sidecar (schema 2) holds the same entries as the JSON one in less than
half the space, plus the source modification time and a SHA-256 of the sidecar
itself. An index opened from it keeps only the identifiers and their order in
memory: `b + 17n` bytes on a 64-bit target for `n` records whose IDs total `b`
bytes. On the publication benchmark's 1-million-read FASTQ, opening it and
fetching 10,000 reads peaked at 35 MiB, where the JSON sidecar needed 384 MiB
(validation record). Use it for
large sources; both encodings return the same records. The
[layout specification](../algorithms/binary-sequence-index-v1.md) lists every
byte. The complete program, its invocation and its output are in the
[binary sidecar recipe](../examples/sequence-index-recipes.md#binary-sidecar-and-quick-check).

`open` and `open_with(..., SourceCheck::Full)` check a binary sidecar's source
exactly as a JSON sidecar's. `open_with(..., SourceCheck::Quick)` is opt-in:

| | Full check (default) | Quick check |
|---|---|---|
| Sidecar | Fully validated, including its checksum | The same |
| Source bytes read | All of them, twice | Up to 32 sampled records |
| Source verified | Whole-file SHA-256, every record reparsed and compared | Length, modification time, and each sampled record's SHA-256, parse, ID and sequence length |
| JSON sidecars | Accepted | Refused: `invalid_index` |
| Opening 1 million reads and fetching 10,000 | 3.1 s, 35 MiB | 0.3 s, 35 MiB |

Samples are the first and the last record and others evenly spaced by record
number; the specification gives the exact positions.

### What a quick check misses

A quick check trusts that the source changed only if its length or modification
time did. It does not detect:

- a change that keeps the length and keeps or restores the modification time,
  for example with `touch -r`, a tool that preserves times, or a write within a
  coarse filesystem time resolution, unless it touches a sampled record;
- a different file of the same length copied over the source with the original
  time, unless a sampled record differs;
- a sidecar rewritten consistently, with a recomputed checksum: the checksum
  detects damage, not deliberate changes.

When such a change is missed:

- `get` still never returns changed bytes: every lookup verifies the record's
  span digest, so a changed record fails with `stale_source`;
- an identifier added by the change is reported missing, and one removed by it
  fails when fetched;
- `len`, `ids` and `source_sha256` describe the indexed content, not the current
  file; after a quick check `source_sha256` is the saved value, not a fresh hash.

Use the quick check for sources that cannot change unnoticed, such as read-only
archived reads that a full check has verified once. A refusal is never a reason
to skip checks: if a quick check refuses a source whose time changed, for example
after copying it, open it with the full check and, if that passes, save a new
binary sidecar from the opened index to record the new time.

## Limits

`IndexOptions` is copied into the index; limits are inclusive and zero is valid:

| Budget | Default | Scope |
|---|---|---|
| `max_records` | 100,000 | Number of retained entries |
| `max_id_bytes` | 16,777,216 | Aggregate UTF-8 bytes of keys, not full headers |
| `max_index_bytes` | 67,108,864 | JSON form of the entries while building; the sidecar file when opening; bytes written when saving |
| `read_limits.max_line_bytes` | 8,388,608 | Physical line including line ending |
| `read_limits.max_record_symbols` | 100,000,000 | Symbols in any one record |
| `read_limits.max_total_bytes` | 1,073,741,824 | Whole uncompressed source |

Build checks entry, aggregate ID and serialized-size budgets before retaining
each new entry. Even an empty index has a nonzero JSON envelope, so a zero-byte
sidecar budget fails. A zero-entry/zero-ID budget allows an empty source if the
remaining limits suffice. Open bounds sidecar bytes before parsing and validates
entry/ID budgets before constructing lookup structures. Whitespace in a manually
reformatted sidecar counts toward its loaded-file byte budget. The JSON form is
always larger than the binary one, so an index built within the budget can be
saved in either encoding. Saving refuses a sidecar over the budget, which can
happen only when an index opened from a binary sidecar is saved as JSON; no file
is left behind.

These budgets are not a precise peak-memory ceiling. Parsing, validation,
duplicate detection and verified reopening may temporarily own additional
metadata. Checked arithmetic rejects overflow.

## Errors and recovery

`IndexError` implements `std::error::Error`, exposes `code()`, and is
non-exhaustive. Existing `BioIoError` variants are not expanded.

| Code | Meaning / response |
|---|---|
| Existing reader codes | Invalid header/residue/quality, truncated input or reader limit; fix source/profile, retaining record/line context |
| `duplicate_record_id` | Two first-token IDs collide; choose a corrected source, not a silent overwrite |
| `unsupported_index_source` | Nonregular or gzip/BGZF source; use a stable plain regular file |
| `invalid_index` | Malformed, incompatible or forged sidecar; rebuild from the intended stable source |
| `stale_source` | Observed source identity/bytes changed; investigate and explicitly rebuild |
| `limit_exceeded` | A reader/index budget is exceeded; choose justified larger limits or smaller input |
| `count_overflow` | Size cannot be represented safely; do not truncate it |
| `output_error` | Existing destination or unavailable temporary name; choose a new trusted path |
| `io_error` | Filesystem/read/write/allocation failure; inspect permissions, paths or storage |

A missing key is not an error: `get` returns `Ok(None)`.

Messages specific to binary sidecars and quick checks:

| Code and message | Meaning / response |
|---|---|
| `invalid_index`: quick source checks need a binary sidecar, which records the source modification time | Open the JSON sidecar with the full check, or save a binary sidecar from an index opened with it |
| `invalid_index`: binary sidecar checksum does not match its contents; length does not match its record and ID counts; shorter than its header and checksum | Damaged or truncated sidecar; rebuild it from the stable source |
| `invalid_index`: sidecar identity, size or timestamps differ from the verified file | The binary sidecar changed after opening; reopen it |
| `stale_source`: source modification time differs from the index | The source may still hold the same bytes, for example after a copy; open it with the full check and, if that passes, save a new binary sidecar |
| `stale_source`: record span bytes differ from the verified SHA-256 | A sampled or fetched record changed, or the sidecar's digest was altered; run the full check and rebuild |

Index budget details include resource, inclusive limit and attempted size.
Duplicate-ID details include the ID and both one-based record ordinals.
Lookup does not silently use cached records after detecting changed source data.
No automatic deletion, repair, index replacement or source modification occurs.

## Compatibility and verification

Functional reference: [Biopython 1.88 SeqIO indexing](https://biopython.org/docs/1.88/api/Bio.SeqIO.html#Bio.SeqIO.index),
consulted 2026-10-16. This implements the source-ordered first-token ID / seek-and-
parse workflow for an explicit FASTA/Sanger FASTQ subset. It is not Python API,
SQLite `index_db`, serialization or performance equivalence. Biopython's wider
formats, BGZF, multi-source databases and custom key functions remain outside
this plain-file `SequenceIndex` profile. BGZF has the separate API linked above.
No Biopython runtime is required or executed for these tests.

JSON sidecars written by earlier versions open unchanged, and `save_new` writes
the same bytes as before; a test compares its output with a sidecar written by
`096a948b`. Binary sidecars and `SourceCheck` were added on 2026-10-22
(decision 034);
earlier versions reject a binary sidecar as invalid JSON.

Fingerprints use pinned RustCrypto `sha2` 0.10.9 with `std` and `force-soft`,
without its assembly backend; seven lockfile packages were added without
upgrading prior dependencies or the pinned Rust 1.85.0 compiler. The project
introduces no unsafe code or FFI. Maintained dependencies retain their own
internal implementations; this is not a claim that all dependency code is unsafe-free.

Acceptance includes independent sequential/indexed equality across all three
molecules, generated wrapping/quality cases, literal raw offsets, real phiX174
data, exact/zero budgets, corrupt schema and forged spans, stale files, special
files and create-only cleanup. Every round trip also reopens a binary sidecar with
both checks. Runtime evidence, scope and phase measurements are recorded in the
validation report; the binary sidecar and
quick check have their own record.
Measured opening/retrieval costs must not be described as a universal speedup.
