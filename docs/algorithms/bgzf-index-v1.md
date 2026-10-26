# Verified BGZF sequence indexes

## Contents

- [Scope and status](#scope-and-status)
- [Source and manifest](#source-and-manifest)
- [Coordinates and empty blocks](#coordinates-and-empty-blocks)
- [Verification and retrieval](#verification-and-retrieval)
- [Rust interface and example](#rust-interface-and-example)
- [Limits and acceptance](#limits-and-acceptance)

[API index](../README.md) · [BGZF contract](bgzf-v1.md) ·
[Plain indexes](../api/sequence-index-v1.md)

## Scope and status

The additive core and public `io::bgzf_index` exports are implemented and
verified; the `index-bgzf`/`get-bgzf-record` CLI and strict-JSON wiring were
removed with the `bioferris` CLI (decision 029). Stage 6 is accepted for the
documented BGZF stream/index profiles; see the
main validation record.
Prior isolated/staging results remain separate checkpoints.
`BgzfSequenceIndex<A>` supports
one stable regular BGZF file with an explicit DNA/RNA/protein alphabet and the
existing FASTA or Sanger FASTQ profile. It does not change `SequenceIndex`, index
ordinary gzip, add a BAM parser, implement SQLite or promise constant-time open.

The functional reference is public
[Biopython 1.88 SeqIO indexing](https://biopython.org/docs/1.88/api/Bio.SeqIO.html#Bio.SeqIO.index).
Retained native public API observations compare seven keys, raw spans,
case/headers, qualities and block tables across three fixtures; see
[provenance](../../seq-bio/tests/fixtures/bgzf-provenance.md#native-sequence-index-observations).
They are not source copying or runtime dependencies.

## Source and manifest

Reuse the existing `IndexAlphabet`, `IndexedFormat`, `IndexOptions`, bounded
FASTA/Sanger FASTQ readers and their private raw-span observations. Keys remain
case-sensitive first header tokens; duplicates fail. Records preserve complete
headers, source case, typed sequence and decoded Phred scores.

A distinct object-only schema-version-1 sidecar identifies BGZF explicitly. It
retains compressed source length/SHA-256, decompressed length, every validated
block boundary and source-ordered entries. Entries contain original decompressed
raw start/length, canonical virtual start/end, sequence length, record/header-line
ordinals and SHA-256 of the complete decompressed raw span. Paths are supplied by
the caller and never followed from a sidecar.

Unknown/duplicate fields, positional arrays, wrong versions/molecules/profiles,
invalid digests, duplicated keys, malformed block tables, out-of-bounds/overlapping
spans and inconsistent source metadata fail. Persisted offsets are not trusted
merely because a block can decode at that address.

## Coordinates and empty blocks

Block metadata maps decompressed absolute byte positions to standard virtual
positions using binary search. At a logical boundary shared by empty members,
the canonical index encoding chooses the earliest physical block beginning at
that logical offset. Ends of 65,536-byte blocks use the next block at offset zero;
65,536 is never stored in the 16-bit component. Physical EOF at offset zero is
also a readable end-position alias after a valid final marker.

The inverse map accepts structurally valid equivalent boundary positions against
the fully verified table. Retrieval compares logical end coordinates, since a
decoder may normalize an empty-member target to its next block. Persisted entry
encodings must still equal the exact canonical rederivation; aliases are not a
permission to forge sidecar offsets.

## Verification and retrieval

Build validates every compressed block, including empty members and final EOF,
and fingerprints all compressed bytes. It parses every record and independently
hashes every full decompressed raw span. Open redoes this work and compares the
entire manifest, including the block table, before exposing any lookup.

Lookups retain a file handle, check observed source identity/size/timestamps,
seek to the verified virtual start, decode a bounded raw span, verify its SHA-256,
parse exactly one record through the existing reader and compare ID/length/end.
Missing IDs return `None` after metadata checks. No partial record is returned.

Sources must remain stable. Metadata/span checks are not an atomic filesystem
snapshot, a mutation lock, authentication or protection against every hostile
concurrent change. The original source path is checked alongside the retained
handle, matching the plain-index stability policy. Sidecar publication reuses the
existing create-only transaction and never replaces an existing destination.

## Rust interface and example

`BgzfSequenceIndex::<Dna|Rna|Protein>::build(source, format, options)` verifies
without creating a sidecar. `save_new(path)` publishes it, and
`open(source, sidecar, options)` fully revalidates it. `get(&mut self, id)` returns
`Result<Option<Record<A>>, BgzfIndexError>`; `ids()` borrows source-ordered keys.
Getters include `len`, `is_empty`, `format`, compressed `source_len`, decompressed
`uncompressed_len`, `source_sha256` and `block_count`.

`virtual_offset(uncompressed_byte_position)` and
`uncompressed_offset(BgzfVirtualOffset)` perform checked binary-search mappings.
These in-memory mappings do not recheck filesystem metadata; lookups and sidecar
publication do. Never add/subtract packed offsets as though they were bytes.

Literal input is the FASTA record `>gene demo\naCgT\n`. Run the complete program
in a new trusted empty directory. It was executed against isolated Rust modules;
those public exports are wired in the main checkout, where the guide programs
and executable recipes also pass.

```rust
use bioferris::{Dna, Record};
use bioferris::io::{LossPolicy, write_record_fasta};
use bioferris::io::bgzf::{BgzfWriteOptions, write_bgzf_new};
use bioferris::io::bgzf_index::{BgzfIndexOptions, BgzfSequenceIndex};
use bioferris::io::index::IndexedFormat;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let record = Record::<Dna>::from_ascii("gene demo", "aCgT", None)?;
    write_bgzf_new("example.bgz", BgzfWriteOptions::default(), |output| {
        write_record_fasta(output, &record, LossPolicy::Reject)
    })?;
    let mut index = BgzfSequenceIndex::<Dna>::build(
        "example.bgz", IndexedFormat::Fasta, BgzfIndexOptions::default(),
    )?;
    println!("IDs: {}", index.ids().collect::<Vec<_>>().join(","));
    let found = index.get("gene")?.ok_or("missing expected record")?;
    println!("Source case: {}", found.source_sequence());
    println!("Normalized: {}", found.sequence());
    println!("Blocks: {}", index.block_count());
    println!("Missing: {}", index.get("absent")?.is_none());
    index.save_new("example.index.json")?;
    let mut reopened = BgzfSequenceIndex::<Dna>::open(
        "example.bgz", "example.index.json", BgzfIndexOptions::default(),
    )?;
    println!("Reopened: {}", reopened.get("gene")?.ok_or("missing expected record")?.source_sequence());
    Ok(())
}
```

Observed stdout, exit 0, empty stderr:

```text
IDs: gene
Source case: aCgT
Normalized: ACGT
Blocks: 2
Missing: true
Reopened: aCgT
```

Block counts include the final empty member. Original case is separate from
normalized computational symbols. Repeating this example in the same directory
refuses `example.bgz`; the first output is not overwritten.

## Limits and acceptance

Existing reader/entry/ID/sidecar budgets remain explicit. Additional budgets bound
compressed source bytes and retained block metadata, including empty blocks.
Checked arithmetic covers all byte totals, source offsets and virtual encodings.
Default decompressed total is the existing 1 GiB reader limit; compressed bytes
default to 1 GiB and blocks to 1,000,000. These are not exact peak-memory claims.

Acceptance includes native and independently constructed BGZF FASTA/FASTQ,
records and headers crossing blocks, CRLF/wrapping/empty FASTA, first-token keys,
metadata/quality equality, duplicate IDs, every exact/zero budget, ordinary gzip
refusal, forged sidecars with correct source digests, empty-member/65,536-byte
boundaries, stale/replaced sources and create-only cleanup. Public exports and
worked examples are wired in the main checkout; their main verification supports
accepted documented scope in the
main validation record. The promoted
[`compressed_io.rs`](../../seq-bio/examples/compressed_io.rs) and its
[recipe checks](../../seq-bio/tests/compression_recipes.rs) implement the
[manual workflows](../api/compressed-io-v1.md).

The main tests, documentation checks and complete workspace harness output with
terminal audit are recorded separately from the isolated checkpoint below.
Final source/evidence review and the local commit accompany this acceptance.

The isolated core checkpoint passed 16 tests, including the seven-record native comparison
and a 5,000-record shared-member regression. Raw-span fingerprinting uses one
sequential decoding pass, rather than decoding the same member for every key.
`BgzfIndexError::{Index, Bgzf}` retains `IndexError` or `BgzfError` without string
erasure. Duplicate keys give `duplicate_record_id`; schema/forgery errors
`invalid_index`; detected changes `stale_source`; budgets `limit_exceeded`;
malformed BGZF retains its precise framing/integrity code. `Error::source`
preserves the underlying error.
