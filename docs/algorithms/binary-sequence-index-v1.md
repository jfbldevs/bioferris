# Binary sequence index sidecar

## Contents

- [Scope](#scope)
- [Layout](#layout)
- [Validation when opening](#validation-when-opening)
- [Full source check](#full-source-check)
- [Quick source check](#quick-source-check)
- [Lookups and memory](#lookups-and-memory)
- [Compatibility](#compatibility)

[Sequence index guide](../api/sequence-index-v1.md#binary-sidecars-and-quick-checks) ·
[Worked recipe](../examples/sequence-index-recipes.md#binary-sidecar-and-quick-check) ·
Decision 034 ·
Validation record

Profile identifier: `BINARY-SEQUENCE-INDEX-v1`, sidecar schema version 2, added
on 2026-10-22.

## Scope

A binary sidecar is a second encoding of the plain-file `SequenceIndex` sidecar.
It holds the same entries as the schema-1 JSON object: per record its first-token
ID, raw byte offset and length, sequence length, header line and span SHA-256,
plus the source length and SHA-256. It adds the source modification time, used
only by the quick source check, and a checksum of the sidecar itself. Writing it:
`SequenceIndex::save_new_as(path, SidecarFormat::Binary)`. `SequenceIndex::open`
and `open_with` recognize it from its first eight bytes.

## Layout

All integers are little-endian. For `n` records whose IDs total `b` UTF-8 bytes,
the file is `120 + 73n + b` bytes long.

| Offset | Bytes | Field |
|---:|---:|---|
| 0 | 8 | Magic `89 42 46 49 44 58 0D 0A` (`\x89BFIDX\r\n`) |
| 8 | 4 | Schema version, `2` |
| 12 | 1 | Molecule: `1` DNA, `2` RNA, `3` protein |
| 13 | 1 | Source format: `1` FASTA, `2` Sanger FASTQ |
| 14 | 1 | Key policy: `1`, `record-id-first-token-v1` |
| 15 | 1 | Reserved, `0` |
| 16 | 8 | Record count `n` |
| 24 | 8 | ID bytes `b`, without line feeds |
| 32 | 8 | Source length in bytes |
| 40 | 8 | Source modification time: signed seconds from the Unix epoch |
| 48 | 4 | Source modification time: nanoseconds, below 1,000,000,000 |
| 52 | 4 | Reserved, `0` |
| 56 | 32 | Source SHA-256 |
| 88 | 64 per record | Entries in source order: offset (8), byte length (8), sequence length (8), header line (8), span SHA-256 (32) |
| 88 + 64n | `b + n` | IDs in source order, each followed by a line feed |
| 88 + 65n + b | 8 per record | Lookup: zero-based entry positions in ascending byte order of their IDs |
| 88 + 73n + b | 32 | SHA-256 of every preceding byte |

The record number is the entry's position plus one, so it is not stored. The
modification time counts nanoseconds forward from the stated second, also before
the epoch. Offsets, lengths, header lines and digests follow the
[JSON sidecar rules](../api/sequence-index-v1.md#persistence-and-verification).

## Validation when opening

Both source checks first validate the whole sidecar, in this order; every failure
is `invalid_index` unless stated:

1. The path is a regular file no longer than `max_index_bytes`
   (`limit_exceeded`, resource `index_bytes`).
2. It is at least 120 bytes. The magic, schema version 2 and key policy 1 match;
   the molecule equals the requested alphabet; the format code is known; the
   reserved bytes are zero; the nanoseconds are below one billion.
3. `n` is within `max_records` and `b` within `max_id_bytes` (`limit_exceeded`,
   resources `index_records` and `index_id_bytes`).
4. The length equals `120 + 73n + b`.
5. Each entry follows the JSON entry rules in source order: a nonzero byte length,
   no overlap with the previous span, an end within the source length, a sequence
   length no larger than the byte length, increasing header lines, and nonempty
   sequences in Sanger FASTQ.
6. The ID section is UTF-8 and holds exactly `n` line-feed-terminated IDs, each
   nonempty and free of whitespace and control characters.
7. Every lookup position is below `n`, and the IDs read in lookup order strictly
   increase. This also rules out duplicate IDs and makes the lookup a
   permutation, so every record is reachable.
8. The final 32 bytes equal the SHA-256 of all bytes before them.
9. The sidecar's identity, size and timestamps did not change while reading.

Only the IDs and lookup positions are kept; entries are read again from the file
when needed.

## Full source check

`SourceCheck::Full`, which `open` uses, then checks the source:

1. The source is a regular, uncompressed file of the saved length, within
   `read_limits.max_total_bytes` (`unsupported_index_source`, `stale_source`,
   `limit_exceeded`).
2. One pass hashes the whole source and each saved span. A different source
   digest is `stale_source`; matching source bytes with a different span digest
   is `invalid_index`.
3. A second pass reparses every record with the saved format and the caller's
   reader limits. Each record's ID, offset, byte length, sequence length and
   header line must equal its entry, and the record count must equal `n`;
   otherwise `invalid_index`. A record that no longer parses is `invalid_index`
   ("source does not match the declared profile"); reader I/O and limit errors
   keep their codes.
4. The source's identity, size and timestamps did not change during the check.

## Quick source check

`SourceCheck::Quick` validates the sidecar as above, then:

1. The source is a regular, uncompressed file of the saved length, within
   `read_limits.max_total_bytes`.
2. Its modification time equals the saved one (`stale_source`, "source
   modification time differs from the index").
3. It verifies `S = min(n, 32)` sampled records at positions
   `floor(k * (n - 1) / (S - 1))` for `k = 0, ..., S - 1` (position 0 when
   `S = 1`): the first and the last record and others evenly spaced by record
   number. Each is read as a lookup reads it: the span's SHA-256 must equal the
   saved one (`stale_source`), it must parse (`invalid_index`), and its ID and
   sequence length must equal the entry (`invalid_index`).

It does not hash the whole source, read the other records or the bytes between
records, or check header line numbers. A change that keeps the length and keeps
or restores the modification time, outside the sampled records, passes it. The
[sequence index guide](../api/sequence-index-v1.md#what-a-quick-check-misses)
lists the consequences.

## Lookups and memory

An opened binary sidecar keeps the source and sidecar files open. In memory it
holds the ID section (`b + n` bytes) and two arrays of `n` machine words: each
ID's end, and the lookup positions. That is `b + 17n` bytes on a 64-bit target.

`get(id)` binary-searches the IDs in memory, reads the record's 64-byte entry from
the sidecar, then reads, hashes and parses the source span as for every index.
Before and after, it checks the identity, size and timestamps of both the source
and the sidecar; a changed sidecar is `invalid_index` ("sidecar identity, size or
timestamps differ from the verified file").

Time to open is linear in the sidecar length for the validation, plus linear in
the source length for a full check, or proportional to the sampled records'
sizes for a quick check.

## Compatibility

- The schema-1 JSON sidecar is unchanged; `save_new` still writes it, with the
  same bytes as before, and `open` reads both encodings.
- Versions before this profile reject a binary sidecar as invalid JSON
  (`invalid_index`).
- Saving is deterministic for a given source content and modification time. Two
  copies of the same content with different times give sidecars that differ only
  in the time fields and the checksum; both pass a full check.
- The layout is BioFerris's own. It is not a samtools `.fai` or `.fqi` index or a
  Biopython `SQLite` index.
