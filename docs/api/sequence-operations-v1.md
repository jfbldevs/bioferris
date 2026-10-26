# Checked sequence slicing and literal search

## Contents

- [Capability and input contract](#capability-and-input-contract)
- [Run a Rust analysis](#run-a-rust-analysis)
- [Coordinates and concatenation](#coordinates-and-concatenation)
- [Literal search and overlap](#literal-search-and-overlap)
- [Record slicing](#record-slicing)
- [Errors and resource limits](#errors-and-resource-limits)
- [Compatibility and verification](#compatibility-and-verification)

[API capability index](../README.md#api-index) · [Quick start](../quick-start.md)

## Capability and input contract

`SequenceExt<A>` is included in `bioferris::prelude::*`. It extends the existing
immutable `Sequence<Dna>`, `Sequence<Rna>` and `Sequence<Protein>` types without
changing their alphabets or alignment behavior.

| Capability | Rust access |
|---|---|
| Checked region extraction | `sequence.slice(start..end)` |
| Same-alphabet concatenation | `first.concat(&second)` |
| First literal match | `sequence.find(pattern)` |
| Literal containment | `sequence.contains(pattern)` |
| Non-overlapping match count | `sequence.count(pattern)` |
| Overlapping match count | `sequence.count_overlapping(pattern)` |
| Lazy match positions | `sequence.find_all(pattern, overlap)` |
| Case/quality-preserving record extraction | `record.slice(start..end)` |

The manual workflow is a Rust program, as shown below. The `slice`, `concat`,
`find` and `count` CLI operations and their JSON requests were removed with the
`bioferris` CLI (decision 029).

Inputs are fully defined, ungapped ASCII sequences. Constructors and search
patterns normalize letter case to uppercase. DNA rejects `U`, RNA rejects `T`,
and proteins reject stop markers such as `*`; whitespace and non-ASCII bytes are
not silently removed. Valid ambiguity symbols remain literal characters.

## Run a Rust analysis

From the repository root, use the pinned Rust toolchain described in the
[quick start](../quick-start.md). This complete example is checked in as
[`sequence_operations.rs`](../../seq-bio/examples/sequence_operations.rs).
Its literal inputs are DNA `atatatn`, the DNA suffix `gg`, RNA `augn`, protein
`MABZX`, and a DNA record whose original text is `aTaTaTn` with seven explicit
Phred scores.

<!-- sequence-operations-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("atatatn")?;
    println!("DNA: {sequence}");
    println!("Slice 1..5: {}", sequence.slice(1..5)?);
    println!("Concatenated: {}", sequence.concat(&dna("gg")?)?);
    println!("First TA: {:?}", sequence.find("ta")?);
    println!("Contains NN: {}", sequence.contains("NN")?);
    println!("ATA, non-overlapping: {}", sequence.count("ATA")?);
    println!("ATA, overlapping: {}", sequence.count_overlapping("ATA")?);
    let matches: Vec<_> = sequence
        .find_all("ATA", MatchOverlap::Overlapping)?
        .collect();
    println!("ATA positions: {matches:?}");
    println!("Empty-pattern count: {}", sequence.count("")?);
    println!(
        "Terminal empty slice: {} symbols",
        sequence.slice(7..7)?.len()
    );
    println!("RNA slice: {}", rna("augn")?.slice(1..3)?);
    println!("Protein BZ: {:?}", protein("MABZX")?.find("bz")?);

    let record = Record::<Dna>::from_ascii(
        "read1 mixed case",
        "aTaTaTn",
        Some(vec![30, 31, 32, 33, 34, 35, 36]),
    )?;
    let region = record.slice(1..5)?;
    println!("Record header: {}", region.header());
    println!("Record source: {}", region.source_sequence());
    println!("Record normalized: {}", region.sequence());
    println!("Record qualities: {:?}", region.qualities());
    println!("Original source: {}", record.source_sequence());

    if let Err(error) = sequence.slice(7..8) {
        println!("Invalid slice: {}", error.code());
    }
    if let Err(error) = sequence.find("AU") {
        println!("Invalid DNA pattern: {}", error.code());
    }
    Ok(())
}
```

Complete invocation:

<!-- sequence-operations-command -->
```bash
cargo bio sequence_operations
```

Complete stdout, exit 0, no stderr:

<!-- sequence-operations-output -->
```text
DNA: ATATATN
Slice 1..5: TATA
Concatenated: ATATATNGG
First TA: Some(1)
Contains NN: false
ATA, non-overlapping: 1
ATA, overlapping: 2
ATA positions: [0, 2]
Empty-pattern count: 8
Terminal empty slice: 0 symbols
RNA slice: UG
Protein BZ: Some(2)
Record header: read1 mixed case
Record source: TaTa
Record normalized: TATA
Record qualities: Some([31, 32, 33, 34])
Original source: aTaTaTn
Invalid slice: invalid_range
Invalid DNA pattern: invalid_sequence
```

Positions are zero-based. The overlapping `ATA` matches share the symbol at
position 2; non-overlapping counting consumes the first match before continuing.
`N` is not a wildcard, so the one terminal `N` does not satisfy pattern `NN`.
The last two lines demonstrate handled errors: the program succeeds because it
reports those deliberately invalid requests rather than propagating them.

## Coordinates and concatenation

Ranges are Rust's zero-based, half-open `start..end`: the start is included and
the end is excluded. Every slice must satisfy `start <= end <= sequence.len()`.
`0..0`, `len..len` and `0..len` are valid, including for an empty sequence.
Out-of-range coordinates are rejected, not clipped. There are no negative
indices, implicit defaults, reversed traversal or strides in this API.

Slicing produces an owned typed sequence. Concatenation produces an owned typed
sequence containing all of the first operand followed by the second. Empty
operands are valid. Neither operation mutates either original sequence.

Concatenation requires identical compile-time alphabets. This complete program
intentionally fails compilation: a DNA value cannot concatenate an RNA value,
even if their present symbols happen to be shared by both alphabets.

```compile_fail
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let first = dna("AC")?;
    let second = rna("AC")?;
    let _combined = first.concat(&second)?;
    Ok(())
}
```

The compiler reports incompatible DNA/RNA types; no biological output is
produced. Use an explicit biological conversion, such as back-transcription,
only when that transformation is actually intended.

## Literal search and overlap

`find` returns `Result<Option<usize>, SequenceOperationError>`: `Some(position)`
for the first match, `None` for a valid pattern that is absent, or `Err` for an
invalid pattern. `contains` returns the corresponding boolean. `count` returns
the non-overlapping count; `count_overlapping` allows consecutive matches to
share symbols. Both counts are `usize` values inside a `Result`.

`find_all` validates and owns the normalized pattern, borrows the sequence and
returns `SequenceMatches<'_, A>`. Choose `MatchOverlap::NonOverlapping` or
`MatchOverlap::Overlapping` explicitly. It yields `usize` start positions in
ascending order and is fused: once exhausted, subsequent calls return `None`.
It does not precompute or store every match. Dropping the iterator stops scanning;
collecting it into a `Vec` is an explicit caller allocation.

Literal equality applies to every accepted ambiguity symbol. For example,
pattern `N` matches a stored `N`, not every nucleotide, and protein `J` does not
match stored `I` or `L`. These are not motif, regular-expression, approximate
alignment or ambiguity-expansion searches. Lowercase input patterns are
normalized before matching, consistently with the sequence constructors.

An empty pattern matches every boundary `0..=len`, regardless of overlap mode.
Therefore `find("")` is `Some(0)`, `contains("")` is `true`, and both counts are
`len + 1`. An empty source has one such boundary. A nonempty pattern longer than
the source has no matches, but is still fully validated before returning.

## Record slicing

`Record<A>::slice(range)` uses the same checked coordinates and returns a new
record. It preserves the complete original header, parsed ID and description.
It slices the original sequence text without changing letter case, slices any
quality vector with the same range, and retains the normalized typed sequence.
The source record remains unchanged. Missing qualities remain `None`; a valid
empty slice of a record with qualities retains `Some(empty_vector)`.

Headers are retained as text, not interpreted as genomic coordinates or
rewritten to describe the new region. This method does not slice annotated
GenBank/EMBL records, shift feature locations, merge records or infer annotation
loss policies. Those transformations have separate acceptance requirements.
An empty record is representable in memory, but the Sanger FASTQ writer still
rejects an empty read under its existing [file profile](file-io-v1.md).

## Errors and resource limits

| Rust error | Stable `code()` | Meaning |
|---|---|---|
| `InvalidRange { start, end, length }` | `invalid_range` | Reversed or out-of-bounds range; supplied coordinates and source length are retained |
| `InvalidSequence(SequenceError)` | `invalid_sequence` | Invalid pattern byte; the inner error reports alphabet, byte and zero-based pattern offset |
| `LengthOverflow` | `length_overflow` | A concatenated length overflows `usize` or exceeds the platform's `isize::MAX` allocation bound |

Errors implement `std::error::Error` and human-readable `Display`. Match typed
variants or stable codes for control flow instead of parsing messages. Search
validation errors are distinct from the valid absence result `None` or count 0.
Alphabet mismatches in `concat` are compile-time errors rather than runtime
`SequenceOperationError` values.

The search implementation uses linear-time prefix preprocessing: for source
length `n` and pattern length `m`, a full scan is O(n + m) with O(m) auxiliary
memory. Match positions are lazy; collecting `k` matches uses another O(k)
space. Slicing allocates a copy of the selected symbols; concatenation allocates
a copy of both inputs. Record slicing also copies its preserved header metadata,
source text and available qualities. Original sequences are not copied merely
to search them.

These are algorithmic bounds, not measured speed or throughput claims. There is
no configurable search-memory cap in this increment. Inputs and allocated
outputs must fit memory; ordinary allocator exhaustion is not converted into a
recoverable biological error. File-reader limits remain separate. The combined
length check is a guarded error contract, not a test that allocates an impossibly
large sequence.

## Compatibility and verification

Reference scope is the publicly documented sequence slicing, concatenation and
literal search/count behavior in the [Biopython 1.88 sequence chapter](https://biopython.org/docs/1.88/Tutorial/chapter_seq_objects.html)
and [Bio.Seq API](https://biopython.org/docs/1.88/api/Bio.Seq.html), consulted
2026-10-15. This is a Rust-native subset, not Python signature or object parity.
Intentional differences include strict nonnegative ranges, alphabet-typed
concatenation, uppercase normalization of both sources and patterns, validation
of all pattern bytes, and `Option<usize>` instead of Python's missing-index
sentinel. No reference implementation source was copied and no Python runtime is
required.

The [Rust verification suite](../../seq-bio/tests/sequence_operations.rs) compares
all short DNA sources over `A/C/N` through length 5 against all patterns through
length 4 using an independent naive byte-search oracle. Property tests cover
DNA, RNA and protein search, slicing, concatenation and quality-preserving record
slicing. Deterministic cases cover every byte as a one-symbol pattern, invalid
offsets, empty patterns, repeated prefixes, iterator exhaustion and record
metadata. The suite also runs the actual Cargo example and compares its source
and complete output with this page; guide rustdoc checks the incompatible-alphabet
example as a compilation failure. The existing public phiX174 CDS fixture checks
integration with translation, and a FASTQ region is written and read back with
unchanged source case and matching quality scores.

Verification establishes the stated contracts. It is not a live Biopython
differential campaign, a performance benchmark or evidence of broader scientific
algorithm fidelity. Existing frozen alignment qualification is unchanged.
Commands, outcomes and retained logs are in the
stage 1 validation report.
