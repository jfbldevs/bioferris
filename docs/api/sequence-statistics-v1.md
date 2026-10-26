# Sequence statistics and bounded k-mers

## Contents

- [Capabilities and defaults](#capabilities-and-defaults)
- [Run a Rust analysis](#run-a-rust-analysis)
- [Summary fields and aggregation](#summary-fields-and-aggregation)
- [GC and ambiguity policies](#gc-and-ambiguity-policies)
- [K-mer windows and resource budgets](#k-mer-windows-and-resource-budgets)
- [Options](#options)
- [Errors and streaming behavior](#errors-and-streaming-behavior)
- [Compatibility and verification](#compatibility-and-verification)

[API capability index](../README.md#api-index) · [Quick start](../quick-start.md)

## Capabilities and defaults

| Capability | Rust access |
|---|---|
| Length, observed-symbol composition and ambiguity count | `sequence.statistics()` |
| Explicit nucleotide GC policy | `sequence.statistics_with(options)` |
| Optional bounded forward k-mer counts | `StatisticsOptions::kmers` |
| Aggregate borrowed sequences | `summarize_sequences(sequences, options)` |
| Streaming FASTA/Sanger FASTQ summaries | `summarize_records(records, options)` |

The `StatisticsExt` trait and summary/configuration types are exported by
`bioferris::prelude::*`. Only the existing typed DNA, RNA and extended stop-free
protein alphabets implement `StatisticsAlphabet`; custom alphabets do not inherit
biological assumptions automatically. Defaults select `GcPolicy::Remove` and no
k-mer calculation. Inputs retain the existing ASCII, case-normalization and
validation rules in the [sequence contract](seq-bio-v1.md).

The manual workflow is a Rust program. The `stats` and `stats-file` CLI
operations and their JSON requests were removed with the `bioferris` CLI
(decision 029).

## Run a Rust analysis

Run from the repository root with its pinned Rust toolchain. The complete
program is checked in as
[`sequence_statistics.rs`](../../seq-bio/examples/sequence_statistics.rs).
Its literal inputs are DNA `acgnsw`, protein `MABZXUO`, an empty sequence iterator,
and this complete FASTA text embedded in the program:

```fasta
>short
G
>long
AAAA
```

The program demonstrates all three GC policies, both k-mer ambiguity policies,
a streaming summary and a handled resource-limit error.

<!-- statistics-program -->
```rust
use std::io::Cursor;

use bioferris::io::FastaReader;
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("acgnsw")?;
    let basic = sequence.statistics()?;
    println!("DNA: {sequence}");
    println!("Length: {}; ambiguous: {}", basic.length, basic.ambiguous);
    println!("Composition: {:?}", basic.composition);
    for (name, policy) in [
        ("remove", GcPolicy::Remove),
        ("ignore", GcPolicy::Ignore),
        ("weighted", GcPolicy::Weighted),
    ] {
        let summary = sequence.statistics_with(StatisticsOptions {
            gc_policy: policy,
            ..StatisticsOptions::default()
        })?;
        if let Some(gc) = summary.gc {
            println!(
                "GC ({name}): {:.6} / {} = {:.6}",
                gc.gc_equivalent, gc.denominator, gc.fraction
            );
        }
    }
    for (name, ambiguity) in [
        ("skip", KmerAmbiguity::Skip),
        ("literal", KmerAmbiguity::Literal),
    ] {
        let summary = sequence.statistics_with(StatisticsOptions {
            kmers: Some(KmerOptions {
                ambiguity,
                ..KmerOptions::new(2)
            }),
            ..StatisticsOptions::default()
        })?;
        if let Some(kmers) = summary.kmers {
            println!("2-mers ({name}): {:?}", kmers.counts);
            println!(
                "Windows: {}; counted: {}; skipped: {}",
                kmers.windows, kmers.counted, kmers.skipped
            );
        }
    }

    let fasta = ">short\nG\n>long\nAAAA\n";
    let records = FastaReader::<_, Dna>::new(Cursor::new(fasta));
    let summary = summarize_records(records, StatisticsOptions::default())?;
    println!(
        "Stream: {} records; {} symbols",
        summary.records, summary.length
    );
    if let Some(gc) = summary.gc {
        println!("Stream GC: {:.6}", gc.fraction);
    }
    let peptide = protein("MABZXUO")?.statistics()?;
    println!(
        "Protein: {} symbols; {} ambiguous; GC absent: {}",
        peptide.length,
        peptide.ambiguous,
        peptide.gc.is_none()
    );
    let empty = summarize_sequences(
        std::iter::empty::<&DnaSequence>(),
        StatisticsOptions::default(),
    )?;
    println!(
        "No records: {}; minimum: {:?}",
        empty.records, empty.minimum_length
    );

    let limited = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_windows: 1,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    if let Err(error) = sequence.statistics_with(limited) {
        println!("Window limit: {}", error.code());
    }
    Ok(())
}
```

Complete invocation:

<!-- statistics-command -->
```bash
cargo bio sequence_statistics
```

Complete stdout, exit 0, no stderr:

<!-- statistics-output -->
```text
DNA: ACGNSW
Length: 6; ambiguous: 3
Composition: {'A': 1, 'C': 1, 'G': 1, 'N': 1, 'S': 1, 'W': 1}
GC (remove): 3.000000 / 5 = 0.600000
GC (ignore): 3.000000 / 6 = 0.500000
GC (weighted): 3.500000 / 6 = 0.583333
2-mers (skip): {"AC": 1, "CG": 1}
Windows: 5; counted: 2; skipped: 3
2-mers (literal): {"AC": 1, "CG": 1, "GN": 1, "NS": 1, "SW": 1}
Windows: 5; counted: 5; skipped: 0
Stream: 2 records; 5 symbols
Stream GC: 0.200000
Protein: 7 symbols; 3 ambiguous; GC absent: true
No records: 0; minimum: None
Window limit: limit_exceeded
```

The three ambiguous DNA symbols are `N`, `S` and `W`, even though `S` and `W`
have definite GC status. Skipping ambiguous 2-mers therefore excludes `GN`, `NS`
and `SW`. The FASTA summary counts one GC-equivalent symbol among five symbols;
averaging the two record fractions would incorrectly give 0.5 instead of 0.2.
The final line is a handled error, so this demonstration exits successfully.

## Summary fields and aggregation

`sequence.statistics()` summarizes one sequence. `statistics_with(options)` adds
explicit policies. `summarize_sequences` accepts an iterator of borrowed
`&Sequence<A>` values. `summarize_records` accepts an iterator of
`Result<Record<A>, BioIoError>`, including the existing streaming readers. Both
aggregate without retaining a collection of input records. Summaries do not
modify source sequences, headers, quality scores or annotations.

| `SequenceStatistics` field | Meaning |
|---|---|
| `records: u64` | Number of input sequences, including empty sequences |
| `length: u64` | Sum of all sequence lengths |
| `minimum_length`, `maximum_length` | `Option<usize>` extrema; `None` for no records |
| `composition: BTreeMap<char, u64>` | Positive counts of observed uppercase symbols, sorted by symbol; absent symbols are omitted |
| `ambiguous: u64` | Count of ambiguous symbols under the molecule-specific rule below |
| `gc: Option<GcContent>` | Nucleotide GC data; `None` for protein |
| `kmers: Option<KmerCounts>` | Requested k-mer calculation, or `None` when disabled |

Zero input records have zero counts and no extrema. One empty sequence has
`records = 1` and both extrema `Some(0)`. Nucleotide summaries with a zero GC
denominator use fraction 0.0, not NaN; protein GC is absent even when empty.
If k-mers are requested for empty input, `kmers` remains `Some` with zero windows
and an empty key map. Summaries do not invent unseen symbols or missing qualities.

GC is calculated from aggregate composition after all input has succeeded.
Quality scores are validated by the input reader but do not weight composition,
GC or k-mer counts. The current API provides an aggregate, not a retained list of
per-record results. To report records separately, iterate the reader and call
`record.sequence().statistics_with(options)` for each successful record.

## GC and ambiguity policies

`GcContent` exposes the chosen `policy`, floating-point `gc_equivalent`, integer
`denominator`, and `fraction` in 0..=1. A percentage is `100 * fraction`;
the API does not multiply it automatically.

| `GcPolicy` | GC numerator | Denominator |
|---|---|---|
| `Remove` (default) | Count G, C and S | Count A, C, G, T/U, S and W |
| `Ignore` | Count G, C and S | All input symbols |
| `Weighted` | Equal-weight concrete IUPAC expansion | All input symbols |

For weighted GC, `N` contributes 1/2, `B` contributes 2/3 and `D` contributes 1/3.
These policies follow the [Biopython 1.88 GC contract](https://biopython.org/docs/1.88/api/Bio.SeqUtils.html#Bio.SeqUtils.gc_fraction)
within BioFerris's validated alphabets. No ambiguity is removed from the source,
composition or total length: only the GC numerator/denominator policy changes.

Ambiguity counting is deliberately separate. Every nucleotide symbol other than
A/C/G/T for DNA or A/C/G/U for RNA is ambiguous, including S/W. Protein B/Z/J/X
are ambiguous; U and O are concrete selenocysteine and pyrrolysine symbols.
Protein GC stays `None` under every `GcPolicy`.

Counts are checked `u64` integers. GC fractions and weighted contributions use
`f64`; normal-sized independent reference comparisons allow absolute error
1e-12. Integer-to-floating conversion can lose precision above 2^53. This is not
a claim of bit-for-bit agreement with every floating-point reference calculation.

## K-mer windows and resource budgets

Set `StatisticsOptions { kmers: Some(KmerOptions::new(k)), ..Default::default() }`.
Operations accept `1 <= k <= 64`; construction alone does not validate `k`.
Invalid options fail before the input iterator is consumed, even for empty input.
Use `options.validate()?` to perform the same configuration check before opening
an input stream. This checks option validity, not whether unseen input will
exceed a valid budget; zero budgets remain valid configurations.

All overlapping forward windows are considered within each record. A record of
length `n` contributes `max(0, n - k + 1)` windows. Windows never bridge two
records. There is no reverse-complement canonicalization or ambiguity expansion.
`KmerAmbiguity::Skip` omits any window containing a molecule-specific ambiguous
symbol. `Literal` keeps those uppercase letters as literal keys. U/O-containing
protein windows are not skipped merely because those residues are uncommon.

`KmerCounts` reports `k`, `ambiguity`, `windows`, `counted`, `skipped` and a
lexicographically sorted `BTreeMap<String, u64>` named `counts`. Successful
summaries satisfy `windows = counted + skipped`, and key counts sum to `counted`.
Repeated occurrences consume another window, but not another distinct-key slot.

| `KmerOptions` budget | Default | Scope |
|---|---:|---|
| `max_windows` | 1,000,000 | All considered windows, including skipped ones |
| `max_distinct` | 100,000 | Distinct keys across the entire summary |
| `max_key_bytes` | 8,388,608 | Stored key text across all distinct keys |

Budgets are inclusive and aggregate across records, not reset per record. Zero is
valid: a zero-window budget allows records shorter than `k`; zero key budgets can
still permit windows that are all skipped. Exceeding a budget is an error, not a
truncated successful count. Increase limits explicitly when the full analysis is
intended and sufficient resources are available.

The window budget is checked against the full next record's possible windows
before its symbols are scanned. A failure reports the attempted aggregate for
that record: for `AAAA`, `k = 2`, and `max_windows = 1`, `observed` is 3.
Distinct-key and key-byte budgets are checked before inserting each new key;
key bytes count text length, not the node overhead of the ordered map.

The implementation reuses an O(k) scratch buffer. Stored key data scale with
`distinct * k`; map-node overhead is additionally bounded by `max_distinct` but
is not included in `max_key_bytes`. Reader buffers and the current record are
separate allocations governed by the existing [file limits](file-io-v1.md#compression-streaming-and-limits).
This is not a fixed process-memory ceiling or a performance benchmark. Allocator
exhaustion is not converted into a recoverable biological error.

## Options

| Field | Values and rule | Default |
|---|---|---|
| `StatisticsOptions::gc_policy` | `GcPolicy::{Remove, Ignore, Weighted}`; no effect for proteins, whose GC is absent | `Remove` |
| `StatisticsOptions::kmers` | `Option<KmerOptions>`; `None` disables k-mers | `None` |
| `KmerOptions::k` | Forward window length; operations require 1..=64 | Set by the caller |
| `KmerOptions::ambiguity` | `KmerAmbiguity::{Skip, Literal}` | `Skip` |
| `KmerOptions::max_windows` | Nonnegative `u64`, including skipped windows | 1,000,000 |
| `KmerOptions::max_distinct` | Nonnegative `usize` | 100,000 |
| `KmerOptions::max_key_bytes` | Nonnegative `usize` | 8,388,608 |

K-mer budgets apply to the entire input, not each record. File-reader defaults
separately cap physical lines at 8 MiB, records at 100 million symbols and total
decompressed input at 1 GiB; configure readers explicitly when needed.

## Errors and streaming behavior

| `StatisticsError` variant | Stable code | Meaning |
|---|---|---|
| `InvalidOptions(String)` | `invalid_statistics_options` | Unsupported configuration, including `k` outside 1..=64 |
| `LimitExceeded { resource, limit, observed }` | `limit_exceeded` | An aggregate k-mer budget would be exceeded |
| `CountOverflow` | `count_overflow` | A count, total length or key-byte total cannot be represented |
| `Input(BioIoError)` | The reader's own code | Failed input, with source error and record/line context retained |

Limit resources are `kmer_windows`, `kmer_distinct` and `kmer_key_bytes`. Match
typed variants or `code()` instead of parsing human-readable messages.
The error in the complete example is `kmer_windows`, limit 1, observed 5.

An input failure or exceeded budget returns `Err`, never a successful partial
summary. Later records are not requested after that failure; buffered readers may
already have read ahead at the byte level. Plain and gzip FASTA/Sanger FASTQ reuse
the established readers and their validation, rather than bypassing malformed
records or inferring quality encodings. A damaged gzip trailer is an input error,
not a reason to accept only the prefix already decoded.

## Compatibility and verification

GC reference: the public [Biopython 1.88 `Bio.SeqUtils` documentation](https://biopython.org/docs/1.88/api/Bio.SeqUtils.html),
consulted 2026-10-15. Inputs remain BioFerris typed sequences: DNA/RNA reject
unsupported characters such as `X`, gaps and whitespace even when a reference
string-based utility may accept more symbols. K-mer ordering, forward orientation,
ambiguity policies and resource budgets are explicit BioFerris contracts, not a
claim of a matching Biopython k-mer API. No reference implementation source was
copied and no external-language runtime is required.

The [independent Rust suite](../../seq-bio/tests/statistics.rs) checks every
three-symbol IUPAC DNA/RNA sequence under every GC policy using concrete-base
sets, plus public GC examples. Generated DNA/RNA/protein record collections are
compared with naive composition/window oracles, including aggregation order and
record boundaries. Deterministic cases cover each budget, skipped windows,
empty input, invalid `k`, streaming failures, gzip/FASTQ and the existing public
phiX174 fixture. Separate core tests exercise arithmetic-overflow guards without
allocating impossible data volumes. The runnable example's complete source and
stdout are checked against this page by an integration test.

These checks establish the documented subset. They are not a live Biopython
differential campaign, a measured throughput claim, or new qualification evidence
for the unchanged alignment algorithms.
Commands, results and retained logs appear in the
stage 2 validation report.
