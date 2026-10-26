# BLAST result input: Rust workflows

## Contents

- [Scope and integration status](#scope-and-integration-status)
- [Write a Rust analysis](#write-a-rust-analysis)
- [Stream large XML reports](#stream-large-xml-reports)
- [Inspect XML reports](#inspect-xml-reports)
- [Errors and recovery](#errors-and-recovery)
- [Compatibility and limits](#compatibility-and-limits)

[API index](../README.md) · [Input contract](../algorithms/blast-input-v1.md)
· [Fixture and reference provenance](../../seq-bio/tests/fixtures/blast-provenance.md)

## Scope and integration status

The Rust API and Cargo example are implemented and verified on Rust 1.85.0.
Stage 5 is accepted for the explicit tabular/XML profiles in the
main validation record. That
record also covers CLI operations, strict JSON requests and discovery, which were
removed with the `bioferris` CLI (decision 029). The complete
program and real Cargo alias below pass the tests. Earlier isolated/staging
results remain separate checkpoints; this is not a package release or universal
SearchIO/API-equivalence claim.

Import reported BLAST tabular outfmt 6/7, legacy XML outfmt 5, or single-file XML2
outfmt 16. This does not run BLAST, contact a search server or reconstruct a
missing alignment. Choose the format explicitly; no molecule/program guessing.

| Source profile | Rust entry points and format option |
|---|---|
| Outfmt 6 rows | `parse_blast_tabular` / `read_blast_tabular`, `BlastTabularFormat::Uncommented` |
| Outfmt 7 comments/rows | Same helpers, `BlastTabularFormat::Commented` |
| Legacy XML | `parse_blast_xml` / `read_blast_xml`, `BlastXmlFormat::Legacy` |
| Single-file namespace-aware XML2 | Same helpers, `BlastXmlFormat::Xml2` |

All APIs are in `bioferris::io::blast`. `parse_*` accepts a buffered reader;
`read_*` accepts a named plain/gzip path, detected by content signature. Both
keep every query until the document ends. The streaming iterators
`BlastTabularReader` and `BlastXmlReader` yield one complete query at a time,
for [large reports](#stream-large-xml-reports); consume them to EOF to validate
the remaining source and compressed trailers.

## Write a Rust analysis

Work from the repository root. The complete runnable program is
[`seq-bio/examples/blast_analysis.rs`](../../seq-bio/examples/blast_analysis.rs).
Its literal tabular input is one standard twelve-column row, with tabs and a
final newline. The subject endpoints are intentionally reversed:

<!-- blast-input -->
```text
q	s	100	4	0	0	1	4	8	5	1e-20	42
```

The program also reads three versioned fixtures from the repository root. The
commented outfmt 7 file
[`blast-example.tab7`](../../seq-bio/tests/fixtures/blast-example.tab7) holds
one query with a hit and one without:

<!-- blast-commented-input -->
```text
# BLASTN 2.15.0+
# Query: q example query
# Database: synthetic
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 1 hits found
q	s	100	4	0	0	1	4	8	5	1e-20	42
# BLASTN 2.15.0+
# Query: absent no matching subject
# Database: synthetic
# 0 hits found
# BLAST processed 2 queries
```

The synthetic XML reports
[`blast-example.xml`](../../seq-bio/tests/fixtures/blast-example.xml) and
[`blast-example-xml2.xml`](../../seq-bio/tests/fixtures/blast-example-xml2.xml)
each contain the two explicitly named no-hit queries `one` and `two`, with
lengths 4 and 5. The legacy file has one report shared by both queries; XML2 has
two distinct reports. These are format demonstrations, not actual search runs.

Complete program, including an omitted-score and a resource-limit boundary:

<!-- blast-rust-program -->
```rust
//! Read reported BLAST coordinates without guessing missing biological values.

use std::io::Cursor;

use bioferris::io::blast::{
    BlastTabularFormat, BlastTabularOptions, BlastXmlFormat, BlastXmlOptions, parse_blast_tabular,
    read_blast_tabular, read_blast_xml,
};

const INPUT: &str = "q\ts\t100\t4\t0\t0\t1\t4\t8\t5\t1e-20\t42\n";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = parse_blast_tabular(Cursor::new(INPUT), BlastTabularOptions::default())?;
    let query = &document.queries()[0];
    let hit = &query.hits()[0];
    let hsp = &hit.hsps()[0];
    let query_coordinates = hsp.query_coordinates().ok_or("missing query coordinates")?;
    let hit_coordinates = hsp.hit_coordinates().ok_or("missing subject coordinates")?;

    println!("Query: {} ({} hit)", query.id(), query.hits().len());
    println!("Query range: {:?}", query_coordinates.range());
    println!(
        "Subject: {} {:?} (reported {}..{})",
        hit.id(),
        hit_coordinates.range(),
        hit_coordinates.reported_start(),
        hit_coordinates.reported_end()
    );
    println!("Reported strand: {:?}", hsp.hit_strand());
    println!(
        "E-value source: {}",
        hsp.field("evalue").ok_or("missing evalue")?
    );

    let minimal = parse_blast_tabular(
        Cursor::new("q\ts\n"),
        BlastTabularOptions {
            fields: Some(vec!["qseqid".into(), "sseqid".into()]),
            ..BlastTabularOptions::default()
        },
    )?;
    println!(
        "Omitted e-value: {:?}",
        minimal.queries()[0].hits()[0].hsps()[0].evalue()
    );

    let mut limited = BlastTabularOptions::default();
    limited.limits.max_queries = 0;
    if let Err(error) = parse_blast_tabular(Cursor::new(INPUT), limited) {
        println!("Rejected zero query budget: {}", error.code());
    }

    let commented = read_blast_tabular(
        "seq-bio/tests/fixtures/blast-example.tab7",
        BlastTabularOptions {
            format: BlastTabularFormat::Commented,
            ..BlastTabularOptions::default()
        },
    )?;
    for query in commented.queries() {
        println!(
            "Commented query {} ({}): {} hit(s)",
            query.id(),
            query.description().unwrap_or("no description"),
            query.hits().len()
        );
    }
    println!(
        "Declared processed queries: {:?}",
        commented.processed_queries()
    );

    for (format, path) in [
        (
            BlastXmlFormat::Legacy,
            "seq-bio/tests/fixtures/blast-example.xml",
        ),
        (
            BlastXmlFormat::Xml2,
            "seq-bio/tests/fixtures/blast-example-xml2.xml",
        ),
    ] {
        let options = BlastXmlOptions {
            format,
            ..BlastXmlOptions::default()
        };
        for query in read_blast_xml(path, options)?.queries() {
            println!(
                "{format:?} XML query {}: length {:?}, {} hit(s)",
                query.id(),
                query.length(),
                query.hits().len()
            );
        }
    }
    Ok(())
}
```

<!-- blast-rust-command -->
```bash
cargo bio blast_analysis
```

Complete observed stdout, empty stderr, exit 0:

<!-- blast-rust-output -->
```text
Query: q (1 hit)
Query range: 0..4
Subject: s 4..8 (reported 8..5)
Reported strand: None
E-value source: 1e-20
Omitted e-value: None
Rejected zero query budget: limit_exceeded
Commented query q (example query): 1 hit(s)
Commented query absent (no matching subject): 0 hit(s)
Declared processed queries: Some(2)
Legacy XML query one: length Some(4), 0 hit(s)
Legacy XML query two: length Some(5), 0 hit(s)
Xml2 XML query one: length Some(4), 0 hit(s)
Xml2 XML query two: length Some(5), 0 hit(s)
```

Coordinates are ascending zero-based half-open ranges. Original positive one-based
endpoints remain separately available; decreasing endpoints alone do not establish
a biological strand. The lexical e-value is retained as `1e-20`. Two-column
input has no reported score, so `None` is preserved instead of manufacturing zero.
The parser returns a typed limit error rather than a silently truncated document.
The commented file keeps its no-hit query and the declared count of processed
queries, and both XML profiles yield the same queries and lengths.

## Stream large XML reports

`parse_blast_xml` and `read_blast_xml` keep every query until the document ends.
On the publication benchmark's largest legacy report, 150 MB holding 1,000
queries with up to 100 hits each, that peaked at 535 MiB. `BlastXmlReader`
yields one complete query at a time instead; on the same report it peaked at
8 MiB, and at 7 MiB on the 175 MB XML2 version, with the same queries, hits and
HSPs (validation record).
Use it for large reports, and the collecting helpers for reports that fit
comfortably in memory.

The complete program
[`seq-bio/examples/blast_stream.rs`](../../seq-bio/examples/blast_stream.rs)
streams four versioned reports: the synthetic legacy example above, and three
unmodified NCBI BLAST 2.9.0+ reports described in the
[fixture provenance](../../seq-bio/tests/fixtures/blast-provenance.md), a blastn
search in legacy XML and XML2 and a blastx search in legacy XML.
`bioferris::io::open_input` opens a plain or gzip file, recognized by its content.

<!-- blast-stream-program -->
```rust
//! Stream BLAST XML reports one complete query at a time, as for files too large
//! to collect into one document.

use bioferris::io::blast::{BlastXmlFormat, BlastXmlOptions, BlastXmlReader};
use bioferris::io::open_input;

const REPORTS: [(BlastXmlFormat, &str); 4] = [
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-example.xml",
    ),
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-nucleotide.xml",
    ),
    (
        BlastXmlFormat::Xml2,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-nucleotide-xml2.xml",
    ),
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-translated.xml",
    ),
];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for (format, path) in REPORTS {
        let name = path.rsplit('/').next().unwrap_or(path);
        let options = BlastXmlOptions {
            format,
            ..BlastXmlOptions::default()
        };
        // Each item is one complete query; the reader keeps no earlier ones.
        for query in BlastXmlReader::new(open_input(path)?, options)? {
            let query = query?;
            let hsps: usize = query.hits().iter().map(|hit| hit.hsps().len()).sum();
            println!(
                "{name}: {} with {} hit(s), {hsps} HSP(s); first hit {}",
                query.id(),
                query.hits().len(),
                query.hits().first().map_or("none", |hit| hit.id())
            );
        }
    }

    // A query larger than its source budget fails instead of being buffered.
    let mut options = BlastXmlOptions::default();
    options.limits.max_query_bytes = 1024;
    let mut reader = BlastXmlReader::new(open_input(REPORTS[1].1)?, options)?;
    if let Some(Err(error)) = reader.next() {
        println!("Query over a 1 KiB budget: {}: {error}", error.code());
    }
    Ok(())
}
```

<!-- blast-stream-command -->
```bash
cargo bio blast_stream
```

Complete observed stdout, empty stderr, exit 0:

<!-- blast-stream-output -->
```text
blast-example.xml: one with 0 hit(s), 0 HSP(s); first hit none
blast-example.xml: two with 0 hit(s), 0 HSP(s); first hit none
blast-ncbi-2.9.0-nucleotide.xml: G26684.1 with 10 hit(s), 14 HSP(s); first hit gi|372099107|ref|NC_000069.6|
blast-ncbi-2.9.0-nucleotide-xml2.xml: G26684.1 with 10 hit(s), 14 HSP(s); first hit gi|372099107|ref|NC_000069.6|
blast-ncbi-2.9.0-translated.xml: AI021773.1 with 10 hit(s), 10 HSP(s); first hit gi|1530504495|emb|VDM03167.1|
Query over a 1 KiB budget: limit_exceeded: BLAST XML: I/O error: BLAST XML source/event/capture byte limit exceeded (record 1, line 49)
```

Each item is one complete query in source order; once the loop moves on, that
query is dropped, so memory holds one query and the reader's own state. Both
encodings of the blastn search yield the same query and hits. The reader's
state is bounded event buffers, document and report metadata, and the IDs of
the completed queries, which it keeps to reject a repeated query ID.

Budgets that matter for large reports, all inclusive, in `BlastXmlOptions`:

| Field | Default | Bounds |
|---|---|---|
| `limits.read.max_total_bytes` | 1 GiB | The whole decoded source |
| `limits.max_queries` | 100,000 | Queries in the document |
| `limits.max_query_bytes` | 64 MiB | One query's source, including its tags |
| `limits.max_query_id_bytes` | 16 MiB | The IDs of completed queries that the reader keeps |
| `limits.max_hits_per_query`, `limits.max_hsps_per_query` | 100,000 and 1,000,000 | Hits and HSPs of one query |
| `max_metadata_bytes` | 16 MiB | Document and report metadata |

Raise the ones a trusted large report needs; they bound source bytes and
counts, not memory. The last line of the output shows a query larger than its
budget: the 2.9.0 blastn query exceeds 1 KiB, so the reader returns
`limit_exceeded` for it and then ends, without a partial query. Consume the
iterator to its end even after the queries of interest: only then are the
closing markup and a gzip trailer verified, and a query already yielded does
not prove that the rest of the file is valid.

## Inspect XML reports

`BlastDocument::source_fields()` and each query's `report_fields()`,
`program()`, `version()` and `target()` retain the report headers, parameters
and query metadata that distinguish the two files. Ordinary `search/Search` XML2
results are supported; iterative/PSI/DELTA and bl2seq result modes fail
explicitly. Known external DTD text never causes a network request; declared
entities are refused.

## Errors and recovery

Errors are `BioIoError` values with stable `code()` identifiers. In the program
above, a zero query budget rejects the first query with `limit_exceeded`.
No partial result is returned. If the larger trusted input is intended,
raise the explicit budget and retry; do not treat a failure as a valid empty
search. Similarly, missing paths/decompression errors use `io_error`, truncated
XML uses `truncated_input`, inconsistent dimensions use `length_mismatch`, and
unsupported semantics use `unsupported_format_feature`. Codes, not message
wording, are the control-flow contract.

## Compatibility and limits

The patched XML engine is pinned to `quick-xml` 0.41.0 on Rust 1.85.0. Each
element has a fixed maximum of 256 namespace declarations; exceeding it is an
`InvalidStructure` error (code `invalid_structure`). Existing
decoded/unescaped attribute whitespace remains preserved. The
[security verification](../testing/security.md) explains the upstream advisories
and boundary regressions; this is not blanket XML/XSD support.

The [technical contract](../algorithms/blast-input-v1.md) defines all defaults,
column mappings, identity selection, XML profiles, metadata preservation, strict
options and error categories. Defaults cap decoded bytes, physical lines,
queries, per-query hits/HSPs/source bytes, completed-ID bytes and XML
event/depth/report metadata. These are source/count limits, not exact memory
measurements. The readers never write files.

Retained Biopython 1.88 public API observations cover eight inputs, eleven queries
and 52 HSPs on CPython 3.12.3/NumPy 2.3.3. All compared explicit source values agree
exactly in that corpus. Unlike NCBIXML, the importer does not synthesize omitted
XML2 positive counts or missing frames from program context. Missing values remain
absent; original IDs are not heuristically rewritten. This is not a claim of
universal SearchIO/API/format equivalence, XML lexical round-trip, experimental
validity or BLAST-engine fidelity. The source data and wheel/fixture SHA-256
values are recorded in the
[provenance](../../seq-bio/tests/fixtures/blast-provenance.md).

The tests in `seq-bio/tests/blast_recipes.rs` compare the complete literal
programs and input fences with their files and run the real
`cargo bio blast_analysis` and `cargo bio blast_stream` aliases against the
outputs above. The earlier staging
run is retained separately. The
main validation record records
accepted Stage 5 scope, toolkit/optimized checks, complete workspace harness
output and terminal audit, documentation and response comparisons. Final
source/evidence review and the local commit accompany this acceptance. Acceptance
does not expand the declared reference compatibility scope.
