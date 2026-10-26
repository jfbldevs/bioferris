# BLAST result input contract

## Contents

- [Scope and delivery state](#scope-and-delivery-state)
- [Shared result model](#shared-result-model)
- [Tabular columns and validation](#tabular-columns-and-validation)
- [Commented query blocks](#commented-query-blocks)
- [Streaming and resource limits](#streaming-and-resource-limits)
- [XML profiles and safety](#xml-profiles-and-safety)
- [Rust interfaces](#rust-interfaces)
- [Verified input/output examples](#verified-inputoutput-examples)
- [Errors and compatibility boundaries](#errors-and-compatibility-boundaries)
- [References and compatibility evidence](#references-and-compatibility-evidence)

[API index](../README.md) · Implementation plan
· [Manual Rust guide](../api/blast-input-v1.md)

## Scope and delivery state

Stage 5 implements tabular BLAST outfmt 6/7, legacy outfmt 5 XML and single-file
outfmt 16 XML2 for the explicit profiles below. The core passed isolated Rust
1.85.0 checks, including retained public-reference observations. **Stage 5 is
accepted for the documented tabular/XML profiles.** The
main validation record records
final checks and scope; the CLI adapter, strict requests and discovery it also
covered were removed with the `bioferris` CLI (decision 029). The earlier isolated observations below remain
their own checkpoint, separate from main-checkout acceptance.
Reading a result never executes a BLAST search.

The functional reference is the public Biopython 1.88 SearchIO API. Primary NCBI
format documentation defines source fields and XML structure. The parser does not
copy reference implementation source, normalize identifiers heuristically, discard
unknown columns or synthesize absent scores.

## Shared result model

`BlastDocument` contains ordered queries, document comments and an optional source
processed-query count. `BlastQuery` contains a chosen identifier with its source
field name, an optional separately preserved comment-header identifier,
description, sequence length, program/version/database, ordered field names,
comments, source metadata and hits. `BlastHit` contains its chosen identifier and
source field, optional length/description, source metadata and ordered HSPs.
`BlastHsp` preserves all source column names/values in order and exposes optional
typed scores, counts, percentages, aligned strings, frames and coordinates.
XML queries additionally share immutable report metadata through an `Arc`-backed
slice returned by `report_fields()`. Shared report metadata is not duplicated
into document-level source fields. A legacy report containing no queries retains
its otherwise orphaned metadata at document level.

Coordinates are ascending zero-based half-open ranges: source endpoints `a,b`
become `[min(a,b)-1,max(a,b))`. Both original positive one-based endpoints remain
available. `Increasing`, `Decreasing` and `Point` describe source coordinate
order; they do not claim a protein strand. A biological strand is present only
when explicitly supplied or implied by a nonzero frame. A zero frame is retained
as zero, not converted to a nucleotide strand. Translated coordinates remain in
the source sequence's units; amino acid alignment strings are not used to divide
nucleotide positions by three.

Missing optional values remain absent. Numeric score getters use finite binary64
values; original source text remains available even when decimal precision is
greater than binary64 or a very small positive e-value underflows to zero. Integer
coordinates/counts use checked unsigned 64-bit values; frames are in -3..3.
No identity percentage, mismatch count, strand or score is invented from other
fields. Query/hit metadata must agree across their HSP rows.

## Tabular columns and validation

`BlastTabularOptions` explicitly selects uncommented outfmt 6 or commented outfmt
7. The default is outfmt 6 with these twelve columns:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

An explicit field list may reorder, omit or extend columns. `std` expands to the
default twelve; duplicate names are errors. Field names use ASCII letters, digits
and underscores, beginning with a letter, with a 256-byte name limit. At least one query identity field
(`qseqid`, `qaccver`, `qacc`, `qgi`) and one subject identity field (`sseqid`,
`saccver`, `sacc`, `sgi`) are required. That order also defines identity selection
when several alternatives are supplied. All original fields remain available.
Coordinate pairs must be requested together; a frame or aligned string may be
requested independently. Unknown custom fields remain raw values with no inferred
type. `sallseqid` alone cannot supply a unique subject identity.

The default delimiter is a tab. A caller may select a non-control printable ASCII
byte, except `#`, to match NCBI `delim` output. Delimiters are not quoted or escaped;
each row must have exactly the selected number of fields. LF/CRLF and a final
line without a newline are supported. Empty physical lines are tolerated.
Nonempty rows with control bytes, invalid UTF-8, empty identities, invalid numeric
values, zero coordinates, incomplete coordinate pairs or contradictory fields
fail with line/query context.

Known numeric fields are validated even when only their source text is exposed.
Percentages must lie in 0..100; counts are nonnegative; alignment/sequence lengths
are positive; e-values and bit scores are finite and nonnegative. Raw scores are
finite. Coordinate bounds are checked against `qlen`/`slen` when present. Explicit
frame signs and subject strands must agree with non-point coordinate direction
and with each other. `frames` must contain exactly two valid frames separated by
`/`, agreeing with separately supplied frame columns.

When `length` is present, individual identity/positive/mismatch/gap counts cannot
exceed it; identities cannot exceed positives. When identities, mismatches and
gaps are all supplied, their sum must equal the reported length. Supplied alignment
strings use ASCII letters, `*` and `-`, preserve case, have equal column lengths
when both are present, and agree with explicit `length`. Rounded percentages are
not forced to equal ratios recomputed from counts. BTOP and taxonomy/title/list
fields are preserved without pretending to reconstruct an alignment or taxonomy.

## Commented query blocks

Outfmt 7 preserves every comment exactly except its physical line ending,
including unknown comments and the final processed-query summary. Header fields
recognize program/version, `Query`, `Database`, `RID`, `Fields` and the declared
number of result rows. No-hit query blocks produce a query with zero hits; absent
scores remain absent. A no-hit query's identifier comes from the source query
header and is labeled accordingly.

`# Fields:` labels select their canonical NCBI column names. Unknown labels need
an explicit caller field mapping with matching column count. Recognized labels
must agree with the same-position explicit name. This check intentionally differs
from Biopython's documented policy of ignoring a supplied field list in commented
mode. Raw labels are retained in comments. Fields and row-count declarations
cannot be repeated or occur after data rows. A data row requires a query header,
field declaration and row-count declaration. Counts describe HSP rows, not distinct
subject IDs. Completed blocks must contain exactly the declared number of rows.
The optional final processed-query summary must agree with the number of blocks.

Rows within a query may revisit a subject; HSPs merge into that subject in source
order and hits retain first-occurrence order. Outfmt 6 queries must be contiguous.
A completed query identifier cannot reappear later; that is an explicit error,
not silent reordering or a second incomplete result. Commented blocks likewise
have unique chosen identifiers. A query header ID can differ from an accession
column; both are retained without reconstructing one from the other.

## Streaming and resource limits

`BlastTabularReader<R: BufRead>` yields complete queries and fuses after the first
error. It retains one current query, one boundary lookahead line, and a bounded
set of completed query IDs. Query comments stay with their query; comments outside
query blocks and a footer remain accessible on the reader. Consume to EOF to
validate the footer, later query groups, decompression trailers and all input.
Previously yielded queries do not imply an unread suffix is valid.

The eager `parse_blast_tabular`/`read_blast_tabular` helpers collect a document,
including otherwise unassigned comments. Named paths reuse content-signature
plain/gzip input. Parser options reuse `ReadLimits`: 8 MiB per physical line and
1 GiB total decoded bytes by default. BLAST adds defaults of 100,000 queries,
100,000 hits per query, 1,000,000 HSPs per query, 64 MiB source bytes per query,
16 MiB retained completed-query identifier bytes and at most 128 columns. Zero
budgets permit only inputs requiring no resource of that kind. These are source,
entry and byte budgets, not exact process-memory limits.

## XML profiles and safety

XML uses the established safe Rust `quick-xml` event parser, the same result model,
and two explicit profiles:

1. Legacy outfmt 5 follows NCBI `NCBI_BlastOutput.mod.dtd`: `BlastOutput` report
   metadata, `BlastOutput_iterations/Iteration`, `Iteration_hits/Hit`, `Hit_hsps/Hsp`,
   parameters, statistics, messages and optional query-level fields.
2. Single-file outfmt 16 follows NCBI `schema_alt/NCBI_BlastOutput2.mod.xsd` with
   namespace `http://www.ncbi.nlm.nih.gov`. `BlastXML2` contains `BlastOutput2`
   reports; each `report/Report` contains `program`, `version`, `reference`,
   `search-target/Target`, `params/Parameters`, and `results/Results`. `Results`
   distinguishes `search/Search`, `iterations/Iteration` and `bl2seq/Search`.
   The implemented profile accepts `search/Search` only, with exactly one Search
   per report and any number of reports. Iterative and bl2seq modes are rejected;
   rounds are never collapsed. The separate `schema/` underscore-tag variant is
   not the same profile.

XML2 `Hit/description/HitDescr` may contain multiple identifiers, accessions,
titles and taxonomy descriptions. Preserve their order instead of merging them
into an invented single title. HSP fields include original query/hit endpoints,
explicit strands/frames, scores, identity/positive/gap counts, aligned strings and
midline; report parameters/statistics and optional messages remain represented.
Generated IDs remain raw unless an explicit separately verified identity policy
is added. XML2 error reports must become explicit errors, not empty successes.

The XML dependency is pinned `quick-xml = 0.41.0`, default features disabled:
its official manifest declares Rust 1.79, MIT, and `memchr` as its sole mandatory
dependency. Its streaming event API and namespace resolution are used without an
application-defined XML tokenizer. No remote DTD/entity resolution, `XInclude`,
network access or processing instruction execution occurs. XML 1.0 UTF-8 and
US-ASCII declarations are supported; literal non-ASCII bytes contradicting an
ASCII declaration fail, while legal numeric character references remain valid.
Unknown/duplicate/out-of-order declaration attributes fail. Text is checked for
legal XML 1.0 characters; predefined entities and legal numeric references are
decoded once. XML2 uses exactly namespace `http://www.ncbi.nlm.nih.gov`; legacy
XML uses no namespace. Namespace declarations and the root XML2 schema-location
attribute are retained. Other attributes and unknown semantic elements fail.
Each element may declare at most 256 namespace bindings; excess declarations fail
with `InvalidStructure` before unbounded resolver allocation. This explicit fixed
safety limit and the patched attribute duplicate checks address
[RUSTSEC-2026-0194](https://rustsec.org/advisories/RUSTSEC-2026-0194.html) and
[RUSTSEC-2026-0195](https://rustsec.org/advisories/RUSTSEC-2026-0195.html).
Retained attribute values remain decoded/unescaped with source whitespace
preserved; the dependency migration does not silently introduce additional
attribute normalization. See the [security verification](../testing/security.md).
The known legacy external DOCTYPE is retained without opening the named resource;
internal entity declarations are refused. No full DTD/XSD validation or guarantee
of original lexical XML serialization is implied.

`BlastXmlReader<R: BufRead>` yields one complete query at a time and fuses after
the first error. It retains a bounded event buffer, query/metadata subtree,
namespace/depth state and completed IDs, not an entire unbounded XML tree.
`parse_blast_xml` and `read_blast_xml` collect only after consuming the complete
document and transport. Supported program metadata names are blastn, blastp,
blastx, tblastn and tblastx, case-insensitively checked but preserved verbatim.
XML2 requires an explicit query-id and at least one HitDescr per hit; a new
identity policy is required for documents without these. Legacy optional query
ID/description/length may fall back to explicitly retained report-level values.
Repeated chosen query IDs are rejected as unsupported noniterative semantics.

Shared source budgets apply to XML too. Extra XML defaults are depth 64, 8 MiB
per source event and 16 MiB total document/report metadata source bytes. Event
and capture limits apply before unrestricted buffering. Query-source limits
include the query opening/closing tags. These are source/count budgets, not an
exact heap-size guarantee. Consuming to EOF is necessary to validate trailing
markup and gzip integrity; an earlier yielded query is not proof of a valid suffix.

XML-specific isolated acceptance tests cover real public v1/v2 data, multiple/no-hit queries,
reverse/translated coordinates, repeated descriptions, absent optional fields,
malformed nesting/attributes/entities, unsupported namespaces/encodings/rounds,
truncation, bounded events/depth and transport errors. The unmodified public
reports and retained observations are described in the
[fixture provenance](../../seq-bio/tests/fixtures/blast-provenance.md).

## Rust interfaces

The public module is `bioferris::io::blast`. Tabular streaming uses
`BlastTabularReader::new(reader, BlastTabularOptions)`; XML streaming uses
`BlastXmlReader::new(reader, BlastXmlOptions)`. Their iterator items are
`Result<BlastQuery, BioIoError>`. Both expose accumulated document comments;
XML exposes document source fields and tabular exposes the processed-query footer.
The eager helpers are `parse_blast_tabular`, `read_blast_tabular`,
`parse_blast_xml` and `read_blast_xml`. Named readers use content-signature
plain/gzip detection, including concatenated gzip members. Inputs are immutable.

## Verified input/output examples

These examples were first executed in the isolated harness with Rust 1.85.0.
Main-checkout guide programs and final checks support the scoped acceptance
recorded in the main validation record.

### Manual Rust program

Literal input has two columns separated by a tab and one final newline:

```text
q	s
```

Complete Rust program:

```rust
use std::io::Cursor;
use bioferris::io::blast::{BlastTabularOptions, parse_blast_tabular};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = "q\ts\n";
    let document = parse_blast_tabular(
        Cursor::new(input),
        BlastTabularOptions {
            fields: Some(vec!["qseqid".into(), "sseqid".into()]),
            ..BlastTabularOptions::default()
        },
    )?;
    let query = &document.queries()[0];
    let hit = &query.hits()[0];
    let hsp = &hit.hsps()[0];
    println!("{}: {} hit", query.id(), query.hits().len());
    println!("{}: evalue={:?}, coordinates={:?}", hit.id(), hsp.evalue(), hsp.query_coordinates());
    Ok(())
}
```

The minimal program above can be run with `cargo run --quiet` in an external
Cargo project depending on BioFerris. The repository's separate, richer
[`blast_analysis.rs`](../../seq-bio/examples/blast_analysis.rs) example runs from
the repository root with `cargo bio blast_analysis`; its different complete
input/output is documented in the [manual guide](../api/blast-input-v1.md).
The observed output below belongs to the complete minimal program above:

```text
q: 1 hit
s: evalue=None, coordinates=None
```

No score or coordinate was reported, so the parser does not invent either.

### Commented result with a hit and a no-hit query

The complete literal fixture is
[`blast-example.tab7`](../../seq-bio/tests/fixtures/blast-example.tab7):

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

The [manual guide's program](../api/blast-input-v1.md#write-a-rust-analysis)
reads this file with `BlastTabularFormat::Commented`; its complete output
includes `Commented query q (example query): 1 hit(s)`,
`Commented query absent (no matching subject): 0 hit(s)` and
`Declared processed queries: Some(2)`. The source subject endpoints 8..5 become
4..8. Direction is decreasing; no biological strand is inferred from tabular
endpoints alone. The no-hit query is retained rather than disappearing.

The same program shows a resource-limit failure: a zero query budget fails with
`limit_exceeded` and returns no partial query list. Raise the explicit budget
only when the larger trusted input is intended; the parser never silently
truncates to zero.

## Errors and compatibility boundaries

The shared typed `BioIoError` carries a stable category, format, message and
one-based record/line where known. Categories include `io_error` for open,
decompression and transport failures; `invalid_header`, `invalid_structure`,
`invalid_sequence`, `length_mismatch`, `truncated_input`,
`unsupported_format_feature` and `limit_exceeded`. Messages are explanatory;
control flow must use codes. No partial eager document is returned on failure.

Observed compatibility is scoped to eight retained inputs, eleven queries and
52 HSPs under Biopython 1.88/CPython 3.12.3/NumPy 2.3.3 on the observed platform.
All compared explicit source coordinates, scores, supplied counts, alignments,
midlines, identities and ordering agree. Numeric comparisons are exact for this
retained corpus; no tolerance or broader fidelity claim is inferred.

Intentional differences and unsupported semantics remain explicit:

- BioFerris keeps generated/raw IDs and source numeric spellings; it does not
  promise Biopython object identity, ID heuristics, warning classes or formatting.
- NCBIXML may synthesize a missing XML2 positive count from identities and absent
  frames from program context. BioFerris leaves omitted counts/frames absent;
  tests verify both the retained observation and this source-preservation policy.
- Unknown custom tabular columns are retained. Unknown XML semantic elements
  fail; broader XML variants cannot silently lose information.
- Commented-mode explicit field mappings are checked rather than ignored.
- Iterative/PSI/DELTA, XML2 bl2seq, external entities, XInclude, arbitrary encodings,
  result writing/indexing, live BLAST execution and universal SearchIO parity
  are not implemented by this importer.

Isolated verification checkpoint (2026-10-18): 37 core tests and 14 adapter tests
pass in debug and release with Rust 1.85.0; strict all-target Clippy including
pedantic, owned-source formatting and strict library rustdoc generation pass.
The complete Rust programs in this contract and the separate API guide each
pass rustdoc execution. Five adapter success/error output pairs are byte-identical
between debug/release. All twelve retained data/request/observation SHA-256 values
match the provenance record. This paragraph records the isolated checkpoint,
not a new main-checkout test result. The
[promoted recipe tests](../../seq-bio/tests/blast_recipes.rs), the CLI tests
(since removed with the CLI) and
[`blast_analysis.rs`](../../seq-bio/examples/blast_analysis.rs) existed alongside
the implemented public wiring. Their main-checkout recipes/Cargo alias,
documentation checks and complete workspace harness output/terminal audit now
support scoped acceptance in the
main validation record.
That main evidence is distinct from the isolated checkpoint. Final source/evidence
review and the local commit accompany this acceptance; no closed scientific campaign
is reopened or inferred from these parser checks.

## References and compatibility evidence

- [NCBI BLAST command options and default tabular fields](https://www.ncbi.nlm.nih.gov/books/NBK279684/).
- [NCBI custom output examples](https://www.ncbi.nlm.nih.gov/books/NBK569862/).
- [Biopython 1.88 BLAST SearchIO public API](https://biopython.org/docs/1.88/api/Bio.SearchIO.BlastIO.html).
- [Biopython 1.88 BLAST XML public API](https://biopython.org/docs/1.88/api/Bio.Blast.NCBIXML.html).
- [Biopython 1.88 SearchIO coordinate conventions](https://biopython.org/docs/1.88/Tutorial/chapter_searchio.html).
- [NCBI legacy XML DTD](https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd).
- [NCBI XML2 schema_alt specification](https://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.mod.xsd).
- [quick-xml 0.41.0 official manifest](https://github.com/tafia/quick-xml/blob/v0.41.0/Cargo.toml).

Scientific examples and independent parser tests establish only their declared
input profile. Exact lexical preservation, normalized semantic agreement and
Biopython object/API equivalence are distinct claims. No live BLAST execution,
universal SearchIO parity or performance result is implied by the isolated
checks or accepted main-checkout scope in this contract.
