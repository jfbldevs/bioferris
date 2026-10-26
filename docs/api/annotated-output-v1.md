# Checked annotated output v1

## Capability index

| Workflow | Rust interface in `io::annotated` |
|---|---|
| Same-profile stream output | `write_genbank`, `write_embl` |
| Create-only file output | `write_genbank_new`, `write_embl_new` |
| Strict annotated slice | `AnnotatedRecord::slice` |
| Annotated DNA/RNA reverse complement | `AnnotatedRecord::reverse_complement` |
| Inspect changes and losses | `AnnotatedTransform::report`, a typed report |

## Contents

- [Run an analysis](#run-an-analysis)
- [Writing contract](#writing-contract)
- [Transformation options](#transformation-options)
- [Locations and biological invariants](#locations-and-biological-invariants)
- [Metadata and qualifiers](#metadata-and-qualifiers)
- [Source classification](#source-classification)
- [Reports and agent use](#reports-and-agent-use)
- [Errors and limits](#errors-and-limits)
- [Compatibility and verification](#compatibility-and-verification)

[API index](../README.md#api-index) · [Annotated input](annotated-input-v1.md) ·
[Complete inputs and outputs](../examples/annotated-output-recipes.md)

## Run an analysis

Write an ordinary Rust program importing `bioferris`. The complete checked-in
[annotated_records.rs](../../seq-bio/examples/annotated_records.rs) reads two
literal fixtures, slices a complete CDS, writes GenBank, reverse-complements
EMBL with explicit loss reporting, and demonstrates two refusals. From the
repository root:

```bash
cargo bio annotated_records
```

The [worked guide](../examples/annotated-output-recipes.md) places the full input,
program, invocation and observed output together. `cargo bio` is ordinary Cargo
example execution.

## Writing contract

Call `write_genbank(&mut output, &record)` or `write_embl(&mut output, &record)`
with any Rust `Write` destination. Return type is `Result<(), BioIoError>`.
`record.format()` returns `AnnotatedFormat::{GenBank, Embl}`. A writer accepts
only its own source profile: cross-format calls fail before producing bytes.
Use `write_genbank_new(path, &record)` / `write_embl_new(path, &record)` for
create-only publication in a trusted directory. They never overwrite an existing
entry, including a dangling symlink. There is no compression or append option.

The writer renders the parsed model rather than returning `raw_record()`. It
preserves identity, description, accessions, declared molecule, source sequence
case, ordered metadata values, features, repeated qualifiers and decoded
qualifier values. Existing quoted/unquoted style and multiline qualifier value
breaks remain significant. Double quotes inside quoted values are escaped.
Unsupported locations can be written unchanged semantically without certifying
their validity; they still cannot be extracted or transformed.

Output uses LF and a final newline, canonical exact locations, 60-symbol rows,
10-symbol groups and regenerated position counters. Metadata indentation can be
retained separately from values. Physical whitespace, sequence wrapping and the
original spelling/wrapping of exact locations are not a byte-round-trip contract.
EMBL SQ lengths and standard A/C/G/T/other counts are recomputed from literal
case-insensitive sequence letters. A U is `other`, not a substituted T. Only
length-only SQ or the supported standard count grammar can be regenerated;
opaque SQ extensions are rejected, not silently deleted. This does not expand
the reader's strict RNA/U storage profile into general INSDC RNA compatibility.

Rendering is bounded and reparsed for supported semantic equality before the
external stream receives bytes. A stream I/O error can nevertheless leave a
prefix; callers own stream rollback/flush. Path helpers publish only after the
complete write succeeds, using the shared same-directory create-only mechanism.
There is no multi-file transaction, hostile-directory safety or power-loss
durability guarantee. Repeated stream calls can form a multi-record document,
but a later failure cannot roll back records already written.

## Transformation options

`AnnotationTransformOptions::new("LOCAL_fragment")?` requires a nonempty suffix
of ASCII letters, digits or underscores after `LOCAL_`, at most 16 total bytes.
The name must differ from the source locus, ID and accessions. A GenBank result
uses this local ID without an accession; modern EMBL requires `SV 1`, so its
record ID is `LOCAL_fragment.1`. This is a local serialization identity, not a
database-assigned accession/version. No old accession is incremented or copied.

GenBank transformations require `.with_date("16-OCT-2026")?`. Dates use exact
`DD-MMM-YYYY`, years 0001–9999, uppercase month abbreviations and calendar-valid
days including leap-year checks. EMBL can omit this option, or records the
supplied date as an explicitly local DT annotation. There is no clock lookup
or fabricated default date. Old database revision dates are not reused.

The default `AnnotationLossPolicy::Reject` refuses unpreservable information.
`.with_loss_policy(AnnotationLossPolicy::Discard)` opts into the specific
reported removals below; it does not enable partial-feature clipping or unknown
location interpretation. Methods return `AnnotatedTransform<A>` with `record()`,
`report()` and consuming `into_record()`. All source data remain unchanged.
Keep the report if provenance/loss matters; `into_record()` deliberately extracts
only the new record. New descriptions identify the operation and source ID.

## Locations and biological invariants

Rust slices are zero-based, half-open and checked: `3..12` selects nine symbols.
Negative indexes, strides, implicit clipping and circular wrap are unsupported.
Ordinary exact features must be entirely inside or entirely outside the selected
region, considering every join leaf. A partially intersected feature always
fails, even under `Discard`; wholly excluded features require that policy.
Nested joins retain written order, overlaps and repeated segments. The result
of a proper slice is linear; a full-length slice retains recognized topology.

Reverse complement supports actual DNA/RNA alphabets, preserves the sequence's
case pattern in reversed positions, and retains recognized topology. It maps
each leaf `[a,b)` to `complement([n-b,n-a))`, recursively preserving join order
and simplifying double complements. Feature vector order does not get sorted
by its new coordinates. For every retained ordinary feature, extracting its
new location from the new sequence yields the same biological sequence as
before transformation. CDS translation/codon-start qualifiers therefore need
no automatic recomputation for these complete features.

One explicit exception is a plain `source` range spanning the entire record:
it denotes coverage. A nonempty slice resizes it to the new extent; reverse
complement leaves it as plain whole-record coverage. An empty slice cannot
retain a nonempty source location and needs explicit discard. Other `source`
expressions receive ordinary exact-feature rules. The report marks coverage
changes separately from ordinary coordinate mapping.

Fuzzy, remote, between-base, `order` and other unsupported locations reject all
transforms, including full-length slices, even with `Discard`. Unknown topology
also fails. A complemented feature on protein/custom alphabets rejects even a
slice or an otherwise permitted whole-feature removal; its biological extraction
is unsupported. There is no remote resolver, automatic CDS translation, frameshift
repair or inference of biological semantics from free text.

## Metadata and qualifiers

Local identity, description and structural header replacement are explicitly
authorized by requesting a derived record and are reported. Old accessions and
versions are removed. GenBank SOURCE/ORGANISM and EMBL OS/OC/OG taxonomy fields
are retained. Unreviewed metadata such as references, reference ranges,
DBLINK/DR (including MD5), comments, keywords and old revision dates either fail
under `Reject` or are removed individually under `Discard`. Structural values
are checked; custom text is not silently treated as a harmless separator.
Nonempty ORIGIN, opaque SQ payloads and extra trailing LOCUS tokens need explicit
discard before regeneration, with old/new payloads retained in the report.
The caller supplies the new LOCUS date, but the source division must be present
and recognized; no biological classification is fabricated. An unchanged
writer refuses opaque SQ outright because it has no loss-policy argument.

For retained features, the reviewed coordinate-independent qualifier names are:

```text
gene gene_synonym locus_tag product protein_id codon_start transl_table
translation pseudo pseudogene organism mol_type strain isolate db_xref
culture_collection plasmid organelle
```

Their values/order remain annotations from the source, not newly verified
biological assertions. In particular, preserving `/translation` does not
validate that it was correct initially. Qualifiers outside this list—including
`note`, `transl_except`, `anticodon`, `rpt_unit_range` and `direction`—fail by
default. Under `Discard`, the **entire feature** is removed and reported. This
prevents removing an exceptional-codon qualifier while leaving a misleading
dependent translation. Expanding this list requires a separate reviewed contract.

## Source classification

Derived headers retain the source division, not a guessed classification based
on the local filename or identifier. In particular, `SYN` means synthetic sequence,
not locally generated output. The supported division tokens are format-specific:

| Source profile | Retained division tokens |
|---|---|
| GenBank | `PRI ROD MAM VRT INV PLN BCT VRL PHG SYN UNA EST PAT STS GSS HTG HTC ENV` |
| EMBL | `PHG ENV FUN HUM INV MAM VRT MUS PLN PRO ROD SYN TGN UNC VRL` |

EMBL also retains its source data class: `STD PAT EST GSS HTC HTG WGS TSA STS`.
The assembly class `CON` is outside this literal-sequence transformation profile.
Missing/unknown classification, codes belonging only to the other format, and
`CON` all produce `unsupported_format_feature`, even with `Discard`. Writers of
unchanged records retain opaque headers without claiming to validate classification.
These field meanings follow the [NCBI field guide](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
and [ENA release 143 manual, sections 3.1–3.2](https://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/doc/usrman.txt).
Retaining a source label is not an independent taxonomic validation.

For a minimal reproducible example, save this complete program as
`seq-bio/examples/classification.rs`. The literal four-base test input labels
its source `PRI`; it is a synthetic test of header handling, not a real accession.

```rust
use bioferris::Dna;
use bioferris::io::annotated::{AnnotationTransformOptions, parse_genbank};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = "LOCUS       ORIGINAL 4 bp DNA linear PRI 15-OCT-2026\n\
ORIGIN\n\
        1 ACGT\n\
//\n";
    let options = AnnotationTransformOptions::new("LOCAL_part")?
        .with_date("16-OCT-2026")?;
    let records = parse_genbank::<Dna>(input)?;
    let derived = records[0].slice(1..3, options.clone())?;
    let header = &derived.record().metadata()[0].value;

    let unknown = parse_genbank::<Dna>(&input.replace("PRI", "ZZZ"))?;
    let error = match unknown[0].slice(1..3, options) {
        Err(error) => error,
        Ok(_) => return Err("unknown classification was unexpectedly accepted".into()),
    };
    let observed = format!("{header}\n{}", error.code());
    assert_eq!(
        observed,
        "LOCAL_part 2 bp DNA linear PRI 16-OCT-2026\nunsupported_format_feature"
    );
    println!("{observed}");
    Ok(())
}
```

From the repository root, after saving the file:

```bash
cargo bio classification --offline --locked
```

Complete stdout, empty stderr, exit 0 (the expected refusal is handled):

```text
LOCAL_part 2 bp DNA linear PRI 16-OCT-2026
unsupported_format_feature
```

The half-open slice contains `CG`, retains `PRI`, and does not become `SYN`.
Changing the source division to `ZZZ` fails instead of guessing. The complete
program and its output assertion are executed by the Markdown Rustdoc checks;
the integration suite also checks both loss policies and both source formats.

## Reports and agent use

The typed `AnnotationTransformReport` exposes source/derived IDs, an
`AnnotationOperation` and ordered `changes`. Source indexes are zero-based:

| Change | Meaning |
|---|---|
| `MetadataChanged` | Source index (or `None` for an addition), key, before/after optional values and reason; `after: None` means removal |
| `LocationChanged` | Original feature index, old/new location text and a source-coverage marker |
| `FeatureRemoved` | Original feature index, key, original location and reason; its qualifiers leave with it |

Inspect these types directly; do not parse a display string. Errors expose
`kind()`, `code()`, `format()` and a diagnostic message. Agents write and run
the same library programs as humans and must explicitly choose a
loss policy, new local identity and applicable date. Raw record/feature/location/
qualifier text on a successful derived record is regenerated coherently by
rendering/reparsing, never left pointing at the original coordinates.

## Errors and limits

| Category | Typical cause / response |
|---|---|
| `invalid_header` | Invalid/reused local identity or invalid/missing required date; supply explicit valid options |
| `length_mismatch` | Reversed/out-of-bounds slice; check the original length |
| `unsupported_format_feature` | Cross-format writing, protein reverse complement, unsupported location/topology/classification or unrepresentable output; choose a supported profile |
| `data_loss` | Partial feature, excluded feature or unreviewed metadata/qualifier; use a different region or review explicit discard; partial clipping stays prohibited |
| `limit_exceeded` | Default line, record or document bounds; do not bypass validation |
| `io_error` | Stream/filesystem failure, existing destination or publication failure; retain the original and choose a new trusted path |

The existing defaults apply: 8 MiB physical lines, 100 million sequence symbols,
1 GiB serialized document and 64 location-nesting levels. APIs are eager, not
streaming parsers or zero-copy transformations. A canonical location must fit a
single bounded physical line even if its input was wrapped. Model, output and
reparsed model can coexist. These are logical input/output bounds, not
peak-process-memory or running-time guarantees. No new configurable annotated
budgets are introduced.

## Compatibility and verification

This is a deliberately bounded clean-room profile, not full Biopython or INSDC
equivalence. [Biopython 1.88](https://biopython.org/docs/1.88/api/Bio.SeqRecord.html)
documents feature-adjusting reverse complement and defaults that omit general
annotations. BioFerris additionally requires local identity and explicit loss
handling; it does not reproduce Python slicing/default-copy behavior.

Format references consulted 2026-10-16: [INSDC Feature Table Definition 11.4](https://www.insdc.org/submitting-standards/feature-table/),
[NCBI GenBank sample/field descriptions](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
and [ENA release 143 manual](https://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/doc/usrman.txt).
Local derived records are analysis artifacts, not database submissions. Actual
checks and limitations are recorded in the validation report.
The [worked recipes](../examples/annotated-output-recipes.md) are executable
contracts, not a live reference-runtime comparison or new scientific campaign.
