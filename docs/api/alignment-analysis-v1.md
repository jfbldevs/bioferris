# Alignment analysis: coordinates, slices, statistics and consensus

## Contents

- [Availability and purpose](#availability-and-purpose)
- [Write a Rust analysis](#write-a-rust-analysis)
- [Understand cells and adapters](#understand-cells-and-adapters)
- [Map and slice source coordinates](#map-and-slice-source-coordinates)
- [Interpret identity and gap statistics](#interpret-identity-and-gap-statistics)
- [Choose consensus policies](#choose-consensus-policies)
- [Handle annotations explicitly](#handle-annotations-explicitly)
- [Handle missing data without guessing](#handle-missing-data-without-guessing)
- [Limits and errors](#limits-and-errors)
- [Compatibility and verification](#compatibility-and-verification)

[API index](../README.md) · [Analysis contract](../algorithms/alignment-analysis-v1.md) ·
[Input format profiles](alignment-input-v1.md) ·
Implementation plan

## Availability and purpose

The Rust API and Cargo example are implemented and verified on Rust 1.85.0.
Stage 4 is accepted for the documented analysis/policy scope; see the
main validation record. That
record also covers CLI operations, strict JSON requests and discovery, which were
removed with the `bioferris` CLI (decision 029). The complete
program and the real Cargo alias below pass the tests. Earlier isolated/staging
evidence remains a separate checkpoint, not a claim of package release or
universal Biopython compatibility.

Use `AlignmentAnalysis<A>` to inspect an existing DNA, RNA or protein alignment.
It never realigns sequences or changes scientific kernels.

| Task | Rust entry point |
|---|---|
| Slice rows/columns without renumbering source positions | `slice_rows`, `slice_columns`, `slice_with` |
| Map a column to source or source to column | `column_to_sequence`, `sequence_to_column` |
| Count gaps, missing cells and literal identity | `statistics` |
| Produce a gapped, policy-qualified consensus | `consensus` |

Keep alignment columns, ungapped positions and original sequence positions
distinct. They coincide only for a complete ungapped forward row with origin 0.

## Write a Rust analysis

This complete program includes all literal FASTA/Stockholm/NEXUS inputs, explicit
policies and a small accepted Smith–Waterman result. Complete runnable source:
[`seq-bio/examples/alignment_analysis.rs`](../../seq-bio/examples/alignment_analysis.rs).

<!-- alignment-analysis-rust-program -->
```rust
//! Checked alignment statistics, source coordinates, consensus and annotation loss.

use bioferris::io::alignment::{AlignmentFormat, parse_alignment};
use bioferris::prelude::*;

const FASTA: &str = ">reference\nAC-GT\n>sample\nA--GT\n>ambiguous\nAN-G-\n";
const STOCKHOLM: &str = "# STOCKHOLM 1.0\n#=GF ID original\n# comment with coordinates 1-4\n#=GS reference DE original reference row\nreference AC-G\nsample ACTG\n#=GC SS_cons <..>\n#=GR reference SS <..>\n//\n";
const NEXUS: &str = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA GAP=~ MISSING=!;\nMATRIX\na A~!C\nb CG!C\n;\nEND;\n";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = parse_alignment::<Dna>(FASTA.as_bytes(), AlignmentFormat::Fasta)?;
    let alignment = AlignmentAnalysis::from_document(&document)?;
    report_statistics(&alignment)?;
    report_coordinates_and_slice(&alignment)?;
    report_consensus(&alignment)?;
    report_missing_data()?;
    report_annotation_loss()?;
    report_local_pairwise_origin()?;
    Ok(())
}

fn report_statistics(alignment: &AlignmentAnalysis<Dna>) -> Result<(), Box<dyn std::error::Error>> {
    let stats = alignment.statistics(AlignmentStatisticsOptions::default())?;
    println!("Rectangle: {} rows x {} columns", stats.rows, stats.columns);
    println!("Cells: {:?}", stats.cells);
    println!("Gap runs: {:?}", stats.gap_runs);
    println!("Pairs: {:?}", stats.pairs);
    println!("Identity: {:?} / {}", stats.identity, stats.denominator);
    let ignored = alignment.statistics(AlignmentStatisticsOptions {
        ambiguity: AlignmentAmbiguityPolicy::Ignore,
        denominator: IdentityDenominator::ResidueAndSingleGap,
    })?;
    println!(
        "Ignore ambiguity, include single gaps: {:?} / {}",
        ignored.identity, ignored.denominator
    );
    Ok(())
}

fn report_coordinates_and_slice(
    alignment: &AlignmentAnalysis<Dna>,
) -> Result<(), Box<dyn std::error::Error>> {
    let oriented =
        alignment.with_row_coordinates(0, RowCoordinates::new(20, 10..14, Strand::Reverse)?)?;
    let positions: Result<Vec<_>, _> = (0..oriented.width())
        .map(|column| oriented.column_to_sequence(0, column))
        .collect();
    println!("Reverse source positions: {:?}", positions?);
    println!(
        "Source position 12 -> column {:?}",
        oriented.sequence_to_column(0, 12)?
    );
    println!(
        "Unaligned source position 0 -> {:?}",
        oriented.sequence_to_column(0, 0)?
    );
    let sliced = oriented.slice_with(0..2, 1..4, AlignmentAnnotationPolicy::Reject)?;
    for row in sliced.alignment().rows() {
        println!(
            "Slice {}: {} {:?}",
            row.identifier(),
            row.aligned(),
            row.coordinates()
        );
    }
    Ok(())
}

fn report_consensus(alignment: &AlignmentAnalysis<Dna>) -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "Default consensus: {}",
        alignment.consensus(ConsensusOptions::default())?.sequence
    );
    let voted = alignment.consensus(ConsensusOptions {
        gaps: ConsensusGapPolicy::Vote,
        ..Default::default()
    })?;
    println!("Gap-voting consensus: {}", voted.sequence);
    for column in &voted.columns {
        println!(
            "Column {}: {:?}, votes {}/{}, winners {:?}",
            column.column,
            column.status,
            column.winning_count,
            column.denominator,
            String::from_utf8(column.winners.clone())?
        );
    }
    for min_support in [0.3, 0.5] {
        let resolved = alignment.consensus(ConsensusOptions {
            gaps: ConsensusGapPolicy::Vote,
            ties: ConsensusTiePolicy::AsciiOrder,
            min_support,
            ..Default::default()
        })?;
        println!(
            "ASCII ties at support {min_support}: {}, column 1 {:?}",
            resolved.sequence, resolved.columns[1].status
        );
    }
    println!(
        "Consensus ignoring ambiguity: {}",
        alignment
            .consensus(ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })?
            .sequence
    );
    report_consensus_rejections(alignment);
    let empty = alignment.slice_rows(0..0)?;
    println!(
        "Zero-row slice: {} x {}, consensus {}",
        empty.sequence_count(),
        empty.width(),
        empty.consensus(ConsensusOptions::default())?.sequence
    );
    Ok(())
}

fn report_consensus_rejections(alignment: &AlignmentAnalysis<Dna>) {
    for (label, options) in [
        (
            "Reject gaps",
            ConsensusOptions {
                gaps: ConsensusGapPolicy::Reject,
                ..Default::default()
            },
        ),
        (
            "Reject ambiguity",
            ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Reject,
                ..Default::default()
            },
        ),
        (
            "Reject ties",
            ConsensusOptions {
                ties: ConsensusTiePolicy::Error,
                ..Default::default()
            },
        ),
        (
            "Invalid support",
            ConsensusOptions {
                min_support: 1.1,
                ..Default::default()
            },
        ),
    ] {
        if let Err(error) = alignment.consensus(options) {
            println!("{label}: {}: {error}", error.code());
        }
    }
}

fn report_missing_data() -> Result<(), Box<dyn std::error::Error>> {
    let nexus = parse_alignment::<Dna>(NEXUS.as_bytes(), AlignmentFormat::Nexus)?;
    let missing = AlignmentAnalysis::from_document(&nexus)?;
    println!(
        "Custom glyphs: {} / {}",
        missing.rows()[0].aligned(),
        missing.rows()[1].aligned()
    );
    if let Err(error) = missing.column_to_sequence(0, 0) {
        println!("Missing origin: {}: {error}", error.code());
    }
    println!(
        "Consensus ignoring missing data: {}",
        missing.consensus(ConsensusOptions::default())?.sequence
    );
    if let Err(error) = missing.consensus(ConsensusOptions {
        missing: ConsensusMissingPolicy::Reject,
        ..Default::default()
    }) {
        println!("Reject missing data: {}: {error}", error.code());
    }
    Ok(())
}

fn report_annotation_loss() -> Result<(), Box<dyn std::error::Error>> {
    let annotated = parse_alignment::<Dna>(STOCKHOLM.as_bytes(), AlignmentFormat::Stockholm)?;
    let annotated = AlignmentAnalysis::from_document(&annotated)?;
    if let Err(error) = annotated.slice_columns(1..3) {
        println!("Strict annotation slice: {}: {error}", error.code());
    }
    let discarded = annotated.slice_with(0..2, 1..3, AlignmentAnnotationPolicy::Discard)?;
    for annotation in discarded.discarded_annotations() {
        println!(
            "Discard: {:?} {:?} {:?} line {} {:?}",
            annotation.scope(),
            annotation.sequence(),
            annotation.tag(),
            annotation.line(),
            annotation.value()
        );
    }
    Ok(())
}

fn report_local_pairwise_origin() -> Result<(), Box<dyn std::error::Error>> {
    let first = dna("TTACGTAA")?;
    let second = dna("GGACGTCC")?;
    let run = Pairwise::smith_waterman().run(&first, &second)?;
    let local = AlignmentAnalysis::from_pairwise(run.first(), run.second(), &run)?;
    println!(
        "Local first-row origin: {:?}",
        local.rows()[0].coordinates()
    );
    println!(
        "Local first column -> source {:?}",
        local.column_to_sequence(0, 0)?
    );
    Ok(())
}
```

Run from the repository root:

<!-- alignment-analysis-rust-command -->
```bash
cargo bio alignment_analysis
```

Complete observed stdout; stderr is empty; exit 0:

<!-- alignment-analysis-rust-output -->
```text
Rectangle: 3 rows x 5 columns
Cells: AlignmentCellCounts { residues: 10, gaps: 5, missing: 0, ambiguous: 1 }
Gap runs: AlignmentGapRuns { leading: 0, internal: 3, trailing: 1 }
Pairs: AlignmentPairCounts { total: 15, missing: 0, excluded_ambiguity: 0, double_gap: 3, single_gap: 4, residue_pairs: 8, identities: 7, mismatches: 1 }
Identity: Some(0.875) / 8
Ignore ambiguity, include single gaps: Some(0.7) / 10
Reverse source positions: [Some(13), Some(12), None, Some(11), Some(10)]
Source position 12 -> column Some(1)
Unaligned source position 0 -> None
Slice reference: C-G Some(RowCoordinates { source_length: 20, range: 11..13, strand: Reverse })
Slice sample: --G Some(RowCoordinates { source_length: 3, range: 1..2, strand: Forward })
Default consensus: ANNGT
Gap-voting consensus: AN-GT
Column 0: Accepted, votes 3/3, winners "A"
Column 1: Tie, votes 1/3, winners "-CN"
Column 2: Gap, votes 3/3, winners "-"
Column 3: Accepted, votes 3/3, winners "G"
Column 4: Accepted, votes 2/3, winners "T"
ASCII ties at support 0.3: A--GT, column 1 TieResolved
ASCII ties at support 0.5: AN-GT, column 1 BelowThreshold
Consensus ignoring ambiguity: ACNGT
Reject gaps: rejected_cell: policy rejects gap byte 45 at row 1, column 1
Reject ambiguity: rejected_cell: policy rejects ambiguity byte 78 at row 2, column 1
Reject ties: consensus_tie: consensus tie at column 1 between ASCII glyphs [67, 78]
Invalid support: invalid_minimum_support: consensus minimum support must be finite and in 0..=1; received 1.1
Zero-row slice: 0 x 5, consensus NNNNN
Custom glyphs: A~!C / CG!C
Missing origin: missing_coordinates: row 0 has no known residue-to-source mapping; missing data is not guessed
Consensus ignoring missing data: NGNC
Reject missing data: rejected_cell: policy rejects missing byte 33 at row 0, column 2
Strict annotation slice: annotation_loss: slice could invalidate 5 annotations; explicit discard permission is required
Discard: File None "ID" line 2 "original"
Discard: Comment None "comment" line 3 " comment with coordinates 1-4"
Discard: Sequence Some("reference") "DE" line 4 "original reference row"
Discard: Column None "SS_cons" line 7 "<..>"
Discard: Residue Some("reference") "SS" line 8 "<..>"
Local first-row origin: Some(RowCoordinates { source_length: 8, range: 2..6, strand: Forward })
Local first column -> source Some(2)
```

The reverse origin maps the first reference residue to source position 13,
not 10. A gap maps to `None`; a valid source position outside the aligned interval
also maps to `None`. Neither is a missing-data error. The slice's positions remain
in the original source interval, so the reference origin is 11..13, not 0..2.

The statistics count seven equal residue pairs among eight eligible residue
pairs, giving 0.875. Including single-gap pairs while excluding ambiguous
evidence gives seven identities among ten pairs, 0.7. Column 2 is all gaps:
default consensus has no evidence there and emits N; gap voting emits a
supported gap. Column 1 is a three-way tie under gap voting, retained as status
`Tie` rather than mislabeled as a supported N. Explicit ASCII resolution at
support 0.3 accepts `-` in that tied column; at support 0.5 it emits unknown with
`BelowThreshold`. The example also demonstrates literal versus ignored
ambiguity, all three cell-rejection policies, error-policy ties and an invalid
support threshold. These decisions are not inferred from the output glyph alone.

The annotated slice refuses by default. Explicit discard returns each removed
annotation's full original scope, target, tag, source line and value, including
the structural annotation and comment. The local pairwise adapter retains origin
2..6 within its full eight-base source, rather than inventing origin zero.

## Understand cells and adapters

`AlignmentAnalysis<A>` is a separate immutable rectangle. It does not weaken the
accepted computational `MultipleAlignment` invariants. It can retain all-gap
rows/columns, zero columns and zero-row slices with their explicit width.

Each cell is one of:

- `AlignmentCell::Residue(Symbol<A>)`: a validated normalized residue, including
  permitted ambiguity symbols.
- `AlignmentCell::Gap(original_byte)`: an explicit alignment gap such as
  `-`, Stockholm `.`, or a declared NEXUS gap glyph.
- `AlignmentCell::Missing(original_byte)`: unknown data, not an assumed residue
  or gap, such as `?` or a declared NEXUS missing glyph.

`rows()` retains source order; each row exposes `identifier()`, `cells()`,
`aligned()` and optional `coordinates()`. `width()`, `sequence_count()`,
`annotations()` and `limits()` inspect the rectangle.

Constructors return `Result` because dimensions and adapters are checked:

| Source | Constructor | Coordinate origin |
|---|---|---|
| Parsed `AlignmentDocument<A>` | `from_document(&document)` | Ungapped positions unless missing data prevents mapping |
| Existing `MultipleAlignment<A>` | `from_multiple(&alignment)` | Ungapped positions |
| Existing `MultipleAlignmentView<A>` result | `from_multiple_result(&result)` | Ungapped positions |
| Existing `PairwiseAlignmentView` plus original sequences | `from_pairwise(&first, &second, &result)` | Full source lengths and traceback ranges, including local offsets |

Every constructor has a `_with_limits` variant. Distance-only/value results have
no traceback and cannot be adapted as if they contained one. Pairwise adapters
render against the supplied originals and compute identity from rendered residue
symbols, not from unchecked Match/Substitution operation labels.
Adapters do not reconstruct missing original genomic coordinates from IDs.

## Map and slice source coordinates

All row/column indices and source positions are zero-based. All ranges are
ascending, half-open. Negative, reversed and out-of-bounds ranges are not Python
slice clipping. `slice_rows(rows)` and `slice_columns(columns)` default to
annotation refusal; `slice_with(rows, columns, policy)` returns an
`AlignmentSlice<A>` with `alignment()`, `discarded_annotations()`, and
`into_parts()`.

`RowCoordinates::new(source_length, range, strand)` declares an explicit source
interval. `with_row_coordinates(row, origin)` returns a new rectangle after
checking that range length equals the represented residues. This is a caller
assertion: cells are **not** reversed or complemented, and names do not infer
strand. The interval must lie inside the full source length.

For origin [s,e), residue k maps to s+k forward or e-1-k reverse. In
`AC-GT` with reverse origin 10..14, columns map to 13,12,None,11,10. Selecting
columns 1..4 gives `C-G` with source interval 11..13 and maps 12,None,11.
A valid source position outside the selected interval maps to no column;
a position at or beyond full source length is an error. Empty slices retain a
valid empty interval at the relevant oriented boundary.

Rows containing missing data have no inferred mapping anywhere, even at a
concrete cell. Attaching an origin while missing cells remain is refused.
Slicing away all missing cells does not magically recover one: the caller may
then attach a justified explicit origin. See the
[missing-data recipe](#handle-missing-data-without-guessing).

## Interpret identity and gap statistics

`statistics(AlignmentStatisticsOptions)` returns exact aggregate counts plus
ordered `row_statistics` and `column_statistics`.

`cells.ambiguous` is a subset of `cells.residues`, not another disjoint cell
state. DNA/RNA ambiguity means non-ACGT/U; protein ambiguity means outside the
canonical 20 residues. Counts always describe source cells, even when a comparison
policy excludes some of them.

Gap runs are leading, internal or trailing. A wholly gapped nonempty row has one
**leading** run only, so these categories partition runs. Missing cells split
runs. Zero-column rows have no runs.

Every unordered row pair in each column is assigned once, in this priority:

1. `missing`: at least one missing cell.
2. `excluded_ambiguity`: a nonmissing pair containing ambiguity under Ignore.
3. `double_gap`: two gaps, never a residue identity.
4. `single_gap`: one residue and one gap.
5. `residue_pairs`: two retained residues, split into identities and mismatches.

`pairs.total` equals the sum of the first five categories.
`residue_pairs = identities + mismatches`. Equality is literal, not IUPAC set
intersection or substitution-matrix similarity.

| Option | Values / meaning |
|---|---|
| `ambiguity` | `Literal` (default): compare ambiguity glyphs literally; `Ignore`: exclude affected evidence; `Reject`: return a typed error |
| `denominator` | `ResiduePairs` (default); `ResidueAndSingleGap`; `AllPairs`, explicitly including missing, excluded and double-gap pairs |

`identity = identities / denominator`. A zero denominator returns `None`, not
NaN or perfect identity. Policy and
denominator accompany each result; percentages from different policies are not
interchangeable. Counts use checked integer arithmetic. The reported fraction is
binary64; converting counts above 2^53 may round, while the integer numerator and
denominator remain available exactly. Statistics ambiguity refusal traverses rows
then columns; consensus refusal uses columns then rows.

## Choose consensus policies

`consensus(ConsensusOptions)` returns a gapped string and one evidence record per
column. It does not return an ungapped typed `Sequence` or an inferred ancestor.

| Independent option | Rust variants | Default |
|---|---|---|
| `gaps` | `Ignore`, `Vote`, `Reject` | `Ignore` |
| `missing` | `Ignore`, `Reject` | `Ignore` |
| `ambiguity` | `Literal`, `Ignore`, `Reject` | `Literal` |
| `ties` | `Unknown`, `AsciiOrder`, `Error` | `Unknown` |
| `min_support` | Finite inclusive threshold in 0..=1 | 0.5 |

Ignored cells leave the voting denominator. Gap voting combines gap glyphs into
one `-` ballot without changing stored rows. A unique accepted gap emits `-`.
The canonical unknown output is N for DNA/RNA or X for proteins. No-evidence
columns emit that unknown even at a zero threshold.

Each column retains `column`, `symbol`, `status`, `winning_count`,
`denominator`, optional `support`, and all maximum-vote `winners` in ASCII
order. These statuses distinguish biologically different outputs:

| Status | Meaning |
|---|---|
| `Accepted` | Unique residue winner reaches minimum support |
| `Gap` | Unique gap winner reaches minimum support |
| `Tie` | Unknown-policy tied votes emit N/X |
| `TieResolved` | Explicit ASCII-order tie resolution reaches support |
| `BelowThreshold` | Winning support is insufficient |
| `NoEvidence` | No voting cells; support is absent |

Ties are handled **before** minimum support. Error/unknown tie policies therefore
remain observable below the threshold. ASCII-order selects the smallest glyph
(including `-`) but still must meet support. Consensus refusal traverses columns
then rows and returns no partial success.

## Handle annotations explicitly

The program's `STOCKHOLM` input carries five annotations: a file `ID`, a
comment, a sequence `DE`, a column `SS_cons` and a residue `SS` line.

Read-only analysis preserves every annotation. An actual identity slice retains
them even with discard permission. Any other annotated slice either fails or
requires explicit all-annotation discard and a complete ordered loss report.
Opaque structural text cannot safely be clipped just because its width matches:
paired secondary-structure symbols can refer outside the slice, row selection
can stale a consensus, and text can contain coordinates or checksums.

In the program's output, selecting columns 1..3 fails with `annotation_loss`
because all five annotations could become stale; the error retains each one's
scope, target, tag, one-based source line and original value. The same selection
with `AlignmentAnnotationPolicy::Discard` succeeds, yields a rectangle without
annotations and returns all five in `discarded_annotations()`. No structural
clipping is claimed.

## Handle missing data without guessing

The program's `NEXUS` input declares `~` as gap and `!` as missing data.
`AlignmentAnalysis` preserves both glyphs and their distinct typed states
(`Custom glyphs: A~!C / CG!C`). The literal A at row 0/column 0 does not justify
guessing its original source coordinate when the row also contains missing data:
mapping row 0 fails with `missing_coordinates`, an error rather than the `None`
returned for a known gap or a valid unaligned source position.

## Limits and errors

`AlignmentAnalysisLimits` defaults are `max_rows=100000`, `max_columns=1000000`
and `max_cells=10000000`; pass other limits through the `_with_limits`
constructors.
The product is checked for overflow and budget compliance before allocating the
analysis copy. The existing eager document parser runs **before** that dimension
check, under its separate 1 GiB decompressed input limit. These options do not
promise an exact process-memory cap or streamed alignment analysis.

The current eager alignment parser has a fixed decompressed-input bound; analysis
limits do not replace or configure that parser bound. Derived slices stay within
their source limits. Empty dimensions and all-gap columns are allowed; a zero-row selection keeps width,
so consensus remains one unknown per retained column.

| Stable code | Meaning / response |
|---|---|
| `limit_exceeded`, `count_overflow` | Explicit dimensions/counts cannot be represented or exceed budget; no partial result |
| `invalid_range`, `invalid_row`, `invalid_column` | Check strict zero-based bounds |
| `invalid_source_position` | Position is at/beyond the declared full source length |
| `invalid_coordinates` | Invalid interval or residue-count mismatch; do not guess a new origin |
| `missing_coordinates` | Missing data or retained unknown mapping prevents coordinate inference |
| `annotation_loss` | Nontrivial annotated slice needs a deliberate discard decision; full affected metadata retained |
| `rejected_cell` | Chosen gap/missing/ambiguity policy refused a concrete row/column/byte |
| `invalid_minimum_support` | Supply a finite threshold in 0..=1 |
| `consensus_tie` | Error-policy tie; tied column/glyphs are retained |
| `invalid_pairwise`, `invalid_cell` | Source/result cannot be represented consistently; inspect typed cause/context |
| Existing I/O/format codes | Filesystem, malformed alignment, decompression or parser-limit failure while reading the input; inspect the stable code and input/format/source context |

`AlignmentAnalysisError` is typed and implements `Error` with `code()`. Use
stable codes for control flow, not message parsing.

## Compatibility and verification

The [analysis contract](../algorithms/alignment-analysis-v1.md) defines BioFerris
policy rather than blanket Biopython slicing/consensus parity. Existing imported
metadata, all-gap columns and missing states are retained; checked analysis does
not broaden the computational `MultipleAlignment` model or rewrite kernels.

Standalone core checks cover 55,980 explicit-pair statistical comparisons,
251,910 independently sorted-ballot consensus-policy comparisons and 2,186
gapped forward/reverse rectangles with inverse maps and composed slices.
Accepted Smith–Waterman and MAFFT results are adapted in core tests.
These are scoped isolated oracle checkpoints, not independently a main-workspace
acceptance result or newly qualified scientific alignment fidelity.

The tests in `seq-bio/tests/alignment_analysis_recipes.rs` check the complete
program fence and the real `cargo bio alignment_analysis` alias against the
output above. The standalone/staging oracle checks above remain separate
checkpoints.

The main validation record
records accepted Stage 4 scope, final toolkit/optimized checks, complete workspace
harness output and terminal audit, documentation and response comparisons.
Final source/evidence review and the local commit accompany this acceptance. This
analysis introduces no Python runtime, network client, file export, performance
claim or newly qualified alignment-engine fidelity.
