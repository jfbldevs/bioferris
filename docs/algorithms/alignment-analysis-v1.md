# Alignment analysis and coordinate contract

## Contents

- [Scope and implementation state](#scope-and-implementation-state)
- [Analysis model and adapters](#analysis-model-and-adapters)
- [Coordinates and slices](#coordinates-and-slices)
- [Annotations](#annotations)
- [Statistics](#statistics)
- [Consensus](#consensus)
- [Acceptance](#acceptance)

[API index](../README.md) · Implementation plan

Complete manual workflows are in the
[analysis guide](../api/alignment-analysis-v1.md). Their isolated/staging checks
are distinct from final main-checkout verification and stage acceptance.

## Scope and implementation state

Stage 4 implemented and verified public Rust/prelude, CLI and strict-JSON
interfaces; the CLI and strict-JSON interfaces were later removed with the
`bioferris` CLI (decision 029). The documented analysis/policy scope is accepted
in the main validation record. It includes
checked row/column slicing, sequence/alignment coordinate mapping, gap/identity
statistics and explicit-policy consensus. References are the
[Biopython 1.88 alignment guide](https://biopython.org/docs/1.88/Tutorial/chapter_align.html)
and [Stockholm format](https://sonnhammer.sbc.su.se/Stockholm.html). The contracts
below define conservative BioFerris behavior rather than blanket slicing parity.

## Analysis model and adapters

Use a separate immutable typed analysis rectangle. The existing computational
`MultipleAlignment` rejects all-gap columns, but valid imported alignments and
row slices may contain them. Analysis must retain residue, gap and missing cells
as distinct states, including source glyphs and ordered imported annotations.
An empty row slice retains an explicit column count; zero columns and all-gap
rows/columns are valid. No accepted computational invariant is weakened.

Adapters accept existing multiple-alignment views, canonical pairwise alignment
views with their two original sequences, and `io::alignment::AlignmentDocument`.
Distance-only results have no traceback and are ineligible. Pairwise comparisons
use actual rendered symbols, not unverified Match/Substitution operation labels.
Existing scientific kernels and parsers are reused without modification.

Explicit limits bound rows, columns and total cells before new model allocation.
Defaults are 100,000 rows, 1,000,000 columns and 10,000,000 cells. Derived slices
remain within their source limits. Overflow is an error, not wrapped dimensions.

## Coordinates and slices

All ranges are zero-based and half-open; reversed/out-of-bounds ranges fail.
Pairwise rows retain complete source lengths and alignment source ranges, including
local origins. MSA/import rows default to ungapped-row positions. A valid source
position outside the aligned range maps to no column; beyond-source positions
are errors. Gaps map to no sequence position. Slices preserve source coordinates
instead of renumbering selected residues from zero.

Explicit row-origin metadata can describe a forward or reverse source interval.
Its length must equal the row's represented residues and be within source length.
For origin `[s,e)`, residue k maps to `s+k` forward or `e-1-k` reverse.
A selected oriented residue interval `[a,b)` maps to `[s+a,s+b)` or
`[e-b,e-a)`, respectively. No orientation is inferred from row names or X-drop
extension direction. Rows with missing cells refuse coordinate mapping because
missing data need not establish a residue/gap assignment. Removing all missing
cells by slicing does not invent an origin: mapping remains unknown, but the
caller may then attach a valid explicit origin to that now-complete row.

For `AC-GT` with reverse origin `[10,14)`, column positions are
13,12,None,11,10; selecting columns1..4 retains12,None,11.

## Annotations

Read-only analysis preserves annotations. A true identity slice may retain them.
Any nontrivial annotated slice defaults to an annotation-loss error. Explicit
discard returns the derived rectangle plus an ordered typed loss report retaining
each removed annotation's scope, target, tag, source line and value.

Width-checked opaque text is not necessarily safe to clip: Stockholm structural
annotations may contain paired symbols, row selection can invalidate CLUSTAL
consensus, and metadata may contain coordinates or checksums. Generic substring
projection must not be presented as biological annotation preservation.

## Statistics

Report exact residue/gap/missing totals by row and column, and leading/internal/
trailing gap runs. A run spanning an entirely gapped row is classified as leading
only, so the run categories partition the total. Missing cells separate gap runs.

For each column, classify every unordered row pair with priority: missing pair,
excluded-ambiguity pair, double-gap pair, single-gap pair, residue pair. Literal
residue equality contributes identity; double gaps never contribute residue
identity. Residue-pair mismatches are residue pairs minus identities. Aggregate
counts can be derived from per-symbol multiplicities rather than allocating a
quadratic pair matrix. Use checked integer arithmetic throughout.

Ambiguity policies are literal (default), ignore or reject. DNA/RNA ambiguity is
outside A/C/G/T/U; canonical protein analysis uses the standard 20 residues.
Denominator policies are residue pairs only (default), residue plus single-gap
pairs, or all column pairs. The last explicitly includes missing/ambiguity/double
gap categories in its denominator. Return no fraction for a zero denominator;
never fabricate NaN or perfect identity for absent evidence.

## Consensus

Return a gapped consensus string plus per-column evidence, not an ungapped
`Sequence`. Options independently select gap handling (ignore/vote/reject),
missing handling (ignore/reject), ambiguity handling (literal/ignore/reject),
tie handling (unknown/ASCII-order/error), and a finite minimum support in 0..1.
Defaults ignore gaps/missing, count literal ambiguity, resolve ties to unknown,
and require inclusive minimum support 0.5. Ignored cells are excluded from the
denominator. Canonical unknown is N for DNA/RNA and X for proteins.

A unique gap winner under gap voting yields `-`. Ties, insufficient support and
no evidence yield the declared unknown; per-column status distinguishes those
cases. Empty evidence never selects the first alphabet letter. Result metadata
retains winning count, evidence denominator and selected policy decisions.
Tie handling precedes the support-threshold check. Unknown/error tie policies
therefore remain observable below the threshold; ASCII-order tie resolution must
still pass the threshold before its selected symbol is accepted. Voting combines
all gap glyphs into `-` without changing the glyphs retained in analysis rows.

## Acceptance

Required evidence includes hand-countable identity/gap/missing/ambiguity cases,
mapping inverses, slice composition, source immutability, empty dimensions,
reverse/local origins and stable errors. Imported Stockholm/NEXUS metadata and
missing cells must survive read-only analysis; transforms must reject or report
loss. Real accepted pairwise/MSA results are adapted without recomputing kernels.
The isolated Rust checkpoint passed 32 internal tests, including 55,980 independent
explicit-pair statistics comparisons, 251,910 sorted-ballot consensus-policy
comparisons and 2,186 forward/reverse gapped rectangles with inverse maps and
composed slices. Real accepted Smith-Waterman and MAFFT results are adapted in
the tests. Strict standalone compilation, Clippy, rustdoc and formatting passed.
These are isolated-core checks, not final main-checkout verification or acceptance.
Seven additional external-caller tests passed against the isolated facade on
2026-10-18. They independently count the published fixture's cells/pairs/gap runs,
check reverse and local source origins, preserve unknown origins after removing
missing cells, compare full ordered annotation-loss reports, distinguish tie
resolution from minimum support, retain explicit zero-row widths and reject
dimension/range violations without changing the input. The
[external-caller tests](../../seq-bio/tests/alignment_analysis_public.rs),
the CLI checks (removed with the CLI),
[recipe checks](../../seq-bio/tests/alignment_analysis_recipes.rs) and
[`alignment_analysis.rs`](../../seq-bio/examples/alignment_analysis.rs) example
were promoted in the main checkout. The removed `alignment-stats`,
`alignment-consensus`, `alignment-slice` and `alignment-map` CLI operations shared
these computations and policies. The main tests/examples, documentation
checks and complete workspace harness output/terminal audit now support scoped
acceptance in the main validation record.
The isolated/staging results remain earlier checkpoints. Final source/evidence
review and the local commit accompany this acceptance; no accepted alignment kernel
or scientific qualification is changed.
