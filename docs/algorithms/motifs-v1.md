# Nucleotide motifs and scoring matrices

## Contents

- [Scope and current state](#scope-and-current-state)
- [IUPAC pattern semantics](#iupac-pattern-semantics)
- [Matrix model](#matrix-model)
- [Scoring and reverse strands](#scoring-and-reverse-strands)
- [Limits and interfaces](#limits-and-interfaces)
- [Verification](#verification)

[API index](../README.md) · Implementation plan

## Scope and current state

Stage 3 includes both IUPAC motif matching and count/probability/log-odds matrices.
The public Rust facade and prelude helpers are implemented and verified; the CLI
operations and strict JSON requests that Stage 3 also accepted were removed with
the `bioferris` CLI (decision 029). Stage 3 is accepted for this
documented matching/matrix scope in the
main validation record.
Earlier isolated numerical and staged-integration results remain separate
checkpoints; acceptance is not universal native/API parity.
Only existing typed DNA/RNA alphabets are selected. Primary references:

- [IUPAC-IUB nucleotide symbols](https://iubmb.qmul.ac.uk/misc/naseq.html).
- [Biopython 1.88 sequence utilities](https://biopython.org/docs/1.88/api/Bio.SeqUtils.html).
- [Biopython 1.88 motif tutorial](https://biopython.org/docs/1.88/Tutorial/chapter_motifs.html).
- [Public matrix API](https://biopython.org/docs/1.88/api/Bio.motifs.matrix.html).

## IUPAC pattern semantics

Compile each pattern symbol into its A/C/G/T set; RNA U uses the T bit while
retaining the existing RNA validator. Empty patterns fail. An input shorter than
the pattern has no windows and succeeds with no matches.

Input ambiguity is explicit:

- `Definite` (default): every input possibility is accepted by the pattern,
  `input_mask & pattern_mask == input_mask`.
- `Possible`: at least one possibility agrees, `input_mask & pattern_mask != 0`.
- `Reject`: any ambiguous source symbol is an error with original-input offset.

Search all overlapping windows. Select forward/reverse/both strands explicitly,
defaulting to both. Reverse matching uses a reverse-complemented pattern against
the original source. Return ascending zero-based half-open original-input ranges
and an explicit strand, sorted by start then forward before reverse. A palindromic
pattern on both strands yields two labeled matches at the same range. Pattern
ambiguity is never expanded into a potentially exponential list of sequences.

For `AGRYN` and pattern `R`, forward definite matches are 0..1, 1..2 and 2..3;
possible also includes 4..5. Reverse definite matches 3..4; reverse possible also
includes 4..5. `AAAA` and forward `AA` match 0..2, 1..3 and 2..4.

## Matrix model

Use explicit column arrays in A/C/G/T order (A/C/G/U for RNA). Counts and
pseudocounts must be finite and nonnegative. Matrices have nonempty equal widths;
weighted counts and unequal column totals are permitted. Construct counts from
canonical equal-length instances or explicit columns. Ambiguous instances fail.

Normalize each column after adding its per-base pseudocounts:
`P(base,column) = (count + pseudocount) / sum_base(count + pseudocount)`.
A zero-total column is an error. Zero original counts with positive pseudocount
totals are valid. Background weights must be finite and strictly positive; retain
their normalized distribution, defaulting to uniform. No zero background or
implicit pseudocount is inferred.

Column and background totals use an independently implemented Neumaier
compensated sum in alphabet order. This reproduces the improved floating-point
`sum` behavior in the pinned CPython 3.12 reference. Adding the four adjusted
counts sequentially can change the last binary64 bits and is not equivalent.
Log-odds cells are evaluated as `ln(P / background) / ln(2)` in binary64, rather
than using a platform `log2` call. P=0 produces negative infinity. Core
models retain that value explicitly; the removed CLI encoded it in JSON as the
string `"-inf"` rather than silently converting nonfinite values into null.

An independent exact example uses instances AC, AC, AC, CA and pseudocount 1:

| Base | Counts | Probabilities | Uniform-background log2 scores |
|---|---|---|---|
| A | 3,1 | 0.5,0.25 | 1,0 |
| C | 1,3 | 0.25,0.5 | 0,1 |
| G | 0,0 | 0.125,0.125 | -1,-1 |
| T/U | 0,0 | 0.125,0.125 | -1,-1 |

## Scoring and reverse strands

Score each complete concrete source window against the built scoring matrix.
PSSM scanning rejects ambiguous input explicitly; IUPAC set matching does not
define a probabilistic score for an ambiguous source. The reference documents
float32 scan results. Each window adds score cells sequentially in binary64 and
converts the final sum once to binary32. A finite, binary32-representable
threshold is also rounded once to binary32 before inclusive `score >= threshold`
comparison. These rules are verified by direct public-API observations, including
cancellation and thresholds between adjacent binary32 values.

The observed Biopython 1.88/NumPy 2.3.3 search raises a zero-dimensional-array
error when there is exactly one source window. BioFerris intentionally supports
that valid case; its score is compared with the reference's `calculate` result,
not with the failing search. Native NaN scores for ambiguous windows are likewise
not reproduced: this profile returns an explicit ambiguity error instead.

Reverse scoring reverses/complements the **built PSSM**, or its equivalent lookup:
`S_rc(base,j) = S(complement(base),width-1-j)`. Do not rebuild reversed counts
with an unchanged asymmetric background. With a uniform one-column probability
and background A=.5,C=.25,G=.125,T=.125, scores are A=-1,C=0,G=1,T=1;
reverse scoring on original A gives the original T score, 1.

For the two-column example and source ACGTAC, forward scores are
2,-1,-2,-1,2; reverse scores are -2,-1,2,-1,-2. An inclusive threshold 2
retains forward0..2, reverse2..4 and forward4..6. Core/JSON coordinates remain
original-input half-open ranges; Python-negative reverse indices are not used.

## Limits and interfaces

Search budgets bound total evaluated windows, symbol comparisons and returned
matches across selected strands. Exceeding a budget fails without truncating a
successful result. Zero result budget permits only no-hit input. Matrix arithmetic
must reject overflow/NaN, and scans must not allocate an unbounded score grid.

The Rust facade exposes typed `IupacMotif<A>`, search options/results,
count/weight/scoring matrices and typed errors. A sealed nucleotide alphabet
trait prevents inventing semantics for arbitrary custom alphabets. Core objects
remain immutable; short constructors and prelude exports support manual programs.
The removed `motif-find`, `motif-matrix` and `motif-scan` CLI operations provided
strict JSON access through the same core.

## Verification

Independent set-expansion tests cover every IUPAC symbol combination and both
strands, overlaps, palindromes, empty/short sources and limits. Matrix tests use
the exact rational/integer example, asymmetric backgrounds, zero-probability
cells and independently calculated probabilities/scores. Numerical observations
must pin Biopython, NumPy, interpreter and platform; exact matrix/scan comparisons
retain their precision and do not silently replace float32 with float64 semantics.
The [360-request retained corpus](../../seq-bio/tests/fixtures/motifs-provenance.md)
matches 5,244 PWM cells, 5,244 PSSM cells, 28,158 forward/reverse scan scores,
9,127 hits and 18 zero-column errors exactly in the standalone Rust harness.
At the isolated checkpoint, 18 internal tests also covered 450 IUPAC
input/pattern policy cases and 225 two-symbol reverse-complement cases. Six
additional external-caller tests passed against the isolated facade on
2026-10-18. They pin the retained corpus
digest/runtime metadata and exercise prelude constructors, literal matrix results,
RNA/DNA parity, asymmetric-background reverse scoring, negative infinity, budgets
and immutable inputs. This is historical isolated evidence, not a report of final
main-checkout test success. The public facade wiring and promoted
[`motif_analysis.rs`](../../seq-bio/examples/motif_analysis.rs),
[external-caller tests](../../seq-bio/tests/motifs.rs) and
[executable recipe checks](../../seq-bio/tests/motif_recipes.rs) now exist in the
main checkout. The [manual guide](../api/motifs-v1.md) records complete
observed workflows and their evidence limits. The main-checkout tests, recipes,
documentation checks and complete workspace harness output/terminal audit now
support scoped acceptance in the
main validation record.
Staged checks remain a preparatory checkpoint; final source/evidence review and
the local commit accompany this acceptance.
