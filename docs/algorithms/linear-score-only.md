# Exact Linear Score-Only Execution

## Contents

- [Scope](#scope)
- [Recurrence and result](#recurrence-and-result)
- [Arithmetic and errors](#arithmetic-and-errors)
- [Complexity and memory](#complexity-and-memory)
- [Validation](#validation)

[Needleman–Wunsch](needleman-wunsch.md) · [Smith–Waterman](smith-waterman.md) ·
Validation record

Profile identifier: `LINEAR-SCORE-ONLY-v1`, added on 2026-10-22.

## Scope

This profile computes only the optimum integer score of the existing linear-gap
objectives: global `NW-global-linear-v1` and local Smith–Waterman. It is an
execution profile, not a new biological objective, and covers every alphabet and
substitution scorer those kernels accept. It returns the same score as the
corresponding complete alignment, without the predecessor matrix, operations or
aligned ranges.

Entry points: `NeedlemanWunsch::score` and `SmithWaterman::score` in
`align-pairwise`, used by `Pairwise::needleman_wunsch().score(...)` and
`Pairwise::smith_waterman().score(...)` in the `align` facade. Before this
profile, the facade's `score()` built the complete traceback; it now uses these
kernels. The traceback kernels themselves are unchanged.

## Recurrence and result

The score-only kernels evaluate the recurrence of the
[global](needleman-wunsch.md#equations) and
[local](smith-waterman.md#published-recurrence-and-specialization) specifications
cell by cell in the same row-major order, over one rolling row of `n + 1`
scores for a second sequence of length `n`:

- global: the boundary row and column hold multiples of the gap score, and the
  result is the last cell of the last row;
- local: every cell is at least `0`, and the result is the largest cell value,
  `0` when no positive cell exists.

The maximum of the three candidates is the same number whichever tie-breaking
order a traceback would use, so the score equals `align(first, second)?.score()`.

## Arithmetic and errors

The same preflight selects the arithmetic: when every possible path score fits in
the 64-bit `Score`, unchecked addition is exact; otherwise each candidate is
checked as in the traceback kernel. A global cell whose mathematical maximum is
not representable, and any overflowing local candidate, return
`AlignmentError::ScoreOverflow` naming the same row and column as `align`.
Dimension overflow is reported only for the score row and the path length. The
traceback matrix is never allocated, so an input whose matrix `align` cannot
allocate can still be scored.

## Complexity and memory

For lengths `m` and `n`:

- time: `Θ(mn)`;
- algorithm-owned memory: one score row of `8(n + 1)` bytes, versus
  `(m + 1)(n + 1)` additional traceback bytes for complete alignment.

For two sequences of 31,623 bases, the row takes 252,992 bytes where the
traceback matrix of complete alignment takes about 1 GB.

## Validation

`align-pairwise/tests/linear_score_only.rs` compares both kernels with their
traceback kernels on random DNA inputs and scoring parameters, through both the
unchecked and checked arithmetic paths, on overflowing inputs and on a
2,000 × 2,500 pair. The validation record
lists the checks and a laboratory comparison on the publication benchmark's
pairwise inputs.
