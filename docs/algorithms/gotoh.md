# Gotoh Global Alignment (`GOTOH-global-affine-v1`)

## Status and exact target

Status: complete as of 2026-08-03 for the exact profile defined here.

The target is the scalar, deterministic, global pairwise algorithm obtained from
equations (1), (4), and (5) of Gotoh (1982), specialized to one affine gap function
shared by deletions, insertions, and terminal gaps.

The public score convention is maximization:

```text
gap(k) = open + (k - 1) * extend,  k >= 1
```

`open <= extend <= 0`. The original paper minimizes distances and writes
`w_k = u*k + v`, with `u >= 0` and `v >= 0`. The exact parameter mapping is:

```text
open   = -(u + v)
extend = -u
```

The target does not include local, semiglobal, overlap, free-end, piecewise-affine,
SIMD, or parallel variants.

## Primary reference

Osamu Gotoh. "An Improved Algorithm for Matching Biological Sequences."
*Journal of Molecular Biology* 162(3), 705-708, 1982.

- DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9)
- PMID: [7166760](https://pubmed.ncbi.nlm.nih.gov/7166760/)
- Publication date: 15 December 1982
- Exact reproduced material: the affine restriction `w_k = u*k + v`, the
  `D/P/Q` recurrences, global boundary conditions, linear-space score observation,
  and the traceback tie guidance in Figure 1.

No source code from an alignment implementation was inspected or used. Biopython
1.87 is used only as an independently executed black-box oracle after the scalar
implementation and independent tests exist. Its documented affine convention is
`open + (k - 1) * extend`, which matches this public API directly.

## Objective

For complete typed sequences

```text
A = a_1 ... a_M
B = b_1 ... b_N
```

and substitution score `s(a_i, b_j)`, return one maximum-score global alignment.
Every symbol in both inputs must be consumed exactly once. A contiguous run of `k`
deletions or insertions contributes `gap(k)` once, rather than charging an identical
linear score to each gap column.

The same implementation operates on DNA, RNA, and protein `Sequence<A>` values.
Ambiguity symbols retain the semantics of the supplied substitution scorer.

## Published recurrence and score conversion

Gotoh minimizes the distance matrix `D`:

```text
D[m,n] = min(
    D[m-1,n-1] + d(a_m,b_n),
    P[m,n],
    Q[m,n]
)

P[m,n] = min(D[m-1,n] + w_1, P[m-1,n] + u)
Q[m,n] = min(D[m,n-1] + w_1, Q[m,n-1] + u)
```

where `w_k = u*k + v`. Negating distances produces the maximization recurrence used
here. `H` corresponds to `-D`, `E` to `-P`, and `F` to `-Q`:

```text
E[i,j] = max(H[i-1,j] + open, E[i-1,j] + extend)
F[i,j] = max(H[i,j-1] + open, F[i,j-1] + extend)
H[i,j] = max(
    H[i-1,j-1] + s(a_i,b_j),
    E[i,j],
    F[i,j]
)
```

`E` ends in a deletion: `a_i` is aligned to a gap. `F` ends in an insertion:
`b_j` is aligned to a gap.

The global boundaries are:

```text
H[0,0] = 0
H[i,0] = gap(i)  for i > 0
H[0,j] = gap(j)  for j > 0

E[i,0] = gap(i)  for i > 0
F[0,j] = gap(j)  for j > 0
E[0,j] and F[i,0] are unreachable for positive j and i
```

Opening `E` from `H` or `F` from `H` preserves the alignment space of the published
`D/P/Q` recurrence. In particular, a deletion run may be followed immediately by an
insertion run if that path is optimal. This behavior is not silently replaced by a
restricted three-state recurrence that requires an intervening aligned symbol.

## Reconstructed pseudocode

```text
INPUT: A[1..M], B[1..N], substitution s, affine scores open and extend
OUTPUT: one optimal score and global operation path

allocate H_row[0..N]
allocate E_row[0..N]
allocate trace[0..M, 0..N]

H_row[0] <- 0
for j <- 1..N:
    H_row[j] <- open + (j - 1) * extend
    trace[0,j].H <- INSERTION
    trace[0,j].F <- OPEN if j = 1 else EXTEND

for i <- 1..M:
    old_diagonal <- H_row[0]
    H_row[0] <- open + (i - 1) * extend
    E_row[0] <- H_row[0]
    trace[i,0].H <- DELETION
    trace[i,0].E <- OPEN if i = 1 else EXTEND

    F <- unreachable
    for j <- 1..N:
        old_H <- H_row[j]

        E_open   <- old_H + open
        E_extend <- E_row[j] + extend, if i > 1
        E_row[j] <- maximum(E_open, E_extend), preferring EXTEND on a tie

        F_open   <- H_row[j-1] + open
        F_extend <- F + extend, if j > 1
        F <- maximum(F_open, F_extend), preferring EXTEND on a tie

        diagonal <- old_diagonal + s(A[i], B[j])
        H_row[j] <- maximum(diagonal, E_row[j], F),
                    preferring DIAGONAL, then DELETION, then INSERTION
        record all selected canonical predecessors in trace[i,j]
        old_diagonal <- old_H

score <- H_row[N]

state <- H
(i,j) <- (M,N)
while i > 0 or j > 0:
    if state = H:
        DIAGONAL: emit match/substitution; i--; j--
        DELETION: state <- E
        INSERTION: state <- F
    else if state = E:
        emit deletion; i--
        state <- E for EXTEND, otherwise H
    else:
        emit insertion; j--
        state <- F for EXTEND, otherwise H

reverse emitted operations
return score and operations over ranges 0..M and 0..N
```

The checked implementation evaluates the same recurrence with exact overflow
detection. A bounded fast path uses ordinary addition only after substitution bounds,
gap scores, and the maximum path length prove that every partial path score fits in
`i64`.

## Data structures

- `AffineGap`: validated `open` and `extend` scores using the public length-one-open
  convention.
- `H_row[N+1]`: one rolling row of best scores.
- `E_row[N+1]`: one rolling row of deletion-state scores.
- `F`: the insertion-state score for the current cell.
- Packed one-byte traceback cell:
  - selected `H` predecessor;
  - whether `E` opened or extended;
  - whether `F` opened or extended.
- `Vec<AlignmentOperation>`: the canonical forward traceback.
- Optional caller-owned `AffineWorkspace`: retains the two rolling rows, checked-path
  validity row, and packed traceback capacity across sequential calls. Returned
  operations remain independently owned.

No score matrix is retained after a row is completed. The traceback matrix is kept
because this phase returns an alignment, not only an optimal score.

## Complexity and memory

For lengths `M` and `N`:

- Time: `Theta(MN)`.
- Rolling score memory: `2 * (N + 1) * sizeof(Score)` bytes.
- Traceback memory: `(M + 1) * (N + 1)` bytes.
- Maximum operation capacity: `(M + N) * sizeof(AlignmentOperation)` bytes.
- Overall memory while constructing a returned alignment: `Theta(MN)`.

This reproduces Gotoh's essential `MN` recurrence. The paper's score-only
linear-memory observation is preserved in the rolling score layout, while deterministic
traceback necessarily adds the direction matrix.

## Parameters and public behavior

- `substitution`: any `SubstitutionScorer<A>`.
- `gap.open_score()`: score for the first column of a gap run.
- `gap.extension_score()`: score for every later column in that same run.
- `gap.score_for_length(k)`: checked score of a complete positive-length gap.
- Inputs: immutable typed sequences with one common alphabet.
- Output: `PairwiseAlignment` consuming both complete sequences.
- Empty/empty result: score zero, empty operations, ranges `0..0` and `0..0`.
- Empty/non-empty result: one gap run, not repeated opening penalties.
- Execution is deterministic and has no global mutable state.
- `align_with_workspace` is result-identical to `align`; it changes only scratch
  lifetime. One mutable workspace can be reused across Gotoh, affine local/end-gap,
  protein, DNA, and RNA calls, but cannot be shared concurrently.
- `AffineWorkspace::reserve(M, N)` fallibly moves worst-case scratch reservation out
  of a hot loop. Capacity may exceed the active-call memory estimate after a larger
  prior alignment.

Positive gap rewards are rejected. To retain the original `u,v >= 0` restriction,
an opening score less severe than its extension score (`open > extend`) is also
rejected.

## Tie-breaking

The paper retains multiple direction bits and tells traceback to continue straight
when possible at a branch point. This library returns one canonical optimum:

1. Within `E` or `F`, extension wins a tie with opening. This is the explicit
   straight-ahead rule.
2. Within `H`, diagonal wins, then deletion (`E`), then insertion (`F`).
3. Boundaries have their only valid gap direction.

The second rule is not specified by Gotoh. It is an explicit deterministic extension
consistent with the existing pairwise API. Scores remain independent of this choice.

## Elements not specified or ambiguous

- The publication uses distance minimization; the API uses similarity maximization.
  The algebraic sign conversion is explicit above.
- Modern libraries disagree whether a length-one gap costs `open` or
  `open + extend`. This API uses `open`; the mapping to the paper is explicit.
- The paper describes a two-pass construction of a multi-bit direction matrix. This
  implementation records both gap-state decisions during the single score pass.
- The publication does not define one total order across diagonal, deletion, and
  insertion ties.
- Exact integer overflow, allocation failure, typed ambiguity symbols, and a stable
  operation-based public result are engineering requirements absent from the paper.

## Components

1. Validated affine scoring primitive and errors in `align-core`.
2. Scalar `H/E/F` recurrence and packed traceback in `align-pairwise`.
3. Independent affine score replay grouped by contiguous operation runs.
4. Manual unit tests for boundaries, internal and terminal gaps, ties, and all three
   biological alphabets.
5. Exhaustive small-input oracle that enumerates complete operation paths without
   using the Gotoh recurrence.
6. Property tests for score replay, symmetry under symmetric scoring, structural
   validity, determinism, and reduction to linear gaps when `open == extend`.
7. Static Biopython 1.87 black-box fixtures generated in an isolated sibling project.
8. Criterion time/memory benchmarks and a mirrored cross-language benchmark.
9. Fresh-versus-reused complete-result equivalence tests and matched AB/BA workspace
   benchmarks against the same release worker.

## Validation and acceptance criteria

This variant may be marked complete only when all of the following hold:

- the implementation matches the equations and public convention above;
- manually calculated DNA, RNA, and protein cases pass;
- empty, identical, fully different, repetitive, initial, internal, and terminal-gap
  cases pass;
- independent replay confirms every reported score and complete-range consumption;
- exhaustive enumeration confirms every optimal score in the bounded test domain;
- property tests are deterministic and use expected values not generated by Gotoh;
- the `open == extend` specialization agrees with the completed linear global
  algorithm;
- static black-box cases agree with pinned Biopython on score and, when completely
  enumerable, contain the returned traceback among the reference optima;
- benchmarks record time and modeled memory on reproducible inputs;
- a mirrored isolated benchmark records Rust/Biopython intervals and relative speed;
- formatting, linting, tests, documentation, dependency policy, and license gates pass;
- the public API contains no placeholders, unsafe code, external engine, FFI, or
  subprocess path;
- limitations and differences remain documented.

Completion evidence:

- all 91 workspace tests pass, including manual, property, regression, and scientific
  cases for the three completed pairwise profiles;
- the independent exhaustive oracle agrees on all 961 binary sequence pairs through
  length four;
- 11 pinned Biopython cases agree on every score; the 10 bounded cases record 13
  reference tracebacks, and each Rust canonical result belongs to its case's set;
- the scientific P02185/P68871 profile independently reproduces score `-53`;
- Criterion, memory, and mirrored Biopython measurements are recorded under
  `docs/benchmarks`;
- formatting, Clippy with warnings denied, documentation, doctests, offline checks,
  and all cargo-deny advisories/bans/licenses/sources checks pass.

## Differences from the 1982 program description

- Scores are maximized and represented exactly as `i64`; the paper minimizes
  non-negative distances.
- One packed direction byte is written in one pass instead of building a multi-path
  three-bit matrix with two orientation-swapped passes.
- Only one canonical optimal alignment is returned.
- The public gap parameters use the modern length-one-open convention.
- The implementation reports allocation, dimension, structure, and score overflow
  failures rather than relying on machine-specific numeric behavior.
- No source or behavior of the historical program is assumed beyond the publication.

## Limitations

- Global mode only.
- One symmetric affine function for both sequences and all positions.
- One optimal traceback; alternative optima are not enumerated.
- Full `Theta(MN)` traceback memory.
- No asymmetric, free-end, local, piecewise-affine, convex, or arbitrary gap costs.
- The scalar implementation remains the correctness baseline. A runtime-dispatched,
  exact AVX2 i16/i32 backend is available through `AffineWorkspace` and documented
  separately in `affine-simd.md`; there is no parallel execution.
