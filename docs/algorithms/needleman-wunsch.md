# Needleman–Wunsch Global Alignment (`NW-global-linear-v1`)

## Status and exact target

Status: complete for `NW-global-linear-v1`; portfolio acceptance reconfirmed on
2026-08-12.

Phase 1 implements the scalar global alignment specialization named
`NW-global-linear-v1`. It consumes both input sequences completely, applies a linear
score once per gap symbol, retains one rolling score row plus the full predecessor
matrix, and returns one deterministic optimal traceback.

The implementation was reconstructed from the publication and the specification in
this document. No source code from another alignment project was consulted or adapted.

## Primary reference

Needleman, S. B., and Wunsch, C. D. (1970). “A general method applicable to the
search for similarities in the amino acid sequence of two proteins.” *Journal of
Molecular Biology* 48(3), 443–453.
[doi:10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836(70)90057-4).
PMID 5420325.

Alphabet references:

- IUPAC–IUB Commission on Biochemical Nomenclature (1972). “A One-Letter Notation
  for Amino Acid Sequences.”
  [doi:10.1351/pac197231040639](https://doi.org/10.1351/pac197231040639).
- Cornish-Bowden, A. (1985). “Nomenclature for incompletely specified bases in
  nucleic acid sequences: recommendations 1984.” *Nucleic Acids Research* 13(9),
  3021–3030. [doi:10.1093/nar/13.9.3021](https://doi.org/10.1093/nar/13.9.3021).
- INSDC Feature Table Definitions, amino-acid abbreviations and symbols,
  <https://www.insdc.org/submitting-standards/feature-table/>.

## Objective

Given sequences

\[
X=x_1x_2\ldots x_m,\qquad Y=y_1y_2\ldots y_n,
\]

a substitution score \(s(x_i,y_j)\), and a non-positive linear gap score \(g\), find
a global alignment of maximum total score. Each alignment column is exactly one of a
symbol pair, a symbol from \(X\) against a gap, or a gap against a symbol from \(Y\).
Gap-gap columns are invalid.

## Fidelity boundary

The 1970 paper presents a general maximum-match path method. It permits arbitrary cell
weights and states that a gap penalty may depend on gap size or direction. Its examples
also discuss paths that begin or end on array borders and a randomization procedure for
estimating match significance.

The paper therefore does not uniquely prescribe all conventions of the modern
three-predecessor recurrence when scores may be negative. This implementation makes
the following specialization explicit:

- strict global boundaries require complete consumption of both sequences;
- a gap of length \(k\) scores \(k g\), not one constant score per gap run;
- only unit diagonal, vertical, and horizontal transitions are needed;
- statistical significance by randomized sequences is not part of alignment;
- the output is one optimal alignment rather than the set of all optimal alignments.

This is not represented as a literal reproduction of the paper's FORTRAN experiment.

## Equations

Let \(F(i,j)\) be the maximum score for globally aligning prefixes \(X_{1..i}\) and
\(Y_{1..j}\). Initialize:

\[
F(0,0)=0,
\]

\[
F(i,0)=F(i-1,0)+g \quad (1\le i\le m),
\]

\[
F(0,j)=F(0,j-1)+g \quad (1\le j\le n).
\]

For \(i>0\) and \(j>0\):

\[
F(i,j)=\max\begin{cases}
F(i-1,j-1)+s(x_i,y_j) & \text{diagonal},\\
F(i-1,j)+g & \text{deletion},\\
F(i,j-1)+g & \text{insertion}.
\end{cases}
\]

The optimal global score is \(F(m,n)\).

## Reconstructed pseudocode

```text
function GLOBAL_LINEAR_NW(X, Y, substitution, gap):
    rows    <- length(X) + 1
    columns <- length(Y) + 1
    allocate score[columns]
    allocate predecessor[rows, columns]

    bounded <- substitution declares valid inclusive bounds AND
               every possible path sum through rows + columns fits in i64

    score[0] <- 0
    for j from 1 through length(Y):
        score[j] <- add(score[j - 1], gap, bounded)
        predecessor[0, j] <- LEFT

    for i from 1 through length(X):
        diagonal_base <- score[0]
        score[0] <- add(diagonal_base, gap, bounded)
        predecessor[i, 0] <- UP
        for j from 1 through length(Y):
            up_base <- score[j]
            diagonal <- add(diagonal_base, substitution(X[i], Y[j]), bounded)
            up       <- add(up_base, gap, bounded)
            left     <- add(score[j - 1], gap, bounded)
            choose maximum using DIAGONAL > UP > LEFT for ties
            in checked mode, reject a positive overflow or an
                all-negative-overflow candidate set
            score[j] <- representable maximum
            predecessor[i, j] <- chosen predecessor
            diagonal_base <- up_base

    i <- length(X)
    j <- length(Y)
    operations <- empty list
    while i > 0 or j > 0:
        follow predecessor[i, j]
        append its operation and decrement consumed index or indices
    reverse operations
    return score[length(Y)] and operations
```

`add(value, increment, bounded)` uses ordinary `i64` addition only after the range
proof succeeds. Otherwise it applies the checked candidate rules documented below.

## Data structures

- `Sequence<A>` stores validated alphabet-indexed symbols.
- A `Vec<i64>` of length \(n+1\) updates one rolling row of \(F\) in place.
- A row-major one-byte predecessor matrix stores start, diagonal, up, or left.
- An operation `Vec` is constructed backwards, reversed once, and moved into the
  result without another element-buffer allocation.
- `PairwiseAlignment` stores the score, complete source ranges, and operations.

## Complexity and memory

For lengths \(m\) and \(n\):

- time: \(\Theta(mn)\);
- rolling score row: \(8(n+1)\) bytes;
- traceback matrix: \((m+1)(n+1)\) bytes;
- temporary traceback: at most \((m+n)\,\mathrm{sizeof}(operation)\) bytes;
- total asymptotic space: \(\Theta(mn)\).

Score-only execution, `NeedlemanWunsch::score`, keeps only the rolling score row:
\(8(n+1)\) bytes. It is specified in [linear score-only execution](linear-score-only.md).

Allocator metadata, allocator rounding, and input storage are not included in the
exact element-buffer estimate exposed by `memory_requirement`.

## Parameters and public behavior

- `SubstitutionScorer<A>` supplies every symbol-pair score and may declare inclusive
  score bounds.
- `MatchMismatch` uses literal symbol equality, including ambiguity symbols.
- `DenseSubstitutionMatrix<A>` accepts one complete row-major score table.
- `LinearGap` accepts zero or negative scores and rejects positive scores.
- Scores are `i64`. When declared scorer bounds prove that every possible path score
  fits, one preflight check removes overflow branches from the scalar cell loop. If no
  bounds are declared or the proof fails, candidate additions use the fully checked
  path. A positively overflowing candidate must dominate every representable score and
  returns an error. A negatively overflowing candidate is discarded if another
  candidate is representable; if all three underflow, the cell returns an error.
- Matrix dimensions and buffer sizes are checked before allocation.

## Tie-breaking

Equal candidates are selected in fixed order:

1. diagonal;
2. deletion, consuming the first sequence and placing a gap in the second;
3. insertion, placing a gap in the first sequence and consuming the second.

This makes repeated runs bit-for-bit deterministic. Swapping inputs preserves the
optimal score for a symmetric scorer but need not preserve the same visual gap
placement.

## Components

1. Typed alphabet and sequence validation in `align-core`.
2. Scoring and linear gap validation in `align-core`.
3. Dimension and memory preflight.
4. Arithmetic-range preflight and scalar rolling-row fill.
5. Deterministic traceback.
6. Checked result construction and rendering.
7. Independent public replay verifier.

## Validation and acceptance criteria

The algorithm is accepted only if all of the following pass:

- manually calculated empty, identical, dissimilar, repetitive, and classic examples;
- leading, internal, and terminal gap cases;
- DNA, RNA, canonical protein, ambiguity, selenocysteine, and pyrrolysine input;
- invalid-input, invalid-scoring, dimension-overflow, and score-overflow checks;
- independent replay of every operation and recomputation of every score;
- exhaustive enumeration of every alignment for all binary sequence pairs through
  length four;
- property tests for validity, score symmetry, reconstruction, and determinism;
- static black-box agreement with Biopython 1.87 for scores and, where bounded, the
  complete optimal traceback and metric set;
- workspace formatting, lint, test, documentation, dependency, and license gates;
- reproducible time and exact element-buffer memory benchmarks.

The internal Criterion record is in `docs/benchmarks/needleman-wunsch.md`. A separate
mirrored end-to-end protocol measured this implementation at 1.42×–4.75× the speed of
Biopython 1.87 for the tested lengths; scope and limitations are documented in
`docs/benchmarks/needleman-wunsch-vs-biopython.md`.

Expected values in manual tests are calculated independently. The exhaustive oracle
recursively enumerates all legal column paths and does not call or reuse the dynamic
programming implementation.

## Differences from the original program

- The implementation supports DNA, RNA, and protein through one typed architecture.
- It uses a strict global boundary condition.
- Gap cost is linear per emitted gap symbol.
- It accepts arbitrary signed integer substitution matrices rather than the paper's
  codon-correspondence experimental tables.
- It does not reproduce randomized significance estimation.
- It reports arithmetic and resource errors instead of overflowing or partially
  returning an alignment.

## Limitations

- Memory is quadratic because traceback retains the complete matrix.
- Only global alignment and linear gaps are implemented.
- Named biological matrices such as PAM and BLOSUM are not bundled in this phase.
- Ambiguity symbols have literal semantics unless the caller supplies a custom matrix.
- Only one optimal traceback is returned.
- SIMD, parallel filling, affine gaps, local/end-free modes, and accelerated
  algorithms are intentionally outside this phase.
- Cargo never invokes an external program. Black-box conformance currently covers one
  explicitly configured Biopython 1.87 profile through static fixtures. The large
  scientific case compares only the optimal score; complete optimal tracebacks and
  metrics are compared for the bounded short cases.
