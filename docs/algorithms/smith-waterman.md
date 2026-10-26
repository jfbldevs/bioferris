# Smith–Waterman Local Alignment (`SW-local-linear-v1`)

## Status and exact target

Status: complete for `SW-local-linear-v1`; portfolio acceptance reconfirmed on
2026-08-12.

Phase 2 implements the scalar local-alignment specialization named
`SW-local-linear-v1`. It finds one maximum-scoring pair of subsequences, applies a
linear score once per gap symbol, retains one rolling score row plus the full
predecessor matrix, and returns one deterministic optimal traceback.

This document was written before implementation and remains the scientific and
behavioral contract. The implementation was reconstructed without consulting,
translating, or adapting source code from another alignment project.

## Primary reference

Smith, T. F., and Waterman, M. S. (1981). “Identification of Common Molecular
Subsequences.” *Journal of Molecular Biology* 147(1), 195–197.
[doi:10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).
PMID 7265238.

The publication defines zero boundary values, a recurrence containing zero, symbol
association, and deletions of arbitrary length, followed by traceback from the maximum
matrix element until a zero element. The authors' publication list also records this
paper at <https://dornsife.usc.edu/msw/research-papers/>.

Alphabet definitions remain those recorded for the common architecture in
`docs/algorithms/needleman-wunsch.md`.

## Objective

Given sequences

\[
X=x_1x_2\ldots x_m,\qquad Y=y_1y_2\ldots y_n,
\]

a substitution score \(s(x_i,y_j)\), and a non-positive linear gap score \(g\), find
contiguous ranges of both source sequences whose gapped alignment has maximum score.
Unaligned prefixes and suffixes contribute no score. The empty alignment with score
zero is always a candidate.

## Published recurrence and specialization

The paper assigns a positive deletion weight \(W_k\) to a deletion of length \(k\)
and gives the recurrence

\[
H(i,j)=\max\left\{
H(i-1,j-1)+s(x_i,y_j),
\max_{1\le k\le i}\{H(i-k,j)-W_k\},
\max_{1\le l\le j}\{H(i,j-l)-W_l\},
0
\right\}.
\]

For this phase, \(W_k=-kg\). Because a length-\(k\) gap then has exactly the sum of
\(k\) unit gap transitions, the recurrence specializes to

\[
H(i,j)=\max\begin{cases}
0,\\
H(i-1,j-1)+s(x_i,y_j) & \text{diagonal},\\
H(i-1,j)+g & \text{deletion},\\
H(i,j-1)+g & \text{insertion}.
\end{cases}
\]

Initialize

\[
H(i,0)=0\quad(0\le i\le m),\qquad
H(0,j)=0\quad(0\le j\le n).
\]

The optimal local score is the maximum value anywhere in \(H\). Traceback begins at
the selected maximum cell and terminates at the first zero predecessor.

## Reconstructed pseudocode

```text
function LOCAL_LINEAR_SW(X, Y, substitution, gap):
    rows    <- length(X) + 1
    columns <- length(Y) + 1
    allocate score[columns] initialized to zero
    allocate predecessor[rows, columns] initialized to STOP

    endpoint <- (0, 0)
    maximum  <- 0
    bounded  <- substitution declares valid inclusive bounds AND
                every possible candidate sum fits in i64

    for i from 1 through length(X):
        diagonal_base <- score[0]
        score[0] <- 0
        for j from 1 through length(Y):
            up_base <- score[j]
            diagonal <- add(diagonal_base, substitution(X[i], Y[j]), bounded)
            up       <- add(up_base, gap, bounded)
            left     <- add(score[j - 1], gap, bounded)

            if every candidate is at most zero:
                score[j] <- 0
                predecessor[i, j] <- STOP
            else:
                choose maximum using DIAGONAL > UP > LEFT for ties
                score[j] <- selected score
                predecessor[i, j] <- selected predecessor

            if score[j] is strictly greater than maximum:
                maximum <- score[j]
                endpoint <- (i, j)

            diagonal_base <- up_base

    traceback from endpoint until STOP, reverse operations,
    and return the resulting half-open source ranges
```

Ordinary `i64` addition is used only after the range proof succeeds. Otherwise each
candidate uses checked arithmetic. A positive overflow makes the mathematical cell
maximum unrepresentable and returns an error. Matrix bases are non-negative, so adding
an `i64` non-positive increment cannot underflow; any non-positive result loses to the
zero candidate.

## Data structures

- `Sequence<A>` stores validated alphabet-indexed symbols.
- A `Vec<i64>` of length \(n+1\) updates one rolling row of \(H\) in place.
- A row-major one-byte predecessor matrix stores stop, diagonal, up, or left.
- The endpoint stores the first row-major occurrence of the global maximum.
- An operation `Vec` is constructed backwards, reversed once, and moved into the
  result.
- `PairwiseAlignment` stores the score, selected half-open source ranges, and forward
  operations.

## Complexity and memory

For lengths \(m\) and \(n\):

- time: \(\Theta(mn)\);
- rolling score row: \(8(n+1)\) bytes;
- traceback matrix: \((m+1)(n+1)\) bytes;
- temporary traceback: at most \((m+n)\,\mathrm{sizeof}(operation)\) bytes;
- total asymptotic space: \(\Theta(mn)\).

Score-only execution, `SmithWaterman::score`, keeps only the rolling score row:
\(8(n+1)\) bytes. It is specified in [linear score-only execution](linear-score-only.md).

The published arbitrary-gap recurrence is more general. Its direct evaluation is not
the algorithm implemented in this phase.

## Parameters and public behavior

- Public constructor: `SmithWaterman::new(substitution, gap)`.
- Public execution: `SmithWaterman::align(first, second)`.
- Independent replay: `verify_local_alignment(...)`.
- `SubstitutionScorer<A>` supplies symbol-pair scores and may declare exact bounds.
- `LinearGap` accepts zero or negative scores and rejects positive scores.
- DNA, RNA, and protein use the same statically typed API.
- An optimum of zero returns score zero, empty operations, and canonical ranges
  `0..0` for both inputs.
- Only one optimal local alignment is returned.

## Tie-breaking

The publication does not specify how equal candidates or equal maximum cells are
resolved. This profile fixes both choices:

1. zero wins whenever the best candidate is zero, terminating a zero-score path;
2. among positive equal candidates: diagonal, deletion (`Up`), insertion (`Left`);
3. among equal positive matrix maxima, retain the first cell in row-major order.

These rules produce deterministic ranges and operations. Swapping inputs preserves
the optimal score for a symmetric scorer but need not preserve selected ranges or gap
placement.

## Elements not specified by the authors

- exact tie-breaking;
- one canonical representation when no positive similarity exists;
- integer width and overflow behavior;
- allocation-failure behavior;
- a stable operation-based API;
- ambiguity-symbol semantics;
- whether only the best or subsequent non-associated alignments must be returned.

This profile resolves each item explicitly. Subsequent-best local alignments mentioned
in the paper are outside this phase.

## Components

1. Existing typed alphabets, sequences, scoring, gaps, results, and errors.
2. Dimension, memory, and arithmetic-range preflight.
3. Scalar rolling-row local dynamic programming.
4. Deterministic maximum selection and traceback-to-zero.
5. Checked result construction and rendering.
6. Independent public local-alignment replay verifier.
7. Manual, exhaustive, property, and black-box validation.
8. Reproducible time and exact element-buffer memory benchmarks.

## Validation and acceptance criteria

The algorithm may be marked complete only if all of the following pass:

- manually calculated empty, identical, dissimilar, repetitive, and internal-match
  examples;
- alignments excluding non-homologous prefixes and suffixes;
- internal gaps and exclusion of penalized terminal gaps from selected local ranges;
- DNA, RNA, canonical protein, ambiguity, selenocysteine, and pyrrolysine input;
- zero-optimum canonicalization, invalid input, dimension overflow, and score overflow;
- independent replay of every operation and recomputation of every score;
- exhaustive enumeration of all alignments of every pair of binary-sequence
  substrings through length four;
- property tests for bounded ranges, score symmetry, reconstruction, and determinism;
- static black-box agreement with an explicitly configured reputable reference;
- reproducible time and exact element-buffer memory benchmarks;
- formatting, lint, test, documentation, dependency, source, and license gates.

Expected values in manual tests must not be generated by this implementation. The
exhaustive oracle must enumerate legal alignment paths independently and must not call
or reuse the dynamic-programming implementation.

## Differences from the 1981 program and example

- The original recurrence accepts arbitrary gap weights \(W_k\); this profile uses
  the explicit linear specialization \(W_k=-kg\).
- The paper's illustrative gap weight \(W_k=1+k/3\) is affine after scaling and is
  intentionally deferred to the Gotoh phase. It is not silently approximated here.
- Scores are signed 64-bit integers, so the paper's fractional example is not used as
  an exact golden test.
- The implementation returns only the best local alignment, not the paper's proposed
  sequence of subsequent non-associated alignments.
- It reports arithmetic and resource errors rather than overflowing or partially
  returning an alignment.

## Limitations

- Traceback memory is quadratic even though score storage is linear in the second
  sequence length.
- Only local alignment with a linear per-symbol gap score is implemented.
- Arbitrary, affine, and convex gap functions are not part of this phase.
- Named biological matrices such as PAM and BLOSUM are not bundled.
- Ambiguity symbols have literal semantics unless the caller supplies a custom matrix.
- Only one optimal traceback is returned.
- SIMD and parallel filling are intentionally deferred until scalar correctness and
  all acceptance gates are established.

Benchmark records are maintained in `docs/benchmarks/smith-waterman.md` and
`docs/benchmarks/smith-waterman-vs-biopython.md`.
