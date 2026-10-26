# Affine X-Drop Seed Extension

## Status and exact target

Status: implemented and validated on 2026-08-03 for profile
`XDROP-forward-affine-v1`.

The target profile is `XDROP-forward-affine-v1`. It reconstructs the dynamic-
programming X-drop extension described explicitly by Zhang, Schwartz, Wagner, and
Miller (2000), with the gapped-extension objective introduced by Gapped BLAST
(1997). It extends alignments from a fixed origin, selects the maximum-scoring
retained endpoint, uses affine gaps, and evaluates cells in anti-diagonal order.

The profile also defines a deterministic convenience operation around one forced
seed pair. Upstream and downstream extensions are independent applications of the
same forward primitive; the seed substitution score is included exactly once.

X-drop is a heuristic. A pruned path can never re-enter the search, so a larger
unrestricted semi-global score may exist. This phase does not implement BLAST word
hits, HSP construction, E-values, the Zhang greedy recurrence, Z-drop, fixed-width
banded or wavefront alignment, SIMD, or parallel execution.

## Scientific references

Primary algorithmic specification:

1. Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller. "A Greedy
   Algorithm for Aligning DNA Sequences." *Journal of Computational Biology*
   7(1-2):203-214, 2000.
   DOI: [10.1089/10665270050081478](https://doi.org/10.1089/10665270050081478).

   Reproduced material: Section 2 and Figure 2's dynamic-programming X-drop
   objective, anti-diagonal schedule, strict `T-X` pruning test, active-range
   contraction, fixed origin, open endpoint, and exhaustion rule. Sections 5.2 and
   5.3 establish that affine gaps and symbol-dependent substitution scores are valid
   dynamic-programming generalizations. The paper's specialized greedy algorithm is
   not reproduced by this profile.

Original and gapped BLAST context:

2. Stephen F. Altschul, Warren Gish, Webb Miller, Eugene W. Myers, and David J.
   Lipman. "Basic Local Alignment Search Tool." *Journal of Molecular Biology*
   215(3):403-410, 1990.
   DOI: [10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2).

   Used for the original score-drop extension concept and seed-and-extend context.

3. Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schaffer, Jinghui Zhang,
   Zheng Zhang, Webb Miller, and David J. Lipman. "Gapped BLAST and PSI-BLAST: a
   New Generation of Protein Database Search Programs." *Nucleic Acids Research*
   25(17):3389-3402, 1997.
   DOI: [10.1093/nar/25.17.3389](https://doi.org/10.1093/nar/25.17.3389).

   Reproduced material: cell-level gapped extension from a seed, independent forward
   and backward traversal, adaptive search-region exhaustion, affine gap scoring,
   and selection of the best score seen.

Affine recurrence:

4. Osamu Gotoh. "An Improved Algorithm for Matching Biological Sequences."
   *Journal of Molecular Biology* 162(3):705-708, 1982.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).

Independent black-box reference:

5. Biotite 1.7.1 `biotite.sequence.align.align_local_gapped`, BSD-3-Clause. Its
   documented downstream/upstream/both seed-extension modes, affine gap tuple, and
   strict X-drop behavior provide independent score observations. Biotite and Python
   remain in the sibling validation laboratory. No reference source is inspected,
   copied, translated, linked, or invoked by Cargo.

## Objective and coordinates

Let

```text
A = a_0 ... a_(m-1)
B = b_0 ... b_(n-1)
```

and let `(i,j)` identify prefixes `A[0..i]` and `B[0..j]`. The primitive begins at
`(0,0)` with score zero. It seeks the maximum-scoring path from the origin to any
retained cell. The selected endpoint is open:

```text
(i*, j*) = arg max retained H[i,j].
```

The empty extension at the origin is a valid result. Consequently the primitive's
score is never negative even when every non-empty extension is unfavorable.

There is no Smith-Waterman zero reset in interior cells. A positive-scoring path
cannot restart after leaving the origin.

## Affine recurrence

The workspace convention for a gap run of length `k >= 1` is

```text
gap(k) = open + (k - 1) * extend,
open <= extend <= 0.
```

This maps Zhang et al.'s Section 5.2 score `I*ind + G*a` to

```text
open   = a + ind
extend = ind.
```

It also maps the Gapped BLAST cost `Go + k*Ge` to signed scores

```text
open   = -(Go + Ge)
extend = -Ge.
```

For deletion state `D`, insertion state `I`, and best state `H`:

```text
D[i,j] = max(H[i-1,j] + open, D[i-1,j] + extend)
I[i,j] = max(H[i,j-1] + open, I[i,j-1] + extend)

H[i,j] = max(
    H[i-1,j-1] + s(a_(i-1), b_(j-1)),
    D[i,j],
    I[i,j]
)
```

with

```text
H[0,0] = 0
D[0,0] = I[0,0] = unreachable.
```

Boundary values are produced by the same recurrence. Thus `H[i,0]` is one deletion
run and `H[0,j]` is one insertion run until X-drop pruning removes it.

## Published X-drop pruning rule

Anti-diagonal `k` contains every real cell with `i+j=k`. Let `T_(k-1)` be the
maximum score retained through anti-diagonal `k-1`. Each candidate on anti-diagonal
`k` is first evaluated from retained predecessors. It is then discarded exactly when

```text
H[i,j] < T_(k-1) - X,
```

where `X >= 0`. Equivalently, a drop strictly greater than `X` is pruned. Equality
survives. The comparison is implemented without subtracting in the score type so an
extreme threshold cannot overflow.

Figure 2 updates `T'` while evaluating the current anti-diagonal but compares cells
against the prior `T`. This profile preserves that ordering: all cells on one
anti-diagonal use the same `T_(k-1)`, and only after that anti-diagonal is complete is
the global maximum advanced to `T_k`.

Zhang et al.'s core Figure 2 uses linear gaps and half-nodes to depend on one prior
anti-diagonal. Symbol-dependent scores can be odd and affine scoring has three states,
so this profile uses the mathematically direct regular-node recurrence above and
retains the two preceding anti-diagonals. It applies exactly the published pruning
predicate to regular DP cells.

### Affine-state pruning resolution

The papers describe retaining or discarding a path-graph cell by its optimal cell
score, but do not specify independent X-drop tests for Gotoh's internal states. This
profile resolves the ambiguity as follows:

- compute `D`, `I`, and `H` for the cell;
- test only `H` against the published `T-X` threshold;
- if the cell fails, discard `H`, `D`, and `I` together;
- if it survives, retain all three reachable affine states.

This follows Gapped BLAST's cell-level language and does not invent a second,
state-specific threshold.

## Active anti-diagonal region

Only retained cells can create successors. With regular nodes, a candidate cell on
anti-diagonal `k` can be reached by:

- insertion or deletion from anti-diagonal `k-1`;
- a diagonal substitution from anti-diagonal `k-2`.

If a retained frontier spans first-prefix indices `[L,U]`, insertion successors use
`[L,U]`, deletion successors use `[L+1,U+1]`, and diagonal successors use
`[L+1,U+1]` from the second-prior frontier. Each successor interval is clipped to
real cells before the surviving intervals are merged; an entirely out-of-matrix
interval cannot widen another legal interval. Legal coordinates satisfy

```text
max(0, k-n) <= i <= min(m, k).
```

Interior unreachable holes remain explicit; only unreachable cells at both edges are
trimmed. Traversal stops before `m+n` only when neither of the two preceding retained
frontiers can generate a real successor. Keeping two frontiers is essential: a direct
diagonal successor may remain possible after the intervening anti-diagonal becomes
empty.

## Seeded extension

`XDropSeed(first_index, second_index)` identifies one forced aligned symbol pair.
The seed indices must be inside their respective sequences. The seed operation is a
match or substitution and contributes

```text
s(first[first_index], second[second_index])
```

exactly once.

Directions are defined as:

- `Downstream`: run the forward primitive on suffixes strictly after the seed;
- `Upstream`: run the forward primitive on reversed prefixes strictly before the
  seed, then reverse the returned operations into biological order;
- `Both`: run both independent extensions and concatenate upstream operations, the
  forced seed operation, and downstream operations.

Adding the seed score as a constant does not change either side's X-drop comparisons.
For `Both`, the selected score is

```text
upstream_extension + seed_substitution + downstream_extension.
```

This is the two-way seed extension described by Gapped BLAST. It does not search for
or move the supplied seed.

## Own scalar pseudocode

```text
xdrop_forward(A, B, scorer, gap, X, trace):
    previous2 <- empty frontier
    previous1 <- frontier k=0 containing origin H=0
    best <- origin
    T <- 0

    for k <- 1 .. len(A)+len(B):
        candidate_range <- legal successor union(previous1, previous2)
        if candidate_range is empty:
            return best with FRONTIERS_EXHAUSTED

        current <- unreachable cells covering candidate_range
        T_next <- T

        for each first prefix i in candidate_range:
            j <- k-i
            compute affine D, I, and H from retained predecessors
            count this real candidate cell as evaluated

            if H is reachable and H >= T-X:
                retain the complete affine cell
                record deterministic traceback when requested
                update selected endpoint
                T_next <- max(T_next, H)

        trim unreachable cells from both frontier edges
        rotate previous2, previous1, current
        T <- T_next

    return best with MATRIX_EXHAUSTED
```

Traceback stores one packed byte for each evaluated candidate plus one metadata entry
per evaluated anti-diagonal. Score-only execution retains only three rolling affine
frontiers.

## Determinism and tie-breaking

The X-drop papers specify the score and pruning rule but do not prescribe traceback
or endpoint ties. This profile fixes them:

1. best-state tie: insertion, deletion, diagonal;
2. gap-state tie: extension before opening;
3. endpoint tie: greater `i+j`, then greater `j`, then greater `i`;
4. the seed pair is always represented by a diagonal match/substitution operation.

Repeated calls with equal validated inputs and configuration must be bit-for-bit
deterministic.

## Public API

- `XDrop<S>::new(substitution, gap, threshold)`;
- immutable scorer, gap, and threshold accessors;
- `score(first, second)` and `align(first, second)` for the published forward
  prefix-extension primitive;
- `score_from_seed(first, second, seed, direction)` and
  `align_from_seed(first, second, seed, direction)`;
- `XDropDirection::{Upstream, Downstream, Both}`;
- score, source-range, evaluated/retained-cell, anti-diagonal, maximum-frontier, and
  termination metadata;
- checked worst-case memory modeling for score and traceback execution.

All sequence APIs are generic over `Alphabet` and `SubstitutionScorer`, so DNA, RNA,
protein, ambiguity symbols, dense matrices, and custom scorers share one algorithm.

## Data structures and memory

- Each rolling frontier owns a contiguous first-prefix interval and affine states
  with explicit reachability. Capacity grows through checked reservations to the
  largest interval actually observed rather than the full theoretical width.
- Three frontiers retain anti-diagonals `k-2`, `k-1`, and `k`.
- Traceback execution appends a packed one-byte trace cell for every evaluated
  candidate and stores the first-prefix lower bound and trace offset for each
  anti-diagonal.
- Returned operations contain at most `m+n` entries for one forward extension.
- All dimensions, byte totals, score additions, and reservations are checked.

Let `W` be the maximum candidate width and `Q` the number of evaluated real candidate
cells:

```text
score auxiliary space  Theta(W)
alignment history      Theta(Q)
returned operations    O(m+n)
```

Worst case has `W = min(m,n)+1` and `Q = (m+1)(n+1)`, so X-drop retains quadratic
worst-case time and traceback space. Favorable extensions terminate with a much
smaller adaptive region.

## Complexity

```text
fill time       Theta(Q)
traceback time  O(first_end + second_end)
score space     Theta(W)
alignment space Theta(Q + first_end + second_end)
```

The two-way seed convenience operation is the sum of its independent upstream and
downstream traversals. It never evaluates a cell crossing through the seed from one
side to the other because the seed is forced.

## Elements omitted or ambiguous in the publications

| Issue | Resolution in `XDROP-forward-affine-v1` |
|---|---|
| The explicit Figure 2 recurrence is linear-gap, while affine gaps are sketched later. | Use the published cell pruning rule around the standard three-state affine recurrence. |
| Figure 2's half-node optimization requires divisible substitution scores. | Use regular nodes and two prior anti-diagonals, preserving generic integer substitution matrices. |
| The prose says score drops by `X`; Figure 2 prunes with strict `< T-X`. | Follow the pseudocode: equality survives. |
| `T'` is updated during an anti-diagonal, but pruning uses `T`. | Every cell uses the best score from completed earlier anti-diagonals. |
| Affine substates have no published independent pruning rule. | Retain or discard the entire cell according to `H`; preserve all states in a retained cell. |
| Traceback ties are unspecified. | Insertion, deletion, diagonal; gap extension before opening. |
| Maximum endpoint ties are unspecified. | Prefer farther total consumption, then greater second-sequence consumption. |
| The seed-selection heuristic is BLAST-pipeline-specific. | Require an explicit seed; do not search for one. |
| Forward/backward result composition is described but not formalized. | Run independent reversed-prefix and suffix extensions and score the seed once. |
| Empty input has no seed. | Forward extension is valid and returns the origin; seeded APIs reject out-of-bounds seeds. |
| Score infinity and overflow are abstract. | Store reachability independently, discard negatively underflowed losing candidates, and report positive optimal overflow. |

## Independently testable components

1. Candidate ranges include insertion, deletion, and second-prior diagonal successors.
2. Strict threshold equality survives and strict excess is pruned.
3. All cells on one anti-diagonal compare against the prior completed maximum.
4. An empty intermediate anti-diagonal still permits a second-prior diagonal
   successor.
5. Wide/no-effective-drop execution matches an independent full affine semi-global
   matrix and exhaustive path enumeration.
6. Every returned operation path independently replays to its score and ranges.
7. Zero-drop identity traversal reaches the complete diagonal even when both
   intervening gap frontiers are pruned; symmetric scoring preserves the best score
   under input transposition.
8. Upstream, downstream, and both-direction seed scores equal independent primitive
   compositions.
9. Static fixtures match Biotite 1.7.1 score observations without Cargo invoking it.
10. Benchmarks consume score, traceback, cell counts, and exact modeled memory.

## Acceptance criteria

This phase may be marked complete only when all of the following hold:

- hand-calculated identity, mismatch, initial/internal/terminal affine gap, empty,
  repetitive, completely different, strict-threshold, and early-exhaustion cases pass;
- DNA, RNA, and protein use the same public implementation;
- an independent non-X-drop oracle confirms the selected score whenever pruning is
  disabled in effect;
- a separate literal X-drop oracle agrees on every evaluated/retained cell and
  endpoint for exhaustive short inputs;
- property tests cover deterministic repetition, range bounds, score replay,
  transposition symmetry, zero-drop identity, seeded composition, and oracle
  equivalence;
- immutable Biotite 1.7.1 fixtures agree on every compatible recorded score;
- Python, Biotite, fixture generation, and live reference calls stay outside Cargo
  and the product workspace;
- Criterion records score, alignment, threshold scaling, early termination, seeded
  directions, scientific protein data, and modeled memory;
- mirrored score and traceback benchmarks use identical sequences, scoring, X-drop,
  seeds, warm-up, repetitions, and consumed outputs;
- formatting, locked workspace tests, examples, Clippy, rustdoc with warnings denied,
  and all `cargo-deny` checks pass;
- production code contains no placeholder, fallback aligner, external engine, FFI,
  subprocess, `unsafe`, or code derived from another implementation.

## Known profile limitations

- X-drop can prune the unrestricted optimum and is not an exact local aligner.
- A large threshold can approach quadratic work and memory for traceback.
- A small threshold can stop at a short low-scoring bridge that a biologically useful
  alignment would cross.
- Forward extension has a fixed origin and open endpoint; it is not global, overlap,
  or arbitrary free-end alignment.
- Seeded extension requires a correct caller-supplied seed and always forces that
  symbol pair into the result.
- Scalar execution is the correctness baseline. SIMD and parallel wavefronts require
  separate post-correctness review and differential tests.
- Biotite chooses a shorter endpoint for one fixture with two score-15 endpoints;
  this profile deliberately chooses the farther endpoint according to its published
  tie policy. All 24 reference scores and every unambiguous range agree.

## Completion evidence

The implementation and public API are in `align-fast/src/xdrop.rs` and are exported
from `align-fast`. There is no fallback aligner, FFI, subprocess, external engine,
`unsafe`, placeholder, or reference-derived source.

Correctness evidence:

- 3 internal invariant tests, 11 manual/exhaustive integration tests, 4 property
  groups at 192 cases each, and 1 static-reference test pass;
- 3,844 short binary/configuration combinations agree cell-for-cell with an
  independent dense literal X-drop oracle;
- 675 tiny input/configuration combinations agree between the X-drop result, an
  unrestricted affine matrix, and enumeration of every operation path;
- every returned path is replayed independently for score, operation labels,
  affine gap runs, consumption, and ranges;
- empty, identity, completely different, repetitive, strict-threshold, early-stop,
  initial/internal/terminal gap, deterministic tie, overflow, all seed directions,
  DNA, RNA, protein, and ambiguity cases are explicit;
- 24 pinned Biotite 1.7.1 black-box cases have zero score disagreements and zero
  unambiguous-range disagreements. The sole equal-score endpoint-tie difference is
  named and required by the test.

Reproducibility and performance evidence:

- the isolated fixture regenerates byte-for-byte from a SHA-256-locked CPython 3.13
  environment; its copied workspace SHA-256 is
  `d839e6c77396d366d8a316f61192c1c9ce59fde55fc71c84b93662809cea5c10`;
- the mirrored seven-repeat protocol shows Rust faster at every length and in both
  operations, from 1.83x to 823.69x for score and 5.83x to 849.55x for traceback;
- Criterion covers score, traceback, thresholds, upstream/downstream/both seeds,
  curated UniProt globins, and checked worst-case memory;
- complete records are
  `x-drop.md` and
  `x-drop-vs-biotite.md`.

Final gates on 2026-08-03 passed: rustfmt check; Clippy across the locked workspace,
all targets, and all features with warnings denied; locked workspace tests across all
targets/features; rustdoc and doctests with warnings denied; the executable X-drop
example; and cargo-deny advisories, bans, licenses, and sources. The isolated
generator's `--check` mode and a fresh hash-locked wheel installation also passed.
