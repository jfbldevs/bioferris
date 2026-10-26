# Exact global affine WFA

## Scope

Status: complete for the exact global affine profile; portfolio acceptance
reconfirmed on 2026-08-12.

`WavefrontAffine` implements exact end-to-end gap-affine Wavefront Alignment. It
minimizes positive penalties over two complete sequences and returns either the exact
penalty summary or a deterministic complete traceback.

The retained gap convention is

```text
gap(k) = gap_opening + k * gap_extension,  k > 0
```

where mismatch and extension penalties are strictly positive and the opening penalty
is nonnegative. The conventional maximizing alignment score is the negative of the
penalty.

## Scientific basis

The implementation follows the exact affine WFA recurrence described by Marco-Sola
et al., *Fast gap-affine pairwise alignment using the wavefront algorithm*,
Bioinformatics 2021, DOI
[`10.1093/bioinformatics/btaa777`](https://doi.org/10.1093/bioinformatics/btaa777).

For penalty `s` and diagonal `k = h - v`, the algorithm stores furthest-reaching
horizontal offsets for match/mismatch (`M`), insertion (`I`), and deletion (`D`)
states. Predecessor wavefronts are selected by their exact penalty increments:

```text
I_s,k = max(M_(s-o-e),k-1 + 1, I_(s-e),k-1 + 1)
D_s,k = max(M_(s-o-e),k+1,     D_(s-e),k+1)
X_s,k = max(M_(s-x),k + 1, I_s,k, D_s,k)
M_s,k = extend_exact_matches(X_s,k)
```

Here `x`, `o`, and `e` are mismatch, opening, and extension penalties. Coordinates
are checked before every transition and exact-match extension. The search terminates
when the target diagonal reaches both sequence ends.

## Public API

- `WavefrontAffinePenalties::new(mismatch, opening, extension)` validates the model;
- `WavefrontAffine::penalty` and `penalty_with_workspace` return only the exact
  optimum through the lean score path;
- `WavefrontAffine::score` returns `WavefrontScore` without traceback;
- `WavefrontAffine::score_with_workspace` reuses rolling and recycled wavefront
  storage;
- `WavefrontAffine::align` returns `WavefrontAlignment` and its
  `PairwiseAlignment`;
- `WavefrontAffine::align_with_workspace` returns the same independently owned
  result while reusing internal buffers;
- `WavefrontAffine::traceback_in_workspace` returns a compact
  `WavefrontTracebackView` whose `WavefrontTraceRun` iterator is in forward
  alignment order;
- `WavefrontAffine::align_in_workspace` returns a `WavefrontAlignmentView` borrowed
  from the workspace, avoiding the final operation-vector copy for repeated calls;
- `WavefrontWorkspace` reports retained bytes and supports explicit release;
- `WavefrontAffine::memory_requirement` returns a checked upper-bound model;
- score summaries expose wavefront, offset, match-extension, width, and modeled-byte
  diagnostics.

The implementation uses sparse weighted-score scheduling, a 64-bit ring scheduler
for small penalty increments, rolling score storage for score-only calls, recycled
wavefront allocations, and full predecessor history only when traceback is
requested. Validated sequences lazily expose compact alphabet indices so equality
extension compares bytes without changing literal symbol semantics. Optimization
of scheduling, equality extension, active intervals, or allocation does not alter
wavefront order or tie-breaking.

For unequal inputs with at most 8,192 symbols on the shorter side, the exact
completion bound also checks whether the inputs differ by one contiguous gap. It
extends their equal prefix and compares the remaining suffix after applying the
net length shift. A successful proof supplies the real `o + k*e` path cost and
can restrict the search to width-one wavefronts. A failed proof falls back to a
valid substitution-plus-terminal-gap alignment. Both alternatives are actual
alignment costs, so pruning remains admissible and the algorithm remains exact.

Compact traceback is run-length encoded inside `WavefrontWorkspace`. A view
borrows that storage until the next mutable workspace use and exposes penalty,
alignment score, input lengths, column count, run count, and a forward iterator.
Materialization remains explicit for consumers that require a per-column vector.

## Correctness contract

For every successful call:

1. the alignment consumes both inputs completely;
2. each operation replays against the corresponding symbols;
3. gap runs use `o + k*e` exactly;
4. the replayed penalty equals the returned penalty;
5. a dense affine dynamic-programming oracle yields the same optimum;
6. empty inputs, identity, substitutions, terminal/internal gaps, long gaps,
   repetitive ties, and penalty-boundary cases remain defined;
7. fresh, reused-owned, compact, and borrowed-workspace APIs return identical
   penalties, and expanding compact runs yields the exact same operation sequence.

Unit, exhaustive, property, and WFA2 2.3.6 black-box tests enforce this contract.
The pinned fixture contains 19 affine DNA/RNA/protein observations and does not
require WFA2 at test time.

## Validation and benchmarks

```text
cargo test -p align-fast --test wavefront
cargo test -p align-fast --test wavefront_properties
cargo test -p align-fast --test wfa2_reference
cargo bench -p align-benchmarks --bench wfa
cargo run --release --locked -p align-benchmarks --example wfa_validation_matrix
```

See the benchmark protocol and
the WFA2 comparison. The closure matrix and all
four retained process rounds are available in
`wfa-gap-asymmetry-divergence-closure-2026-08-11.json`.
