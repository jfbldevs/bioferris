# Affine Local and Explicit End-Gap Alignment

## Reproduction targets and status

This document specifies two completed scalar profiles:

- `AFFINE-ENDGAP-v1`: the Gotoh affine recurrence with four explicit free-end
  permissions; and
- `SMITH-WATERMAN-GOTOH-local-affine-v1`: Smith-Waterman local restart applied
  to the three-state Gotoh recurrence.

Status: complete and validated on 2026-08-05.

These profiles close the current pairwise boundary-mode gap without relabeling a
linear algorithm as affine. They share the existing typed alphabet, substitution,
gap, result, and error model for protein, DNA, and RNA. The already completed
`GOTOH-global-affine-v1`, `ENDGAP-linear-v1`, and `SW-local-linear-v1` profiles
remain independent reduction oracles.

## Scientific references

Primary sources:

1. Temple F. Smith and Michael S. Waterman. "Identification of Common Molecular
   Subsequences." *Journal of Molecular Biology* 147(1), 195-197, 1981.
   DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).
   The zero-restart recurrence and maximum-over-all-cells local endpoint are
   normative for the local profile.
2. Osamu Gotoh. "An Improved Algorithm for Matching Biological Sequences."
   *Journal of Molecular Biology* 162(3), 705-708, 1982.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
   Equations (1), (4), and (5), the affine restriction, alternative boundary
   initialization, and straight-ahead gap traceback guidance are normative.
3. Saul B. Needleman and Christian D. Wunsch. "A General Method Applicable to
   the Search for Similarities in the Amino Acid Sequence of Two Proteins."
   *Journal of Molecular Biology* 48(3), 443-453, 1970.
   DOI: [10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836(70)90057-4).
   Prefix dynamic programming and complete-sequence global boundaries are the
   base case recovered when no end is free.
4. Peter H. Sellers. "On the Theory and Computation of Evolutionary Distances."
   *SIAM Journal on Applied Mathematics* 26(4), 787-793, 1974.
   DOI: [10.1137/0126070](https://doi.org/10.1137/0126070).
5. Michael S. Waterman, Temple F. Smith, and William A. Beyer. "Some Biological
   Sequence Metrics." *Advances in Mathematics* 20(3), 367-387, 1976.
   DOI: [10.1016/0001-8708(76)90202-4](https://doi.org/10.1016/0001-8708(76)90202-4).

Later clarification and comparison sources:

6. Stephen F. Altschul and Bruce W. Erickson. "Optimal Sequence Alignment Using
   Affine Gap Costs." *Bulletin of Mathematical Biology* 48(5-6), 603-616,
   1986. DOI:
   [10.1007/BF02462326](https://doi.org/10.1007/BF02462326).
   This is used to audit affine-state traceback completeness, not as a source of
   implementation code.
7. Jeff Daily. "Parasail: SIMD C Library for Global, Semi-Global, and Local
   Pairwise Sequence Alignments." *BMC Bioinformatics* 17, 81, 2016.
   DOI: [10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z).
   Its explicit global/semiglobal/local terminology supports the public naming.
8. Biopython 1.87 `Bio.Align.PairwiseAligner` documentation. This is used only
   for isolated black-box observations after independent Rust oracles exist.

There is no universal definition of "semiglobal" or "overlap" across tools.
This specification therefore defines four booleans and treats names only as
presets. No reference implementation source may be inspected or adapted.

## Common score model

For typed sequences `A=a_1...a_M` and `B=b_1...b_N`, substitution scorer `s`,
and affine gap parameters:

```text
gap(k) = open + (k - 1) * extend,  k >= 1
open <= extend <= 0
```

`E[i,j]` ends with `a_i` aligned to a gap (deletion), `F[i,j]` ends with
`b_j` aligned to a gap (insertion), and `H[i,j]` is the best state:

```text
E[i,j] = max(H[i-1,j] + open, E[i-1,j] + extend)
F[i,j] = max(H[i,j-1] + open, F[i,j-1] + extend)
H[i,j] = max(H[i-1,j-1] + s(a_i,b_j), E[i,j], F[i,j])
```

Opening from `H` permits adjacent opposite gap runs, as the frozen global Gotoh
profile does. No transition is silently removed.

## `AFFINE-ENDGAP-v1`

### Objective

`FreeEndGaps(FS, FE, SS, SE)` controls whether a prefix or suffix of the first or
second sequence may be omitted without score. The returned half-open ranges
identify the scored core. Free overhangs are not emitted as zero-score operations.

Presets retain their existing algebraic meanings:

```text
global                       = (false, false, false, false)
semiglobal                   = (true,  true,  true,  true )
fit_first_in_second          = (false, false, true,  true )
fit_second_in_first          = (true,  true,  false, false)
first_suffix_second_prefix   = (true,  false, false, true )
second_suffix_first_prefix   = (false, true,  true,  false)
```

The suffix-prefix presets are the two directed overlap orientations. This
profile does not require a nonempty overlap: with non-positive evidence, an empty
boundary-to-boundary core is a valid mathematical optimum and must remain visible.

### Boundary states

Use an explicit unreachable state, never a numeric sentinel. At the origin:

```text
H[0,0] = 0
E[0,0] = F[0,0] = unreachable
```

For the first-sequence boundary:

```text
if FS:
    H[i,0] = 0;      E[i,0] = unreachable
else:
    H[i,0] = E[i,0] = gap(i)
F[i,0] = unreachable
```

For the second-sequence boundary:

```text
if SS:
    H[0,j] = 0;      F[0,j] = unreachable
else:
    H[0,j] = F[0,j] = gap(j)
E[0,j] = unreachable
```

Endpoint candidates are:

```text
always: H[M,N]
if FE:  H[i,N] for 0 <= i <= M
if SE:  H[M,j] for 0 <= j <= N
```

The endpoint must lie on the last row or last column. Both suffixes cannot be
discarded from an arbitrary interior cell; that would be local alignment.
Traceback may stop at `(i,0)` only when `FS`, or at `(0,j)` only when `SS`.

### Deterministic ties

1. A gap extension wins an equal-score gap opening.
2. `H` ties prefer diagonal, deletion, then insertion.
3. Endpoint ties prefer greater `row + column`, then greater `row`.
4. Traceback emits one forward operation stream after reversing the recovered
   predecessor path.

Rules 1 and 2 reproduce the completed global Gotoh profile. Rule 3 makes free-end
endpoint selection deterministic instead of leaving it to matrix iteration order.

## `SMITH-WATERMAN-GOTOH-local-affine-v1`

### Objective and recurrence

The local profile selects arbitrary ranges from both sequences and adds the
Smith-Waterman restart candidate:

```text
H[0,j] = H[i,0] = 0
E[0,j], F[0,j], E[i,0], F[i,0] = unreachable

E[i,j] = max(H[i-1,j] + open, E[i-1,j] + extend)
F[i,j] = max(H[i,j-1] + open, F[i,j-1] + extend)
H[i,j] = max(0, H[i-1,j-1] + s(a_i,b_j), E[i,j], F[i,j])
```

Gap-state values are not clamped independently; only `H` is restarted. The
endpoint is the maximum `H[i,j]` over the whole matrix. Traceback stops at the
`Start` predecessor of a zero `H` state.

### Deterministic ties and zero optimum

1. Gap extension wins an equal opening.
2. For a positive cell, diagonal wins deletion, which wins insertion.
3. A maximum not greater than zero selects `Start`, never a zero-score operation.
4. Equal positive endpoint maxima retain their first row-major occurrence.
5. A zero optimum returns score zero, empty operations, and canonical ranges
   `0..0` for both sequences.

Scores and endpoint selection reduce exactly to `SW-local-linear-v1` when
`open == extend`. Equal-score operation order can differ because this affine profile
retains extension-before-opening state ties while linear traceback has no gap state.

## Pseudocode

```text
AFFINE_ALIGN(A, B, scorer, gap, boundary):
    checked_allocate rolling H row, rolling E row, E-valid flags
    checked_allocate packed traceback matrix
    initialize H/E/F boundaries from boundary
    endpoint <- required/global or local origin

    for i in 1..M:
        diagonal_base <- previous H[i-1,0]
        initialize H[i,0], E[i,0], and row-local F state
        for j in 1..N:
            old_H <- H[i-1,j]
            E <- choose_extend_then_open(old_H + open, old_E + extend)
            F <- choose_extend_then_open(H[i,j-1] + open, old_F + extend)
            diagonal <- diagonal_base + s(a_i,b_j)

            if local:
                H <- choose(Start=0, diagonal, E, F)
                observe every H endpoint using strict improvement
            else:
                H <- choose(diagonal, E, F)
                observe last-column endpoint

            store H predecessor and both gap predecessors
            diagonal_base <- old_H

    if end-gap:
        observe H[M,N], configured last-column candidates,
        and configured last-row candidates

    traceback H/E/F states until the permitted boundary or local Start
    reverse operations
    construct checked PairwiseAlignment ranges
    independently replay score, ranges, operations, and affine gap runs
```

The bounded arithmetic path is enabled only when scorer bounds prove every
candidate fits `Score`. Otherwise every addition is checked, negative-overflowed
unreachable candidates are discarded, and positive overflow is a typed error.

## Data structures and memory

- one rolling `H` row;
- one rolling deletion-state `E` row;
- one deletion-validity row for checked arithmetic;
- one scalar row-local insertion state `F`;
- one packed byte per matrix cell containing the selected `H` predecessor and
  the open/extend choice for both gaps;
- one fallibly reserved operation vector;
- explicit endpoint coordinates and score;
- typed `Sequence<A>`, `AffineGap`, `SubstitutionScorer<A>`, and
  `PairwiseAlignment` values.
- optional caller-owned `AffineWorkspace` shared sequentially with global Gotoh;
  it retains score, checked-validity, and packed-trace capacities while each returned
  operation vector remains independent.

No full score matrix, unsafe code, FFI, external runtime, or new production
dependency is permitted.

## Complexity

For `M=|A|` and `N=|B|`:

- time: `O(MN)`;
- rolling score storage: `O(N)`;
- traceback storage: `O(MN)` bytes;
- returned operations: `O(M+N)` worst case.

All dimensions, element counts, byte counts, and score additions are checked.

## Public API contract

The phase adds:

```text
AffineEndGapAligner<S>::new(substitution, gap, FreeEndGaps)
AffineLocalAligner<S>::new(substitution, gap)
AffineEndGapAligner<S>::align_with_workspace(..., &mut AffineWorkspace)
AffineLocalAligner<S>::align_with_workspace(..., &mut AffineWorkspace)
verify_affine_end_gap_alignment(...)
verify_local_affine_alignment(...)
```

Both aligners expose scoring accessors, exact memory estimates, stable profile
identifiers, and one generic implementation for protein, DNA, and RNA. Existing
`Gotoh` and `SmithWaterman` APIs remain source-compatible.
The ordinary `align` methods allocate an ephemeral workspace and delegate to the
same implementation, so the recurrence, boundaries, endpoint ordering, and
traceback are not duplicated.

## Published omissions and explicit resolutions

- **Semiglobal terminology:** resolved through four booleans; labels are only
  documented presets.
- **Empty overlap:** permitted and represented honestly; a future nonempty-overlap
  policy must be a separate profile.
- **Affine traceback ties:** extension before opening, then diagonal/deletion/
  insertion, as frozen above. When `open == extend`, optimal scores and endpoint
  selection reduce to the linear profiles, but a repeated affine state can choose
  a different equal-score operation order than a memoryless linear traceback.
- **Local affine composition:** explicitly named Smith-Waterman-Gotoh because no
  claim is made that either paper alone specifies this complete API and tie order.
- **Terminal affine gaps:** a scored terminal run uses the same
  `open + (k-1)*extend` function as an internal run. A configured free overhang is
  excluded from the scored range instead.
- **Adjacent opposite gaps:** retained because opening transitions originate at
  `H`; they are not silently forbidden.
- **Ambiguity symbols:** interpreted only by the supplied scorer.
- **Extreme arithmetic:** no `Score::MIN` sentinel and no saturating arithmetic.

## Independently testable components

1. Boundary initialization for all 16 `FreeEndGaps` configurations.
2. Affine gap-state open/extend transitions and tie packing.
3. Endpoint selection for last-row, last-column, and local all-cell domains.
4. Three-state traceback and range recovery.
5. Independent affine end-gap verifier.
6. Independent affine local verifier with positive-prefix checks.
7. Reduction to global Gotoh.
8. Reduction of the local specialization to Smith-Waterman when `open == extend`.
9. Checked overflow and exact memory accounting.
10. Static black-box fixtures and process benchmarks.

## Acceptance criteria

The profiles may be marked complete only when:

- manual protein, DNA, and RNA cases cover empty inputs, identity, complete
  disagreement, repeated symbols, leading/internal/trailing affine gaps,
  fit-both-directions, both directed overlaps, all-free semiglobal, local restart,
  and adjacent opposite gaps;
- all 16 free-end configurations and local alignment agree with an independent
  exhaustive path oracle on bounded binary sequences;
- `FreeEndGaps::global()` agrees exactly with the existing Gotoh score, ranges,
  and operation stream;
- free-end modes agree with an independent exhaustive affine recurrence, and the
  local `open == extend` specialization agrees with Smith-Waterman; affine
  extension-before-opening ties are explicitly tested and documented;
- property tests cover score symmetry under sequence/mode mirroring, determinism,
  range consumption, score upper bounds, replay, and operation tampering;
- independently calculated overflow and memory cases pass;
- reviewed Biopython 1.87 fixtures compare exact scores and compatible normalized
  ranges/operations for DNA, RNA, and protein without executing Python in Cargo;
- an isolated mirrored in-process comparison against Biopython and Criterion
  benchmarks include local, fitting, semiglobal scientific protein, length, and
  global-reduction workloads; overlap shares the same fill kernel and has separate
  exhaustive correctness coverage;
- documentation records every opaque-reference tie difference rather than fitting
  the implementation to it;
- debug and release tests, examples, formatting, workspace Clippy and rustdoc with
  warnings denied, dependency policy, source hygiene, and `git diff --check` pass;
- production code contains no placeholder, mock, unchecked fallback, unsafe code,
  external engine, subprocess, FFI, or unimplemented route.

No expected result may be generated by either new Rust aligner.

The initial pre-implementation draft requested fresh-process time/RSS measurements.
That protocol was revised explicitly before completion: for sub-millisecond pairwise
calls, process startup measures CPython and loader cost rather than the alignment
kernel, while RSS is dominated by the language runtime and allocator. The accepted
protocol instead times repeated complete calls in isolated runtimes and publishes an
exact Rust owned-buffer model. The decision and both harnesses are retained in
`docs/benchmarks/affine-alignment-modes-vs-biopython.md`.

## Implementation and validation record

Production implementation:

- `align-pairwise/src/affine_modes.rs` contains the shared scalar rolling-row fill,
  explicit local/four-end boundaries, packed three-state traceback, deterministic
  endpoints, checked arithmetic, and public aligners;
- `align-pairwise/src/verification.rs` independently replays selected source ranges,
  operation kinds, affine gap runs, local positive prefixes, and reported scores;
- no new production dependency, unsafe block, FFI, subprocess, external engine, or
  runtime sidecar was introduced.

Validation evidence:

- 15 manual acceptance tests cover all required alphabets and edge classes, gap-run
  placement, every empty-side permission, arithmetic extremes, custom matrices,
  deterministic repetition, adjacent opposite gaps, and deliberate result tampering;
- an independent complete-path enumerator checks all 16 four-end modes on every
  binary pair through length three: 3,600 mode/input combinations;
- the same independent enumerator checks every local substring path for every binary
  pair through length four: 961 input pairs;
- six property tests run 256 generated cases each by default for replay,
  determinism, symmetry, exact global Gotoh reduction, and linear score/endpoint
  reductions; the recorded release stress run used 5,000 cases per property (30,000
  generated cases total);
- 27 static Biopython 1.87 cases cover all 16 end configurations, local affine
  alignment, DNA, RNA, protein, internal affine runs, ambiguity, repetition, and two
  scientific score-only cases; 25 bounded cases retain all 31 distinct compatible
  optimal reference paths, with zero disagreements;
- independent score and traceback verifiers are invoked by native, exhaustive,
  property, and black-box tests.
- workspace-specific tests reuse one scratch object across changing shapes, modes,
  orientations, all supported alphabets, bounded/checked arithmetic, and a failed
  overflow call; 192 property cases compare every complete result with a fresh call.

Fixture SHA-256 is
`334557f3e00f6ca3964f179e461b4ff12c2916d2b8cee34b194e79c57465cf18`.
The external generator SHA-256 is
`dfe257d5e201a9f629e4de91e19d7dd161e1857b9f3d84594b9e628ad8ec8a96`,
and its case specification SHA-256 is
`cbda5eec88fa9364b941957814d933a7181f14ee626ea381438d4be8475cb9c8`.

Measured time and exact memory are recorded in
`docs/benchmarks/affine-alignment-modes.md`. Under the mirrored complete-call
protocol, Rust was 1.21×–2.80× faster than Biopython for affine fitting and
2.56×–3.91× faster for affine local alignment. These performance observations are
not used as correctness evidence.

## Remaining limitations

- single-threaded; the exact scalar baseline now has an optional runtime-dispatched
  AVX2 i16/i32 backend, but no non-x86 SIMD implementation;
- one symmetric affine gap function for both directions;
- one canonical optimum rather than all ties;
- no dual-affine, convex, PSSM-specific, GPU, or linear-space traceback;
- the scalar/SIMD selector is calibrated by shape and numeric safety, but there is no
  cross-algorithm regime selector;
- empty overlap is allowed by the overlap presets.
