# Exact Affine Score-Only Execution

Status: implementation and bounded validation passed on 2026-08-07 UTC;
release-candidate stratification and T3 remain open.

Profile identifier: `AFFINE-SCORE-ONLY-v1`.

## Scope

This profile computes only the optimum integer score of the existing exact affine
objectives. It covers global Gotoh, Smith--Waterman--Gotoh local alignment, and all
16 combinations of the four independent free-end permissions. It supports protein,
DNA, and RNA through the same typed sequence and substitution-scoring contracts as
complete affine alignment.

This is an execution profile, not a new biological objective. It must return the
same score as the corresponding complete alignment while omitting traceback
construction, operation allocation, aligned ranges, and rendered rows.

## Scientific references

1. Osamu Gotoh. "An improved algorithm for matching biological sequences."
   *Journal of Molecular Biology* 162(3), 1982, 705--708.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
2. Thomas F. Smith and Michael S. Waterman. "Identification of common molecular
   subsequences." *Journal of Molecular Biology* 147(1), 1981, 195--197.
   DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).
3. Jeff Daily, Ananth Kalyanaraman, Sriram Krishnamoorthy, and Bin Ren. "On the
   Impact of Widening Vector Registers on Sequence Alignment." *45th International
   Conference on Parallel Processing*, 2016, 506--515.
   DOI: [10.1109/ICPP.2016.65](https://doi.org/10.1109/ICPP.2016.65).
4. Jeff Daily. "Parasail: SIMD C library for global, semi-global, and local pairwise
   sequence alignments." *BMC Bioinformatics* 17, 81, 2016.
   DOI: [10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z).

Only papers, equations, published descriptions, and independently observable
black-box behavior informed this implementation. Source code from Parasail,
Biopython, or another alignment implementation was not an implementation source.

## Objective and equations

Let `H[i,j]` be the best score for the two prefixes, `D[i,j]` a path whose last
column deletes a symbol from the second sequence, and `I[i,j]` a path whose last
column inserts a symbol into the first sequence. Let `o <= e <= 0` be the library's
length-one gap-open and later-symbol extension scores.

```text
D[i,j] = max(H[i-1,j] + o, D[i-1,j] + e)
I[i,j] = max(H[i,j-1] + o, I[i,j-1] + e)
H[i,j] = max(H[i-1,j-1] + substitution(first[i], second[j]),
             D[i,j],
             I[i,j])
```

Local mode adds `0` to the final maximum. Global and free-end boundary values are
exactly those specified in `gotoh.md` and `affine-alignment-modes.md`.

For explicit free ends, the score is selected from:

```text
H[m,n]                         always
max H[i,n], 0 <= i <= m        when first_end is free
max H[m,j], 0 <= j <= n        when second_end is free
```

Start permissions alter only the first column or first row. Local alignment returns
the maximum of zero and every interior `H` value.

## Scalar pseudocode

```text
initialize H[0..n] from second_start
initialize rolling D[0..n] as unavailable
best_local = 0
best_last_column = H[n]

for i in 1..m:
    diagonal = old H[0]
    update H[0] from first_start
    I = unavailable

    for j in 1..n:
        old_h = H[j]
        D[j] = max(old_h + open, D[j] + extend when available)
        I = max(H[j-1] + open, I + extend when available)
        H[j] = max(diagonal + substitution(i,j), D[j], I)
        if local: H[j] = max(0, H[j])
        best_local = max(best_local, H[j])
        diagonal = old_h

    best_last_column = max(best_last_column, H[n])

select the objective score from H, best_local, and best_last_column
```

The bounded path uses ordinary signed-64-bit additions only after a conservative
whole-path proof. A scorer without valid bounds uses checked additions plus explicit
gap-state representability flags. Positive overflow is an error. A negatively
overflowing candidate is unavailable, matching the complete-alignment contract.

## SIMD pseudocode

The implementation reuses `AFFINE-SIMD-PREFIX-SCAN-v1`. One second-sequence column
is processed at a time and the first-sequence axis occupies vector lanes.

```text
prepare alphabet-indexed substitution profile
prepare each second-symbol profile offset
initialize one H column and one I column
precompute every prefix-scan extension delta and inactive sentinel

for each symbol in second:
    for each full padded vector block in first:
        diagonal = shifted(old_H) + profile
        I = max(old_H + open, old_I + extend)
        Q = max(diagonal, I, zero when local)
        D = affine_prefix_max(Q, carry, open, extend)
        H = max(Q, D)
        update objective-score accumulators

return selected optimum score
```

The score-only specialization does not compute trace codes or write a traceback
matrix. Runtime ISA dispatch, signed-16/signed-32 preflight proofs, and strict SIMD
errors remain identical to the complete backend. Automatic crossover gates are
calibrated independently because removing traceback changes both scalar and SIMD
costs.

The global path is selected at dispatch time and monomorphized separately from
local and free-end execution. It therefore omits endpoint scans, local restarts,
intermediate terminal-row writes, and free-end branches. Its rolling columns store
only rows `1..=m`; the top boundary is carried as a scalar, while the padded tail is
processed as a complete vector. Safe `chunks_exact_mut` traversal proves each
score/profile block in bounds without `unsafe`. Prefix-scan constants and
second-symbol profile offsets are computed once before the recurrence rather than
inside every vector block or column.

The frozen score-only `Auto` gate requires at least 4,096 interior cells. The i16
kernel additionally requires `first_len >= 10`; the i32 kernel permits
`first_len == 8` or `first_len >= 11`. Strict SIMD remains available below those
shape gates for diagnostics. These thresholds affect execution choice only, never
the objective or returned score.

## Public API contract

The three existing affine aligners expose:

```text
score(first, second) -> Result<Score, AlignmentError>
score_with_workspace(first, second, workspace) -> Result<Score, AlignmentError>
```

`Gotoh`, `AffineEndGapAligner`, and `AffineLocalAligner` retain their current
configuration types. `AffineWorkspace` retains the existing
`AffineKernelPolicy::{Auto, Scalar, Simd}` selection and reports the completed
backend through `last_backend()`.

`AffineWorkspace::reserve_score` pre-reserves scalar score scratch, while
`reserve_simd_score::<A>` pre-reserves alphabet-specific SIMD score scratch.
`AffineScoreMemoryEstimate` and `simd_score_memory_requirement::<A>` expose the
active element-buffer model before allocation.

Score-only memory estimates and reservation methods must describe only active score
buffers, checked-state flags, SIMD profiles, and cached symbol indices. They must
exclude traceback cells and operation vectors. Reusing a workspace may retain
inactive capacity from earlier complete alignments; such retained capacity is not an
active score-only requirement.

## Data structures

- Scalar: one rolling `H` row, one rolling `D` row, and one byte per column only
  for checked-path gap-state representability.
- SIMD: padded boundary-free `H` and `I` columns, an alphabet-indexed substitution
  profile, second-sequence profile offsets, precomputed scan constants, and a
  signed-64-bit terminal-score row.
- Result: one signed-64-bit `Score`.

No traceback matrix, operation vector, endpoint coordinates, or aligned sequence is
part of this profile.

## Complexity

For lengths `m` and `n`:

- scalar time: `O(mn)`;
- scalar active auxiliary space: `O(n)`;
- SIMD time: `O(mn / L + n log L)` for lane count `L`, excluding profile setup;
- SIMD active auxiliary space: `O(m * alphabet_width + m + n)`;
- result space: `O(1)`.

The SIMD profile is intentionally retained across calls in a caller-owned workspace,
but it is rebuilt when the input or scorer changes.

## Parameters and determinism

All substitution scores, affine gap scores, free-end flags, typed alphabet rules,
and literal ambiguity semantics are unchanged. Score-only output has no traceback
tie to expose. The optimum integer score must nevertheless equal the deterministic
complete alignment's score for every accepted input.

The API is deterministic for a fixed input, scorer, gap model, policy, architecture,
and crate version. Scalar and SIMD backends must return the same integer score.

## Published ambiguities and explicit resolutions

- Gotoh's paper is expressed as a distance minimization. This library uses exact
  score maximization and the already documented algebraic conversion.
- Gap-open conventions differ across tools. The profile retains
  `open + (length - 1) * extend`; adapters must map parameters explicitly.
- Papers do not define software overflow behavior. This implementation uses a
  conservative fast-path proof and checked signed-64-bit fallback.
- Publications do not define automatic SIMD crossover thresholds. Thresholds are
  implementation policy, frozen only after repeated causal benchmarks, and never
  affect the mathematical result.
- A score-only result cannot identify an endpoint or path. Callers requiring either
  must use `align` or `align_with_workspace`; the library will not synthesize them.

## Acceptance criteria

The profile is complete only when all of the following hold:

1. Manual global, local, semiglobal, overlap, and empty-boundary scores pass.
2. Score-only equals complete alignment for DNA, RNA, protein, cached matrices, and
   all 16 free-end configurations.
3. Exhaustive short-input and property-based scalar/SIMD comparisons pass.
4. Bounded and checked arithmetic produce the same scores or typed overflow errors
   as complete alignment.
5. A fresh score-only workspace does not allocate or resize traceback storage.
6. Active memory estimates match the buffers used by scalar and SIMD execution.
7. Prepared black-box preflight and postflight agree with pinned Biopython and
   Parasail cases at matched score-only scope.
8. Criterion reports scalar, strict SIMD, and automatic score-only timings separately
   from complete traceback timings.
9. Length, identity, alphabet, matrix, gap, and boundary strata retain raw paired
   samples and confidence intervals in the isolated AVC laboratory.
10. Automatic SIMD gates select only repeatedly measured winning shapes.
11. A sealed, uninspected holdout is opened only after all preceding gates pass and
    is reported as a separate release-candidate subgate; quarantined records cannot
    participate.

## Acceptance status

Criteria 1--8 and 10 pass. The dedicated release test target has 11 passing tests,
including 8,100 exhaustive scalar score/complete-result comparisons, 7,056 strict
SIMD comparisons on non-empty inputs, 192 generated three-objective cases, 96
forced-i32 generated cases, every free-end bit pattern, checked overflow behavior,
and zero traceback capacity in fresh score-only workspaces.

Criterion 7 is backed by nine final prepared external observations: six against
Biopython 1.87 and three against Parasail 1.3.4. Every preflight score, timed result,
and postflight score agreed. Raw paired measurements and content hashes are recorded
in
`../validation/avc-v1-affine-score-only-checkpoint.md`.
Criterion 9 remains partial because the final campaigns cover only six frozen
development cases rather than the full declared performance surface.

Criterion 11 is open. No holdout case was executed, but one source-pair record was
accidentally displayed during inspection and is therefore quarantined. Neither that
pair nor its currently materialized derived case was used for tuning, testing, or
measurement. They must be excluded and replaced before any T3 claim.

## Known limitations

- The only planned vector backend is runtime-dispatched AVX2 on x86/x86-64.
- ARM NEON/SVE, WebAssembly SIMD, GPU, banding, X/Z-drop, PSSMs, two-piece affine
  gaps, and inter-alignment batching are outside this profile.
- Performance superiority is never inferred from correctness or from a small set of
  development cases.
- On the three frozen global development cases, the optimized Rust kernel remains
  behind Parasail 1.3.4 `scan_32`: reference/Rust ratios are 0.9041 for RNA, 0.9441
  for DNA, and 0.7759 for protein. Closing that remaining gap requires a separately
  documented kernel change rather than weakening arithmetic or validation.
- CPU time, peak RSS, pinned affinity, recorded power mode, broad performance
  stratification, and a resealed T3 holdout remain open.
