# Exact Affine SIMD Backend

Status: complete and validated on 2026-08-07 UTC for the exact profile and
platform scope defined here. All acceptance evidence is pinned in the final
checkpoint linked below.

Profile identifier: `AFFINE-SIMD-PREFIX-SCAN-v1`.

## Scope

This is an execution backend for the existing exact three-state affine recurrence.
It accelerates global Gotoh, affine local alignment, and all 16 combinations of the
four independent free-end flags. It does not define a new biological objective and
does not replace any named algorithm. Protein, DNA, and RNA use the same typed
sequence and substitution-scoring APIs as the scalar implementation.

This profile returns complete alignments with traceback. The separately documented
[`AFFINE-SCORE-ONLY-v1`](affine-score-only.md) specialization reuses the exact
prefix recurrence without trace codes and has its own crossover and benchmark
evidence.

## Scientific references

1. Osamu Gotoh. "An improved algorithm for matching biological sequences."
   *Journal of Molecular Biology* 162(3), 1982, 705--708.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
2. Jeff Daily, Ananth Kalyanaraman, Sriram Krishnamoorthy, and Bin Ren. "On the
   Impact of Widening Vector Registers on Sequence Alignment." *45th International
   Conference on Parallel Processing*, 2016, 506--515.
   DOI: [10.1109/ICPP.2016.65](https://doi.org/10.1109/ICPP.2016.65).
3. Jeff Daily. "Parasail: SIMD C library for global, semi-global, and local pairwise
   sequence alignments." *BMC Bioinformatics* 17, 81, 2016.
   DOI: [10.1186/s12859-016-0930-z](https://doi.org/10.1186/s12859-016-0930-z).
4. Thomas F. Smith and Michael S. Waterman. "Identification of common molecular
   subsequences." *Journal of Molecular Biology* 147(1), 1981, 195--197.
   DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).

Only papers, equations, published descriptions, and independently observable
black-box behavior were used. No source code from Parasail or any other alignment
implementation was inspected, copied, translated, or adapted.

## Exact objective

Let `H[i,j]` be the best score for two prefixes, `D[i,j]` a path ending in a
deletion from the second sequence, and `I[i,j]` a path ending in an insertion into
the first sequence. Let `o` be this library's gap-open score and `e` its extension
score. A gap of length one scores `o`; each later symbol scores `e`.

```text
D[i,j] = max(H[i-1,j] + o, D[i-1,j] + e)
I[i,j] = max(H[i,j-1] + o, I[i,j-1] + e)
S[i,j] = H[i-1,j-1] + substitution(first[i], second[j])
H[i,j] = max(S[i,j], D[i,j], I[i,j])
```

Local mode adds zero to the last maximum. Boundary initialization and endpoint
selection are exactly those in `gotoh.md` and `affine-alignment-modes.md`.

Daily et al. express opening as a term followed by extension. The exact mapping is
`Gext = e` and `Gopen_paper = o - e`; substitution yields this library's convention
without charging an extra extension to a one-symbol gap.

## Max-plus prefix transformation

One second-sequence column is evaluated at a time. For fixed `j`, compute the
horizontal insertion and a temporary value excluding the vertical deletion:

```text
I[i] = max(H_previous[i] + o, I_previous[i] + e)
Q[i] = max(H_previous[i-1] + substitution(i,j), I[i])
```

Local mode also includes zero in `Q[i]`. Because the validated affine contract has
`o <= e`, opening a deletion from an `H` value already ending in deletion cannot
beat extending it. The remaining recurrence is therefore exactly:

```text
D[i] = max(Q[i-1] + o, D[i-1] + e)
H[i] = max(Q[i], D[i])
```

For a vector of `p` consecutive rows, initialize:

```text
seed[0] = max(H_carry + o, D_carry + e)
seed[l] = Q[l-1] + o                         for l > 0
```

At the first scored row only the open candidate is valid. An inclusive max-plus
scan applies offsets `1, 2, 4, ...`:

```text
D_next[l] = max(D_current[l], D_current[l-offset] + offset * e)
```

Only lanes at or above the offset participate. After `ceil(log2(p))` stages, every
active lane is identical to the scalar deletion recurrence. The balanced
Hillis--Steele scan and consecutive-row layout are explicit reconstruction choices;
the 2016 paper describes parallel scan and a concrete striped layout.

## Pseudocode

```text
prove all partial scores and scan arithmetic fit i16 or i32
detect AVX2; otherwise return strict error or use scalar fallback
pad score/profile backing to a complete vector
build a substitution profile indexed by second symbol and first position
initialize the column-zero boundary

for column j from 1 through N:
    initialize the row-zero boundary
    diagonal_carry = old H[0]
    deletion_carry is unavailable at row zero

    for each padded SIMD block of first-sequence rows:
        old_H = load previous-column H
        old_I = load previous-column I
        diagonal = shift_in(diagonal_carry, old_H) + profile[j]
        I = max(old_H + open, old_I + extend)
        Q = max(diagonal, I [, zero in local mode])
        D = max-plus-prefix-scan(seed(Q, H carry, D carry), extend)
        H = max(Q, D)

        reproduce scalar predecessor choices as vector masks
        store H and I for the next column
        copy trace bytes for active lanes only
        observe local endpoints in canonical row-major order
        advance scalar carries from the last active lane

    record H[M,j] for free-end endpoint selection

select the endpoint with the established scalar rules
trace back through the column-major trace matrix
return the ordinary independently owned PairwiseAlignment
```

Padding replaces a separate scalar tail. Inactive lanes can receive arbitrary
intermediate values, but dependencies flow only from lower to higher lane indices;
they cannot affect an active lane, endpoint observation, or traceback byte.

## Data structures and public API

- `AffineKernelPolicy::Auto` selects SIMD only when architecture, numeric proof,
  and calibrated shape gates pass. `Scalar` forces the established implementation;
  `Simd` requires this backend and otherwise returns `SimdUnavailable`.
- `AffineBackend` reports `Scalar`, `SimdI16 { lanes }`, or
  `SimdI32 { lanes }` for the most recent successful workspace call.
- `AffineWorkspace` retains two rolling score columns, a terminal-score row, the
  transposed substitution profile, second-sequence alphabet indices, and a separate
  column-major trace matrix.
- The backing score/profile vectors contain padded `i32` elements. The 16-lane
  signed-16-bit path safely reinterprets that storage with `pulp::bytemuck`; the
  8-lane signed-32-bit path uses it directly.
- `AffineWorkspace::simd_memory_requirement`, `reserve_simd`, and the SIMD capacity
  accessors expose the active padded model and caller-owned allocation lifetime.
- Public `Score` remains `i64`. Terminal values are widened exactly. Returned
  operation vectors are independently owned and never alias workspace storage.

## Runtime dispatch and safety

The workspace keeps `unsafe_code = "forbid"`. It pins `pulp = 0.22.3` with default
features disabled and only `std` and `x86-v3` enabled. On x86/x86-64,
`pulp::x86::V3::try_new` performs runtime AVX2 feature detection and
`Simd::vectorize` enters the specialized implementation. Other architectures and
x86 processors without this feature set use the exact scalar path under `Auto`.

Project code contains no `unsafe` block and no raw-pointer access. It invokes only
Pulp's safe vector and intrinsic-wrapper methods. AVX-512 and non-x86 vector
implementations are intentionally outside this profile.

Pulp is MIT licensed. Its resolved normal/build graph is `bytemuck`, `cfg-if`,
`libm`, `num-complex`, `num-traits`, `autocfg`, `paste`, `pulp-wasm-simd-flag`,
`reborrow`, `raw-cpuid`, `bitflags`, and `version_check`. The locked graph must pass
the workspace's ban, license, and source audit and Rust 1.85.0 MSRV gate. The one
informational `paste` maintenance advisory is explicitly reviewed in
`docs/dependencies.md`; it is not a vulnerability waiver.

## Lane-width proof and fallback

The preferred kernel uses 16 signed-16-bit lanes. A conservative preflight proves
that every substitution value, affine boundary, complete or partial path score, and
prefix-scan offset is representable in `i16`. If only the analogous `i32` proof
passes, the 8-lane kernel is selected. If the `i32` proof cannot be established,
`Auto` uses the ordinary bounded-`i64` or checked scalar recurrence and `Simd`
returns a typed error.

No SIMD score arithmetic saturates or clips. Unknown scorer bounds therefore cannot
enter SIMD even if sampled values appear small. Empty inputs contain no vector
recurrence cells and follow the same scalar/strict policy.

## Automatic dispatch gate

Calibration used forced scalar and forced SIMD complete calls with one reused
workspace and identical inputs. The frozen `Auto` gate is deliberately conservative:

- at least 4,096 interior DP cells are required;
- signed-16-bit SIMD requires `first_len >= 7`;
- signed-32-bit SIMD requires `first_len` in `6..=8` or `first_len >= 11`.

At the 4,096-cell boundary, i16 width 6 was effectively parity/slightly slower and
i32 widths 5 and 9--10 did not show a stable win. Strict SIMD remains available for
diagnostics below the gate, but `Auto` does not select it. These thresholds are an
implementation policy, not a scientific property, and must be recalibrated before a
different SIMD architecture is enabled.

## Deterministic tie-breaking

Vector masks reproduce the scalar decisions, not merely an optimal score:

1. gap extension wins equality with gap opening;
2. deletion wins equality with insertion;
3. diagonal wins equality with either gap state;
4. local zero wins every non-positive candidate and records `Start`;
5. a local endpoint is the earliest positive maximum in row-major order;
6. free-end endpoint rank is score, consumed symbols, then row;
7. operation classification still depends on actual symbol equality in traceback.

## Complexity and memory

For sequence lengths `M` and `N`, alphabet width `A`, and `p` lanes:

- time is `O(MN)` scalar work-equivalent, implemented as approximately
  `O(N * ceil(M/p) * log p)` vector scan stages plus `O(AM)` profile construction;
- rolling score storage is `O(M + p)` because the first dimension is padded;
- the substitution profile is `O(A(M + p))`;
- complete traceback remains `O(MN)` one-byte cells;
- traceback and returned operations are `O(M + N)`.

The public memory estimate includes active padding and the maximum operation vector.
Allocator metadata, allocator rounding, inputs, and excess capacity retained from a
larger previous call are excluded and reported separately through capacity accessors.

## Published omissions and explicit resolutions

| Omission or mismatch | Resolution in this profile |
|---|---|
| The scan paper omits traceback details. | Reconstruct the established scalar predecessor comparisons from final vectors and store the same four trace fields. |
| The concrete published layout is striped and uses two passes. | Use consecutive-row blocks and an associative max-plus prefix scan; preserve the equations and document the layout difference. |
| The papers do not define four independent free-end flags. | Reuse the frozen scalar boundary/endpoint contract and verify all 16 masks. |
| Paper gap notation charges extension after opening. | Map `Gopen_paper = o - e` and `Gext = e`. |
| The papers do not define a safe narrow-lane dispatch proof. | Conservatively bound every partial path and prefix operation for i16, then i32; otherwise use exact scalar arithmetic. |
| A partial vector is implementation-specific. | Pad vector backing, keep dataflow one-directional, and expose only active lanes to trace and endpoint logic. |
| CPU-dependent selection can obscure reproducibility. | Expose forced scalar/strict SIMD policies and report the selected backend and lane count. |

## Acceptance criteria

The profile is complete only when all of the following hold:

1. Complete score, ranges, and operations equal forced-scalar results.
2. Exhaustive short inputs cover global, local, and all 16 free-end masks.
3. Property tests cross 8/16-lane and padded boundaries, change scores, gaps,
   workspace shapes, modes, alphabets, and scorer implementations.
4. Independent global, local, and free-end verifiers replay every tested result.
5. i16/i32 eligibility, unknown bounds, overflow boundaries, empty input,
   unavailable ISA, and strict errors are tested.
6. Forced scalar, automatic, and strict SIMD calls are deterministic and never
   mutate previously returned results.
7. Rust 1.85.0 and all formatting, lint, test, rustdoc, dependency-license,
   advisory, ban, and source gates pass.
8. Criterion measures crossover lengths, vector/tail widths, all affine modes,
   generic and dense scorers, protein/DNA/RNA, and fresh/reused workspaces.
9. Same-binary scalar/SIMD and external Biopython/Parasail campaigns use prepared
   complete calls, 31 alternating AB/BA samples, and at least 0.25 seconds per
   measured batch.
10. Every external case passes preflight, independent replay, postflight, duration,
    digest, and stderr gates.
11. A claimed win requires a 95% bootstrap lower confidence bound above one; a
    losing shape is not selected by `Auto`.
12. No reference source is used and no external runtime enters the product graph.

The final evidence record is
`../validation/avc-v1-affine-simd-checkpoint.md`.

## Limitations

- Only AVX2 on x86/x86-64 is implemented; there is no AVX-512, Arm, WASM, GPU, or
  inter-alignment batching backend.
- Score-only execution is a separate public profile; dual-affine gaps, PSSMs,
  banding, X-drop, and Z-drop fusion remain outside this profile.
- Complete traceback remains `O(MN)` and can become memory-bandwidth limited.
- The automatic gate is calibrated on one host and is conservative; it is not a
  regime-wide performance guarantee.
