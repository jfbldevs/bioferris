# Myers Bit-Vector Edit Distance (`MYERS-bit-vector-v1`)

## Reproduced publication and scope

This profile reconstructs the bit-vector dynamic-programming recurrence in Gene
Myers's 1999 Journal of the ACM paper. It implements two boundary conditions over
the same published recurrence:

1. **Approximate suffix search**, the problem stated in the paper: for every text
   column, report whether some text suffix ending at that column is within a given
   unit edit distance of the complete pattern.
2. **Global edit distance**, obtained by supplying the global top-boundary
   horizontal difference to the paper's general block-advance operation.

The public profile is deliberately narrower than a general biological scoring
aligner. Every substitution, insertion, and deletion costs one; equal symbols cost
zero. DNA, RNA, and protein symbols use literal equality. In particular, IUB
ambiguity symbols and `X` are ordinary symbols rather than wildcards.

The scalar implementation uses `u64` blocks and supports patterns longer than one
machine word. Long global calls add an independently implemented exact adaptive
band; they do not use SIMD, FFI, subprocesses, source code from another
implementation, or a fallback dynamic-programming algorithm.

## Scientific references

Primary source:

1. Gene Myers. "A Fast Bit-Vector Algorithm for Approximate String Matching Based
   on Dynamic Programming." *Journal of the ACM* 46(3), 395–415, 1999.
   DOI: [10.1145/316542.316550](https://doi.org/10.1145/316542.316550).

   Reproduced material: the unit-cost edit-distance recurrence; vertical and
   horizontal difference representation; `Pv`, `Mv`, `Ph`, and `Mh` bit-vectors;
   the `Eq`, `Xv`, and `Xh` equations; the single-word algorithm; the general
   `Advance_Block` interface with horizontal input and output; and the
   `ceil(m / w)` multiword complexity described in Sections 3 and 4.

Independent black-box reference:

2. Martin Šošić and Mile Šikić. "Edlib: a C/C++ Library for Fast, Exact Sequence
   Alignment Using Edit Distance." *Bioinformatics* 33(9), 1394–1395, 2017.
   DOI: [10.1093/bioinformatics/btw753](https://doi.org/10.1093/bioinformatics/btw753).

   Edlib is based on Myers's bit-vector algorithm, documents an exact Ukkonen-style
   band extension, and exposes global edit distance and paths. The paper and its
   supplement inform the global band conditions; implementation source is not
   consulted. A pinned Python wheel is called only from the isolated sibling
   validation laboratory. Edlib is not a build, test, runtime, or source dependency
   of this Rust workspace.

Needleman–Wunsch with unit costs is also used as a mathematically independent test
oracle for small generated cases. That comparison establishes recurrence
equivalence; it does not imply that Needleman–Wunsch is executed by this profile.

## Definitions and boundary conditions

Let pattern `P = p_1 ... p_m`, text `T = t_1 ... t_n`, and `C[i,j]` be the minimum
number of edits required at dynamic-programming cell `(i,j)`. The recurrence is

```text
C[i,j] = min(
    C[i-1,j-1] + [p_i != t_j],
    C[i-1,j]   + 1,
    C[i,j-1]   + 1
)
```

where `[condition]` is one when true and zero otherwise. Both supported modes use

```text
C[i,0] = i.
```

They differ only in the top row:

```text
approximate suffix search: C[0,j] = 0
global edit distance:      C[0,j] = j
```

Consequently, the horizontal difference entering the top of each text column is
zero for approximate search and `+1` for global distance. The original search
boundary makes any text prefix free, so `C[m,j]` is the best distance between the
complete pattern and a suffix ending at text position `j`. The global boundary
makes `C[m,n]` the Levenshtein distance between both complete sequences.

Public search endpoints are Rust-style exclusive offsets `1..=n`. The artificial
column zero is not reported. For an empty pattern, every scanned text endpoint has
distance zero; an empty text therefore produces no search hits.

## Bit-vector equations

Define signed differences

```text
Dv[i,j] = C[i,j] - C[i-1,j]
Dh[i,j] = C[i,j] - C[i,j-1].
```

The unit-cost recurrence guarantees that every difference is in `{-1, 0, +1}`.
For the current text column, `Pv` and `Mv` mark rows whose vertical differences are
`+1` and `-1`; `Ph` and `Mh` do the same for horizontal differences. Bit zero is
the first pattern row. `Peq[s]` marks pattern rows whose literal symbol equals `s`.

For one block, with the current symbol mask `Eq`, Myers gives

```text
Xv = Eq | Mv
Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq
Ph = Mv | ~(Xh | Pv)
Mh = Pv & Xh
```

Unsigned wrapping addition is intentional. Complements and updates are masked to
the active block width, so unused high bits in the last block cannot become state.
The horizontal difference leaving the last active row is `+1` when its `Ph` bit is
set, `-1` when its `Mh` bit is set, and zero otherwise.

The unshifted horizontal vectors are then shifted one row. The incoming horizontal
difference is injected into bit zero: `Ph[0] = 1` for `+1`, `Mh[0] = 1` for `-1`,
and both are zero for `0`. As required by the paper's block operation, an incoming
`-1` also sets `Eq[0]` before evaluating the equations. The next vertical state is

```text
Pv = shifted_Mh | ~(Xv | shifted_Ph)
Mv = shifted_Ph & Xv.
```

Initial state is `Pv = 1` at every active pattern row, `Mv = 0`, and bottom score
`m`, matching column zero.

## Multiword block decomposition

The pattern is partitioned from top to bottom into `ceil(m / 64)` blocks. The full
frontier advances every block for every text symbol in top-to-bottom order. The top
boundary supplies the first horizontal input; each block's horizontal output
becomes the next block's input. The final block output updates the bottom score for
the text column.

```text
advance_column(symbol, top_input):
    horizontal = top_input
    for block in pattern blocks from first to last:
        Eq = Peq[symbol, block]
        (Pv[block], Mv[block], horizontal, Ph, Mh) =
            advance_block(Pv[block], Mv[block], Eq, horizontal, width(block))
    bottom_score += horizontal
```

Pattern preprocessing requires `O(|Σ| ceil(m/64))` words. Scoring a text of length
`n` requires `O(n ceil(m/64))` time and `O(ceil(m/64))` mutable working words, in
the worst case, in agreement with the arbitrary-length construction stated by
Myers. The preprocessed pattern can be reused for multiple texts.

After the scalar implementation passed all correctness tests, a one-block
specialization placed `Pv` and `Mv` in scalar locals for patterns of at most 64
symbols. The equations, masks, boundary input, score updates, and results are
identical; only the otherwise one-element frontier allocation and loop are removed.

## Exact global band extension

Global patterns spanning at least eight words use an exact block band. At each text
column, lower-bound inequalities discard leading and trailing blocks that cannot
belong to any global alignment within `max_distance`. A block immediately below the
active range is initialized from the preceding bottom score only when it can still
enter the band. No surviving cell is approximated or pruned heuristically.

`global_distance_bounded(text, max_distance)` exposes this execution directly and
returns `None` exactly when the true distance exceeds the bound. Automatic
`global_distance` begins with

```text
max_distance = max(32, max(m,n) / 18, abs(m-n))
```

capped at `max(m,n)`. After a failed attempt, twice the current limit is rounded
upward to a power of two and capped again. This keeps the first low-divergence band
narrow while avoiding unnecessary intermediate retries. Automatic execution
remains exact for every divergence. Its best case visits a narrow block interval;
its worst case approaches the full frontier and may repeat work while widening.

For ordinary allocatable dimensions, the bounded kernel selects overflow-proven
arithmetic once before scanning columns. A separate cold checked path preserves
defined behavior for theoretical extreme dimensions. Complete 64-bit blocks are
also specialized once so their no-op masks and partial-word branches are absent
from the hot loop.

Patterns below eight words retain the full-frontier path because the setup and
lower-bound checks cost more than they save at that scale. The public bounded API
still works at every length.

## Global traceback extension

The 1999 paper reports search endpoints and scores; it does not specify an edit
script or traceback algorithm. Returning a `PairwiseAlignment` is therefore a
documented library extension, not a claim about unpublished Myers behavior.

The optimized extension stores two decision masks per retained block: one marks
valid diagonal predecessors and one marks valid deletions. Insertion is the
remaining valid predecessor. Equality positions always take the valid zero-cost
diagonal. For mismatches, the diagonal mask is derived from the published vertical
and horizontal difference updates while the column is advanced. These decisions
are sufficient to recover one exact path without retaining four full difference
planes or recomputing a conventional matrix.

For long patterns, traceback first obtains the exact distance and reruns the bounded
kernel with that known limit. Only active block decision masks are appended, with a
compact `(first_block, history_end)` record for each column. Short patterns retain
the same two-mask representation for the full frontier. The deterministic choice
at every cell is:

1. diagonal (`Match` for equal symbols, otherwise `Substitution`);
2. up (`Deletion`);
3. left (`Insertion`).

The backward walk takes `O(m+n)` time after scoring. Full-frontier history has a
conservative `2 n ceil(m/64)`-word upper bound; the long banded implementation
retains only active records plus `O(n)` column metadata. History is intentionally
separated from distance-only and search APIs, which retain linear frontier space in
the number of blocks.

The returned alignment score is the negative edit distance. This maps the fixed
cost model to the workspace convention `match = 0`, `mismatch = -1`, and
`gap = -1`; the dedicated distance APIs continue to return a nonnegative `usize`.

## Public components

- `MyersPattern<A>` owns a validated pattern and its `Peq` table and may be reused.
- `global_distance(text)` returns exact global unit edit distance without traceback
  history.
- `global_distance_bounded(text, max_distance)` returns the exact distance when it
  is within the inclusive limit and `None` otherwise.
- `global_distance_with_workspace` and
  `global_distance_bounded_with_workspace` reuse frontier and band-score buffers.
- `find_all(text, max_distance)` returns every exclusive text endpoint whose best
  pattern-to-suffix distance is no greater than the threshold.
- `align_global(text)` returns one deterministic optimal global operation path and
  the negated distance score.
- `align_global_with_workspace` reuses frontier and compact history memory while
  returning an independently owned alignment.
- `MyersWorkspace` exposes retained-byte reporting and explicit buffer release.
- `MyersMatch` records an exclusive endpoint and its exact distance.
- `MyersMemoryEstimate` exposes modeled persistent, distance-working, traceback,
  and output buffer sizes for reproducible memory benchmarks.
- `verify_unit_edit_alignment` independently replays an alignment, checks literal
  operation semantics and global consumption, and compares the recomputed edit
  count with the reported negative score.

Construction and execution return typed errors for dimension arithmetic,
allocation failure, score conversion, invalid history invariants, and invalid
alignment structure. No unchecked indexing or `unsafe` code is used.

## Parameters and determinism

Approximate search has an inclusive maximum distance that filters reported
endpoints but still advances its exact search frontier. Bounded global distance has
an independent inclusive limit that safely restricts the global frontier.
Automatic global distance and traceback choose and widen exact limits internally.
There are no scoring parameters because the scientific profile fixes unit costs.

Pattern order and text order are never exchanged automatically. Peq layout,
block traversal, hit order, and traceback tie-breaking are deterministic. Search
hits are returned in increasing exclusive endpoint order.

## Published omissions and explicit resolutions

| Issue | Resolution in `MYERS-bit-vector-v1` |
|---|---|
| The main paper problem has a free top row, while users also expect global edit distance. | Expose both. Search uses horizontal input `0`; global distance uses `+1` through the published general block interface. |
| The paper does not publish traceback. | Derive exact predecessor-decision masks from the published difference updates and document this as an extension. |
| Full multiword traversal wastes work for similar long sequences. | Apply the exact global band inequalities described by the Edlib paper/supplement, expose a bounded API, and widen automatically without changing the optimum. |
| Word size is machine-dependent in the analysis. | Fix blocks to `u64` for reproducible layout and behavior. |
| A partial final word has unspecified inactive bits. | Mask every complement and state update to the active width. |
| Biological ambiguity semantics are outside the string algorithm. | Use literal validated-symbol equality; no wildcard interpretation. |
| Search coordinates can be reported with different conventions. | Return exclusive Rust offsets and omit artificial column zero. |
| Multiple optimal global paths may exist. | Prefer diagonal, then deletion, then insertion during derived traceback. |
| Empty-pattern endpoint reporting is not an API-level concern in the paper. | Report distance-zero hits for scanned endpoints only; an empty text has no endpoint. |

## Independently testable implementation stages

1. Build and inspect literal `Peq` masks, including partial and multiple blocks.
2. Advance one block for all three horizontal inputs and verify its output/state.
3. Compare single- and multi-block global distances with hand-computed cases and
   an independent row dynamic program.
4. Verify approximate-search endpoint scores against the search boundary matrix.
5. Compare bounded and automatically widened distance with the complete frontier at
   word boundaries and on diverse long generated inputs.
6. Derive traceback from compact decision masks and replay it independently.
7. Compare static results with the pinned Edlib black box.
8. Benchmark preprocessing, reused-pattern/workspace distance, known and automatic
   limits, divergence, search, and traceback separately.

## Acceptance criteria

This phase may be marked complete only when all of the following hold:

- hand-calculated empty, identical, substitution, insertion, deletion, terminal
  gap, internal gap, completely different, and repetitive cases pass;
- DNA, RNA, and protein cases demonstrate literal alphabet-generic behavior;
- exhaustive short inputs agree with an independent scalar edit-distance matrix;
- property tests cover symmetry of global distance, identity, length lower bound,
  upper bound, triangle inequality, valid replay, and deterministic repetition;
- pattern lengths around block boundaries (`63`, `64`, `65`, `127`, `128`, and
  `129`) agree with the independent oracle;
- search endpoints and distances agree with an independently initialized
  approximate-search matrix;
- returned global operations replay to the reported score independently;
- insufficient bounded limits return `None`, sufficient limits return the exact
  full-frontier score, and long banded paths equal the independent oracle;
- committed, immutable fixtures generated by the isolated Edlib laboratory agree
  in distance and admit the Rust path as an optimal edit script;
- deterministic Criterion benchmarks record time and modeled memory for synthetic
  and curated scientific inputs;
- a mirrored Rust/Edlib benchmark uses equivalent inputs, preprocessing policy,
  requested result kind, warm-up, repetitions, and result consumption;
- formatting, locked tests, examples, Clippy, documentation, license policy, and
  dependency audit gates pass with no placeholders or unfinished routes.

## Known profile limitations

- Costs are fixed to unit Levenshtein distance; substitution matrices, arbitrary
  linear scores, and affine gaps are outside this algorithm.
- Ambiguity symbols are literal rather than degenerate sets.
- Traceback history remains `O(n ceil(m/64))` in the worst case even though active
  low-distance bands are usually smaller and distance/search use linear frontier
  memory.
- This scalar implementation has no SIMD, threading, or instruction-set-specific
  dispatch.
- At high divergence, automatic widening can revisit columns before approaching
  full-frontier work; callers with a trustworthy threshold should use the bounded
  API directly.
- Approximate search returns endpoints and distances, not start coordinates or one
  alignment for every hit; the original bit-vector frontier does not retain that
  information.

## Completion evidence (updated 2026-08-10)

The versioned profile and its public API satisfy this phase's acceptance criteria:

- the complete `align-fast` unit, integration, property, and doc-test set passes,
  including bounded/full-frontier differential tests over word boundaries and 160
  diverse long generated pairs;
- exhaustive binary DNA pairs through length four and generated DNA through length
  139 agree with independent row dynamic programming, while explicit cases cross
  63/64/65 and 127/128/129 word boundaries;
- all returned paths pass the independent unit-edit verifier and repeated calls are
  deterministic, including 512-, 513-, and 1,024-symbol banded tracebacks;
- 20 pinned Edlib 1.3.9.post1 cases agree on every global distance and all 812
  per-endpoint search scores; the fixture covers DNA, RNA, protein, empty inputs,
  ambiguity symbols, word boundaries, and the curated scientific protein pair;
- `cargo test --workspace --all-targets --locked` passes the complete workspace,
  examples, and benchmark smoke cases;
- release Criterion time and exact buffer-memory measurements are recorded in
  `../benchmarks/myers.md`;
- the counterbalanced one-shot Rust/Edlib comparison, with Rust distance and
  traceback wins at all five tested lengths in both rounds, is recorded in
  `../benchmarks/myers-vs-edlib.md`;
- formatting, Clippy with warnings denied, rustdoc with warnings denied, and
  `cargo-deny` advisories, bans, licenses, and sources all pass;
- production code contains no `unsafe`, external engine, FFI, subprocess, fallback
  algorithm, placeholder, mock, or unfinished route.

The remaining items in **Known profile limitations** are documented boundaries, not
unfinished code paths. Affine costs and SIMD belong to separate scientific
profiles.
