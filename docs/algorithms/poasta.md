# POASTA: Exact Accelerated Gap-Affine Partial Order Alignment

## Reproduction target and status

This document was created before the accelerated implementation and froze the
scientific and engineering contract for
`POASTA-2025-global-affine-ASTAR-DFS-superbubble-v1` before production code was
written.

Status: complete and validated on 2026-08-05.

The target is the exact global, gap-affine POASTA algorithm published by van
Dijk et al. (2025). It aligns one complete query to one complete source-to-sink
path through a character-labeled partial-order DAG. It must preserve all three
characteristic accelerations from the publication:

1. sparse A* traversal of the three-state affine alignment graph using the
   minimum-remaining-gap heuristic;
2. depth-first traversal of previously unvisited zero-cost exact matches; and
3. safe pruning informed by a superbubble index and implicitly opened gaps.

The already implemented full-matrix profile
`POASTA-2025-supplement-global-affine-DP-v1` is the independent semantic oracle.
An implementation containing only A*, only the heuristic, or only exact-match
DFS must not be called the frozen POASTA profile.

The first production profile is scalar and single-threaded. SIMD and parallel
graph construction are later optimizations, not substitutes for any published
step.

## Primary references and frozen versions

1. van Dijk LR, Manson AL, Earl AM, Garimella KV, Abeel T. *Fast and exact
   gap-affine partial order alignment with POASTA*. Bioinformatics.
   2025;41(1):btae757. Published 3 January 2025; corrected typeset record
   22 January 2025.
   DOI: [10.1093/bioinformatics/btae757](https://doi.org/10.1093/bioinformatics/btae757).
   The main paper is normative for the A* traversal, exact-match DFS, and
   superbubble pruning.
2. The official supplemental material distributed with the article, dated
   19 November 2024. Supplemental Equation 2 and pages 5-7 are normative for
   the affine alignment graph, heuristic equations, multiple-exit-position
   bounds, and the affine pruning rule.
3. Gartner F, Muller L, Stadler PF. *Superbubbles revisited*. Algorithms for
   Molecular Biology. 2018;13:16.
   DOI: [10.1186/s13015-018-0134-3](https://doi.org/10.1186/s13015-018-0134-3).
   Algorithm 2 is normative for linear-time superbubble detection in a DAG.
4. Gotoh O. *An improved algorithm for matching biological sequences*. Journal
   of Molecular Biology. 1982;162(3):705-708.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
   This is the affine pairwise foundation generalized by the first reference.
5. Lee C, Grasso C, Sharlow MF. *Multiple sequence alignment using partial order
   graphs*. Bioinformatics. 2002;18(3):452-464.
   DOI: [10.1093/bioinformatics/18.3.452](https://doi.org/10.1093/bioinformatics/18.3.452).
   This defines the PO-MSA graph semantics and progressive graph incorporation.

The retrieved supplemental artifacts are pinned by:

```text
btae757_supplementary_data.zip
SHA-256 d7639047ecbec4215d6420d53bdf935283f85223f9de9bd0d81484aeccb80d6d

034a0_Poasta_paper_OUP_Supplemental__Revision.pdf
SHA-256 1ccf7287107cd98fab1d526cbff152ca47436524c7d5e1071f5d6ce948dbefde
```

No source code from POASTA, SPOA, abPOA, or another aligner may be inspected,
translated, adapted, or used in production. Pinned executables may be invoked
only by the isolated black-box validation project.

## Objective and graph model

Let `G=(V,E)` be a validated character-labeled DAG and let
`Q=q_1...q_m` be a typed protein, DNA, or RNA sequence. Add an implicit virtual
source `nu` with an edge to every root and an implicit virtual sink `tau` with an
incoming edge from every sink. For an empty graph, add `nu -> tau`.

POASTA finds a minimum-cost path from `<nu,0,M>` to `<tau,m,M>` in the
alignment product graph. Equivalently, the public API returns the negated cost
as a maximum score and one deterministic alignment of the complete query to a
complete source-to-sink graph path.

Nodes outside the selected graph path are alternatives and cost nothing.
Leading, internal, and trailing gaps are penalized. The graph is immutable for a
single sequence-to-graph call. A separate iterative adapter applies the existing
Lee graph-incorporation rules after each verified alignment.

## Published cost model

POASTA minimizes non-negative cost:

```text
C = N_m * Delta_m
  + N_x * Delta_x
  + N_o * Delta_o
  + N_g * Delta_e
```

where `N_o` is the number of gap runs and `N_g` is the number of symbols in all
gaps. The frozen profile requires:

```text
Delta_m = 0
Delta_x >= 0
Delta_o >= 0
Delta_e >= 0
Delta_o > Delta_e
gap_cost(k) = Delta_o + k * Delta_e, k > 0
```

The paper benchmark uses `Delta_x=4`, `Delta_o=6`, and `Delta_e=2`.
`PoastaAffineCosts` is the sole scoring type accepted by this accelerated
profile. Arbitrary substitution matrices are supported by the full-matrix
oracle but are not silently attributed to POASTA because exact-match DFS and
the heuristic proof depend on zero-cost literal matches and non-negative other
costs.

The workspace score is exactly `-C`. The existing `AffineGap` conversion is:

```text
match_score     = 0
mismatch_score  = -Delta_x
open_score      = -(Delta_o + Delta_e)
extension_score = -Delta_e
```

## Published affine alignment graph

A search state is `(u,i,S)` where `u` is the last consumed graph node or `nu`,
`i` is the consumed query-prefix length in `0..=m`, and `S` is `M`, `I`, or
`D`. The real-node transitions reconstructed from Supplemental Equation 2 are:

```text
M(u,i) -> I(u,i+1)                         Delta_o + Delta_e
I(u,i) -> I(u,i+1)                         Delta_e

M(u,i) -> D(v,i), for each graph edge u->v Delta_o + Delta_e
D(u,i) -> D(v,i), for each graph edge u->v Delta_e

M(u,i) -> M(v,i+1), for each edge u->v     Delta(label(v), q[i])
I(u,i) -> M(u,i)                           0
D(u,i) -> M(u,i)                           0
```

Edges from `nu` to roots use the same transitions. A zero-cost edge connects
`M(sink,m)` to `M(tau,m)`. Insertions cannot originate at `tau`, and deletions
cannot consume `tau`. The goal for a non-empty graph is therefore any
`M(sink,m)`; for an empty graph it is `M(nu,m)`.

Every edge has non-negative cost, so Dijkstra and A* are exact. Gap closing is a
real zero-cost transition. Adjacent insertion and deletion runs are therefore
legal and must not be removed from the search space.

## Minimum-remaining-gap A* heuristic

For cursor `u`, let `d_min(u,tau)` and `d_max(u,tau)` be the minimum and maximum
edge counts to `tau`. Define the remaining real graph symbols and query symbols:

```text
g_min = d_min(u,tau) - 1
g_max = d_max(u,tau) - 1
q_rem = m - i
```

The minimum unavoidable number of gap symbols is:

```text
N_gap = q_rem - g_max  if g_max < q_rem
        g_min - q_rem  if g_min > q_rem
        0              otherwise.
```

Supplemental Definitions 2-4 give:

```text
h(M) = 0                              if N_gap = 0
       Delta_o + N_gap * Delta_e      otherwise

h(I) = 0                              if N_gap = 0
       N_gap * Delta_e                if g_max < q_rem
       Delta_o + N_gap * Delta_e      otherwise

h(D) = 0                              if N_gap = 0
       N_gap * Delta_e                if g_min > q_rem
       Delta_o + N_gap * Delta_e      otherwise.
```

The state-specific cases avoid charging a second opening cost when all
unavoidable remaining gaps can extend the currently open gap. Mismatches are
assigned a lower bound of zero. Thus `h` never exceeds the true remaining cost.
Minimum and maximum distances are computed once in reverse topological order in
`O(|V|+|E|)` time and `O(|V|)` memory, including implicit sentinels.

## Sparse A* baseline

The sparse baseline retains only discovered alignment states. A packed key combines
the query prefix and affine state, and each compact record stores the best exact
`u64` path cost, canonical predecessor key/tie code, and expansion cost. A bucket
frontier handles priorities through 65,536; larger priorities use an overflow heap.
Both paths are ordered by increasing `f=g+h`, then deterministic state keys. Stale
entries are ignored.

A lower `g` reopens a previously expanded state. This makes exactness independent
of heuristic consistency; admissibility alone is sufficient. Equal-cost arrivals
may improve traceback metadata under the frozen tie rules without requiring a
second expansion because descendants refer to predecessor state keys rather than
copying traceback paths.

After finding a goal of cost `C*`, search continues while the smallest frontier
`f <= C*`. This settles equal-cost endpoint candidates and canonical predecessor
ties. It terminates when the frontier is empty or its minimum `f > C*`.

## Exact-match DFS

When an `M(u,i)` state is selected from the A* frontier, POASTA starts a DFS:

1. inspect every deterministic graph successor `v`;
2. if `i < m`, `label(v)==q[i]`, and `M(v,i+1)` has not been visited at an
   equal or lower cost, relax the zero-cost match, push it onto the DFS stack,
   and suppress only the deletion alternative adjacent to that successor;
3. enqueue a mismatch and its adjacent deletion normally;
4. for an exact target already reached at an equal or lower cost, enqueue its
   adjacent deletion normally;
5. enqueue insertion when at least one real successor was a mismatch or an
   already-reached exact target, when there was no newly visitable exact target,
   or when the graph path has ended while query symbols remain;
6. process all branches in postorder-style DFS until the stack is empty, while
   retaining the ordinary queue for costly alternatives.

The publication's correctness argument is that a newly unvisited zero-cost match
cannot be improved by a path containing only non-negative-cost alternatives.
Once the exact target is already visited, that argument no longer permits
discarding alternatives, so its deletion and the insertion alternative are
queued normally. Suppression is per successor rather than per source state: a
different mismatching branch must retain both gap directions even when another
branch matches exactly. This resolution is required by the finite alignment
graph; suppressing every gap from a source having any exact successor changes
the optimum on a three-symbol counterexample. `I` and `D` states and their
zero-cost closing edges remain ordinary A* work.

The implementation keeps a feature-isolation mode that disables DFS. It is an
exact A* baseline for differential tests, not a separately advertised algorithm.

## Superbubble index

The graph is augmented by `nu` and `tau` for preprocessing. Superbubbles are
identified with Gartner et al. Algorithm 2, not with a guessed degree-only
pattern:

1. construct a deterministic reverse DFS postorder `pi`;
2. compute `OutParent(v)=min(pi(p))` over predecessors and
   `OutChild(v)=max(pi(c))` over successors, using negative/positive infinity at
   the augmented endpoints;
3. scan positions in reverse with the published exit-candidate stack and
   `outParentMap` interval minima;
4. report `(s,t)` exactly when conditions F1 and F2 and the nearest-entrance
   minimality condition are satisfied.

For each node contained in each reported bubble, the index stores its bubble
exit plus the shortest and longest path length from that node to the exit.
Membership includes the entrance and interior and excludes the exit; nested
bubbles are retained. Distances count graph edges, matching the number of
zero-cost diagonal transitions required to reach the exit.

The detector itself is `O(|V|+|E|)`. Materializing every node-to-containing-
bubble entry is output-sensitive; its explicit bound is `O(|V|+|E|+K)`, where
`K` is the number of stored membership entries. Distances are computed within
each bubble in reverse topological order. This may revisit edges for nested
bubbles and has worst-case `O(sum_b(|V_b|+|E_b|))` preprocessing time. That cost
must be reported rather than mislabeled as globally linear.

## Superbubble-informed pruning

For every bubble exit `t`, retain an ordered map
`J_t^M = {j -> C(t,j,M)}` of settled match-state query positions and their
lowest costs. When a candidate `(v,i,S)` of cost `C` belongs to bubble `(s,t)`,
its best-case zero-cost match paths reach exit positions:

```text
j_min = i + d_min(v,t)
j_max = i + d_max(v,t).
```

If either bound exceeds `m`, it is clipped to the valid query domain; if the
resulting interval is empty, that bubble cannot prune the candidate.

For a reached exit state `(t,j,M)` with cost `C_j`, an implicit gap supplies an
upper bound at exit coordinate `x`:

```text
U_j(x) = C_j                                      if x = j
         C_j + Delta_o + (x-j)*Delta_e            if x > j
         C_j + Delta_o + (j-x)*Delta_e            if x < j and a downstream
                                                    path of at least j-x nodes
                                                    exists after t
         infinity                                  otherwise.

U(x) = min_j U_j(x).
```

The `x>j` case is an implicitly opened insertion at `t`. The `x<j` case is an
implicitly opened deletion after `t`, traced backward through the same number of
best-case match edges as described in the supplement. Feasibility is checked
against the maximum remaining path length after `t`.

Because each implicit-gap cost is piecewise linear, the supplement proves that
it is sufficient to compare `C` at:

```text
{j_min, j_max}
union
{j-1, j, j+1 | j is a reached position in [j_min,j_max]}
```

after retaining only valid positions in the interval. The candidate is pruned
only if a finite `U(x) <= C` for every checked position. Missing evidence always
keeps the candidate.

Only reached `M` states are recorded, matching the published B-tree reduction.
For an `I` or `D` candidate, the supplement adds a mandatory guard: do not prune
if any reached match state satisfies

```text
C < C(t,j,M) + Delta_o.
```

This preserves the advantage of extending the candidate's already-open gap.
All arithmetic is checked and exact. Pruning may remove equal-cost alternatives
but may never change the optimum score.

The phrase "popped from the A* queue" in the supplement conflicts with the main
paper's illustrated interaction between DFS matches and reached bubble exits.
This profile records an exit when an `M` state is settled for expansion from
either frontier. That is the only explicit resolution: it broadens available
upper-bound evidence but does not change the proof because every recorded state
has a realized path and exact cost. Differential tests must separately validate
queue-only and both-frontier recording before the latter remains enabled.

## Reconstructed pseudocode

```text
POASTA_ALIGN(graph, query, costs, options):
    VALIDATE(graph, costs)
    metadata <- MIN_MAX_REMAINING_PATH_LENGTHS(graph + nu + tau)
    bubbles  <- SUPERBUBBLE_INDEX(graph + nu + tau)

    records[<nu,0,M>] <- {g=0, predecessor=Start}
    PUSH_ASTAR(<nu,0,M>, f=h(<nu,0,M>))
    best_goal <- none

    while astar_queue is not empty:
        item <- POP_MIN_VALID_ITEM()
        if best_goal exists and item.f > best_goal.cost:
            break
        if item is stale or was expanded at an equal/lower g:
            continue

        if exact_match_dfs and item.state is M:
            PUSH_DFS(item.state)
            while dfs_stack is not empty:
                state <- DFS_TOP()
                if not entered:
                    if PRUNABLE(state, bubble_index, reached_exits):
                        mark pruned; POP_DFS(); continue
                    mark expanded at g
                    RECORD_EXIT_IF_MATCH_STATE(state)
                    for successor in stable order:
                        if successor is a newly visitable exact match:
                            RELAX_TO_DFS(successor)
                        else:
                            QUEUE_DIAGONAL_IF_APPLICABLE(successor)
                            QUEUE_DELETION(successor)
                            require_insertion <- true
                    if no new exact match or require_insertion:
                        QUEUE_INSERTION(state)
                else if unprocessed DFS child exists:
                    process that child
                else:
                    POP_DFS()
        else:
            if PRUNABLE(item.state, bubble_index, reached_exits):
                mark pruned; continue
            mark expanded at g
            RECORD_EXIT_IF_MATCH_STATE(item.state)
            RELAX_ALL_PUBLISHED_OUTGOING_EDGES(item.state)

        for every settled M(sink,m):
            update best_goal by cost, then canonical sink ID

    TRACE_BACK_FROM(best_goal) through sparse predecessor records
    NEGATE_AND_CHECK_PUBLIC_SCORE(best_goal.cost)
    INDEPENDENTLY_REPLAY_SCORE_AND_PATH()
    return alignment, search statistics, and memory evidence
```

Feature-isolation runs use the same implementation with:

```text
A*:       DFS=false, pruning=false
A*+DFS:   DFS=true,  pruning=false
POASTA:   DFS=true,  pruning=true
```

## Data structures

- existing typed and validated `PoaGraph<A>` with stable node identifiers;
- implicit source and sink cursors, never exposed as biological nodes;
- checked min/max remaining path-length vectors;
- deterministic reverse DFS postorder and Gartner candidate stack;
- output-sensitive node-to-bubble membership entries;
- one sorted reached-position vector per bubble exit, augmented lazily with
  insertion-prefix minima and a deletion range-minimum tree;
- sparse state records keyed by an injective checked-width encoding of
  `(cursor, query_prefix, affine_state)`;
- priority-indexed buckets through 65,536 plus a reversed binary-heap overflow
  frontier for larger exact priorities;
- a reused, fallibly grown DFS stack;
- compact predecessor state plus a one-byte tie transition for traceback;
- exact counters for discovered, queued, reopened, expanded, DFS-expanded,
  bubble-tested, bubble-pruned, stale, and peak-live states;
- checked modeled memory evidence for all algorithm-owned buffers.

Production code has no external engine, subprocess, FFI, unsafe block, global
mutable cache, or new dependency.

## Complexity and memory

Let `m=|Q|`, `n=|V|`, `e=|E|`, `R` be discovered alignment states, `X` expanded
states, and `K` node-to-bubble membership entries.

- graph validation and path-length metadata: `O(n+e)`;
- superbubble detection: `O(n+e)`;
- bubble materialization: output-sensitive as documented above;
- A* frontier operations: amortized constant time in the bounded bucket range and
  `O(log R)` for overflow priorities, with `O(R log R)` worst case;
- state relaxation: `O(m*e + m*n)` worst case, the same finite alignment graph
  bound as full dynamic programming;
- sparse state and traceback storage: `O(R)`;
- path metadata: `O(n)`;
- bubble index and reached exits: `O(K + B*m)` memory in the worst case, where
  `B` is the number of distinct bubble exits; each upper-bound query uses two
  binary searches, an `O(1)` insertion-prefix lookup, and an `O(log m)` deletion
  range-minimum query after a lazy `O(P)` rebuild for `P` reached positions;
- final traceback: `O(n+m)`.

The algorithm remains exact but has no universal subquadratic bound: an
uninformative heuristic or dissimilar inputs can discover most alignment states.
Measured state counts, modeled bytes, elapsed time, and peak RSS are mandatory
benchmark outputs.

## Public API and deterministic behavior

The public kernel is `Poasta`, constructed from `PoastaAffineCosts`. It borrows
an immutable `PoaGraph<A>` and `Sequence<A>` and returns an owned
`PoastaAlignment`. One generic implementation supports protein, DNA, and RNA.
No alphabet is guessed at runtime.

The default constructor enables all published accelerations on the sparse engine. A
clearly named options object may disable DFS or pruning for verification and
benchmarking, or select the exact dense affine DAG backend. Such a run reports its
active features/backend and does not claim the frozen profile identifier unless the
published sparse configuration is active.

`PoastaBackend::Auto` is an exact policy seam. The current validated policy first
attempts an admissible banded relaxed-linear certificate and then an exact affine
certificate. If neither proves a complete optimum, it continues with the same exact
sparse A* and match-DFS engine. The progressive adapter reuses the heuristic
allocation and row metadata between records. Its compact exact-affine band computes
the integer sublevel interval of the prefix/suffix path-length distance in closed
form, using `u128` translated coordinates to preserve checked `usize` behavior. The
exact cost buffers retain their initialized high-water length between alternating
diagonal and affine certificates. When both certificates run in one search step,
they also share one cached graph-prefix length-bound table. Every kernel still
initializes and accesses only its active logical prefix. Dense execution is never
selected silently; callers may request it explicitly for validation or workload
studies.

The publication fixes optimal costs but not a complete traceback order. This
profile freezes:

1. lower exact `g` always wins;
2. equal-cost `I` chooses extension before opening;
3. equal-cost `D` chooses extension before opening, then smaller predecessor ID;
4. equal-cost `M` chooses diagonal, then deletion close, then insertion close;
5. equal diagonal candidates choose the smaller predecessor ID;
6. equal goal costs choose the smaller sink ID;
7. graph successors, roots, and queue state keys use ascending stable IDs;
8. DFS pushes successors in reverse stable order so they are processed ascending;
9. a pruned equal-cost alternative cannot replace an already realized canonical
   state.

Rules 2-6 match the completed full-matrix oracle for alternatives retained in
the sparse search. Published DFS may deliberately suppress an equal-score path
before the later tie is observable, so the accelerated profile claims the same
optimum and deterministic output, not the same canonical traceback as the
full-matrix oracle. Rules 7-9 are explicit engineering resolutions absent from
the paper.

## Numeric and allocation policy

Sparse path costs and `f=g+h` values use `u64`. Inputs are non-negative `i64`
costs, and every addition and multiplication is checked. Any overflow is reported
as a typed error with a widened `u128` diagnostic value. Every optimum representable
by the public `i64` score contract necessarily fits in `u64`; conversion to the
negated public score is still checked explicitly.

Dimensions and modeled byte counts use checked `usize`. Major vectors, maps,
queues, and traceback outputs reserve fallibly before growth. No saturation,
wrapping, infinity sentinel inside a valid numeric cost, panic-based input
validation, or partial result is permitted.

## Independently testable implementation stages

1. Published-cost validation and score conversion.
2. Augmented graph cursors and min/max remaining distances.
3. Sparse ordinary A* without DFS or pruning.
4. Canonical sparse traceback and independent replay.
5. Exact-match DFS in isolation.
6. Reverse DFS postorder and Gartner Algorithm 2 superbubble detection.
7. Bubble memberships and min/max exit distances.
8. Reached-exit maps and implicit insertion/deletion upper bounds.
9. Affine-state pruning guard.
10. Full sparse POASTA kernel with all accelerations and compact frontier/state
    storage.
11. Exact dense affine DAG backend and differential agreement with sparse/scalar
    execution.
12. Iterative graph incorporation through existing Lee semantics.
13. Benchmarks, isolated black-box comparison, and scientific MSA metrics.

Every stage after item 3 must be differentially runnable against the previous
stage and the full-matrix oracle.

## Validation and acceptance criteria

The profile may be marked complete only when all of the following are proven:

- manual DNA, RNA, and protein cases cover identity, mismatch, empty inputs,
  leading/internal/trailing gaps, adjacent opposite gaps, repeated symbols,
  branches, recombination, multiple roots/sinks, and nested bubbles;
- every returned score and operation stream passes the independent affine replay;
- ordinary A*, A*+DFS, and full POASTA agree with the scalar oracle on exhaustive
  bounded linear sequences and independently generated small DAGs;
- random and property-generated DAGs exercise all valid cost tuples, alphabet
  types, deterministic repetition, node relabeling, and branch-order changes;
- a brute-force superbubble predicate agrees with the Gartner detector for every
  bounded generated DAG;
- enabling DFS never changes score, and tests demonstrate that it suppresses gap
  alternatives only for a newly unvisited exact target;
- enabling pruning never changes score, with explicit tests for no reached exit,
  multiple reached positions, implicit insertion, implicit deletion, boundary
  positions, nested bubbles, and the affine `I`/`D` guard;
- state-count tests demonstrate strictly fewer discovered/queued states for a
  published-style long match and strictly fewer expansions for a dominated
  bubble fixture, without using timing as a correctness assertion;
- the iterative adapter passes all existing ten SPOA 4.1.5 global-affine fixtures
  and adds static POASTA black-box observations when a pinned compatible release
  is available in the isolated sibling project;
- black-box score equality is reported separately from traceback, graph, SP, and
  TC equality because optimal ties can change progressive graph construction;
- scientific DNA, RNA, and protein datasets report SP and TC where references
  exist, with provenance and hashes;
- Criterion and isolated process benchmarks cover linear high identity, linear
  low identity, redundant graphs, branch-heavy graphs, nested bubbles, and graph
  growth, reporting time, modeled memory, peak RSS, and all feature-isolation
  modes;
- no expected value is generated by this implementation;
- repeated debug and release tests, Clippy with warnings denied, rustdoc with
  warnings denied, formatting, `git diff --check`, MSRV, and dependency-policy
  gates pass for the complete workspace;
- public code and documentation contain no placeholder, `TODO`, mock, hidden
  fallback to full-matrix DP, unsafe block, external runtime, FFI, or subprocess;
- limitations and every observed difference from the opaque reference remain
  documented.

## Completion evidence

The frozen profile satisfies the acceptance criteria above:

- all four feature modes agree with the scalar full-matrix oracle on exhaustive
  binary linear pairs and property-generated DAGs; an additional run completed
  5,000 cases per POASTA property without failure;
- the Gartner Algorithm 2 detector agrees with an independent formal predicate
  on every one of the 3,481 valid ordered DAGs through six vertices;
- 74 crate tests, two SPOA integration tests, and one POA 2.0 integration test
  pass in both debug and release configurations;
- all ten static SPOA 4.1.5 cases have exact scoreable objectives, graph
  cardinalities, rectangular alignments, SP/TC `1.0`, and source recovery;
- six BAliBASE protein cases, one Rfam RNA case, and one Ensembl transcript case
  are deterministic, source-preserving, and recorded with SP/TC, graph, search,
  memory, provenance, and file hashes;
- Criterion covers high and low identity, redundant, branch-heavy, nested-bubble,
  feature-isolation, and progressive-growth workloads; a separate 31-pair native
  process protocol records validated comparison and peak RSS against SPOA;
- Rust 1.85.0 workspace tests and doctests, workspace Clippy for all targets and
  features, warnings-denied rustdoc, formatting, source hygiene, `git diff
  --check`, and all cargo-deny advisory, ban, license, and source gates pass;
- the executable example completes and reports verified alignments for the common
  typed API.

The reviewed benchmark record, including negative performance evidence, is
`../benchmarks/poasta.md`. No criterion above depends
on an expected value generated by this implementation.

## Published omissions and explicit resolutions

- **Complete ties:** resolved by the deterministic order above; only the optimal
  score is claimed publication-equivalent. Exact-match DFS can suppress an
  equal-cost gap-first path, so its canonical alignment may differ from both the
  scalar oracle and an opaque reference without constituting a score error.
- **Empty graph:** the implicit `nu -> tau` edge yields one insertion run, matching
  the completed scalar oracle.
- **Disconnected DAG components:** `nu` connects all roots and `tau` all sinks;
  one complete component path is selected globally.
- **DFS visited definition:** a target discovered at equal or lower cost is treated
  as visited, so gap alternatives are retained; only a strictly new or improved
  zero-cost exact target permits suppression, and only for its own adjacent
  deletion. A mismatch or already-reached sibling also retains insertion from
  the shared source state.
- **DFS exit recording:** settled `M` states from either frontier are recorded, as
  justified above, and isolated differential tests are mandatory.
- **Bubble distance overflow:** all lengths and endpoints are checked; overflow is
  a typed error rather than a truncated index.
- **Implicit deletion feasibility:** it is limited by the maximum remaining real
  path length after the bubble exit; infeasible bounds are ignored.
- **Nested bubbles:** every containing bubble is tested. A state is pruned when any
  one bubble supplies complete dominance evidence.
- **Range beyond the query:** only valid exit positions in `0..=m` are compared;
  an empty feasible interval cannot prune.
- **Duplicate reached position:** retain its lowest realized match-state cost.
- **Repeated iterative validation:** a public stand-alone call validates its DAG
  once. The transactional iterative adapter reuses the topological order returned
  by the just-completed graph update and consumes the verifier report returned by
  the accelerated kernel, avoiding duplicate validation and replay without
  weakening either public boundary.
- **Memory definition:** both modeled owned-buffer bytes and process peak RSS are
  reported; neither is presented as the other.
- **Progressive update:** the accelerator changes only sequence-to-graph search.
  Existing Lee node-fusion and path-incorporation semantics remain unchanged.

## Differences from the published program

The profile is reconstructed from publications rather than program source. Known
API-level differences are intentional:

- biological alphabets and sequences are statically typed;
- sentinels are implicit and stable biological node IDs are never renumbered;
- exact numeric, allocation, traceback, and verification failures are typed;
- a deterministic total tie order is public;
- search statistics and modeled memory are part of the result evidence;
- feature-isolation modes exist for falsifiable differential validation;
- the first implementation is scalar and single-threaded.

Any later behavioral difference discovered by opaque black-box testing must be
recorded here; it must not be silently patched by deriving behavior from source.

## Limitations

- global sequence-to-DAG alignment only;
- one affine cost function shared by both gap directions;
- literal zero-cost matches and one non-negative mismatch cost;
- one canonical optimum, not enumeration of all optimal paths;
- acyclic, residue-labeled POA graphs only;
- no arbitrary string-labeled or cyclic sequence graph;
- no local, semiglobal, overlap, free-end, convex, or dual-affine mode;
- no SIMD, GPU, or parallel traversal in the initial profile;
- worst-case search still visits the full alignment graph;
- superbubble preprocessing can be output-heavy for deeply nested bubbles;
- high-diversity inputs may favor contiguous full-matrix SIMD DP, as reported by
  the paper; any automatic selector requires a separately specified and tested
  policy.
