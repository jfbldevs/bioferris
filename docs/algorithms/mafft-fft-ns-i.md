# MAFFT FFT-NS-i

## Reproduction target

This document specifies the immutable scientific profile
`MAFFT-2002-FFT-NS-i-v2`.

The target is the iterative-refinement method introduced by Katoh, Misawa,
Kuma, and Miyata in 2002. Its initial alignment is the project's separately
versioned `MAFFT-2002-FFT-NS-2-v1` profile. Refinement reconstructs the
best-first, tree-dependent procedure cited by that paper from Hirosawa et al.
and uses Gotoh's factorized weighted sum-of-pairs (WSP) objective.

This is not a byte-for-byte reconstruction of an undocumented historical or
current MAFFT executable. In particular, current MAFFT has changed matrices,
gap parameters, distance corrections, tree construction, and refinement
engineering since 2002. The production crate contains no MAFFT source,
translated source, FFI, subprocess invocation, or executable dependency.
MAFFT 7.526 remains an external comparator in an isolated validation
laboratory; the v2 stopping rule was also checked independently against the
released controller source.

## Primary references

1. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic Acids
   Research. 2002;30(14):3059-3066.
   DOI: [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
   This is the primary MAFFT source. It defines FFT-NS-i as an FFT-NS-2
   alignment followed by repeated group division and realignment, using
   tree-dependent restricted partitioning, until no better score is found.
2. Hirosawa M, Totoki Y, Hoshida M, Ishikawa M. *Comprehensive study on
   iterative algorithms of multiple sequence alignment*. Computer
   Applications in the Biosciences. 1995;11(1):13-18.
   DOI: [10.1093/bioinformatics/11.1.13](https://doi.org/10.1093/bioinformatics/11.1.13).
   This paper supplies the best-first iteration and the `2N-3` tree-edge
   partitions. It states that a guide tree is constructed from the current MSA
   at the beginning of each cycle.
3. Gotoh O. *A weighting system and algorithm for aligning many
   phylogenetically related sequences*. Computer Applications in the
   Biosciences. 1995;11(5):543-551.
   DOI: [10.1093/bioinformatics/11.5.543](https://doi.org/10.1093/bioinformatics/11.5.543).
   This paper defines WSP, edge-product pair weights, the three-way weighting
   approximation, and rigorous factorization across a tree-edge split.
4. Gotoh O. *Optimal alignment between groups of sequences and its application
   to multiple sequence alignment*. Computer Applications in the Biosciences.
   1993;9(3):361-370.
   DOI: [10.1093/bioinformatics/9.3.361](https://doi.org/10.1093/bioinformatics/9.3.361).
   This paper supplies the group-to-group dynamic-programming basis.
5. Berger MP, Munson PJ. *A novel randomized iterative strategy for aligning
   multiple protein sequences*. Computer Applications in the Biosciences.
   1991;7(4):479-484.
   DOI: [10.1093/bioinformatics/7.4.479](https://doi.org/10.1093/bioinformatics/7.4.479).
   This is an earlier iterative realignment source cited by Katoh et al.

Inherited FFT, matrix, tuple, guide-tree, and gap-scoring references are listed
in [mafft-fft-ns-1.md](mafft-fft-ns-1.md) and
[mafft-fft-ns-2.md](mafft-fft-ns-2.md).

Official operational evidence:

6. The official MAFFT algorithm page defines FFT-NS-i as FFT-NS-2 followed by
   WSP iterative refinement, stopping when WSP no longer improves or when the
   configured cycle limit is reached:
   <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.
7. The official manual identifies
   `mafft --retree 2 --maxiterate 1000 INPUT` as the operational FFT-NS-i
   selection:
   <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.
8. The official source archive exposes the released best-first controller in
   `core/tditeration.c`. After retaining a strictly improving candidate, it
   compares that candidate's score with earlier accepted-score history and
   stops on exact equality. The source archive is published at
   <https://mafft.cbrc.jp/alignment/software/source.html>. MAFFT 7.526's pinned
   `dvtditr` binary contains the corresponding `Oscillating?`, `Converged`, and
   cycle-limit diagnostics.

Items 6-8 identify the current black-box command and its 1,000-cycle bound.
They do not override the 2002 scientific scoring profile.

## Objective and characteristic pipeline

FFT-NS-i corrects progressive-alignment decisions that FFT-NS-2 cannot undo.
Its characteristic pipeline is:

1. construct a complete FFT-NS-2 initial alignment;
2. infer a guide tree from the current MSA at the start of a refinement cycle;
3. derive factorized phylogenetic pair weights from that tree;
4. enumerate every unique split induced by an edge of the unrooted tree;
5. remove columns that are all gaps within either resulting subgroup;
6. independently realign each pair of subgroup profiles with the published FFT
   profile aligner;
7. score every candidate with the same cycle's full WSP objective;
8. replace the current MSA with the strictly highest-scoring candidate;
9. repeat until no candidate strictly improves WSP, an accepted best-candidate
   score repeats an earlier best-first history score, or 1,000 accepted cycles
   have been completed.

Stopping after step 1 would silently substitute FFT-NS-2. Testing only one
split would substitute a greedy or random refinement. Rebuilding the two
subalignments from their ungapped sequences would discard their fixed internal
homology and would also be a different algorithm.

## Best-first oscillation stopping

Tree reconstruction changes pair weights between cycles, so strict improvement
under each cycle's fixed objective does not prove monotonic improvement under
one global objective. A deterministic alignment can therefore alternate among
states while every individual cycle reports a positive gain.

The released MAFFT best-first controller prevents this by retaining accepted
candidate scores. After accepting and retaining a candidate, it stops if that
score is exactly equal to an earlier accepted score. The C controller starts
the comparison at history entry one rather than zero; this profile preserves
that boundary and exact `f64` equality. The repeated candidate remains the
returned alignment. No tolerance, hash approximation, or unaccepted candidate
participates in the decision.

Version 1 omitted this native stopping rule and could consume the entire
1,000-cycle limit on an oscillation. Version 2 adds the explicit
`FftNsiTermination::Oscillating` state and the released stopping semantics. The
change is versioned because it can alter the number of audited cycles and, when
the configured cycle limit does not land on the repeated state, the returned
alignment.

## Initial alignment

The initial MSA is exactly `MAFFT-2002-FFT-NS-2-v1`, including its complete
FFT-NS-1 evidence, aligned Poisson distance matrix, second UPGMA tree, and
second progressive pass. A direct call to FFT-NS-i with a zero refinement
limit is not part of this immutable profile; the published profile always
permits refinement. Tests may exercise internal components independently.

## Cycle guide tree

Hirosawa et al. require a guide tree constructed from the current MSA at the
beginning of every best-first cycle. This profile computes the same explicit
pairwise-deletion Poisson distance used by `MAFFT-2002-FFT-NS-2-v1`:

```text
C_ij = { c | row i and row j both contain residues at column c }
q_ij = identical residues in C_ij / |C_ij|
D_ij = -ln(max(q_ij, 10^-12)).
```

Two empty rows have distance zero. A non-empty pair with no comparable
residues has the maximum finite distance. Exact typed-symbol equality defines
identity. Deterministic size-weighted UPGMA constructs the cycle tree.

The rooted UPGMA representation has a degree-two root. For partitioning and
Gotoh weights it is converted to an unrooted binary tree by suppressing that
artificial root and replacing its two incident edges with one edge whose
length is their sum. The resulting tree has exactly `2N-3` edges for `N >= 2`.
For two sequences it contains one edge; for a singleton it contains none.

## Tree-dependent partitions

Deleting each unrooted edge partitions the leaf indices into non-empty sets
`L` and `R`. Complementary descriptions of the same split are canonicalized:
the side with fewer leaves is stored first; equal-size sides are ordered
lexicographically. Splits are then evaluated in lexicographic order of their
canonical leaf lists. The order cannot change the selected numerical maximum,
but it supplies a deterministic tie rule.

For each subgroup, columns containing only gaps in that subgroup are removed
before realignment. Every retained row otherwise remains unchanged, so the
homology relationships internal to each subgroup are fixed during this
candidate operation.

## Factorized Gotoh pair weights

For an unrooted tree `T`, Gotoh constrains the weight of the leaf-to-leaf path
`p=(j,k)` to factor over its edges:

```text
w_jk = product over e in path(j,k) of u_e.
```

This makes every pair weight crossing a cut edge `e=(x,y)` separable:

```text
w_jk = [product over path(j,x) of u] * u_e
       * [product over path(y,k) of u],  j in L, k in R.
```

The common positive cut factor does not affect the profile-DP optimum.
Therefore the two subgroup profile weights are the normalized products from
each leaf to its endpoint of the cut edge.

### Three-way approximation

At an internal node, each incident direction is replaced by an effective
branch length. If the direction reaches a leaf through an edge of length `l`,
its effective length is `l`. If it reaches another internal node whose two
remaining effective branches are `r` and `s`, it is

```text
l + (r parallel s),       r parallel s = r*s/(r+s).
```

This is Gotoh's synthetic-resistance construction. For effective branch
lengths `a`, `b`, and `c`, define

```text
S    = ab + ac + bc
p_ab = c(a+b)/(2S)
p_ac = b(a+c)/(2S)
p_bc = a(b+c)/(2S).
```

The local factors that exactly factor these three pair weights are

```text
t_a = sqrt(p_ab * p_ac / p_bc)
t_b = sqrt(p_ab * p_bc / p_ac)
t_c = sqrt(p_ac * p_bc / p_ab).
```

The adjustable equalization exponent is applied as `t_direction^F`. A terminal
tree edge receives the factor at its one internal endpoint. An internal edge
receives the product of its two flanking local factors. Pair weights are then
calculated as path products and normalized so that

```text
sum_(j<k) w_jk = N(N-1)/2.
```

Thus equal weights would all be one. This profile fixes `F=1.15`, the value
used for the reported alignment in Gotoh's paper. The paper reports broad
behavior over `1.0 <= F <= 1.3` but Katoh et al. do not print the MAFFT 2002
executable's value. Choosing and versioning `1.15` is therefore an explicit
resolution, not a claim about hidden source behavior.

UPGMA may produce zero-length edges, whereas Gotoh's derivation assumes
strictly positive branch lengths. Before effective-resistance calculation,
each edge is floored to `10^-12`. If every inferred edge is zero, all pair
weights are exactly one. Products are accumulated in logarithmic space and
normalized with a maximum shift to prevent underflow.

## Weighted sum-of-pairs objective

For an MSA `A` and its induced pairwise alignments `A_jk`,

```text
WSP(A) = sum_(j<k) w_jk * S(A_jk).
```

Pairwise alignment extraction first removes columns that are gap-gap for that
pair. Its score is

```text
S(A_jk) = sum over residue-residue columns of M_hat(a,b)
          - S_op * number_of_maximal_residue-gap_runs.
```

The orientation of a gap run is significant: a gap in `j` followed by a gap
in `k` begins a second run. The inherited 2002 defaults are `S_op=2.4` and
`S_a=0.06`; proteins use normalized JTT 200-PAM and DNA/RNA use normalized
K80 200-PAM. This induced-pair definition is consistent with a one-row
application of the publication's extension-free profile recurrence.

Katoh et al. do not publish an executable formula for scoring terminal gaps,
gap-gap columns, ambiguity symbols, or adjacent existing and newly introduced
gaps in the refinement objective. The rules above are the versioned resolution
used both by candidate acceptance and by an independent test oracle.

The guide tree and weights are fixed while all candidates in one cycle are
compared. They are reconstructed only after a candidate is accepted. Therefore
`score_before` and `score_after` within a cycle are directly comparable, while
raw WSP values from different cycles need not be monotone because their
phylogenetic weights may change.

### Exact contribution reuse

Realigning a tree-edge split cannot change either subgroup's induced pairwise
alignments. At the start of a cycle, the implementation therefore caches every
weighted pair contribution from the current MSA. For a candidate it recomputes
only pairs crossing that candidate's split and reuses cached contributions for
pairs on the same side.

This does not use a delta approximation and does not change floating-point
selection. Candidate WSP is accumulated in the original row-major pair order,
choosing either the cached or recomputed contribution at each position. It is
therefore bit-for-bit equal to a complete direct evaluation after subgroup
homology has been validated. Debug assertions and a separate direct oracle
test this invariant.

## Candidate profile realignment

For a cut edge, each subgroup profile receives the cut-factorized leaf weights
described above. The two profiles are realigned with the same deterministic FFT
correlation, homologous-segment chaining, normalized substitution matrix, and
position-specific extension-free dynamic programming used by FFT-NS-1 and
FFT-NS-2. The FFT heuristic therefore remains a defining part of every
refinement candidate.

The candidate's profile operations are replayed independently before it may be
scored. Its rows are restored to original input order and any global all-gap
column is removed. The final candidate must consume every original typed symbol
exactly once and preserve every within-subgroup induced pairwise alignment.

## Best-first acceptance and stopping

All candidates in a cycle are generated from the same current MSA. The
candidate with the greatest finite WSP is selected only if

```text
candidate_WSP > current_WSP.
```

No epsilon or relative tolerance weakens this strict published improvement
condition. Every WSP is checked finite before ordinary IEEE-754 `>` comparison.
Equal scores retain the lexicographically earlier canonical split but are never
accepted over an equal-scoring current MSA.

One accepted best candidate is one refinement cycle. The algorithm stops after
a full cycle with no strict improvement or after 1,000 accepted cycles. The
limit matches the official FFT-NS-i operational maximum and prevents an
unbounded public call if finite-precision reweighting visits a long sequence of
strictly increasing states.

## Deterministic candidate parallelism

The scientific default remains one worker. `with_thread_count(T)` is an
execution setting, not a scientific parameter: it partitions the canonical
candidate list into at most `T` fixed contiguous, nonempty ranges. Each scoped
worker evaluates its range in order and retains only its first maximum. Worker
winners are then reduced in range order with the same strict `>` comparison as
the sequential loop.

No candidate score reduction is parallelized internally. WSP terms retain their
original row-major addition order, and the global candidate reduction retains
canonical split order. Consequently thread scheduling cannot change any score
bit, selected split, cycle decision, merge evidence, or returned alignment.
Zero workers is rejected; values greater than the candidate count are capped by
the available `2N-3` independent tasks.

## Opt-in performance diagnostics

`align_profiled` retains the shared coarse complete-call phases.
`align_detailed_profiled` additionally measures aggregate candidate worker work
for projection, profile preparation, FFT/anchor search, profile DP, profile
assembly, checked alignment materialization, validation, WSP scoring, and
accounting. Candidate work is summed over candidates, so with multiple workers
it is not a second complete-call wall clock. A nested profile-DP scope further
separates workspace initialization, recurrence fill, endpoint/traceback,
accounting, operation-score replay, and residual orchestration, and reports the
number of non-empty cells evaluated.

The ordinary `align` path and coarse `align_profiled` path instantiate the
candidate evaluator without inner clocks. Timing values never participate in
scores, tie resolution, cycle acceptance, output construction, or the profile
identifier.

## Tie-breaking

The papers do not publish all tree, candidate, FFT, or DP ties. This profile
uses the following total deterministic order:

1. aligned distances are evaluated in input-index order;
2. UPGMA uses distance, then minimum descendant indices;
3. unrooted splits use smaller cardinality, then lexicographic leaf indices;
4. candidate selection uses larger WSP, then canonical split order;
5. FFT, segment, profile-DP, and traceback ties use the complete ordering of
   `MAFFT-2002-FFT-NS-1-v1`;
6. rows are always returned in original input order.

No random number, hash iteration order, thread scheduling, locale, clock, or
external executable affects a library result.

## Pseudocode

```text
FFT_NS_I(records):
    validate non-empty collection and unique identifiers
    initial <- FFT_NS_2(records)
    current <- initial.alignment
    cycles  <- []

    repeat up to 1000 accepted cycles:
        distances <- aligned_pairwise_deletion_poisson(current)
        rooted     <- deterministic_UPGMA(distances)
        unrooted   <- suppress_degree_two_root(rooted)
        edge_logs  <- Gotoh_three_way_factors(unrooted, F=1.15)
        pair_w     <- normalized_path_products(edge_logs)
        pair_cache <- weighted pair contributions of current in pair order
        before     <- SUM_IN_PAIR_ORDER(pair_cache)
        best       <- NONE

        splits <- canonical_edge_splits(unrooted)
        ranges <- FIXED_CONTIGUOUS_RANGES(splits, min(T, length(splits)))
        range_winners <- array(length(ranges))

        for range in ranges, optionally on scoped workers:
            range_best <- NONE
            for split in range:
                left  <- project current to split.left; remove all-gap columns
                right <- project current to split.right; remove all-gap columns
                wl    <- normalized leaf-to-cut products(edge_logs, split.left)
                wr    <- normalized leaf-to-cut products(edge_logs, split.right)

                candidate <- FFT_PROFILE_ALIGN(left, right, wl, wr)
                independently replay candidate profile operations
                restore original row order; validate exact symbol consumption
                score <- SUM_IN_PAIR_ORDER(
                    recomputed contribution for pairs crossing split,
                    cached contribution for pairs within either subgroup
                )

                if score is greater than range_best score, retaining split order on ties:
                    range_best <- candidate
            range_winners[range.index] <- range_best

        best <- REDUCE_IN_CANONICAL_ORDER_WITH_STRICT_GREATER(range_winners)

        if best.score <= before:
            record convergence cycle and stop

        current <- best.alignment
        record accepted cycle, split, before, after, and merge evidence

    reconstruct final tree, factors, pair weights, and final WSP from current
    return initial, cycles, final evidence, current, and memory diagnostics
```

## Data structures

- the complete immutable `FftNs2Alignment` audit record;
- a temporary cycle-specific aligned distance matrix and a retained validated
  rooted guide tree;
- an internal compact unrooted adjacency tree with positive effective lengths;
- logarithmic edge factors and a dense symmetric pair-weight matrix;
- a dense cycle cache of fixed weighted pair contributions;
- canonical bipartitions with leaf-to-cut profile weights;
- subgroup `WorkingProfile` rows preserving typed optional symbols;
- existing FFT signals, segment chains, profile scores, and compact traceback;
- one retained accepted-candidate record per cycle plus the terminal rejected
  cycle summary;
- one best candidate per active scoped worker and a conservative concurrent
  worker high-water estimate;
- a final `MultipleAlignment` and explicit major-buffer memory accounting.

## Complexity

Let `N` be the sequence count, `W` the current alignment width, `P(W)` the
cost of one FFT-assisted profile/profile realignment, and `C` the number of
accepted cycles.

- aligned distance reconstruction per attempted cycle: `O(N^2 W)`;
- scalar UPGMA per attempted cycle: `O(N^3)`;
- unrooting, effective lengths, and edge factors: `O(N^2)` in the clear scalar
  implementation (`O(N)` is possible with rerooting messages);
- dense pair weights and one WSP evaluation: `O(N^2 W)`;
- candidate count per attempted cycle: `2N-3`;
- candidate profile realignment: `O(N * P(W))` over all `2N-3` splits;
- candidate WSP scoring across all tree edges:
  `O(N^3 + W * sum_e |L_e||R_e|)`, where the sum is over edge splits;
  this is typically `O(N^2 W log N)` for a balanced tree and remains
  `O(N^3 W)` in the worst unbalanced case;
- dense pair-weight memory: `O(N^2)`;
- cached weighted pair-contribution memory: `O(N^2)`;
- one profile traceback: `O(W^2)` worst case;
- with `T` active candidate workers, worst-case concurrent candidate storage is
  `O(TW^2 + TNW)` and ideal candidate wall time is divided by at most
  `min(T, 2N-3)`; no linear speedup is assumed;
- returned initial/final alignments and cycle evidence: `O(NW + CN)` plus
  retained merge diagnostics.

The one-worker mode remains the auditable correctness baseline. Exact reuse,
sparse profile frequencies, rolling score rows, packed traceback, and optional
candidate parallelism preserve its arithmetic and decision order. Linear-time
tree messages and SIMD remain deferred.

## Published omissions and explicit resolutions

| Omission or ambiguity | `v2` resolution |
| --- | --- |
| Which Hirosawa search variant FFT-NS-i uses | Best-first over all `2N-3` tree-edge splits, the cited procedure whose stopping rule matches Katoh's text. |
| Tree reconstruction frequency | From the current MSA at the start of every cycle, as stated by Hirosawa et al. |
| Aligned distance details | FFT-NS-2's pairwise-deletion Poisson distance with a `10^-12` floor. |
| Rooted versus unrooted edges | Suppress the degree-two UPGMA root and sum its two edge lengths. |
| Pair-weight method | Gotoh three-way edge-product approximation. |
| Equalization exponent used by MAFFT 2002 | Fixed and versioned at Gotoh's reported `F=1.15`. |
| Zero tree branches | Floor to `10^-12`; an entirely zero tree uses equal weights. |
| Subgroup treatment | Preserve internal alignment and remove subgroup-all-gap columns. |
| WSP treatment of gap-gap columns | Delete them from the induced pairwise alignment. |
| WSP gap cost | One `S_op` per maximal orientation-specific residue-gap run, including terminal runs. |
| Candidate improvement tolerance | Strict finite `>` with no epsilon. |
| Equal candidate scores | Canonical split order. |
| Best-first oscillation | After retaining a strictly improving candidate, stop when its exact score repeats any accepted history entry after entry zero. |
| Maximum cycles | 1,000 accepted cycles. |
| Candidate worker count | One by default; an explicit positive count uses fixed contiguous ranges and canonical reduction without changing results. |
| Singleton input | Return FFT-NS-2 unchanged with a converged zero-candidate terminal cycle. |

Every resolution is part of the profile identifier. Changing any one requires
a new profile version and new validation evidence.

## Validation and acceptance criteria

FFT-NS-i can be marked complete only when all of the following hold:

1. this document contains the references, equations, original pseudocode, data
   structures, complexity, parameters, ties, omissions, and limitations;
2. the retained initial audit is bit-for-bit equal to direct FFT-NS-2;
3. hand-computed trees cover root suppression and all `2N-3` unique splits;
4. hand-computed three-leaf and four-leaf cases verify effective resistance,
   edge-product factorization, normalization, symmetry, and positivity;
5. a separately written direct path oracle verifies every pair weight;
6. hand-computed WSP cases cover matches, substitutions, leading, internal,
   terminal, adjacent opposite-orientation gaps, and gap-gap deletion;
7. every candidate operation stream is independently replayed and every WSP is
   recomputed by a separate test oracle;
8. accepted cycles strictly improve their own fixed-weight objective, and
   termination is auditable as no-improvement convergence, an exact repeated
   best-first score, or the documented cycle limit;
9. within-subgroup pairwise homologies and all original typed sequences are
   preserved exactly;
10. unit and property tests cover protein, DNA, RNA, empty sequences,
    singleton, identical, divergent, repetitive, ambiguity-rich, leading-gap,
    internal-gap, and terminal-gap cases;
11. output and all exposed diagnostics are deterministic, and one-worker versus
    multi-worker execution is bit-identical;
12. an inert fixture generated by official MAFFT 7.526 checks exact trivial
    cases and SP/TC agreement on nontrivial cases;
13. a fixed external BAliBASE 3 subset records SP and TC without vendoring the
    scientific dataset;
14. Criterion reports sequence-length, sequence-count, convergence, direct
    FFT-NS-2 versus FFT-NS-i costs, and candidate-thread scaling;
15. matched fresh-process runs compare a Rust worker with official MAFFT 7.526
    and record binary, package, harness, and input hashes;
16. formatting, all-target Clippy, workspace tests, release tests, Rustdoc, and
    dependency license/source/advisory checks pass;
17. the public API and example are English, safe, deterministic, and contain no
    unfinished route, mock, placeholder, or runtime reference dependency.

## Acceptance evidence

The predecessor `MAFFT-2002-FFT-NS-i-v1` profile satisfied the original
criteria on 2026-08-04:

- component/unit tests plus integration/property/reference suites
  cover the equations, tree splits, factorization, WSP, all three alphabets,
  boundary cases, exact source replay, exact cached-versus-direct scoring,
  strict improvement, convergence, and deterministic parallel reduction; the
  three property tests each execute 64 generated cases;
- the complete workspace all-target/all-feature test matrix and the
  `align-mafft` release all-target matrix pass with zero failures;
- the pinned nine-case MAFFT 7.526 fixture passes mean SP/TC gates with
  `0.830670290`/`0.723853659` and three exact alignments;
- the six-family BAliBASE subset records macro SP/TC
  `0.888317779`/`0.564488115` and micro SP/TC
  `0.851233015`/`0.782608696`;
- Criterion covers length, sequence count, an FFT path, proteins, direct
  FFT-NS-2 cost, a genuinely accepted refinement, and a 1-versus-16-worker
  sweep; 16 workers give a measured `4.76x` in-process speedup at 16-by-256;
- at matched 16-worker settings, static Rust is 62.80x, 24.10x, and 3.84x
  faster than official MAFFT at 4-by-64, 8-by-128, and 16-by-256;
- formatting, Clippy with warnings denied, Rustdoc with warnings denied, and
  Cargo dependency advisory/license/source/bans checks pass.

Version 2 adds a minimized five-sequence oscillation fixture, exact history
boundary tests, deterministic serial/parallel termination checks, and a frozen
post-holdout scaling campaign. On its pathological 20-by-352 tuning row, the
native stopping rule reduces accepted cycles from 1,000 to eight and fresh
process time by 99.14% while retaining the predecessor's exact output SHA-256.
The official MAFFT workflow terminates normally on the same input. WSL timings
remain development diagnostics pending the two-host bare-metal campaign.

Detailed measurements and immutable evidence are in
mafft-fft-ns-i.md and
mafft-fft-ns-i-vs-reference.md.

## Black-box reference boundary

The isolated comparator is the official Linux amd64 MAFFT 7.526 package. The
pinned command is:

```text
mafft --retree 2 --maxiterate 1000 --thread T --inputorder --quiet INPUT
```

The package is opaque: its source is neither downloaded nor inspected. Output
FASTA and timing observations are inert development evidence and are never
loaded by product code. Exact alignment equality is required only for
unambiguous trivial cases. Nontrivial comparison uses valid symbol recovery,
SP, TC, determinism, and runtime because the current executable differs from
the 2002 profile in documented and undocumented ways.

## Known limitations of this profile

- It reconstructs the published 2002 scientific pipeline, not current MAFFT
  defaults or every behavior of an unavailable historical executable.
- Katoh et al. omit the exact equalization exponent; `F=1.15` is an explicit
  value reported by Gotoh, not source-derived evidence about MAFFT.
- The aligned-distance, zero-edge, ambiguity, induced-pair gap, terminal-gap,
  and tie rules include explicit totalizations of publication omissions.
- Cross-partition WSP scoring remains scalar within each candidate and can
  dominate realignment for many or highly divergent sequences.
- Tree effective lengths are computed independently per directed edge in the
  first scalar version rather than by an optimized rerooting pass.
- UPGMA, FFT, profile DP, and each candidate's WSP kernel remain scalar. Only
  independent candidates are optionally parallel, using newly scoped workers
  in every attempted cycle; small inputs can therefore be slower with many
  workers. SIMD is not implemented.
- FFT partitioning remains heuristic; a candidate profile DP can still occupy
  a full quadratic rectangle.
- The result retains the complete FFT-NS-2 audit plus refinement evidence, so
  its memory use is higher than an executable that discards intermediate data.
- Parallel execution retains one winning candidate per worker and increases
  peak memory in exchange for throughput.
