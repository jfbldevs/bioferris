# MAFFT G-INS-i

## Reproduction target

This document specifies the immutable scientific profile
`MAFFT-2005-G-INS-i-v1`.

The target is the ancestral G-INS-i strategy introduced in MAFFT 5.3 by
Katoh, Kuma, Toh, and Miyata in 2005. It is a global-consistency multiple
sequence aligner: all sequence pairs are globally aligned, their gap-free
segments form a primary consistency library, that library contributes to every
progressive profile alignment, and the resulting G-INS-1 alignment is refined
against a fixed `WSP+I` objective.

This profile is deliberately versioned. Current MAFFT documentation identifies
the 2005 description as an ancestral version, and current executables have
accumulated parameter and engineering changes. This is therefore not a claim
of byte-for-byte equivalence with either an unavailable MAFFT 5.3 executable or
MAFFT 7.526. No MAFFT source code, translated code, FFI, subprocess, or runtime
executable is used by the Rust implementation. MAFFT 7.526 is permitted only as
an opaque comparator in the isolated validation laboratory.

## Primary references

1. Katoh K, Kuma K, Toh H, Miyata T. *MAFFT version 5: improvement in
   accuracy of multiple sequence alignment*. Nucleic Acids Research.
   2005;33(2):511-518.
   DOI: [10.1093/nar/gki198](https://doi.org/10.1093/nar/gki198).
   This is the primary definition of MAFFT 5.3, G-INS-i, the primary
   consistency library, the importance equations, `W^I=2.7`, and `WSP+I`.
2. Katoh K, Kuma K, Miyata T, Toh H. *Improvement in the Accuracy of
   Multiple Sequence Alignment Program MAFFT*. Genome Informatics.
   2005;16(1):22-33.
   DOI: [10.11234/gi1990.16.22](https://doi.org/10.11234/gi1990.16.22).
   This publication by the authors restates the complete [GL]-INS-i pipeline
   and equations in a more explicit algorithmic sequence.
3. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic
   Acids Research. 2002;30(14):3059-3066.
   DOI: [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
   This defines the FFT-assisted global pair alignment, normalized matrix,
   position-specific extension-free gap cost, and profile DP inherited by the
   2005 method.
4. Thompson JD, Higgins DG, Gibson TJ. *CLUSTAL W: improving the sensitivity
   of progressive multiple sequence alignment through sequence weighting,
   position-specific gap penalties and weight matrix choice*. Nucleic Acids
   Research. 1994;22(22):4673-4680.
   DOI: [10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
   This supplies the progressive-stage tree weights cited by Katoh et al.
5. Gotoh O. *A weighting system and algorithm for aligning many
   phylogenetically related sequences*. Computer Applications in the
   Biosciences. 1995;11(5):543-551.
   DOI: [10.1093/bioinformatics/11.5.543](https://doi.org/10.1093/bioinformatics/11.5.543).
   This supplies the iterative-stage factorized pair weights and WSP score.
6. Hirosawa M, Totoki Y, Hoshida M, Ishikawa M. *Comprehensive study on
   iterative algorithms of multiple sequence alignment*. Computer
   Applications in the Biosciences. 1995;11(1):13-18.
   DOI: [10.1093/bioinformatics/11.1.13](https://doi.org/10.1093/bioinformatics/11.1.13).
   This defines tree-dependent restricted partitioning.
7. Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein
   blocks*. Proceedings of the National Academy of Sciences.
   1992;89(22):10915-10919.
   DOI: [10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915).
   This is the source of BLOSUM62.
8. Kimura M. *A simple method for estimating evolutionary rates of base
   substitutions through comparative studies of nucleotide sequences*.
   Journal of Molecular Evolution. 1980;16:111-120.
   DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581).
   This supplies the nucleotide substitution model.
9. Needleman SB, Wunsch CD. *A general method applicable to the search for
   similarities in the amino acid sequence of two proteins*. Journal of
   Molecular Biology. 1970;48(3):443-453.
   DOI: [10.1016/0022-2836(70)90057-4](https://doi.org/10.1016/0022-2836(70)90057-4).

Official operational evidence:

10. The official algorithm page describes G-INS-i as a global method for
    sequences whose complete regions are alignable and defines the modern
    opaque command `mafft --globalpair --maxiterate 1000 INPUT`:
    <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.
11. The official MAFFT site explicitly labels the 2005 paper as describing
    ancestral G-INS-i, L-INS-i, and E-INS-i strategies:
    <https://mafft.cbrc.jp/alignment/software/>.

Items 10-11 establish the reference boundary and intended use. They do not
replace the version 5.3 scientific profile with current executable behavior.

## Objective and characteristic pipeline

G-INS-i targets sets of globally homologous protein, DNA, or RNA sequences of
similar domain organization. Its characteristic pipeline is:

1. compute every unordered global pair alignment with the 2002 FFT
   approximation and the 2005 scoring profile;
2. derive a normalized score distance for every pair;
3. build one mixed-linkage UPGMA-like guide tree;
4. compute Thompson/CLUSTALW sequence weights from that tree;
5. split every pair alignment into maximal gap-free segments;
6. derive residue frequencies, segment importance values, and the sparse
   primary importance matrix `I`;
7. progressively align profiles along the fixed guide tree with full scalar
   profile DP and the combined substitution-plus-importance matrix, producing
   G-INS-1;
8. derive Gotoh pair weights and all tree-edge partitions from the same fixed
   tree;
9. repeatedly realign the two fixed subalignments of every partition and
   accept only the strictly best improvement in `WSP+I`;
10. stop when a complete cycle has no strict improvement or after 1,000
    accepted cycles.

Replacing the all-pairs stage with tuple counts, rebuilding the guide tree,
using only WSP, applying library extension, enabling FFT in the progressive
stage, or returning G-INS-1 without refinement would be different algorithms.

## Versioned parameters

| Parameter | Value | Scientific role |
| --- | ---: | --- |
| Protein matrix | BLOSUM62 | Best of the five matrices in the 2005 study |
| Nucleotide matrix | K80, transition/transversion ratio 2 | Published MAFFT nucleotide model |
| `S^op` | `1.53` | Extension-free gap-opening cost |
| `S^a` | `0.123` | Normalized-matrix offset |
| `W^I` | `2.7` | Importance-library multiplier |
| Mixed-linkage coefficient `s` | `0.1` | MAFFT UPGMA-like linkage |
| FFT correlation peaks | `20` | Inherited 2002 pairwise approximation |
| Sliding window | `30` | Inherited 2002 pairwise approximation |
| Segment threshold | `0.7` per column | Inherited 2002 pairwise approximation |
| Maximum FFT segment length | `150` | Inherited 2002 pairwise approximation |
| Gotoh equalization exponent | `1.15` | Explicit resolution already used by the versioned FFT-NS-i profile |
| Zero branch floor | `10^-12` | Finite totalization of zero UPGMA branches |
| Maximum accepted cycles | `1,000` | Published/official iterative bound |

MAFFT 5.3 predates the later nucleotide-specific threefold gap-penalty change.
Consequently DNA and RNA use the same `S^op` and `S^a` values as proteins in
this versioned profile. RNA `U` and DNA `T` represent the same canonical K80
state while remaining different typed alphabets at the public API.

## Normalized scoring matrices

For raw substitution scores `M_ab` and background frequencies `f_a`, define

```text
average1 = sum_a     f_a     M_aa
average2 = sum_a,b   f_a f_b M_ab

M_hat(a,b) = (M_ab - average2) / (average1 - average2) + S^a.
```

The expected score for identical background-distributed residues is
`1 + S^a`; the expected independent-pair score is `S^a`. Protein scores use
the published BLOSUM62 integer matrix and its associated background
frequencies. Nucleotide scores use the K80 200-PAM construction already
documented for FFT-NS-1, with transition/transversion ratio 2 and equal base
frequencies.

Ambiguity characters are expectations over their canonical states. Protein
`B`, `Z`, `J`, and `X` are equal mixtures; `U` maps to cysteine and `O` to
lysine because the 20-state 2005 matrix has no separate rows. IUPAC nucleotide
ambiguities are equal mixtures over their represented bases. These rules are
an explicit typed-alphabet extension, not a claim about undocumented MAFFT 5.3
parsing.

## All-pairs global alignment and distance

Each unordered pair is aligned as two one-row profiles using FFT correlation,
homologous-segment chaining, and the extension-free profile DP specified by
`MAFFT-2002-FFT-NS-1-v1`, but with the 2005 matrix and gap parameters. The
pair score `t_st` is independently replayed from the resulting profile path.
Self score `t_ss` is the ungapped sum of `M_hat(a,a)`.

The 2005 paper says that the pairwise scores replace shared tuple counts in the
earlier distance equation. This profile therefore defines

```text
D_st = clamp(1 - t_st / min(t_ss, t_tt), 0, 1).
```

Two empty sequences have distance zero. If exactly one sequence is empty, or
the denominator is non-positive, distinct sequences have distance one. The
clamp is required by the public finite `[0,1]` distance invariant and
totalizes negative pair scores and finite-rounding overshoot.

## Mixed-linkage guide tree

At every agglomerative step the closest two active clusters `L` and `R` are
joined. For any third cluster `C`, MAFFT's modified linkage is

```text
D_PC = s * (D_LC + D_RC) / 2
       + (1 - s) * min(D_LC, D_RC),       s = 0.1.
```

This is neither size-weighted ordinary UPGMA nor pure single linkage. Merge
height is half the selected distance. A negative child branch caused by a
non-ultrametric mixed-linkage update is stored as zero while preserving the
parent topology; Gotoh weighting subsequently applies its positive floor.

Tree ties use selected distance and then the lexicographically ordered minimum
input indices of the two clusters. Children are stored with the smaller
minimum input index first. The tree is built once and is not reconstructed
from the progressive alignment or during refinement, as explicitly reported
for [GHF]-INS-i.

## Primary consistency library

### Gap-free segments

Every pair alignment is divided into maximal consecutive runs of
residue-residue columns. Segment `n` for sequences `s` and `t` stores:

```text
P(s,t,n) = ordered original residue-index pairs in the run
L(s,t,n) = |P(s,t,n)|
S(s,t,n) = sum over (p,q) in P(s,t,n) of M_hat(s[p], t[q]).
```

Gap-opening costs separate segments but are not included in `S`, because the
paper defines `S` as the alignment score of the gap-free segment itself.

### Residue frequency

Let `w_t` be the normalized Thompson/CLUSTALW weight of sequence `t`. The
frequency with which residue `p` of sequence `s` participates in the primary
library is

```text
f(s,p) = sum over n,t,q where (s,t,p,q) in P(s,t,n) of w_t.
```

For a global pair alignment, a residue contributes at most once for any other
sequence, but residues aligned to gaps do not contribute for that pair.

### Segment importance and sparse matrix

```text
E(s,t,n) = S(s,t,n)
           * sum over (i,j) in P(s,t,n)
               [f(s,i) + f(t,j)] / [2 L(s,t,n)]

I(s,t,p,q) = sum over n containing (p,q) of E(s,t,n),
             or 0 if no segment contains the residue pair.
```

G-INS-i has one non-overlapping global alignment per pair, so a residue pair
normally belongs to at most one segment. The sum is retained in the data model
because it is the published definition and will also support later L/E-INS-i
libraries without changing the primitive.

No T-Coffee/ProbCons library extension is performed. The official MAFFT
description explicitly distinguishes its primary-only consistency score from
three-sequence library extension.

## Consistency-aware profile alignment

For profile groups `G_i` and `G_j`, columns `p` and `q` receive

```text
H(G_i,G_j,p,q) =
    sum over s in G_i, t in G_j of w_st *
        [M_hat(A(s,p), A(t,q)) + W^I I(s,t,r_s(p),r_t(q))].
```

`r_s(p)` is the original ungapped residue index represented by sequence `s` at
profile column `p`; a gap has no residue index. Progressive `w_st` factors as
the product of normalized Thompson sequence weights inside the two child
profiles. Iterative `w_st` factors into Gotoh leaf-to-cut weights.

The implementation materializes one temporary dense `H` bonus rectangle by
scattering the sparse importance entries through residue-to-column maps. Base
substitution expectations remain in the existing transformed sparse profile
score. The full score is independently replayed after traceback.

The 2005 table explicitly disables the FFT approximation in the progressive
alignment stage. Accordingly every G-INS-1 profile merge evaluates the full
scalar DP rectangle. Only the preceding all-pairs global alignments use FFT.

## Extension-free profile gap cost

For a gap consuming profile columns from boundary `x` to boundary `i`, the
inherited cost is

```text
G(x,i) = S^op * [1 - (g_start(x) + g_end(i))/2].
```

`g_start` and `g_end` are normalized weighted fractions of rows whose existing
gap starts or ends at the corresponding boundary. There is no separate affine
extension state: the random-match offset `S^a` makes a long gap comparable to
an unrelated residue run. Leading, internal, and terminal gaps all use the
same boundary rule.

For equal DP scores, diagonal is preferred, followed by a left-profile run,
then a right-profile run. Equal run prefixes retain the latest origin, which
prefers a shorter newly introduced gap. These total rules are inherited from
the independently validated scalar profile primitive.

## G-INS-1 progressive alignment

The fixed guide tree is traversed in deterministic postorder. A leaf is a
one-row ungapped profile with its Thompson weight. Each branch aligns its two
child profiles with the full consistency-aware DP described above. Within each
child, previously established residue homologies are fixed. Rows are restored
to original input order after the root merge. This initial result is G-INS-1,
not the final G-INS-i result.

## Fixed-tree iterative refinement

The degree-two root of the fixed tree is suppressed to produce its `2N-3`
unrooted edges. Deleting each edge yields one canonical bipartition. The side
with fewer leaves is stored first; equal-size sides are ordered
lexicographically. Partitions are evaluated in that canonical order.

Gotoh's three-way factorization, effective-resistance construction,
`F=1.15`, logarithmic path accumulation, and mean-one pair normalization are
identical to the separately documented FFT-NS-i primitive. Unlike FFT-NS-i,
the G-INS-i tree and its pair weights remain fixed for the entire search.

For one partition:

1. project the current MSA onto each side;
2. remove columns that are all gaps within that side;
3. preserve both projected subalignments exactly;
4. assign normalized Gotoh leaf-to-cut weights;
5. realign the two profiles with full consistency-aware DP;
6. verify all original sequences and both fixed subalignments;
7. evaluate the complete candidate with the fixed `WSP+I` objective.

All candidates in a cycle are generated from the same current alignment. The
strictly highest scoring candidate is accepted only when it improves the
current objective. One accepted candidate is one cycle.

## WSP+I objective

For an MSA `A`, fixed Gotoh pair weights `w_st`, and induced pair alignment
`A_st`, define

```text
pair_score(A_st) =
    sum of M_hat over residue-residue columns
    - S^op * number of maximal orientation-specific residue-gap runs

importance(A_st) =
    sum of I(s,t,p,q) for residue pairs (p,q) aligned in one MSA column

WSP+I(A) = sum over s<t of w_st *
           [pair_score(A_st) + W^I * importance(A_st)].
```

Gap-gap columns are absent from the induced pair alignment. A gap in `s`
followed by a gap in `t` starts a different orientation-specific run. Terminal
gap runs are penalized. The same `W^I=2.7` used in group DP multiplies the
objective importance term; this keeps candidate generation and acceptance on
one versioned scale.

Within-side induced pair alignments cannot change for a valid candidate. Their
weighted contributions are cached, while crossing-pair contributions are
recomputed in original row-major pair order. A direct evaluator verifies that
the exact cache and full objective agree.

The 2005 executable anchored highly conserved regions during refinement to
reduce CPU time. Its complete anchoring threshold and boundary rules are not
published. This scalar profile therefore evaluates the full profile rectangle
and anchors no conserved region. The omission removes a performance
restriction; it does not remove pairwise evidence, change `H`, change the
fixed subalignments, or weaken candidate acceptance.

## Stopping and determinism

The best candidate is accepted only under ordinary finite IEEE-754

```text
candidate_score > current_score.
```

There is no epsilon. Equal candidate scores retain the first canonical
partition and an equal score never replaces the current alignment. Refinement
stops after a full non-improving cycle or after 1,000 accepted cycles.

Input order determines all remaining ties. No random number, hash iteration
order, locale, clock, thread schedule, external executable, or filesystem state
can affect a library result.

## Pseudocode

```text
G_INS_I(records):
    validate non-empty record collection and unique identifiers
    scorer <- normalized BLOSUM62 or K80 with Sop=1.53, Sa=0.123

    raw_pairs <- []
    distance_values <- zero N-by-N matrix
    for each unordered pair (s,t) in input order:
        pair_alignment <- FFT_GLOBAL_ALIGN(records[s], records[t], scorer)
        verify pair score by replay
        segments <- MAXIMAL_GAP_FREE_SEGMENTS(pair_alignment)
        t_st <- pair alignment score
        D_st <- TOTALIZED_NORMALIZED_SCORE_DISTANCE(t_st, t_ss, t_tt)
        store raw pair, segments, score, and distance

    tree <- MIXED_LINKAGE_TREE(distance_values, s=0.1)
    thompson_weights <- CLUSTALW_PATH_WEIGHTS(tree)

    frequency <- zero value for every original residue
    for each raw pair segment and residue pair (p,q):
        frequency[s,p] += thompson_weights[t]
        frequency[t,q] += thompson_weights[s]

    primary_library <- empty sparse pair library
    for each raw segment:
        E <- segment score times its mean endpoint frequency
        add E for every residue pair in the segment

    initial <- PROGRESSIVE_POSTORDER(tree):
        leaf: one-row profile with Thompson weight
        branch(left,right):
            H <- BASE_PROFILE_SCORE + 2.7 * SCATTER(primary_library)
            FULL_PROFILE_DP(left, right, H, Sop=1.53)
    current <- initial

    gotoh <- FIXED_GOTOH_WEIGHTING(tree, F=1.15)
    objective_cache <- EXACT_WSP_PLUS_I_CACHE(current, gotoh, primary_library)

    repeat up to 1000 accepted cycles:
        before <- objective_cache.total
        best <- NONE
        for split in CANONICAL_TREE_EDGE_SPLITS(gotoh):
            left, right <- PROJECT_FIXED_SUBALIGNMENTS(current, split)
            H <- BASE_PROFILE_SCORE + 2.7 * SCATTER(primary_library)
            candidate <- FULL_PROFILE_DP(left, right, H, Sop=1.53)
            validate source sequences and fixed subalignments
            score <- EXACT_CACHED_WSP_PLUS_I(candidate, split)
            if score is strictly greater than best.score:
                best <- candidate

        if best is NONE or best.score <= before:
            record convergence and stop

        current <- best.alignment
        objective_cache <- exact contributions for current
        record accepted cycle

    verify final objective directly
    return alignment, pair evidence, library statistics, fixed tree,
           progressive merges, refinement evidence, and memory diagnostics
```

## Data structures

- a dense symmetric `[0,1]` pair-score distance matrix;
- one raw pair record at library-construction time, containing its global
  score and maximal gap-free segments;
- sparse ordered residue-pair importance entries for every unordered sequence
  pair;
- per-residue temporary Thompson frequency arrays;
- one fixed mixed-linkage `GuideTree`;
- normalized Thompson weights for progressive alignment;
- an unrooted fixed tree, canonical partitions, Gotoh path factors, and a dense
  symmetric pair-weight matrix for refinement;
- typed `WorkingProfile` rows whose optional cells preserve original residue
  order and original input index;
- transformed sparse base profile scores plus one temporary dense consistency
  bonus rectangle;
- rolling DP score rows and packed traceback cells;
- an exact per-pair `WSP+I` contribution cache;
- auditable pair, progressive-merge, cycle, termination, and memory summaries;
- the final typed `MultipleAlignment` in original row order.

## Complexity

Let `N` be the sequence count, `L` the maximum ungapped length, `W` a current
MSA width, `C` the number of accepted cycles, and `P(a,b)` the FFT-assisted
pair-alignment cost for lengths `a` and `b`.

- all-pairs stage: `sum_(s<t) P(L_s,L_t)`, worst-case `O(N^2 L^2)`;
- raw segment and primary-library construction: `O(N^2 L)` time and space in
  the normal one-global-alignment-per-pair case;
- scalar mixed-linkage tree: `O(N^3)` time and `O(N^2)` space;
- progressive alignment: `N-1` full profile rectangles, worst-case
  `O(N W^2)` time with `O(W^2 + NW)` working space;
- fixed Gotoh factors and pair weights: `O(N^2)` clear scalar work and space;
- one candidate profile realignment: `O(W^2 + K_cross)` to scatter its
  crossing-library entries, plus candidate validation;
- candidate objective scoring: proportional to crossing induced-pair work,
  with same-side contributions reused exactly;
- one attempted cycle evaluates `2N-3` candidates, with worst-case
  `O(N^3 W + N W^2)` scalar time;
- total refinement is the attempted-cycle cost times at most `C+1`;
- returned alignment and evidence: `O(NW + CN + N^2)` excluding the retained
  `O(N^2 L)` primary library.

These are conservative bounds. The correctness baseline is intentionally
scalar. Sparse scattering, contribution reuse, rolling score rows, and packed
traceback are clarity-preserving optimizations; SIMD and parallel candidate
evaluation may be added only after the scalar profile satisfies every gate.

## Published omissions and explicit resolutions

| Omission or ambiguity | `v1` resolution |
| --- | --- |
| Exact historical executable | Reproduce MAFFT 5.3 from the 2005 publications; MAFFT 7.526 is black-box evidence only. |
| Pair score distance details | Earlier MAFFT equation with pair alignment scores, clamped to `[0,1]`; empty/non-positive denominator is explicitly totalized. |
| Modified UPGMA linkage | Published mixed linkage with `s=0.1`, not ordinary size-weighted UPGMA. |
| Protein background frequencies | Frequencies associated with the published BLOSUM62 profile. |
| DNA/RNA gap change in later MAFFT | Not applied; v5.3 uses `1.53/0.123` for all supported molecule types. |
| Segment score | Sum of normalized substitution scores inside the maximal gap-free run; surrounding gap cost excluded. |
| Progressive pair-weight factorization | Product of normalized Thompson leaf weights within the two profiles. |
| Iterative equalization exponent | `F=1.15`, the explicit versioned resolution shared with FFT-NS-i. |
| Guide-tree reconstruction | None after the all-pairs tree, matching the 2005 report. |
| Importance multiplier in final objective | Apply the same published `W^I=2.7` as the profile matrix. |
| Conserved-region anchors | Omitted because thresholds/boundaries are unpublished; evaluate the full rectangle. |
| Library extension | None; primary pairwise library only. |
| Gap-gap and terminal scoring | Delete pairwise gap-gap columns; penalize every maximal orientation-specific gap run, including terminal runs. |
| Candidate search | Best-first over every canonical fixed-tree edge split. |
| Candidate acceptance | Strict finite `>` with no epsilon. |
| Equal candidates | First canonical partition. |
| Singleton | Return the one-row G-INS-1 alignment and a converged zero-partition refinement record. |

Every resolution is part of the profile identifier. Changing one requires a
new profile version and independent evidence.

## Validation and acceptance criteria

G-INS-i can be marked complete only when all of the following hold:

1. this document contains references, equations, original pseudocode, data
   structures, complexity, parameters, ties, omissions, and limitations;
2. hand-computed tests verify BLOSUM62/K80 normalization and ambiguity
   expectations;
3. every pair path and score is independently replayed and every recovered row
   equals its original typed sequence;
4. hand-computed distance and mixed-linkage cases verify the `s=0.1` equation,
   ties, clamping, empty sequences, and topology;
5. hand-computed pair alignments verify maximal gap-free segmentation,
   `f`, `E`, sparse `I`, symmetry, and absence of library extension;
6. hand-computed profile cells verify base-plus-importance scattering with
   Thompson and Gotoh weights;
7. the progressive stage uses no FFT partition and preserves every child
   subalignment;
8. direct and cached `WSP+I` evaluators agree exactly on matches,
   substitutions, leading/internal/terminal gaps, opposite gap orientations,
   and gap-gap columns;
9. every refinement candidate preserves both fixed subalignments and every
   original sequence, and every accepted cycle strictly improves one fixed
   objective;
10. tests cover protein, DNA, RNA, empty, singleton, identical, divergent,
    repetitive, ambiguity-rich, leading-gap, internal-gap, and terminal-gap
    inputs;
11. property tests independently verify rectangularity, exact source
    consumption, finite evidence, objective monotonicity, and determinism;
12. inert fixtures from official MAFFT 7.526 verify exact trivial cases and
    SP/TC agreement on nontrivial cases without entering production code;
13. a fixed external BAliBASE subset records SP and TC without using this
    implementation to generate expected results;
14. Criterion benchmarks record complete-pipeline sequence-length,
    sequence-count, protein, DNA, and RNA costs, while checked memory evidence
    separately accounts for pair-library, progressive, and refinement phases;
15. matched fresh-process runs compare the Rust worker with opaque MAFFT 7.526
    using mirrored inputs, correctness checks, warm-up, sampling, binary and
    harness hashes, raw intervals, time, and peak memory;
16. formatting, warnings-denied Clippy, workspace/release tests, Rustdoc, and
    transitive dependency license/source/advisory checks pass;
17. the public API and examples are safe, deterministic, English-only, and
    contain no unfinished route, placeholder, mock, or runtime reference
    dependency.

## Recorded validation evidence

The scalar `v1` profile completed its declared validation laboratory on
2026-08-05:

- the `align-mafft` suite passes 104 unit, integration, reference, property,
  and documentation tests; the property tests exercise 96 generated cases in
  addition to their fixed cases;
- five inert MAFFT 7.526 cases produce mean SP `0.974657120` and mean TC
  `0.921038961`; the two unambiguous identical cases are exactly equal;
- the six-family BAliBASE 3 subset produces macro SP/TC
  `0.902059636/0.703589360` and micro SP/TC
  `0.888139574/0.828063241`;
- Criterion covers length scaling, sequence-count scaling, protein, DNA, and
  RNA, and `GinsiMemoryEstimate` uses checked phase-specific arithmetic;
- the reproducible Windows process laboratory records valid output, SP/TC,
  warm-ups, nine samples, raw ranges, artifact hashes, modeled memory, and
  sampled process working sets;
- debug and release workspace tests, warnings-denied Clippy, formatting,
  Rustdoc, and `cargo-deny check` pass; the dependency audit reports
  advisories, bans, licenses, and sources as `ok`.

Full timing, BAliBASE, memory, artifact identity, and reproduction details are
in `docs/benchmarks/mafft-g-ins-i.md`.

Within the explicitly versioned resolutions and limitations in this document,
`MAFFT-2005-G-INS-i-v1` is complete as a deterministic scalar scientific
profile. This status does not apply to modern MAFFT behavior or to later
G-INS-i optimization profiles.

## Black-box reference boundary

The isolated comparator is the official signed Windows MAFFT 7.526 package
(`v7.526`, 2024-04-26), SHA-256
`7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7`.
The pinned operational command is

```text
mafft.bat --globalpair --maxiterate 1000 --thread 1 --inputorder --quiet
          --out OUTPUT INPUT
```

The package is opaque: its source is neither downloaded nor inspected. The
signed archive, executable, inputs, outputs, and orchestration stay in the
sibling validation laboratory and never enter the library at runtime. Exact
alignment equality is required only for unambiguous trivial cases. Nontrivial
comparison uses valid source recovery, determinism, SP, TC, score diagnostics,
runtime, and scoped memory observations because the modern executable is
documented to differ from this ancestral profile.

## Known limitations of this profile

- It reconstructs the published MAFFT 5.3 strategy, not all current MAFFT
  behavior or an unavailable historical binary.
- The publications do not completely define empty-sequence distance,
  BLOSUM ambiguity rows, every DP tie, negative mixed-linkage branch handling,
  or exact terminal-gap objective semantics; each is versioned above.
- The v5.3 nucleotide parameters predate later RNA-specific tuning and should
  not be represented as modern MAFFT defaults.
- Conserved-region anchoring is intentionally absent because its complete
  historical rule is unpublished; refinement therefore uses more CPU and may
  explore candidates the historical heuristic excluded.
- The all-pairs primary library requires quadratic pair count and normally
  `O(N^2 L)` retained evidence, so this profile is intended for the small and
  medium globally homologous sets described by the authors, not HomFam-scale
  collections.
- Full progressive and refinement rectangles remain quadratic in profile
  width in the scalar baseline.
- BLOSUM62 is appropriate for the protein profile studied by the authors, but
  no single substitution model is optimal for every evolutionary distance or
  compositional bias.
- G-INS-i assumes the complete input region is globally alignable. Sequences
  with unrelated flanks, multiple isolated motifs, rearrangements, or major
  domain-architecture differences require later L-INS-i or E-INS-i profiles.
