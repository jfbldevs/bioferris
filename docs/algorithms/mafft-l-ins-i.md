# MAFFT L-INS-i

## Reproduction target

This document specifies the immutable scientific profile
`MAFFT-2005-L-INS-i-v1`.

The target is the ancestral L-INS-i strategy reported by Katoh, Kuma, Miyata,
and Toh in *Genome Informatics* in 2005. It is a local-consistency multiple
sequence aligner: every sequence pair is aligned with a local affine-gap
Smith-Waterman algorithm, only the locally aligned gap-free segments enter the
primary consistency library, the library contributes to every progressive
profile alignment, and the resulting L-INS-1 alignment is refined against one
fixed `WSP+I` objective.

This target must not be confused with H-INS-i or F-INS-i in the earlier MAFFT
5.3 NAR paper. Those experimental strategies delegated local pair alignment to
FASTA and are therefore outside this library's no-wrapper and no-external-engine
boundary. The later authors' publication names L-INS-i directly, cites
Smith-Waterman for its all-pairs stage, and the authors' subsequent review
identifies an affine gap cost. This profile reconstructs that algorithm from
the publications and documented parameters without using FASTA, MAFFT source,
translated code, FFI, subprocesses, or runtime executables.

The profile is deliberately versioned. Current MAFFT has accumulated parameter
and engineering changes, so this is not a claim of byte-for-byte equivalence
with an unavailable historical executable or MAFFT 7.526. MAFFT 7.526 is used
only as an opaque comparator in the isolated validation laboratory.

## Primary references

1. Katoh K, Kuma K, Miyata T, Toh H. *Improvement in the Accuracy of Multiple
   Sequence Alignment Program MAFFT*. Genome Informatics. 2005;16(1):22-33.
   DOI: [10.11234/gi1990.16.22](https://doi.org/10.11234/gi1990.16.22).
   This is the primary L-INS-i publication. It specifies the [GL]-INS-i
   pipeline, local all-pairs stage, primary gap-free-segment library, importance
   equations, fixed guide tree, progressive alignment, and `WSP+I` refinement.
2. Katoh K, Kuma K, Toh H, Miyata T. *MAFFT version 5: improvement in accuracy
   of multiple sequence alignment*. Nucleic Acids Research.
   2005;33(2):511-518.
   DOI: [10.1093/nar/gki198](https://doi.org/10.1093/nar/gki198).
   This defines the immediate ancestral consistency framework and explains the
   distinction between global and local pair evidence. Its H/F-INS-i external
   FASTA routes are historical context, not implementation sources.
3. Smith TF, Waterman MS. *Identification of common molecular subsequences*.
   Journal of Molecular Biology. 1981;147(1):195-197.
   DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).
   This defines the zero-restarted local dynamic-programming principle.
4. Gotoh O. *An improved algorithm for matching biological sequences*.
   Journal of Molecular Biology. 1982;162(3):705-708.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
   This supplies the affine-gap state decomposition used in the local
   pairwise stage.
5. Katoh K, Toh H. *Recent developments in the MAFFT multiple sequence
   alignment program*. Briefings in Bioinformatics. 2008;9(4):286-298.
   DOI: [10.1093/bib/bbn013](https://doi.org/10.1093/bib/bbn013).
   This authors' review explicitly distinguishes L-INS-i as local pair
   alignment with an affine gap cost and explains its terminal-flank use case.
6. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic
   Acids Research. 2002;30(14):3059-3066.
   DOI: [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
   This defines the normalized substitution matrices, pair-score distance,
   modified linkage, position-specific extension-free group gap cost, and
   profile DP inherited by [GL]-INS-i.
7. Thompson JD, Higgins DG, Gibson TJ. *CLUSTAL W: improving the sensitivity
   of progressive multiple sequence alignment through sequence weighting,
   position-specific gap penalties and weight matrix choice*. Nucleic Acids
   Research. 1994;22(22):4673-4680.
   DOI: [10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
   This supplies the progressive-stage tree weights cited by Katoh et al.
8. Gotoh O. *A weighting system and algorithm for aligning many
   phylogenetically related sequences*. Computer Applications in the
   Biosciences. 1995;11(5):543-551.
   DOI: [10.1093/bioinformatics/11.5.543](https://doi.org/10.1093/bioinformatics/11.5.543).
   This supplies iterative pair weights and the WSP objective.
9. Hirosawa M, Totoki Y, Hoshida M, Ishikawa M. *Comprehensive study on
   iterative algorithms of multiple sequence alignment*. Computer
   Applications in the Biosciences. 1995;11(1):13-18.
   DOI: [10.1093/bioinformatics/11.1.13](https://doi.org/10.1093/bioinformatics/11.1.13).
   This defines tree-dependent restricted partitioning.
10. Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein
    blocks*. Proceedings of the National Academy of Sciences.
    1992;89(22):10915-10919.
    DOI: [10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915).
11. Kimura M. *A simple method for estimating evolutionary rates of base
    substitutions through comparative studies of nucleotide sequences*.
    Journal of Molecular Evolution. 1980;16:111-120.
    DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581).

Official operational evidence:

12. The official algorithm page describes L-INS-i as suitable for one
    alignable domain surrounded by unrelated terminal flanks and gives the
    modern opaque command `mafft --localpair --maxiterate 1000 INPUT`:
    <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.
13. The official MAFFT manual records `--lop=-2.00`, `--lep=0.1`,
    `--lexp=-0.1`, `--op=1.53`, `--ep=0.123`, and `--weighti=2.7`:
    <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.

Items 12-13 establish documented parameters, intended use, and the opaque
reference boundary. They do not replace the 2005 scientific target with
current executable behavior.

## Objective and characteristic pipeline

L-INS-i targets protein, DNA, or RNA sequences with one homologous domain and
possibly unrelated N- and C-terminal extensions. Its characteristic pipeline
is:

1. compute every unordered pair with scalar local affine-gap
   Smith-Waterman/Gotoh DP;
2. derive a normalized score distance from each maximum local score;
3. build one mixed-linkage UPGMA-like guide tree;
4. compute Thompson/CLUSTALW sequence weights from that tree;
5. divide each local traceback into maximal gap-free residue-pair segments;
6. derive residue frequencies, segment importance values, and the sparse
   primary importance matrix `I` from those local segments only;
7. progressively align complete input profiles along the fixed tree with full
   scalar profile DP and the combined substitution-plus-importance matrix,
   producing L-INS-1;
8. derive Gotoh pair weights and all tree-edge partitions from the same fixed
   tree;
9. repeatedly realign the two fixed subalignments of every partition and
   accept only the strictly best improvement in `WSP+I`;
10. stop when a complete cycle has no strict improvement or after 1,000
    accepted cycles.

The local stage does not truncate the final MSA. Every input residue is present
exactly once in L-INS-1 and the final L-INS-i alignment. Locality controls the
distance and consistency evidence: residues outside the selected local domain
do not enter `P`, `S`, `f`, `E`, or `I`.

Replacing the local affine all-pairs stage with global alignment, tuple counts,
linear gaps, multiple disjoint local hits, generalized-affine gaps, external
FASTA, or an executable call would be a different algorithm. Rebuilding the
guide tree, using only WSP, performing T-Coffee library extension, or returning
L-INS-1 without refinement would also be different algorithms.

## Versioned parameters

| Parameter | Value | Scientific role |
| --- | ---: | --- |
| Protein matrix | BLOSUM62 | Matrix selected in the 2005 study |
| Nucleotide matrix | K80, transition/transversion ratio 2 | Published MAFFT nucleotide model |
| Group `S^op` | `1.53` | Extension-free progressive/refinement gap-opening cost |
| Group `S^a` | `0.123` | Progressive/refinement normalized-matrix offset |
| Local `lop` | `-2.0` | Signed affine gap-open contribution |
| Local `lep` | `0.1` | Positive offset penalty subtracted from the local normalized matrix |
| Local `lexp` | `-0.1` | Signed affine per-extension contribution |
| `W^I` | `2.7` | Primary-library importance multiplier |
| Mixed-linkage coefficient `s` | `0.1` | MAFFT UPGMA-like linkage |
| Gotoh equalization exponent | `1.15` | Shared fixed-tree iterative weighting resolution |
| Zero branch floor | `10^-12` | Finite totalization of zero tree branches |
| Maximum accepted cycles | `1,000` | Published/official iterative bound |

MAFFT 5 predates the later threefold nucleotide pair-gap change. Therefore
protein, DNA, and RNA use the same signed local gap contributions in this
ancestral profile. RNA `U` and DNA `T` represent the same canonical K80 state
while remaining different typed alphabets in the public API.

## Normalized scoring matrices

For raw substitution scores `M_ab` and background frequencies `pi_a`, define

```text
average1 = sum_a     pi_a      M_aa
average2 = sum_a,b   pi_a pi_b M_ab

M_hat_delta(a,b) = (M_ab - average2) / (average1 - average2) + delta.
```

The complete-profile stages use `delta=0.123`. The local all-pairs stage uses
the documented positive `lep=0.1` as a penalty and therefore sets
`delta=-0.1`. Protein scores use the published BLOSUM62 integer matrix and
associated background frequencies. Nucleotide scores use the K80 200-PAM
construction with transition/transversion ratio 2 and equal canonical
frequencies.

The 2005 paper names Smith-Waterman but does not print its local substitution
recurrence, while the official manual calls `lep` an offset and publishes a
positive command-line value. Its sign was therefore resolved experimentally
at the allowed opaque-reference boundary. A sweep of MAFFT 7.526's public
`--lep` option over `-1.0..1.0`, holding every other option fixed, showed that
increasing `lep` progressively removed adversarial flank homologies and
increased the local-score guide-tree distances. Negative values caused the
same flanks to be over-aligned. This monotonic behavior is consistent with
`M_hat_-lep` and contradicts `M_hat_+lep`. The probe recorded only final FASTA
and documented `--treeout` output; source and undocumented intermediate files
were not inspected.

Ambiguity symbols are expectations over canonical states. Protein `B`, `Z`,
`J`, and `X` are equal mixtures; `U` maps to cysteine and `O` to lysine because
the historical 20-state matrix has no separate rows. IUPAC nucleotide
ambiguities are equal mixtures over represented bases. This is an explicit
typed-alphabet extension, not a claim about undocumented MAFFT 5 parsing.

## Local affine pair alignment

Let `X[i,j]` be the best positive local score ending with residue `i-1` from
the first sequence aligned to a gap, `Y[i,j]` the symmetric state, and
`H[i,j]` the best local score ending at boundary `(i,j)`. With signed
`g_o=-2.0` and `g_e=-0.1`, this profile defines:

```text
X[i,j] = max(H[i-1,j] + g_o, X[i-1,j] + g_e)
Y[i,j] = max(H[i,j-1] + g_o, Y[i,j-1] + g_e)

H[i,j] = max(
    0,
    H[i-1,j-1] + M_hat_-0.1(a[i-1], b[j-1]),
    X[i,j],
    Y[i,j]
).
```

Boundary `H` cells are zero and boundary gap states are unreachable. A
one-residue gap contributes `g_o`; each additional residue contributes `g_e`.
This open-only first-residue convention is an explicit resolution because the
official manual publishes the three parameter values but not the recurrence.
Direct transitions between opposite gap states are not allowed.

The selected endpoint is the first strict maximum in row-major `(i,j)` order.
Within `H`, positive ties prefer diagonal, then `X`, then `Y`. Within a gap
state, a tie prefers opening from `H` over extension, making the newest gap
shorter. A zero score always selects restart. Traceback stops at the selected
zero boundary and reports half-open source ranges plus forward operations.
When no positive score exists, the canonical result is score zero, ranges
`0..0`, and no operations.

The local score is independently replayed from the reported source ranges and
operations. Replay verifies substitution scores, affine gap runs, exact range
consumption, positive local prefixes, the terminal score, and source bounds;
it does not read the DP buffers.

## Pair-score distance

For pair score `t_st` and ungapped self scores calculated with the local
matrix,

```text
t_ss = sum_p M_hat_-0.1(s[p], s[p])
t_tt = sum_q M_hat_-0.1(t[q], t[q])

D_st = clamp(1 - t_st / min(t_ss, t_tt), 0, 1).
```

Two empty sequences have distance zero. If exactly one sequence is empty, the
normalization denominator is non-positive, or no positive local hit exists,
distinct sequences have distance one. The clamp totalizes negative,
non-positive, and finite-rounding edge cases while preserving the public
finite `[0,1]` distance invariant.

Because the denominator covers the complete shorter sequence while `t_st`
may cover only one local domain, long unrelated extensions can increase this
distance. The publications do not define an alternative local-only
normalization, so the inherited MAFFT pair-score equation is retained and the
choice is versioned.

## Primary local consistency library

### Gap-free segments

Only the selected maximum local traceback is used. It is divided into maximal
consecutive residue-residue runs. Segment `n` for sequences `s` and `t` stores

```text
P(s,t,n) = ordered original residue-index pairs in the run
L(s,t,n) = |P(s,t,n)|
S(s,t,n) = sum over (p,q) in P(s,t,n)
               M_hat_-0.1(s[p], t[q]).
```

Internal affine gap columns separate segments. Unaligned terminal prefixes
and suffixes do not form segments. Gap contributions are excluded from `S`
because `S` is the score of the gap-free segment. Original residue indices,
not local-slice offsets, are stored.

### Residue frequency, segment importance, and sparse matrix

Let `w_t` be the normalized Thompson/CLUSTALW weight of sequence `t`:

```text
f(s,p) = sum over n,t,q where (s,t,p,q) in P(s,t,n) of w_t.

E(s,t,n) = S(s,t,n)
           * sum over (i,j) in P(s,t,n)
               [f(s,i) + f(t,j)] / [2 L(s,t,n)].

I(s,t,p,q) = sum over n containing (p,q) of E(s,t,n),
             or 0 otherwise.
```

The data model retains the published sum even though this profile selects one
non-overlapping local traceback per pair and a residue pair normally occurs in
at most one segment. No T-Coffee/ProbCons three-sequence library extension is
performed.

## Mixed-linkage guide tree

At each agglomerative step the closest active clusters `L` and `R` are joined.
For a third cluster `C`, the fixed update is

```text
D_PC = s * (D_LC + D_RC) / 2
       + (1 - s) * min(D_LC, D_RC),       s = 0.1.
```

Merge height is half the selected distance. A negative child branch caused by
a non-ultrametric update is stored as zero; Gotoh weighting later applies its
positive floor. Distance ties use the lexicographically ordered minimum input
indices of the two clusters, and children store the smaller minimum index
first. The tree is built once and never reconstructed.

## Consistency-aware complete-profile alignment

For complete profile groups `G_i` and `G_j`, columns `p` and `q` receive

```text
H(G_i,G_j,p,q) =
    sum over s in G_i, t in G_j of w_st *
      [M_hat_0.123(A(s,p), A(t,q))
       + W^I I(s,t,r_s(p),r_t(q))].
```

`r_s(p)` is the original residue index at profile column `p`; a gap has no
index. Progressive `w_st` factors from normalized Thompson leaf weights.
Iterative `w_st` factors from Gotoh leaf-to-cut weights. The implementation
scatters the sparse local importance entries through residue-to-column maps
into one temporary dense bonus rectangle; base substitution expectations stay
in the transformed sparse profile score.

Every progressive and refinement merge evaluates the full scalar DP
rectangle. The 2005 executable's complete conserved-anchor rule is not
published, so this profile uses no anchor restriction. FFT is not used in the
local pair, progressive, or refinement stages.

For a gap consuming profile columns from boundary `x` to `i`, the inherited
complete-profile cost is extension-free:

```text
G(x,i) = S^op * [1 - (g_start(x) + g_end(i))/2],  S^op=1.53.
```

`g_start` and `g_end` are normalized weighted fractions of rows whose existing
gap starts or ends at the boundary. Leading, internal, and terminal gaps use
the same rule. Equal complete-profile scores prefer diagonal, a left-profile
gap run, then a right-profile gap run; equal run prefixes retain the latest
origin and therefore the shorter newly introduced gap.

## L-INS-1 and fixed-tree iterative refinement

The fixed guide tree is traversed in deterministic postorder. Leaves are full
ungapped input sequences, including locally unaligned flanks. Each branch
aligns its two child profiles with the complete consistency-aware DP above and
preserves all established within-child homologies. Rows return to original
input order at the root. This result is L-INS-1.

The degree-two root is suppressed to produce `2N-3` unrooted edges. Deleting
each edge yields one canonical bipartition; the smaller side is stored first,
with lexical ordering for equal sizes. For each partition:

1. project the current MSA onto both sides;
2. remove columns that are all gaps within a side;
3. preserve both projected subalignments exactly;
4. assign normalized Gotoh leaf-to-cut weights with equalization exponent
   `1.15`;
5. globally realign the two complete profiles using local-library bonuses;
6. verify every original sequence and both fixed subalignments;
7. evaluate the complete candidate against the fixed `WSP+I` objective.

All candidates in a cycle start from the same current MSA. Only the strictly
highest improving candidate is accepted, so one accepted candidate is one
cycle.

## WSP+I objective

For complete MSA `A`, fixed Gotoh pair weights `w_st`, and its induced pair
alignment `A_st`, define

```text
pair_score(A_st) =
    sum of M_hat_0.123 over residue-residue columns
    - 1.53 * number of maximal orientation-specific residue-gap runs

importance(A_st) =
    sum of I(s,t,p,q) for residue pairs aligned in one MSA column

WSP+I(A) = sum over s<t of w_st *
           [pair_score(A_st) + 2.7 * importance(A_st)].
```

Gap-gap columns are removed from induced pair alignments. A gap in `s`
followed by a gap in `t` starts a different orientation-specific run, and
terminal gap runs are penalized. Same-side pair contributions are cached while
crossing contributions are recomputed in original row-major pair order. An
independent direct evaluator must agree with the cache before returning.

Candidate acceptance uses ordinary finite IEEE-754 strict comparison:

```text
candidate_score > current_score.
```

There is no epsilon. Equal candidates retain the first canonical partition;
an equal score never replaces the current alignment. Refinement stops after a
complete non-improving attempt or 1,000 accepted cycles.

## Pseudocode

```text
L_INS_I(records):
    validate non-empty records and unique identifiers
    local_scorer <- normalized BLOSUM62 or K80 with signed offset -0.1
    group_scorer <- normalized BLOSUM62 or K80 with offset 0.123

    raw_pairs <- []
    distance_values <- zero N-by-N matrix
    for each unordered pair (s,t) in input order:
        local <- LOCAL_AFFINE_SMITH_WATERMAN(
            records[s], records[t], local_scorer,
            gap_open=-2.0, gap_extend=-0.1)
        INDEPENDENTLY_REPLAY_LOCAL(local)
        segments <- MAXIMAL_LOCAL_GAP_FREE_SEGMENTS(local)
        D_st <- TOTALIZED_NORMALIZED_SCORE_DISTANCE(
            local.score, self_score[s], self_score[t])
        store ranges, score, segments, distance, and memory evidence

    tree <- MIXED_LINKAGE_TREE(distance_values, s=0.1)
    thompson_weights <- CLUSTALW_PATH_WEIGHTS(tree)

    frequency <- zero value for every original residue
    for each local segment and residue pair (p,q):
        frequency[s,p] += thompson_weights[t]
        frequency[t,q] += thompson_weights[s]

    primary_library <- empty sparse pair library
    for each local segment:
        E <- segment score times its mean endpoint frequency
        add E for every original residue pair in the segment

    initial <- PROGRESSIVE_POSTORDER(tree):
        leaf: one full ungapped record with Thompson weight
        branch(left,right):
            H <- GROUP_BASE_SCORE + 2.7 * SCATTER(primary_library)
            FULL_PROFILE_DP(left, right, H, group_gap_open=1.53)
    current <- initial

    gotoh <- FIXED_GOTOH_WEIGHTING(tree, F=1.15)
    objective_cache <- EXACT_WSP_PLUS_I_CACHE(current, gotoh, library)

    repeat up to 1000 accepted cycles:
        before <- objective_cache.total
        best <- NONE
        for split in CANONICAL_TREE_EDGE_SPLITS(gotoh):
            left, right <- PROJECT_FIXED_SUBALIGNMENTS(current, split)
            H <- GROUP_BASE_SCORE + 2.7 * SCATTER(primary_library)
            candidate <- FULL_PROFILE_DP(left, right, H, 1.53)
            verify sources and both fixed subalignments
            score <- EXACT_CACHED_WSP_PLUS_I(candidate, split)
            if score is strictly greater than best.score:
                best <- candidate

        if best is NONE or best.score <= before:
            record convergence and stop

        current <- best.alignment
        objective_cache <- exact contributions for current
        record accepted cycle

    verify the final objective directly
    return the full alignment, local ranges and pair evidence, sparse library,
           fixed tree, progressive merges, refinement evidence, and memory
```

## Data structures

- rolling `H`, `X`, and `Y` score rows plus packed local traceback states;
- one temporary local operation path and explicit source ranges;
- a dense symmetric finite pair-score distance matrix;
- raw local gap-free segments with original residue-index pairs;
- sparse ordered importance entries for every unordered sequence pair;
- per-residue Thompson frequency arrays;
- one fixed mixed-linkage guide tree and normalized Thompson weights;
- an unrooted fixed tree, canonical partitions, Gotoh path factors, and dense
  pair weights;
- typed complete `WorkingProfile` rows that preserve original residue order;
- transformed sparse base profile scores plus a dense consistency bonus
  rectangle;
- rolling complete-profile score rows and packed traceback cells;
- an exact per-pair `WSP+I` contribution cache;
- auditable pair, merge, cycle, termination, and memory summaries;
- the final typed `MultipleAlignment` in original input order.

## Complexity

Let `N` be the sequence count, `L` the maximum ungapped length, `W` a current
MSA width, and `C` the number of accepted refinement cycles.

- all-pairs local affine stage: `O(N^2 L^2)` time and `O(L^2)` peak traceback
  space per sequential pair, plus `O(L)` rolling scores;
- raw local segments and retained primary library: `O(N^2 L)` normal-case
  time and space for one local traceback per pair;
- scalar mixed-linkage tree: `O(N^3)` time and `O(N^2)` space;
- progressive complete-profile stage: `N-1` full rectangles, worst-case
  `O(N W^2)` time and `O(W^2 + NW)` working space;
- fixed Gotoh factors and pair weights: `O(N^2)` clear scalar work and space;
- one refinement candidate: `O(W^2 + K_cross)` profile work plus validation
  and crossing induced-pair scoring;
- one attempted cycle: `2N-3` candidates, conservatively
  `O(N^3 W + N W^2)` total work; independent candidates can be divided among
  `P <= 2N-3` workers without changing their canonical reduction order;
- total refinement: attempted-cycle cost times at most `C+1`;
- returned alignment and evidence: `O(NW + CN + N^2)` excluding the retained
  `O(N^2 L)` primary library.

The local and per-candidate DP kernels remain scalar. Fixed-tree refinement
candidates can be evaluated concurrently through `with_thread_count` or the
facade's `.threads(count)`. Each candidate keeps its original partition index and
the reduction applies the scalar score/tie rule, so the scientific profile and
bit-exact output do not depend on worker count. SIMD, parallel pair evaluation,
banding, and anchor restrictions remain future optimizations and may not alter
this profile's results.

## Published omissions and explicit resolutions

| Omission or ambiguity | `v1` resolution |
| --- | --- |
| H/F-INS-i versus L-INS-i | Target the 2005 Genome Informatics L-INS-i plus the authors' affine clarification; use no FASTA or external engine. |
| Exact historical executable | Reproduce the publications; MAFFT 7.526 is opaque evidence only. |
| Local affine recurrence | Signed `lop=-2.0`, `lexp=-0.1`; a one-residue gap pays open only. |
| Local offset semantics | Treat public `lep=0.1` as a positive penalty and use signed matrix offset `-0.1`; this is supported by an opaque public-option sweep. |
| Local endpoint ties | First strict row-major maximum. |
| Local state ties | Diagonal, first-sequence gap, second-sequence gap; gap-open beats equal extension. |
| Zero-score local result | Canonical empty ranges `0..0`, no operations, no library segment. |
| Multiple local hits | Retain one deterministic maximum traceback, as Smith-Waterman returns one optimum; no Waterman-Eggert enumeration. |
| Segment score | Sum local-offset substitution scores in the maximal gap-free run; exclude gap cost. |
| Local flank handling | Flanks are absent from pair evidence but remain exactly once in every complete MSA. |
| Pair score distance | Inherited normalized score equation using complete local self scores, clamped to `[0,1]`. |
| DNA/RNA historical parameters | Use the ancestral common local and group parameters; do not apply later threefold nucleotide tuning. |
| Modified UPGMA | Mixed linkage with `s=0.1`. |
| Progressive pair weights | Product of normalized Thompson leaf weights within child profiles. |
| Iterative equalization | Fixed `F=1.15`, shared with the versioned FFT-NS-i/G-INS-i resolution. |
| Guide-tree reconstruction | None after the all-pairs tree. |
| Importance multiplier | Apply `W^I=2.7` in both profile DP and final objective. |
| Conserved anchors | Omit unpublished threshold/boundary rules and evaluate full rectangles. |
| Library extension | None; primary local pair library only. |
| Complete-MSA terminal gaps | Penalize maximal orientation-specific runs, including terminal runs. |
| Candidate selection | Strict best-first over every canonical fixed-tree edge split. |
| Singleton | Return the one-row L-INS-1 alignment and a converged zero-partition record. |

Changing any resolution requires a new profile identifier and independent
validation evidence.

## Validation and acceptance criteria

L-INS-i may be marked complete only when all of the following hold:

1. this document contains exact references, equations, original pseudocode,
   data structures, complexity, parameters, tie rules, omissions, and limits;
2. hand-computed tests verify both matrix offsets and all local affine states;
3. an independent exhaustive oracle verifies local optimum scores for all
   short DNA pairs, including every affine path and every possible local
   source range;
4. local traceback replay independently verifies score, ranges, positive
   prefixes, gap runs, and source consumption;
5. tests cover empty/no-hit, identical, divergent, repetitive, ambiguous,
   leading/internal/terminal-gap, and long-flank cases;
6. tests prove that terminal flanks are absent from `I` but retained exactly
   once in L-INS-1 and final L-INS-i output;
7. hand-computed tests verify local segmentation, original residue indices,
   `f`, `E`, sparse `I`, symmetry, and lack of library extension;
8. distance, mixed-linkage, Thompson weighting, profile bonus, and fixed-tree
   tests cover all shared [GL]-INS-i equations;
9. every progressive merge uses full complete profiles, no FFT, and preserves
   both child subalignments;
10. direct and cached `WSP+I` evaluators agree on substitutions and every gap
    orientation, including complete terminal flanks;
11. every refinement candidate preserves both fixed subalignments and every
    source sequence, and accepted cycles strictly improve one fixed objective;
12. protein, DNA, and RNA pass the same deterministic public pipeline;
13. property tests verify rectangularity, exact source recovery, finite
    evidence, local optimality on bounded inputs, objective monotonicity,
    symmetry metrics, and repeated-run determinism;
14. inert fixtures generated only by official opaque MAFFT 7.526 verify exact
    trivial cases and SP/TC agreement for nontrivial local/flank cases;
15. a fixed BAliBASE 3 subset including Ref4 records SP and TC without using
    this implementation to generate expected data;
16. Criterion records length, sequence-count, molecule-type, and flank-scaling
    costs, while checked memory evidence accounts separately for local-pair,
    library, progressive, and refinement phases;
17. matched fresh-process runs compare the isolated Rust worker with opaque
    MAFFT 7.526 using mirrored inputs, warm-up, raw samples, artifact hashes,
    correctness metrics, time, and scoped peak-memory observations;
18. formatting, warnings-denied Clippy, debug/release workspace tests,
    Rustdoc, and transitive advisory/license/source checks pass;
19. the public API and examples are deterministic, safe, English-only, and
    contain no placeholder, unfinished path, mock, or runtime dependency.

## Validation record

The `MAFFT-2005-L-INS-i-v1` profile satisfied the criteria above on
2026-08-05:

- the hand-derived unit, exhaustive-oracle, property, integration, and
  reference suites passed in both debug and release workspace builds;
- formatting, warnings-denied Clippy over every target and feature, Rustdoc
  with denied warnings, and the advisory/license/duplicate/source audit all
  passed;
- seven inert MAFFT 7.526 cases produced six exact alignments; the non-exact
  globin case scored SP `0.972696246` and TC `0.943396226`, while aggregate
  SP/TC were `0.996099464` / `0.991913747`. The fixture SHA-256 is
  `4c1273b81a396c19670094ba03c3a63ca3c8613b042345d435cfd838d45ea55a`;
- the independent seven-point `--lep` sweep resolved the local matrix offset
  sign without source inspection. Its inert fixture SHA-256 is
  `2b2e50e9ee9e5340a70985f9a2c43c018c0bcfe879796ddeb2d0ad8aa197b480`
  and the isolated raw record SHA-256 is
  `9bccaa57c0e623f929f1675fed4520020ea3f22af79ba8952eac473182884f8c`;
- the selected BAliBASE 3 R1-R5 subset recorded Rust macro SP/TC
  `0.870354493` / `0.590635999` and micro SP/TC `0.850494883` /
  `0.780632411`; both the ancestral result and the separate opaque MAFFT
  7.526 result are retained in the
  benchmark record;
- Criterion recorded deterministic length, count, terminal-flank, protein,
  DNA, and RNA workloads. Matched fresh-process checks recovered every source
  and achieved SP=TC=`1.0` on both comparison datasets; the complete protocol,
  raw ranges, hashes, and memory scopes are in the same benchmark record;
- targeted source scans found no unfinished implementation primitive, unsafe
  block, FFI declaration, process launch, or runtime reference dependency in
  the L-INS-i library path.

Accordingly, this implementation is complete for the explicitly
versioned `MAFFT-2005-L-INS-i-v1` profile, subject to the limitations below.
This status does not imply completion of E-INS-i, other MAFFT strategies, or
the wider algorithm catalog.

## Black-box reference boundary

The reference laboratory must use the official signed MAFFT 7.526 Windows
package already pinned for G-INS-i, without downloading or inspecting source.
The operational command is

```text
mafft.bat --localpair --maxiterate 1000 --thread 1 --inputorder --quiet
          --out OUTPUT INPUT
```

Inputs, executable, outputs, orchestration, and raw measurements remain in the
sibling validation repository. Only inert fixtures and reproducible aggregate
records may enter this repository. Exact equality is required for
unambiguous trivial cases; nontrivial evidence uses valid source recovery,
determinism, SP, TC, score diagnostics, runtime, and memory because the modern
executable is not the historical target.

## Known limitations of this profile

- It reconstructs the ancestral published strategy, not all current MAFFT
  behavior or an unavailable historical binary.
- The local recurrence, open-versus-extension convention, endpoint ties,
  ambiguity symbols, and exact local-distance normalization are incompletely
  specified by the authors; every choice is versioned above.
- Only one maximum local domain per pair is retained. L-INS-i assumes one
  alignable domain and is not the E-INS-i solution for several separated
  motifs or long internal unalignable regions.
- Long unrelated flanks remain in the final full-length MSA and may be aligned
  by the complete-profile stages even though they contribute no local library
  evidence. This is characteristic of L-INS-i, not local-output semantics.
- The ancestral nucleotide parameters predate later RNA-specific tuning and
  must not be represented as modern MAFFT defaults.
- Conserved-region anchoring is absent because its complete historical rule is
  unpublished; refinement may therefore use more work and explore
  candidates excluded by the historical heuristic.
- All-pairs local traceback and the retained library scale quadratically in
  sequence count. This profile targets the small and medium sets described by
  the authors, not HomFam-scale collections.
- The local, progressive, and individual refinement rectangles are scalar and
  quadratic in sequence or profile width. Refinement is parallel only across
  independent fixed-tree candidates; no SIMD or intra-candidate parallelism is
  claimed here.
