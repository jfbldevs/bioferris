# MAFFT E-INS-i

## Reproduction target

This document specifies the immutable scientific profile
`MAFFT-2008-E-INS-i-v1`.

The target is the E-INS-i strategy explicitly described by Katoh and Toh in
their 2008 review of MAFFT version 6. E-INS-i is a consistency-aware iterative
multiple sequence aligner for sequences that share several conserved motifs in
the same order but contain long, unrelated internal regions. Its all-pairs
stage uses Altschul's local generalized-affine model so that unrelated residues
between conserved blocks are consumed without declaring them homologous. The
remaining guide-tree, primary-library, complete-profile progressive, and
fixed-tree `WSP+I` refinement stages are the published G-INS-i framework.

The 2008 authors' review states that E-INS-i was added after the MAFFT version
5 publication. Consequently, this profile does not label the 2005 G-INS-i
paper itself as an E-INS-i publication. It uses that paper only for the shared
consistency pipeline and the 2008 review as the defining E-INS-i source.

This profile deliberately reproduces the ancestral parameter route. MAFFT
changed E-INS-i parameters in version 7.243 in June 2015; current MAFFT exposes
the earlier behavior through `--oldgenafpair`. The later parameterization is a
different future profile. No MAFFT source code, translated code, FFI,
subprocess, or runtime executable is part of this implementation.

## Primary references

1. Katoh K, Toh H. *Recent developments in the MAFFT multiple sequence
   alignment program*. Briefings in Bioinformatics. 2008;9(4):286-298.
   DOI: [10.1093/bib/bbn013](https://doi.org/10.1093/bib/bbn013).
   This is the first publication to describe E-INS-i explicitly. It identifies
   E-INS-i as a G-INS-i variant, specifies local generalized-affine all-pairs
   alignment, explains that internal unalignable regions remain unaligned in
   that stage, and states that the final MSA is nevertheless full length.
2. Altschul SF. *Generalized affine gap costs for protein sequence alignment*.
   Proteins. 1998;32(1):88-96.
   DOI:
   [10.1002/(SICI)1097-0134(19980701)32:1&lt;88::AID-PROT10&gt;3.0.CO;2-J](https://doi.org/10.1002/(SICI)1097-0134(19980701)32%3A1%3C88%3A%3AAID-PROT10%3E3.0.CO%3B2-J).
   This defines a gap containing residues from both sequences, its
   three-parameter cost, the path-graph diagonal skip, local use, and `O(mn)`
   dynamic-programming bound.
3. Katoh K, Kuma K, Miyata T, Toh H. *Improvement in the Accuracy of Multiple
   Sequence Alignment Program MAFFT*. Genome Informatics.
   2005;16(1):22-33.
   DOI: [10.11234/gi1990.16.22](https://doi.org/10.11234/gi1990.16.22).
   This defines the shared all-pairs segment library, residue frequency,
   importance equations, complete-profile progressive stage, and `WSP+I`
   iterative refinement.
4. Katoh K, Kuma K, Toh H, Miyata T. *MAFFT version 5: improvement in accuracy
   of multiple sequence alignment*. Nucleic Acids Research.
   2005;33(2):511-518.
   DOI: [10.1093/nar/gki198](https://doi.org/10.1093/nar/gki198).
   This is the immediate G-INS-i consistency-framework publication and
   parameter context; its H/F-INS-i external routes are not used.
5. Smith TF, Waterman MS. *Identification of common molecular subsequences*.
   Journal of Molecular Biology. 1981;147(1):195-197.
   DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5).
6. Gotoh O. *An improved algorithm for matching biological sequences*.
   Journal of Molecular Biology. 1982;162(3):705-708.
   DOI: [10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9).
7. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic
   Acids Research. 2002;30(14):3059-3066.
   DOI: [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
8. Thompson JD, Higgins DG, Gibson TJ. *CLUSTAL W: improving the sensitivity
   of progressive multiple sequence alignment through sequence weighting,
   position-specific gap penalties and weight matrix choice*. Nucleic Acids
   Research. 1994;22(22):4673-4680.
   DOI: [10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
9. Gotoh O. *A weighting system and algorithm for aligning many
   phylogenetically related sequences*. Computer Applications in the
   Biosciences. 1995;11(5):543-551.
   DOI: [10.1093/bioinformatics/11.5.543](https://doi.org/10.1093/bioinformatics/11.5.543).
10. Hirosawa M, Totoki Y, Hoshida M, Ishikawa M. *Comprehensive study on
    iterative algorithms of multiple sequence alignment*. Computer
    Applications in the Biosciences. 1995;11(1):13-18.
    DOI: [10.1093/bioinformatics/11.1.13](https://doi.org/10.1093/bioinformatics/11.1.13).
11. Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein
    blocks*. Proceedings of the National Academy of Sciences.
    1992;89(22):10915-10919.
    DOI: [10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915).
12. Kimura M. *A simple method for estimating evolutionary rates of base
    substitutions through comparative studies of nucleotide sequences*.
    Journal of Molecular Evolution. 1980;16:111-120.
    DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581).

Official operational evidence:

13. The official algorithm page describes E-INS-i, its ordered-motif
    assumption, the June 2015 parameter change, and `--oldgenafpair`:
    <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.
14. The official 2007 manual records the ancestral command and parameters:
    <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.

Items 13-14 establish public parameter names and intended behavior. They do
not replace the papers or make current executable output the scientific
definition.

## Objective and characteristic pipeline

E-INS-i targets protein, DNA, or RNA sequences with multiple homologous blocks
in a conserved order, separated by regions for which residue correspondence is
not supported. Its characteristic pipeline is:

1. align every unordered sequence pair with scalar local DP containing
   ordinary affine indels and one generalized unaligned-region state;
2. derive a normalized score distance from each generalized local optimum;
3. build one mixed-linkage UPGMA-like guide tree;
4. compute Thompson/CLUSTALW sequence weights from that tree;
5. retain only aligned residue-residue runs as the primary consistency
   segments; ordinary gaps and generalized skipped residues separate runs;
6. derive residue frequencies, segment importance values, and the sparse
   primary importance matrix `I`;
7. progressively align every complete input sequence along the fixed tree with
   the E-INS-i complete-profile score, producing E-INS-1;
8. derive fixed Gotoh weights and every canonical tree-edge partition;
9. realign fixed subalignments and accept only the strictly best `WSP+I`
   improvement per cycle;
10. stop at convergence or after 1,000 accepted cycles.

Only the all-pairs stage leaves internal residues unaligned. E-INS-1 and the
final E-INS-i MSA contain every source residue exactly once. A generalized
skip-pair operation consumes one residue from each sequence for cost purposes,
but explicitly does not assert a homologous residue pair and cannot enter the
consistency library.

Replacing the generalized-affine stage with one Smith-Waterman hit, multiple
independent Smith-Waterman hits, Waterman-Eggert enumeration, DIALIGN fragments,
an HMM, or an external program would be a different algorithm. Aligning skipped
residue pairs into one output column would also silently change their meaning.

## Versioned parameters

| Parameter | Value | Source and role |
| --- | ---: | --- |
| Profile identifier | `MAFFT-2008-E-INS-i-v1` | This specification |
| Maximum accepted refinement cycles | 1,000 | Published/official E-INS-i command |
| Importance multiplier `W^I` | 2.7 | 2005 consistency paper/manual |
| Mixed-linkage arithmetic fraction | 0.1 | Inherited published MAFFT linkage |
| Complete-profile gap opening | 1.53 | Ancestral `--op` |
| Complete-profile matrix offset | 0.0 | E-INS-i recommendation `--ep 0` |
| Local ordinary gap opening | -2.0 | Ancestral signed `--lop` |
| Local matrix offset penalty | 0.1 | Public `--lep`, subtracted from the normalized matrix |
| Local ordinary/unpaired extension | -0.1 | Ancestral signed `--lexp` |
| Generalized skip opening | -6.0 | Ancestral signed `--LOP` |
| Generalized paired-skip extension | 0.0 | Ancestral signed `--LEXP` |
| Gotoh equalization exponent | 1.15 | Shared versioned [GLE]-INS-i resolution |
| Zero-branch floor | `1e-12` | Shared finite weighting resolution |

Protein pair scores use normalized BLOSUM62. DNA and RNA use the same
normalized Kimura two-parameter matrix with transition/transversion ratio two.
Ambiguity symbols are canonical expectations. The pair matrix is
`M_hat_-0.1`; the complete-profile and final-objective matrix is `M_hat_0.0`.
The ancestral common protein/nucleotide penalties are retained instead of
claiming later RNA tuning.

The signed generalized values correspond to Altschul's non-negative cost
triple as follows:

```text
a = 6.0, b = 0.1, c = 0.0
generalized_gap_score(k1,k2) =
    -a - b * |k1-k2| - c * min(k1,k2).
```

This path competes with the cheaper ordinary affine model for one-sided indels.
An ordinary run of length `k >= 1` scores `-2.0 - 0.1*(k-1)` under MAFFT's
versioned open-only convention. In contrast, every generalized skipped residue
or residue pair pays its extension contribution in addition to the one
`-6.0` opening, matching Altschul's convention that the opening does not
include the first skipped item.

## Scalar generalized-affine local recurrence

Let `x[1..m]` and `y[1..n]` be sequences and let
`S(i,j)=M_hat_-0.1(x_i,y_j)`. Four states are used:

- `H(i,j)`: best local path ending by aligning `x_i` with `y_j`;
- `X(i,j)`: path inside an ordinary first-sequence-only gap;
- `Y(i,j)`: path inside an ordinary second-sequence-only gap;
- `Z(i,j)`: path inside one generalized unaligned region.

Unavailable boundary states are negative infinity. With
`o=-2.0`, `e=-0.1`, `O=-6.0`, and `E=0.0`:

```text
H(i,j) = S(i,j) + max(
    0,
    H(i-1,j-1),
    X(i-1,j-1),
    Y(i-1,j-1),
    Z(i-1,j-1))

X(i,j) = max(
    H(i-1,j) + o,
    X(i-1,j) + e)

Y(i,j) = max(
    H(i,j-1) + o,
    Y(i,j-1) + e)

Z(i,j) = max(
    H(i-1,j-1) + O + E,   # open with one unaligned pair
    H(i-1,j)   + O + e,   # open with one unpaired x residue
    H(i,j-1)   + O + e,   # open with one unpaired y residue
    Z(i-1,j-1) + E,       # extend with one unaligned pair
    Z(i-1,j)   + e,       # extend with one unpaired x residue
    Z(i,j-1)   + e).      # extend with one unpaired y residue
```

The optimum is the maximum positive `H` cell. A path cannot end in `X`, `Y`,
or `Z`, because trailing unaligned residues do not improve a local score.
Likewise, no skipped prefix is needed: `H` can restart at any residue pair.
The generalized recurrence evaluates `O(mn)` cells and exactly realizes
Altschul's gap cost while retaining a distinct cheaper state for ordinary
biological indels.

### Traceback and tie-breaking

The first strict row-major maximum is the endpoint. A non-positive predecessor
of `H` selects restart. Equal positive `H` predecessors prefer direct aligned
continuation, ordinary first-only gap, ordinary second-only gap, then
generalized skip. Ordinary gap ties prefer opening over extension. Generalized
state ties prefer all opening moves over all extension moves; within either
class they prefer paired, first-only, then second-only consumption. These rules
are implementation-profile choices because the publications omit them.

Traceback reports half-open source ranges and six operation kinds:

```text
Aligned
OrdinaryFirstOnly
OrdinarySecondOnly
SkipPair
SkipFirstOnly
SkipSecondOnly
```

`SkipPair` consumes both residues but creates no homology. For the temporary
two-row profile it is serialized deterministically as a first-only column then
a second-only column. Thus no downstream generic profile routine can mistake
the pair for an aligned residue pair. Only `Aligned` operations form maximal
gap-free library segments. An independent replay checks state legality, exact
range consumption, one opening per run, every extension contribution, local
prefix positivity, and the terminal score without reading DP buffers.

When no positive aligned cell exists, the canonical result has score zero,
ranges `0..0`, no operations, and no consistency segment.

## Pair-score distance and primary library

For generalized pair score `t_st`, self scores use complete ungapped sequences
and the same local matrix:

```text
t_ss = sum_p M_hat_-0.1(s[p],s[p])
t_tt = sum_q M_hat_-0.1(t[q],t[q])
D_st = clamp(1 - t_st / min(t_ss,t_tt), 0, 1).
```

Empty and non-positive normalizations use the explicit totalization defined by
the L-INS-i profile. A generalized path may connect several conserved blocks,
so `t_st` includes their substitution scores, ordinary gaps, and generalized
skip costs. Skipped terminal residues remain outside the local range.

Every maximal run of `Aligned` operations becomes segment `n`:

```text
P(s,t,n) = ordered original residue-index pairs
L(s,t,n) = |P(s,t,n)|
S(s,t,n) = sum over (p,q) in P(s,t,n)
               M_hat_-0.1(s[p],t[q]).
```

The published 2005 equations are then applied unchanged:

```text
f(s,p) = sum over n,t,q where (s,t,p,q) in P(s,t,n) of w_t

E(s,t,n) = S(s,t,n)
           * sum over (i,j) in P(s,t,n)
               [f(s,i) + f(t,j)] / [2 L(s,t,n)]

I(s,t,p,q) = sum over n containing (p,q) of E(s,t,n), or 0.
```

Skipped residues contribute to neither `P`, `S`, `f`, `E`, nor `I`. There is no
three-sequence library extension.

## Shared tree, progressive, and refinement stages

The fixed mixed-linkage update is

```text
D_PC = 0.1 * (D_LC + D_RC) / 2 + 0.9 * min(D_LC,D_RC).
```

Tree construction, deterministic ties, Thompson weights, root suppression,
Gotoh weights, and canonical partitions are identical to the documented
G-INS-i/L-INS-i implementation.

For complete profile groups `G_i` and `G_j`, columns `p` and `q` receive

```text
H_group(G_i,G_j,p,q) =
    sum over s in G_i, t in G_j of w_st *
      [M_hat_0.0(A(s,p),A(t,q))
       + 2.7 * I(s,t,r_s(p),r_t(q))].
```

The extension-free position-specific complete-profile gap opening is `1.53`.
Full scalar DP rectangles are used because the complete historical anchor rule
is unpublished. Every input residue, including all generalized skipped
regions, enters E-INS-1 and all refinement candidates.

The final fixed objective is

```text
WSP+I(A) = sum over s<t of w_st *
  [pair_score_0.0(A_st) + 2.7 * importance(A_st)],
```

where `pair_score_0.0` uses `M_hat_0.0`, one `1.53` penalty per maximal
orientation-specific residue-gap run, and includes terminal runs. Candidate
selection, exact caching, direct replay, strict `>` acceptance, and convergence
rules are shared with the [GLE]-INS-i profiles.

## Pseudocode

```text
E_INS_I(records):
    validate non-empty records and unique identifiers
    pair_scorer  <- normalized BLOSUM62 or K80 with offset -0.1
    group_scorer <- normalized BLOSUM62 or K80 with offset 0.0

    for every unordered pair (s,t) in lexical input order:
        generalized <- LOCAL_GENERALIZED_AFFINE(
            records[s], records[t], pair_scorer,
            ordinary_open=-2.0,
            ordinary_extension=-0.1,
            skip_open=-6.0,
            paired_skip_extension=0.0)
        INDEPENDENTLY_REPLAY_GENERALIZED(generalized)
        segments <- MAXIMAL_ALIGNED_RUNS_ONLY(generalized)
        distance[s,t] <- TOTALIZED_NORMALIZED_SCORE_DISTANCE(
            generalized.score, self_score[s], self_score[t])

    tree <- MIXED_LINKAGE_TREE(distance, arithmetic_fraction=0.1)
    thompson <- CLUSTALW_PATH_WEIGHTS(tree)
    library <- PRIMARY_IMPORTANCE_LIBRARY(segments, thompson)

    initial <- PROGRESSIVE_POSTORDER(tree):
        align complete child profiles with
        M_hat_0.0 + 2.7 * library bonus and gap_open=1.53

    gotoh <- FIXED_GOTOH_WEIGHTING(tree, equalization=1.15)
    cache <- EXACT_WSP_PLUS_I_CACHE(initial, gotoh, library, M_hat_0.0)
    current <- initial

    repeat until convergence or 1000 accepted cycles:
        evaluate every canonical fixed-tree split from current
        preserve both projected subalignments and every source residue
        select the first strictly best WSP+I improvement
        if none improves: stop
        current <- selected candidate
        cache <- independently verified candidate contributions

    verify final cache with a direct WSP+I evaluator
    return complete MSA and pair/tree/library/refinement/memory evidence
```

## Data structures

- four rolling generalized-local score rows and packed per-state traceback;
- typed six-way operation path plus two half-open source ranges;
- explicit generalized-run statistics: openings, paired skips, and unpaired
  residues;
- temporary full-source pair profile that never co-columns a skipped pair;
- raw aligned segments with original residue indices;
- dense finite distance matrix and one fixed mixed-linkage guide tree;
- sparse ordered primary importance library and residue frequencies;
- complete typed progressive/refinement profiles and transformed score tables;
- exact per-pair `WSP+I` contribution cache;
- auditable pair, merge, cycle, termination, and checked memory records;
- final typed `MultipleAlignment` in original input order.

## Complexity

Let `N` be sequence count, `L` maximum ungapped length, `W` an MSA width, and
`C` accepted refinement cycles.

- one generalized local pair: `O(L^2)` time, `O(L)` rolling scores, and
  `O(L^2)` packed traceback;
- all generalized pairs: `O(N^2 L^2)` time, with sequential pair workspace
  and retained `O(N^2 L)` normal-case segment evidence;
- mixed-linkage tree: `O(N^3)` scalar time and `O(N^2)` space;
- progressive stage: `N-1` full profile rectangles, worst-case `O(NW^2)` time;
- one attempted refinement cycle: the shared conservative
  `O(N^3 W + N W^2)` bound over `2N-3` candidates;
- total refinement: attempted-cycle cost times at most `C+1`;
- returned evidence: `O(NW + CN + N^2)` excluding the primary library.

The correctness profile is scalar. SIMD, pair parallelism, candidate
parallelism, banding, anchor restrictions, and linear-space traceback are
future optimizations and may not alter results.

## Published omissions and explicit resolutions

| Omission or ambiguity | `v1` resolution |
| --- | --- |
| First explicit E-INS-i version | Target the 2008 version-6 description, not the 2005 G/H/F paper or post-7.243 route. |
| Exact MAFFT generalized recurrence | Reconstruct Altschul's four-move path with separate ordinary `X/Y` and generalized `Z` states as written above. |
| Mapping of `lop/lexp/LOP/LEXP` | Ordinary gaps use `lop/lexp`; generalized unequal and paired skip moves use `lexp` and `LEXP`, respectively. |
| Opening convention | Ordinary gaps use the shared MAFFT open-only convention; generalized gaps pay `LOP` plus the first extension, as Altschul states. |
| `LEXP=0` behavior | A paired skipped residue is free after the one generalized opening; it is consumed but never homologized. |
| Public `lep=0.1` sign | Subtract the positive penalty, supported by the independent L-INS-i option sweep and reused by the same documented local parameter. |
| Complete-profile `ep` | Use the specifically recommended E-INS-i value zero. |
| Local endpoint and state ties | First strict row-major endpoint and the deterministic ordering above. |
| Generalized skip serialization | Emit two one-sided columns for `SkipPair`; never place the residues in one homology column. |
| Multiple conserved motifs | One optimal generalized local path may contain multiple aligned segments in conserved order. No independent-hit enumeration. |
| Pair distance | Use the inherited complete-self-score normalization and explicit clamp. |
| Progressive/refinement stages | Reuse the exact shared G-INS-i pipeline with group offset zero; all residues are full length. |
| Conserved anchors | Omit unpublished restrictions and evaluate full rectangles. |
| Nucleotide parameters | Retain the ancestral common values; document that these are not modern RNA defaults. |
| Modern comparator | Invoke MAFFT 7.526 opaquely with `--oldgenafpair`; current default E-INS-i is not the target. |

Changing any resolution requires a new profile identifier and independent
validation evidence.

## Validation and acceptance criteria

E-INS-i may be marked complete only when all of the following hold:

1. this specification contains references, equations, original pseudocode,
   data structures, complexity, parameters, ties, omissions, and limits;
2. hand-computed tests cover every transition and both opening conventions;
3. an independent exhaustive path oracle matches all short DNA scores and
   validates every legal local range and generalized skip shape;
4. independent replay verifies the score, state transitions, source ranges,
   exact consumption, gap openings, extensions, and positive local prefixes;
5. tests cover empty/no-hit, identical, divergent, repetitive, ambiguous,
   leading/internal/terminal indels, and unequal long internal regions;
6. tests prove that two or more ordered conserved motifs can enter one pair
   path while intervening skipped residues never enter the library;
7. source recovery proves that every skipped and terminal residue appears
   exactly once in E-INS-1 and final E-INS-i output;
8. hand-computed tests verify segment boundaries, original indices, `f`, `E`,
   sparse `I`, symmetry, and absence of library extension;
9. shared distance, tree, weighting, profile, objective, and fixed-subalignment
   invariants remain green with group offset zero;
10. protein, DNA, and RNA pass one deterministic public API;
11. property tests cover rectangularity, source recovery, finite evidence,
    bounded-input optimality, score replay, objective monotonicity, and repeated
    determinism;
12. an inert official-option sweep records the observable ancestral parameter
    signs and the switch between affine and generalized behavior;
13. opaque MAFFT 7.526 `--oldgenafpair` fixtures gate exact trivial cases and
    SP/TC agreement on nontrivial multiple-motif cases;
14. an appropriate scientific dataset, including long-insertion/internal-gap
    cases, records independent SP and TC without self-generated expectations;
15. Criterion records length, count, molecule, internal-region length, and
    internal-length-imbalance scaling plus checked memory;
16. matched fresh-process comparisons with opaque MAFFT record source hashes,
    raw samples, correctness metrics, runtime, and scoped memory;
17. formatting, warnings-denied Clippy, debug/release workspace tests,
    warnings-denied Rustdoc, and transitive dependency audits pass;
18. the public API and examples are safe, deterministic, English-only, and
    contain no unfinished path or runtime external dependency.

## Implementation and validation record

Status: **complete for `MAFFT-2008-E-INS-i-v1`** as of 2026-08-05. This status
applies only to the immutable ancestral profile defined in this document. It
does not claim reproduction of MAFFT's post-7.243 parameterization or any
undocumented implementation detail.

The implementation is split into a scalar generalized-affine pair engine,
the typed E-INS-i public facade, and the shared versioned [GLE]-INS-i
consistency pipeline. Twenty E-INS-i-specific test functions cover the
generalized recurrence, exhaustive short-path optimality, independent replay,
full-source recovery, profile invariants, supported alphabets, deterministic
properties, opaque parameter observations, and MAFFT comparison fixtures.
Shared tree, weighting, profile, objective, and refinement tests also run for
the E-INS-i route.

The pinned eight-case MAFFT 7.526 `--oldgenafpair` suite has mean
SP `0.930925786`, mean TC `0.826272888`, and four exact alignments. Its inert
fixture SHA-256 is
`20d0be6db389901b27c9482bfe7ba7aa73d75ac51556375dc103e7172e2e0af8`.
The independent ancestral option-sweep fixture SHA-256 is
`f06fc949998ea77dbd2c91ebcaa46d862f2ac91c53d64f4e06194ae05da40991`.
Both generators reproduced their pinned records byte for byte in the isolated
validation laboratory.

The selected BAliBASE 3 R1-R5 protein subset records macro SP/TC
`0.870354493 / 0.590635999` and micro SP/TC
`0.850494883 / 0.780632411`. The official opaque MAFFT comparator records
macro SP/TC `0.918443969 / 0.701347660` and micro SP/TC
`0.904680423 / 0.835968379` on the same references. These are measurements,
not identity requirements. Full Criterion, matched-process, memory, artifact,
and dataset details are recorded in
`docs/benchmarks/mafft-e-ins-i.md`.

The following completion gates passed on 2026-08-05:

```text
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --release --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked
cargo deny check advisories bans licenses sources
cargo run -p align-mafft --example einsi --locked
cargo bench -p align-benchmarks --bench mafft_einsi --locked -- --noplot
```

No runtime dependency, FFI, subprocess, sidecar, external binary, translated
source, placeholder, mock, or unfinished execution path is used. External
programs remain confined to the sibling black-box validation laboratory.

## Black-box reference boundary

The pinned official MAFFT 7.526 Linux and Windows packages may be invoked only
in the isolated sibling validation laboratory. The ancestral comparison command
is

```text
mafft --genafpair --oldgenafpair --ep 0 --maxiterate 1000
      --thread 1 --inputorder --quiet INPUT
```

The source is not downloaded or inspected. Executables, inputs, raw output,
and orchestration remain outside this repository. Only inert fixtures and
reproducible aggregate records enter the project. Exact equality is required
for unambiguous cases; nontrivial cases use source recovery, determinism, SP,
TC, score diagnostics, time, and scoped memory because even the old-parameter
route runs inside a much newer engineered executable.

## Known limitations of this profile

- It reconstructs the published ancestral strategy and explicit Altschul
  model, not undocumented implementation details or post-7.243 E-INS-i.
- The MAFFT authors publish parameter names but omit the combined ordinary and
  generalized recurrence, traceback ties, exact local normalization edge
  cases, and internal serialization. Every resolution is versioned above.
- E-INS-i assumes that conserved motifs occur in the same order. It does not
  model rearrangements, inversions, circular permutations, duplications, or
  genome-scale synteny.
- A skipped pair is a cost-model device, not evidence of homology. The
  deterministic two-column serialization is intentionally not a biological
  alignment of those residues.
- The final MSA is full length. Complete-profile progressive/refinement stages
  can place residues from pairwise-unrelated regions into common columns even
  though those columns receive no primary-library support.
- The ancestral nucleotide parameters predate later RNA-specific tuning.
- All-pairs DP and retained segment evidence scale quadratically in sequence
  count; scalar traceback and profile rectangles scale quadratically in width.
- No SIMD, parallel all-pairs calculation, generalized banding, or historical
  conserved-anchor heuristic is claimed.
