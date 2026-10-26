# MAFFT FFT-NS-1

## Reproduction target

This crate reproduces the single-pass progressive strategy described by Katoh,
Misawa, Kuma, and Miyata in 2002. The versioned scientific contract is named
`MAFFT-2002-FFT-NS-1-v1`.

The contract deliberately does not combine the 2002 method with defaults from a
later MAFFT release. In particular, MAFFT 7.526 defaults to BLOSUM62,
`--op 1.53`, and `--ep 0.123`, corrects the length dependence of the 6-mer
distance, and documents a modified UPGMA procedure. Those later choices are
useful for black-box comparison, but they are not the algorithm specified by
the original publication.

No MAFFT source code was downloaded or inspected. The implementation is derived
from the publications, their displayed equations, the official manual and
algorithm description, and independent tests. The official MAFFT 7.526 binary
is used only in the isolated reference-validation laboratory.

## Primary references

1. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic Acids
   Research. 2002;30(14):3059-3066.
   [doi:10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
2. Jones DT, Taylor WR, Thornton JM. *The rapid generation of mutation data
   matrices from protein sequences*. Computer Applications in the Biosciences.
   1992;8(3):275-282.
   [doi:10.1093/bioinformatics/8.3.275](https://doi.org/10.1093/bioinformatics/8.3.275).
3. Grantham R. *Amino acid difference formula to help explain protein
   evolution*. Science. 1974;185(4154):862-864.
   [doi:10.1126/science.185.4154.862](https://doi.org/10.1126/science.185.4154.862).
4. Kimura M. *A simple method for estimating evolutionary rates of base
   substitutions through comparative studies of nucleotide sequences*. Journal
   of Molecular Evolution. 1980;16:111-120.
   [doi:10.1007/BF01731581](https://doi.org/10.1007/BF01731581).
5. Thompson JD, Higgins DG, Gibson TJ. *CLUSTAL W: improving the sensitivity of
   progressive multiple sequence alignment through sequence weighting,
   position-specific gap penalties and weight matrix choice*. Nucleic Acids
   Research. 1994;22(22):4673-4680.
   [doi:10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
6. Cooley JW, Tukey JW. *An algorithm for the machine calculation of complex
   Fourier series*. Mathematics of Computation. 1965;19(90):297-301.
   [doi:10.1090/S0025-5718-1965-0178586-1](https://doi.org/10.1090/S0025-5718-1965-0178586-1).
7. Katoh K. *MAFFT manual*. The documented FFT-NS-1 invocation is
   `mafft --retree 1 --maxiterate 0`; the same manual identifies shared 6-mer
   distance and FFT profile alignment.
   <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.
8. MAFFT project. *Details of the algorithms*. This post-publication record is
   used only to identify later changes, not to replace the 2002 equations.
   <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.

The 2002 article has no separate supplementary algorithm file. The displayed
equations and prose in the article are therefore the highest-priority complete
description. The official MAFFT algorithm page and manual are secondary records
used to identify the command corresponding to FFT-NS-1 and to enumerate later
behavioral changes.

## Objective

FFT-NS-1 is a fast progressive multiple-sequence aligner. It has three defining
stages:

1. estimate every pairwise distance by shared 6-tuples;
2. construct one UPGMA guide tree;
3. progressively align profiles in tree order, using Fourier correlation to
   locate homologous blocks and divide each profile-profile dynamic-programming
   matrix.

There is exactly one guide-tree construction and one progressive pass. A second
tree/pass would be FFT-NS-2. Iterative realignment would be FFT-NS-i. Neither is
part of this phase.

## Alphabet model

The public implementation supports `Protein`, `Dna`, and `Rna` through the same
typed API.

### Protein 6-tuple compression

Canonical residues are compressed to the six Dayhoff classes used by the
publication:

| Code | Residues |
| --- | --- |
| 0 | `A G P S T` |
| 1 | `D E N Q` |
| 2 | `H K R` |
| 3 | `I L M V` |
| 4 | `F W Y` |
| 5 | `C` |

`B`, `Z`, and `J` map unambiguously to classes 1, 1, and 3. `U` uses the
cysteine class and `O` uses the lysine class. `X` invalidates a tuple window.
These extended-symbol rules are explicit project resolutions; the 2002 paper
only specifies the 20 canonical residues.

### Nucleic-acid tuples

DNA and RNA use four canonical states. `T` and `U` occupy the same state.
Windows containing an IUPAC ambiguity symbol are omitted from the tuple
histogram. Ambiguity symbols remain supported in profile scoring through equal
mixtures of their compatible canonical bases.

## Rough distance matrix

For sequence `i`, let `C_i(q)` be the number of occurrences of compressed tuple
`q`. The number of shared tuples is the multiset intersection

```text
T_ij = sum_q min(C_i(q), C_j(q)).
```

The self-count is `T_ii = sum_q C_i(q)`. Equation 8 of Katoh et al. is

```text
D_ij = 1 - T_ij / min(T_ii, T_jj).
```

This profile uses the uncorrected 2002 distance. It does not apply the empirical
length correction introduced in MAFFT 6.

The paper assumes ordinary biological sequences long enough to contain a
6-tuple. For a collection containing a non-empty sequence shorter than six, the
project uses one common tuple length

```text
k = min(6, shortest non-empty input length).
```

An empty sequence has no tuples and has distance one from a non-empty sequence.
Two tuple-empty byte-identical sequences have distance zero; other tuple-empty
pairs have distance one. This is a documented totalization, not a claim about
the historical executable.

## Guide tree and sequence weights

Standard size-weighted UPGMA is applied to the rough distance matrix. If clusters
`A` and `B` merge, their distance to another cluster `C` is

```text
d(A union B, C) = (|A| d(A,C) + |B| d(B,C)) / (|A| + |B|).
```

The parent height is `d(A,B)/2`; child-edge lengths are parent height minus
child height. Cluster pairs with equal distance are ordered by the pair of
minimum input indices. Children are stored with the smaller minimum input index
on the left.

The 2002 article states that progressive sequence weights are calculated in the
same manner as CLUSTALW. For leaf `s`, this implementation uses the CLUSTALW
tree-path allocation

```text
w_s = sum_(edge e on root-to-s path) max(length(e), 0) / leaves_below(e).
```

Weights are normalized within each profile before correlation, match scoring,
or gap-frequency calculation. If the complete tree has zero length, equal
weights are used. A zero individual weight is floored to machine-stable positive
mass before normalization so duplicate sequences cannot disappear from a
profile.

## Fourier homology search

### Protein signals

For canonical amino acid `a`, Grantham volume `v(a)` and polarity `p(a)` are
standardized over the 20 canonical amino acids:

```text
v_hat(a) = (v(a) - mean(v)) / sd(v)
p_hat(a) = (p(a) - mean(p)) / sd(p).
```

A profile column is the weighted sum of its rows' two-component signals. Gaps
contribute zero. Extended residues contribute the equal mixture of their
documented canonical analogues.

### Nucleotide signals

A nucleotide profile column has four weighted frequency channels, one each for
`A`, `C`, `G`, and `T/U`. IUPAC symbols distribute their weight equally across
compatible bases; gaps contribute zero.

### Cross-correlation

For signal channel `r`, profiles `A` and `B`, and lag `k`, the paper defines

```text
c_r(k) = sum_n A_r(n) B_r(n + k),
c(k)   = sum_r c_r(k).
```

The implementation evaluates this linear cross-correlation as a convolution of
`reverse(A_r)` and `B_r`. Each convolution uses the exactly pinned RustFFT 6.4.1
complex transform with zero padding to the next power of two at least
`len(A) + len(B) - 1`. Plans are cached by transform length and processing uses
per-call buffers and scratch, so concurrent candidates share immutable plans but
not mutable workspaces. RustFFT exposes runtime-dispatched SIMD through a safe API
and retains a scalar fallback; it supplies no alignment logic.

The 20 highest correlation lags are retained. Correlation ties prefer the
smaller absolute lag and then the smaller signed lag. This ordering is not
specified by the paper.

### Homologous segments

At each retained lag, a 30-column sliding window is scored with the same
profile-profile homology score used by dynamic programming. A window is accepted
when its score exceeds `30 * 0.7`. Consecutive accepted windows on the same lag
are merged. Merged segments longer than 150 columns are split into consecutive
pieces of at most 150 columns, as specified in the paper.

Candidate segments from all lags are ordered by a two-dimensional weighted
interval-chain dynamic program. A predecessor must end no later than the next
segment starts in both profiles. The highest total segment score is selected;
ties prefer more covered columns and then lexicographically earlier coordinates.
This is the project's explicit formulation of the paper's segment-level DP,
whose boundary and tie details are not published.

The center column pair of every selected segment is fixed as a homologous pair.
The remaining profile-profile matrix is divided at those centers and each
rectangle is aligned independently. Fixing the centers makes the published
matrix division executable and prevents a partition from crossing the detected
homologous block.

## Similarity matrices

Katoh et al. normalize a raw log-odds matrix `M` as

```text
average1 = sum_a f_a M_aa
average2 = sum_a sum_b f_a f_b M_ab
M_hat_ab = (M_ab - average2) / (average1 - average2) + S_a.
```

Consequently, the expected score of two identical random draws is `1 + S_a`
and the expected score of independent random draws is `S_a`. The paper defaults
are `S_a = 0.06` and gap opening `S_op = 2.4`.

### Protein matrix

The JTT one-PAM mutation-probability matrix and residue frequencies are taken
from Tables II and III of Jones et al. Each probability row is renormalized
after the paper's five-decimal rounding, and matrix exponentiation produces the
200-PAM transition matrix. For target residue `a` and source residue `b`, the raw
score is

```text
M_ab = ln(P_200(a | b) / f_a).
```

The rounded published one-PAM table introduces very small reversibility error,
so the two orientations are averaged before Katoh normalization. The logarithm
base and a common multiplicative factor cancel in that normalization.

### Nucleotide matrix

DNA and RNA use a Kimura two-parameter transition matrix at 200 PAM, equal base
frequencies, and transition/transversion rate ratio two. With instantaneous
transition rate `alpha`, per-transversion rate `beta`, `alpha/beta = 2`, and
expected substitutions per site `(alpha + 2 beta)t = 2`,

```text
P_same       = 1/4 + 1/4 exp(-4 beta t)
                     + 1/2 exp(-2 (alpha + beta)t)
P_transition = 1/4 + 1/4 exp(-4 beta t)
                     - 1/2 exp(-2 (alpha + beta)t)
P_transvert   = 1/4 - 1/4 exp(-4 beta t).
```

Raw scores are `ln(P_ab / 0.25)` and are normalized by the same Katoh equation.
The phrase “transition/transversion ratio of 2.0” does not state whether the
ratio is instantaneous or an observed finite-time count. This profile fixes the
standard K80 instantaneous-rate interpretation above.

Ambiguity-symbol scores are expectations over equal canonical mixtures. Protein
`B`, `Z`, and `J` use their two canonical residues; `X` uses all 20; `U` uses
cysteine; and `O` uses lysine. This behavior is not specified in the 2002 paper.

## Profile homology score

For normalized profile weights `w_n` and `w_m`, columns `i` and `j` have score

```text
H(i,j) = sum_n sum_m w_n w_m M_hat[A(n,i), B(m,j)],
```

where a residue-gap or gap-gap pair contributes zero. Equivalently, the
implementation multiplies the two weighted residue-frequency vectors by the
normalized matrix.

## Position-specific, extension-free gap cost

The publication's recurrence is

```text
P(i,j) = H(i,j) + max(
    P(i-1,j-1),
    max_x(P(x,j-1) - G_1(i,x)),
    max_y(P(i-1,y) - G_2(j,y))
).
```

There is no explicit extension cost. For a new gap spanning an existing profile
from boundary `x` to boundary `i`,

```text
G_1(i,x) = S_op * (1 - (g_start(x) + g_end(i)) / 2).
```

`g_start` is the normalized weight of rows whose existing gap begins at the
left boundary, and `g_end` is the normalized weight of rows whose existing gap
ends at the right boundary. `G_2` is symmetric. Terminal outside cells are
treated as residues solely for detecting an existing terminal-gap boundary.

Because the start and end terms are separable, prefix maxima evaluate the
published all-gap-length recurrence in `O(L_1 L_2)` time. The implementation
stores the selected origin for traceback. The independent cubic recurrence is
retained in tests as an oracle.

Global terminal gaps pay one opening cost and no extension cost. This boundary
condition is consistent with a global Needleman-Wunsch alignment, but the paper
does not publish its initialization. Empty/empty scores zero; aligning a
non-empty profile with an empty one introduces exactly one gap opening, adjusted
by existing terminal-gap evidence.

## Tie-breaking

The historical executable's tie rules are not published. This implementation
uses one total deterministic policy:

1. higher score;
2. diagonal profile transition;
3. gap in the right profile;
4. gap in the left profile;
5. for equal transitions, the shortest newly introduced gap;
6. lower original sequence indices for tree, segment, and output ordering.

FFT-lag and segment-chain ties use the rules stated in their sections. No hash
iteration order, thread scheduling, wall clock, or random number affects output.

## Pseudocode

```text
FFT_NS_1(records):
    validate non-empty collection and unique identifiers
    k <- min(6, shortest non-empty input length), or 1 if all are empty

    for each sequence:
        tuples <- sorted compressed k-tuples, omitting invalid windows
    for each sequence pair i,j:
        T_ij <- multiset_intersection(tuples_i, tuples_j)
        D_ij <- published normalized shared-tuple distance

    tree <- deterministic size-weighted UPGMA(D)
    weights <- CLUSTALW path weights(tree)
    scorer <- normalized JTT200 or K80-200 matrix

    ALIGN_NODE(node):
        if node is leaf:
            return one-row ungapped profile with its tree weight

        left  <- ALIGN_NODE(node.left)
        right <- ALIGN_NODE(node.right)

        signals_left  <- weighted profile signals(left)
        signals_right <- weighted profile signals(right)
        correlations <- FFT_cross_correlations(signals_left, signals_right)
        lags <- best_20(correlations)

        segments <- sliding_windows(lags, width=30, threshold=0.7)
        segments <- merge_successive_and_split_at_150(segments)
        chain <- maximum_weight_monotone_segment_chain(segments)

        operations <- []
        for each selected segment center in order:
            operations += PROFILE_DP(rectangle before center)
            operations += DIAGONAL(center)
        operations += PROFILE_DP(trailing rectangle)

        return apply_operations(left, right, operations)

    profile <- ALIGN_NODE(tree.root)
    reorder profile rows to original input order
    independently validate rectangularity, symbol consumption, and score replay
    return alignment, guide tree, distances, weights, score, merge diagnostics
```

## Data structures

- sorted integer tuple-code arrays for sparse multiset intersection;
- a dense symmetric `DistanceMatrix`;
- deterministic binary `GuideTree` nodes with UPGMA heights;
- internal weighted profile rows containing typed optional symbols;
- dense per-column weighted residue frequencies;
- complex transform buffers, per-call scratch, and immutable cached RustFFT plans;
- homologous segment records and predecessor indices;
- rolling DP score rows plus one compact traceback record per cell;
- public merge diagnostics recording FFT candidates, selected segments, and
  profile dimensions.

## Complexity

Let `K` be the number of sequences, `L` a representative ungapped length, `W`
an intermediate profile width, and `R` the alphabet size.

- tuple extraction: `O(KL log L)` for sorted tuple codes;
- all-pairs tuple distance: `O(K^2 L)`;
- simple UPGMA implementation: `O(K^3)` time and `O(K^2)` memory;
- one profile FFT search: `O(C W log W)`, where `C=2` for protein and `C=4`
  for nucleic acid;
- one partitioned profile DP: at most `O(W_left W_right R^2)` to form homology
  scores and `O(W_left W_right)` for the separable recurrence;
- traceback memory: `O(W_left W_right)` for the largest partition;
- final alignment storage: `O(KW_final)`.

The FFT partition often reduces the DP area, but the worst case remains the full
quadratic profile matrix. UPGMA is intentionally scalar and simple in this first
correct implementation; a heap-based tree builder is a later optimization only
if benchmarks justify it.

## Published omissions and explicit resolutions

| Omission or ambiguity | `v1` resolution |
| --- | --- |
| FFT padding, transform convention, and equal-peak order | Linear convolution, power-of-two zero padding through pinned RustFFT 6.4.1, deterministic lag order described above |
| Exact segment-DP state and ties | Maximum-weight monotone 2D interval chain |
| Meaning of a matrix “division” point | Fix each selected segment's center pair and align intervening rectangles |
| DP boundary initialization and terminal gaps | Global alignment; one extension-free opening per terminal run |
| Equal DP scores | Published project-wide tie order above |
| Short or empty sequences in 6-tuple distance | Common shorter `k` and explicit empty rules |
| Protein ambiguity symbols | Equal canonical mixtures and documented biochemical analogues |
| Nucleotide IUPAC symbols | Equal compatible-base mixtures; ambiguous tuple windows omitted |
| K80 transition/transversion convention | Instantaneous `alpha/beta = 2` |
| Rounded JTT table reversibility | Average the two log-odds orientations |
| Zero-length guide trees | Equal sequence weights |
| Threading and SIMD | Alignment reductions retain deterministic order; the generic FFT primitive uses runtime-dispatched SIMD with a scalar fallback |

These resolutions are part of the versioned contract. They must not be changed
silently; a behavior change requires a new profile version.

## Validation and acceptance evidence

The implementation satisfies the following completion gates:

1. hand-calculated tuple intersections, UPGMA merges, normalized-matrix moments,
   profile scores, extension-free gap runs, and segment chains;
2. direct `O(NM)` correlation versus the radix-2 FFT for exhaustive short signal
   vectors;
3. optimized profile DP versus an independent cubic recurrence for exhaustive
   short score grids and property-generated profiles;
4. independent replay of every returned alignment's symbol consumption and
   objective score;
5. empty, singleton, identical, dissimilar, repetitive, leading-gap,
   internal-gap, and terminal-gap cases for protein, DNA, and RNA;
6. property tests for rectangularity, permutation-stable input row order,
   deterministic repeatability, reverse pair distance, and ungapped recovery;
7. checked-in black-box observations from an official MAFFT 7.526 executable
   invoked as `mafft --retree 1 --maxiterate 0 --thread 1 --quiet` in the
   isolated reference laboratory;
8. exact agreement for deliberately unambiguous trivial cases and documented SP
   and TC thresholds for nontrivial cases, recognizing the published version
   differences listed above;
9. BAliBASE protein subset SP/TC measurements generated independently of this
   implementation's expected outputs;
10. reproducible Criterion time/memory-shape benchmarks plus process-level
    timing against the pinned official executable;
11. formatting, Clippy with warnings denied, all workspace tests, rustdoc with
    warnings denied, dependency license/advisory/source policy, and source-hygiene
    checks.

The crate contributes 37 deterministic tests: 27 component/oracle tests, seven
public integration cases, two property tests with 128 generated inputs each,
and one static reference suite. Radix-2 correlation is checked against direct
correlation for exhaustive short binary signals. The optimized extension-free
profile recurrence is checked against an independently written cubic oracle,
and every progressive merge score is replayed from its returned operations.

The eight-case official MAFFT 7.526 fixture has mean SP `0.733667546`, mean TC
`0.602408143`, and three byte-identical alignments. Its checked-in SHA-256 is
`583c98c19ddbd1ef41b6c95d1f33fa3ef4f75d00577bebd86c5effc3bc8a5215`.
On the predeclared six-family BAliBASE 3 subset, macro SP/TC are
`0.877536523`/`0.557288444` and micro SP/TC are
`0.852642174`/`0.782608696`. Criterion, memory, and matched process results are
recorded in `../benchmarks/mafft-fft-ns-1.md`
and
`../benchmarks/mafft-fft-ns-1-vs-reference.md`.

The MAFFT executable does not expose its internal progressive objective score,
so black-box validation compares valid symbol consumption, exact alignment text
where meaningful, SP, TC, and deterministic output. Internal scores are instead
verified by a separate replay implementation and hand-derived fixtures.

## Known limitations of this profile

- It reproduces the 2002 scientific method, not every undocumented behavior of
  the historical MAFFT 1 executable or the current 7.526 executable.
- It does not include MAFFT 6's 6-mer length correction or the later modified
  UPGMA fragment heuristic.
- It does not include a second progressive pass, iterative refinement,
  consistency scoring, structural information, direction adjustment, or
  add-to-alignment modes.
- Fourier arithmetic is scalar `f64`; SIMD and parallel pairwise distances are
  deferred post-correctness optimizations.
- The published JTT one-PAM table is rounded to five decimal places. The exact
  unrounded historical matrix was obtainable from the authors but is not printed
  in the paper; `v1` uses the printed scientific record.
