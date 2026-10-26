# MAFFT FFT-NS-2

## Reproduction target

Status: complete for `MAFFT-2002-FFT-NS-2-v1`; portfolio acceptance reconfirmed on
2026-08-12.

This document specifies the immutable scientific profile
`MAFFT-2002-FFT-NS-2-v1`.

The target is the two-pass progressive method introduced by Katoh, Misawa,
Kuma, and Miyata in 2002. It is not a reconstruction of the current MAFFT
defaults. The first pass is exactly the project's separately versioned
`MAFFT-2002-FFT-NS-1-v1` profile. A second guide tree is inferred from that
first multiple alignment and the original, ungapped input sequences are then
progressively aligned again.

The production crate contains no MAFFT source, translated source, FFI,
subprocess invocation, or executable dependency. MAFFT 7.526 is used only as
an opaque black-box comparator in an external validation laboratory.

## Primary references

1. Katoh K, Misawa K, Kuma K, Miyata T. *MAFFT: a novel method for rapid
   multiple sequence alignment based on fast Fourier transform*. Nucleic Acids
   Research. 2002;30(14):3059-3066.
   DOI: [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436).
   This is the primary algorithm source. Its FFT-NS-2 section states that the
   input sequences are realigned along a guide tree inferred from the
   FFT-NS-1 alignment.
2. Jones DT, Taylor WR, Thornton JM. *The rapid generation of mutation data
   matrices from protein sequences*. Computer Applications in the Biosciences.
   1992;8(3):275-282.
   DOI: [10.1093/bioinformatics/8.3.275](https://doi.org/10.1093/bioinformatics/8.3.275).
3. Grantham R. *Amino acid difference formula to help explain protein
   evolution*. Science. 1974;185(4154):862-864.
   DOI: [10.1126/science.185.4154.862](https://doi.org/10.1126/science.185.4154.862).
4. Kimura M. *A simple method for estimating evolutionary rates of base
   substitutions through comparative studies of nucleotide sequences*.
   Journal of Molecular Evolution. 1980;16:111-120.
   DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581).
5. Thompson JD, Higgins DG, Gibson TJ. *CLUSTAL W: improving the sensitivity of
   progressive multiple sequence alignment through sequence weighting,
   position-specific gap penalties and weight matrix choice*. Nucleic Acids
   Research. 1994;22(22):4673-4680.
   DOI: [10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
6. Grishin NV. *Estimation of the number of amino acid substitutions per site
   when the substitution rate varies among sites*. Journal of Molecular
   Evolution. 1995;41:675-679.
   DOI: [10.1007/BF00175826](https://doi.org/10.1007/BF00175826).
   This paper states the standard Poisson relation `q = exp(-d)` between the
   unchanged fraction and substitutions per site.
7. Cooley JW, Tukey JW. *An algorithm for the machine calculation of complex
   Fourier series*. Mathematics of Computation. 1965;19(90):297-301.
   DOI: [10.1090/S0025-5718-1965-0178586-1](https://doi.org/10.1090/S0025-5718-1965-0178586-1).

Official historical and operational evidence:

8. The official MAFFT algorithm page describes FFT-NS-2 as FFT-NS-1 followed
   by recomputation of the guide tree from the first MSA and a second
   progressive alignment:
   <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.
9. The official MAFFT change log records that an early second phase used the
   number of substitutions with Poisson correction and that a later release
   replaced it with a pairwise alignment score:
   <https://mafft.cbrc.jp/alignment/software/changelog.html>.
10. The official manual identifies `mafft --retree 2 --maxiterate 0 INPUT` as
    the operational FFT-NS-2 selection:
    <https://mafft.cbrc.jp/alignment/software/manual/manual.html>.

Items 8-10 resolve version history and black-box invocation only. They do not
override the 2002 paper's scoring parameters.

## Objective

FFT-NS-2 reduces progressive-alignment error caused by the very rough shared
tuple tree used by FFT-NS-1. Its characteristic pipeline is:

1. construct the rough shared-tuple distance matrix;
2. construct the first UPGMA guide tree;
3. perform the complete FFT-NS-1 progressive alignment;
4. estimate a new evolutionary distance matrix from that MSA;
5. construct a second UPGMA guide tree and new sequence weights;
6. discard the first alignment's inserted gaps as alignment state;
7. progressively realign the original sequences along the second tree.

Step 7 always runs, even when both tree topologies are equal. Stopping after
step 3 would silently substitute FFT-NS-1. Reusing the gapped rows from step 3
as fixed profiles would also be a different method: the paper says the input
sequences are realigned.

FFT-NS-2 does not perform tree-dependent bipartition refinement and does not
accept or reject individual realignments by an objective score. Those are
FFT-NS-i behaviors and are outside this phase.

## First progressive pass

The complete first pass is the public `MAFFT-2002-FFT-NS-1-v1` algorithm
specified in [mafft-fft-ns-1.md](mafft-fft-ns-1.md). In summary:

- protein sequences use compressed Dayhoff-class tuples and nucleotide
  sequences use canonical-base tuples;
- the initial distance is based on the multiset count of shared tuples;
- deterministic size-weighted UPGMA supplies the first guide tree;
- CLUSTALW root-to-leaf path weights supply sequence weights;
- the pinned runtime-dispatched FFT primitive locates candidate homologous regions;
- segment chaining partitions profile/profile dynamic programming;
- the normalized JTT 200-PAM or K80 200-PAM matrix and the publication's
  position-specific, extension-free gap recurrence score profile merges.

The first-pass result exposed by FFT-NS-2 must be bit-for-bit identical to a
direct call to `MafftFftNs1` on the same typed records.

## Distance inferred from the first alignment

### Published statement and historical evidence

Katoh et al. state only that the second guide tree is inferred from the
FFT-NS-1 alignment. They do not publish an equation, a treatment of gaps, a
saturation rule, or behavior for ambiguity symbols and empty sequences.

The official change log supplies important version evidence: the original
second phase used a number of substitutions with Poisson correction; a later
release changed this to a pairwise alignment score. Because this profile
targets 2002, it uses the former family of distance estimates rather than the
later score-based behavior.

### Explicit reproducible resolution

Let `X_i(c)` be row `i` at column `c` of the first MSA. For a pair `(i,j)`, the
set of comparable columns uses pairwise deletion:

```text
C_ij = { c | X_i(c) is a residue and X_j(c) is a residue }
n_ij = |C_ij|
m_ij = sum over c in C_ij of I[X_i(c) = X_j(c)].
```

The unchanged fraction is

```text
q_ij = m_ij / n_ij                         when n_ij > 0
q_ij = 1                                  when both rows are empty
q_ij = 0                                  otherwise.
```

The Poisson-corrected evolutionary distance is

```text
epsilon = 10^-12
D2_ii   = 0
D2_ij   = -ln(max(q_ij, epsilon)), i != j.
```

The Grishin relation `q = exp(-d)` gives the inverse `d = -ln(q)`. The floor
makes completely non-identical or non-overlapping pairs finite, with a maximum
distance of `-ln(10^-12)`. This saturation value is an implementation
resolution, not a parameter printed by Katoh et al.

Pairwise deletion is used because the Poisson model estimates substitutions,
not insertion/deletion events. A residue-gap or gap-gap column therefore does
not increment either `n_ij` or `m_ij`. A non-empty pair with no overlapping
residues is maximally distant. Two empty rows are identical by definition.

Typed symbols are identical only when their encoded symbols are equal. For
example, `N` does not count as an identity with `A`, and protein `B` does not
count as an identity with `D`. The 2002 paper does not specify ambiguity-code
probabilities for guide-tree inference; exact equality avoids silently
crediting an uncertain match.

The result exposes the unscaled finite matrix. Distances are not clamped to
`[0,1]`: they represent corrected expected substitutions per site. The custom
matrix type validates square shape, symmetry, a zero diagonal, nonnegative
finite values, and checked allocation sizes.

## Second guide tree and weights

The second matrix is clustered with the same deterministic, size-weighted
UPGMA definition used by this profile's first pass. If clusters `A` and `B`
merge, their distance to an active cluster `C` is

```text
D_(A union B),C = (|A| D_A,C + |B| D_B,C) / (|A| + |B|).
```

The merge height is half the selected inter-cluster distance. Branch lengths
are parent height minus child height. The second pass recomputes normalized
CLUSTALW path weights from this new tree; it does not retain the first-pass
weights for scoring.

The 2002 paper names UPGMA and does not describe the mixed-linkage procedure
documented for later MAFFT versions. Accordingly, this versioned profile uses
ordinary size-weighted UPGMA in both passes.

## Second progressive pass

The second pass starts from one ungapped leaf profile per original
`SequenceRecord`. It uses the same alphabet-specific scorer and the same
published constants as the first pass:

```text
tuple length                         6
FFT peaks inspected                 20
homology window                     30
homology threshold per column       0.7
maximum retained segment length     150
gap opening S_op                    2.4
random-match offset S_a             0.06
```

For normalized profile weights `w_n` and `w_m`, the profile homology score is

```text
H(i,j) = sum_n sum_m w_n w_m M_hat[A(n,i), B(m,j)],
```

where residue-gap and gap-gap pairs contribute zero. The normalized matrix is

```text
average1 = sum_a f_a M_aa
average2 = sum_a sum_b f_a f_b M_ab
M_hat_ab = (M_ab - average2) / (average1 - average2) + S_a.
```

The extension-free profile recurrence is

```text
P(i,j) = H(i,j) + max(
    P(i-1,j-1),
    max_x(P(x,j-1) - G_1(i,x)),
    max_y(P(i-1,y) - G_2(j,y))
).
```

The implementation uses the same separable-prefix optimization, traceback,
FFT correlation, segment selection, and score replay already specified and
tested for FFT-NS-1. Sharing these project-owned primitives does not collapse
the two pipelines: FFT-NS-2 retains and exposes both passes and always executes
the second one.

## Tie-breaking

The historical executable's ties are not published. This profile defines a
total deterministic order:

1. pairwise distances are evaluated in increasing input-index order;
2. UPGMA selects the smallest distance, then the lexicographically smallest
   pair of minimum descendant indices;
3. the child containing the lower minimum input index is the left child;
4. zero-distance branches use equal positive weights if all path weights
   collapse to zero;
5. FFT lag, homologous segment, profile-DP, and traceback ties use the complete
   ordering specified by `MAFFT-2002-FFT-NS-1-v1`;
6. returned rows always follow original input order.

No random number, hash iteration order, thread scheduling, locale, clock, or
external program can affect a library result.

## Pseudocode

```text
FFT_NS_2(records):
    validate non-empty collection and unique identifiers

    first <- FFT_NS_1(records)

    for every pair of first.alignment rows i,j:
        comparable <- 0
        identical  <- 0
        for each alignment column c:
            if both cells are residues:
                comparable += 1
                identical += I[cell_i == cell_j]

        if i == j or both rows are empty:
            D2[i,j] <- 0
        else if comparable == 0:
            D2[i,j] <- -ln(epsilon)
        else:
            q <- identical / comparable
            D2[i,j] <- -ln(max(q, epsilon))

    tree2    <- deterministic size-weighted UPGMA(D2)
    weights2 <- CLUSTALW_path_weights(tree2)
    scorer   <- normalized JTT200 or K80-200 matrix

    second <- PROGRESSIVE_FFT_PASS(
        original_ungapped_records,
        tree2,
        weights2,
        scorer,
        paper_2002_parameters
    )

    independently replay every second-pass profile merge score
    validate rectangularity, input order, and exact symbol consumption
    return first, D2, tree2, weights2, second, and memory diagnostics
```

## Data structures

- the complete immutable `FftNs1Alignment` audit record for pass one;
- a dense symmetric `FftNs2DistanceMatrix` of Poisson-corrected values;
- a second validated binary `GuideTree`;
- a second normalized sequence-weight vector;
- the same internal weighted profile, FFT, segment-chain, and compact
  traceback structures as FFT-NS-1;
- separate first- and second-pass merge diagnostics;
- a final `MultipleAlignment` in original input order;
- explicit major-buffer memory accounting for both passes.

## Complexity

Let `K` be the sequence count, `L` a representative ungapped length, `W1` the
first alignment width, and `P(K,L)` the time of one FFT-assisted progressive
pass as detailed in the FFT-NS-1 specification.

- first rough-distance and progressive pass: the FFT-NS-1 costs;
- aligned Poisson distances: `O(K^2 W1)` time and `O(K^2)` memory;
- second scalar UPGMA: `O(K^3)` time and `O(K^2)` memory;
- second progressive pass: another `P(K,L)`;
- returned first and final alignments: `O(K(W1 + W2))` typed cells in this
  audit-oriented API.

The overall time is

```text
O(FFT_NS_1 + K^2 W1 + K^3 + second_progressive_pass).
```

The FFT partitions often reduce profile-DP work, but either progressive pass
can still require a full quadratic profile rectangle. The initial scalar tree
implementation deliberately favors auditability over an asymptotically faster
heap or nearest-neighbor structure.

## Components and independent checks

Implementation is divided into independently testable components:

1. the already completed FFT-NS-1 pass;
2. aligned-column identity counts;
3. Poisson correction and saturation;
4. finite evolutionary-distance matrix validation;
5. deterministic second UPGMA tree;
6. second path-weight calculation;
7. reusable deterministic FFT progressive pass;
8. final alignment and score replay;
9. phase-aware memory accounting.

Tests must not obtain expected distance values, trees, paths, or scores by
calling the same production routine under test. Hand calculations and simple
cubic/direct oracles are used where applicable.

## Published omissions and explicit resolutions

| Omission | Resolution in this profile |
|---|---|
| Exact second-distance equation | Historical Poisson correction, `d = -ln(q)`. |
| Columns containing gaps | Pairwise deletion. |
| No comparable residues | Maximum finite distance, except two empty rows have distance zero. |
| Infinite distance at zero identity | Identity floor `10^-12`. |
| Ambiguous symbol identity | Exact typed-symbol equality only. |
| Tree linkage details | Ordinary size-weighted UPGMA from the 2002 paper. |
| Whether first-pass gaps remain fixed | They do not; original sequences are realigned from leaf profiles. |
| Whether to skip an unchanged second tree | Never skip the second progressive pass. |
| Tree and DP ties | Stable total orders documented above. |
| Empty and singleton collections | Empty collection is rejected; singleton is returned through two auditable zero-merge passes. |

These decisions are part of the stable profile identifier. A future profile
that reproduces later score-based distances or mixed linkage requires a new
identifier and separate evidence.

## Validation and acceptance criteria

FFT-NS-2 can be marked complete only when all of the following hold:

1. this scientific specification contains references, equations, original
   pseudocode, data structures, complexity, parameters, tie rules, omissions,
   limitations, and acceptance evidence;
2. the first-pass audit is bit-for-bit equal to a direct FFT-NS-1 run;
3. hand-computed aligned distances cover identities, substitutions, internal
   and terminal gaps, ambiguity symbols, empty rows, no overlap, and complete
   divergence;
4. an independent test oracle recomputes every matrix entry directly from the
   exposed first alignment;
5. the second tree is built from the second matrix and the second weights are
   independently normalized and finite;
6. both identical-tree and changed-tree examples prove that the second pass
   executes;
7. every returned row consumes exactly its corresponding original typed
   sequence and the alignment contains no all-gap column;
8. unit and property tests cover protein, DNA, RNA, empty, singleton,
   identical, completely different, repetitive, ambiguity-rich, leading-gap,
   internal-gap, and terminal-gap inputs;
9. output, diagnostics, distances, trees, weights, and scores are deterministic;
10. each profile merge score has independent operation replay evidence;
11. an inert fixture generated from the pinned official MAFFT 7.526 binary is
    checked for exact trivial cases and SP/TC agreement on nontrivial cases;
12. a fixed external BAliBASE 3 subset records SP and TC metrics without
    vendoring the scientific dataset;
13. Criterion reports length scaling, sequence-count scaling, a long FFT path,
    globin fragments, and direct FFT-NS-1 versus FFT-NS-2 cost;
14. matched fresh-process measurements compare the Rust worker and official
    MAFFT 7.526 with binary, package, harness, and input hashes;
15. formatting, all-target Clippy, workspace tests, release tests, Rustdoc, and
    dependency license/source/advisory checks all pass;
16. the public API, example, benchmark protocol, profile differences, and
    limitations are documented in English;
17. the implementation contains no unfinished path and no runtime dependency
    on a reference aligner.

## Black-box reference boundary

The isolated comparator is the official Linux amd64 MAFFT 7.526 package,
identified by package SHA-256 in the generated fixture. The command is:

```text
mafft --retree 2 --maxiterate 0 --thread 1 --inputorder --quiet INPUT
```

The package is treated as opaque. No source is downloaded or inspected. Its
output is inert test data; it is neither executed nor read by product code.

Exact equality is required only for cases whose optimal alignment is
unambiguous, such as identical inputs. Other cases use SP and TC because MAFFT
7.526 differs materially from the 2002 target:

- BLOSUM62 instead of JTT 200 PAM by default;
- later gap and offset parameters;
- post-2002 shared-tuple length correction;
- mixed-linkage rather than the paper's ordinary UPGMA;
- later pairwise-score distance rather than the historical Poisson distance;
- numerous later FFT and profile-alignment revisions.

## Known limitations of this profile

- It reconstructs the 2002 method, not byte-for-byte behavior of any released
  MAFFT executable.
- Pairwise deletion, ambiguity handling, zero-identity saturation, and empty
  row behavior are explicit resolutions of publication omissions.
- The JTT 1-PAM source table is printed to five decimal places; the first-pass
  profile renormalizes its rows before exponentiation.
- Both tree construction and all-pairs aligned-distance calculation are scalar.
- FFT and profile DP are scalar and single-threaded; SIMD and parallelism are
  deferred until the scalar profile passes every acceptance gate.
- UPGMA is cubic in sequence count in this first auditable implementation.
- The audit result retains the first and final MSA, increasing memory compared
  with a streaming executable that discards the first MSA.
- FFT partitioning is heuristic, and worst-case profile DP remains quadratic.
- The checked-in BAliBASE evidence uses a fixed representative subset, not the
  complete benchmark distribution.
