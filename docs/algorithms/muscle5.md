# MUSCLE 5 Multiple Sequence Alignment

## Reproduced profile

This phase targets `MUSCLE-5.1.0-Nature-2022-default-v2`, an independent,
deterministic reconstruction of the MUSCLE 5 method published by Edgar in 2022.
Status: implemented. Version 1 passed the 2026-08-04 release gate, and its
validation evidence and benchmark records stay pinned to it. Version 2
(2026-10-22) changes only how refinement scores the current alignment; see the
[version history](#version-history).
The pinned profile includes all sequence-alignment modes described as the main
MUSCLE 5 pipelines:

- `Align`, implemented as Parallel Perturbed ProbCons (PPP): pair-HMM posterior
  probabilities, two probabilistic-consistency transformations, expected-error
  UPGMA, posterior-sum progressive profile alignment, and 100 randomized
  partition-refinement attempts;
- `Super5`: decreasing-length input ordering, UCLUST-EE redundancy filtering,
  coarse clustering and subdivision, PPP cluster alignment, consensus-tree
  construction, PPP-pa/PPP-ps cluster merging, and redundant-row insertion;
- stratified ensembles, with four guide-tree permutations for each perturbation
  seed; and
- diversified ensembles, with a new perturbation seed per replicate while the
  four guide-tree permutations cycle in published order.

The implementation accepts typed protein, DNA, and RNA sequences. It is safe
Rust, performs no external I/O, and has no runtime dependency on MUSCLE,
ProbCons, Python, or another aligner. It does not call, bind, translate, inspect,
or derive from MUSCLE or ProbCons source code. Correct scalar kernels were
validated first; independent sequence-pair work is now dispatched through
deterministic scoped threads while every per-pair recurrence remains scalar.

This implementation targets only the versioned MUSCLE 5 pipelines; it is not a
generic progressive alignment renamed MUSCLE or an
approximation that omits consistency, expected-accuracy scoring, randomized
refinement, guide-tree permutations, or the characteristic Super5 clustering
pipeline.

## Scientific sources

Primary MUSCLE 5 source:

1. Robert C. Edgar. "MUSCLE v5 enables improved estimates of phylogenetic tree
   confidence by ensemble bootstrapping." bioRxiv, 2021. DOI:
   [10.1101/2021.06.20.449169](https://doi.org/10.1101/2021.06.20.449169).
   This preprint introduces PPP, Super5, perturbed HMMs, guide-tree
   permutations, and alignment ensembles.
2. Robert C. Edgar. "MUSCLE5: High-accuracy alignment ensembles enable unbiased
   assessments of sequence homology and phylogeny." *Nature Communications*
   13:6968, 2022. DOI:
   [10.1038/s41467-022-34630-w](https://doi.org/10.1038/s41467-022-34630-w).
   The published article and its 22-page supplementary information are the
   normative algorithm description for this profile.

Original component methods:

3. Chuong B. Do, Mahathi S. P. Mahabhashyam, Michael Brudno, and Serafim
   Batzoglou. "ProbCons: Probabilistic consistency-based multiple sequence
   alignment." *Genome Research* 15(2), 330-340, 2005. DOI:
   [10.1101/gr.2821705](https://doi.org/10.1101/gr.2821705). This defines
   posterior expected-accuracy alignment, probabilistic consistency, the
   expected-accuracy guide tree, progressive profile scoring, and randomized
   iterative refinement on which PPP is based.
4. Ian Holmes and Richard Durbin. "Dynamic programming alignment accuracy."
   *Journal of Computational Biology* 5(3), 493-504, 1998. DOI:
   [10.1089/cmb.1998.5.493](https://doi.org/10.1089/cmb.1998.5.493). This
   establishes posterior decoding and maximum expected-accuracy alignment.
5. Robert C. Edgar. "Search and clustering orders of magnitude faster than
   BLAST." *Bioinformatics* 26(19), 2460-2461, 2010. DOI:
   [10.1093/bioinformatics/btq461](https://doi.org/10.1093/bioinformatics/btq461).
   This defines U-sort/UCLUST, which MUSCLE 5 adapts with expected error as the
   acceptance distance.
6. Steven Henikoff and Jorja G. Henikoff. "Amino acid substitution matrices
   from protein blocks." *Proceedings of the National Academy of Sciences*
   89(22), 10915-10919, 1992. DOI:
   [10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915).
   This is the original source for BLOSUM62.
7. Yi-Kuo Yu, John C. Wootton, and Stephen F. Altschul. "The compositional
   adjustment of amino acid substitution matrices." *Proceedings of the
   National Academy of Sciences* 100(26), 15688-15693, 2003. DOI:
   [10.1073/pnas.2533904100](https://doi.org/10.1073/pnas.2533904100). The NCBI
   public-domain matrix-frequency record accompanying this work supplies the
   full-precision BLOSUM62 joint and marginal probabilities named, but not
   printed, by the MUSCLE 5 paper.

Version and manual source:

8. Robert C. Edgar, official *MUSCLE v5 manual*: the `align`, `super5`, HMM,
   HMM-perturbation, guide-tree-permutation, stratified-ensemble, and
   diversified-ensemble pages. These pages specify the documented command
   defaults: two consistency iterations, 100 refinement iterations, perturbation
   amplitude 0.25, four guide-tree permutations, 16 stratified replicates, and
   100 diversified replicates.

Independent black-box reference:

9. The unmodified official `muscle5.1.linux_intel64` release executable is run
   only in the isolated sibling validation laboratory. It reports
   `muscle 5.1.linux64 [12f0e2]`, built 2022-01-13. Its SHA-256 is
   `2bd82928a5cc72b5086b60928a3a820a0ce066ca9fbc0fecd1543eebca981e76`.
   Its source is not downloaded or inspected. Reviewed observations may be
   stored as immutable test fixtures, but Cargo never invokes the executable.

Scientific validation datasets:

10. Julie D. Thompson et al. "BAliBASE 3.0: latest developments of the multiple
    sequence alignment benchmark." *Proteins* 61(1), 127-136, 2005. DOI:
    [10.1002/prot.20527](https://doi.org/10.1002/prot.20527).
11. Robert C. Edgar. "Quality measures for protein alignment benchmarks."
    *Nucleic Acids Research* 38(7), 2145-2153, 2010. DOI:
    [10.1093/nar/gkp1196](https://doi.org/10.1093/nar/gkp1196). This defines the
    PREFAB Q-score family used by the MUSCLE publications.
12. Van Walle et al., SABmark, and Sievers et al., HomFam, are reserved for the
    shared cross-family validation suite after stable redistribution provenance
    is pinned. Their omission is validation scope, not an algorithmic shortcut.

## Version boundary

The scientific description is the 2022 Nature profile. The black-box executable
is exactly release 5.1.0, the closest official release artifact associated with
that publication. Later MUSCLE 5 releases add behavior, including structural
alignment, that is outside this profile and outside the current project scope.

The public API restores stable input row order. Release 5.1.0 normally emits
guide-tree order. Black-box comparison therefore matches rows by unique
identifier before comparing residue-pair sets, SP, TC, and aligned strings.

## Inputs, outputs, and invariants

Input is a nonempty ordered collection of uniquely named, typed, ungapped
sequences. Individual sequences may be empty. Every input in one call has the
same compile-time alphabet: `Protein`, `Dna`, or `Rna`.

Every returned replicate contains:

- a validated MSA in original input order;
- mode, perturbation seed, and guide-tree permutation;
- the pairwise expected-error matrix and deterministic guide tree;
- counts and checksums for raw and consistency-transformed posteriors;
- progressive merge and refinement evidence;
- Super5 cluster, representative, and sampled-pair evidence when applicable;
- an independently recomputable final posterior objective; and
- explicit algorithm-owned resident and peak-memory estimates.

Removing gaps from output row `i` must reproduce input sequence `i` exactly.
Identifiers are unchanged, all rows have equal width, no column is all gaps, and
equal parameters and seeds produce byte-identical output and metadata.

## Pair-HMM

For sequences `x=x_1...x_m` and `y=y_1...y_n`, the five emitting states are:

```text
M   emits (x_i,y_j)
Ix  emits x_i and a gap in y       short gap
Iy  emits a gap in x and y_j       short gap
Jx  emits x_i and a gap in y       long gap
Jy  emits a gap in x and y_j       long gap
```

The symmetric published default transitions are:

```text
S -> M=.62, Ix=.02, Iy=.02, Jx=.18, Jy=.18  (printed weights)
M -> M=.96, Ix=.012, Iy=.012, Jx=.008, Jy=.008
I -> I=.40, M=.60
J -> J=.90, M=.10
M -> E=.62, Ix/Iy -> E=.02, Jx/Jy -> E=.18
```

Transitions between opposite insertion states are absent. Start and end weights
are symmetric as stated by the paper. The rounded start weights printed in
Supplementary Figure S1 sum to 1.02, despite the Methods requirement that every
exclusive set sum to one. This profile divides that complete set by 1.02,
retaining the printed ratios exactly. A common start/end scale cancels from all
posteriors, so this explicit normalization resolves the contradiction without
changing an alignment probability. The other printed outgoing sets already sum
to one.

Let `q(a,b)` be a match emission and `p(a)` an insertion emission. With
state-specific coordinate increments `Delta(s)`, the log-domain form of the
forward recurrence is

```text
F_s(i,j) = log_emit_s(i,j)
           + logsumexp_r(F_r(i-Delta_i(s),j-Delta_j(s)) + log t_rs).
```

The start weight replaces the predecessor sum at the first emitted state. The
backward recurrence sums over the next emitted state and ends with the published
state-to-end weight. The log partition function is

```text
log Z = logsumexp_s(F_s(m,n) + log t_sE).
```

The posterior match probability is

```text
P_xy(i,j) = exp(F_M(i,j) + B_M(i,j) - log Z).
```

The implementation has two mathematically equivalent numeric evaluators. It
first scores one concrete valid state path. If that path has log probability at
least `-600`, then `Z >= exp(-600)` and ordinary `f64` probability-domain
forward/backward is safe for every contribution that can affect a posterior at
machine precision. Posterior evaluation uses `(F_M/Z) * B_M`, avoiding a
premature product underflow. A forward or backward value small enough to
underflow is then less than `exp(-145)` relative to the proven partition lower
bound. For harder or longer pairs, the implementation selects the full
log-domain recurrence automatically. Its stable two-term operation is
`max(a,b) + log1p(exp(min(a,b)-max(a,b)))`, with `-infinity` for unreachable
states. A unit oracle evaluates the same pairs in both domains and requires
partition and posterior agreement within `2e-12`.

Both evaluators expand the published five-state topology explicitly, so absent
transitions are not tested inside every cell. Neither evaluator substitutes
Viterbi scores for posterior probabilities. Final posterior roundoff is clamped
to `[0,1]` before sparsification.

### Protein emissions

Canonical protein `q(a,b)` is the full-precision symmetric BLOSUM62 joint
frequency table named by the MUSCLE 5 paper. `p(a)=sum_b q(a,b)` uses the
corresponding NCBI marginal frequencies in canonical order
`ARNDCQEGHILKMFPSTWYV`. The constants are scientific data, not code from an
alignment implementation.

Protein ambiguities are probability mixtures: `B={D,N}`, `Z={E,Q}`,
`J={I,L}`, `X` is the BLOSUM62 marginal mixture over 20 residues, `U={C}`, and
`O={K}`. For mixtures `u` and `v`, emissions are

```text
q(u,v) = sum_a sum_b u_a v_b q(a,b)
p(u)   = sum_a u_a p(a).
```

### Nucleotide emissions and publication ambiguity

The paper states that MUSCLE 5 adds nucleotide support but describes all match
emissions as BLOSUM62 and publishes no DNA/RNA emission table. BLOSUM62 cannot
define nucleotide emissions. The official manual does not resolve this omission,
and the reference CLI does not expose posterior probabilities. This profile
therefore makes the missing choice explicit rather than silently claiming an
author value:

```text
canonical p(A)=p(C)=p(G)=p(T/U)=0.25
q(a,a)=0.22 for each canonical base
q(a,b)=0.01 for each ordered a != b
```

The matrix sums to one and every row/column sums to 0.25. It is a symmetric,
identity-biased coupled model with a transition/transversion ratio of one. DNA
and RNA share probabilities while preserving `T` and `U`. IUB symbols distribute
uniformly over represented canonical bases. This extension is versioned as
`symmetric-nucleotide-joint-v1`; black-box alignment agreement is reported
separately and never presented as proof that the hidden reference matrix is
identical.

## HMM perturbations

For nonzero integer seed `s`, every probability is independently perturbed:

```text
P_k* = (1 + alpha delta_k) P_k,       delta_k ~ Uniform[-1,1]
P_k' = P_k* / sum_(i in same event set) P_i*,   alpha=.25.
```

The exclusive sets are start/end transitions, transitions from `M`, from each
short-gap state, from each long-gap state, match emissions, and insertion
emissions. Symmetry-equivalent probabilities share one perturbation so sequence
swap and reversal symmetry remain intact. Seed zero means exactly no HMM
perturbation. The paper does not identify its pseudo-random generator. This
profile uses SplitMix64 with a published stage-domain derivation; its output is
deterministic across machines but is not expected to duplicate the private
random stream of the reference executable.

## Sparse posterior representation

A posterior matrix stores rows of `(column, probability)` entries in ascending
column order. Raw values below `0.01`, the cutoff published by ProbCons, are
discarded after posterior construction and after each consistency round. Reverse
orientation is represented by transpose, not a second HMM calculation. Identity
posteriors used during consistency are exact diagonal matrices.

The checksum is computed over sequence indices, dimensions, coordinates, and
IEEE-754 probability bits in stable order. It is evidence, not an expected test
oracle.

Unordered pair jobs are enumerated lexicographically. Workloads of at least
65,536 pair-lattice cells use up to 24 available scoped workers, divided into
fixed contiguous job chunks. Results are joined and installed in original job
order. Scheduling therefore cannot change floating addition order within a
posterior or change checksums. The sequential and parallel collectors are
required to agree bit-for-bit. The memory estimate reports the sum of concurrent
worker high-water scratch, not merely one worker's matrix.

Each worker owns one reusable Forward, Backward, and dense-posterior workspace.
The buffers are cleared and reused for subsequent jobs in that worker's fixed
chunk; no probability, traversal order, cutoff, or traceback decision changes.
Consistency workers similarly reuse one dense accumulator. This removes allocator
traffic while the reported concurrent high-water remains the sum of the largest
logical workspace assigned to each active worker.

## Expected-accuracy pair alignment and expected error

For a candidate alignment `A`, expected accuracy is

```text
EA(A|x,y) = sum_((i,j) in matched(A)) P_xy(i,j).
```

Maximum expected accuracy is found with gap score zero:

```text
E(i,j) = max(
    E(i-1,j-1) + P_xy(i,j),
    E(i-1,j),
    E(i,j-1)).
```

Ties prefer a residue-residue match, then consuming from the first sequence,
then consuming from the second. The expected error used as a distance is the
published mean posterior error over alignment columns:

```text
EE(x,y) = 1 - EA(A*|x,y) / |A*|.
```

The empty/empty distance is zero; empty/nonempty is one. The returned pair path
is independently scoreable from the stored posterior.

## Probabilistic consistency

For sequence set `S` of size `N`, one consistency transformation is

```text
P'_xy(i,j) = (1/N) sum_(z in S) sum_k P_xz(i,k) P_zy(k,j).
```

When `z=x` or `z=y`, the corresponding matrix is an identity, so the original
`P_xy` contributes twice. Sparse row products accumulate in coordinate order.
Values are divided by `N`, clamped only for floating roundoff, filtered at 0.01,
and stored in stable row/column order. The published default applies two rounds.
Independent output-pair jobs use the same fixed chunking and ordered join as raw
posteriors; intermediates and sparse entries are still accumulated sequentially in
their documented order, so parallelism is bit-identical to one-worker execution.

## Expected-error UPGMA guide tree

Pairwise `EE` is the initial cluster distance. At each step the closest active
clusters merge. For every other cluster `C`, arithmetic average linkage is

```text
d(A union B,C) = (|A| d(A,C) + |B| d(B,C)) / (|A|+|B|).
```

This size-weighted recurrence is the standard UPGMA form and equals the mean of
all leaf-pair expected errors. The ProbCons prose gives the equal-size special
case. A merge tie is resolved by the lexicographically smallest pair of sorted
leaf-index vectors. Children are ordered by their minimum input index. Parent
height is half the merge distance; negative limbs from inversions are clamped to
zero and recorded.

## Guide-tree permutations

For the complete guide tree, its leaves are divided into three connected
parts `a,b,c`: the unrooted subtree is first split as close as possible to a
`1:2` leaf ratio without cutting an existing descendant group, then the larger
side is split as close as possible to `1:1`. Stable leaf-vector order breaks
equal-size choices. The four published resolutions are:

```text
none = original rooted topology
abc  = ((a,b),c)
acb  = ((a,c),b)
bca  = ((b,c),a)
```

Parts and child nodes remain internally unchanged. Subtrees with fewer than
three leaves are unchanged. Permutations alter progressive order, never the
posterior probabilities.

## Progressive PPP profile alignment

A working profile retains stable input row indices, aligned typed cells, and for
each non-gap cell its original one-based sequence position. For profiles `X,Y`,
the score for aligning columns `i,j` is the published PPP-pa sum

```text
M_ij = sum_(x in X) sum_(y in Y) P_xy(pos_x(i),pos_y(j)),
```

where any gap cell contributes zero. A maximum-score global DP uses this match
score and zero gap score. Ties prefer matching profile columns, then consuming
from `X`, then consuming from `Y`. Traceback inserts profile-wide gaps, removes
no existing column, and preserves internal column order.

Original-position-to-profile-column maps turn each sparse posterior entry directly
into one column-score update; absent posterior cells are never probed. The DP keeps
two rolling `f64` score rows and one byte of traceback per lattice cell rather than
a dense `f64` score lattice. The optimized build uses the DP optimum as the
candidate cross-partition score because both are the same PPP-pa sum; checked test
builds replay it independently before any refinement decision.

The guide tree is traversed postorder. A leaf becomes a one-row profile; each
branch aligns its two child profiles. Final rows are reordered to stable input
order and all-gap columns, which should be impossible by construction, are
rejected.

## Randomized partition refinement

The published default makes 100 attempts. At attempt `r`, a deterministic random
bit partitions rows into two nonempty groups; if one side is empty, the lowest
stable index is moved to it. All-gap columns are removed independently from each
subprofile, and the two profiles are realigned using PPP-pa. Because internal
columns are fixed and only cross-profile posterior terms change, the candidate
must not decrease the cross-partition posterior objective.

The current MSA, restricted to the two subprofiles, is itself one path through
the PPP-pa lattice. Its cross-partition objective is therefore computed by
replaying that path over the same column scores the DP uses, adding them in path
order exactly as the DP does. Rounded floating-point addition is monotone, so the
DP optimum can never fall below this value, and a candidate identical to the
current MSA scores exactly the same. The candidate is accepted only when its
fixed-precision objective is strictly larger; equality retains the current MSA,
and a smaller value is an internal-invariant error. Every attempt records both
scores, partition, candidate width, and decision.

The paper does not publish the RNG or its stream interaction with HMM
perturbation. This profile derives independent SplitMix64 streams for HMM,
refinement, Super5 clustering, and PPP-ps sampling, preventing implementation
order from changing scientific output.

A rejected partition is memoized while the current profile is unchanged because
the posterior library, extracted subprofiles, DP path, and score are then a pure
function of that partition. The complete 100-draw schedule and evidence remain
visible. Any accepted candidate changes the profile and immediately clears the
cache; no attempt is dropped or reordered.

## Super5

Super5 follows the published divide-and-conquer sequence:

1. Sort inputs by decreasing ungapped length, then stable input index.
2. Run UCLUST-EE redundancy filtering with threshold `EE < 0.01`; a sequence is
   deferred after it matches a representative. Candidate representatives are
   ranked by shared canonical `k`-mer counts and stable index. At most 16 failed
   posterior/EE comparisons are made before starting a new representative.
3. Cluster representatives with `EE < 0.30` by the same U-sort rule.
4. Subdivide a cluster larger than 500 with `EE < 0.10`.
5. Randomly and deterministically subdivide any remaining group larger than 500.
6. Align every coarse cluster with PPP.
7. Construct its consensus from the most frequent non-gap symbol per column;
   canonical alphabet order breaks a residue tie, and a column is deleted when
   gap is the strict or tied majority.
8. Compute all consensus-pair EE values and a deliberately biased UPGMA tree.
   The publication says "biased" but gives no equation. This profile uses
   cluster-size-weighted average linkage, exposing the sizes and every merge.
9. Merge cluster MSAs by PPP-pa when `|X||Y| <= 2000`; otherwise PPP-ps samples
   exactly 2000 distinct cross-profile sequence pairs without replacement and
   applies the same column equation to that subset.
10. Reinsert deferred redundant sequences through the stored representative
    pair path, transitively if needed. A representative gap propagates to all
    attached rows; an inserted residue creates a new MSA column only at its
    anchored representative interval.

The Super5 expected-error matrix and guide-tree leaves describe coarse-cluster
consensuses, not original input sequences. Cross-cluster PPP-pa/PPP-ps obtains raw
pair posteriors directly for the selected original sequence pairs. A global
ProbCons transformation cannot be defined when Super5 deliberately omits most
pair matrices; the paper gives no alternative restricted-consistency equation, so
these cross matrices are explicitly untransformed rather than silently inventing
one.

Consequently, the final Super5 objective covers the posterior pairs actually
available to the divide-and-conquer pipeline: intra-cluster pairs, selected
cross-cluster PPP-pa/PPP-ps pairs, and representative/deferred pairs. Public
`evaluated_objective_pair_count` and `possible_objective_pair_count` distinguish
that audited subset from all unordered original pairs. Align mode always evaluates
all possible pairs.

Within one run, a pair posterior depends only on its two sequences, in order,
and on the run's HMM and cutoff. The UCLUST comparisons of steps 2 to 4 keep the
posteriors they compute, and the later stages, the cluster PPP libraries, the
consensus library and the PPP-pa/PPP-ps cross pairs, take an ordered pair from
that store when it is there and compute it otherwise. A reverse pair is still
computed: its pair HMM runs in another order, so its posterior need not equal
the transpose bit for bit. Results, checksums and evidence are those of
computing every pair (decision 035).

The UCLUST paper and MUSCLE supplement omit the exact nucleotide/protein `k`,
hash collision policy, candidate tie order, biased-UPGMA formula, random
subdivision stream, PPP-ps sample count when the threshold is crossed, and
transitive insertion tie behavior. The rules above are explicit resolutions,
not claims about hidden executable internals. Protein uses canonical 3-mers;
nucleotides use canonical 6-mers; an ambiguity breaks a word.

## Ensembles

A stratified ensemble with seed count `K` emits, for every seed
`s=0..K-1`, permutations `none,abc,acb,bca` in that order. The documented
default `K=4` therefore has 16 replicates named `none.s`, `abc.s`, `acb.s`, and
`bca.s`.

A diversified ensemble of `R` replicates emits seed `r` with permutation
`[none,abc,acb,bca][r mod 4]`. The documented default is `R=100`. Replicates are
returned in generation order, and their alignments are independently owned.

Release 5.1.0 has a CLI edge case where `-stratified -replicates 1` writes one
plain AFA in a probe. The typed API follows the published ensemble definition
and therefore returns four permutations for one stratified seed.

## Fixed decisions and numerical ties

Posterior probabilities remain `f64`; decisions that compare complete alignment
objectives quantize once at `1e-9` units. Quantizing alone does not make such a
comparison safe: two sums of the same terms, added in different orders, can
differ in their last bits and still round onto adjacent `1e-9` steps. Refinement
therefore takes both compared values from one column-score matrix, each summed
along its own path in path order, as described under refinement above. Other
comparisons use total ordering after rejecting NaN. The deterministic priority
order is:

```text
DP: diagonal > consume left > consume right
tree merge: lexicographically smaller sorted leaf vectors
profile consensus: larger count > canonical alphabet order
cluster candidate: larger shared-kmer count > smaller stable index
refinement: strict improvement only
```

Checked size arithmetic rejects dimensions whose allocation or worst-case score
count would overflow `usize` or fixed objective units.

## Data structures

```text
HmmParameters       normalized transition and emission event sets
PreparedHmm         probability/log tables plus alphabet emission lookup
SparsePosterior     dimensions plus sorted sparse rows
PosteriorLibrary    triangular pair map with transposed access
ExpectedErrorMatrix symmetric dense pair distances
Muscle5Tree         deterministic binary nodes, members, heights, limbs
WorkingProfile      row indices, aligned cells, original residue positions
ProgressiveStep     children, score, operations, width, posterior pair count
RefinementAttempt   partition, before/after objective, decision
Super5Evidence      representatives, clusters, subdivisions, PPP-pa/ps merges
PairHmmWorkspace    reusable per-worker Forward, Backward, and dense buffers
Muscle5Alignment    MSA plus complete auditable evidence
Muscle5Ensemble     ordered, named alignment replicates
```

## Pseudocode

```text
PPP(records, seed, permutation, parameters):
    validate(records, parameters)
    H = published_hmm_for_alphabet()
    if seed != 0: H = perturb_and_normalize(H, seed, alpha=.25)

    P = ordered_parallel_pair_posteriors(records, H, cutoff=.01)
    repeat consistency_iterations times (default 2):
        P = ordered_parallel_consistency(P, identity_posteriors, cutoff=.01)

    for each sequence pair x,y:
        path_xy = maximum_expected_accuracy(P_xy)
        EE_xy = 1 - posterior_sum(path_xy) / path_width(path_xy)

    T = upgma(EE)
    T = permute(T, permutation)
    A = progressive_profile_alignment(records, T, P)

    repeat refinement_iterations times (default 100):
        (X,Y) = deterministic_random_nonempty_partition(A.rows)
        (candidate, S) = ppp_profile_align(strip_gap_columns(X), strip_gap_columns(Y), P)
        current = path_sum(A as a path over the column scores S)
        if fixed(candidate.score) > fixed(current): A = candidate

    verify_replay_and_rectangularity(A, records)
    return A and evidence

SUPER5(records, seed, permutation, parameters):
    ordered = decreasing_length(records)
    (representatives, deferred) = uclust_ee(ordered, .01, max_rejects=16)
    coarse = uclust_ee(representatives, .30, max_rejects=16)
    coarse = subdivide_large(coarse, .10, then_random, max_size=500)
    cluster_msas = map(cluster -> PPP(cluster), coarse)
    consensuses = map(consensus, cluster_msas)
    T = size_weighted_upgma(all_pair_EE(consensuses))
    Q = direct_untransformed_posteriors(required_original_pairs)
    A = progressive_cluster_merge(T,
          ppp_pa(Q) if cross_pairs <= 2000 else ppp_ps(Q, sample=2000))
    A = reinsert_redundant_rows(A, deferred)
    verify_replay_and_rectangularity(A, records)
    return A and evidence

STRATIFIED(records, K=4):
    return [pipeline(records, seed, permutation)
            for seed in 0..K
            for permutation in [none,abc,acb,bca]]

DIVERSIFIED(records, R=100):
    return [pipeline(records, r, permutations[r mod 4]) for r in 0..R]
```

## Complexity

Let `N` be sequence count, `L` a typical sequence length, `W` a typical profile
width, and `K` the average retained posterior entries per row.

- Pair-HMM forward/backward time is `O(N^2 L^2)` and dense temporary space is
  `O(L^2)` per pair; sparse resident posterior space is `O(N^2 L K)`.
- One sparse consistency iteration is at most `O(N^3 L K^2)` time and
  `O(N^2 L K)` resident space. Dense worst case is `O(N^3 L^3)` time and
  `O(N^2 L^2)` space, matching the ProbCons analysis.
- Expected-accuracy pair decoding is `O(N^2 L^2)` time and `O(L^2)` traceback
  space in the scalar correctness profile.
- UPGMA is `O(N^2 log N)` time with a stale-safe priority queue and `O(N^2)`
  space; a scan implementation is `O(N^3)` and is not used.
- Progressive alignment scatters retained sparse posterior entries into a dense
  column-score matrix, then runs the DP. Worst-case time remains `O(N^2 W^2)`
  across the tree. Working space is `O(W^2 + NW)`: `O(W^2)` column scores and
  byte traceback, `O(W)` rolling scores, and `O(NW)` profile cells/maps.
- One randomized refinement is `O(N^2 W^2)` in the worst balanced partition;
  100 rounds multiply this term by 100.
- Super5 bounds each PPP cluster at 500 rows after subdivision. U-sort candidate
  ranking is `O(N^2)` in the worst case but limits expensive EE rejections to 16
  per query; consensus alignment is quadratic in cluster count and consensus
  width. PPP-ps bounds a cluster-profile column score to 2000 sequence-pair
  contributions. The store of UCLUST posteriors holds one sparse matrix per
  comparison, `O(C L K)` for `C` comparisons, at most 16 per query.

Every result distinguishes resident posterior storage, peak pair-HMM scratch,
peak profile-DP scratch, returned MSA bytes, and an estimated algorithm-owned
peak. Process RSS is measured separately by the benchmark worker. In Super5, the
pair-HMM scratch estimate describes computing every requested pair, so with
reused posteriors the actual scratch can only be lower; the store of UCLUST
posteriors is not part of the estimate.

## Component test plan

Implementation is divided into independently testable components:

1. default/perturbed HMM normalization and symmetry;
2. canonical and ambiguous emissions;
3. guarded probability-domain Forward/Backward against the independent log-domain
   evaluator, plus independently enumerated tiny paths;
4. posterior transpose, sparsification, and checksum stability;
5. expected-accuracy DP and EE calculation;
6. one and two consistency transformations against direct dense multiplication;
7. deterministic UPGMA and four tree permutations;
8. PPP-pa profile-column scores and traceback;
9. strict randomized refinement decisions;
10. U-sort, clustering thresholds, subdivision, consensus, PPP-ps sampling, and
    redundant-row insertion;
11. stratified and diversified replicate schedules; and
12. public pipeline invariants, score verification, datasets, and benchmarks.

The test matrix also forces sequential and parallel collectors on identical jobs
and requires identical pair order, probability bits, checksums, and memory-accounting
semantics.

## Acceptance criteria

The algorithm is complete only when all of the following hold:

- this contract exists before implementation and every departure from an
  underspecified publication detail remains visible;
- `Align`, `Super5`, stratified ensembles, and diversified ensembles compile for
  protein, DNA, and RNA and preserve every ungapped input;
- hand-computed tests cover HMM path probability, posterior probability,
  expected accuracy, EE, consistency, guide trees, permutations, profile score,
  refinement, Super5 clustering, consensus, sampling, and reinsertion;
- empty, singleton, identical, disjoint, terminal-gap, internal-gap, repetitive,
  ambiguity-symbol, and duplicate-identifier cases are covered;
- an independent verifier recomputes posterior path scores, final profile
  objective, row replay, and every accepted refinement without production
  traceback or objective code;
- property tests cover rectangularity, replay, determinism, input-order
  stability, posterior bounds, HMM normalization, transpose symmetry, tree leaf
  coverage, strict refinement monotonicity, and ensemble schedules;
- immutable black-box fixtures generated by the pinned MUSCLE 5.1.0 executable
  cover `Align`, `Super5`, all guide permutations, perturbed seeds, proteins,
  DNA, RNA, ties, repeats, and divergent families;
- reference comparisons report exact aligned-string agreement plus residue-pair
  SP and TC, so biologically equivalent row orders are not counted as failures;
- BAliBASE subset SP/TC results and fixture provenance are recorded;
- release benchmarks report wall time, algorithm-owned memory, and process peak
  memory for Rust and explicitly configured MUSCLE 5.1.0 on identical inputs;
- `cargo test --workspace`, strict Clippy, formatting, documentation,
  dependency-license checks, and source-hygiene scans pass; and
- the public API is documented, has an example for all three molecule types,
  contains no placeholder, `TODO`, mock, or unimplemented path, and invokes no
  external aligner.

No expected alignment may be generated by the Rust implementation under test.

## Recorded validation evidence

This evidence was recorded for version 1 and is unchanged. The checks of
version 2 are in its validation record.

- All hand-derived, exhaustive, property, independent-replay, and static-reference
  tests pass for the three supported alphabets and both top-level modes.
- Against 15 pinned official MUSCLE 5.1 cases, 9 are byte-identical; mean SP is
  `0.889090956` and mean TC is `0.779340659`. Every per-case and aggregate gate
  passes.
- On the predeclared six-family BAliBASE 3 subset, macro SP/TC are
  `0.9132516705990265`/`0.6491193278793675`; micro SP/TC are
  `0.8934071464519376`/`0.8300395256916996`.
- Under the matched static Linux process protocol, Rust is faster in all five
  predeclared Align/Super5 workloads, with median speedups from 1.13x to 4.36x.

Exact fixtures, dataset hashes, measurement ranges, memory scopes, and
reproduction commands are recorded in
`../validation/black-box.md`,
`../benchmarks/muscle5.md`, and
`../benchmarks/muscle5-vs-reference.md`.

## Known scientific and compatibility limits

- The nucleotide emission table is absent from the paper and manual. The
  explicit symmetric model above is scientifically coherent but not claimed to
  reproduce an undisclosed reference table.
- The paper does not identify its PRNG. Fixed seeds reproduce this Rust profile,
  not necessarily byte-identical perturbed reference replicates.
- Several Super5 operational details are absent; the deterministic resolutions
  above preserve every published pipeline stage but can change cluster boundaries
  relative to the executable.
- Stable input order intentionally differs from default reference presentation.
- Pair and consistency jobs use deterministic scoped threads; their numerical
  kernels and profile DP remain scalar and contain no SIMD. Measured wins on the
  pinned process corpus are evidence for those workloads, not a universal speed
  claim across hardware, sequence distributions, or large datasets.
- Structural and genomic alignment remain outside scope.

These limitations must be quantified by black-box and scientific-dataset results;
none permits replacing PPP or Super5 with a simpler alignment algorithm.

## Version history

| Version | Date | Refinement compares |
|---|---|---|
| `MUSCLE-5.1.0-Nature-2022-default-v1` | 2026-08-04 release gate | The DP optimum with an independent row-pair sum of the current objective |
| `MUSCLE-5.1.0-Nature-2022-default-v2` | 2026-10-22 | The DP optimum with the current path summed over the same column scores |

Parameters, posteriors, guide trees, progressive alignment and Super5 are
unchanged. The two versions can decide differently only when the current
alignment's row-pair sum and its path sum round onto different `1e-9` steps.

On 2026-10-23, Super5 began reusing the pair posteriors of its UCLUST passes
within a run. That changes no output, so the profile identifier stays
`MUSCLE-5.1.0-Nature-2022-default-v2`; the
validation record holds the checks.

In the publication benchmark,
version 1 stopped with `PPP refinement candidate decreased its partition
objective` on 61 of 257 `muscle5` cases and 40 of 264 `muscle5-super5` cases.
For `muscle5`, none of 93 cases under 2,000 residues failed, but 45 of 63 above
10,000 did. In instrumented runs of the three smallest failed cases, every failing
attempt realigned a partition into the current MSA itself: its row-pair objective
was exactly the current one, while the DP's column-by-column sum of the same terms
was 1.4 to 5.5 × 10⁻¹² lower and rounded onto the `1e-9` step below. The
validation record holds those runs and
the checks of version 2; the rationale is
decision 030.
