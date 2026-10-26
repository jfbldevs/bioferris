# Clustal Omega Progressive HMM Multiple Alignment

## Reproduced profiles and status

Status: complete for the two versioned default profiles defined below. Optional
features outside those profiles remain explicitly out of scope under "Known profile
boundaries."

This phase reconstructs the default unaligned-sequence pipeline under two explicit,
molecule-specific profiles:

- `CLUSTALO-1.0.2-default-protein-v1` reproduces the protein-only release described
  by Sievers et al. in 2011; and
- `CLUSTALO-1.2.4-default-nucleotide-v1` is the companion DNA/RNA profile. Explicit
  nucleotide support first appeared in version 1.1.0; version 1.2.4 is selected as
  the stable public reference release for black-box validation.

The split is mandatory scientific provenance. Version 1.0.2 accepted nucleotide
letters but interpreted them as amino acids, so attributing a DNA or RNA method to
the 2011 profile would be false. DNA and RNA use one nucleotide model; `T` and `U`
occupy the same model state while typed output preserves the input molecule.

Both profiles preserve the characteristic published pipeline:

1. one-tuple distances for unaligned sequences;
2. the modified mBed embedding and bisecting clustering scheme;
3. real within-cluster distance matrices and UPGMA guide trees;
4. conversion of every sequence or intermediate alignment into a profile HMM;
5. progressive HMM-HMM alignment in guide-tree order;
6. posterior-probability Maximum Accuracy (MAC) alignment when it fits the memory
   budget; and
7. the five-state Viterbi HMM-HMM fallback when MAC exceeds that budget.

The scalar equations are the correctness contract. A validated optional worker
count parallelizes independent complete-matrix distance pairs and independent
progressive subtrees. Collection order, guide topology, merge trace, scores, and
alignment are identical to the one-worker execution.

This is an independent reconstruction from publications, supplementary material,
published equations, author papers, and official manuals. It does not call, link,
translate, inspect, or derive from Clustal Omega, HHalign, MUSCLE, or other aligner
source code.

The profile covers the ordinary command-line operation on a nonempty collection of
unaligned sequences. External Profile Alignment, HMMER file import, two pre-aligned
profile input, and optional guide-tree/HMM iteration are separately observable
Clustal Omega features, not hidden branches of this default profile. They are listed
as extensions rather than silently approximated here.

## Scientific sources

Primary Clustal Omega source:

1. Fabian Sievers, Andreas Wilm, David Dineen, Toby J. Gibson, Kevin Karplus,
   Weizhong Li, Rodrigo Lopez, Hamish McWilliam, Michael Remmert, Johannes Soding,
   Julie D. Thompson, and Desmond G. Higgins. "Fast, scalable generation of
   high-quality protein multiple sequence alignments using Clustal Omega."
   *Molecular Systems Biology* 7, 539, 2011.
   DOI: [10.1038/msb.2011.75](https://doi.org/10.1038/msb.2011.75).
2. The official supplementary information to source 1. It specifies the modified
   mBed construction, cluster limit of 100, UPGMA subtrees, profile-HMM alignment,
   MAC/Viterbi switch, terminal-core handling, EPA, iteration, benchmark versions,
   and datasets.

Primary component sources:

3. Gordon Blackshields, Fabian Sievers, Weifeng Shi, Andreas Wilm, and Desmond G.
   Higgins. "Sequence embedding for fast construction of guide trees for multiple
   sequence alignment." *Algorithms for Molecular Biology* 5, 21, 2010.
   DOI: [10.1186/1748-7188-5-21](https://doi.org/10.1186/1748-7188-5-21).
4. Johannes Soding. "Protein homology detection by HMM-HMM comparison."
   *Bioinformatics* 21(7), 951-960, 2005.
   DOI: [10.1093/bioinformatics/bti125](https://doi.org/10.1093/bioinformatics/bti125).
5. Ian Holmes and Richard Durbin. "Dynamic Programming Alignment Accuracy."
   *Journal of Computational Biology* 5(3), 493-504, 1998.
   DOI: [10.1089/cmb.1998.5.493](https://doi.org/10.1089/cmb.1998.5.493).
6. Andreas Biegert and Johannes Soding. "De novo identification of highly diverged
   protein repeats by probabilistic consistency." *Bioinformatics* 24(6), 807-814,
   2008. DOI: [10.1093/bioinformatics/btn039](https://doi.org/10.1093/bioinformatics/btn039).
7. W. John Wilbur and David J. Lipman. "Rapid similarity searches of nucleic acid
   and protein data banks." *Proceedings of the National Academy of Sciences*
   80(3), 726-730, 1983.
   DOI: [10.1073/pnas.80.3.726](https://doi.org/10.1073/pnas.80.3.726).
8. Julie D. Thompson, Desmond G. Higgins, and Toby J. Gibson. "CLUSTAL W:
   improving the sensitivity of progressive multiple sequence alignment through
   sequence weighting, position-specific gap penalties and weight matrix choice."
   *Nucleic Acids Research* 22(22), 4673-4680, 1994.
   DOI: [10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
9. Walter M. Fitch and Emanuel Margoliash. "Construction of phylogenetic trees."
   *Science* 155(3760), 279-284, 1967.
   DOI: [10.1126/science.155.3760.279](https://doi.org/10.1126/science.155.3760.279).
10. Robert R. Sokal and Charles D. Michener. "A statistical method for evaluating
   systematic relationships." *University of Kansas Science Bulletin* 38,
   1409-1438, 1958. This is the original UPGMA source.
11. Steven Henikoff and Jorja G. Henikoff. "Position-based sequence weights."
   *Journal of Molecular Biology* 243(4), 574-578, 1994.
   DOI: [10.1006/jmbi.1994.1644](https://doi.org/10.1006/jmbi.1994.1644).
12. Gaston H. Gonnet, Mark A. Cohen, and Steven A. Benner. "Exhaustive matching of
    the entire protein sequence database." *Science* 256(5062), 1443-1445, 1992.
    DOI: [10.1126/science.1604319](https://doi.org/10.1126/science.1604319).

Version, implementation-behavior, and manual sources:

13. Fabian Sievers and Desmond G. Higgins. "Clustal Omega for making accurate
    alignments of many protein sequences." *Protein Science* 27(1), 135-145, 2018.
    DOI: [10.1002/pro.3290](https://doi.org/10.1002/pro.3290). This paper records
    the 1.0.2 protein-only boundary, the 1.1.0 nucleotide addition, the nucleotide
    identity matrix, BRaliBase II calibration, MAC memory factors, and deterministic
    thread behavior.
14. Kieran Boyce, Fabian Sievers, and Desmond G. Higgins. "Instability in
    progressive multiple sequence alignment algorithms." *Algorithms for Molecular
    Biology* 10, 26, 2015.
    DOI: [10.1186/s13015-015-0057-1](https://doi.org/10.1186/s13015-015-0057-1).
    This source documents the prevalence and effect of tied k-tuple distances.
15. The archived official Clustal Omega 1.0.2 README. It specifies mBed as the
    default, real matrices for clusters of at most 100, UPGMA, input/profile modes,
    iteration behavior, and the initial-versus-aligned distance distinction.
16. The HHalign manual corresponding to the published substitution-matrix
    pseudocount mode. It records the diversity-dependent admixture, transition
    priors, free terminal gaps, and MAC threshold. Manual values resolve parameters
    omitted by papers; each such use is identified below.

Published immutable Gonnet scores are taken from the NCBI BLAST matrix distribution,
not from an aligner implementation. The reviewed matrix is
`https://ftp.ncbi.nih.gov/blast/matrices/GONNET`; its SHA-256 is recorded beside the
generated numeric table in the implementation.

Independent black-box validation uses the Clustal Omega 1.2.4 executable only in
the isolated sibling validation laboratory. The executable is not a Cargo
dependency and is never called by library tests. Static observed alignments,
distance matrices, guide trees, and timing records are allowed fixtures; reference
source is not inspected.

## Inputs and outputs

Input is an ordered, nonempty collection of uniquely named, typed, ungapped
sequences. All records in one invocation have the same alphabet. Empty individual
sequences are accepted. The result contains:

- a validated rectangular multiple alignment in original input order;
- the deterministic rooted UPGMA guide tree;
- the guide-tree method and, when a full matrix was materialized, its pair distances;
- mBed seed and cluster metadata when embedding was used;
- every progressive merge, selected HMM-HMM engine, core interval, score, and memory
  observation; and
- peak concurrent workspace metadata separated from returned-alignment storage.

Removing gaps from every output row must reproduce its corresponding input exactly.
Every row has the same width and no output column consists entirely of gaps.

## Stage 1: one-tuple distance

The supplement explicitly identifies the guide distance as the Wilbur--Lipman
k-tuple measure used by ClustalW 1.83/2. The original mBed experiments used protein
tuples as long as two; the reproduced command profiles and all 1.2.4 observations
use `k=1`. A one-tuple hit has score one, a mismatch has score zero, starting an
internal gap event costs `g=3`, extending that event costs zero, and unmatched
terminal material is free.

The Wilbur--Lipman objective is specified first as an exact, full-window,
three-state scalar recurrence. This representation is mathematically equivalent for
`k=1`, is easy to audit exhaustively, and remains the optimization oracle. For
symbols `x_i,y_j`, let `M`, `X`, and `Y` end in a paired column, a gap in `y`, and a
gap in `x`, respectively:

```text
M(i,j) = max(M(i-1,j-1), X(i-1,j-1), Y(i-1,j-1)) + [x_i = y_j]
X(i,j) = max(M(i-1,j)-g, Y(i-1,j)-g, X(i-1,j))
Y(i,j) = max(M(i,j-1)-g, X(i,j-1)-g, Y(i,j-1))

M(i,0) = M(0,j) = 0
X(i,0), Y(0,j) = unreachable
S(x,y) = max(0, M, X, Y over the last row or last column)
d(x,y) = 1 - S(x,y) / min(|x|, |y|).
```

For two empty sequences `d=0`; when exactly one sequence is empty `d=1`.
Distances are quantized to `DISTANCE_SCALE = 1,000,000,000` by round-to-nearest
with upward half ties. Typed symbols are compared by exact alphabet value; DNA `T`
and RNA `U` occur in different typed invocations and are never mixed implicitly.
IUB ambiguity and extended protein symbols are exact observable symbols rather than
undocumented wildcards.

Production uses the equivalent sparse hit-chain evaluator. A chain of increasing
one-tuple hits changes alignment diagonal exactly when it opens one internal gap
event. For a hit `(i,j)` on diagonal `delta=i-j`,

```text
H(i,j) = 1 + max(0,
                   best earlier hit on diagonal delta,
                   best earlier hit with i'<i and j'<j minus g).
```

Hits from one `i` row are committed only after the row is evaluated, enforcing both
strict inequalities. A Fenwick prefix-maximum tree supplies the dominance query;
an indexed diagonal table supplies the no-gap transition. If `r` is the number of
equal-symbol pairs, this uses `O(r log(min(|x|,|y|)))` time and linear auxiliary
memory. Worst-case single-symbol repeats remain quadratic-logarithmic; the dense
rolling oracle is `O(|x|*|y|)` time and `O(min(|x|,|y|))` memory.

The sparse and dense scores are bit-identical for every pair of binary sequences
through length four and four gap penalties. A separate alignment-path enumeration
checks the public distance and all boundary conventions. The five static Clustal
Omega 1.2.4 cases reproduce every printed distance to the reference's six-decimal
precision.

## Stage 2: modified mBed and deterministic bisecting k-means

For `N` records, the Clustal Omega supplement changes the original mBed seed count
from `(log2 N)^2` to `log2 N`. This profile uses

```text
t = max(1, ceil(log2 N)).
```

Records are stably sorted by `(length,input_index)`. Initial seeds are sampled at a
constant stride over that order, as published for mBed; duplicate biological
sequences are removed while retaining the earlier sampled index. If deduplication
leaves no seed, the smallest input index is retained. Each sequence becomes

```text
F(s) = [d(s,R1), d(s,R2), ..., d(s,Rt)].
```

Squared Euclidean distance drives k-means assignment. Euclidean distances used in
trees are obtained by deterministic integer square root and normalized by
`sqrt(t)`, so topology does not depend on a platform `libm` rounding choice.

The publications name bisecting k-means and k-means++ but do not publish random
seeds, empty-cluster repair, stopping rules, or tie-breaking. The reproducible
resolution is:

1. use a public fixed `u64` seed, with a sub-seed derived from the sorted member
   indices for every split;
2. select the first centroid uniformly with SplitMix64 and the second by the
   k-means++ squared-distance distribution;
3. assign to the nearer centroid, preferring centroid zero on an exact tie;
4. if a side is empty, move the point farthest from the occupied centroid, then the
   smallest input index;
5. recompute each integer barycentre with nearest rounding;
6. stop when assignments do not change, or after 100 iterations; and
7. if all vectors coincide, split the stable member order in half.

An oversized cluster is repeatedly bisected until every leaf cluster contains at
most 100 sequences. The next split is the eligible cluster with largest
within-cluster sum of squared errors, then smallest input index. The fixed seed is a
documented reconstruction parameter, not a claim about an unpublished historical
seed.

For every leaf cluster, real one-tuple distances for all member pairs are computed
and a UPGMA subtree is built. Cluster barycentres are joined by a second UPGMA tree
using vector distances, and each barycentre leaf is replaced by its real subtree.
For at most 100 records the default path is therefore an effective full distance
matrix and one UPGMA tree, as stated by the official README. Forcing full mode skips
embedding and constructs the same kind of real matrix directly.

## Stage 3: deterministic Clustal Omega guide clustering

The supplement calls this stage UPGMA and states that Clustal Omega uses the fast
MUSCLE routine, but it does not publish that routine's update and tie rules. Source
code inspection is prohibited. Static 1.2.4 guide trees distinguish conventional
leaf-count-weighted UPGMA from the observable equal-child-node update below; the
reconstruction records the latter as black-box-resolved behavior rather than
mislabeling it as the original Sokal--Michener equation.

For active nodes `A` and `B`, select the smallest distance. Exact ties use the
lexicographically smallest pair of minimum original input indices. If the selected
distance is `d(A,B)`, the new node height is

```text
h(U) = max(ceil(d(A,B)/2), h(A), h(B)).
```

The `max` is an explicit safeguard for quantized, non-ultrametric observations.
Child limbs are `h(U)-h(child)` and are never negative. For another active node `C`,
each child node has equal weight regardless of how many leaves it contains:

```text
d(U,C) = round_nearest((d(A,C) + d(B,C)) / 2).
```

All topology decisions and updates use checked fixed-point integers. Children
retain the deterministic active-slot order of the selected pair; clade comparisons
ignore left/right presentation. The public tree converts limb lengths to `f64` only
after topology is fixed. All five fixture topologies and every printed branch length
agree with Clustal Omega 1.2.4 within its Newick precision. A unit case that produces
a different third merge under leaf-weighted UPGMA protects this inferred behavior.

## Stage 4: profile-HMM construction

Every non-all-gap alignment column becomes one match state. A single sequence is
therefore a profile HMM with one match state per residue. Intermediate alignment
gaps define delete paths. Insert states are present at every boundary and receive
published transition priors even when the current alignment has no observed insert
column.

### Position-based sequence weights

For column `i`, let `r_i` be the number of distinct observed residue categories and
`n_i(a)` the number of rows containing category `a`. The unnormalized Henikoff
weight of sequence `s` is

```text
w_s = sum_(i where s_i is a residue) 1 / (r_i * n_i(s_i)).
```

Weights are normalized to sum to the row count. A row with no observed residue gets
the smallest positive weight assigned to an observed row; if all rows are empty,
weights are uniform. This empty-row rule is not specified by the papers and is made
explicit for total API behavior.

### Emissions and effective diversity

Weighted canonical counts in column `i` are normalized to observed frequencies
`f_i(a)`. Ambiguity symbols divide their weight uniformly among represented
canonical residues. Effective diversity is the entropy number

```text
Neff(i) = 2 ^ (-sum_a f_i(a) * log2(f_i(a))).
```

The substitution-matrix pseudocount prior is

```text
g_i(a) = sum_b Q(a|b) * f_i(b)
tau(i) = pca / (1 + ((Neff(i)-1)/pcb)^pcc)
e_i(a) = (1-tau(i))*f_i(a) + tau(i)*g_i(a).
```

The protein profile uses the published/manual `pca=1.0`, `pcb=1.5`, and `pcc=1.0`.
For protein,

```text
Q(a|b) = background(a) * 2^(GONNET(a,b)/2) / Z_b,
```

where the divisor two interprets the published integer table in half-bit units and
`Z_b` normalizes the column. This conversion is an explicit reconstruction of the
paper's "substitution matrix method" because it does not publish the table-to-prior
normalization.

For DNA/RNA, Sievers and Higgins specify an identity substitution matrix and state
that nucleotide pseudocount and gap parameters were calibrated on BRaliBase II, but
do not publish the calibrated values. The same normalized odds construction is
applied to a four-by-four matrix with diagonal one and off-diagonal zero:

```text
Q(a|b) = (1/4) * 2^([a=b]/2) / Z_b.
```

Canonical nucleotide background is uniform. `pcb=1.5` and `pcc=1.0` are retained;
`pca=0.55` is an explicit black-box calibration against 1.2.4, not a claimed value
from the paper. The public builder exposes it, the profile identifier versions it,
and tests isolate matrix semantics from that calibrated admixture. All emission
rows are renormalized and checked to be finite and positive.

### Transition probabilities

Weighted paths supply observed `M->M`, `M->D`, `D->M`, and `D->D` counts. Gaps before
a row's first residue or after its last residue are terminal overhangs and do not
create delete transitions. There are no observed insert transitions in an ordinary
rectangular profile, but insert states remain reachable through priors. The manual
values are

```text
prior(M->I) = prior(M->D) = gapd * 0.0286
prior(M->M) = 1 - prior(M->I) - prior(M->D)
prior(I->I) = prior(D->D) = gape / (gape - 1 + 1/0.75)
prior(I->M) = prior(D->M) = 1 - prior(I->I)
```

with `gapb=1.0`, `gapd=0.15`, and `gape=1.0`. An outgoing family is the normalized
sum of its weighted observations and `gapb` times its prior. Transitions are stored
at the HHalign boundary index used by the published recurrences, with a prior-only
trailing boundary. Unsupported `I->D` and `D->I` transitions are excluded as in
Soding (2005).

## Stage 5: HMM-HMM scoring and Viterbi alignment

For match columns `q_i` and `p_j`, the log-sum-of-odds column score from Soding is

```text
Saa(i,j) = log2(sum_a q_i(a)*p_j(a)/background(a)) + shift.
```

The default reconstruction shift is `-0.03` bits. The 2005 paper reports `-0.1`,
while later HHalign manuals report `-0.03`; Clustal Omega publications omit the
embedded value. The value is therefore a versioned public parameter and a black-box
calibration target, not an unqualified historical claim.

There are five pair states:

```text
MM, MI, IM, DG, GD.
```

`II` and `DD` are excluded. A state may extend itself; opening or closing a side
state passes through `MM`. Let `Tq_i(X,Y)` and `Tp_j(X,Y)` be profile transition
probabilities. Representative Viterbi recurrences are

```text
VMM(i,j) = Saa(i,j) + max(
    VMM(i-1,j-1) + log2(Tq(M,M)*Tp(M,M)),
    VMI(i-1,j-1) + log2(Tq(M,M)*Tp(I,M)),
    VIM(i-1,j-1) + log2(Tq(I,M)*Tp(M,M)),
    VDG(i-1,j-1) + log2(Tq(D,M)*Tp(M,M)),
    VGD(i-1,j-1) + log2(Tq(M,M)*Tp(D,M)))

VMI(i,j) = max(
    VMM(i-1,j) + log2(Tq_(i-1)(M,M)*Tp_j(M,I)),
    VMI(i-1,j) + log2(Tq_(i-1)(M,M)*Tp_j(I,I)))

VDG(i,j) = max(
    VMM(i-1,j) + log2(Tq_(i-1)(M,D)),
    VDG(i-1,j) + log2(Tq_(i-1)(D,D))).
```

`VIM` and `VGD` are symmetric, using `Tq_i` with `Tp_(j-1)` and `Tp_(j-1)`,
respectively. Closing any side state into `MM(i,j)` uses the two profiles' boundary
indices `i-1,j-1`, exactly as in Soding's equations 5--7. The four gap-transition log
terms are multiplied by the published/manual factors `0.6` for open and `0.6` for
extension. `MM(i,0)` and `MM(0,j)` are zero and other boundary states are
unreachable. The endpoint is the best `MM` cell in the last row or last column.
This is Soding's semi-global core; terminal gaps cost zero.

For entry to `MM`, Viterbi ties use the fixed predecessor order `MM, MI, IM, DG,
GD`. A side-state opening/extension tie prefers extension. Equal endpoints prefer
the larger consumed core and then the smaller left coordinate. These deterministic
rules resolve omissions in the publications. An independent enumerator visits
every legal five-state path for tiny profile pairs and checks both the optimum score
and the exact returned core traceback.

## Stage 6: Forward/Backward posterior probabilities and MAC

When the configured memory budget admits MAC, max is replaced by sum in the five
state recurrences. Let `F_X(i,j)` and `B_X(i,j)` be their mathematical log-space
Forward and Backward values and `Z` the log-sum over legal semi-global `MM`
endpoints. The match posterior is

```text
P(i,j) = exp2(F_MM(i,j) + B_MM(i,j) - Z).
```

The production evaluator computes the equivalent probability-space recurrences in
contiguous row order. Every Forward row has an explicit base-two scale `s_i` and
stores `F'_X(i,j)=F_X(i,j)*2^(-s_i)`; dependencies from the preceding row are
multiplied by the exact scale ratio. Backward uses the same construction in reverse
with two rolling rows. Rescaling occurs only outside high/low water bounds and uses
an integral power of two, avoiding both per-cell transcendental functions and silent
underflow. Endpoint partition terms and posterior terms reintroduce their recorded
scales. This is an arithmetic evaluation strategy, not an approximate HMM or a
different recurrence.

An independent five-state log-sum-exp implementation is compiled for tests. Across
900 ordered pairs of binary DNA profiles through length four under both side-tie
orientations (1,800 evaluations), plus eight indel/repeat/general evaluations, it
produces the same core
interval and operation sequence; scores and posterior bounds agree within
`1e-10`. Values are clamped only for final roundoff to `[0,1]`; a non-finite
internal value is an error. The MAC recurrence extended to local HMM-HMM alignment
by Biegert and Soding is

```text
A(i,j) = max(
    0,
    A(i-1,j-1) + P(i,j) - T,
    A(i-1,j)   - T/2,
    A(i,j-1)   - T/2).
```

The documented default greediness threshold is `T=0.35`. The best cell in the whole
matrix ends the local core. Traceback stops at zero. MAC ties prefer a posterior
match, then extension on the side with smaller minimum input index, then the other
side, then stop.

Forward values for all five states, one precomputed emission-weight matrix, the
match posterior matrix, and one-byte MAC trace are retained. Backward and MAC scores
use two rolling rows. The resident interior-cell model is therefore 57 bytes
(40 Forward, 8 emission, 8 posterior, and 1 trace), plus row scales, transition
weights, and rolling buffers. Automatic mode chooses MAC only when this checked
workspace is no larger than the configured limit. Otherwise it selects the
five-state Viterbi recurrence; the algorithm is never silently replaced by ordinary
affine profile alignment.

## Terminal core reattachment and progressive merging

Both engines align an overlapping HMM core. Every original profile column outside
that core must be restored. Left prefixes are right-aligned in
`max(left_prefix,right_prefix)` free columns; right suffixes are left-aligned in
`max(left_suffix,right_suffix)` free columns. This is the explicit reattachment rule
for an omission in the supplement. If MAC finds no positive core, complete profiles
are right-aligned in free terminal columns. That empty-core convention is a
documented deterministic resolution supported by the non-repetitive 1.2.4 fixtures.

`MM` consumes one column from each profile. `MI` and `DG` consume one column from
the first profile and insert a gap into every row of the second; `IM` and `GD` do the
opposite. The relative order of every column within each child profile remains
fixed. After every merge, an independent replay verifies sequence consumption,
operation legality, and absence of all-gap columns.

The guide tree is traversed in deterministic postorder until the root profile is
formed. With more than one worker, two nontrivial child subtrees may be evaluated
concurrently with a fixed worker-budget split. Their traces are concatenated in
left-postorder, right-postorder, parent order regardless of completion timing.
Complete guide-distance pairs use the same indexed deterministic collection rule.
Rows are restored to original input order only after the final merge.

## Own pseudocode

```text
function CLUSTAL_OMEGA(records, parameters):
    validate records and molecule/profile agreement

    if forced_full:
        D <- deterministic parallel map of all ONE_TUPLE_DISTANCE pairs
        tree <- UPGMA(D)
        guide_evidence <- FULL(D)
    else if count(records) <= 100:
        D <- deterministic parallel map of all ONE_TUPLE_DISTANCE pairs
        tree <- UPGMA(D)
        guide_evidence <- EFFECTIVE_FULL_MBED(D)
    else:
        seeds <- SYSTEMATIC_SEEDS(records, ceil(log2(N)))
        vectors <- distances of every record to seeds
        clusters <- DETERMINISTIC_BISECTING_KMEANS(vectors, maximum_size=100)
        for cluster in clusters:
            subtree[cluster] <- UPGMA(real one-tuple distances within cluster)
        supertree <- UPGMA(Euclidean distances between cluster barycentres)
        tree <- GRAFT(subtree, supertree)
        guide_evidence <- MBED(seeds, clusters)

    function ALIGN_NODE(node, worker_budget):
        if node is a leaf:
            return one-row profile for its input sequence

        if both children are branches and worker_budget > 1:
            left, right <- evaluate child subtrees concurrently with split budget
        else:
            left <- ALIGN_NODE(node.left, worker_budget)
            right <- ALIGN_NODE(node.right, worker_budget)
        q <- BUILD_PROFILE_HMM(left)
        p <- BUILD_PROFILE_HMM(right)

        if MAC_MEMORY(q,p) <= memory_limit:
            forward <- FIVE_STATE_FORWARD(q,p)
            posterior <- FIVE_STATE_BACKWARD_POSTERIOR(q,p,forward)
            core <- MAC_TRACEBACK(posterior, threshold=0.35)
            engine <- MAC
        else:
            core <- FIVE_STATE_VITERBI_TRACEBACK(q,p)
            engine <- VITERBI

        path <- REATTACH_FREE_TERMINAL_COLUMNS(core, left, right)
        merged <- MERGE_FIXED_CHILD_COLUMNS(left, right, path)
        REPLAY_AND_VALIDATE(merged, path)
        concatenate left trace, right trace, then parent step
        record concurrent memory high-water observation
        return merged

    root_profile <- ALIGN_NODE(tree.root, parameters.thread_count)
    restore original row order
    validate ungapped rows, rectangularity, and columns
    return alignment, tree, guide evidence, trace, memory metadata
```

## Parameters

The stable defaults and their provenance are:

| Parameter | Protein 1.0.2 | Nucleotide 1.2.4 | Provenance |
|---|---:|---:|---|
| Tuple size | 1 | 1 | authors' later analysis; nucleotide companion choice |
| One-tuple gap event | 3 | 3 | Wilbur--Lipman/Clustal family plus 1.2.4 black box |
| Guide mode | modified mBed | modified mBed | paper and official README |
| Maximum leaf cluster | 100 | 100 | supplement |
| Seed count | `ceil(log2 N)` | `ceil(log2 N)` | supplement |
| HMM residue matrix | Gonnet | identity | 2005 and 2018 papers |
| Emission `pca,pcb,pcc` | `1,1.5,1` | `0.55,1.5,1` | HHalign manual; explicit nucleotide black-box calibration |
| Transition `gapb,gapd,gape` | `1,0.15,1` | `1,0.15,1` | HHalign manual; RNA values remain calibration targets |
| Column shift | `-0.03` bits | `-0.03` bits | later manual; explicit ambiguity resolution |
| Gap path factors | `0.6` each | `0.6` each | HHalign manual |
| MAC threshold | `0.35` | `0.35` | HHalign manual |
| MAC memory limit | 2 GiB | 8 GiB | 1.0.2 and 1.2.3+ descriptions |
| Worker count | 1 | 1 | deterministic scalar default; explicit performance parameter |

The public builder permits force-full guide trees, force-MAC or force-Viterbi,
memory limit, MAC threshold, pseudocount strengths, and deterministic k-means seed
to be changed only after finite/range validation. Worker count is a checked positive
integer. The aligner exposes effective parameters, and each result retains the
stable scientific profile identifier plus all observable guide/merge decisions.

## Complexity

Let `N` be sequence count, `L` maximum raw length, `t=ceil(log2 N)`, and let
progressive merge `r` align HMM lengths `x_r,y_r`.

- A pair with `r` one-tuple hits uses `O(r log L)` time and `O(L+sigma)` auxiliary
  memory; `r=O(L^2)` in the repetitive worst case. Full distance construction is
  `O(sum_pairs r_pair log L)` time plus `O(N^2)` matrix space. The retained dense
  test oracle is `O(L^2)` time and `O(L)` memory per pair.
- mBed embedding evaluates only `N*t` such sparse pairs and retains `O(N*t)` vector
  space. Bisecting k-means is `O(I*N*t*log N)` for at most `I=100` Lloyd iterations
  per split. Real matrices inside bounded clusters evaluate at most `O(N*100)`
  additional pairs and use `O(100^2)` transient matrix space.
- This clarity-first UPGMA is `O(N^3)` in forced-full mode and
  `O(C^3 + N*100^2)` in mBed mode for `C` leaf clusters; its matrix space is
  quadratic in the active cluster count. A later optimized nearest-neighbor UPGMA
  must preserve the scalar topology oracle.
- Profile construction at merge `r` is
  `O(rows_r*(x_r+y_r) + (x_r+y_r)*sigma^2)` because substitution priors are dense.
- Viterbi is `O(sum_r x_r*y_r*sigma)` time and uses a full five-state score and
  traceback matrix plus a cached column-score matrix, hence `O(x_r*y_r)` space.
- MAC has the same asymptotic time and `O(x_r*y_r)` space with a 57-byte resident
  interior-cell model plus `O(x_r+y_r)` scaled/rolling buffers.
- With `P` workers, complete guide-distance work is partitioned deterministically.
  Progressive work has the same total complexity; its span follows the guide-tree
  dependency path and cannot parallelize a parent merge before both children finish.
  Concurrent peak workspace can rise by as much as the sum of active subtree peaks.
- The returned alignment occupies `O(N*W)` for final width `W`.

## Test and validation boundaries

The implementation is divided into independently testable components:

1. one-tuple recurrence, free boundaries, gap-event score, and exact fixed-point
   distances;
2. fixed-point UPGMA topology, heights, branch lengths, and tie order;
3. systematic seeds, embedding vectors, deterministic k-means++, cluster cap, and
   subtree grafting;
4. Henikoff weights and ambiguity projections;
5. emission pseudocounts, transition counts, priors, and normalization;
6. five Viterbi state recurrences, independent exhaustive optimum, and traceback;
7. Forward/Backward partition consistency, row scaling, and posterior bounds;
8. exhaustive scaled-versus-log-space MAC differential checks and threshold behavior;
9. terminal reattachment and profile merge replay;
10. complete protein, DNA, and RNA pipelines;
11. static Clustal Omega 1.2.4 observations from the isolated laboratory;
12. bit-identical one-worker/eight-worker pipeline execution; and
13. scientific SP/TC measurements plus time and memory benchmarks.

Expected values for hand cases are calculated independently on paper or by a small
oracle that implements the mathematical definition through a structurally different
enumeration. The production implementation never generates its own expected test
fixtures.

## Acceptance criteria

This phase is complete only when all of the following are evidenced:

- this specification contains the exact profiles, references, equations, own
  pseudocode, structures, parameters, tie rules, omissions, and explicit resolutions;
- protein, DNA, and RNA APIs are safe, deterministic, and preserve typed input;
- empty rows, identical and wholly different records, repeats, ambiguity symbols,
  and leading/internal/trailing gaps are covered;
- unit tests independently verify every component listed above;
- property tests verify determinism, row consumption, input-order preservation,
  rectangularity, no all-gap columns, normalized HMM probabilities, finite scores,
  and bounded posteriors;
- score and traceback are replayed independently;
- black-box cases compare alignments, one-tuple distance matrices, guide trees, SP,
  and TC where the reference exposes them;
- BAliBASE and HomFam subsets are evaluated, with a BRaliBase II nucleotide subset
  when licensing and stable retrieval permit redistribution or pinned local use;
- Criterion reports increasing lengths and sequence counts, both HMM engines, mBed,
  complete process time, and explicit peak memory;
- the isolated process benchmark compares Rust with Clustal Omega 1.2.4 under
  equivalent molecule type, guide mode, thread count, and output work;
- formatting, Clippy, all tests, rustdoc warnings, dependency licenses, advisories,
  bans, sources, fixture hashes, and local links pass;
- no external runtime, FFI, unsafe code, copied implementation, placeholder, mock,
  or unimplemented execution path exists; and
- all observed reference differences and performance limitations are documented.

## Known profile boundaries

- Default unaligned-sequence alignment is reproduced. EPA, HMMER parsing,
  pre-aligned profile input, user-provided guide trees, and iteration require their
  own documented profiles because their inputs and scientific contracts differ.
- The original 1.0.2 executable is protein-only. Nucleotide comparison against
  1.2.4 validates the later companion profile, not a nonexistent 2011 feature.
- Published descriptions omit several embedded HHalign values. Those values remain
  explicit, versioned reconstruction choices even if black-box calibration improves
  agreement.
- One-tuple distance uses residue order but remains a low-resolution identity-only
  heuristic and produces many ties. The repetitive DNA fixture has several
  score-equivalent placements; the deterministic Rust placement differs from the
  1.2.4 presentation while preserving a valid alignment.
- The dense one-tuple, log-space MAC, and exhaustive Viterbi evaluators remain
  test-only scalar correctness anchors. Production uses sparse hit chaining,
  row-scaled MAC, deterministic complete-distance/subtree parallelism, contiguous
  emissions, and cached transition/column weights. SIMD emission scores,
  cache-blocked DP, linear-space traceback, a repetitive-input fallback, parallel
  mBed embedding, and faster guide clustering remain optimization targets.
- Progressive alignment remains greedy: columns fixed in an earlier child profile
  cannot be rearranged at a later merge.
- The scientific dataset record covers six BAliBASE families and one HomFam family,
  not either complete corpus. A stably retrievable, license-reviewed BRaliBase II
  nucleotide subset was not acquired in this phase; nucleotide calibration is
  therefore supported by equations, component tests, and 1.2.4 black-box cases.
- One-worker process timing remains slower than the packaged 1.2.4 executable, and
  the smallest 16-worker case remains startup-bound. Rust crosses over on the two
  larger equal-thread process cases; exact measurements and OS caveats are recorded
  separately.
