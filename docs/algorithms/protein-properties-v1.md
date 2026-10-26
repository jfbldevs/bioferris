# Protein property model and numerical contract

## Contents

- [Scope and implementation state](#scope-and-implementation-state)
- [Mass and composition properties](#mass-and-composition-properties)
- [Charge and isoelectric point](#charge-and-isoelectric-point)
- [Reference observations](#reference-observations)
- [Interfaces and acceptance](#interfaces-and-acceptance)

[API index](../README.md) · [Complete workflows](../api/protein-properties-v1.md)
· Six-proposal plan

## Scope and implementation state

Stage 2 accepted public Rust/CLI/JSON interfaces and verified examples for the
documented models; integration evidence
records completed gates. The CLI/JSON interfaces were later removed with the
`bioferris` CLI (decision 029). All five requested properties are
in scope: average molecular weight, charge at a specified pH, isoelectric point,
aromaticity and Kyte-Doolittle GRAVY. The model is a nonempty, linear,
unmodified peptide containing only the canonical 20 amino acids. The existing
broader protein alphabet is retained, with explicit property-analysis errors for
ambiguous/noncanonical residues. Nonfinite pH and pH outside inclusive 0..14 are
explicit BioFerris errors. No sequence symbols are dropped or averaged.

## Mass and composition properties

Use average free-amino-acid masses in Da and subtract `(length - 1) * 18.0153`
for peptide-bond water loss. Tables and numerical reduction order must not be
mixed with differently rounded residue-mass tables. Sources:

- [Young et al., Science 2018 supplementary Table S1](https://ora.ox.ac.uk/objects/uuid%3Ab6403cfa-56ae-4675-981b-cc77f69e7dc6/files/m3613772303ddf27d4bea7f8a94b033c2).
- [NIST water molecular weight](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=CAL).
- [Biopython 1.88 ProtParam API](https://biopython.org/docs/1.88/api/Bio.SeqUtils.ProtParam.html).

```text
A 89.0932  C 121.1582 D 133.1027 E 147.1293 F 165.1891
G 75.0666  H 155.1546 I 131.1729 K 146.1876 L 131.1729
M 149.2113 N 132.1179 P 115.1305 Q 146.1445 R 174.2010
S 105.0926 T 119.1192 V 117.1463 W 204.2252 Y 181.1885
```

Aromaticity is the sum of F, W and Y frequencies. GRAVY averages the sequence's
Kyte-Doolittle values, from [ExPASy ProtScale](https://web.expasy.org/protscale/pscale/Hphob.Doolittle.html):

```text
A 1.8 C 2.5 D -3.5 E -3.5 F 2.8 G -0.4 H -3.2 I 4.5 K -3.9 L 3.8
M 1.9 N -3.5 P -1.6 Q -3.5 R -4.5 S -0.8 T -0.7 V 4.2 W -0.9 Y -1.3
```

Full-precision floating-point reduction order matters. An independently authored
Neumaier compensated scalar sum accumulates mass and hydropathy values in source
sequence order. Aromaticity adds the separately divided F, W and Y fractions,
in that order; combining the counts before division can change the last bits.
These reductions match the retained observations exactly on the capture platform.

## Charge and isoelectric point

Use the independent-ionization model: each basic group contributes
`count / (1 + 10^(pH - pKa))`; each acidic group subtracts
`count / (1 + 10^(pKa - pH))`. A linear peptide has one free N-terminus and one
free C-terminus. The published Bjellqvist baseline is N-terminal 7.50,
C-terminal 3.55, C 9.00, D 4.05, E 4.45, H 5.98, K 10.00, R 12.00 and Y 10.00; see
[Kozlowski 2016, Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC5075173/).

The terminal-adjusted model uses N-terminal A 7.59 / E 7.70 / M 7.00 / P 8.36 /
S 6.93 / T 6.82 / V 7.44 and C-terminal D 4.55 / E 4.75.
The latter corrections apply to the free C-terminal
group, not the acidic residue's side chain. Other termini retain the baseline
values. Ordinary binary64 addition combines positive contributions in
N-terminal/K/R/H order and negative contributions in C-terminal/D/E/C/Y order,
then subtracts the totals. Compensating this charge sum changes the reference
bits and is deliberately not used. Placement and arithmetic order are checked
against every single residue, every ordered residue pair and longer cases.
Each contribution evaluates the reciprocal ionization fraction first, then
multiplies by the binary64 count; dividing the count by the denominator instead
is mathematically equivalent but can round differently.

Biopython's [public pI API](https://biopython.org/docs/1.88/api/Bio.SeqUtils.IsoelectricPoint.html)
specifies bisection with initial pH 7.775 and bounds 4.05..12. The first probe is
the specified 7.775, not the arithmetic midpoint of those bounds. The verified
algorithm returns the current probe when interval width is at most 0.0001.
Otherwise, positive charge moves the lower bound to the current probe,
nonpositive charge moves the upper bound to that probe, and the next probe is
the midpoint `(lower + upper) / 2`.
The published bracket does not contain the root for every peptide, so the result
must expose whether it is bracketed or approaches the lower/upper bound. A
boundary-limited estimate must not be described as a zero-charge solution.

These are theoretical unmodified-peptide estimates, not measurements or
post-translational/structural corrections. See the
[ExPASy model explanation](https://web.expasy.org/compute_pi/pi_tool-doc.html).

## Reference observations

An isolated official Biopython 1.88 wheel was used for permitted public-API
observations under the existing CPython 3.12.3 interpreter. The wheel digest and
isolation procedure are recorded in the [CDS provenance](../../seq-bio/tests/fixtures/cds-provenance.md#biopython-behavioral-snapshot).
Implementation source is not inspected or copied. Observations cover all 20 single
residues, all 400 ordered residue pairs, published examples and longer composition
and boundary cases, retaining IEEE 754 binary64 bits at 13 explicit pH values.
The [retained fixture and provenance](../../seq-bio/tests/fixtures/protein-properties-provenance.md)
pin all 7,344 binary64 values and input sequences. Independent temporary debug
and optimized Rust harnesses match them exactly. The public Rust API, the CLI
(since removed) and the recipe checks also passed, as did the unrestricted
workspace regression of 1,447 tests across 206 harnesses at that increment. No Python process is required to run the shipped toolkit or its
Rust tests.

JSON serialization uses round-trippable decimals. Exact deserialization requires
the existing `serde_json` dependency's `float_roundtrip` feature: the default
parser can change the last bits of a valid decimal. Typed request pH and parsed
numeric responses must retain the original binary64 value. No package-version
change is needed for this requirement.

The retained comparisons state the reference version, runtime, platform,
models, inputs and numerical fields. Rounded examples, model agreement and exact
binary64 agreement are distinct evidence; none implies universal floating-point
identity on all platforms or sequences.

## Interfaces and acceptance

The Rust API is a borrowed validated `ProteinAnalysis` with methods
`molecular_weight`, `charge_at_ph`, `isoelectric_point`, `aromaticity` and `gravy`,
plus an aggregate `ProteinProperties` result and explicit pH options. Short
sequence extension methods provide the common path. Typed errors identify empty
input, unsupported residues with offsets and invalid pH. A pI estimate exposes
its bracket/boundary status.

The removed CLI/strict JSON interface exposed the aggregate properties, selected
models and pH. Acceptance required independent model checks, exact
retained-reference comparisons where claimed, boundaries, malformed inputs,
Rust/CLI equality, complete examples, documentation updates and the plan's
integration gates. Monoisotopic masses, modified proteins
and other hydropathy/pKa models are separate additions.
