# Protein properties: Rust workflows

## Contents

- [Scope and status](#scope-and-status)
- [Write a Rust analysis](#write-a-rust-analysis)
- [Errors and recovery](#errors-and-recovery)
- [Models, precision and limitations](#models-precision-and-limitations)

[API index](../README.md) · [Numerical contract](../algorithms/protein-properties-v1.md)
· [Reference provenance](../../seq-bio/tests/fixtures/protein-properties-provenance.md)

## Scope and status

The Rust implementation and its executable example are implemented and verified
for the documented models; see the
validation record. The CLI and JSON
adapters that record describes were removed with the `bioferris` CLI
(decision 029). The retained
reference corpus matches all 7,344 numerical values exactly on the observed
platform. That is a scoped compatibility result, not universal Biopython parity.

Analyze a nonempty, canonical-20, unmodified linear peptide for:

| Property | Result and unit | Model |
|---|---|---|
| Molecular weight | Average molecular mass in Da | Free amino acid masses minus peptide-bond water loss |
| Charge at pH | Net elementary charges; fractional estimates are normal | Terminal-adjusted Bjellqvist independent ionization |
| Isoelectric point | Bounded pH estimate with explicit boundary status | Specified bisection in 4.05..12 |
| Aromaticity | F/W/Y fraction, mathematically 0..1 | Separate phenylalanine, tryptophan and tyrosine frequencies |
| GRAVY | Mean hydropathy value | Kyte-Doolittle |

Only charge depends on the selected pH; its default is 7.0. A pI boundary marker
means the model root is outside the search interval, not that the returned value
has zero charge. No modifications, ambiguity averaging or alternative models are
inferred. Standard protein construction still accepts the broader historical
alphabet; this analysis explicitly rejects unsupported residues.

## Write a Rust analysis

Work from the repository root. Save this complete program as
[`seq-bio/examples/protein_analysis.rs`](../../seq-bio/examples/protein_analysis.rs).
Literal inputs are `INGAR`, `PETER`, `DD` and the intentionally unsupported `AX`.

<!-- protein-properties-rust-program -->
```rust
//! Canonical peptide properties with explicit pH, pI boundaries and input errors.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let peptide = protein("INGAR")?;
    let properties = peptide.protein_properties()?;

    println!("Protein: {peptide}");
    println!("Mass (Da): {:.4}", properties.molecular_weight());
    println!("Charge at pH 7: {:.6}", properties.charge());
    println!("Isoelectric point: {:.6}", properties.isoelectric_point());
    println!("Aromaticity: {:.6}", properties.aromaticity());
    println!("GRAVY: {:.6}", properties.gravy());

    let options = ProteinPropertiesOptions { ph: 4.53 };
    let custom = protein("PETER")?.protein_properties_with(options)?;
    println!("PETER charge at pH 4.53: {:.6}", custom.charge());

    let acidic = protein("DD")?.protein_properties()?;
    println!("DD pI boundary: {:?}", acidic.pi_estimate().boundary());

    if let Err(error) = protein("AX")?.protein_properties() {
        println!("Rejected AX: {error}");
    }

    let out_of_range = ProteinPropertiesOptions { ph: 15.0 };
    if let Err(error) = protein("INGAR")?.protein_properties_with(out_of_range) {
        println!("Rejected pH 15 ({}): {error}", error.code());
    }

    Ok(())
}
```

<!-- protein-properties-rust-command -->
```bash
cargo bio protein_analysis
```

Complete stdout, exit 0, empty stderr:

<!-- protein-properties-rust-output -->
```text
Protein: INGAR
Mass (Da): 529.5904
Charge at pH 7: 0.760092
Isoelectric point: 9.750021
Aromaticity: 0.000000
GRAVY: -0.420000
PETER charge at pH 4.53: 0.002796
DD pI boundary: Lower
Rejected AX: protein property analysis requires canonical amino acids; found X at residue 1
Rejected pH 15 (invalid_ph): protein charge pH must be finite and in 0..=14 (received 15)
```

The program deliberately formats human-readable numbers to four/six decimal
places; the returned `f64` values are not rounded. `DD` demonstrates the lower
pI boundary, while `AX` and pH 15 demonstrate refusals; see
[Errors and recovery](#errors-and-recovery).
The original sequences are immutable. For repeated individual queries,
`analyze_protein()` returns a borrowed `ProteinAnalysis`; its methods expose the
five calculations and the complete pI estimate without cloning the sequence.

## Errors and recovery

Analysis returns `Err(ProteinPropertyError)` instead of repairing input. Each
variant has a stable identifier from `code()`; the program above prints two of
them in its last lines.

| Variant | `code()` | Cause | Recovery |
|---|---|---|---|
| `EmptySequence` | `empty_protein` | No residues | Supply a peptide |
| `UnsupportedResidue { offset, byte }` | `unsupported_protein_residue` | B, Z, J, X, U or O at the zero-based `offset` | Choose the sequence or model deliberately; never delete X or substitute an average residue |
| `InvalidPh { ph }` | `invalid_ph` | A pH that is not finite or lies outside 0..=14, including NaN and infinities | Select a supported pH; values are never clamped |

`AX` is a valid general protein but not a canonical-20 property input: the error
names X at residue 1. Bytes outside the protein alphabet fail earlier, when
`protein()` builds the sequence.

## Models, precision and limitations

The [numerical contract](../algorithms/protein-properties-v1.md) pins all mass,
hydropathy and pKa constants, evaluation order, pI initialization and stopping
rule. Mass/GRAVY use compensated source-order sums; charge uses ordered ordinary
sums and reciprocal-first contributions; aromaticity adds separately divided
frequencies. Results are binary64 values; the program above rounds only for
display.

The observed reference is Biopython 1.88 / CPython 3.12.3 on Linux x86-64. Its
432 sequences include all single residues, all ordered pairs, a published
152-residue protein and longer stress/boundary inputs, with 13 pH values each.
The fixture and provenance are retained for offline Rust bit comparisons. Other
platforms/interpreters, arbitrary sequences or additional numerical models are
not covered by a universal exact-fidelity claim.

This is not monoisotopic mass, circular-peptide mass, modification analysis,
protein structure prediction, disulfide modeling, an extinction coefficient or
another hydropathy/pKa scale. pI boundary detection adds explicit BioFerris
metadata; it does not extend the reference's numerical bracket. There is no
performance claim. The validation record holds the checks that passed for this
scope.
