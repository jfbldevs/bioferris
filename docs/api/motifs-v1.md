# Nucleotide motifs: matching and scoring matrices

## Contents

- [Availability and purpose](#availability-and-purpose)
- [Write a Rust analysis](#write-a-rust-analysis)
- [Choose matching semantics](#choose-matching-semantics)
- [Build and score matrices](#build-and-score-matrices)
- [Defaults and resource limits](#defaults-and-resource-limits)
- [Errors and recovery](#errors-and-recovery)
- [Reference compatibility and verification](#reference-compatibility-and-verification)

[API index](../README.md) · [Numerical contract](../algorithms/motifs-v1.md) ·
Implementation plan

## Availability and purpose

The Rust API and Cargo example are implemented and verified on Rust 1.85.0.
Stage 3 is accepted for the documented matching/matrix profiles and retained
reference corpus; see the main validation record.
That record also covers CLI operations, strict JSON requests and capability
discovery, which were removed with the `bioferris` CLI
(decision 029). The complete program and the real `cargo bio` alias
below are checked by the tests. Scoped acceptance does not establish a package
release or universal native/API equivalence.

This module provides two different analyses for DNA and RNA:

| Question | Rust entry point |
|---|---|
| Does an IUPAC pattern match here? | `dna_motif` / `rna_motif`, `IupacMotif::find_with` |
| What count, probability and log-odds model do these sites define? | `MotifCounts` → `MotifWeights` → `MotifScores` |
| Which concrete windows reach my score threshold? | `MotifScores::find_with` |

These are not regexes, protein motifs, motif discovery, binding-site prediction,
or a significance test. A score cutoff is supplied by the caller; no p-value,
false-discovery rate or learned background is inferred.

## Write a Rust analysis

Manual use means writing an ordinary Rust program, as with other BioFerris
analyses. All input sequences, instances, custom columns and policies appear in
this complete example. DNA is normalized to uppercase; RNA uses U, not T.

Complete runnable example:
[`seq-bio/examples/motif_analysis.rs`](../../seq-bio/examples/motif_analysis.rs):

<!-- motifs-rust-program -->
```rust
//! IUPAC matching and explicit nucleotide count, probability and scoring workflows.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let source = dna("AGRYN")?;
    let motif = dna_motif("R")?;
    for ambiguity in [MotifInputAmbiguity::Definite, MotifInputAmbiguity::Possible] {
        let hits = motif.find_with(
            &source,
            MotifSearchOptions {
                ambiguity,
                ..Default::default()
            },
        )?;
        println!("{ambiguity:?}:");
        for hit in hits {
            println!("  {:?} {:?}", hit.strand(), hit.range());
        }
    }
    if let Err(error) = motif.find_with(
        &source,
        MotifSearchOptions {
            ambiguity: MotifInputAmbiguity::Reject,
            ..Default::default()
        },
    ) {
        println!("Reject: {}: {error}", error.code());
    }
    let no_windows = MotifSearchOptions {
        max_windows: 0,
        ..Default::default()
    };
    if let Err(error) = dna_motif("A")?.find_with(&dna("A")?, no_windows) {
        println!("Budget: {}: {error}", error.code());
    }

    let instances = [dna("AC")?, dna("AC")?, dna("AC")?, dna("CA")?];
    let counts = MotifCounts::from_instances(&instances)?;
    let weights = counts.normalize([1.0; 4])?;
    let scores = weights.log_odds_uniform()?;
    println!("Counts A/C/G/T: {:?}", counts.columns());
    println!("PWM A/C/G/T: {:?}", weights.columns());
    println!("PSSM A/C/G/T: {:?}", scores.columns());
    let concrete = dna("ACGTAC")?;
    println!("Forward scores: {:?}", scores.calculate(&concrete, 100)?);
    for hit in scores.find_with(
        &concrete,
        MatrixSearchOptions {
            threshold: 2.0,
            ..Default::default()
        },
    )? {
        println!(
            "Scored hit: {:?} {:?} {}",
            hit.strand(),
            hit.range(),
            hit.score()
        );
    }

    let custom = MotifCounts::<Dna>::from_columns(vec![[1.0; 4]])?
        .normalize_default()?
        .log_odds([4.0, 2.0, 1.0, 1.0])?;
    println!("Normalized background: {:?}", custom.background());
    println!("Asymmetric PSSM: {:?}", custom.columns());
    for hit in custom.find_with(
        &dna("A")?,
        MatrixSearchOptions {
            threshold: 1.0,
            ..Default::default()
        },
    )? {
        println!(
            "Asymmetric hit: {:?} {:?} {}",
            hit.strand(),
            hit.range(),
            hit.score()
        );
    }

    let sparse = MotifCounts::from_instances(&[dna("A")?])?
        .normalize_default()?
        .log_odds_uniform()?;
    println!(
        "Zero-probability score for C: {}",
        sparse.score_window(&dna("C")?)?
    );
    let cancellation = MotifScores::<Dna>::from_columns(
        vec![
            [1.0, 0.0, 0.0, 0.0],
            [2_f64.powi(-25), 0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 0.0],
        ],
        [0.25; 4],
    )?;
    println!(
        "Float64 accumulation, final float32: {}",
        cancellation.score_window(&dna("AAA")?)?
    );
    Ok(())
}
```

From the repository root:

<!-- motifs-rust-command -->
```bash
cargo bio motif_analysis
```

Complete observed program stdout; stderr is empty; exit 0:

<!-- motifs-rust-output -->
```text
Definite:
  Forward 0..1
  Forward 1..2
  Forward 2..3
  Reverse 3..4
Possible:
  Forward 0..1
  Forward 1..2
  Forward 2..3
  Reverse 3..4
  Forward 4..5
  Reverse 4..5
Reject: ambiguous_motif_input: motif workflow requires concrete bases; found R at nucleotide 2
Budget: limit_exceeded: motif windows count 2 exceeds limit 0
Counts A/C/G/T: [[3.0, 1.0, 0.0, 0.0], [1.0, 3.0, 0.0, 0.0]]
PWM A/C/G/T: [[0.5, 0.25, 0.125, 0.125], [0.25, 0.5, 0.125, 0.125]]
PSSM A/C/G/T: [[1.0, 0.0, -1.0, -1.0], [0.0, 1.0, -1.0, -1.0]]
Forward scores: [2.0, -1.0, -2.0, -1.0, 2.0]
Scored hit: Forward 0..2 2
Scored hit: Reverse 2..4 2
Scored hit: Forward 4..6 2
Normalized background: [0.5, 0.25, 0.125, 0.125]
Asymmetric PSSM: [[-1.0, 0.0, 1.0, 1.0]]
Asymmetric hit: Reverse 0..1 1
Zero-probability score for C: -inf
Float64 accumulation, final float32: 0.000000029802322
```

Every reported interval is zero-based, half-open and refers to the original
input, including reverse hits. `Definite` accepts the source R at position 2
because every possibility in R is allowed by the pattern R; it excludes N at
position 4. `Possible` includes N on both strands. `Reject` refuses the first
ambiguous source base at its original offset; the example handles the error
without changing the sequence. Pattern A over source A has two windows under the
default both-strand search, even though only the forward strand matches, so a
zero-window budget refuses the search before scanning; no partial hits are
returned.

The A/C/G/T arrays are **position-major**: one four-value array per motif
column. Their last values are U for RNA. Positive log odds favor the motif
relative to the supplied background; negative values favor that background.
A probability of zero has the meaningful score `-inf`, not a missing value.

## Choose matching semantics

`DnaMotif` and `RnaMotif` are aliases of `IupacMotif<Dna/Rna>`.
`new` validates and owns a nonempty pattern; `pattern()` returns its normalized
text and `width()` its number of positions. `find_in(&sequence)` uses the defaults;
`find_with(&sequence, MotifSearchOptions { ... })` makes policies explicit.

Each IUPAC symbol denotes a set of concrete bases. Source ambiguity has its own
policy, independent of ambiguity in the pattern:

| Policy | Required relation | Example: source N, pattern R |
|---|---|---|
| `Definite` (default) | Source set is a subset of pattern set | No |
| `Possible` | Source and pattern sets intersect | Yes |
| `Reject` | Source must be concrete everywhere | Error |

`StrandSelection::{Forward, Reverse, Both}` defaults to `Both`. Reverse matching
uses the reverse complement of the **pattern** against the original sequence.
All overlaps remain. Results sort by original start, forward before reverse.
A reverse-complement palindrome therefore has two labeled hits at one span;
these labels are not duplicate positions to discard implicitly.

For literal matching without IUPAC set interpretation, use the separate
[sequence operations](sequence-operations-v1.md). An empty IUPAC motif is an
error, unlike that API's empty literal pattern. A source shorter than the motif
has no windows and no hits; reject-ambiguity still validates the whole source.

## Build and score matrices

- `MotifCounts::from_instances(&[Sequence<A>])` requires at least one nonempty,
  equal-width concrete DNA/RNA instance. Gaps and ambiguous instances are refused.
- `MotifCounts::<A>::from_columns(Vec<[f64; 4]>)` accepts nonnegative finite
  weighted counts and unequal column totals, in A/C/G/T or A/C/G/U order.
- `normalize([a,c,g,t_or_u])` adds four nonnegative finite pseudocounts to every
  column and divides by that column's total. `normalize_default()` adds zero.
  A zero-total adjusted column is an error; it is not a uniform distribution.
- `log_odds(background)` normalizes four strictly positive finite background
  weights and computes `ln(probability/background)/ln(2)`.
  `log_odds_uniform()` selects uniform background.
- `columns()` and `width()` inspect each immutable model;
  `MotifScores::background()` returns its normalized background.
- `score_window(&sequence)` requires exactly the matrix width and returns one
  `f32`; `calculate(&sequence, max_windows)` returns all forward scores.
  `find_with(&sequence, MatrixSearchOptions { ... })` returns threshold hits
  with range, strand and `f32` score.

Column/background sums use Neumaier compensation in binary64. Window sums are
sequential binary64 and round once to binary32. A finite representable threshold
also rounds once to binary32 before the inclusive comparison. The example's
cancellation test retains 2^-25; repeatedly rounding partial sums to binary32
would erase it. Returned scores are exact binary32 values, not independently
recomputed binary64 scans.

Reverse scanning complements and reverses the **built PSSM**. It does not rebuild
a PWM against an unchanged asymmetric background. In the example, A scores -1
forward and +1 reverse under background 0.5/0.25/0.125/0.125. This valid one-window
scan is supported. Both scoring and matrix-instance input reject ambiguous
bases explicitly rather than inventing ambiguity probabilities.

Advanced callers can construct `MotifScores::from_columns(columns, background)`
for an already-defined PSSM. Finite scores and negative infinity are accepted;
NaN and positive infinity are errors. This constructor does not rederive scores
from counts.

## Defaults and resource limits

| Option | Default | Applies to |
|---|---:|---|
| `strands` | `both` | Pattern search and score scan |
| `ambiguity` | `definite` | Pattern search only |
| `threshold` | 0 | Score scan only |
| `pseudocounts` | [0,0,0,0] | Matrices and score scans |
| `background` | [0.25,0.25,0.25,0.25] | Matrices and score scans |
| `max_matches` | 100000 | All selected strands together |
| `max_windows` | 10000000 | All selected strands together |
| `max_comparisons` | 100000000 | Windows × motif width, all selected strands |

Search budgets are fields of `MotifSearchOptions` and `MatrixSearchOptions`.
Finite fractional counts and pseudocounts are allowed. Zero pseudocounts can yield
negative infinity; all-zero adjusted columns cannot define a probability model.

The constructors validate matrices but impose no instance-count or model-width
limit; Rust callers own their input allocation. Search through
`find_with` preflights complete worst-case windows × width, even if early
mismatches would do less work. The convenience `calculate` method limits only
its forward output windows; it has no separate comparison or hit-count budget.
Limits are not process-memory guarantees.
A zero hit budget permits only no-hit results. Exceeding any search budget fails
the whole operation; no truncated success or partial hits are returned.
The existing FASTA/decompression reader also applies its documented
[file-input limits](file-io-v1.md).

## Errors and recovery

| Stable code | Meaning / appropriate caller action |
|---|---|
| `invalid_sequence` | Source/instance violates its explicit DNA/RNA alphabet; inspect input and offset |
| `invalid_motif_pattern` | Invalid pattern byte, including T in RNA; correct the chosen molecule or input |
| `empty_motif_pattern` | Empty IUPAC pattern; provide a biological width |
| `ambiguous_motif_input` | Reject mode, count instances or PSSM scan encountered ambiguity; choose another justified analysis, not silent base substitution |
| `invalid_motif_matrix` | Nonfinite/negative count or pseudocount, invalid background, zero total, or invalid score arithmetic |
| `motif_dimension_mismatch` | Unequal instance widths or wrong exact scoring-window length |
| `limit_exceeded` | A search exceeds an explicit budget; inspect the required count before deliberately changing limits |

`MotifError` is typed, implements `Error`, exposes `code()` and retains
variant-specific fields. Reading instance sequences from files reports the
[file-input](file-io-v1.md) errors instead.
A valid ambiguous DNA sequence is not necessarily valid for PSSM scoring.

## Reference compatibility and verification

The [algorithm contract](../algorithms/motifs-v1.md) and
[retained-observation provenance](../../seq-bio/tests/fixtures/motifs-provenance.md)
pin Biopython 1.88, NumPy 2.3.3 and CPython 3.12.3. The standalone core agrees on
5,244 PWM cells, 5,244 PSSM cells, 28,158 forward/reverse scores, 9,127 selected
hits and 18 zero-total errors in the 360-request corpus, under the recorded
platform and arithmetic profile. These counts are not a universal fidelity claim.

Intentional differences are explicit: original-input half-open reverse ranges
rather than Python-negative indices; typed ambiguity/zero-total errors rather
than NaN/ZeroDivisionError; and valid one-window scans rather than the observed
native zero-dimensional-array search failure. RNA uses the same set/matrix
semantics with molecule-specific input validation.

The documentation-driven tests in `seq-bio/tests/motif_recipes.rs` check the
complete program fence and the real `cargo bio motif_analysis` alias against the
output above. Earlier standalone and staging results remain separately scoped
numerical/integration checkpoints, not substitutes for these checks.

Stage 3 is accepted for this documented scope in the
main validation record, which
records toolkit/optimized checks, complete workspace harness output and its
terminal audit, documentation and response comparisons. Final source/evidence
review and the local commit accompany this acceptance. Neither that acceptance nor
the retained numerical corpus implies universal parity or a new scientific
qualification.
