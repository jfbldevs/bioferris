# Motif matrix reference fixture provenance

## Contents

- [Purpose and boundaries](#purpose-and-boundaries)
- [Pinned reference environment](#pinned-reference-environment)
- [Input generation and captured values](#input-generation-and-captured-values)
- [Observed numerical rules](#observed-numerical-rules)
- [Reproduction](#reproduction)
- [Capture program](#capture-program)
- [Verification and limitations](#verification-and-limitations)

[API index](../../../docs/README.md) ·
[Motif contract](../../../docs/algorithms/motifs-v1.md)

## Purpose and boundaries

`motifs-biopython-1.88.json` contains permitted public-API observations for the
stage 3 nucleotide matrix workflow. Inputs are independently constructed synthetic
counts and concrete DNA, not reference implementation source or a frozen
scientific corpus. No reference source was inspected or copied. Biopython and
NumPy are observation tools only; BioFerris computations and tests use Rust.

These observations qualify a declared numerical profile. They do not establish
universal equivalence, gene regulation predictions, transcription-factor binding
accuracy, protein-motif support or acceptance of the complete stage 3 feature.
At capture time the Rust module was compiled and tested through a separate
temporary harness, pending the parent increment's integration and verification.

## Pinned reference environment

Observed on 2026-10-17:

| Component | Value |
|---|---|
| Biopython | 1.88 |
| NumPy | 2.3.3 |
| CPython | 3.12.3 |
| Platform | Linux-6.18.33.2-microsoft-standard-WSL2-x86_64-with-glibc2.39 |
| Byte order | Little endian |
| Observation environment | `OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`; isolated `PYTHONPATH`; `python3 -s` |
| Rust checker | rustc 1.85.0 (4d91de4e4 2025-02-17), edition 2024 |

The two pre-existing isolated wheels were hashed directly:

```text
9384c30460832926e765319040d732d1252a713a05af1dd4cc8b19cb1d590439  biopython-1.88-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl
d9192da52b9745f7f0766531dcfa978b7763916f158bb63bdb8a1eca0068ab20  numpy-2.3.3-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl
```

No global installation, runtime upgrade or dependency change was performed for
this capture.

## Input generation and captured values

The deterministic Cartesian product contains **360 matrix requests**:

- Widths: 1, 2, 3, 4, 5 and 8.
- Counts: positive integers, fractional counts, sparse counts, zero-total columns
  and mixed decimal scales.
- Pseudocounts: scalar 0, scalar 0.5, per-base 0.1/0.2/0.3/0.4 and per-base
  1/2/3/4.
- Background weights: uniform 1/1/1/1, symmetric 0.3/0.2/0.2/0.3 and asymmetric
  scaled weights 3/1/2/4.

Every successful model scans the same literal 44-base sequence:

```text
ACGTACGTTGCAAAAACCCCGGGGTTTTACGTCAGTGCATACGT
```

Case order cycles through threshold zero, an observed forward score, the next
binary64 number above that score and the next binary32 number above it. A model
without a finite forward score uses zero as the anchor.

The file records original inputs alongside their exact binary64 bits, normalized
probability and score matrix bits, and all forward and reverse-complement PSSM
scan scores as binary32 bits. Hexadecimal strings are IEEE 754 big-endian
representations, independent of host byte order. Input bit strings avoid
mistaking a JSON parser's decimal rounding for an algorithm discrepancy.
Negative infinity is retained as `fff0000000000000` in binary64 or `ff800000`
in binary32; no invalid JSON infinities or silently substituted nulls occur.

Native search observations preserve each returned Python position. Reverse
positions map to `original_start = input_length + native_position`; the width
determines the half-open end. Hits are then sorted by original start and forward
before reverse, matching BioFerris's explicit ordering.

Of the 360 requests, 342 succeed and 18 raise native `ZeroDivisionError` for
zero-total columns with zero pseudocounts. BioFerris returns a typed invalid-matrix
error for those cases. Four additional direct-PSSM probes cover cancellation,
one-window scanning and two threshold-rounding boundaries.

## Observed numerical rules

The public matrix API is the count/probability/score functional reference:
[Biopython matrix API](https://biopython.org/docs/latest/api/Bio.motifs.matrix.html)
(the page identifies version 1.88). The versioned
[motif tutorial](https://biopython.org/docs/1.88/Tutorial/chapter_motifs.html)
documents pseudocounts, log-odds scoring, float32 scans and reverse positions.
The versioned matrix API URL was unavailable through the browser during this
capture; the resolved public API page and the installed package both identify
1.88.

Public observations established the following arithmetic profile:

- Column and background totals use Neumaier compensated summation. Python 3.12's
  public [release notes](https://docs.python.org/3.12/whatsnew/3.12.html#other-language-changes)
  identify that algorithm for floating-point `sum`.
- Log-odds evaluation is `ln(probability/background)/ln(2)`; using a direct
  `log2` can differ in the final binary64 bit.
- A scan adds the selected binary64 cells sequentially, then rounds once to
  binary32.
- Native threshold comparisons round the threshold to binary32 and are inclusive.
- Reverse scanning complements the already-built PSSM. Rebuilding a probability
  model against an unchanged asymmetric background is a different computation.

An initial Rust comparison found 615 binary64 probability/score-cell differences,
while every binary32 scan and selected native hit already matched. A public
arithmetic probe isolated the normalization issue: for adjusted counts
0.6/0.825/1.05/0.9, sequential summation gives
`0x1.affffffffffffp+1`, whereas the reference's denominator is
`0x1.b000000000000p+1`. Independently implementing Neumaier summation removed
all 615 differences without changing the captured expected values.

The fixture also preserves an intentional difference: native `search` fails for
exactly one scoring window under NumPy 2.3.3 because its calculation is
zero-dimensional. BioFerris accepts the valid one-window input. The captured
native error is evidence of that difference, not an expected BioFerris failure.

## Reproduction

The full observation program follows. Save that exact code as `capture.py`.
Use an already prepared isolated site containing the pinned wheel contents;
this command does not install packages. Set `REFERENCE_SITE` to that site and
write to a new temporary output, rather than overwriting accepted expectations:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$REFERENCE_SITE" \
  python3 -s capture.py > recaptured.json
cmp seq-bio/tests/fixtures/motifs-biopython-1.88.json recaptured.json
sha256sum capture.py recaptured.json
```

On the recorded platform, a second independent execution was byte-identical.
The executed program and fixture hashes are:

```text
dbdbc364307a014af25ccb2a1bbced0955107dee1a4d1cb5b87baba93a9ae411  capture.py
75e37c85cd6370b2b17a9773e5226a459071f2156d6b2a3285d7321f25f5c824  motifs-biopython-1.88.json
```

The fixture contains platform metadata, so another platform can legitimately
produce a different file even before numerical comparisons. Investigate and
record differences rather than automatically refreshing expectations.

## Capture program

```python
"""Observe public reference APIs; never inspect their implementation source."""

import itertools
import json
import math
import platform
import struct
import sys

import Bio
import numpy as np
from Bio.Seq import Seq
from Bio.motifs.matrix import FrequencyPositionMatrix, PositionSpecificScoringMatrix

ALPHABET = "ACGT"
WIDTHS = [1, 2, 3, 4, 5, 8]
SOURCE = "ACGTACGTTGCAAAAACCCCGGGGTTTTACGTCAGTGCATACGT"
PSEUDOCOUNTS = [
    ("zero", 0.0),
    ("half", 0.5),
    ("fractional", {base: (i + 1) / 10 for i, base in enumerate(ALPHABET)}),
    ("integer", {base: i + 1.0 for i, base in enumerate(ALPHABET)}),
]
BACKGROUNDS = [
    ("uniform", [1.0, 1.0, 1.0, 1.0]),
    ("symmetric_gc", [0.3, 0.2, 0.2, 0.3]),
    ("asymmetric_scaled", [3.0, 1.0, 2.0, 4.0]),
]
FAMILIES = ["integer", "fractional", "sparse", "zero_column", "mixed_scale"]


def bits64(value):
    return struct.pack(">d", float(value)).hex()


def bits32(value):
    return struct.pack(">f", float(value)).hex()


def row_values(columns):
    return {base: [column[i] for column in columns] for i, base in enumerate(ALPHABET)}


def column_bits(matrix, width):
    return [[bits64(matrix[base][j]) for base in ALPHABET] for j in range(width)]


def count_columns(family, width):
    result = []
    for j in range(width):
        if family == "integer":
            column = [float(((j + 1) * (b + 2) + j * b) % 11 + 1) for b in range(4)]
        elif family == "fractional":
            column = [(j + 1) * (b + 1) / 10 + ((j + b) % 3) / 8 for b in range(4)]
        elif family == "sparse":
            column = [float(j + 1) if b == j % 4 else 0.0 for b in range(4)]
        elif family == "zero_column":
            column = [0.0 if j % 3 == 0 else float(j + b + 1) for b in range(4)]
        elif family == "mixed_scale":
            column = [((j + b) % 7 + 1) * 10.0 ** ((j * 3 + b) % 5 - 2) for b in range(4)]
        else:
            raise AssertionError(family)
        result.append(column)
    return result


def original_hits(pssm, sequence, threshold):
    hits = []
    for position, score in pssm.search(Seq(sequence), threshold=threshold, both=True):
        position = int(position)
        start = position if position >= 0 else len(sequence) + position
        hits.append({
            "native_position": position,
            "start": start,
            "end": start + pssm.length,
            "strand": "forward" if position >= 0 else "reverse",
            "score_f32_bits": bits32(score),
        })
    return sorted(hits, key=lambda hit: (hit["start"], hit["strand"] != "forward"))


def capture_matrix(width, family, pc_name, pc, bg_name, background, index):
    columns = count_columns(family, width)
    per_base_pc = [float(pc[base]) for base in ALPHABET] if isinstance(pc, dict) else [pc] * 4
    case = {
        "id": f"w{width}-{family}-{pc_name}-{bg_name}",
        "width": width,
        "counts": columns,
        "counts_f64_bits": [[bits64(value) for value in column] for column in columns],
        "pseudocounts": per_base_pc,
        "pseudocounts_f64_bits": [bits64(value) for value in per_base_pc],
        "background": background,
        "background_f64_bits": [bits64(value) for value in background],
        "sequence": SOURCE,
    }
    try:
        pwm = FrequencyPositionMatrix(ALPHABET, row_values(columns)).normalize(pseudocounts=pc)
        pssm = pwm.log_odds(dict(zip(ALPHABET, background)))
    except Exception as error:
        case["reference_error"] = {"type": type(error).__name__, "message": str(error)}
        return case
    forward = np.asarray(pssm.calculate(Seq(SOURCE)), dtype=np.float32).reshape(-1)
    reverse = np.asarray(pssm.reverse_complement().calculate(Seq(SOURCE)), dtype=np.float32).reshape(-1)
    finite_scores = [float(value) for value in forward if math.isfinite(float(value))]
    anchor = finite_scores[len(finite_scores) // 2] if finite_scores else 0.0
    threshold_mode = index % 4
    if threshold_mode == 0:
        threshold, threshold_kind = 0.0, "zero"
    elif threshold_mode == 1:
        threshold, threshold_kind = anchor, "observed_score"
    elif threshold_mode == 2:
        threshold, threshold_kind = math.nextafter(anchor, math.inf), "next_f64_after_observed_score"
    else:
        threshold = float(np.nextafter(np.float32(anchor), np.float32(math.inf)))
        threshold_kind = "next_f32_after_observed_score"
    case.update({
        "weights_f64_bits": column_bits(pwm, width),
        "scores_f64_bits": column_bits(pssm, width),
        "forward_scores_f32_bits": [bits32(value) for value in forward],
        "reverse_scores_f32_bits": [bits32(value) for value in reverse],
        "threshold": threshold,
        "threshold_f64_bits": bits64(threshold),
        "threshold_kind": threshold_kind,
        "native_hits": original_hits(pssm, SOURCE, threshold),
    })
    return case


def public_probe(columns, sequence, threshold):
    pssm = PositionSpecificScoringMatrix(ALPHABET, row_values(columns))
    result = {
        "columns": columns,
        "sequence": sequence,
        "threshold": threshold,
        "forward_scores_f32_bits": [bits32(value) for value in np.asarray(pssm.calculate(Seq(sequence))).reshape(-1)],
    }
    try:
        result["native_hits"] = original_hits(pssm, sequence, threshold)
    except Exception as error:
        result["reference_search_error"] = {"type": type(error).__name__, "message": str(error)}
    return result


cases = []
for width, family, (pc_name, pc), (bg_name, background) in itertools.product(WIDTHS, FAMILIES, PSEUDOCOUNTS, BACKGROUNDS):
    cases.append(capture_matrix(width, family, pc_name, pc, bg_name, background, len(cases)))

single = float(np.float32(0.1))
probes = {
    "double_accumulation": public_probe([[1.0] * 4, [2.0 ** -25] * 4, [-1.0] * 4], "AAAA", 0.0),
    "one_window_native_search": public_probe([[0.1] * 4], "A", 0.0),
    "threshold_above_score_but_same_float32": public_probe([[0.1] * 4], "AA", single + 1e-12),
    "threshold_next_float32": public_probe([[0.1] * 4], "AA", float(np.nextafter(np.float32(single), np.float32(math.inf)))),
}
document = {
    "schema_version": 1,
    "profile": "biopython-1.88-motif-matrices-v1",
    "reference": {
        "biopython": Bio.__version__,
        "numpy": np.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "byteorder": sys.byteorder,
        "biopython_wheel_sha256": "9384c30460832926e765319040d732d1252a713a05af1dd4cc8b19cb1d590439",
        "numpy_wheel_sha256": "d9192da52b9745f7f0766531dcfa978b7763916f158bb63bdb8a1eca0068ab20",
    },
    "generation": {
        "widths": WIDTHS,
        "count_families": FAMILIES,
        "pseudocount_profiles": [name for name, _ in PSEUDOCOUNTS],
        "background_profiles": [name for name, _ in BACKGROUNDS],
        "threshold_modes": ["zero", "observed_score", "next_f64_after_observed_score", "next_f32_after_observed_score"],
        "matrix_order": "column-major, A/C/G/T within each column",
        "binary_encoding": "IEEE 754 big-endian hexadecimal; 16 digits for binary64, 8 digits for binary32",
        "native_hit_order": "sorted original start then forward before reverse; native_position preserved",
    },
    "cases": cases,
    "probes": probes,
}
json.dump(document, sys.stdout, indent=2, allow_nan=False)
sys.stdout.write("\n")
```

## Verification and limitations

The independently authored Rust checker read the stored input bits and compared
actual results against the captured output bits. Its final report was:

```text
cases=360 expected_errors=18 weight_cells=5244 score_cells=5244 scan_scores=28158 native_hits=9127 matrix_mismatches=0 scan_mismatches=0 hit_mismatches=0
```

The permanent internal Rust regressions additionally exercise the four direct
numerical probes and the intentional one-window difference. All 18 module tests
passed in the standalone harness, along with formatting and strict Clippy
`all`/`pedantic` checks. The remaining tests cover independently specified IUPAC
sets, reverse complements, overlaps, coordinates, RNA, empty/short inputs,
invalid values and resource limits. No Cargo target or shared module export was
changed by this disjoint draft.

This synthetic corpus is finite. Exact bit agreement is established only for its
recorded environment, inputs and operations. It does not prove matching across
platform math libraries, interpreter versions, every possible float, unsupported
backgrounds, ambiguous matrix scan input, score distributions, motif file
formats or database workflows. Source ambiguity matching has its own explicit
set semantics and independent Rust tests; it is not inferred from PSSM scans.
