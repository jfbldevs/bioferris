# Edlib 1.3.9.post1 Fixtures

## Myers bit-vector profile

`myers.json` is static black-box evidence generated outside this Rust workspace by
the Python API of Edlib 1.3.9.post1. It is test data, not implementation source or a
build, runtime, or distribution dependency.

Recorded environment and configuration:

- CPython 3.13.14 on Windows 11 build 10.0.26200;
- Edlib 1.3.9.post1, installed from the PyPI source distribution;
- global API: `edlib.align(pattern, text, mode="NW", task="path", k=-1)`;
- search API: reversed pattern and reversed text prefix with `mode="SHW"`,
  `task="distance"`, and `k=-1`;
- literal equality with no additional equalities;
- fixture SHA-256:
  `bde4c1abfa84444c94a003c1ea534709d30a4cadab62c6ce1d4ebf2b83acd4c3`;
- generator SHA-256:
  `de8b72921a14cb6d256a602f845db0d4e84ea3ceb6db3b5c24271dc698efd507`;
- case specification SHA-256:
  `0770b14f9b2366b736583e658d92788da70fe851f41850b6ab0f98b0f5189e49`;
- PyPI source-distribution SHA-256:
  `b0fb6e85882cab02208ccd6daa46f80cb9ff1d05764e91bf22920a01d7a6fbfa`.

Reversing both the pattern and a text prefix maps suffixes ending at that prefix's
last symbol to prefixes of the reversed text. Edlib's `SHW` boundary makes the
remaining reversed-text suffix free, so its distance equals Myers's required
per-endpoint approximate-suffix-search score.

The suite contains 20 DNA, RNA, and protein cases. It covers empty inputs, identity,
all unit edit types, initial/internal/terminal gaps, completely different and
repetitive strings, literal ambiguity symbols, pattern word boundaries 63/64/65 and
127/128/129, and the curated UniProt P02185/P68871 protein pair. Every reference path
is replayed by the external generator before export. Rust must reproduce every
global distance, every per-endpoint search score, and a separately verified optimal
path; it need not duplicate Edlib's tie-breaking.

The external generator and pinned environment are maintained in the sibling
`rust-sequence-analysis-reference-validation` laboratory. Cargo never invokes them.
See `docs/validation/black-box.md` for the isolation and regeneration contract.
