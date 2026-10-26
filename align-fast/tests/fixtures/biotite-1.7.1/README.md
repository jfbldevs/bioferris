# Biotite X-drop black-box fixture

`xdrop.json` contains observations from the public
`biotite.sequence.align.align_local_gapped()` API in Biotite 1.7.1
(BSD-3-Clause). The package is installed only in the sibling
`rust-sequence-analysis-reference-validation/xdrop` laboratory from a pinned,
SHA-256-verified Windows CPython 3.13 wheel and transitive lock file.

The generator treats Biotite as an opaque black box. No Biotite source was
inspected, copied, translated, linked, invoked by Cargo, or placed in this
workspace. The checked-in JSON is inert data containing scores, ranges, traces,
rendered alignments, package provenance, and input/generator hashes.

Normative conformance compares every affine seeded-extension score and all
unambiguous consumed ranges. One fixture (`dna_long_affine_gap`) has two endpoints
with score 15: Biotite reports the shorter prefix, while Rust applies its documented
farthest-endpoint tie policy and consumes the complete region. The test requires
that exact, named divergence and rejects any other range difference. Rust
tracebacks are checked by an independent operation replay because equally scoring
paths may legitimately differ. RNA support is covered by native tests; Biotite's
nucleotide sequence type does not accept `U` and is therefore not used as an RNA
oracle.
