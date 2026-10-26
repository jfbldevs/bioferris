# Biopython 1.87 pairwise fixtures

This directory contains inert black-box observations for the retained pairwise
portfolio:

- `needleman-wunsch.json` — exact global linear alignment;
- `smith-waterman.json` — exact local linear alignment;
- `gotoh.json` — exact global affine alignment;
- `affine-modes.json` — Gotoh local and free-end modes.

The fixtures were generated outside this workspace with
`Bio.Align.PairwiseAligner` from Biopython 1.87 on CPython 3.13.14. Each JSON file
records its own case and generator hashes. The external generator independently
replays reference paths; Rust comparisons require equal optima and valid replay but
do not require identical tie-breaking.

Biopython, Python, and NumPy are not build or runtime dependencies. Cargo consumes
only these static JSON observations.
