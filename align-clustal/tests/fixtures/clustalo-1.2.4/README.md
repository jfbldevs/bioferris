# Clustal Omega 1.2.4 observations

`default.json` is a static black-box observation generated in the isolated sibling
laboratory `rust-sequence-analysis-reference-validation/clustalo`.

- Reference: Ubuntu `clustalo` 1.2.4 (`1.2.4-8build1`)
- Executable SHA-256: `e695ae9d1fdc9b27540eaf041740b3dafc6e38670cbadee8afe3f6d9114062ea`
- Settings: explicit molecule type, full distance matrix, input output order, one
  thread, no iteration
- Observation date: 2026-08-03

The executable is never called by Cargo. The fixture contains only inputs and
observable outputs; no reference implementation source was inspected or copied.
The sibling generator's `--check` mode reproduces `default.json` byte-for-byte.
