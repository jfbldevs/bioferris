# Preserved analysis/interoperability preparation snapshots

## Source index

| Preserved preparation | Current compiled target |
|---|---|
| `motif_analysis.rs` | [Motif example](../examples/motif_analysis.rs) |
| `alignment_analysis.rs` | [Alignment-analysis example](../examples/alignment_analysis.rs) |
| `blast_analysis.rs` | [BLAST example](../examples/blast_analysis.rs) |
| `compressed_io.rs` | [Compressed-I/O example](../examples/compressed_io.rs) |
| `motif_commands.rs` | Motif CLI tests; removed with the CLI |
| `alignment_analysis_commands.rs` | Alignment-analysis CLI tests; removed with the CLI |
| `motifs.rs` | [Public motif tests](../tests/motifs.rs) |
| `alignment_analysis_public.rs` | [Public alignment-analysis tests](../tests/alignment_analysis_public.rs) |
| `motif_recipes.rs` | [Motif documentation checks](../tests/motif_recipes.rs) |
| `alignment_analysis_recipes.rs` | [Alignment documentation checks](../tests/alignment_analysis_recipes.rs) |
| `blast_commands.rs` | BLAST CLI tests; removed with the CLI |
| `blast_recipes.rs` | [BLAST documentation checks](../tests/blast_recipes.rs) |
| `compression_cli.rs` | Compression CLI tests; removed with the CLI |
| `compression_recipes.rs` | [Compression documentation checks](../tests/compression_recipes.rs) |
| `analysis_request_boundaries.rs` | Cross-feature request/discovery tests; removed with the CLI |

## Purpose and limits

These files preserve the exact pre-promotion sources pinned by the earlier
preparation manifest.
They are historical evidence, not Cargo targets, alternate public implementations
or current usage instructions. Production examples/tests live at the linked
destinations. Main integration later strengthened the discovery checks in
`tests/analysis_request_boundaries.rs`, so that file was intentionally no longer
byte-identical to this earlier snapshot. The five CLI and request test targets
were deleted with the `bioferris` CLI (decision 029); their snapshots
here stay unchanged as historical evidence.

The current validation record
distinguishes standalone/staging checks from final main-checkout gates. Do not
edit these retained snapshots to match later production changes or claim their
old test results validate new source. No external runtime or frozen scientific
campaign is required to use the compiled targets.
