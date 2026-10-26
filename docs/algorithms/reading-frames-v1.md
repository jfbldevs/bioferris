# Six-frame translation and ORF scanning specification

## Contents

- [Scope and references](#scope-and-references)
- [Six-frame translation](#six-frame-translation)
- [ORF scanning](#orf-scanning)
- [Coordinates and ordering](#coordinates-and-ordering)
- [Acceptance criteria](#acceptance-criteria)

## Scope and references

This additive increment extends typed DNA/RNA, not alignment kernels. Biological
assignments and initiation sets reuse the NCBI codes of the shared translation
engine: codes 1, 2 and 11 when it was written, every current NCBI table since
2026-10-22. ORF discovery refuses tables 27, 28 and 31, whose codons that can end
translation also code amino acids. No new dependency, runtime or automatic
annotation interpretation is introduced.

Public references consulted 2026-10-16:

- [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi):
  assignments and initiation sets for the three supported codes.
- [Biopython 1.88 ORF cookbook](https://biopython.org/docs/1.88/Tutorial/chapter_cookbook.html#identifying-open-reading-frames):
  six-frame stop-delimited candidate regions and strand-aware coordinates.

ORF policies are an explicit BioFerris contract, not a claim that Biopython has
one universal ORF caller. Public descriptions guide independent Rust code;
reference implementation source is not copied. Candidates are not gene predictions.

## Six-frame translation

`ReadingFramesExt::translate_six_frames()` returns `[FrameTranslation; 6]` for
table 1; `translate_six_frames_with(GeneticCode)` selects table 1, 2 or 11.
Both methods are infallible on a validated sequence (apart from ordinary Rust
allocation failure). Each result exposes `strand()`, `reading_frame()`,
`genetic_code()`, `range()`, `peptide()`, `skipped_leading_bases()` and
`ignored_trailing_bases()`.

Translate complete codons in all three offsets of the input and its reverse
complement. Retain stop markers; do not recode initiators. Use the existing IUPAC
residue resolution, including B/Z/J/X. Incomplete trailing codons are explicitly
excluded and their count reported. This new workflow has that fixed policy;
existing single-frame translation still rejects partial codons by default.
Short/empty inputs still produce six entries. The actual leading count is
`min(frame_offset, input_length)`; a frame beyond the length has an empty range,
empty peptide and zero trailing bases. The selected frame offset remains 0/1/2.

## ORF scanning

`ReadingFramesExt::find_orfs()` uses `OrfOptions::default()`;
`find_orfs_with(OrfOptions)` returns `Result<Vec<OpenReadingFrame>, OrfError>`.
The result exposes `strand()`, `reading_frame()`, `genetic_code()`, `range()`,
`stop_range()`, `peptide()` and checked `to_protein()`.

| Option | Choices | Default |
|---|---|---|
| `genetic_code` | Existing `GeneticCode` | Table 1 |
| `strands` | `StrandSelection::{Both, Forward, Reverse}` | Both |
| `start_policy` | `OrfStartPolicy::{Atg, GeneticCode, AnySense}` | Atg |
| `stop_policy` | `OrfStopPolicy::{Require, AllowPartial}` | Require |
| `min_amino_acids` | Positive `usize`, excluding stop | 1 |
| `max_orfs` | Maximum retained results; zero permits only an empty result | 100,000 |

The v1 scanner accepts only concrete A/C/G/T DNA or A/C/G/U RNA. It rejects the
first ambiguous symbol in the original input with its byte and zero-based
offset, even in a skipped or incomplete region. Six-frame translation does not
have this restriction. Invalid symbols/gaps are rejected by sequence construction.

Scan each selected strand's three frames in translation direction. A stop resets
the current candidate. Start-based modes retain the earliest accepted start in
each stop-delimited interval, not every nested start. `Atg` accepts ATG/AUG;
`GeneticCode` accepts the selected NCBI initiation set. Both recode the first
residue to M. `AnySense` begins at the first complete non-stop codon after the
frame boundary or preceding stop, and never recodes that residue.

`Require` emits only candidates with a terminal stop. `AllowPartial` also emits
the last candidate at the sequence edge, using only its complete codons and a
missing `stop_range()`. It does not infer a missing start in start-based modes.
Empty candidates are never emitted. Minimum length is inclusive. One/two trailing
bases never become a residue or part of an ORF range. `min_amino_acids == 0` is
`OrfError::InvalidMinimumLength`; excess output is `ResultLimitExceeded { limit }`,
not truncated success. Ambiguity is `AmbiguousBase { offset, byte }`. Errors
return no partial result. Option validation precedes ambiguity validation.

## Coordinates and ordering

`Strand::{Forward, Reverse}` names the orientation. Reading frames are offsets
0/1/2 into that oriented sequence, not signed numbers or genomic modulo values.
All public ranges are zero-based, half-open, ascending coordinates into the
**original input**. Map oriented reverse interval `s, e)` in length `n` to
`[n-e, n-s)`. An ORF's main range excludes its stop codon; `stop_range()` stores
that separate three-base interval. Peptides always read in translation direction.

Order is forward frames 0/1/2, then reverse frames 0/1/2. Within a frame, candidates
follow translation direction (therefore decreasing genomic positions on reverse).
Strand selection filters this order; minimum length filters without reordering.

For `ATGAAATAA`, the forward frame-0 ORF is range `0..6`, stop `6..9`, peptide
`MK`. For `TTATTTCAT`, its reverse counterpart is range `3..9`, stop `0..3`,
peptide `MK`. These are specification expectations, verified by acceptance tests.

Work and retained peptide storage are linear in input length plus result count:
no all-nested-start enumeration, circular wrapping or repeated suffix translation.
`max_orfs` bounds the result count, not exact process memory. Rust callers own
input-size control.

## Acceptance criteria

- Independently enumerate concrete codons for codes 1/2/11 and verify strand,
  frame, stop, initiation, short-input and boundary behavior.
- Compare scanning against a separately structured test oracle and verify DNA/RNA
  equivalence, coordinate extraction and deterministic ordering.
- Verify selected public real-record inputs without scientific campaign reuse.
- Keep existing translation and alignment behavior unchanged.
- Verify real CLI/strict JSON requests, errors, discovery and literal guide output
  (the CLI was removed later ([decision 029)).
- Complete formatting, strict Clippy/docs, Rust regression and release checks;
  record actual results and limits before claiming implementation complete.
