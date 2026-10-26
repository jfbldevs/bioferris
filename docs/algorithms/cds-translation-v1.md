# Annotation-aware CDS translation contract

## Contents

- [Purpose and references](#purpose-and-references)
- [Feature interpretation](#feature-interpretation)
- [Reference comparison](#reference-comparison)
- [Errors and limits](#errors-and-limits)
- [Interfaces](#interfaces)

[API index](../README.md) · Implementation plan

## Purpose and references

This contract defines stage 1 of the authorized six-proposal plan. Implementation
and verification are in progress; this contract alone is not completion evidence.
It connects existing exact local feature extraction to the existing complete-CDS
translation engine. Public references, consulted 2026-10-17:

- [Biopython 1.88 feature translation](https://biopython.org/docs/1.88/api/Bio.SeqFeature.html#Bio.SeqFeature.SeqFeature.translate).
- [INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/).
- [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

No reference implementation source is used. The accepted concrete/IUPAC codon
resolver and the NCBI tables remain the shared computational engine; tables other
than 1, 2 and 11 and `/ribosomal_slippage` were added on 2026-10-22
(decision 033).

## Feature interpretation

Only the case-sensitive feature key `CDS` is accepted. Extraction uses exact
local points/ranges, ordered `join` expressions and `complement`, including the
existing explicitly joined circular-origin representation. Unsupported locations
retain the existing extraction error; topology never creates implicit wrapping.

After extraction into biological coding orientation:

- `/codon_start` is one ASCII decimal integer 1, 2 or 3, default 1. Its zero-based
  offset applies to the extracted sequence, after all joins/complements.
- `/transl_table` is one ASCII decimal supported NCBI table ID. It takes precedence
  over the caller's `default_genetic_code`, whose default is table 1.
- Duplicate, missing-value, nondecimal and out-of-range interpretation qualifiers
  are errors. Surrounding ASCII whitespace is permitted; internal whitespace is not.
- `/pseudo`, `/pseudogene`, `/transl_except`, `/exception` and `/trans_splicing`
  always fail explicitly. They are not bypassed by disabling reference
  comparison. Other qualifiers remain metadata.
- `/ribosomal_slippage` must be a flag, at most once. The location is extracted
  and translated exactly as annotated: a join usually encodes the frameshift,
  with overlapping parts repeating a base and separated parts skipping bases. The
  result reports the flag. A product read without the shift it needs fails the
  complete-CDS rules or the reference check.
- The existing complete-CDS rules require a valid table-specific initiator,
  complete codons, a last codon that can end translation and no internal definite
  stop. In tables 27, 28 and 31 that last codon may be one that otherwise codes an
  amino acid. The initiator becomes M; the terminal stop is excluded. Existing
  IUPAC resolution applies.

The result retains the extracted feature length, selected table and codon start,
peptide and reference-comparison status. Translation error offsets are relative
to the extracted coding-oriented feature, not genomic coordinates. The source
location remains available separately. Extraction never mutates the parent record.

## Reference comparison

`CdsReferencePolicy` has three explicit choices:

- `CheckIfPresent` (default): compare `/translation` if present; its absence is
  represented by `Absent`.
- `RequireMatch`: require `/translation` and exact peptide equality.
- `Ignore`: do not interpret `/translation`; return `Ignored` and no retained
  reference string. Complete-CDS and exception checks still apply.

When comparison is enabled, require exactly one valued `/translation`, remove
ASCII whitespace introduced by wrapping and require an uppercase, nonempty
stop-free protein string accepted by the existing protein alphabet. No stop,
case, residue or mismatch correction is inferred. Equality includes every residue
and length. A mismatch errors with the reference, computed peptide and zero-based
first differing residue (the shorter length for a strict prefix). A successful
comparison has status `Matched` and retains the normalized reference.

Biopython's documented feature translation provides extraction/table/frame
semantics; automatic comparison and conservative qualifier refusals are additional
BioFerris validation policies, not claims of identical exception behavior.

## Errors and limits

Typed errors distinguish a non-CDS key, unsupported biological qualifier,
duplicate/invalid qualifier, unsupported table, extraction failure, complete-CDS
failure, missing reference and reference mismatch. No partially accepted peptide
is returned on error. Existing annotation parsing/extraction limits apply.

The removed CLI processed whole files: it preserved record and feature source
order, numbered features one-based among all features, aborted on the first
failure and applied a global `max_features` budget of 100,000 CDS features.
Rust callers translate features one at a time and own any such budget.

## Interfaces

Rust adds `Feature::translate_cds` and `Feature::translate_cds_with`, with typed
`CdsTranslationOptions`, `CdsTranslation`, `CdsReferencePolicy`,
`CdsReferenceStatus` and `CdsTranslationError`. The methods accept the parent DNA
or RNA sequence and delegate to existing extraction and translation.

The `translate-cds-features` CLI/JSON adapter over the same computation was
removed with the `bioferris` CLI (decision 029).
Public usage examples and observed output are added after execution.
