# Annotated public-data fixture provenance

## Contents

- [Source snapshots](#source-snapshots)
- [Independent expectations](#independent-expectations)
- [Molecule regression and remaining limitation](#molecule-regression-and-remaining-limitation)
- [Integrity and offline execution](#integrity-and-offline-execution)

[Annotated input guide](../../../docs/api/annotated-input-v1.md)

## Source snapshots

Retrieved 2026-10-15 from primary NCBI/ENA public sequence-data endpoints. Files
contain complete response bodies, including source metadata, line wrapping and
terminal newlines. No molecule relabeling, T/U replacement, annotation trimming
or coordinate editing was performed. The saved bytes were compared directly
with the downloaded response bodies before recording their checksums.

| Fixture | Accession/version | Primary source | Purpose |
|---|---|---|---|
| [annotated-dq067560.1.gb](annotated-dq067560.1.gb) | DQ067560.1 | [NCBI GenBank response](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=DQ067560.1&rettype=gb&retmode=text) | Positive genomic-DNA import and exact exon extraction |
| [annotated-dq067560.1.embl](annotated-dq067560.1.embl) | DQ067560.1 | [ENA EMBL response](https://www.ebi.ac.uk/ena/browser/api/embl/DQ067560.1?download=false) | Same accession in modern EMBL |
| [annotated-x56734.1.gb](annotated-x56734.1.gb) | X56734.1 | [NCBI GenBank response](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X56734.1&rettype=gb&retmode=text) | Explicit RNA/T-storage compatibility gap |
| [annotated-m77789.2.embl](annotated-m77789.2.embl) | M77789.2 | [ENA EMBL response](https://www.ebi.ac.uk/ena/browser/api/embl/M77789.2?download=false) | Positive `other DNA` molecule-profile regression |

These are public biological data, not reference implementation source. The test
does not claim that the two databases expose identical metadata; their shared
sequence and the stated feature facts are compared explicitly. Accession versions
pin sequence identity but do not guarantee that future server metadata will stay
byte-identical, which is why these offline response snapshots are retained.

## Independent expectations

DQ067560.1 is the *Salmo trutta* p53 gene fragment described by the source.
GenBank LOCUS and EMBL ID/SQ declare 534 bases; EMBL SQ independently records
134 A, 126 C, 134 G and 140 T. Both feature tables identify exon 5 as 1..187
and exon 6 as 422..534. The test asserts these source facts without generating
expected values from the production parser.

The complete 113-base exon 6 expectation below was independently transcribed
from sequence row 421, omitting its first base (position 421), followed by row
481. File coordinates are one-based inclusive:

```text
GTCCTGCCCCGCGAGGTCACCTGGTCAGAGTTGAGGGGAACCAGCGAGCAGAGTACATGGAGGATGGTAACACTCTGAGACAGAGTGTGCTCGTCCCCTATGAGCCTCCTCAG
```

The source's CDS is `join(<1..187,422..>534)`, with protein identifier
`AAY68035.1`. Its fuzzy ends are preserved, and extraction returns
`unsupported_format_feature`; the parser does not guess a complete CDS or
translate the record using an implicit genetic code.

## Molecule regression and remaining limitation

- X56734.1 is 1,859 bases, declared mRNA, with T in the stored sequence. A DNA
  request fails the molecule check; an RNA request fails strict RNA symbol
  validation. This confirms a documented compatibility limitation, not invalid
  reference data.
- M77789.2 is the 2,686-base pUC19 cloning vector. ENA declares `other DNA`;
  SQ independently records 659 A, 685 C, 675 G and 667 T. The test checks all
  counts, successful typed DNA import, retained source text and rejection of an
  RNA request. The initial allowlist rejection was fixed narrowly: only the
  exact case-insensitive modern-EMBL `other DNA` token was added. Tests ensure
  unrelated molecule tokens and strand-prefixed `ss-other DNA` remain unsupported.

The mapping follows the [ENA ID-line contract](https://ftp.ebi.ac.uk/pub/databases/ena/sequence/release/doc/usrman.txt),
which uses the INSDC molecule vocabulary, and the [INSDC /mol_type definition](https://www.insdc.org/submitting-standards/feature-table/),
which classifies `other DNA` as a synthetic molecule. These primary specifications
were checked on 2026-10-15. Neither source residues nor fixture bytes were changed;
the refinement adds no RNA/T normalization or arbitrary molecule inference.

## Integrity and offline execution

Four source fixtures total **17,541 bytes**:

```text
204edf8e8bc5dab6098702d4bbcb39f41ada5393665a129e2d70b9844fc3fd72  annotated-dq067560.1.gb
ed019b8b2dd1d6e930c219852aa86a796d077f0c66bd2ac3ae61f17da6255347  annotated-dq067560.1.embl
fdb9e4a45838464eacd8f2694a17ae26d4bd9c6ffbf0335c8e32c0fb6174b3fe  annotated-x56734.1.gb
dbbad17ccf47e48cb16d12630fb22b023c7add9deb4fbccb32957d5c804610b9  annotated-m77789.2.embl
```

Run from the repository root:

```bash
cargo test -p bioferris --test annotated_import --offline --locked
```

Observed result after adding these fixtures: **17 passed, 0 failed, 0 ignored**.
The integration tests use `include_str!` and perform no network requests.
No Python/Biopython process, scientific qualification campaign, performance
measurement or universal compatibility claim is involved.
