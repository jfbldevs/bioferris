# Annotated CDS validation fixture provenance

## Contents

- [Public DNA record](#public-dna-record)
- [Independent biological expectations](#independent-biological-expectations)
- [Synthetic interface fixtures](#synthetic-interface-fixtures)
- [Biopython behavioral snapshot](#biopython-behavioral-snapshot)
- [Offline verification](#offline-verification)

[API index](../../../docs/README.md) · [CDS contract](../../../docs/algorithms/cds-translation-v1.md)

## Public DNA record

`annotated-nc_001422.1.gb` is the complete unmodified NCBI GenBank response for
RefSeq accession **NC_001422.1**, Escherichia phage phiX174, downloaded 2026-10-17:

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=gb&retmode=text

The LOCUS declares 5,386 bases, ss-DNA and circular topology. Metadata, feature
qualifiers, sequence case/wrapping and terminal newline are retained. No molecule
relabeling, T/U replacement, annotation rewriting or coordinate conversion was
applied to the response. This is public biological data, not implementation source.
The accession version identifies the sequence; retained bytes and their checksum
also pin the metadata snapshot.
The narrow Git whitespace exemption preserves NCBI's padded `ORIGIN` line and
terminal blank line; these primary-source bytes are not reformatted for linting.

```text
0a8099e207633f8c3ba998bfe2a198518e1de356c7956691101b2512324b807c  annotated-nc_001422.1.gb
```

## Independent biological expectations

The source's E/endolysin CDS is **568..843** in one-based inclusive file
coordinates, protein identifier **NP_040709.1**, 276 bases including the terminal
stop and 91 translated residues. It uses `/codon_start=1` and `/transl_table=11`.
Its nucleotide sequence is independently available in the existing
`phix174-e.fasta` snapshot downloaded on 2026-10-15. Tests can compare extraction
against that separate response and the computed protein against the retained
source `/translation`, without generating an expected protein from BioFerris.

The full genome also contains explicitly joined origin-spanning CDS locations.
Those joins are interpreted in their written biological order; circular metadata
does not cause implicit wrapping. Every accepted complete CDS must match its own
source `/translation` under the declared table. This fixture does not establish
universal annotation interpretation or live Biopython differential equivalence.

The earlier DQ067560.1 fuzzy-CDS and X56734.1 RNA/T-storage snapshots remain
unchanged. Their documented refusals remain part of the supported-profile tests;
they are not relabeled or converted to manufacture positive cases.

## Synthetic interface fixtures

`cds-request.json` requests checked translation of the existing synthetic
`annotated-example.gb`; that CDS yields `MA` and has the matching `/translation`.
`cds-mismatch.gb` deliberately labels the same computed `MA` protein as `MV` to
exercise exact mismatch errors and explicit reference-ignore recovery. It is a
local synthetic record, not a database accession or scientific qualification case.

## Biopython behavioral snapshot

`cds-biopython-1.88.json` records direct observed output from **Biopython 1.88**
on CPython **3.12.3**, Linux x86-64, collected 2026-10-17. The unchanged complete
DNA snapshot above was parsed with the public `Bio.GenBank.FeatureParser` API;
each `CDS` was translated using `feature.translate(record.seq)` with its source
qualifiers and default CDS checks. All 11 peptides exactly matched their source
reference proteins. The retained fields are source feature ordinal, record ID,
reference location representation, selected table/codon start, peptide and exact
source-protein match. Rust tests compare the product/table/frame/ordinal fields;
the reference's location display is not the BioFerris source-location format.

The official wheel was downloaded using [PyPI release metadata](https://pypi.org/pypi/biopython/1.88/json),
verified, and unpacked into a task-specific temporary directory. Only public API
behavior was observed; implementation source was not inspected or copied.
There was no global package installation, interpreter upgrade or new BioFerris
runtime dependency. The reference needs no access to private project data.

```text
Wheel: biopython-1.88-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl
SHA256: 9384c30460832926e765319040d732d1252a713a05af1dd4cc8b19cb1d590439
Source: https://files.pythonhosted.org/packages/17/c2/613ca5f729b55b2f7a9c17a3cd59d3f053ee2816da2d19933acc6b8741ec/biopython-1.88-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl
Snapshot SHA256: 0c493b52a3b5a6f809fd982b6df8f30c3cf3b31e22925aa39111c3aeaaae5420
```

An initial system import confirmed Biopython was absent. An initial attempt to
use the higher-level `SeqIO` entry point in the isolated package failed because
that entry point imports optional NumPy-dependent alignment code. Switching to
the public GenBank parser avoided that unrelated optional dependency; no NumPy
package was installed and no biological expectations were changed.

## Offline verification

Core, CLI and documentation-recipe tests consume checked-in fixtures offline.
No test contacts NCBI or requires an interpreter. Both snapshot checksums are
checked explicitly by Rust tests and the new validation record. The documented
feature subset, commands and observed results are maintained alongside the code.
