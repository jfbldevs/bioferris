# Sequence fixture provenance

## Contents

- [Analysis and interoperability fixtures](#analysis-and-interoperability-fixtures)
- [Existing sequence and record fixtures](#existing-sequence-and-record-fixtures)
- [Real sequence](#real-sequence)
- [NCBI genetic codes](#ncbi-genetic-codes)
- [Streamed GenBank records](#streamed-genbank-records)

## Analysis and interoperability fixtures

The four later analysis/interoperability stages are implemented and verified for
their documented scope; the integration validation
records the checks and evidence limitations. Fixture presence alone is not proof
of acceptance. Automated checks use Rust and retained offline data; reference
packages were observation tools, not application/runtime dependencies.

`motif-instances-dna.fa`, `motif-request.json` and
`motif-recipe-columns-request.json` are synthetic interface examples for the
[complete motif workflows](../../../docs/api/motifs-v1.md): four AC/AC/AC/CA
instances, explicit pseudocounts and a sparse one-column model with `-inf` scores.
`motifs-biopython-1.88.json` instead retains real public-API numerical observations;
its [provenance, precision rules and checksums](motifs-provenance.md) identify
Biopython/NumPy/interpreter versions and the exact corpus/platform scope.

`alignment-analysis.fa`, `alignment-analysis.sto`, `alignment-analysis.nex`
and `alignment-analysis-request.json` are independently authored examples for
[alignment analysis](../../../docs/api/alignment-analysis-v1.md). They retain
all-gap columns, ambiguous residues, custom gap/missing glyphs, reverse origins
and five ordered Stockholm annotations. Their expected policies and complete
loss reports are not new scientific alignment qualification.

`blast-example.tab6`, `blast-example.tab7`, `blast-minimal.tsv`,
`blast-example.xml`, `blast-example-xml2.xml` and `blast-request.json` are
synthetic [BLAST result-reader examples](../../../docs/api/blast-input-v1.md),
not executed searches. The NCBI cookbook excerpt and four unmodified public
NCBI BLAST 2.9.0+ XML reports have separate
[source URLs, transformation notes and checksums](blast-provenance.md).
`blast-biopython-1.88.json` records permitted public-API observations, not invented
expected outputs or copied parser implementation.

`compression-example.fa` and `compression-index-request.json` are synthetic
[compressed-I/O examples](../../../docs/api/compressed-io-v1.md). The request
uses a BGZF source/sidecar created by the preceding create-only recipe.
The public BAM byte fixture, native-generated BGZF streams and retained
`bgzf-*-1.88.json` observations have separate
[provenance and checksums](bgzf-provenance.md). These cover decoded-byte,
block/virtual-position and sequence-index interoperability; compressed-byte
identity is not required, and byte-level BAM coverage is not a BAM parser.

## Existing sequence and record fixtures

`protein-properties-request.json` is a synthetic INGAR analysis request. Its
complete input/result appears in the [protein-property guide](../../../docs/api/protein-properties-v1.md).
The separate `protein-properties-biopython-1.88.json` retains 7,344 observed
reference values and all 432 literal inputs; its [provenance and checksum](protein-properties-provenance.md)
distinguish real public-API observations from synthetic interface examples.

`cds-request.json` and `cds-mismatch.gb` are synthetic annotation-aware translation
examples. The unmodified public `annotated-nc_001422.1.gb` snapshot supplies
complete DNA CDS features, including exact source `/translation` proteins and
origin-spanning joins. See the separate [CDS provenance and checksum](cds-provenance.md).

`orf-request.json` is a synthetic inline request for reverse-strand ORF discovery
on DNA `TTATTTCAT`. The [reading-frame guide](../../../docs/api/reading-frames-v1.md)
shows the complete input and result. Tests separately use the public genetic-code
tables and independent codon/segmented ORF oracles, plus the existing phiX174 and
pUC19 fixtures below; this small request is not scientific qualification data.

`io-example.fastq`, `io-genes.fasta` and `io-request.json` are synthetic, locally
authored file-I/O documentation fixtures. Their complete inputs and observed
outputs appear in [File I/O v1](../../../docs/api/file-io-v1.md) and are checked
by Rust integration tests. They are not reference-fidelity or performance data.

`sequence-request.json`, `statistics-example.fasta` and `statistics-request.json`
are likewise synthetic local fixtures. The [sequence operations](../../../docs/api/sequence-operations-v1.md)
and [statistics](../../../docs/api/sequence-statistics-v1.md) guides show their
complete inputs, runnable commands and observed outputs. The statistics FASTA
deliberately uses unequal record lengths and ambiguous symbols to distinguish
count-first GC aggregation from averaging fractions or joining record boundaries.

`configured-translation-input.fasta` and `configured-translation-request.json`
are synthetic local examples for [explicit genetic codes and frames](../../../docs/api/genetic-codes-v1.md).
The two FASTA records select table 11/frame 1, exercise GTG/CTG initiation and
produce MA/MP in CDS mode. The request selects the first record inline. These
fixtures demonstrate the interface; independent codon tables and expansion tests
establish the specified biological assignments rather than these examples alone.

`alignment-example.fasta` and `alignment-request.json` are synthetic local
interface fixtures: DNA `ACGT`/`AGT` with full identifiers `reference`/`deletion`,
and a Myers value-only request. Their complete inputs, invocations and observed
outputs appear in the alignment CLI recipes.
They check adapters against existing Rust facades, not new reference fidelity.

`index-example.fastq` and `index-request.json` are synthetic local fixtures for
[persistent record indexing](../../../docs/examples/sequence-index-recipes.md).
The FASTQ preserves mixed sequence case, a repeated `+` header and a quality line
containing `@+` (Phred scores 31/10). The request retrieves `read2` from a sidecar
created by the preceding recipe; no prebuilt index or external executable is
needed. Separate generated, corrupt-input, source-change and real phiX174 tests
establish the index contract beyond this demonstration.

`annotated-example.gb` and `annotated-example.embl` are synthetic local
documentation fixtures for checked annotated output. The complete 15-symbol DNA
contains a nine-symbol CDS; the EMBL version also has an opaque database link
and an embedded-coordinate tRNA qualifier to demonstrate explicit refusal/loss.
They are not database accession records or scientific qualification fixtures.

The other byte-exact `annotated-*` public NCBI/ENA snapshots have separate
[provenance, expected biological facts and checksums](annotated-provenance.md).
Their tests run offline and retain documented unsupported cases. These data
intentionally retain original trailing spaces/newlines; `.gitattributes` prevents
line-ending conversion of the raw GenBank/EMBL snapshots.

## Real sequence

`phix174-e.fasta`: Escherichia phage phiX174 RefSeq **NC_001422.1**,
one-based inclusive coordinates **568..843**, plus strand, CDS `phiX174p08`
(endolysin, gene E), 276 nucleotides. Downloaded 2026-10-15 from NCBI:

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&seq_start=568&seq_stop=843&rettype=fasta&retmode=text

Expected 91-residue translation is the independently annotated protein
**NP_040709.1**, recorded in the public GenBank `/translation` qualifier:

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=gb&retmode=text

The record declares bacterial table **11**. This particular CDS starts with ATG
and uses only codons with identical assignments in tables 1 and 11, so its stated
translation also validates the table 1 subset. Configured-translation tests now
also select table 11 explicitly. This record alone does not exercise additional
start codons; separate exhaustive codon/ambiguity tests cover those initiation
sets. Tests read the local fixture offline; no test requires NCBI or an external
interpreter. Historical table-1 validation remains unchanged.

Other test reference data are the 64 assignments and initiation sets of NCBI
genetic codes 1, 2 and 11 (page revision 2024-09-23, consulted 2026-10-15)
and the INSDC/DDBJ IUPAC base/amino-acid ambiguity definitions:

- https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
- https://www.ddbj.nig.ac.jp/ddbj/code-e.html
- https://biopython.org/docs/1.88/api/Bio.Seq.html

The grouped-codon and string-set oracle does not call the production engine to
generate expected values. It tests specification-derived ambiguity rules, not a
live differential run against Biopython.

## NCBI genetic codes

`ncbi-gc.prt` is NCBI's genetic code table file, version 4.6, downloaded
unchanged on 2026-10-22 from
https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt (SHA-256
`2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04`). It defines
the 27 current tables: their amino acids (`ncbieaa`) and their initiators and
codons that can end translation (`sncbieaa`). The production tables were
transcribed from it; `ncbi_genetic_codes.rs` parses it with its own code and
compares every table, codon, initiator and terminator with the library. NCBI
data are a public United States government work.

## Streamed GenBank records

`annotated-stream.gb` is a synthetic three-record GenBank file written for the
[streaming reader guide](../../../docs/api/annotated-input-v1.md#stream-records-from-large-files):
two DNA records around one that declares RNA, so a DNA reader reports the middle
record as an error and can continue past it.

