# BLAST result fixture provenance

## Contents

- [NCBI cookbook tabular rows](#ncbi-cookbook-tabular-rows)
- [Synthetic cases](#synthetic-cases)
- [XML evidence](#xml-evidence)
- [Retained public API observations](#retained-public-api-observations)
- [Checksums](#checksums)

## NCBI cookbook tabular rows

`blast-ncbi-cookbook.tsv` contains the first three result rows of NCBI's public
custom-separator example, retrieved 2026-10-17 from
[Display BLAST search results with custom output format](https://www.ncbi.nlm.nih.gov/books/NBK569862/).
The source example runs `blastn -db ecoli -outfmt '6 delim=@'`. Its `@` delimiters
were mechanically represented as tabs, with all 36 source field values unchanged.
These are published search-result data, not implementation source. This is a
three-row excerpt; no claim is made that it contains every HSP of that query.

Independent checks use the literal source scores and endpoints. In particular,
row 2 has decreasing subject coordinates 6928..6821, normalized to 6820..6928.
Its e-value source text `7.69e-30` remains available alongside its binary64 value.

## Synthetic cases

Unit tests create minimal explicitly labeled synthetic rows/comment blocks to
exercise custom column selection, omitted values, no-hit queries, grouping,
translated/reversed coordinates, malformed data, transport errors and resource
budgets. Their results are independent format expectations, not observed BLAST
executions or copied Biopython outputs.

## XML evidence

The following complete, unmodified public NCBI BLAST 2.9.0+ report data were
retrieved 2026-10-18 from Biopython's version-tagged public test-data directory.
Only biological result files were downloaded and inspected, not parser or test
implementation source. The filenames are renamed locally for discoverability;
their contents are unchanged.

| Local fixture | Original source |
|---|---|
| `blast-ncbi-2.9.0-nucleotide.xml` | [xml_2900_blastn_001.xml](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/Blast/xml_2900_blastn_001.xml) |
| `blast-ncbi-2.9.0-nucleotide-xml2.xml` | [xml_2900_blastn_001_v2.xml](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/Blast/xml_2900_blastn_001_v2.xml) |
| `blast-ncbi-2.9.0-translated.xml` | [xml_2900_blastx_001.xml](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/Blast/xml_2900_blastx_001.xml) |
| `blast-ncbi-2.9.0-translated-xml2.xml` | [xml_2900_blastx_001_v2.xml](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/Blast/xml_2900_blastx_001_v2.xml) |

Each report has one query and ten subject hits; some subjects have multiple HSPs.
Together they retain 48 HSPs. The nucleotide first hit has subject endpoints
101449177..101449144 and explicit reverse strand/frame information. The translated
first hit has query endpoints 20..343, frame +2 and 108 alignment columns: source
coordinates must remain nucleotide units rather than being divided by three.
XML2 descriptions retain separate accession/title/taxonomy fields. Synthetic
tests add repeated HitDescr entries, multiple reports, missing fields, no hits,
malformed XML, transport failures and explicit budgets.

`blast-example.xml`, `blast-example-xml2.xml`, `blast-example.tab7` and
`blast-minimal.tsv` are independently authored synthetic format demonstrations,
not observed search runs. Legacy/XML2 examples each contain two no-hit queries;
the commented tabular example contains one HSP and one no-hit query. The minimal
two-column input has no scores or coordinates to infer.
`blast-example.tab6` is the same one-row literal standard input used by the Rust
program. `blast-request.json` is the complete strict request selecting the two
columns in `blast-minimal.tsv`; it is request data, not a search result.

## Retained public API observations

`blast-biopython-1.88.json` is the unchanged captured output of permitted public
API calls on eight inputs: the four complete public reports, both synthetic XML
examples, the three-row cookbook excerpt and the commented tabular example.
It contains eleven queries and 52 HSPs. The reference environment was:

- Existing system CPython 3.12.3 on x86_64 Linux; no runtime installed/upgraded.
- Biopython 1.88 official CPython 3.12 manylinux wheel,
  [PyPI distribution](https://files.pythonhosted.org/packages/17/c2/613ca5f729b55b2f7a9c17a3cd59d3f053ee2816da2d19933acc6b8741ec/biopython-1.88-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl),
  SHA-256 `9384c30460832926e765319040d732d1252a713a05af1dd4cc8b19cb1d590439`.
- NumPy 2.3.3 official CPython 3.12 manylinux wheel,
  [PyPI distribution](https://files.pythonhosted.org/packages/51/64/7de3c91e821a2debf77c92962ea3fe6ac2bc45d0778c1cbe15d4fce2fd94/numpy-2.3.3-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl),
  SHA-256 `d9192da52b9745f7f0766531dcfa978b7763916f158bb63bdb8a1eca0068ab20`.

Wheel hashes were verified before extraction to a task-specific temporary
`site` directory. Nothing was globally installed. No reference implementation
source was opened or copied. The reference interpreter and wheels are not
required by the application or its automated Rust regression tests.

The capture invoked `Bio.Blast.NCBIXML.parse` for legacy/XML2 and
`Bio.SearchIO.parse(path, "blast-tab", comments=False|True)` for tabular input.
Documented public attributes supplied query/hit identities and lengths, ordered
HSPs, source endpoints, scores/counts, frames and aligned strings. Endpoint
normalization was the explicitly documented transformation
`[min(start,end)-1, max(start,end)]`. The environment was invoked with
`OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1` and the temporary directory
as `PYTHONPATH`; the capture used `json.dumps(..., indent=2, ensure_ascii=False,
allow_nan=False)`. The retained JSON header records exact package/interpreter
versions and each case's public API/profile. Source numeric strings, arbitrary
metadata and absent optional values are separately checked against the input,
not invented by the reference projection.

Rust tests compare every retained comparable explicit coordinate, score, count,
alignment, midline and object order exactly. The initial comparison identified
that NCBIXML fills an absent XML2 positive count from identity and synthesizes
some frames. The fixture was retained unchanged. Tests explicitly verify that
the source did not supply these values and that BioFerris leaves them absent;
source-preserving absence is a documented intentional difference, not a claim
that all Biopython object fields match. No claims of universal compatibility,
performance, reference formatting, XML lexical round-trip or live search fidelity
follow from these observations.

## Checksums

SHA-256 of the retained source/observation fixtures:

```text
589c872d76943cf8b7d43af297063a73fd2fd69e10ec17d2f901df495f814088  blast-ncbi-2.9.0-nucleotide.xml
3679c61cc2952643dea078e45e283afcb917466ad9dd3b6e7c9012373432a353  blast-ncbi-2.9.0-nucleotide-xml2.xml
59f1da333b3de576b70194af95367a2cfe05563c94ac8ec5c3c9d9f799c19087  blast-ncbi-2.9.0-translated.xml
c30798903c086013000d754e9f73907ecc529a1055cfd191c0b08469f2a314b9  blast-ncbi-2.9.0-translated-xml2.xml
3f3133e993d7b851c0e58e9b3a6aa283a3427bf596473e20a05bc733842b6893  blast-example.xml
f82978c2613e856ed53b244af3a5a3475c4d67d4af0782d6e6495b2543d6aa0e  blast-example-xml2.xml
90522bd152f43f72f32fe57cdae5aa8adde631ca2dee5b5f960a188278179794  blast-example.tab7
bb87b844d2e496d02d2efd8dc7540be50d52028a9019434cc754715c746d52cb  blast-example.tab6
a537c7045b1f0e4be7567d82c74917b5be1ad61bdbd941be6b3f88924c924a02  blast-request.json
a4f4dc70741258c7c1d96fa4391f026f4a061e44f859ad3b985b75790cbf339c  blast-minimal.tsv
7d3406924235ab1a5cc913b31eb340d8d3f742aafc997a3bf3f720acdee33f60  blast-ncbi-cookbook.tsv
a59a7d3692b42c6eed5f08d84c9bd97252cce59b2af2c4036ac70bfaa439d4eb  blast-biopython-1.88.json
```

The [input contract](../../../docs/algorithms/blast-input-v1.md) records current
integration state, API/CLI schema, safety boundaries and verified complete
input/output examples. Isolated draft checks do not substitute for parent
workspace integration and acceptance gates.
The [manual and agent guide](../../../docs/api/blast-input-v1.md) reuses that
contract and has a promoted Cargo example plus documentation-driven recipes.
