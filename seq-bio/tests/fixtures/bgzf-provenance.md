# BGZF fixture provenance

## Contents

- [Public BAM byte stream](#public-bam-byte-stream)
- [Biopython-generated stream](#biopython-generated-stream)
- [Observed public API results](#observed-public-api-results)
- [Native sequence-index observations](#native-sequence-index-observations)
- [Scope of evidence](#scope-of-evidence)

## Public BAM byte stream

`bgzf-biopython-ex1.bam` is the unmodified public
[Biopython 1.88 SamBam/ex1.bam test fixture](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/SamBam/ex1.bam),
downloaded on 2026-10-17. Upstream path: `Tests/SamBam/ex1.bam`, release tag
`biopython-188`. It is retained as external test data under its upstream terms;
see the [Biopython licensing information](https://biopython.org/wiki/License).

- Compressed length: 124,584 bytes.
- SHA-256: `c2da01460656f9af6d58aef84c5b548f227db739c155ff6d5e4242f32236f017`.
- Public `Bio.bgzf` observations: 456,614 decompressed bytes, CRC32
  `1043922902`, prefix `42 41 4d 01` and eight blocks including the final marker.

The complete public block tuples are recorded in `bgzf-biopython-1.88.json` and
checked as literal expected values in the Rust tests. Importing this byte stream
does not add a BAM biological parser or establish algorithm execution provenance.

## Biopython-generated stream

`bgzf-biopython-generated.bgz` was generated through the permitted public
`Bio.bgzf.BgzfWriter` API in the existing isolated Biopython 1.88 environment.
No reference implementation source was read or copied. Compression level: 6.
Input is exactly this public data construction, with no biological interpretation:

```text
For zero-based index i in 0..140000: emit the byte (i * 17) modulo 256.
Then append these five bytes: 00 00 3e 0d 0a.
```

The writer received the first 65,536 input bytes and then the remaining bytes,
followed by explicit close/finalization. The 140,005 input bytes have CRC32
`3774126258`. The compressed fixture is 1,623 bytes with SHA-256
`ce044dea1f863f80d6e72bbb9270e623fb58d7e5346c384f552f359b01fd4390`.

## Observed public API results

`bgzf-biopython-1.88.json` contains actual executions of `BgzfBlocks`, binary
`BgzfReader.read/tell/seek`, `BgzfWriter.write/tell/close`,
`make_virtual_offset` and `split_virtual_offset`. It records complete block
tuples, selected boundary positions, component overflow errors and missing-EOF
behavior. SHA-256: `472588c6ffacf3b6d370a7910236c278d96406f055f012da4b3df0ced471d721`.

`bgzf-rust-to-biopython-1.88.json` records the reverse interoperability check:
the new Rust writer encoded the same 140,005-byte input at level 6 with a
65,536-byte maximum candidate size. A separate Biopython 1.88 binary reader
returned exactly the original input, confirmed by direct byte equality and CRC32.
Its SHA-256 is `775680a60c5e9fe7aa5a3651ca2122693b7ca5a4a31b2e9ddf6faca406b355a4`.

The Rust stream's block sizes were 603, 603, 376 and 28 bytes; native Biopython's
were 605, 605, 385 and 28 bytes. Both decoded to the same input. This is concrete
evidence that compressed-byte identity is not part of the current writer contract.

The temporary observation scripts used Python solely as an operator for permitted
public reference APIs, not as application code or runtime dependencies. New
automated implementation tests are Rust and use retained offline observations.

## Native sequence-index observations

On 2026-10-18 the existing CPython 3.12.3 executed public Biopython 1.88 APIs with
NumPy 2.3.3 from an isolated temporary package directory. No runtime upgrade,
global package installation or reference-source inspection was performed. Only
official pinned wheels were unpacked:

| Official package artifact | SHA-256 |
| --- | --- |
| [Biopython 1.88 cp312 manylinux wheel](https://files.pythonhosted.org/packages/17/c2/613ca5f729b55b2f7a9c17a3cd59d3f053ee2816da2d19933acc6b8741ec/biopython-1.88-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl) | `9384c30460832926e765319040d732d1252a713a05af1dd4cc8b19cb1d590439` |
| [NumPy 2.3.3 cp312 manylinux wheel](https://files.pythonhosted.org/packages/51/64/7de3c91e821a2debf77c92962ea3fe6ac2bc45d0778c1cbe15d4fce2fd94/numpy-2.3.3-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl) | `d9192da52b9745f7f0766531dcfa978b7763916f158bb63bdb8a1eca0068ab20` |

`Bio.bgzf.BgzfWriter(..., "wb", compresslevel=6)` received the first 65,536 bytes,
then the remaining bytes, followed by its context-manager close. For the small
inputs, the second write was empty. Generated binary fixtures:

| Fixture | SHA-256 |
| --- | --- |
| `bgzf-native-index.fasta.bgz` | `9d56cef4c81cd4793b3c862ad965f4a328fdb506c6efc9f3a908d7495be81160` |
| `bgzf-native-index.fastq.bgz` | `fb75c2075025134dd42341248aded088e46b666b2297f239d643bae612ad0fe6` |
| `bgzf-native-boundary.fasta.bgz` | `3de3be62684f48d613f3ed559709f2b8aef09b839254ff4c1c62949cce137f0a` |

Exact small FASTA input:

```text
>z description
aCgT
NN
>a empty

>m other
TT
```

Exact wrapped Sanger FASTQ input:

```text
@r description
aCg
TN
+r description
!I~
@+
@q
AA
+
~~
```

The third input is `b">a\n" + b"A" * 65532 + b"\n>b\nTT\n"`: its second header
begins at decompressed byte 65,536. All inputs end in one LF byte, with no CRs.
The observation JSON retains literal small inputs and all input SHA-256 values.

Public `Bio.bgzf.BgzfBlocks`, `Bio.SeqIO.index(path, "fasta"|"fastq")`, mapping
iteration/access and `.get_raw(id)` produced block tuples, ordered IDs, full
headers, sequence-case digests, Sanger qualities and raw-span lengths/digests.
Biopython's complete `SeqRecord.description` maps to BioFerris's complete
`Record.header`, not its separate description suffix. This documents a semantic
mapping, not API identity.

`bgzf-index-biopython-1.88.json` SHA-256:
`376c24835a9a9697eb750f41470f031210c97471af6883dfbc7872ed2e4c343c`.
An offline Rust test checks all three compressed hashes and block tables, all
seven records' headers, lengths, case-sensitive sequence SHA-256, qualities and
raw-span lengths/SHA-256, then saves/reopens each sidecar. Automated tests never
invoke Python or download packages.

`compression-example.fa` is independently authored synthetic documentation input,
not upstream scientific data. Its bytes are `>r description\naCgT\n`.

## Scope of evidence

The isolated draft core passed 28 Rust tests, including an independent bitwise
CRC oracle, directly constructed RFC 1951 stored blocks, 16,065 single-byte
mutations, all compression levels, block/offset boundaries, failed I/O and typed
error propagation. These tests and reference observations qualify the documented
format profile only. Broader gzip headers, compressed-byte identity, full native
API parity, create-only publication, CLI/JSON or universal sequence-index parity
are not established by the original byte fixtures. The new index observations
establish the seven retained records and documented mappings only. Separate
malformed/boundary tests cover BioFerris's stricter policies. Stage 6 now has
public Rust/CLI/strict-request wiring and promoted examples; final main-checkout
acceptance is tracked in the integration validation.
