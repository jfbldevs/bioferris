# SPOA 4.1.5 black-box fixture for POASTA

`global-affine-iterative-v1.json` is an inert artifact generated in the sibling
`rust-sequence-analysis-reference-validation/poa` laboratory. It records ten
DNA, RNA, and protein observations from the documented SPOA 4.1.5 command-line
interface at commit `167bdcc41a25e397cbd4487b6338989266d815d8`.

The external command selects global alignment (`-l 1`), match score `0`,
mismatch score `-4`, and equal first/second affine functions with open score
`-8` and extension score `-2`. Under the published POASTA convention this is
exactly mismatch cost `4`, opening cost `6`, and per-symbol extension cost `2`:
a positive gap of length `k` costs `6 + 2k`. Both FASTA MSA (`-r 1`) and GFA
graph (`-r 3`) outputs are retained as hashes and normalized observations.

The reference executable was built in Release mode with CMake 4.4.2 and g++
13.3.0. Its SHA-256 is
`672d8fe684c063a6154797a14e953aa8a41b88abe26a370f5c6f3bf090813a01`.
The fixture SHA-256 is
`9ac76f111ed7b294fbbf2e4f0c82e9d256044d59c16cd5ade218cd3e2b0b2687`.
The case specification and generator SHA-256 values are recorded inside the
fixture.

`align-poa/tests/poasta_reference.rs` consumes the fixture through `PoastaPoa` only.
No SPOA source, executable, native library, Python runtime, subprocess, or live
reference environment participates in Cargo builds or tests. The Rust test
deserializes only this reviewed JSON. Reference source was not inspected or used as
implementation guidance.
