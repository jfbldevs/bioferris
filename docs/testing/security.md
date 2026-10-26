# Dependency and XML safety verification

## Contents

- [Finding and remediation](#finding-and-remediation)
- [Reproducible dependency check](#reproducible-dependency-check)
- [XML boundary regression](#xml-boundary-regression)
- [Interrupted-read regression](#interrupted-read-regression)
- [Policy and limitations](#policy-and-limitations)

[Testing strategy](README.md) · [Documentation index](../README.md)
· [BLAST contract](../algorithms/blast-input-v1.md) · [File-I/O guide](../api/file-io-v1.md)

## Finding and remediation

The initial current-advisory audit of the former `quick-xml` 0.38.3 pin returned
exit 1: RUSTSEC-2026-0194 and RUSTSEC-2026-0195. The former concerns quadratic
duplicate-attribute checks; the latter concerns unbounded namespace-declaration
allocation inside `NsReader`, which this project's XML adapter uses. These are
availability findings, not evidence of a demonstrated compromise of user data.
See the primary [attribute advisory](https://rustsec.org/advisories/RUSTSEC-2026-0194.html)
and [namespace advisory](https://rustsec.org/advisories/RUSTSEC-2026-0195.html).

The dependency is now exactly 0.41.0, the minimal patched release, with default
features disabled. Its [official manifest](https://github.com/tafia/quick-xml/blob/v0.41.0/Cargo.toml)
declares Rust 1.79; the project remains on 1.85.0. Only this registry package
changed in `Cargo.lock`. Public resolver APIs replace removed methods. The adapter
explicitly sets 256 namespace declarations per element. This equals the 0.41.0
upstream default; setting it keeps the limit fixed if a later release changes
that default. Former decoded/unescaped source-attribute whitespace is retained
rather than silently applying the new attribute normalization. No advisory is
newly ignored and no checks are disabled.

The same audit found that shared line parsing treated a transient
`std::io::ErrorKind::Interrupted` read as a hard I/O error. It is now retried; see
the [interrupted-read regression](#interrupted-read-regression).

## Reproducible dependency check

Input: the current root `Cargo.toml`, `Cargo.lock` and `deny.toml`. The policy
checks all features and reads advisory data only from
`target/qa-tools/advisory-db`, not a global cache. The observed Linux auditor is
the official prebuilt `cargo-deny` 0.20.2 release, unpacked under
`target/qa-tools`; it does not change the product compiler. Its binary SHA-256 is
`b329e25933d01c36dd7c47d84ea5716694f9b7caf53a5003d45674703a8ed54a`. The
CI setup action contains its exact
download invocation and archive SHA-256 verification. The observed archive hash:

```text
9f12ed4c49936e09b48bf862b595cde2fe64fcbd9d74dfacac6131ca824c8d5f
```

With that auditor and advisory snapshot present, run from the repository root:

```bash
./target/qa-tools/cargo-deny-0.20.2-x86_64-unknown-linux-musl/cargo-deny --locked --offline check
```

Observed complete stdout, exit 0; stderr empty:

<!-- qa-security-output -->
```text
advisories ok, bans ok, licenses ok, sources ok
```

Observed advisory database commit:
`999ec89d5d3f10f20be14750ec8c2d8b44e0ac9b`, dated 2026-10-17. This is a dated
result, not an assertion that future advisory checks must pass. A missing auditor
or advisory directory makes the command fail; nothing is downloaded or installed
by it in offline mode.

For durable execution, set `BIOFERRIS_QA_DENY_BIN` to that executable's absolute
path, then run the `dependencies` QA profile. It records the actual tool/scan
statuses, binary hash, database commit and age. A separate successful network
fetch is required in the current run, followed by an offline policy scan. Its
completion receipt must be no older than one hour; an old publisher commit date
alone does not make a successfully refreshed snapshot stale. Missing tools/data,
failed/stale fetch receipts, version mismatch,
new advisories or failed policy checks fail closed. No automatic installation,
package upgrade or exception change is performed by the runner.

## XML boundary regression

The [focused Rust regression](../../seq-bio/tests/xml_security_regressions.rs)
starts from the committed legacy and XML2 fixtures
[`blast-example.xml`](../../seq-bio/tests/fixtures/blast-example.xml) and
[`blast-example-xml2.xml`](../../seq-bio/tests/fixtures/blast-example-xml2.xml).
It inserts distinct declarations into the root element, for example:

```xml
<BlastOutput xmlns:n0="urn:n0" xmlns:n1="urn:n1" ... xmlns:n256="urn:n256">
```

Exactly 256 declarations per element are accepted with the same parsed queries
as the unmodified fixture. XML2's default namespace counts toward that total.
The 257th declaration is rejected. Observed error for the legacy fixture, printed
as `kind: {:?}` of `BioIoError::kind()` followed by the error's `Display` text:

<!-- security-namespace-limit-output -->
```text
kind: InvalidStructure
BLAST XML: start tag declares more than 256 namespace bindings; raise the limit with NamespaceResolver::set_max_declarations_per_element (record 1, line 2)
```

The error kind is the stable contract. The message text is produced by the
dependency and names an upstream method; BioFerris does not expose a way to raise
this limit. Duplicate declarations of one prefix are also rejected as
`InvalidStructure`.

Programs branch on `BioIoError::code()`, which returns the stable code
`invalid_structure`; `record()` and `line()` return `Some(1)` and `Some(2)`.
No partial result is returned. Keep at most 256 declarations per element; the
limit cannot be raised through BioFerris. The
`documented_namespace_limit_error_matches_the_library` test generates the same
257-declaration input and checks the complete report above against the library.

Run the focused regression tests from the repository root:

```bash
cargo test -p bioferris --test xml_security_regressions
```

Observed test-binary stdout, exit 0 (Cargo's build progress goes to stderr):

```text
running 4 tests
test blast_xml_duplicate_root_attributes_are_rejected ... ok
test interrupted_line_reads_retry_without_losing_a_partial_line ... ok
test documented_namespace_limit_error_matches_the_library ... ok
test blast_xml_namespace_declarations_are_capped_at_256_per_element ... ok

test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

Test order and timing can vary between runs. Against the previous source with
`quick-xml` 0.38.3, the cap test failed because 257 declarations were accepted,
and the interrupted-read test failed with `Err(Io)` at byte 0. The duplicate check
passed on both versions: it guards rejection, not the quadratic cost fixed upstream.

The broader [Rust robustness suite](../../seq-bio/tests/robust_streaming.rs) also
checks generated legacy/XML2 inputs at the 256/257 boundary under one-byte reads,
and 512 unsupported ordinary attributes. A separately owned test subprocess has a
ten-second watchdog; no hostile unbounded XML is run. Literal and escaped
tabs/newlines in retained schema-location values remain unchanged. See
[robustness examples and observed results](robustness.md).

## Interrupted-read regression

The same test file wraps finite inputs in a reader that returns one byte per call
and fails once with `Interrupted` at a chosen offset, behind a four-byte
`BufReader`. Every offset is tried for this FASTA input, shown with Rust byte-string escapes:

```text
>sample mixed case\r\nAc\r\ngTN\r\n>second\nRY\n
```

The equivalent FASTQ input and the committed
[`blast-example.tab6`](../../seq-bio/tests/fixtures/blast-example.tab6) fixture
are also checked at every offset. Each result must equal parsing without
interruption, for example two FASTA records, `sample` (`ACGTN`) and
`second` (`RY`). Other read errors are unchanged: they stay typed `Io` errors and
fuse iteration.

## Policy and limitations

The existing `paste` informational unmaintained advisory exception remains
explicit, with owner/review date in `deny.toml`; it is not an exception for these
vulnerabilities. Scanner network access is development tooling only: biological
XML parsing never resolves external DTDs/entities or performs network requests.

These finite checks do not establish complete XML compatibility or a universal
CPU/RAM bound. A permanently blocking or endlessly interrupting custom reader
still requires a caller-managed deadline. Only the Linux audit was run locally; other OS/architecture CI lanes are
configured but not executed.
The existing native/public XML snapshots and all historical evidence are unchanged.
