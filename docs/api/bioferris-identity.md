# BioFerris identity and Rust naming

## Contents

- [Public identity](#public-identity)
- [Run an analysis](#run-an-analysis)
- [Migration and compatibility](#migration-and-compatibility)
- [Rust naming rules](#rust-naming-rules)
- [Verification](#verification)

[API capability index](../README.md#api-index) · [Quick start](../quick-start.md)

## Public identity

- Project: **BioFerris**.
- Cargo package and Rust import: `bioferris`. The `bioferris` executable was
  removed on 2026-10-21 (decision 029).
- Tagline: **Rust-native bioinformatics, from sequences to workflows.**

BioFerris is the toolkit name, not a new biological algorithm or a claim of full
Biopython compatibility. The package is local and unpublished (`publish = false`);
adopting the name does not reserve a registry name, domain or trademark.

## Run an analysis

Run from the root of this repository.
The literal input is coding DNA `ATGGCTTAA`. This complete program is checked in
as [`my_analysis.rs`](../../seq-bio/examples/my_analysis.rs):

<!-- identity-program -->
```rust
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATGGCTTAA")?;

    println!("DNA: {sequence}");
    println!("Complement: {}", sequence.complement()?);
    println!("Reverse complement: {}", sequence.reverse_complement()?);
    println!("RNA: {}", sequence.transcribe()?);
    println!("Protein: {}", sequence.translate_cds()?);

    Ok(())
}
```

Invocation:

```bash
cargo bio my_analysis
```

Complete stdout, exit 0:

<!-- identity-output -->
```text
DNA: ATGGCTTAA
Complement: TACCGAATT
Reverse complement: TTAAGCCAT
RNA: AUGGCUUAA
Protein: MA
```

The CDS translation uses the existing standard code and omits its terminal stop.
Invalid input still returns the same typed biological errors; renaming does not
relax alphabet or translation validation. See the [sequence contract](seq-bio-v1.md).

AI agents use the same interface as people: a Rust program that imports the
library. The command-line agent contract is kept only as a
record of the removed CLI.

## Migration and compatibility

| Surface | Previous spelling | Current spelling |
|---|---|---|
| Cargo package | `seq-bio` | `bioferris` |
| Rust import | `use seq_bio::prelude::*;` | `use bioferris::prelude::*;` |
| CLI binary | `seq-bio` | `bioferris`, itself removed on 2026-10-21 |
| Package selection | `cargo build -p seq-bio` | `cargo build -p bioferris` |
| HTML API entry point | `target/doc/seq_bio/index.html` | `target/doc/bioferris/index.html` |
| Example shortcut | `cargo bio my_analysis` | Unchanged |

The source directory is still `seq-bio/`, and the active checkout path inside the
workspace is unchanged. On 2026-10-20 the enclosing local workspace folder, which
is not a Git repository, was renamed from `RustBio` to `BioFerris`;
repository-relative paths are unaffected. Existing source links, fixture paths
and versioned guide filenames remain valid. A path dependency now uses the key
`bioferris` with the same path to `seq-bio/`. Historical reports, copyright notices,
raw data, schema identifiers and archived commands are not rewritten.

This is a public naming migration, not a compatibility alias: old package/import
names and the old CLI target are no longer built. Update client manifests,
imports and launch commands. Previously built `seq-bio` executables may still
exist in a local `target/` directory; they are stale artifacts, not aliases.
The rename does not delete those files or publish/install any executable globally.

The rename left biological method names, results, defaults, exit statuses and
JSON schemas unchanged; the CLI's help advertised `bioferris`. The identifier
`seq-bio-original-biological-envelope-v1` still names the result-import profile
for responses of the removed CLI. Archived verification must
be reproduced at its recorded commit, not presented as evidence for renamed code.

## Rust naming rules

| Item | Convention | Example |
|---|---|---|
| Variables and fields | `snake_case` | `input_sequence`, `record_count` |
| Functions and methods | `snake_case` | `reverse_complement`, `translate_cds` |
| Modules and new Rust files | `snake_case` | `result_json`, `file_commands.rs` |
| Structs, enums, traits, aliases, variants | `UpperCamelCase` | `DnaSequence`, `TranslationMode`, `Cds` |
| Constants and statics | `SCREAMING_SNAKE_CASE` | `MAX_RECORD_BYTES` |
| Type parameters and lifetimes | Concise uppercase types; lowercase lifetimes | `A`, `T`, `'a`, `'src` |

Use acronym casing such as `DnaSequence`, not `DNASequence`. Rust models data and
behavior through structs/enums, traits and `impl` blocks rather than classes.
Prefer meaningful names, conventional getters such as `sequence()`, and `as_`,
`to_` or `into_` conversion names that reflect borrowing, allocation and ownership.
These conventions follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/naming.html).

The library denies Rust's `nonstandard_style` lint group. Strict
Clippy catches naming warnings in tests and examples as well. These checks enforce
casing, not every semantic naming choice; code review still checks clarity.
Serialized field names and external format tokens obey their own contracts and
must not be changed solely to resemble Rust identifiers.

## Verification

The branding test checks the Cargo package name and description; its checks of
the removed executable's help and archived responses ran until the CLI removal.
Documentation-driven tests compare the program and its complete stdout with the
existing verified example. See the
rename validation report for actual checks.
