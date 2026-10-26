# mafft-core

> The MAFFT alignment engine, in pure Rust. Byte-identical to C MAFFT 7.526.

> **Acknowledgment.** rust-MAFFT is a port of [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> by **Kazutaka Katoh** and colleagues (CBRC). The science is theirs.
> If you use this crate in published work, please cite
> [Katoh & Standley 2013](https://doi.org/10.1093/molbev/mst010);
> full guidance (including mode-specific references) at
> [the project citation page](https://luksgrin.github.io/rust-MAFFT/citation/).

`MafftEngine`, `AlignmentMode`, the progressive-alignment driver
(`merge_step_cached` / `Falign`), iterative refinement (segmented
NW/SW + oscillation detection), and the `--add` / `--addfragments`
machinery for grafting new sequences into an existing MSA. Output
matches the C reference byte-for-byte across BAliBASE 3 (1930/1930)
for every supported mode.

## Install

```sh
cargo add mafft-core
```

Most users want the top-level [`mafft`](https://crates.io/crates/mafft)
crate, which re-exports the public surface here and is more stable
across internal reorganisations. Depend on `mafft-core` directly when
you want a smaller dependency footprint or you're embedding the
engine without the I/O layer.

## Example

```rust,no_run
use mafft_core::{MafftEngine, AlignmentMode};
use mafft_types::{Sequence, SequenceSet, SeqType};

let input = SequenceSet {
    sequences: vec![
        Sequence { name: "a".into(), data: b"ACDEFGHIK".to_vec() },
        Sequence { name: "b".into(), data: b"ACDEFHIK".to_vec() },
    ],
    seq_type: SeqType::Protein,
};
let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
println!("aligned {} sequences to width {}", msa.sequences.len(), msa.width());
```

## Documentation

API reference: [docs.rs/mafft-core](https://docs.rs/mafft-core)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/) ·
[Byte-identity design notes](https://luksgrin.github.io/rust-MAFFT/architecture/byte-identity/)

## License

MIT for the Rust port; the alignment algorithms are derived from
upstream MAFFT under BSD-3-Clause. See `LICENSE-MIT` and `LICENSE-BSD`
in the workspace root.
