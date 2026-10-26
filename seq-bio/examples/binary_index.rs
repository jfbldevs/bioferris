//! Save JSON and binary sidecars, then reopen the binary one with full and quick
//! source checks.

use std::io;
use std::path::PathBuf;

use bioferris::io::index::{
    IndexOptions, IndexedFormat, SequenceIndex, SidecarFormat, SourceCheck,
};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [source, json, binary] = paths.as_slice() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "usage: cargo bio binary_index -- SOURCE.fastq NEW_INDEX.json NEW_INDEX.bfidx",
        )
        .into());
    };

    let options = IndexOptions::default();
    let index = SequenceIndex::<Dna>::build(source, IndexedFormat::FastqSanger, options)?;
    index.save_new(json)?;
    index.save_new_as(binary, SidecarFormat::Binary)?;
    drop(index);
    println!("JSON sidecar: {} bytes", std::fs::metadata(json)?.len());
    println!("Binary sidecar: {} bytes", std::fs::metadata(binary)?.len());

    for check in [SourceCheck::Full, SourceCheck::Quick] {
        let mut index = SequenceIndex::<Dna>::open_with(source, binary, options, check)?;
        let record = index.get("read2")?.ok_or("missing read2")?;
        println!(
            "{check:?} check: {} records; read2 {} {:?}",
            index.len(),
            record.source_sequence(),
            record.qualities()
        );
    }
    if let Err(error) = SequenceIndex::<Dna>::open_with(source, json, options, SourceCheck::Quick) {
        println!("Quick check of the JSON sidecar: {}: {error}", error.code());
    }
    Ok(())
}
