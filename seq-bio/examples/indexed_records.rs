//! Build, persist, reopen and query a stable Sanger FASTQ source from Rust.

use std::io;
use std::path::PathBuf;

use bioferris::io::index::{IndexOptions, IndexedFormat, SequenceIndex};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    let [source, destination] = paths.as_slice() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "usage: cargo bio indexed_records -- SOURCE.fastq NEW_INDEX.json",
        )
        .into());
    };

    let options = IndexOptions::default();
    let index = SequenceIndex::<Dna>::build(source, IndexedFormat::FastqSanger, options)?;
    index.save_new(destination)?;
    if let Err(error) = index.save_new(destination) {
        println!("Existing destination refused: {}", error.code());
    }
    drop(index);

    let mut index = SequenceIndex::<Dna>::open(source, destination, options)?;
    println!("Records: {}", index.len());
    println!("IDs: {}", index.ids().collect::<Vec<_>>().join(", "));
    if let Some(record) = index.get("read2")? {
        println!("Header: {}", record.header());
        println!("Source sequence: {}", record.source_sequence());
        println!("Normalized sequence: {}", record.sequence());
        println!("Phred scores: {:?}", record.qualities());
    }
    println!("Missing ID: {}", index.get("absent")?.is_none());
    Ok(())
}
