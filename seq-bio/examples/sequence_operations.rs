//! Checked slicing, literal search and case/quality-preserving record slicing.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("atatatn")?;
    println!("DNA: {sequence}");
    println!("Slice 1..5: {}", sequence.slice(1..5)?);
    println!("Concatenated: {}", sequence.concat(&dna("gg")?)?);
    println!("First TA: {:?}", sequence.find("ta")?);
    println!("Contains NN: {}", sequence.contains("NN")?);
    println!("ATA, non-overlapping: {}", sequence.count("ATA")?);
    println!("ATA, overlapping: {}", sequence.count_overlapping("ATA")?);
    let matches: Vec<_> = sequence
        .find_all("ATA", MatchOverlap::Overlapping)?
        .collect();
    println!("ATA positions: {matches:?}");
    println!("Empty-pattern count: {}", sequence.count("")?);
    println!(
        "Terminal empty slice: {} symbols",
        sequence.slice(7..7)?.len()
    );
    println!("RNA slice: {}", rna("augn")?.slice(1..3)?);
    println!("Protein BZ: {:?}", protein("MABZX")?.find("bz")?);

    let record = Record::<Dna>::from_ascii(
        "read1 mixed case",
        "aTaTaTn",
        Some(vec![30, 31, 32, 33, 34, 35, 36]),
    )?;
    let region = record.slice(1..5)?;
    println!("Record header: {}", region.header());
    println!("Record source: {}", region.source_sequence());
    println!("Record normalized: {}", region.sequence());
    println!("Record qualities: {:?}", region.qualities());
    println!("Original source: {}", record.source_sequence());

    if let Err(error) = sequence.slice(7..8) {
        println!("Invalid slice: {}", error.code());
    }
    if let Err(error) = sequence.find("AU") {
        println!("Invalid DNA pattern: {}", error.code());
    }
    Ok(())
}
