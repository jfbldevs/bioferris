//! A first sequence-analysis program using the Rust library directly.

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
