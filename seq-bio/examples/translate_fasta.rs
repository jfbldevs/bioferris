//! Translate complete DNA coding sequences from a FASTA file, offline.
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args_os()
        .nth(1)
        .ok_or("usage: translate_fasta <coding-sequences.fasta>")?;
    for record in read_fasta::<Dna>(path)? {
        let product = record.sequence().translate_standard(TranslationOptions {
            mode: TranslationMode::Cds,
            ..TranslationOptions::default()
        })?;
        // Validation is explicit; a full translation retaining '*' cannot be
        // silently passed to the historical stop-free protein alignment API.
        let protein = product.to_protein()?;
        println!(">{}\n{protein}", record.identifier());
    }
    Ok(())
}
