//! Handle invalid input: match a typed error, or request an explicit policy.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // An optional argument replaces the default input; `?` ends the program on
    // input that is not DNA.
    let input = std::env::args().nth(1).unwrap_or_else(|| "ATGA".to_owned());
    let sequence = dna(&input)?;

    match sequence.translate_cds() {
        Ok(protein) => println!("CDS {sequence}: {protein}"),
        Err(error @ TranslationError::PartialCodon { .. }) => {
            println!("Rejected {sequence}: {error}");
        }
        Err(error) => return Err(error.into()),
    }

    let translation = sequence.translate_standard(TranslationOptions {
        partial_codons: PartialCodonPolicy::Ignore,
        ..TranslationOptions::default()
    })?;
    println!(
        "Remainder ignored: {translation} (ignored bases: {})",
        translation.ignored_trailing_bases()
    );
    Ok(())
}
