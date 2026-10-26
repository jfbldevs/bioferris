//! Runs the complete Rust program shown in the quick-start guide.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let coding = dna("ATGGCTTAA")?;
    assert_eq!(coding.reverse_complement()?.to_string(), "TTAAGCCAT");
    assert_eq!(coding.transcribe()?.to_string(), "AUGGCUUAA");

    let product = coding.translate_cds()?;
    assert_eq!(product.as_str(), "MA");
    let translated_protein = product.to_protein()?;
    assert_eq!(translated_protein, protein("MA")?);
    println!("{product}");

    let target = dna("GATTACA")?;
    let query = dna("GCATGCT")?;
    let run = Pairwise::gotoh()
        .scores(2, -1)
        .gaps(-5, -1)
        .run(&target, &query)?;
    println!("{}", run.report());
    Ok(())
}
