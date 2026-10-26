//! Explicit NCBI genetic codes, forward frames and translation metadata.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATATGAAGA")?;
    println!("DNA: {sequence}");
    for genetic_code in [
        GeneticCode::Standard,
        GeneticCode::VertebrateMitochondrial,
        GeneticCode::BacterialArchaealPlastid,
    ] {
        let translation = sequence.translate_with(TranslationConfig {
            genetic_code,
            ..TranslationConfig::default()
        })?;
        println!(
            "Table {} ({}): {translation}",
            genetic_code.id(),
            genetic_code.name()
        );
    }

    let cds_config = TranslationConfig {
        genetic_code: GeneticCode::BacterialArchaealPlastid,
        mode: TranslationMode::Cds,
        ..TranslationConfig::default()
    };
    println!(
        "Table 11 CDS GTGAAATAA: {}",
        dna("GTGAAATAA")?.translate_with(cds_config)?
    );
    println!(
        "Table 11 RNA GUGAAAUAA: {}",
        rna("GUGAAAUAA")?.translate_with(cds_config)?
    );

    let framed = dna("CCATGAAAT")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::Two,
        partial_codons: PartialCodonPolicy::Ignore,
        ..TranslationConfig::default()
    })?;
    println!(
        "Frame {} CCATGAAAT: {framed}; skipped: {}; ignored: {}",
        framed.reading_frame().offset(),
        framed.skipped_leading_bases(),
        framed.ignored_trailing_bases()
    );
    println!("Stop-free protein: {}", framed.to_protein()?);

    let stopped = dna("AATGTGATAA")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::One,
        mode: TranslationMode::Cds,
        ..TranslationConfig::default()
    });
    match stopped {
        Err(TranslationError::InternalStop { offset }) => {
            println!("CDS error: internal stop at original nucleotide {offset}");
        }
        Err(error) => return Err(error.into()),
        Ok(_) => return Err("expected an internal-stop error".into()),
    }

    let empty = dna("AC")?.translate_with(TranslationConfig {
        reading_frame: ReadingFrame::Two,
        ..TranslationConfig::default()
    })?;
    println!(
        "Frame 2 AC: length {}; skipped: {}",
        empty.len(),
        empty.skipped_leading_bases()
    );
    Ok(())
}
