//! Translate all six frames and discover strand-aware, stop-terminated ORFs.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("ATGAAATAA")?;
    println!("DNA: {sequence}");
    for frame in sequence.translate_six_frames() {
        println!(
            "{:?} frame {} {:?}: {} (leading {}, trailing {})",
            frame.strand(),
            frame.reading_frame().offset(),
            frame.range(),
            frame.peptide(),
            frame.skipped_leading_bases(),
            frame.ignored_trailing_bases()
        );
    }
    for input in ["ATGAAATAA", "TTATTTCAT"] {
        for orf in dna(input)?.find_orfs()? {
            println!(
                "ORF {input}: {:?} frame {} {:?}, stop {:?}, protein {}",
                orf.strand(),
                orf.reading_frame().offset(),
                orf.range(),
                orf.stop_range(),
                orf.to_protein()?
            );
        }
    }
    match dna("ATGNNNTAA")?.find_orfs() {
        Err(error) => println!("Ambiguous input: {error}"),
        Ok(_) => return Err("expected ambiguous ORF input to be rejected".into()),
    }

    let bacterial = OrfOptions {
        genetic_code: GeneticCode::BacterialArchaealPlastid,
        strands: StrandSelection::Forward,
        start_policy: OrfStartPolicy::GeneticCode,
        ..OrfOptions::default()
    };
    for orf in rna("GUGAAAUAA")?.find_orfs_with(bacterial)? {
        println!(
            "GUG start, table 11: {:?}, stop {:?}, protein {}",
            orf.range(),
            orf.stop_range(),
            orf.to_protein()?
        );
    }
    let partial = OrfOptions {
        strands: StrandSelection::Forward,
        stop_policy: OrfStopPolicy::AllowPartial,
        ..OrfOptions::default()
    };
    for orf in dna("ATGAA")?.find_orfs_with(partial)? {
        println!(
            "Unterminated ATGAA: {:?}, stop {:?}, protein {}",
            orf.range(),
            orf.stop_range(),
            orf.to_protein()?
        );
    }
    let no_results = OrfOptions {
        max_orfs: 0,
        ..OrfOptions::default()
    };
    if let Err(error) = dna("ATGAAATAA")?.find_orfs_with(no_results) {
        println!("Zero result budget: {error}");
    }
    Ok(())
}
