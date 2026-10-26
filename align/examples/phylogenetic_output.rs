//! Writes one MSA in the supported downstream-analysis formats.

use std::path::PathBuf;

use align::prelude::*;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let Some(output_directory) = std::env::args_os().nth(1).map(PathBuf::from) else {
        eprintln!("usage: phylogenetic_output <output-directory>");
        return Ok(());
    };
    std::fs::create_dir_all(&output_directory)?;

    let records = [
        SequenceRecord::new("reference", DnaSequence::try_from("ACGTACGT")?)?,
        SequenceRecord::new("allele", DnaSequence::try_from("ACGTTCGT")?)?,
        SequenceRecord::new("deletion", DnaSequence::try_from("ACGTCGT")?)?,
    ];
    let run = Multiple::mafft(MafftStrategy::Linsi).run(&records)?;

    run.write_text(output_directory.join("alignment.txt"))?;
    run.write_json(output_directory.join("alignment.json"))?;
    run.write_fasta(output_directory.join("alignment.fasta"))?;
    run.write_clustal(output_directory.join("alignment.aln"))?;
    run.write_phylip(output_directory.join("alignment.phy"))?;
    run.write_nexus(output_directory.join("alignment.nex"))?;
    run.write_stockholm(output_directory.join("alignment.sto"))?;

    let strict_interleaved = PhylipOptions::new()
        .with_layout(PhylipLayout::Interleaved)
        .with_name_style(PhylipNameStyle::Strict)
        .with_block_width(60)?;
    run.write_phylip_with(
        output_directory.join("alignment-strict-interleaved.phy"),
        strict_interleaved,
    )?;

    println!("wrote alignment outputs to {}", output_directory.display());
    Ok(())
}
