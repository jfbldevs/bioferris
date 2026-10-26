//! Aligns short DNA, RNA, and protein examples with one common API.

use std::error::Error;

use align_core::{DnaSequence, LinearGap, MatchMismatch, ProteinSequence, RnaSequence};
use align_pairwise::{NeedlemanWunsch, verify_global_alignment};

fn main() -> Result<(), Box<dyn Error>> {
    let substitution = MatchMismatch::new(2, -1);
    let gap = LinearGap::new(-2)?;
    let aligner = NeedlemanWunsch::new(substitution, gap);

    let dna_first = DnaSequence::from_ascii(b"GATTACA")?;
    let dna_second = DnaSequence::from_ascii(b"GCATGCT")?;
    let dna_alignment = aligner.align(&dna_first, &dna_second)?;
    verify_global_alignment(&dna_first, &dna_second, &dna_alignment, &substitution, gap)?;
    let rendered = dna_alignment.render(&dna_first, &dna_second)?;
    println!("DNA score: {}", dna_alignment.score());
    println!("{}", rendered.first());
    println!("{}", rendered.second());

    let rna_first = RnaSequence::from_ascii(b"GCAUAC")?;
    let rna_second = RnaSequence::from_ascii(b"GCUUAC")?;
    let rna_alignment = aligner.align(&rna_first, &rna_second)?;
    verify_global_alignment(&rna_first, &rna_second, &rna_alignment, &substitution, gap)?;
    println!("RNA score: {}", rna_alignment.score());

    let protein_first = ProteinSequence::from_ascii(b"PLEASANTLY")?;
    let protein_second = ProteinSequence::from_ascii(b"MEANLY")?;
    let protein_alignment = aligner.align(&protein_first, &protein_second)?;
    verify_global_alignment(
        &protein_first,
        &protein_second,
        &protein_alignment,
        &substitution,
        gap,
    )?;
    println!("Protein score: {}", protein_alignment.score());

    Ok(())
}
