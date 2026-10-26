//! Aligns DNA, RNA, and protein globally with one affine-gap API.

use std::error::Error;

use align_core::{AffineGap, DnaSequence, MatchMismatch, ProteinSequence, RnaSequence};
use align_pairwise::{Gotoh, verify_global_affine_alignment};

fn main() -> Result<(), Box<dyn Error>> {
    let substitution = MatchMismatch::new(2, -3);
    let gap = AffineGap::new(-5, -1)?;
    let aligner = Gotoh::new(substitution, gap);

    let dna_first = DnaSequence::from_ascii(b"AAACCCCBBB")?;
    let dna_second = DnaSequence::from_ascii(b"AAABBB")?;
    let dna_alignment = aligner.align(&dna_first, &dna_second)?;
    verify_global_affine_alignment(&dna_first, &dna_second, &dna_alignment, &substitution, gap)?;
    let rendered = dna_alignment.render(&dna_first, &dna_second)?;
    println!("DNA score: {}", dna_alignment.score());
    println!("{}", rendered.first());
    println!("{}", rendered.second());

    let rna_first = RnaSequence::from_ascii(b"AUGCU")?;
    let rna_second = RnaSequence::from_ascii(b"AUCU")?;
    let rna_alignment = aligner.align(&rna_first, &rna_second)?;
    verify_global_affine_alignment(&rna_first, &rna_second, &rna_alignment, &substitution, gap)?;
    println!("RNA score: {}", rna_alignment.score());

    let protein_first = ProteinSequence::from_ascii(b"MKTLL")?;
    let protein_second = ProteinSequence::from_ascii(b"MKLL")?;
    let protein_alignment = aligner.align(&protein_first, &protein_second)?;
    verify_global_affine_alignment(
        &protein_first,
        &protein_second,
        &protein_alignment,
        &substitution,
        gap,
    )?;
    println!("Protein score: {}", protein_alignment.score());

    Ok(())
}
