//! Runs fitting and local affine alignment through the shared typed API.

use align_core::{AffineGap, DnaSequence, MatchMismatch, ProteinSequence};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, FreeEndGaps, verify_affine_end_gap_alignment,
    verify_local_affine_alignment,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let substitution = MatchMismatch::new(3, -4);
    let gap = AffineGap::new(-5, -1)?;

    let query = DnaSequence::from_ascii(b"ACGT")?;
    let target = DnaSequence::from_ascii(b"TTACGTGG")?;
    let free_ends = FreeEndGaps::fit_first_in_second();
    let fitting = AffineEndGapAligner::new(substitution, gap, free_ends);
    let fitted = fitting.align(&query, &target)?;
    verify_affine_end_gap_alignment(
        &query,
        &target,
        &fitted,
        fitting.substitution(),
        fitting.gap(),
        fitting.free_ends(),
    )?;
    let rendered = fitted.render(&query, &target)?;
    println!("fitting score: {}", fitted.score());
    println!(
        "ranges: {:?} versus {:?}",
        fitted.first_range(),
        fitted.second_range()
    );
    println!("{}\n{}", rendered.first(), rendered.second());

    let first_protein = ProteinSequence::from_ascii(b"QQMTEYKAA")?;
    let second_protein = ProteinSequence::from_ascii(b"GGMTEYKPP")?;
    let local = AffineLocalAligner::new(substitution, gap);
    let conserved = local.align(&first_protein, &second_protein)?;
    verify_local_affine_alignment(
        &first_protein,
        &second_protein,
        &conserved,
        local.substitution(),
        local.gap(),
    )?;
    println!("local protein score: {}", conserved.score());
    println!(
        "ranges: {:?} versus {:?}",
        conserved.first_range(),
        conserved.second_range()
    );
    Ok(())
}
