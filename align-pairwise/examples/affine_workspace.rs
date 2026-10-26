//! Reuses affine scratch storage across independent complete alignments.

use align_core::{AffineGap, DnaSequence, MatchMismatch, ProteinSequence};
use align_pairwise::{AffineLocalAligner, AffineWorkspace, Gotoh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let scoring = MatchMismatch::new(3, -4);
    let gap = AffineGap::new(-5, -1)?;
    let global = Gotoh::new(scoring, gap);
    let local = AffineLocalAligner::new(scoring, gap);
    let mut workspace = AffineWorkspace::new();

    // Optional pre-reservation moves scratch allocation outside a known hot loop.
    workspace.reserve(64, 64)?;

    let first_dna = DnaSequence::from_ascii(b"ACGTTGCA")?;
    let second_dna = DnaSequence::from_ascii(b"ACTTGCA")?;
    let global_alignment = global.align_with_workspace(&first_dna, &second_dna, &mut workspace)?;

    // The same workspace can be reused by another affine mode and alphabet. The
    // previously returned alignment remains independently owned.
    let first_protein = ProteinSequence::from_ascii(b"QQMTEYKLVVVG")?;
    let second_protein = ProteinSequence::from_ascii(b"GGMTEYKLVVVGPP")?;
    let local_alignment =
        local.align_with_workspace(&first_protein, &second_protein, &mut workspace)?;

    println!("global DNA score: {}", global_alignment.score());
    println!("local protein score: {}", local_alignment.score());
    println!(
        "retained scratch: {} score columns, {} traceback cells",
        workspace.score_row_capacity(),
        workspace.traceback_capacity()
    );
    Ok(())
}
