//! Computes exact affine scores without constructing traceback state.

use align_core::{AffineGap, Dna, DnaSequence, MatchMismatch};
use align_pairwise::{AffineWorkspace, Gotoh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let first = DnaSequence::from_ascii(b"ACGTTGCA")?;
    let second = DnaSequence::from_ascii(b"ACTTGCA")?;
    let aligner = Gotoh::new(MatchMismatch::new(3, -4), AffineGap::new(-5, -1)?);

    let mut workspace = AffineWorkspace::new();
    workspace.reserve_score(first.len(), second.len())?;
    workspace.reserve_simd_score::<Dna>(first.len(), second.len())?;

    let score = aligner.score_with_workspace(&first, &second, &mut workspace)?;
    println!("exact affine score: {score}");
    println!("selected backend: {:?}", workspace.last_backend());
    println!("traceback capacity: {}", workspace.traceback_capacity());
    Ok(())
}
