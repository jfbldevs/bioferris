//! Exact global affine wavefront alignment.

use align_core::{DnaSequence, ProteinSequence};
use align_fast::{WavefrontAffine, WavefrontAffinePenalties};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let first = DnaSequence::from_ascii(b"GATACA")?;
    let second = DnaSequence::from_ascii(b"GAGATA")?;

    let penalties = WavefrontAffinePenalties::new(4, 6, 2)?;
    let affine = WavefrontAffine::new(penalties).align(&first, &second)?;
    let affine_rendered = affine.alignment().render(&first, &second)?;
    println!("affine penalty: {}", affine.penalty());
    println!("{}", affine_rendered.first());
    println!("{}", affine_rendered.second());
    println!("affine work: {:?}", affine.summary());

    let protein = ProteinSequence::from_ascii(b"MTEYKLVVVG")?;
    let protein_variant = ProteinSequence::from_ascii(b"MTEFKLVVVG")?;
    println!(
        "protein affine penalty: {}",
        WavefrontAffine::new(WavefrontAffinePenalties::new(5, 8, 2)?)
            .score(&protein, &protein_variant)?
            .penalty()
    );
    Ok(())
}
