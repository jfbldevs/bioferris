//! Minimal two-way affine X-drop seed-extension example.

use align_core::{AffineGap, DnaSequence, MatchMismatch};
use align_fast::{XDrop, XDropDirection, XDropSeed};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let first = DnaSequence::from_ascii(b"TTACGTAA")?;
    let second = DnaSequence::from_ascii(b"GGACGTCC")?;
    let aligner = XDrop::new(MatchMismatch::new(2, -3), AffineGap::new(-6, -1)?, 2)?;
    let seed = XDropSeed::new(3, 3);

    let result = aligner.align_from_seed(&first, &second, seed, XDropDirection::Both)?;
    let rendered = result.alignment().render(&first, &second)?;
    println!("score: {}", result.score());
    println!(
        "ranges: {:?} versus {:?}",
        result.first_range(),
        result.second_range()
    );
    println!("{}", rendered.first());
    println!("{}", rendered.second());
    println!("metadata: {:?}", result.summary());
    Ok(())
}
