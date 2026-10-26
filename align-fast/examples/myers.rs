//! Reusable Myers bit-vector distance, search, and global alignment.

use align_core::DnaSequence;
use align_fast::{MyersPattern, verify_unit_edit_alignment};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pattern = DnaSequence::from_ascii(b"ACGTAC")?;
    let text = DnaSequence::from_ascii(b"TTACGTCG")?;
    let myers = MyersPattern::new(&pattern)?;

    let distance = myers.global_distance(&text)?;
    let hits = myers.find_all(&text, 1)?;
    let alignment = myers.align_global(&text)?;
    let verified = verify_unit_edit_alignment(&pattern, &text, &alignment)?;
    let rendered = alignment.render(&pattern, &text)?;

    println!("global distance: {distance}");
    println!("search hits (exclusive end, distance):");
    for hit in hits {
        println!("  ({}, {})", hit.end(), hit.distance());
    }
    println!("{}", rendered.first());
    println!("{}", rendered.second());
    println!("independently verified edits: {}", verified.distance());

    Ok(())
}
