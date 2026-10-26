//! Aligns three proteins with the versioned default Clustal Omega profile.

use align_clustal::{ClustalOmega, OmegaGuideTreeMethod};
use align_core::{Protein, ProteinSequence};
use align_msa::SequenceRecord;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let records = [
        SequenceRecord::new(
            "human",
            ProteinSequence::from_ascii(b"MKTAYIAKQRQISFVKSHFSRQ")?,
        )?,
        SequenceRecord::new(
            "mouse",
            ProteinSequence::from_ascii(b"MKTAYIAKQKQISFVKSHFSRQ")?,
        )?,
        SequenceRecord::new(
            "fish",
            ProteinSequence::from_ascii(b"MKTAYIAKQRQISFVKSFSRQ")?,
        )?,
    ];

    let result = ClustalOmega::<Protein>::published().align(&records)?;
    println!("profile: {}", result.profile_identifier());
    println!("guide: {:?}", result.guide_evidence().method());
    assert_eq!(
        result.guide_evidence().method(),
        OmegaGuideTreeMethod::MbedEffectiveFull
    );
    for row in result.alignment().rows() {
        println!("{:<8} {}", row.identifier(), row.to_ascii_with_gaps());
    }
    for (index, step) in result.progressive_steps().iter().enumerate() {
        println!(
            "merge {index}: {:?}, score {:.3}, core {:?}/{:?}, posterior {:?}",
            step.engine(),
            step.score(),
            step.left_core(),
            step.right_core(),
            step.posterior_bounds()
        );
    }
    println!(
        "peak workspace estimate: {} bytes",
        result.memory_estimate().estimated_peak_bytes()
    );
    Ok(())
}
