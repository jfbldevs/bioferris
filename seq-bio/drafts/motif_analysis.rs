//! IUPAC matching and explicit nucleotide count, probability and scoring workflows.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let source = dna("AGRYN")?;
    let motif = dna_motif("R")?;
    for ambiguity in [MotifInputAmbiguity::Definite, MotifInputAmbiguity::Possible] {
        let hits = motif.find_with(
            &source,
            MotifSearchOptions {
                ambiguity,
                ..Default::default()
            },
        )?;
        println!("{ambiguity:?}:");
        for hit in hits {
            println!("  {:?} {:?}", hit.strand(), hit.range());
        }
    }
    if let Err(error) = motif.find_with(
        &source,
        MotifSearchOptions {
            ambiguity: MotifInputAmbiguity::Reject,
            ..Default::default()
        },
    ) {
        println!("Reject: {}: {error}", error.code());
    }

    let instances = [dna("AC")?, dna("AC")?, dna("AC")?, dna("CA")?];
    let counts = MotifCounts::from_instances(&instances)?;
    let weights = counts.normalize([1.0; 4])?;
    let scores = weights.log_odds_uniform()?;
    println!("Counts A/C/G/T: {:?}", counts.columns());
    println!("PWM A/C/G/T: {:?}", weights.columns());
    println!("PSSM A/C/G/T: {:?}", scores.columns());
    let concrete = dna("ACGTAC")?;
    println!("Forward scores: {:?}", scores.calculate(&concrete, 100)?);
    for hit in scores.find_with(
        &concrete,
        MatrixSearchOptions {
            threshold: 2.0,
            ..Default::default()
        },
    )? {
        println!(
            "Scored hit: {:?} {:?} {}",
            hit.strand(),
            hit.range(),
            hit.score()
        );
    }

    let custom = MotifCounts::<Dna>::from_columns(vec![[1.0; 4]])?
        .normalize_default()?
        .log_odds([4.0, 2.0, 1.0, 1.0])?;
    println!("Normalized background: {:?}", custom.background());
    println!("Asymmetric PSSM: {:?}", custom.columns());
    for hit in custom.find_with(
        &dna("A")?,
        MatrixSearchOptions {
            threshold: 1.0,
            ..Default::default()
        },
    )? {
        println!(
            "Asymmetric hit: {:?} {:?} {}",
            hit.strand(),
            hit.range(),
            hit.score()
        );
    }

    let sparse = MotifCounts::from_instances(&[dna("A")?])?
        .normalize_default()?
        .log_odds_uniform()?;
    println!(
        "Zero-probability score for C: {}",
        sparse.score_window(&dna("C")?)?
    );
    let cancellation = MotifScores::<Dna>::from_columns(
        vec![
            [1.0, 0.0, 0.0, 0.0],
            [2_f64.powi(-25), 0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0, 0.0],
        ],
        [0.25; 4],
    )?;
    println!(
        "Float64 accumulation, final float32: {}",
        cancellation.score_window(&dna("AAA")?)?
    );
    Ok(())
}
