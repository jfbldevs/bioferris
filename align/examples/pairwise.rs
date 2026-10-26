//! Runs global affine DNA alignment through the ergonomic Rust facade.

use align::prelude::*;

fn main() -> align::Result<()> {
    let target = DnaSequence::try_from("GATTACA")?;
    let query = DnaSequence::try_from("GCATGCT")?;

    let aligner = Pairwise::gotoh().scores(2, -1).gaps(-5, -1);
    let run = aligner.run(&target, &query)?;

    println!("{}", run.report().labels("Target", "Query"));
    Ok(())
}
