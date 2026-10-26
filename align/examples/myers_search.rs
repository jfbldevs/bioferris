//! Prints exact Myers approximate-search endpoints as a readable table.

use align::prelude::*;

fn main() -> align::Result<()> {
    let pattern = DnaSequence::try_from("GATTACA")?;
    let text = DnaSequence::try_from("TTTGATTACATTTGACTACA")?;
    let run = Pairwise::myers().run_find_all(&pattern, &text, 2)?;
    println!("{}", run.report());
    Ok(())
}
