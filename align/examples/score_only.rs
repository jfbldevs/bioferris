//! Prints explicit objective summaries for every pairwise score-only path.

use align::prelude::*;

fn main() -> align::Result<()> {
    let first = DnaSequence::try_from("GATTACA")?;
    let second = DnaSequence::try_from("GCATGCT")?;

    println!(
        "{}\n",
        Pairwise::needleman_wunsch()
            .run_score(&first, &second)?
            .report()
    );
    println!(
        "{}\n",
        Pairwise::gotoh().run_score(&first, &second)?.report()
    );
    println!(
        "{}\n",
        Pairwise::myers().run_distance(&first, &second)?.report()
    );
    println!(
        "{}\n",
        Pairwise::wfa()
            .run_penalty(&first, &second)?
            .report()
            .verbosity(Verbosity::Detailed)
    );
    println!(
        "{}",
        Pairwise::xdrop()
            .run_score(&first, &second)?
            .report()
            .verbosity(Verbosity::Detailed)
    );
    Ok(())
}
