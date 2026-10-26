//! Prints compact summaries for all seven pairwise algorithm families.

use align::prelude::*;

fn main() -> align::Result<()> {
    let first = DnaSequence::try_from("GATTACA")?;
    let second = DnaSequence::try_from("GCATGCT")?;

    println!(
        "{}\n",
        Pairwise::needleman_wunsch()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Pairwise::smith_waterman()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Pairwise::gotoh()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Pairwise::myers()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Pairwise::wfa()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Pairwise::adaptive_affine()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}",
        Pairwise::xdrop()
            .run(&first, &second)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    Ok(())
}
