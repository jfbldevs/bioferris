//! Executable input/output recipes for the existing alignment facade.
//!
//! Run with `cargo run -p bioferris --example alignment_recipes -- <recipe>`.
//! Each recipe computes its result; tests compare it with the documented stdout.

use std::error::Error;
use std::fmt::Write as _;
use std::io::{self, Write as _};

use bioferris::align::PairwiseMetric;
use bioferris::prelude::*;

type RecipeResult<T> = Result<T, Box<dyn Error>>;

const DNA_FASTA: &[u8] = b">reference\nacgt acgt\n>allele\nACGTTCGT\n>deletion\nACGTCGT\n";
const PROTEIN_FASTA: &[u8] = b">human\nMKTAYIAKQRQ\n>mouse\nMKTAYIAKQKQ\n>fish\nMKTAYIAQRQ\n";
const USAGE: &str =
    "usage: alignment_recipes <pairwise|edit-search|extension|multiple|graph|formats>";

fn main() -> RecipeResult<()> {
    let mut arguments = std::env::args().skip(1);
    let recipe = arguments
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, USAGE))?;
    if arguments.next().is_some() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE).into());
    }
    let output = recipe_output(&recipe)?;
    io::stdout().lock().write_all(output.as_bytes())?;
    Ok(())
}

fn recipe_output(recipe: &str) -> RecipeResult<String> {
    match recipe {
        "pairwise" => pairwise(),
        "edit-search" => edit_search(),
        "extension" => extension(),
        "multiple" => multiple(),
        "graph" => graph(),
        "formats" => formats(),
        _ => Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE).into()),
    }
}

fn append_pairwise<A: Alphabet, R: PairwiseAlignmentView>(
    output: &mut String,
    run: &PairwiseRun<'_, A, R>,
) -> RecipeResult<()> {
    let (objective, value) = match run.info().metric() {
        PairwiseMetric::Score(value) => ("score", value.to_string()),
        PairwiseMetric::Distance(value) => ("distance", value.to_string()),
        PairwiseMetric::Penalty(value) => ("penalty", value.to_string()),
    };
    writeln!(output, "{}: {objective}={value}", run.info().algorithm())?;
    writeln!(
        output,
        "target={:?} query={:?} cigar={}",
        run.alignment().first_range(),
        run.alignment().second_range(),
        run.cigar()
    )?;
    let rows = run.alignment().render(run.first(), run.second())?;
    writeln!(output, "target {}", rows.first())?;
    writeln!(output, "query  {}", rows.second())?;
    Ok(())
}

fn pairwise() -> RecipeResult<String> {
    let target = dna("GATTACA")?;
    let query = dna("GCATGCT")?;
    let mut output = String::new();
    append_pairwise(
        &mut output,
        &Pairwise::needleman_wunsch()
            .scores(2, -1)
            .gap(-2)
            .run(&target, &query)?,
    )?;
    append_pairwise(
        &mut output,
        &Pairwise::smith_waterman()
            .scores(2, -1)
            .gap(-2)
            .run(&target, &query)?,
    )?;
    append_pairwise(
        &mut output,
        &Pairwise::gotoh()
            .scores(2, -1)
            .gaps(-5, -1)
            .global()
            .run(&target, &query)?,
    )?;
    writeln!(output, "[affine gap example]")?;
    append_pairwise(
        &mut output,
        &Pairwise::gotoh()
            .scores(2, -1)
            .gaps(-5, -1)
            .global()
            .run(&dna("AACCCGG")?, &dna("AAGG")?)?,
    )?;
    Ok(output)
}

fn edit_search() -> RecipeResult<String> {
    let pattern = dna("GATTACA")?;
    let comparison = dna("GACTACA")?;
    let text = dna("TTTGATTACATTTGACTACA")?;
    let mut output = String::new();
    writeln!(
        output,
        "global_distance={}",
        Pairwise::myers().distance(&pattern, &comparison)?
    )?;
    append_pairwise(&mut output, &Pairwise::myers().run(&pattern, &comparison)?)?;
    writeln!(output, "search: maximum_distance=1")?;
    for hit in Pairwise::myers().find_all(&pattern, &text, 1)? {
        writeln!(output, "text_end={} distance={}", hit.end(), hit.distance())?;
    }
    append_pairwise(
        &mut output,
        &Pairwise::wfa()
            .penalties(4, 6, 2)
            .run(&pattern, &comparison)?,
    )?;
    append_pairwise(
        &mut output,
        &Pairwise::adaptive_affine()
            .penalties(4, 6, 2)
            .run(&pattern, &comparison)?,
    )?;
    Ok(output)
}

fn extension() -> RecipeResult<String> {
    let target = dna("TTGATTACAGG")?;
    let query = dna("CCGACTACATT")?;
    let run = Pairwise::xdrop()
        .scores(2, -3)
        .gaps(-6, -1)
        .threshold(8)
        .run_from_seed(&target, &query, XDropSeed::new(5, 5), XDropDirection::Both)?;
    let mut output = String::new();
    append_pairwise(&mut output, &run)?;
    Ok(output)
}

fn append_multiple<A, R>(output: &mut String, run: &MultipleRun<A, R>) -> RecipeResult<()>
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    writeln!(output, "{}", run.info().algorithm())?;
    writeln!(output, "profile={}", run.info().profile())?;
    writeln!(
        output,
        "rows={} columns={}",
        run.alignment().sequence_count(),
        run.alignment().width()
    )?;
    for row in run.alignment().rows() {
        writeln!(output, "{} {}", row.identifier(), row.to_ascii_with_gaps())?;
    }
    Ok(())
}

fn multiple() -> RecipeResult<String> {
    let records = parse_fasta::<Protein>(PROTEIN_FASTA)?;
    let mut output = String::new();
    append_multiple(
        &mut output,
        &Multiple::clustal_omega().threads(1).run(&records)?,
    )?;
    append_multiple(&mut output, &Multiple::muscle5().run(&records)?)?;
    append_multiple(
        &mut output,
        &Multiple::mafft(MafftStrategy::FftNs1).run(&records)?,
    )?;
    Ok(output)
}

fn graph() -> RecipeResult<String> {
    let records = parse_fasta::<Dna>(DNA_FASTA)?;
    let run = Graph::poasta().costs(4, 6, 2).run(&records)?;
    let mut output = String::new();
    append_multiple(&mut output, &run)?;
    writeln!(output, "[GFA]")?;
    let mut gfa = Vec::new();
    bioferris::io::write_poa_gfa_to(&mut gfa, run.result().graph())?;
    output.push_str(&String::from_utf8(gfa)?);
    Ok(output)
}

fn formats() -> RecipeResult<String> {
    let records = parse_fasta::<Dna>(DNA_FASTA)?;
    let mut fasta = Vec::new();
    bioferris::io::write_fasta_to(&mut fasta, &records)?;
    let reparsed = parse_fasta::<Dna>(fasta.as_slice())?;
    let mut output = String::new();
    writeln!(
        output,
        "records_equal_after_round_trip={}",
        records == reparsed
    )?;
    writeln!(output, "[normalized FASTA]")?;
    output.push_str(std::str::from_utf8(&fasta)?);

    let run = Multiple::mafft(MafftStrategy::FftNs1).run(&records)?;
    let mut exports = Vec::new();
    writeln!(exports, "[PHYLIP]")?;
    bioferris::io::write_phylip_to(&mut exports, run.alignment())?;
    writeln!(exports, "[NEXUS]")?;
    bioferris::io::write_nexus_to(&mut exports, run.alignment())?;
    writeln!(exports, "[Stockholm]")?;
    bioferris::io::write_stockholm_to(&mut exports, run.alignment())?;
    writeln!(exports, "[Myers distance JSON]")?;
    let pattern = dna("GATTACA")?;
    let comparison = dna("GACTACA")?;
    Pairwise::myers()
        .run_distance(&pattern, &comparison)?
        .write_json_to(&mut exports)?;
    output.push_str(std::str::from_utf8(&exports)?);
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GUIDE: &str = include_str!("../../docs/examples/alignment-recipes.md");

    fn assert_documented_output(recipe: &str) -> RecipeResult<()> {
        let marker = format!("<!-- output: {recipe} -->\n```text\n");
        let (_, remainder) = GUIDE
            .split_once(&marker)
            .ok_or_else(|| format!("missing output block for {recipe}"))?;
        let (expected, _) = remainder
            .split_once("\n```")
            .ok_or_else(|| format!("unclosed output block for {recipe}"))?;
        assert_eq!(recipe_output(recipe)?, format!("{expected}\n"), "{recipe}");
        Ok(())
    }

    #[test]
    fn pairwise_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("pairwise")
    }

    #[test]
    fn edit_search_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("edit-search")
    }

    #[test]
    fn seeded_extension_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("extension")
    }

    #[test]
    fn multiple_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("multiple")
    }

    #[test]
    fn graph_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("graph")
    }

    #[test]
    fn format_output_matches_documentation() -> RecipeResult<()> {
        assert_documented_output("formats")
    }
}
