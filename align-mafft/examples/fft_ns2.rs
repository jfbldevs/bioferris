//! Demonstrates the scientific FFT-NS-2 profile on protein, DNA, and RNA input.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{MafftAlphabet, MafftFftNs2};
use align_msa::{MultipleAlignment, SequenceRecord};

fn print_alignment<A: Alphabet>(label: &str, alignment: &MultipleAlignment<A>) {
    println!("{label}");
    for row in alignment.rows() {
        println!(">{}\n{}", row.identifier(), row.to_ascii_with_gaps());
    }
}

fn run<A: MafftAlphabet>(
    label: &str,
    records: &[SequenceRecord<A>],
) -> Result<(), Box<dyn std::error::Error>> {
    let result = MafftFftNs2::<A>::new().align(records)?;
    print_alignment(label, result.alignment());
    println!(
        "first merges={}, second merges={}, second tree changed={}, root score={:?}\n",
        result.first_pass().merges().len(),
        result.second_merges().len(),
        result.guide_tree_changed(),
        result.second_root_merge_score()
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proteins = vec![
        SequenceRecord::<Protein>::new(
            "human_beta",
            ProteinSequence::from_ascii("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLST")?,
        )?,
        SequenceRecord::<Protein>::new(
            "human_alpha",
            ProteinSequence::from_ascii("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS")?,
        )?,
        SequenceRecord::<Protein>::new(
            "whale_myoglobin",
            ProteinSequence::from_ascii("MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHL")?,
        )?,
    ];
    run("Protein FFT-NS-2", &proteins)?;

    let dna = vec![
        SequenceRecord::<Dna>::new("dna_a", DnaSequence::from_ascii("TACGTACGT")?)?,
        SequenceRecord::<Dna>::new("dna_b", DnaSequence::from_ascii("ACGTCGT")?)?,
        SequenceRecord::<Dna>::new("dna_c", DnaSequence::from_ascii("ACGTACG")?)?,
    ];
    run("DNA FFT-NS-2", &dna)?;

    let rna = vec![
        SequenceRecord::<Rna>::new("rna_a", RnaSequence::from_ascii("UACGUACGU")?)?,
        SequenceRecord::<Rna>::new("rna_b", RnaSequence::from_ascii("ACGUCGU")?)?,
        SequenceRecord::<Rna>::new("rna_c", RnaSequence::from_ascii("ACGUACG")?)?,
    ];
    run("RNA FFT-NS-2", &rna)?;
    Ok(())
}
