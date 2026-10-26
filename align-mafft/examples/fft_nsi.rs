//! Demonstrates the scientific FFT-NS-i profile on protein, DNA, and RNA input.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{FftNsiTermination, MafftAlphabet, MafftFftNsi};
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
    let result = MafftFftNsi::<A>::new().align(records)?;
    print_alignment(label, result.alignment());
    println!(
        "attempted cycles={}, accepted cycles={}, termination={:?}, changed={}, final WSP={:.6}\n",
        result.cycles().len(),
        result.accepted_cycle_count(),
        result.termination(),
        result.changed_initial_alignment(),
        result.final_score()
    );
    assert!(matches!(
        result.termination(),
        FftNsiTermination::Converged
            | FftNsiTermination::Oscillating
            | FftNsiTermination::CycleLimit
    ));
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
    run("Protein FFT-NS-i", &proteins)?;

    let dna = vec![
        SequenceRecord::<Dna>::new("dna_a", DnaSequence::from_ascii("TACGTACGT")?)?,
        SequenceRecord::<Dna>::new("dna_b", DnaSequence::from_ascii("ACGTCGT")?)?,
        SequenceRecord::<Dna>::new("dna_c", DnaSequence::from_ascii("ACGTACG")?)?,
        SequenceRecord::<Dna>::new("dna_d", DnaSequence::from_ascii("TACGTTACGT")?)?,
    ];
    run("DNA FFT-NS-i", &dna)?;

    let rna = vec![
        SequenceRecord::<Rna>::new("rna_a", RnaSequence::from_ascii("UACGUACGU")?)?,
        SequenceRecord::<Rna>::new("rna_b", RnaSequence::from_ascii("ACGUCGU")?)?,
        SequenceRecord::<Rna>::new("rna_c", RnaSequence::from_ascii("ACGUACG")?)?,
        SequenceRecord::<Rna>::new("rna_d", RnaSequence::from_ascii("UACGUUACGU")?)?,
    ];
    run("RNA FFT-NS-i", &rna)?;
    Ok(())
}
