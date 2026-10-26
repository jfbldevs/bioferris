//! Demonstrates ancestral G-INS-i on protein, DNA, and RNA input.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{GinsiTermination, MafftAlphabet, MafftGinsi};
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
    let result = MafftGinsi::<A>::new().align(records)?;
    print_alignment(label, result.alignment());
    println!(
        "pairs={}, library entries={}, attempted cycles={}, accepted cycles={}, termination={:?}, final WSP+I={:.6}, estimated peak={} bytes\n",
        result.pair_evidence().len(),
        result.library_statistics().entry_count(),
        result.cycles().len(),
        result.accepted_cycle_count(),
        result.termination(),
        result.final_score(),
        result.memory_estimate().estimated_peak_bytes()
    );
    assert!(matches!(
        result.termination(),
        GinsiTermination::Converged | GinsiTermination::CycleLimit
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
    run("Protein G-INS-i", &proteins)?;

    let dna = vec![
        SequenceRecord::<Dna>::new("dna_a", DnaSequence::from_ascii("TACGTACGT")?)?,
        SequenceRecord::<Dna>::new("dna_b", DnaSequence::from_ascii("ACGTCGT")?)?,
        SequenceRecord::<Dna>::new("dna_c", DnaSequence::from_ascii("ACGTACG")?)?,
        SequenceRecord::<Dna>::new("dna_d", DnaSequence::from_ascii("TACGTTACGT")?)?,
    ];
    run("DNA G-INS-i", &dna)?;

    let rna = vec![
        SequenceRecord::<Rna>::new("rna_a", RnaSequence::from_ascii("UACGUACGU")?)?,
        SequenceRecord::<Rna>::new("rna_b", RnaSequence::from_ascii("ACGUCGU")?)?,
        SequenceRecord::<Rna>::new("rna_c", RnaSequence::from_ascii("ACGUACG")?)?,
        SequenceRecord::<Rna>::new("rna_d", RnaSequence::from_ascii("UACGUUACGU")?)?,
    ];
    run("RNA G-INS-i", &rna)?;
    Ok(())
}
