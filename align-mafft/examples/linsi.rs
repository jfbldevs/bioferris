//! Demonstrates ancestral L-INS-i on protein, DNA, and RNA with terminal flanks.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{LinsiTermination, MafftAlphabet, MafftLinsi};
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
    let result = MafftLinsi::<A>::new().align(records)?;
    print_alignment(label, result.alignment());
    for pair in result.pair_evidence() {
        println!(
            "pair {}-{}: ranges {:?}/{:?}, score={:.6}, distance={:.6}, segments={}",
            records[pair.first_sequence()].identifier(),
            records[pair.second_sequence()].identifier(),
            pair.first_aligned_range(),
            pair.second_aligned_range(),
            pair.alignment_score(),
            pair.distance(),
            pair.segment_count()
        );
    }
    println!(
        "library entries={}, attempted cycles={}, accepted cycles={}, termination={:?}, final WSP+I={:.6}, estimated peak={} bytes\n",
        result.library_statistics().entry_count(),
        result.cycles().len(),
        result.accepted_cycle_count(),
        result.termination(),
        result.final_score(),
        result.memory_estimate().estimated_peak_bytes()
    );
    assert!(matches!(
        result.termination(),
        LinsiTermination::Converged | LinsiTermination::CycleLimit
    ));
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proteins = vec![
        SequenceRecord::<Protein>::new(
            "protein_a",
            ProteinSequence::from_ascii("WWWWMKTAYIAKQRQISFVKSHFSRQWWWW")?,
        )?,
        SequenceRecord::<Protein>::new(
            "protein_b",
            ProteinSequence::from_ascii("DDMKTAYIAKQRTISFVKSHFSDD")?,
        )?,
        SequenceRecord::<Protein>::new(
            "protein_c",
            ProteinSequence::from_ascii("CCCCMKTAYIAQRQISFVKSYFSRQCCCC")?,
        )?,
    ];
    run("Protein L-INS-i", &proteins)?;

    let dna = vec![
        SequenceRecord::<Dna>::new("dna_a", DnaSequence::from_ascii("TTTTACGTACGTAAAA")?)?,
        SequenceRecord::<Dna>::new("dna_b", DnaSequence::from_ascii("GGACGTCGTCC")?)?,
        SequenceRecord::<Dna>::new("dna_c", DnaSequence::from_ascii("CCCCACGTACGGGG")?)?,
    ];
    run("DNA L-INS-i", &dna)?;

    let rna = vec![
        SequenceRecord::<Rna>::new("rna_a", RnaSequence::from_ascii("UUUUACGUACGUAAAA")?)?,
        SequenceRecord::<Rna>::new("rna_b", RnaSequence::from_ascii("GGACGUCGUCC")?)?,
        SequenceRecord::<Rna>::new("rna_c", RnaSequence::from_ascii("CCCCACGUACGGGG")?)?,
    ];
    run("RNA L-INS-i", &rna)?;
    Ok(())
}
