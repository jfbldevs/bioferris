//! Demonstrates ancestral E-INS-i on protein, DNA, and RNA with internal regions.

use align_core::{Alphabet, Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use align_mafft::{EinsiTermination, MafftAlphabet, MafftEinsi};
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
    let result = MafftEinsi::<A>::new().align(records)?;
    print_alignment(label, result.alignment());
    for pair in result.pair_evidence() {
        println!(
            "pair {}-{}: ranges {:?}/{:?}, score={:.6}, distance={:.6}, segments={}, skip openings={}, paired skipped={}, unpaired skipped={}",
            records[pair.first_sequence()].identifier(),
            records[pair.second_sequence()].identifier(),
            pair.first_aligned_range(),
            pair.second_aligned_range(),
            pair.alignment_score(),
            pair.distance(),
            pair.segment_count(),
            pair.generalized_skip_opening_count(),
            pair.generalized_skipped_pair_count(),
            pair.generalized_unpaired_residue_count(),
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
        EinsiTermination::Converged | EinsiTermination::CycleLimit
    ));
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proteins = vec![
        SequenceRecord::<Protein>::new(
            "protein_a",
            ProteinSequence::from_ascii("ACDEFGHIKLMNPQRSTVWYWWWWWWWWWWACDEFGHIKLMNPQRSTVWY")?,
        )?,
        SequenceRecord::<Protein>::new(
            "protein_b",
            ProteinSequence::from_ascii(
                "ACDEFGHIKLMNPQRSTVWYDDDDDDDDDDDDDDDDDDACDEFGHIKLMNPQRSTVWY",
            )?,
        )?,
        SequenceRecord::<Protein>::new(
            "protein_c",
            ProteinSequence::from_ascii("ACDEFGHIKLMNPQRSTVWYCCCCCCCCCCCCCCACDEFGHIKLMNPQRSTVWY")?,
        )?,
    ];
    run("Protein E-INS-i", &proteins)?;

    let dna = vec![
        SequenceRecord::<Dna>::new(
            "dna_a",
            DnaSequence::from_ascii("ACGTACGTACGTAAAAAAAAACGTACGTACGT")?,
        )?,
        SequenceRecord::<Dna>::new(
            "dna_b",
            DnaSequence::from_ascii("ACGTACGTACGTCCCCCCCCCCCCCCCCCCACGTACGTACGT")?,
        )?,
        SequenceRecord::<Dna>::new(
            "dna_c",
            DnaSequence::from_ascii("ACGTACGTACGTGGGGGGGGGGGGGACGTACGTACGT")?,
        )?,
    ];
    run("DNA E-INS-i", &dna)?;

    let rna = vec![
        SequenceRecord::<Rna>::new(
            "rna_a",
            RnaSequence::from_ascii("ACGUACGUACGUAAAAAAAAACGUACGUACGU")?,
        )?,
        SequenceRecord::<Rna>::new(
            "rna_b",
            RnaSequence::from_ascii("ACGUACGUACGUCCCCCCCCCCCCCCCCCCACGUACGUACGU")?,
        )?,
        SequenceRecord::<Rna>::new(
            "rna_c",
            RnaSequence::from_ascii("ACGUACGUACGUGGGGGGGGGGGGGACGUACGUACGU")?,
        )?,
    ];
    run("RNA E-INS-i", &rna)?;
    Ok(())
}
