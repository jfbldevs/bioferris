//! Align protein, DNA, and RNA records with the common MUSCLE 5 API.

use align_core::{DnaSequence, ProteinSequence, RnaSequence};
use align_msa::SequenceRecord;
use align_muscle::{Muscle5, Muscle5Parameters};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proteins = vec![
        SequenceRecord::new(
            "human_beta_globin",
            ProteinSequence::from_ascii("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLST")?,
        )?,
        SequenceRecord::new(
            "human_alpha_globin",
            ProteinSequence::from_ascii("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS")?,
        )?,
        SequenceRecord::new(
            "whale_myoglobin",
            ProteinSequence::from_ascii("MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHL")?,
        )?,
    ];
    let protein_result = Muscle5::published_align().align(&proteins)?;
    println!("# protein / Align");
    for row in protein_result.alignment().rows() {
        println!(">{}\n{}", row.identifier(), row.to_ascii_with_gaps());
    }

    let dna = vec![
        SequenceRecord::new("dna_a", DnaSequence::from_ascii("TACGTACGT")?)?,
        SequenceRecord::new("dna_b", DnaSequence::from_ascii("ACGTCGT")?)?,
        SequenceRecord::new("dna_c", DnaSequence::from_ascii("ACGTACG")?)?,
    ];
    let dna_result = Muscle5::published_align().align(&dna)?;
    println!("# DNA / Align");
    for row in dna_result.alignment().rows() {
        println!(">{}\n{}", row.identifier(), row.to_ascii_with_gaps());
    }

    let rna = vec![
        SequenceRecord::new("rna_a", RnaSequence::from_ascii("UACGUACGU")?)?,
        SequenceRecord::new("rna_b", RnaSequence::from_ascii("ACGUCGU")?)?,
        SequenceRecord::new("rna_c", RnaSequence::from_ascii("ACGUACG")?)?,
    ];
    let super5 =
        Muscle5::with_parameters(Muscle5Parameters::published_super5().with_perturbation_seed(0))?;
    let rna_result = super5.align(&rna)?;
    println!("# RNA / Super5");
    for row in rna_result.alignment().rows() {
        println!(">{}\n{}", row.identifier(), row.to_ascii_with_gaps());
    }
    eprintln!(
        "protein_objective={} dna_objective={} rna_objective={} rna_clusters={} estimated_peak_bytes={}",
        protein_result.final_objective_score(),
        dna_result.final_objective_score(),
        rna_result.final_objective_score(),
        rna_result
            .super5_evidence()
            .map_or(0, |evidence| evidence.clusters().len()),
        rna_result.memory_estimate().estimated_peak_bytes(),
    );
    Ok(())
}
