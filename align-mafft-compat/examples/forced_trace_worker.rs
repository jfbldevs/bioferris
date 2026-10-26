//! Original-record read-only diagnostic view of the existing forced-model engine.
//! No native buffers, alignments, guide trees or scores are accepted as inputs.
use align_core::{Alphabet, Dna, Protein, Rna, Sequence as TypedSequence, Symbol};
use mafft_7526_engine::{AlignmentMode, MafftEngine, SeqType, Sequence, SequenceSet};
use serde_json::{Value, json};
use std::error::Error;
use std::io::{self, Read};

fn sequences<A: Alphabet>(request: &Value) -> Result<Vec<Sequence>, Box<dyn Error>> {
    request["records"]
        .as_array()
        .ok_or("original records required")?
        .iter()
        .map(|r| {
            let source =
                TypedSequence::<A>::from_ascii(r["sequence"].as_str().ok_or("sequence required")?)?;
            Ok(Sequence {
                name: r["id"].as_str().ok_or("identifier required")?.to_owned(),
                data: source.iter().map(Symbol::as_ascii).collect(),
            })
        })
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let request: Value = serde_json::from_str(&input)?;
    let (mut sequences, seq_type) = match request["alphabet"].as_str() {
        Some("dna") => (sequences::<Dna>(&request)?, SeqType::Dna),
        Some("rna") => (sequences::<Rna>(&request)?, SeqType::Rna),
        Some("protein") => (sequences::<Protein>(&request)?, SeqType::Protein),
        _ => return Err("explicit molecule required".into()),
    };
    let normalize = request
        .get("normalize_nucleotide_case")
        .map_or(Ok(false), |value| {
            value
                .as_bool()
                .ok_or("boolean nucleotide normalization required")
        })?;
    if normalize {
        if !seq_type.is_nucleotide() {
            return Err("nucleotide normalization needs a nucleotide model".into());
        }
        for sequence in &mut sequences {
            sequence.data.make_ascii_lowercase();
        }
    }
    if sequences.is_empty() || sequences.iter().any(|r| r.data.is_empty()) {
        return Err("nonempty original source required".into());
    }
    let iterations = usize::try_from(request.get("iterations").map_or(Ok(1000), |v| {
        v.as_u64().ok_or("unsigned iterations required")
    })?)?;
    if iterations > 1000 {
        return Err("diagnostic iterations outside controlled range".into());
    }
    let mode = match request["profile"].as_str() {
        Some("mafft-linsi") => AlignmentMode::LInsi { iterations },
        Some("mafft-einsi") => AlignmentMode::EInsi { iterations },
        _ => return Err("explicit supported profile required".into()),
    };
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    let result = pool.install(|| {
        MafftEngine::new(mode)
            .with_reorder(false)
            .align(&SequenceSet {
                sequences,
                seq_type,
            })
    });
    let distances = result.distance_matrix.as_ref().map(|dm| {
        (0..dm.nseq)
            .map(|i| {
                (0..dm.nseq)
                    .map(|j| dm.get(i, j).to_bits())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    });
    let guide = result.guide_tree.as_ref().map(|tree| {
        tree.steps
            .iter()
            .map(|s| {
                json!({"left":s.left,"right":s.right,"left_length_bits":s.left_length.to_bits(),
               "right_length_bits":s.right_length.to_bits()})
            })
            .collect::<Vec<_>>()
    });
    println!(
        "{}",
        json!({"qualified":false,"iterations":iterations,
        "rows":result.names.iter().zip(&result.sequences).map(|(name,seq)|
            json!({"id":name,"aligned":String::from_utf8(seq.clone()).expect("validated ASCII")})).collect::<Vec<_>>(),
        "distances_bits":distances,"guide":guide,"score_bits":result.score.to_bits(),
        "steps":result.step_trace.iter().map(|s|json!({"left_size":s.clus1,"right_size":s.clus2,
            "width":s.width,"score_bits":s.score.to_bits()})).collect::<Vec<_>>() })
    );
    Ok(())
}
