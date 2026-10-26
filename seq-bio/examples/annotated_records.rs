//! Checked annotated writing and derivation with explicit identity and loss policy.

use bioferris::io::annotated::{
    AnnotationChange, AnnotationLossPolicy, AnnotationTransformOptions, parse_embl, parse_genbank,
    write_embl, write_genbank, write_genbank_new,
};
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = std::fs::read_to_string("seq-bio/tests/fixtures/annotated-example.gb")?;
    let source = parse_genbank::<Dna>(&input)?;
    let source = &source[0];
    let options = AnnotationTransformOptions::new("LOCAL_cds")?.with_date("16-OCT-2026")?;
    let sliced = source.slice(3..12, options)?;
    let cds = sliced.record().features()[1].extract(sliced.record().record().sequence())?;
    println!("CDS extraction: {cds}");
    println!("Protein: {}", cds.translate_cds()?);
    println!(
        "Original unchanged: {} {} symbols",
        source.record().id(),
        source.record().sequence().len()
    );
    let mut bytes = Vec::new();
    write_genbank(&mut bytes, sliced.record())?;
    println!("GenBank slice:");
    print!("{}", String::from_utf8(bytes)?);

    let partial = source.slice(
        4..12,
        AnnotationTransformOptions::new("LOCAL_partial")?.with_date("16-OCT-2026")?,
    );
    match partial {
        Ok(_) => return Err("partial CDS unexpectedly accepted".into()),
        Err(error) => println!("Partial CDS refused: {}", error.code()),
    }

    let input = std::fs::read_to_string("seq-bio/tests/fixtures/annotated-example.embl")?;
    let entries = parse_embl::<Dna>(&input)?;
    let entry = &entries[0];
    match entry.reverse_complement(AnnotationTransformOptions::new("LOCAL_rc")?) {
        Ok(_) => return Err("unreviewed metadata unexpectedly accepted".into()),
        Err(error) => println!("Default metadata policy: {}", error.code()),
    }
    let reversed = entry.reverse_complement(
        AnnotationTransformOptions::new("LOCAL_rc")?
            .with_loss_policy(AnnotationLossPolicy::Discard),
    )?;
    for change in &reversed.report().changes {
        match change {
            AnnotationChange::MetadataChanged {
                key, after: None, ..
            } if key == "DR" => println!("Removed metadata: {key}"),
            AnnotationChange::FeatureRemoved { index, key, .. } => {
                println!("Removed feature: {index} {key}");
            }
            _ => {}
        }
    }
    let cds = reversed.record().features()[1].extract(reversed.record().record().sequence())?;
    println!("Reverse-complement CDS extraction: {cds}");
    let mut bytes = Vec::new();
    write_embl(&mut bytes, reversed.record())?;
    println!("EMBL reverse complement:");
    print!("{}", String::from_utf8(bytes)?);

    let mut arguments = std::env::args_os().skip(1);
    if let Some(path) = arguments.next() {
        if arguments.next().is_some() {
            return Err("usage: cargo bio annotated_records -- [NEW_GENBANK_PATH]".into());
        }
        match write_genbank_new(path, sliced.record()) {
            Ok(()) => println!("Saved new GenBank file."),
            Err(error) => {
                eprintln!("Save refused: {}", error.code());
                std::process::exit(1);
            }
        }
    }
    Ok(())
}
