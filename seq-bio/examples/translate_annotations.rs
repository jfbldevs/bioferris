//! Translate an annotated CDS and check its source protein with explicit policies.

use bioferris::{io::annotated::parse_genbank, prelude::*};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const INPUT: &str = r#"LOCUS       DEMO 15 bp DNA circular SYN 15-OCT-2026
DEFINITION  Original coding record.
ACCESSION   DEMO
VERSION     DEMO.2
FEATURES             Location/Qualifiers
     source          1..15
                     /organism="Synthetic construct"
     CDS             4..12
                     /gene="demo"
                     /codon_start=1
                     /transl_table=1
                     /translation="MA"
ORIGIN
        1 aaaATGGCTT AAccc
//
"#;
    let records = parse_genbank::<Dna>(INPUT)?;
    let record = records.first().ok_or("the example requires one record")?;
    let feature = record
        .features()
        .iter()
        .find(|feature| feature.key == "CDS")
        .ok_or("the example requires a CDS")?;
    let product = feature.translate_cds(record.record().sequence())?;
    println!("Record: {}", record.record().id());
    println!("Location: {}", feature.location.raw());
    println!(
        "Extracted DNA: {}",
        feature.extract(record.record().sequence())?
    );
    println!("Genetic code: {}", product.genetic_code().id());
    println!("Codon start: {}", product.codon_start());
    println!("Protein: {}", product.peptide());
    println!("Reference: {:?}", product.reference_status());

    let mut mismatch = feature.clone();
    let reference = mismatch
        .qualifiers
        .iter_mut()
        .find(|qualifier| qualifier.name == "translation")
        .ok_or("the example requires a reference translation")?;
    reference.value = Some("MV".to_owned());
    match mismatch.translate_cds(record.record().sequence()) {
        Err(error) => println!("Mismatch: {error}"),
        Ok(_) => return Err("a mismatching reference must fail".into()),
    }
    let recovered = mismatch.translate_cds_with(
        record.record().sequence(),
        CdsTranslationOptions {
            reference_policy: CdsReferencePolicy::Ignore,
            ..CdsTranslationOptions::default()
        },
    )?;
    println!(
        "Explicitly ignored reference: {:?}",
        recovered.reference_status()
    );
    println!("Computed protein: {}", recovered.peptide());
    Ok(())
}
