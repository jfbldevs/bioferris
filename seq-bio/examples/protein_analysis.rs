//! Canonical peptide properties with explicit pH, pI boundaries and input errors.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let peptide = protein("INGAR")?;
    let properties = peptide.protein_properties()?;

    println!("Protein: {peptide}");
    println!("Mass (Da): {:.4}", properties.molecular_weight());
    println!("Charge at pH 7: {:.6}", properties.charge());
    println!("Isoelectric point: {:.6}", properties.isoelectric_point());
    println!("Aromaticity: {:.6}", properties.aromaticity());
    println!("GRAVY: {:.6}", properties.gravy());

    let options = ProteinPropertiesOptions { ph: 4.53 };
    let custom = protein("PETER")?.protein_properties_with(options)?;
    println!("PETER charge at pH 4.53: {:.6}", custom.charge());

    let acidic = protein("DD")?.protein_properties()?;
    println!("DD pI boundary: {:?}", acidic.pi_estimate().boundary());

    if let Err(error) = protein("AX")?.protein_properties() {
        println!("Rejected AX: {error}");
    }

    let out_of_range = ProteinPropertiesOptions { ph: 15.0 };
    if let Err(error) = protein("INGAR")?.protein_properties_with(out_of_range) {
        println!("Rejected pH 15 ({}): {error}", error.code());
    }

    Ok(())
}
