//! Simple manual Rust workflows preserve all explicit typed contracts.
use bioferris::prelude::*;

#[test]
fn short_constructors_and_translation_shortcuts_match_explicit_options() {
    let sequence = dna("ttgaaataactg").unwrap();
    assert_eq!(sequence.to_string(), "TTGAAATAACTG");
    assert_eq!(sequence.translate().unwrap().as_str(), "LK*L");
    assert_eq!(sequence.translate_to_stop().unwrap().as_str(), "LK");
    let complete = rna("uugaaauaa").unwrap();
    assert_eq!(
        complete.translate_cds().unwrap().to_protein().unwrap(),
        protein("MK").unwrap()
    );
    for mode in [
        TranslationMode::Full,
        TranslationMode::ToStop,
        TranslationMode::Cds,
    ] {
        let simple = match mode {
            TranslationMode::Full => complete.translate(),
            TranslationMode::ToStop => complete.translate_to_stop(),
            TranslationMode::Cds => complete.translate_cds(),
        };
        assert_eq!(
            simple,
            complete.translate_standard(TranslationOptions {
                mode,
                ..TranslationOptions::default()
            })
        );
    }
}

#[test]
fn simple_entry_points_keep_validation_and_actionable_errors() {
    assert!(dna("ACGU").is_err());
    assert!(rna("ACTG").is_err());
    assert!(protein("M*").is_err());
    let partial = dna("ATGA").unwrap();
    assert_eq!(
        partial.translate().unwrap_err(),
        TranslationError::PartialCodon {
            length: 4,
            remainder: 1
        }
    );
    assert!(partial.translate_cds().is_err());
    assert!(partial.translate_to_stop().is_err());
    assert_eq!(
        dna(b"acgt").unwrap().reverse_complement().unwrap(),
        dna("ACGT").unwrap()
    );
}

#[test]
fn standard_result_is_not_shadowed_in_human_workflows() -> Result<(), Box<dyn std::error::Error>> {
    let product = dna("ATGGCTTAA")?.translate_cds()?.to_protein()?;
    let expected = protein("MA")?;
    let run = Pairwise::gotoh().run(&product, &expected)?;
    assert_eq!(run.result().first_range(), 0..2);
    Ok(())
}
