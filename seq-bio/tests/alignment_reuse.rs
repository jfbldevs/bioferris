//! Integration checks for unchanged alignment types, records and file formats.

use bioferris::prelude::*;

#[test]
fn transformed_fasta_records_use_the_unchanged_multiple_alignment_api() {
    let original = bioferris::io::parse_fasta::<Dna>(
        b">reference full header\nACGTACGT\n>allele\nACGTCGT\n>third\nACGTTCGT\n".as_slice(),
    )
    .unwrap();
    assert_eq!(original[0].identifier(), "reference full header");
    let transformed: Vec<_> = original
        .iter()
        .map(|record| {
            SequenceRecord::new(
                record.identifier(),
                record.sequence().reverse_complement().unwrap(),
            )
            .unwrap()
        })
        .collect();
    let result = Multiple::mafft(MafftStrategy::FftNs1)
        .align(&transformed)
        .unwrap();
    let direct = bioferris::align::Multiple::mafft(MafftStrategy::FftNs1)
        .align(&transformed)
        .unwrap();
    assert_eq!(result.alignment(), direct.alignment());
    assert_eq!(result.alignment().sequence_count(), 3);

    let mut fasta = Vec::new();
    bioferris::io::write_fasta_to(&mut fasta, &transformed).unwrap();
    assert_eq!(
        bioferris::io::parse_fasta::<Dna>(fasta.as_slice()).unwrap(),
        transformed
    );
}

#[test]
fn translation_products_use_the_existing_pairwise_types_and_results() {
    let proteins: Vec<_> = ["ATGGCTAAATAG", "ATGGTTAAATAG"]
        .into_iter()
        .map(|input| {
            DnaSequence::try_from(input)
                .unwrap()
                .translate_standard(TranslationOptions {
                    mode: TranslationMode::Cds,
                    ..TranslationOptions::default()
                })
                .unwrap()
                .to_protein()
                .unwrap()
        })
        .collect();
    let aligner = Pairwise::needleman_wunsch().scores(2, -1).gap(-2);
    let result = aligner.align(&proteins[0], &proteins[1]).unwrap();
    assert_eq!(result.score(), 3);
    assert_eq!(result.first_range(), 0..3);
    assert_eq!(result.second_range(), 0..3);
    assert_eq!(
        result,
        bioferris::align::Pairwise::needleman_wunsch()
            .scores(2, -1)
            .gap(-2)
            .align(&proteins[0], &proteins[1])
            .unwrap()
    );
}

#[test]
fn graph_api_is_available_without_a_new_graph_implementation() {
    let records = bioferris::io::parse_fasta::<Dna>(b">a\nACGT\n>b\nAGT\n".as_slice()).unwrap();
    let graph = Graph::poasta().run(&records).unwrap();
    assert_eq!(graph.result().alignment().sequence_count(), 2);
}
