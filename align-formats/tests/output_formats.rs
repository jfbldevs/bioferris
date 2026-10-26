//! Deterministic interchange-format serialization tests.

use std::io::Cursor;

use align_core::{
    AlignmentOperation, Dna, DnaSequence, PairwiseAlignment, ProteinSequence, RnaSequence,
};
use align_fast::MyersPattern;
use align_formats::{
    FormatError, PhylipLayout, PhylipNameStyle, PhylipOptions, write_clustal_to,
    write_myers_tsv_to, write_nexus_to, write_pairwise_fasta_to, write_phylip_to,
    write_phylip_with_options_to, write_poa_dot_to, write_poa_gfa_to, write_stockholm_to,
    write_text_report_to,
};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};
use align_poa::{PoaGraph, PoastaAffineCosts, PoastaPoa};

fn records() -> Vec<SequenceRecord<Dna>> {
    [("reference one", "ACGT"), ("allele", "AGT")]
        .into_iter()
        .map(|(identifier, sequence)| {
            SequenceRecord::new(identifier, DnaSequence::try_from(sequence).unwrap()).unwrap()
        })
        .collect()
}

fn graph() -> PoaGraph<Dna> {
    PoastaPoa::new(PoastaAffineCosts::new(4, 6, 2).unwrap())
        .align(&records())
        .unwrap()
        .into_graph_and_alignment()
        .0
}

fn dna_alignment(identifiers: [&str; 2]) -> MultipleAlignment<Dna> {
    let first = DnaSequence::try_from("ACGT").unwrap();
    let second = DnaSequence::try_from("AGT").unwrap();
    MultipleAlignment::try_new(vec![
        AlignedRow::new(identifiers[0], first.iter().map(Some).collect()).unwrap(),
        AlignedRow::new(
            identifiers[1],
            vec![second.get(0), None, second.get(1), second.get(2)],
        )
        .unwrap(),
    ])
    .unwrap()
}

#[test]
fn pairwise_fasta_is_aligned_and_validates_identifiers() {
    let first = DnaSequence::try_from("ACGT").unwrap();
    let second = DnaSequence::try_from("AGT").unwrap();
    let alignment = PairwiseAlignment::try_new(
        1,
        0..4,
        0..3,
        vec![
            AlignmentOperation::Match,
            AlignmentOperation::Deletion,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ],
    )
    .unwrap();
    let mut output = Vec::new();
    write_pairwise_fasta_to(&mut output, "target", &first, "query", &second, &alignment).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        ">target\nACGT\n>query\nA-GT\n"
    );
    assert!(matches!(
        write_pairwise_fasta_to(Vec::new(), "", &first, "query", &second, &alignment),
        Err(FormatError::InvalidIdentifier { .. })
    ));
}

#[test]
fn aligned_fasta_wraps_at_eighty_columns_without_losing_symbols() {
    let source = "A".repeat(81);
    let first = DnaSequence::try_from(source.as_str()).unwrap();
    let second = DnaSequence::try_from(source.as_str()).unwrap();
    let alignment =
        PairwiseAlignment::try_new(81, 0..81, 0..81, vec![AlignmentOperation::Match; 81]).unwrap();
    let mut output = Vec::new();
    write_pairwise_fasta_to(&mut output, "first", &first, "second", &second, &alignment).unwrap();
    let output = String::from_utf8(output).unwrap();
    let lines = output.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], ">first");
    assert_eq!(lines[1].len(), 80);
    assert_eq!(lines[2], "A");
    assert_eq!(lines[3], ">second");
    assert_eq!(lines[4].len(), 80);
    assert_eq!(lines[5], "A");
}

#[test]
fn clustal_has_header_rows_and_literal_conservation() {
    let result = PoastaPoa::new(PoastaAffineCosts::new(4, 6, 2).unwrap())
        .align(&records())
        .unwrap();
    let mut output = Vec::new();
    write_clustal_to(&mut output, result.alignment()).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "CLUSTAL multiple sequence alignment\n\nreference one  ACGT\nallele         A-GT\n               * **\n\n"
    );
}

#[test]
fn clustal_wraps_rectangular_alignments_at_sixty_columns() {
    let source = "A".repeat(61);
    let long_records = [
        SequenceRecord::new("first", DnaSequence::try_from(source.as_str()).unwrap()).unwrap(),
        SequenceRecord::new("second", DnaSequence::try_from(source.as_str()).unwrap()).unwrap(),
    ];
    let result = PoastaPoa::new(PoastaAffineCosts::new(4, 6, 2).unwrap())
        .align(&long_records)
        .unwrap();
    let mut output = Vec::new();
    write_clustal_to(&mut output, result.alignment()).unwrap();
    let output = String::from_utf8(output).unwrap();
    assert!(output.contains(&format!("first   {}\n", "A".repeat(60))));
    assert!(output.contains("first   A\nsecond  A\n"));
}

#[test]
fn myers_tsv_uses_native_exclusive_endpoints() {
    let pattern = DnaSequence::try_from("A").unwrap();
    let text = DnaSequence::try_from("AAAA").unwrap();
    let matches = MyersPattern::new(&pattern)
        .unwrap()
        .find_all(&text, 0)
        .unwrap();
    let mut output = Vec::new();
    write_myers_tsv_to(&mut output, &matches).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "rank\ttext_end\tdistance\n1\t1\t0\n2\t2\t0\n3\t3\t0\n4\t4\t0\n"
    );
}

#[test]
fn poasta_gfa_and_dot_preserve_stable_graph_identity() {
    let graph = graph();
    let mut gfa = Vec::new();
    write_poa_gfa_to(&mut gfa, &graph).unwrap();
    let gfa = String::from_utf8(gfa).unwrap();
    assert!(gfa.starts_with("H\tVN:Z:1.0\nS\tn0\tA\tCL:i:0\n"));
    assert!(gfa.contains("P\treference%20one\t"));
    assert!(gfa.contains("RC:i:"));

    let mut dot = Vec::new();
    write_poa_dot_to(&mut dot, &graph).unwrap();
    let dot = String::from_utf8(dot).unwrap();
    assert!(dot.starts_with("digraph POASTA {\n  rankdir=LR;\n"));
    assert!(dot.contains("n0 [label=\"n0\\nA\\nclass 0\"]"));
    assert!(dot.ends_with("}\n"));
}

#[test]
fn text_reports_stream_and_end_in_exactly_one_newline() {
    for report in ["", "report", "report\n", "report\n\n\n"] {
        let mut output = Cursor::new(Vec::new());
        write_text_report_to(&mut output, &report).unwrap();
        let expected = if report.is_empty() { "\n" } else { "report\n" };
        assert_eq!(String::from_utf8(output.into_inner()).unwrap(), expected);
    }
}

#[test]
fn text_report_writer_preserves_internal_blank_lines_while_normalizing_the_end() {
    let mut output = Cursor::new(Vec::new());
    write_text_report_to(&mut output, &"header\n\nbody\n\n").unwrap();
    assert_eq!(
        String::from_utf8(output.into_inner()).unwrap(),
        "header\n\nbody\n"
    );
}

#[test]
fn relaxed_sequential_phylip_has_exact_dimensions_and_rows() {
    let alignment = dna_alignment(["taxon_a", "taxon_b"]);
    let mut output = Vec::new();
    write_phylip_to(&mut output, &alignment).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "2 4\ntaxon_a ACGT\ntaxon_b A-GT\n"
    );
}

#[test]
fn strict_interleaved_phylip_uses_fixed_names_and_equal_blocks() {
    let alignment = dna_alignment(["taxon_a", "taxon_b"]);
    let options = PhylipOptions::new()
        .with_layout(PhylipLayout::Interleaved)
        .with_name_style(PhylipNameStyle::Strict)
        .with_block_width(2)
        .unwrap();
    let mut output = Vec::new();
    write_phylip_with_options_to(&mut output, &alignment, options).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "2 4\ntaxon_a   AC\ntaxon_b   A-\n\n          GT\n          GT\n"
    );
}

#[test]
fn phylip_rejects_lossy_identifiers_and_zero_blocks() {
    let long = dna_alignment(["eleven_chars", "taxon_b"]);
    let strict = PhylipOptions::new().with_name_style(PhylipNameStyle::Strict);
    assert!(matches!(
        write_phylip_with_options_to(Vec::new(), &long, strict),
        Err(FormatError::InvalidFormatIdentifier { .. })
    ));

    let whitespace = dna_alignment(["taxon a", "taxon_b"]);
    assert!(matches!(
        write_phylip_to(Vec::new(), &whitespace),
        Err(FormatError::InvalidFormatIdentifier { .. })
    ));
    assert!(matches!(
        PhylipOptions::new().with_block_width(0),
        Err(FormatError::InvalidLineWidth {
            format: "PHYLIP",
            value: 0
        })
    ));
}

#[test]
fn nexus_quotes_taxa_and_selects_the_typed_datatype() {
    let dna = dna_alignment(["sample one", "O'Brien"]);
    let mut output = Vec::new();
    write_nexus_to(&mut output, &dna).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "#NEXUS\n\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=4;\n  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n  MATRIX\n    'sample one' ACGT\n    'O''Brien' A-GT\n  ;\nEND;\n"
    );

    let rna = MultipleAlignment::try_new(vec![
        AlignedRow::new(
            "rna",
            RnaSequence::try_from("ACGU")
                .unwrap()
                .iter()
                .map(Some)
                .collect(),
        )
        .unwrap(),
    ])
    .unwrap();
    let protein = MultipleAlignment::try_new(vec![
        AlignedRow::new(
            "protein",
            ProteinSequence::try_from("MKT")
                .unwrap()
                .iter()
                .map(Some)
                .collect(),
        )
        .unwrap(),
    ])
    .unwrap();
    let mut rna_output = Vec::new();
    let mut protein_output = Vec::new();
    write_nexus_to(&mut rna_output, &rna).unwrap();
    write_nexus_to(&mut protein_output, &protein).unwrap();
    assert!(
        String::from_utf8(rna_output)
            .unwrap()
            .contains("DATATYPE=RNA")
    );
    assert!(
        String::from_utf8(protein_output)
            .unwrap()
            .contains("DATATYPE=PROTEIN")
    );
}

#[test]
fn stockholm_has_the_required_envelope_and_rejects_whitespace_names() {
    let alignment = dna_alignment(["taxon_a", "b"]);
    let mut output = Vec::new();
    write_stockholm_to(&mut output, &alignment).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "# STOCKHOLM 1.0\ntaxon_a ACGT\nb       A-GT\n//\n"
    );

    let invalid = dna_alignment(["taxon a", "taxon_b"]);
    assert!(matches!(
        write_stockholm_to(Vec::new(), &invalid),
        Err(FormatError::InvalidFormatIdentifier { .. })
    ));

    for reserved in ["#=GS", "//"] {
        let invalid = dna_alignment([reserved, "taxon_b"]);
        assert!(matches!(
            write_stockholm_to(Vec::new(), &invalid),
            Err(FormatError::InvalidFormatIdentifier { .. })
        ));
    }
}
