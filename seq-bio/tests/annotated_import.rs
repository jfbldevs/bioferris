//! Independent annotated-file fixtures, extraction boundaries and documented results.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use bioferris::io::annotated::{
    Location, LocationExpression, parse_embl, parse_genbank, read_embl, read_genbank,
};
use bioferris::io::{IoErrorKind, ReadLimits};
use bioferris::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};
use flate2::{Compression, write::GzEncoder};

const GENBANK: &str = r#"LOCUS       SYNTHETIC 18 bp DNA linear SYN 15-OCT-2026
DEFINITION  Example with a wrapped
            description.
ACCESSION   TEST0001 ALT0002
VERSION     TEST0001.2
DBLINK      BioProject: PRJNA000000
COMMENT     Unknown metadata is retained.
FEATURES             Location/Qualifiers
     source          1..18
                     /organism="Synthetic construct"
     CDS             complement(join(1..3,
                     10..12))
                     /gene="demo"
                     /note="first ""quoted""
                     second /literal"
                     /db_xref="DB:1"
                     /db_xref="DB:2"
                     /codon_start=2
                     /transl_table=11
                     /pseudo
     misc_feature    <1..>18
     misc_feature    OTHER.1:1..3
ORIGIN
        1 aTggcttaac ccgggttt
//
"#;

const EMBL: &str = r#"ID   TEST0001; SV 2; linear; genomic DNA; STD; SYN; 18 BP.
XX
AC   TEST0001; ALT0002;
DE   Example with a wrapped
DE   description.
DR   TEST; retained.
FH   Key             Location/Qualifiers
FH
FT   source          1..18
FT                   /organism="Synthetic construct"
FT   CDS             complement(join(1..3,10..12))
FT                   /gene="demo"
FT                   /transl_table=11
XX
SQ   Sequence 18 BP; 3 A; 4 C; 5 G; 6 T; 0 other;
     aTggcttaac ccgggttt 18
//
"#;

fn genbank(sequence: &str, length: usize, molecule: &str) -> String {
    format!(
        "LOCUS       TEST {length} bp {molecule} linear SYN 15-OCT-2026\nORIGIN\n        1 {sequence}\n//\n"
    )
}

fn embl(sequence: &str, length: usize, molecule: &str) -> String {
    format!(
        "ID   TEST; SV 1; linear; {molecule}; STD; SYN; {length} BP.\nSQ   Sequence {length} BP;\n     {sequence} {length}\n//\n"
    )
}

#[test]
fn genbank_retains_metadata_repeated_qualifiers_and_case() {
    let entries = parse_genbank::<Dna>(GENBANK).unwrap();
    assert_eq!(entries.len(), 1);
    let entry = &entries[0];
    assert_eq!(entry.locus(), "SYNTHETIC");
    assert_eq!(entry.declared_molecule(), "DNA");
    assert_eq!(entry.record().id(), "TEST0001.2");
    assert_eq!(
        entry.record().description(),
        Some("Example with a wrapped description.")
    );
    assert_eq!(entry.accessions(), ["TEST0001", "ALT0002"]);
    assert_eq!(entry.record().source_sequence(), "aTggcttaacccgggttt");
    assert_eq!(entry.record().sequence().to_string(), "ATGGCTTAACCCGGGTTT");
    assert!(entry.record().qualities().is_none());
    assert_eq!(entry.raw_record(), GENBANK);
    assert!(
        entry
            .metadata()
            .iter()
            .any(|field| field.key == "DBLINK" && field.value == "BioProject: PRJNA000000")
    );
    assert!(
        entry
            .metadata()
            .iter()
            .any(|field| field.key == "DEFINITION"
                && field.continuation
                && field.raw_line == "            description.")
    );
    let cds = &entry.features()[1];
    assert_eq!(cds.key, "CDS");
    assert_eq!(cds.location.raw(), "complement(join(1..3,\n10..12))");
    assert_eq!(
        cds.extract(entry.record().sequence()).unwrap().to_string(),
        "GGGCAT"
    );
    assert_eq!(
        cds.qualifiers[1].value.as_deref(),
        Some("first \"quoted\"\nsecond /literal")
    );
    assert_eq!(cds.qualifiers[2].name, "db_xref");
    assert_eq!(cds.qualifiers[2].value.as_deref(), Some("DB:1"));
    assert_eq!(cds.qualifiers[3].value.as_deref(), Some("DB:2"));
    assert_eq!(cds.qualifiers[4].value.as_deref(), Some("2"));
    assert_eq!(cds.qualifiers[5].value.as_deref(), Some("11"));
    assert_eq!(cds.qualifiers[6].value, None);
    assert_eq!(cds.qualifiers[6].raw, "/pseudo");
    assert_eq!(cds.raw_lines.len(), 10);
}

#[test]
fn modern_embl_preserves_metadata_and_matches_independent_genbank_sequence() {
    let entries = parse_embl::<Dna>(EMBL).unwrap();
    let entry = &entries[0];
    assert_eq!(entry.locus(), "TEST0001");
    assert_eq!(entry.declared_molecule(), "genomic DNA");
    assert_eq!(entry.record().id(), "TEST0001.2");
    assert_eq!(
        entry.record().description(),
        Some("Example with a wrapped description.")
    );
    assert_eq!(entry.accessions(), ["TEST0001", "ALT0002"]);
    assert_eq!(entry.raw_record(), EMBL);
    assert!(
        entry
            .metadata()
            .iter()
            .any(|field| field.key == "DR" && field.value == "TEST; retained.")
    );
    assert_eq!(entry.features().len(), 2);
    assert_eq!(
        entry.features()[1]
            .extract(entry.record().sequence())
            .unwrap()
            .to_string(),
        "GGGCAT"
    );
    assert_eq!(
        entry.record().sequence(),
        parse_genbank::<Dna>(GENBANK).unwrap()[0]
            .record()
            .sequence()
    );
}

#[test]
fn exact_points_ranges_and_complemented_joins_have_checked_coordinates() {
    let dna = DnaSequence::from_ascii("ATGGCTTAACCCGGGTTT").unwrap();
    assert_eq!(
        Location::parse("1..3").unwrap().expression(),
        &LocationExpression::Range(0..3)
    );
    assert_eq!(
        Location::parse("18").unwrap().expression(),
        &LocationExpression::Range(17..18)
    );
    for (text, expected) in [
        ("1", "A"),
        ("18", "T"),
        ("join(1..3,10..12)", "ATGCCC"),
        ("complement(join(1..3,10..12))", "GGGCAT"),
        ("join(complement(10..12),complement(1..3))", "GGGCAT"),
        ("join(16..18,1..3)", "TTTATG"),
        ("complement(complement(1..3))", "ATG"),
        ("join(1..3,1..3)", "ATGATG"),
    ] {
        assert_eq!(
            Location::parse(text)
                .unwrap()
                .extract(&dna)
                .unwrap()
                .to_string(),
            expected,
            "{text}"
        );
    }
    let rna = RnaSequence::from_ascii("AUGCUU").unwrap();
    assert_eq!(
        Location::parse("complement(1..3)")
            .unwrap()
            .extract(&rna)
            .unwrap()
            .to_string(),
        "CAU"
    );
    let protein = ProteinSequence::from_ascii("MAW").unwrap();
    assert_eq!(
        Location::parse("join(1,3)")
            .unwrap()
            .extract(&protein)
            .unwrap()
            .to_string(),
        "MW"
    );
    assert_eq!(
        Location::parse("complement(1..3)")
            .unwrap()
            .extract(&protein)
            .unwrap_err()
            .kind(),
        IoErrorKind::Unsupported
    );
}

#[test]
fn every_exact_subrange_and_complement_matches_independent_ascii_oracle() {
    let text = "ACGTRYSWKMBDHVN";
    let sequence = DnaSequence::from_ascii(text).unwrap();
    let complements = b"TGCAYRSWMKVHDBN";
    for start in 0..text.len() {
        for end in start + 1..=text.len() {
            let location = format!("{}..{end}", start + 1);
            assert_eq!(
                Location::parse(&location)
                    .unwrap()
                    .extract(&sequence)
                    .unwrap()
                    .to_string(),
                text[start..end]
            );
            let expected: String = complements[start..end]
                .iter()
                .rev()
                .map(|byte| char::from(*byte))
                .collect();
            assert_eq!(
                Location::parse(&format!("complement({location})"))
                    .unwrap()
                    .extract(&sequence)
                    .unwrap()
                    .to_string(),
                expected
            );
        }
    }
}

#[test]
fn fuzzy_remote_between_base_and_order_locations_are_retained_not_extracted() {
    let sequence = DnaSequence::from_ascii("ATGC").unwrap();
    for text in [
        "<1..>4",
        "OTHER.2:1..3",
        "2^3",
        "order(1,4)",
        "one-of(1,2)",
        "1.3",
        "join(1,OTHER.2:3)",
    ] {
        let location = Location::parse(text).unwrap();
        assert_eq!(location.raw(), text);
        let error = location.extract(&sequence).unwrap_err();
        assert_eq!(error.kind(), IoErrorKind::Unsupported, "{text}");
        assert_eq!(error.format(), "insdc-location");
    }
    let entries = parse_genbank::<Dna>(GENBANK).unwrap();
    for feature in &entries[0].features()[2..] {
        assert_eq!(
            feature
                .extract(entries[0].record().sequence())
                .unwrap_err()
                .code(),
            "unsupported_format_feature"
        );
    }
}

#[test]
fn invalid_exact_locations_and_excessive_nesting_are_rejected() {
    for text in [
        "",
        "0",
        "0..3",
        "3..2",
        "1..",
        "..3",
        "1...3",
        "1..2..3",
        "join(1)",
        "join(1,)",
        "complement(1,2)",
        "join(1,2",
        "1)",
        "(1",
        "1..999999999999999999999999999999999999999",
    ] {
        assert_eq!(
            Location::parse(text).unwrap_err().kind(),
            IoErrorKind::InvalidStructure,
            "{text}"
        );
    }
    let deeply_nested = format!("{}1{}", "complement(".repeat(65), ")".repeat(65));
    assert_eq!(
        Location::parse(&deeply_nested).unwrap_err().kind(),
        IoErrorKind::LimitExceeded
    );
    let sequence = DnaSequence::from_ascii("ATG").unwrap();
    assert_eq!(
        Location::parse("1..4")
            .unwrap()
            .extract(&sequence)
            .unwrap_err()
            .kind(),
        IoErrorKind::LengthMismatch
    );
    let too_long_feature = GENBANK.replace("source          1..18", "source          1..19");
    assert_eq!(
        parse_genbank::<Dna>(&too_long_feature).unwrap_err().kind(),
        IoErrorKind::LengthMismatch
    );
}

#[test]
fn multi_record_crlf_no_final_newline_and_empty_documents_are_explicit() {
    let first = genbank("atg", 3, "DNA").replace('\n', "\r\n");
    let second = genbank("CC", 2, "DNA");
    let input = format!("\r\n{first}\n{}", second.trim_end());
    let entries = parse_genbank::<Dna>(&input).unwrap();
    assert_eq!(entries.len(), 2);
    assert_eq!(entries[0].raw_record(), first);
    assert_eq!(entries[1].raw_record(), second.trim_end());
    assert_eq!(entries[1].record().sequence().to_string(), "CC");
    assert!(parse_genbank::<Dna>(" \n\r\n").unwrap().is_empty());
    assert!(parse_embl::<Dna>("").unwrap().is_empty());
    assert_eq!(
        parse_embl::<Dna>(&format!("{EMBL}{EMBL}")).unwrap().len(),
        2
    );
}

#[test]
fn truncation_and_nested_record_headers_fail_with_record_line_context() {
    let truncated = genbank("ATG", 3, "DNA").replace("//\n", "");
    let combined = format!("{}\n{truncated}", genbank("ATG", 3, "DNA"));
    let error = parse_genbank::<Dna>(&combined).unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::Truncated);
    assert_eq!(error.record(), Some(2));
    assert_eq!(error.line(), Some(6));
    assert_eq!(error.format(), "genbank");
    assert_eq!(
        parse_genbank::<Dna>(&format!("{truncated}{GENBANK}"))
            .unwrap_err()
            .kind(),
        IoErrorKind::Truncated
    );
    assert_eq!(
        parse_embl::<Dna>(&EMBL.replace("//\n", ""))
            .unwrap_err()
            .kind(),
        IoErrorKind::Truncated
    );
    assert_eq!(
        parse_genbank::<Dna>(">TEST\nATG\n").unwrap_err().kind(),
        IoErrorKind::InvalidHeader
    );
}

#[test]
fn lengths_counters_and_literal_sequence_are_checked() {
    assert_eq!(
        parse_genbank::<Dna>(&genbank("ATG", 4, "DNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::LengthMismatch
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("ATG", 3, "DNA").replace("1 ATG", "2 ATG"))
            .unwrap_err()
            .kind(),
        IoErrorKind::LengthMismatch
    );
    assert_eq!(
        parse_embl::<Dna>(&embl("ATG", 3, "DNA").replace("ATG 3", "ATG 2"))
            .unwrap_err()
            .kind(),
        IoErrorKind::LengthMismatch
    );
    assert_eq!(
        parse_embl::<Dna>(&embl("ATG", 3, "DNA").replace("Sequence 3", "Sequence 4"))
            .unwrap_err()
            .kind(),
        IoErrorKind::LengthMismatch
    );
    let contig =
        "LOCUS       TEST 3 bp DNA linear SYN 15-OCT-2026\nCONTIG      join(OTHER.1:1..3)\n//\n";
    assert_eq!(
        parse_genbank::<Dna>(contig).unwrap_err().kind(),
        IoErrorKind::Unsupported
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("A-G", 3, "DNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("A1G", 3, "DNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("ATG", 3, "DNA").replace("1 ATG", "ATG"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
}

#[test]
fn typed_molecule_contract_never_converts_rna_storage_t_into_u() {
    let rna = genbank("aug", 3, "mRNA");
    assert_eq!(
        parse_genbank::<Rna>(&rna).unwrap()[0]
            .record()
            .sequence()
            .to_string(),
        "AUG"
    );
    assert_eq!(
        parse_genbank::<Dna>(&rna).unwrap_err().kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_genbank::<Rna>(&genbank("atg", 3, "mRNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_embl::<Rna>(&embl("aug", 3, "mRNA")).unwrap()[0]
            .record()
            .sequence()
            .to_string(),
        "AUG"
    );
    assert_eq!(
        parse_embl::<Rna>(&embl("atg", 3, "mRNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("ATG", 3, "notDNA"))
            .unwrap_err()
            .kind(),
        IoErrorKind::Unsupported
    );
    assert_eq!(
        parse_genbank::<Dna>(&genbank("ATG", 3, "ds-DNA")).unwrap()[0].declared_molecule(),
        "ds-DNA"
    );
    let protein = "LOCUS       PEP 3 aa linear SYN 15-OCT-2026\nORIGIN\n        1 mAw\n//\n";
    assert_eq!(
        parse_genbank::<Protein>(protein).unwrap()[0]
            .record()
            .sequence()
            .to_string(),
        "MAW"
    );
    assert_eq!(
        parse_genbank::<Dna>(protein).unwrap_err().kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_genbank::<Protein>(&genbank("MAW", 3, "protein"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidHeader
    );
    assert_eq!(
        parse_embl::<Protein>(&embl("MAW", 3, "protein"))
            .unwrap_err()
            .kind(),
        IoErrorKind::Unsupported
    );
    assert_eq!(
        parse_genbank::<Protein>(&protein.replace("mAw", "MA*"))
            .unwrap_err()
            .kind(),
        IoErrorKind::InvalidSequence
    );
}

#[test]
fn malformed_headers_features_and_qualifiers_fail_without_guessing() {
    let base = genbank("ATG", 3, "DNA");
    for bad in [
        base.replace("3 bp", "x bp"),
        base.replace("3 bp", "3 nt"),
        base.replace("3 bp DNA", "3 bp"),
    ] {
        assert!(parse_genbank::<Dna>(&bad).is_err());
    }
    let old_embl = "ID   TEST standard; DNA; SYN; 3 BP.\nSQ   Sequence 3 BP;\n     atg 3\n//\n";
    assert_eq!(
        parse_embl::<Dna>(old_embl).unwrap_err().kind(),
        IoErrorKind::Unsupported
    );
    let bad_qualifier = GENBANK.replace("/gene=\"demo\"", "/gene=\"demo\"junk");
    assert_eq!(
        parse_genbank::<Dna>(&bad_qualifier).unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
    let unterminated = "LOCUS       TEST 3 bp DNA linear SYN 15-OCT-2026\nFEATURES             Location/Qualifiers\n     CDS             1..3\n                     /note=\"missing end\nORIGIN\n        1 ATG\n//\n";
    assert_eq!(
        parse_genbank::<Dna>(unterminated).unwrap_err().kind(),
        IoErrorKind::Truncated
    );
    let orphan = GENBANK.replace(
        "     source          1..18",
        "                     /gene=\"orphan\"",
    );
    assert_eq!(
        parse_genbank::<Dna>(&orphan).unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
    let versions = GENBANK.replace(
        "VERSION     TEST0001.2",
        "VERSION     TEST0001.2\nVERSION     TEST0001.3",
    );
    assert_eq!(
        parse_genbank::<Dna>(&versions).unwrap_err().kind(),
        IoErrorKind::InvalidHeader
    );
    let unknown_key = GENBANK.replace("/gene=\"demo\"", "/new-key=\"retained\"");
    assert_eq!(
        parse_genbank::<Dna>(&unknown_key).unwrap()[0].features()[1].qualifiers[0].name,
        "new-key"
    );
    let invalid_key = GENBANK.replace("/gene=\"demo\"", "/123=\"invalid\"");
    assert_eq!(
        parse_genbank::<Dna>(&invalid_key).unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
}

#[test]
fn default_line_and_declared_record_limits_are_enforced() {
    let too_large = genbank("A", ReadLimits::default().max_record_symbols + 1, "DNA");
    assert_eq!(
        parse_genbank::<Dna>(&too_large).unwrap_err().kind(),
        IoErrorKind::LimitExceeded
    );
    let long_line = format!(
        "LOCUS       TEST 1 bp DNA\nCOMMENT     {}\nORIGIN\n        1 A\n//\n",
        "x".repeat(ReadLimits::default().max_line_bytes)
    );
    let error = parse_genbank::<Dna>(&long_line).unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::LimitExceeded);
    assert_eq!(error.line(), Some(2));
}

struct TempInput(PathBuf);

impl TempInput {
    fn new(bytes: &[u8]) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let serial = NEXT.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "seq-bio-annotated-{}-{serial}.data",
            std::process::id()
        ));
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        file.write_all(bytes).unwrap();
        Self(path)
    }
}

impl Drop for TempInput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[test]
fn path_readers_accept_plain_and_gzip_and_surface_late_corruption() {
    let plain = TempInput::new(GENBANK.as_bytes());
    assert_eq!(
        read_genbank::<Dna>(&plain.0).unwrap()[0].raw_record(),
        GENBANK
    );
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(EMBL.as_bytes()).unwrap();
    let mut compressed = encoder.finish().unwrap();
    let valid = TempInput::new(&compressed);
    assert_eq!(read_embl::<Dna>(&valid.0).unwrap()[0].raw_record(), EMBL);
    let crc_index = compressed.len() - 8;
    compressed[crc_index] ^= 1;
    let corrupt = TempInput::new(&compressed);
    assert_eq!(
        read_embl::<Dna>(&corrupt.0).unwrap_err().kind(),
        IoErrorKind::Io
    );
    let invalid_utf8 = TempInput::new(&[0xff]);
    assert_eq!(
        read_genbank::<Dna>(&invalid_utf8.0).unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
}

fn doc_block<'a>(document: &'a str, marker: &str, language: &str) -> &'a str {
    let marker = format!("<!-- {marker} -->");
    let remainder = document.split_once(&marker).unwrap().1;
    let fence = char::from(96).to_string().repeat(3);
    let opening = format!("{fence}{language}\n");
    remainder
        .split_once(&opening)
        .unwrap()
        .1
        .split_once(&fence)
        .unwrap()
        .0
}

#[test]
fn documented_genbank_embl_and_unsupported_location_outputs_match() {
    let docs = include_str!("../../docs/api/annotated-input-v1.md");
    let genbank_input = doc_block(docs, "annotated-genbank-input", "text");
    let embl_input = doc_block(docs, "annotated-embl-input", "text");
    let program = doc_block(docs, "annotated-program", "rust");
    assert!(program.contains(&format!("const GENBANK: &str = r#\"{genbank_input}\"#;")));
    assert!(program.contains(&format!("const EMBL: &str = r#\"{embl_input}\"#;")));
    let gb = parse_genbank::<Dna>(genbank_input).unwrap();
    let em = parse_embl::<Dna>(embl_input).unwrap();
    let feature = &gb[0].features()[0];
    let error = Location::parse("<1..9")
        .unwrap()
        .extract(gb[0].record().sequence())
        .unwrap_err();
    let output = format!(
        "GenBank: {} {}\nEMBL: {} {}\nCDS: {}\nReverse feature: {}\nBoundary: {}\n",
        gb[0].record().id(),
        gb[0].record().sequence(),
        em[0].record().id(),
        em[0].record().sequence(),
        feature.extract(gb[0].record().sequence()).unwrap(),
        Location::parse("complement(join(1..3,4..6))")
            .unwrap()
            .extract(gb[0].record().sequence())
            .unwrap(),
        error.code(),
    );
    assert_eq!(output, doc_block(docs, "annotated-output", "text"));
}

#[test]
fn public_ncbi_ena_dq067560_records_match_declared_biology_and_exact_exon() {
    const EXON_SIX: &str = "GTCCTGCCCCGCGAGGTCACCTGGTCAGAGTTGAGGGGAACCAGCGAGCAGAGTACATGGAGGATGGTAACACTCTGAGACAGAGTGTGCTCGTCCCCTATGAGCCTCCTCAG";
    let genbank = include_str!("fixtures/annotated-dq067560.1.gb");
    let embl = include_str!("fixtures/annotated-dq067560.1.embl");
    let gb = parse_genbank::<Dna>(genbank).unwrap();
    let em = parse_embl::<Dna>(embl).unwrap();
    assert_eq!(gb.len(), 1);
    assert_eq!(em.len(), 1);
    assert_eq!(gb[0].record().sequence(), em[0].record().sequence());
    assert_eq!(
        gb[0].raw_record(),
        genbank.trim_end_matches('\n').to_owned() + "\n"
    );
    assert_eq!(em[0].raw_record(), embl);
    // These expectations come from the source ID/SQ and feature table, not from
    // a parser snapshot. The literal exon was transcribed from rows 421/481.
    for entry in [&gb[0], &em[0]] {
        assert_eq!(entry.record().id(), "DQ067560.1");
        assert_eq!(entry.record().sequence().len(), 534);
        assert_eq!(
            entry.record().description(),
            Some("Salmo trutta p53 tumor suppressor (p53) gene, exons 5, 6 and partial cds.")
        );
        assert_eq!(entry.features().len(), 6);
        let sequence = entry.record().sequence().to_string();
        assert!(
            sequence.starts_with("TACTCGCCAGACCTGAACAAGTTGTTCTGCCAATTGGCAAAGACTTGTCCAGTTCAGATC")
        );
        for (symbol, expected) in [(b'A', 134), (b'C', 126), (b'G', 134), (b'T', 140)] {
            assert_eq!(
                sequence.bytes().filter(|byte| *byte == symbol).count(),
                expected
            );
        }
        let exon_six = entry
            .features()
            .iter()
            .find(|feature| {
                feature.key == "exon"
                    && feature.qualifiers.iter().any(|qualifier| {
                        qualifier.name == "number" && qualifier.value.as_deref() == Some("6")
                    })
            })
            .unwrap();
        assert_eq!(
            exon_six.location.expression(),
            &LocationExpression::Range(421..534)
        );
        assert_eq!(
            exon_six
                .extract(entry.record().sequence())
                .unwrap()
                .to_string(),
            EXON_SIX
        );
        let cds = entry
            .features()
            .iter()
            .find(|feature| feature.key == "CDS")
            .unwrap();
        assert_eq!(cds.location.raw(), "join(<1..187,422..>534)");
        assert!(
            cds.qualifiers
                .iter()
                .any(|qualifier| qualifier.name == "protein_id"
                    && qualifier.value.as_deref() == Some("AAY68035.1"))
        );
        assert_eq!(
            cds.extract(entry.record().sequence()).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
    }
}

#[test]
fn public_ncbi_x56734_rna_storage_limitation_is_an_explicit_rejection() {
    let input = include_str!("fixtures/annotated-x56734.1.gb");
    assert!(input.starts_with("LOCUS       X56734                  1859 bp    mRNA"));
    let wrong_type = parse_genbank::<Dna>(input).unwrap_err();
    assert_eq!(wrong_type.kind(), IoErrorKind::InvalidSequence);
    assert_eq!(wrong_type.line(), Some(1));
    assert!(wrong_type.message().contains("declared molecule mRNA"));
    let storage = parse_genbank::<Rna>(input).unwrap_err();
    assert_eq!(storage.kind(), IoErrorKind::InvalidSequence);
    assert!(storage.message().contains("outside the RNA alphabet"));
}

#[test]
fn public_ena_puc19_other_dna_molecule_maps_to_dna_without_source_rewriting() {
    let input = include_str!("fixtures/annotated-m77789.2.embl");
    assert!(input.starts_with("ID   M77789; SV 2; circular; other DNA; STD; SYN; 2686 BP."));
    let records = parse_embl::<Dna>(input).unwrap();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert_eq!(record.record().id(), "M77789.2");
    assert_eq!(record.declared_molecule(), "other DNA");
    assert_eq!(record.raw_record(), input);
    assert_eq!(record.record().sequence().len(), 2686);
    assert_eq!(
        record.record().description(),
        Some("Cloning vector pUC19, complete sequence.")
    );
    let sequence = record.record().sequence().to_string();
    // These counts are independently declared in the ENA SQ header.
    for (symbol, expected) in [(b'A', 659), (b'C', 685), (b'G', 675), (b'T', 667)] {
        assert_eq!(
            sequence.bytes().filter(|byte| *byte == symbol).count(),
            expected
        );
    }
    assert_eq!(
        parse_embl::<Rna>(input).unwrap_err().kind(),
        IoErrorKind::InvalidSequence
    );
    assert_eq!(
        parse_embl::<Dna>(&input.replacen("; other DNA;", "; OTHER DNA;", 1)).unwrap()[0]
            .declared_molecule(),
        "OTHER DNA"
    );
    for unsupported in ["other RNA", "unassigned DNA", "not DNA", "ss-other DNA"] {
        let altered = input.replacen("; other DNA;", &format!("; {unsupported};"), 1);
        assert_eq!(
            parse_embl::<Dna>(&altered).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
    }
}
