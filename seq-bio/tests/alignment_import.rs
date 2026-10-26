//! Independent alignment fixtures, loss checks and existing-writer compatibility.

use std::error::Error;
use std::fmt::Write as _;
use std::io::{self, BufReader, Read};

use bioferris::io::alignment::{
    AlignmentDocument, AlignmentFormat, AnnotationScope, parse_alignment, read_alignment,
};
use bioferris::io::{
    IoErrorKind, write_alignment_fasta_to, write_clustal_to, write_nexus_to,
    write_phylip_with_options_to, write_stockholm_to,
};
use bioferris::prelude::*;

type TestResult = Result<(), Box<dyn Error>>;

const ALIGNED_FASTA: &str = ">reference\nAC-G\n>sample\nACTG\n";

fn dna_document(input: &str, format: AlignmentFormat) -> AlignmentDocument<Dna> {
    parse_alignment(input.as_bytes(), format).expect("valid independent alignment fixture")
}

fn assert_rows(document: &AlignmentDocument<Dna>, expected: &[(&str, &str)]) {
    let actual: Vec<_> = document
        .rows()
        .iter()
        .map(|row| (row.identifier(), row.aligned()))
        .collect();
    assert_eq!(actual, expected);
}

#[test]
fn documented_fasta_input_output() -> TestResult {
    let document = parse_alignment::<Dna>(ALIGNED_FASTA.as_bytes(), AlignmentFormat::Fasta)?;
    let mut actual = String::new();
    writeln!(
        actual,
        "{} rows, {} columns",
        document.rows().len(),
        document.width()
    )?;
    for row in document.rows() {
        writeln!(actual, "{}: {}", row.identifier(), row.aligned())?;
    }
    let alignment = document.to_multiple_alignment()?;
    assert_eq!(alignment.rows()[0].to_ascii_with_gaps(), "AC-G");
    let guide = include_str!("../../docs/api/alignment-input-v1.md");
    let section = guide
        .split_once("<!-- alignment-fasta-recipe -->")
        .expect("recipe marker")
        .1;
    let section = section
        .split_once("<!-- /alignment-fasta-recipe -->")
        .expect("recipe end")
        .0;
    let input = section
        .split_once("```fasta\n")
        .expect("literal input")
        .1
        .split_once("```")
        .expect("input end")
        .0;
    let output = section
        .split_once("```text\n")
        .expect("literal output")
        .1
        .split_once("```")
        .expect("output end")
        .0;
    assert_eq!(input, ALIGNED_FASTA);
    assert_eq!(actual, output);
    print!("{actual}");
    Ok(())
}

fn documented_output(marker: &str) -> &str {
    let guide = include_str!("../../docs/api/alignment-input-v1.md");
    let section = guide
        .split_once(&format!("<!-- {marker} -->"))
        .expect("output marker")
        .1
        .split_once(&format!("<!-- /{marker} -->"))
        .expect("output end marker")
        .0;
    section
        .split_once("```text\n")
        .expect("text block")
        .1
        .split_once("```")
        .expect("block end")
        .0
}

#[test]
fn documented_format_and_loss_outputs() -> TestResult {
    let inputs = [
        (
            "CLUSTAL",
            AlignmentFormat::Clustal,
            "CLUSTAL W\n\na AC-G\nb ACTG\n",
        ),
        (
            "PHYLIP relaxed sequential",
            AlignmentFormat::Phylip {
                names: PhylipNameStyle::Relaxed,
                layout: PhylipLayout::Sequential,
            },
            "2 4\na AC\n-G\nb AC\nTG\n",
        ),
        (
            "PHYLIP strict interleaved",
            AlignmentFormat::Phylip {
                names: PhylipNameStyle::Strict,
                layout: PhylipLayout::Interleaved,
            },
            "2 4\na         AC\nb         AC\n\n-G\nTG\n",
        ),
        (
            "Stockholm",
            AlignmentFormat::Stockholm,
            "# STOCKHOLM 1.0\na AC-G\nb ACTG\n//\n",
        ),
        (
            "NEXUS",
            AlignmentFormat::Nexus,
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA;\nMATRIX\na AC-G\nb ACTG\n;\nEND;\n",
        ),
    ];
    let mut actual = String::new();
    for (label, format, input) in inputs {
        let document = parse_alignment::<Dna>(input.as_bytes(), format)?;
        let alignment = document.to_multiple_alignment()?;
        assert_eq!(alignment.rows()[0].to_ascii_with_gaps(), "AC-G");
        assert_eq!(alignment.rows()[1].to_ascii_with_gaps(), "ACTG");
        writeln!(
            actual,
            "{label}: {} rows, {} columns",
            document.rows().len(),
            document.width()
        )?;
    }
    assert_eq!(actual, documented_output("alignment-formats-output"));
    print!("{actual}");
    actual.clear();
    let input = b"# STOCKHOLM 1.0\n#=GF ID example\na AC.G\nb ACTG\n//\n";
    let document = parse_alignment::<Dna>(&input[..], AlignmentFormat::Stockholm)?;
    let error = document.to_multiple_alignment().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::DataLoss);
    writeln!(actual, "Strict conversion: {}", error.code())?;
    writeln!(
        actual,
        "Retained annotation: {}={}",
        document.annotations()[0].tag(),
        document.annotations()[0].value()
    )?;
    let alignment = document.to_multiple_alignment_discarding_annotations()?;
    writeln!(
        actual,
        "Explicit conversion: {}",
        alignment.rows()[0].to_ascii_with_gaps()
    )?;
    writeln!(actual, "Unchanged source: {}", document.rows()[0].aligned())?;
    let missing = b">a\nAC?G\n>b\nACTG\n";
    let document = parse_alignment::<Dna>(&missing[..], AlignmentFormat::Fasta)?;
    let error = document
        .to_multiple_alignment_discarding_annotations()
        .unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::DataLoss);
    writeln!(actual, "Missing-cell conversion: {}", error.code())?;
    assert_eq!(actual, documented_output("alignment-loss-output"));
    print!("{actual}");
    Ok(())
}

#[test]
fn fasta_preserves_headers_missing_data_and_all_gap_columns() {
    let document = dna_document(
        ">first full header\r\na c-?\r\n>second\r\nAC-N\r\n",
        AlignmentFormat::Fasta,
    );
    assert_rows(
        &document,
        &[("first full header", "AC-?"), ("second", "AC-N")],
    );
    assert_eq!(document.missing_symbol(), Some(b'?'));
    assert_eq!(document.gap_symbols(), b"-");
    assert_eq!(
        document.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    let all_gap = dna_document(">first\nA-C\n>second\nG-T\n", AlignmentFormat::Fasta);
    assert_eq!(all_gap.width(), 3);
    assert_eq!(
        all_gap.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
    assert_rows(&all_gap, &[("first", "A-C"), ("second", "G-T")]);
}

#[test]
fn fasta_rejects_malformed_nonrectangular_and_invalid_alphabets() {
    for (input, kind) in [
        ("", IoErrorKind::InvalidStructure),
        ("ACG\n", IoErrorKind::InvalidHeader),
        (">\nAC\n", IoErrorKind::InvalidHeader),
        (">a\nAC\n>a\nAC\n", IoErrorKind::InvalidStructure),
        (">a\nAC\n>b\nA\n", IoErrorKind::LengthMismatch),
        (">a\nAU\n", IoErrorKind::InvalidSequence),
        (">a\nA.\n", IoErrorKind::InvalidSequence),
        (">a\n;comment\nAC\n", IoErrorKind::Unsupported),
    ] {
        let error = parse_alignment::<Dna>(input.as_bytes(), AlignmentFormat::Fasta).unwrap_err();
        assert_eq!(error.kind(), kind, "input {input:?}: {error}");
    }
    let protein = parse_alignment::<Protein>(&b">a\nM*\n"[..], AlignmentFormat::Fasta).unwrap_err();
    assert_eq!(protein.kind(), IoErrorKind::InvalidSequence);
    let invalid =
        parse_alignment::<Dna>(&b">a\nAC\n>b\nAU\n"[..], AlignmentFormat::Fasta).unwrap_err();
    assert_eq!(invalid.record(), Some(2));
    assert_eq!(invalid.line(), Some(3));
}

#[test]
fn clustal_checks_blocks_counts_and_preserves_consensus() {
    let input = "CLUSTAL W (2.1) multiple sequence alignment\n\nalpha  AC- 2\nbeta   ACT 3\n       ** \n\nalpha  GT 4\nbeta   G- 4\n       * \n";
    let document = dna_document(input, AlignmentFormat::Clustal);
    assert_rows(&document, &[("alpha", "AC-GT"), ("beta", "ACTG-")]);
    assert_eq!(document.annotations().len(), 1);
    let annotation = &document.annotations()[0];
    assert_eq!(annotation.scope(), AnnotationScope::Column);
    assert_eq!(annotation.tag(), "consensus");
    assert_eq!(annotation.value(), "** * ");
    assert_eq!(
        document.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    assert_eq!(
        document
            .to_multiple_alignment_discarding_annotations()
            .unwrap()
            .width(),
        5
    );
    for replacement in [
        input.replace("beta   G- 4", "gamma  G- 4"),
        input.replace("alpha  GT 4", "alpha  GT 5"),
        input.replace("beta   G- 4", "beta   G 4"),
        input.replace("alpha  AC- 2", "alpha  AC- extra"),
        input.replace("\nbeta   G- 4", ""),
    ] {
        assert!(parse_alignment::<Dna>(replacement.as_bytes(), AlignmentFormat::Clustal).is_err());
    }
}

#[test]
fn all_four_phylip_profiles_have_independent_fixtures() {
    for (names, layout, input) in [
        (
            PhylipNameStyle::Relaxed,
            PhylipLayout::Sequential,
            "2 6\nalpha AC-\nGTA\nbeta ACTG\nTA\n",
        ),
        (
            PhylipNameStyle::Strict,
            PhylipLayout::Sequential,
            "2 6\nalpha     AC-\nGTA\nbeta      ACTG\nTA\n",
        ),
        (
            PhylipNameStyle::Relaxed,
            PhylipLayout::Interleaved,
            "2 6\nalpha AC-\nbeta ACT\n\nG TA\nG TA\n",
        ),
        (
            PhylipNameStyle::Strict,
            PhylipLayout::Interleaved,
            "2 6\nalpha     AC-\nbeta      ACT\n\n          G TA\n          G TA\n",
        ),
    ] {
        let document = dna_document(input, AlignmentFormat::Phylip { names, layout });
        assert_rows(&document, &[("alpha", "AC-GTA"), ("beta", "ACTGTA")]);
        assert_eq!(document.to_multiple_alignment().unwrap().width(), 6);
    }
    let strict = dna_document(
        "2 4\nH. sapiensAC-G\nmouse     ACTG\n",
        AlignmentFormat::Phylip {
            names: PhylipNameStyle::Strict,
            layout: PhylipLayout::Sequential,
        },
    );
    assert_eq!(strict.rows()[0].identifier(), "H. sapiens");
}

#[test]
fn phylip_dimension_truncation_options_and_trailing_data_are_checked() {
    let format = AlignmentFormat::Phylip {
        names: PhylipNameStyle::Relaxed,
        layout: PhylipLayout::Interleaved,
    };
    for input in [
        "2 4\na AC\nb ACT\nG\nG\n",
        "2 4\na AC\nb AC\nG\n",
        "2 4\na AC-GT\nb ACTGT\n",
        "2 4 W\na AC-G\nb ACTG\n",
        "2 4\na AC-G\nb ACTG\n1\n",
        "2 4\na AC.G\nb ACTG\n",
        "2 4\na AC3G\nb ACTG\n",
        "2 184467440737095516160\na A\nb A\n",
    ] {
        assert!(
            parse_alignment::<Dna>(input.as_bytes(), format).is_err(),
            "{input:?}"
        );
    }
    let empty = dna_document(
        "2 0\na\nb\n",
        AlignmentFormat::Phylip {
            names: PhylipNameStyle::Relaxed,
            layout: PhylipLayout::Sequential,
        },
    );
    assert_eq!(empty.width(), 0);
    let missing = dna_document("1 2\na A?\n", format);
    assert_eq!(
        missing.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
}

#[test]
fn stockholm_preserves_all_four_annotation_scopes_and_comments() {
    let input = "# STOCKHOLM 1.0\n#=GF ID example\n#=GF CC first paragraph\n#=GF CC second paragraph\n#=GS alpha DE full description\n# ordinary comment\nalpha AC.G\nbeta ACTG\n#=GR alpha SS HH.H\n#=GC SS_cons HH.H\n//\n";
    let document = dna_document(input, AlignmentFormat::Stockholm);
    assert_rows(&document, &[("alpha", "AC.G"), ("beta", "ACTG")]);
    let annotations = document.annotations();
    assert_eq!(annotations.len(), 7);
    assert_eq!(annotations[1].value(), "first paragraph");
    assert_eq!(annotations[2].value(), "second paragraph");
    assert_eq!(annotations[3].sequence(), Some("alpha"));
    assert_eq!(annotations[3].scope(), AnnotationScope::Sequence);
    assert_eq!(annotations[4].scope(), AnnotationScope::Comment);
    assert_eq!(annotations[5].scope(), AnnotationScope::Residue);
    assert_eq!(annotations[6].scope(), AnnotationScope::Column);
    assert_eq!(
        document.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    let converted = document
        .to_multiple_alignment_discarding_annotations()
        .unwrap();
    assert_eq!(converted.rows()[0].to_ascii_with_gaps(), "AC-G");
    assert_eq!(document.rows()[0].aligned(), "AC.G");
}

#[test]
fn stockholm_wrapping_and_annotation_dimensions_are_checked() {
    let wrapped = "# STOCKHOLM 1.0\na AC\nb AC\n#=GR a PP 99\n#=GC SS_cons HH\n\na .G\nb TG\n#=GR a PP .8\n#=GC SS_cons .H\n//\n";
    let document = dna_document(wrapped, AlignmentFormat::Stockholm);
    assert_rows(&document, &[("a", "AC.G"), ("b", "ACTG")]);
    assert_eq!(document.annotations()[0].value(), "99.8");
    assert_eq!(document.annotations()[1].value(), "HH.H");
    for invalid in [
        wrapped.replace("//\n", ""),
        wrapped.replace("//\n", "//\n# STOCKHOLM 1.0\n"),
        wrapped.replace("#=GR a PP .8", "#=GR a PP 8"),
        wrapped.replace("#=GR a PP 99", "#=GR absent PP 99"),
        wrapped.replace("#=GC SS_cons HH", "#=XX SS_cons HH"),
        wrapped.replace("#=GC SS_cons HH", "#=GC SS_cons H H"),
    ] {
        assert!(
            parse_alignment::<Dna>(invalid.as_bytes(), AlignmentFormat::Stockholm).is_err(),
            "{invalid:?}"
        );
    }
}

const NEXUS: &str = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA MISSING=? GAP=-;\nMATRIX\n'first taxon' AC-G\nsecond_taxon ACTG\n;\nEND;\n";

#[test]
fn nexus_quotes_underscores_datatypes_and_custom_glyphs_are_explicit() {
    let document = dna_document(NEXUS, AlignmentFormat::Nexus);
    assert_rows(
        &document,
        &[("first taxon", "AC-G"), ("second taxon", "ACTG")],
    );
    assert_eq!(document.to_multiple_alignment().unwrap().width(), 4);
    let custom = NEXUS
        .replace("GAP=-", "GAP=.")
        .replace("MISSING=?", "MISSING=@")
        .replace("AC-G", "AC.G")
        .replace("ACTG", "AC@G")
        .replace("'first taxon'", "'a;[name]''s underscore_'");
    let document = dna_document(&custom, AlignmentFormat::Nexus);
    assert_rows(
        &document,
        &[("a;[name]'s underscore_", "AC.G"), ("second taxon", "AC@G")],
    );
    assert_eq!(document.gap_symbols(), b".");
    assert_eq!(document.missing_symbol(), Some(b'@'));
    assert_eq!(
        document.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    // Use an independently authored RNA fixture, not DNA coercion at import.
    let rna =
        "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=4; FORMAT DATATYPE=RNA; MATRIX a ACUG; END;\n";
    assert_eq!(
        parse_alignment::<Rna>(rna.as_bytes(), AlignmentFormat::Nexus)
            .unwrap()
            .rows()[0]
            .aligned(),
        "ACUG"
    );
    assert!(parse_alignment::<Dna>(rna.as_bytes(), AlignmentFormat::Nexus).is_err());
    let protein = "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=4; FORMAT DATATYPE=PROTEIN; MATRIX p MKTW; END;\n";
    assert!(parse_alignment::<Protein>(protein.as_bytes(), AlignmentFormat::Nexus).is_ok());
}

#[test]
fn nexus_comments_are_retained_and_unsupported_constructs_are_not_skipped() {
    let document = dna_document(
        &NEXUS.replace("AC-G", "AC[nested [comment]]-G"),
        AlignmentFormat::Nexus,
    );
    assert_eq!(document.annotations()[0].value(), "nested [comment]");
    assert_eq!(
        document.to_multiple_alignment().unwrap_err().kind(),
        IoErrorKind::DataLoss
    );
    for invalid in [
        NEXUS.replace("END;", "END; BEGIN TREES; TREE t=(a,b); END;"),
        NEXUS.replace("DATATYPE=DNA", "DATATYPE=DNA INTERLEAVE=YES"),
        NEXUS.replace("DATATYPE=DNA", "DATATYPE=DNA MATCHCHAR=."),
        NEXUS.replace("NTAX=2", "NTAX=3"),
        NEXUS.replace("NCHAR=4", "NCHAR=5"),
        NEXUS.replace("NCHAR=4", "NCHAR=4 NCHAR=4"),
        NEXUS.replace("MISSING=?", "MISSING=-"),
        NEXUS.replace("GAP=-", "GAP=A"),
        NEXUS.replace("END;", "END"),
        NEXUS.replace("END;", "[unterminated END;"),
        NEXUS.replace("'first taxon'", "'unterminated taxon"),
    ] {
        assert!(
            parse_alignment::<Dna>(invalid.as_bytes(), AlignmentFormat::Nexus).is_err(),
            "{invalid:?}"
        );
    }
}

#[test]
fn nexus_diagnostics_identify_original_source_lines() {
    let input = format!("\n\n{}", NEXUS.replace("ACTG", "ACT"));
    let error = parse_alignment::<Dna>(input.as_bytes(), AlignmentFormat::Nexus).unwrap_err();
    assert_eq!(error.record(), Some(2));
    assert_eq!(error.line(), Some(9));
    let document = dna_document(NEXUS, AlignmentFormat::Nexus);
    assert_eq!(document.rows()[0].line(), 6);
    assert_eq!(document.rows()[1].line(), 7);
}

#[test]
fn all_existing_alignment_writers_round_trip_supported_semantics() -> TestResult {
    let alignment = dna_document(ALIGNED_FASTA, AlignmentFormat::Fasta).to_multiple_alignment()?;
    for format in [
        AlignmentFormat::Fasta,
        AlignmentFormat::Clustal,
        AlignmentFormat::Stockholm,
        AlignmentFormat::Nexus,
    ] {
        let mut bytes = Vec::new();
        match format {
            AlignmentFormat::Fasta => write_alignment_fasta_to(&mut bytes, &alignment)?,
            AlignmentFormat::Clustal => write_clustal_to(&mut bytes, &alignment)?,
            AlignmentFormat::Stockholm => write_stockholm_to(&mut bytes, &alignment)?,
            AlignmentFormat::Nexus => write_nexus_to(&mut bytes, &alignment)?,
            AlignmentFormat::Phylip { .. } => unreachable!(),
        }
        let document = parse_alignment::<Dna>(&bytes[..], format)?;
        assert_eq!(
            document.to_multiple_alignment_discarding_annotations()?,
            alignment
        );
    }
    for names in [PhylipNameStyle::Strict, PhylipNameStyle::Relaxed] {
        for layout in [PhylipLayout::Sequential, PhylipLayout::Interleaved] {
            let mut bytes = Vec::new();
            write_phylip_with_options_to(
                &mut bytes,
                &alignment,
                PhylipOptions::new()
                    .with_layout(layout)
                    .with_name_style(names)
                    .with_block_width(2)?,
            )?;
            let document =
                parse_alignment::<Dna>(&bytes[..], AlignmentFormat::Phylip { names, layout })?;
            assert_eq!(document.to_multiple_alignment()?, alignment);
        }
    }
    Ok(())
}

#[test]
fn io_and_non_utf8_errors_are_returned_without_panics() {
    struct Failing;
    impl Read for Failing {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("injected read failure"))
        }
    }
    let error =
        parse_alignment::<Dna>(BufReader::new(Failing), AlignmentFormat::Fasta).unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::Io);
    assert!(parse_alignment::<Dna>(&b">name\n\xff\n"[..], AlignmentFormat::Fasta).is_err());
    let strict = AlignmentFormat::Phylip {
        names: PhylipNameStyle::Strict,
        layout: PhylipLayout::Sequential,
    };
    assert!(parse_alignment::<Dna>("1 1\n123456789éA\n".as_bytes(), strict).is_err());
}

#[test]
fn path_reader_detects_gzip_and_checks_integrity() -> TestResult {
    use std::io::Write as _;

    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(ALIGNED_FASTA.as_bytes())?;
    let gzip = encoder.finish()?;
    let mut corrupt = gzip.clone();
    let trailer = corrupt.len() - 8;
    corrupt[trailer] ^= 1;
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)?
        .as_nanos();
    for (index, (bytes, valid)) in [
        (ALIGNED_FASTA.as_bytes().to_vec(), true),
        (gzip, true),
        (corrupt, false),
    ]
    .into_iter()
    .enumerate()
    {
        let path = std::env::temp_dir().join(format!(
            "seq-bio-alignment-{}-{stamp}-{index}.data",
            std::process::id()
        ));
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)?;
        file.write_all(&bytes)?;
        drop(file);
        let result = read_alignment::<Dna>(&path, AlignmentFormat::Fasta);
        std::fs::remove_file(&path)?;
        if valid {
            assert_rows(&result?, &[("reference", "AC-G"), ("sample", "ACTG")]);
        } else {
            assert_eq!(result.unwrap_err().kind(), IoErrorKind::Io);
        }
        assert_eq!(
            read_alignment::<Dna>(&path, AlignmentFormat::Fasta)
                .unwrap_err()
                .kind(),
            IoErrorKind::Io
        );
    }
    Ok(())
}

proptest::proptest! {
    #![proptest_config(proptest::test_runner::Config::with_cases(128))]
    #[test]
    fn arbitrary_utf8_is_handled_without_panics(characters in proptest::collection::vec(proptest::char::any(), 0..512)) {
        let input: String = characters.into_iter().collect();
        for format in [
            AlignmentFormat::Fasta,
            AlignmentFormat::Clustal,
            AlignmentFormat::Phylip { names: PhylipNameStyle::Strict, layout: PhylipLayout::Sequential },
            AlignmentFormat::Phylip { names: PhylipNameStyle::Strict, layout: PhylipLayout::Interleaved },
            AlignmentFormat::Phylip { names: PhylipNameStyle::Relaxed, layout: PhylipLayout::Sequential },
            AlignmentFormat::Phylip { names: PhylipNameStyle::Relaxed, layout: PhylipLayout::Interleaved },
            AlignmentFormat::Stockholm,
            AlignmentFormat::Nexus,
        ] {
            let _ = parse_alignment::<Dna>(input.as_bytes(), format);
            let framed = match format {
                AlignmentFormat::Fasta => format!(">row\n{input}\n"),
                AlignmentFormat::Clustal => format!("CLUSTAL W\n\nrow {input}\n"),
                AlignmentFormat::Phylip { .. } => format!("1 512\nrow       {input}\n"),
                AlignmentFormat::Stockholm => format!("# STOCKHOLM 1.0\nrow {input}\n//\n"),
                AlignmentFormat::Nexus => format!("#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=512; FORMAT DATATYPE=DNA; MATRIX row {input}; END;\n"),
            };
            let _ = parse_alignment::<Dna>(framed.as_bytes(), format);
        }
    }
}
