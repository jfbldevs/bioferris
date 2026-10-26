//! Independent biological, extraction and qualifier checks for annotated CDS translation.

use std::error::Error;

use bioferris::io::IoErrorKind;
use bioferris::io::annotated::{Feature, Location, Qualifier, parse_embl, parse_genbank};
use bioferris::prelude::*;
use proptest::prelude::*;
use sha2::{Digest, Sha256};

const CODES: [GeneticCode; 3] = [
    GeneticCode::Standard,
    GeneticCode::VertebrateMitochondrial,
    GeneticCode::BacterialArchaealPlastid,
];

// Biological initiation and termination sets from NCBI's genetic-code page,
// revision 2024-09-23, consulted 2026-10-17. The explicit codon sets are
// independent of the production engine's indexed tables and ambiguity masks.
// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
fn is_start(codon: &str, code: GeneticCode) -> bool {
    let starts: &[&str] = match code {
        GeneticCode::Standard => &["TTG", "CTG", "ATG"],
        GeneticCode::VertebrateMitochondrial => &["ATT", "ATC", "ATA", "ATG", "GTG"],
        GeneticCode::BacterialArchaealPlastid => &["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
        // The other tables are checked against NCBI gc.prt in ncbi_genetic_codes.rs.
        other => panic!("no independent initiation set for table {}", other.id()),
    };
    starts.contains(&codon)
}

fn is_stop(codon: &str, code: GeneticCode) -> bool {
    let stops: &[&str] = match code {
        GeneticCode::Standard | GeneticCode::BacterialArchaealPlastid => &["TAA", "TAG", "TGA"],
        GeneticCode::VertebrateMitochondrial => &["TAA", "TAG", "AGA", "AGG"],
        other => panic!("no independent stop set for table {}", other.id()),
    };
    stops.contains(&codon)
}

fn feature(location: &str, qualifiers: &[(&str, Option<&str>)]) -> Feature {
    Feature {
        key: "CDS".to_owned(),
        location: Location::parse(location).unwrap(),
        qualifiers: qualifiers
            .iter()
            .map(|(name, value)| Qualifier {
                name: (*name).to_owned(),
                value: value.map(str::to_owned),
                raw: String::new(),
            })
            .collect(),
        raw_lines: Vec::new(),
    }
}

fn complete_feature(sequence: &str, qualifiers: &[(&str, Option<&str>)]) -> Feature {
    feature(&format!("1..{}", sequence.len()), qualifiers)
}

fn options(code: GeneticCode, policy: CdsReferencePolicy) -> CdsTranslationOptions {
    CdsTranslationOptions {
        default_genetic_code: code,
        reference_policy: policy,
    }
}

fn assert_invalid_qualifier(error: CdsTranslationError, expected: &str) {
    match error {
        CdsTranslationError::InvalidQualifier { name, reason } => {
            assert_eq!(name, expected);
            assert!(!reason.is_empty());
        }
        other => panic!("expected invalid /{expected}, got {other:?}"),
    }
}

// This deliberately does not call BioFerris's reverse-complement operation.
fn reverse_complement(input: &str) -> String {
    input
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!("property generator uses concrete DNA"),
        })
        .collect()
}

#[test]
fn default_result_retains_complete_cds_metadata_without_a_reference() {
    let sequence = dna("ATGGCTTAA").unwrap();
    let cds = complete_feature("ATGGCTTAA", &[]);
    let result = cds.translate_cds(&sequence).unwrap();
    assert_eq!(result.peptide(), "MA");
    assert_eq!(result.translation().as_str(), "MA");
    assert_eq!(result.translation().to_protein().unwrap().to_string(), "MA");
    assert_eq!(result.feature_length(), 9);
    assert_eq!(result.codon_start(), 1);
    assert_eq!(result.genetic_code(), GeneticCode::Standard);
    assert_eq!(result.translation().reading_frame(), ReadingFrame::Zero);
    assert_eq!(result.translation().skipped_leading_bases(), 0);
    assert_eq!(result.translation().ignored_trailing_bases(), 0);
    assert_eq!(result.reference_status(), CdsReferenceStatus::Absent);
    assert_eq!(result.reference_translation(), None);
    assert_eq!(
        CdsTranslationOptions::default(),
        options(GeneticCode::Standard, CdsReferencePolicy::CheckIfPresent)
    );

    // These names must also be available at the annotated module boundary.
    let _: bioferris::io::annotated::CdsTranslationOptions = CdsTranslationOptions::default();
    let _: bioferris::CdsTranslation = result;
}

#[test]
fn all_concrete_initiators_and_candidate_stops_follow_the_selected_table() {
    let mut comparisons = 0;
    for code in CODES {
        for first in [b'A', b'C', b'G', b'T'] {
            for second in [b'A', b'C', b'G', b'T'] {
                for third in [b'A', b'C', b'G', b'T'] {
                    let start = String::from_utf8(vec![first, second, third]).unwrap();
                    for stop in ["TAA", "TAG", "TGA", "AGA", "AGG"] {
                        let input = format!("{start}GCT{stop}");
                        let cds = complete_feature(&input, &[]);
                        let config = options(code, CdsReferencePolicy::CheckIfPresent);
                        for result in [
                            cds.translate_cds_with(&dna(&input).unwrap(), config),
                            cds.translate_cds_with(&rna(input.replace('T', "U")).unwrap(), config),
                        ] {
                            if is_start(&start, code) && is_stop(stop, code) {
                                let result = result.unwrap_or_else(|error| {
                                    panic!("table {}, {input}: {error}", code.id())
                                });
                                assert_eq!(result.peptide(), "MA", "table {}, {input}", code.id());
                                assert_eq!(result.genetic_code(), code);
                                assert_eq!(result.feature_length(), 9);
                            } else {
                                assert!(result.is_err(), "table {}, {input}", code.id());
                            }
                            comparisons += 1;
                        }
                    }
                }
            }
        }
    }
    assert_eq!(comparisons, 1_920);
}

#[test]
fn table_qualifier_overrides_the_fallback_and_preserves_internal_assignments() {
    let input = "GTGATATGATAA";
    let sequence = dna(input).unwrap();
    let mitochondrial = complete_feature(input, &[("transl_table", Some("2"))]);
    let result = mitochondrial.translate_cds(&sequence).unwrap();
    assert_eq!(result.peptide(), "MMW");
    assert_eq!(result.genetic_code(), GeneticCode::VertebrateMitochondrial);
    assert_eq!(
        result.translation().genetic_code(),
        GeneticCode::VertebrateMitochondrial
    );

    let bacterial = complete_feature("GTGATATAA", &[("transl_table", Some("11"))]);
    let result = bacterial
        .translate_cds_with(
            &dna("GTGATATAA").unwrap(),
            options(
                GeneticCode::VertebrateMitochondrial,
                CdsReferencePolicy::CheckIfPresent,
            ),
        )
        .unwrap();
    assert_eq!(result.peptide(), "MI");
    assert_eq!(result.genetic_code(), GeneticCode::BacterialArchaealPlastid);

    let fallback = complete_feature(input, &[])
        .translate_cds_with(
            &sequence,
            options(
                GeneticCode::VertebrateMitochondrial,
                CdsReferencePolicy::CheckIfPresent,
            ),
        )
        .unwrap();
    assert_eq!(fallback.peptide(), "MMW");
}

#[test]
fn exact_locations_and_written_join_order_produce_the_known_coding_sequence() {
    for (parent, location) in [
        ("CCCATGAAATAAGGG", "4..12"),
        ("TAACCCATGAAAGGG", "join(7..12,1..3)"),
        ("TTACCCTTTCAT", "complement(join(1..3,7..12))"),
        ("TTACCCTTTCAT", "join(complement(7..12),complement(1..3))"),
        ("ATGAAATAA", "join(1,2,3..9)"),
        ("ATGAAATAA", "complement(complement(1..9))"),
    ] {
        let source = dna(parent).unwrap();
        let cds = feature(location, &[]);
        assert_eq!(
            cds.translate_cds(&source).unwrap().peptide(),
            "MK",
            "{location}"
        );
        assert_eq!(source.to_string(), parent);
    }
}

#[test]
fn codon_start_is_applied_after_the_entire_join_and_reverse_complement() {
    for (parent, location, codon_start, offset) in [
        ("CATGAAATAA", "1..10", "2", ReadingFrame::One),
        ("CCATGAAATAA", "1..11", "3", ReadingFrame::Two),
        (
            "TTACCCTTTCATG",
            "complement(join(1..3,7..13))",
            "2",
            ReadingFrame::One,
        ),
        (
            "TTACCCTTTCATGG",
            "complement(join(1..3,7..14))",
            "3",
            ReadingFrame::Two,
        ),
    ] {
        let cds = feature(location, &[("codon_start", Some(codon_start))]);
        let result = cds.translate_cds(&dna(parent).unwrap()).unwrap();
        assert_eq!(result.peptide(), "MK", "{location}");
        assert_eq!(result.codon_start().to_string(), codon_start);
        assert_eq!(result.translation().reading_frame(), offset);
        assert_eq!(
            result.translation().skipped_leading_bases(),
            offset.offset()
        );
        assert_eq!(result.feature_length(), 9 + offset.offset());
    }
}

#[test]
fn rna_complemented_join_and_alternative_initiation_are_supported() {
    let source = rna("UUACCCUUUCACG").unwrap();
    let cds = feature(
        "complement(join(1..3,7..13))",
        &[("codon_start", Some("2")), ("transl_table", Some("11"))],
    );
    let result = cds.translate_cds(&source).unwrap();
    assert_eq!(result.peptide(), "MK");
    assert_eq!(source.to_string(), "UUACCCUUUCACG");
    assert_eq!(result.genetic_code(), GeneticCode::BacterialArchaealPlastid);
}

#[test]
fn ambiguous_cds_cases_have_independently_known_residues_and_boundary_rules() {
    for (input, code, expected) in [
        ("HTGG CNTAA", GeneticCode::Standard, "MA"),
        ("ATGRAYTAA", GeneticCode::Standard, "MB"),
        ("ATGSARTAR", GeneticCode::Standard, "MZ"),
        ("ATGWTATAR", GeneticCode::Standard, "MJ"),
        ("ATGTANTAA", GeneticCode::Standard, "MX"),
        ("ATNTGATAG", GeneticCode::VertebrateMitochondrial, "MW"),
        ("RTGTGRTAA", GeneticCode::VertebrateMitochondrial, "MW"),
        ("NTGG CNTGA", GeneticCode::BacterialArchaealPlastid, "MA"),
    ] {
        let input = input.replace(' ', "");
        let cds = complete_feature(&input, &[("translation", Some(expected))]);
        let result = cds
            .translate_cds_with(
                &dna(&input).unwrap(),
                options(code, CdsReferencePolicy::RequireMatch),
            )
            .unwrap();
        assert_eq!(result.peptide(), expected, "{input}");
        assert_eq!(result.reference_status(), CdsReferenceStatus::Matched);
    }
    for input in ["NTGGCTTAA", "ATGGCTTAN", "ATGTARGCTTAA"] {
        assert!(
            complete_feature(input, &[])
                .translate_cds(&dna(input).unwrap())
                .is_err(),
            "{input}"
        );
    }
}

#[test]
fn complete_cds_errors_remain_typed_and_report_extracted_feature_offsets() {
    for (input, expected) in [
        ("ATG", TranslationError::CdsTooShort { length: 3 }),
        (
            "ATGAAATA",
            TranslationError::PartialCodon {
                length: 8,
                remainder: 2,
            },
        ),
        (
            "AAAGCTTAA",
            TranslationError::InvalidStartCodon { codon: *b"AAA" },
        ),
        (
            "ATGGCTAAA",
            TranslationError::MissingTerminalStop { codon: *b"AAA" },
        ),
    ] {
        let error = complete_feature(input, &[])
            .translate_cds(&dna(input).unwrap())
            .unwrap_err();
        match error {
            CdsTranslationError::Translation(actual) => assert_eq!(actual, expected),
            other => panic!("expected a translation error for {input}, got {other:?}"),
        }
    }
    let coding_oriented = "CATGTAAGCTTAA";
    let parent = format!("CCC{}GGG", reverse_complement(coding_oriented));
    let cds = feature("complement(4..16)", &[("codon_start", Some("2"))]);
    match cds.translate_cds(&dna(&parent).unwrap()).unwrap_err() {
        CdsTranslationError::Translation(TranslationError::InternalStop { offset }) => {
            assert_eq!(offset, 4);
        }
        other => panic!("expected extracted feature offset 4, got {other:?}"),
    }
}

#[test]
fn only_the_exact_case_sensitive_cds_key_is_translated() {
    let sequence = dna("ATGGCTTAA").unwrap();
    for key in ["cds", "Cds", "gene", "CDS ", ""] {
        let mut cds = complete_feature("ATGGCTTAA", &[]);
        cds.key = key.to_owned();
        match cds.translate_cds(&sequence).unwrap_err() {
            CdsTranslationError::NotCds { key: actual } => assert_eq!(actual, key),
            other => panic!("unexpected error for {key:?}: {other:?}"),
        }
    }
}

#[test]
fn numeric_qualifiers_accept_outer_ascii_whitespace_and_decimal_zero_padding() {
    let cds = complete_feature(
        "CCGTGAAATAA",
        &[
            ("codon_start", Some(" \t003\r\n")),
            ("transl_table", Some("\t011 \r")),
        ],
    );
    let result = cds.translate_cds(&dna("CCGTGAAATAA").unwrap()).unwrap();
    assert_eq!(result.codon_start(), 3);
    assert_eq!(result.genetic_code(), GeneticCode::BacterialArchaealPlastid);
    assert_eq!(result.peptide(), "MK");
}

#[test]
fn numeric_qualifiers_reject_flags_nondecimal_internal_space_and_overflow() {
    let source = dna("ATGGCTTAA").unwrap();
    for name in ["codon_start", "transl_table"] {
        for value in [
            None,
            Some(""),
            Some(" \t"),
            Some("+1"),
            Some("-1"),
            Some("1.0"),
            Some("1 1"),
            Some("1\n1"),
            Some("１"),
            Some("1\u{a0}"),
            Some("65536"),
            Some("999999999999999999999999999999999999"),
        ] {
            let error = complete_feature("ATGGCTTAA", &[(name, value)])
                .translate_cds(&source)
                .unwrap_err();
            assert_invalid_qualifier(error, name);
        }
    }
    for value in ["0", "4", "255", "256"] {
        let error = complete_feature("ATGGCTTAA", &[("codon_start", Some(value))])
            .translate_cds(&source)
            .unwrap_err();
        assert_invalid_qualifier(error, "codon_start");
    }
}

#[test]
fn unsupported_numeric_tables_are_not_replaced_with_a_fallback() {
    let source = dna("ATGGCTTAA").unwrap();
    for id in [0, 7, 8, 17, 20, 34, 65535] {
        let value = id.to_string();
        let error = complete_feature("ATGGCTTAA", &[("transl_table", Some(&value))])
            .translate_cds_with(
                &source,
                options(
                    GeneticCode::BacterialArchaealPlastid,
                    CdsReferencePolicy::Ignore,
                ),
            )
            .unwrap_err();
        match error {
            CdsTranslationError::Translation(TranslationError::UnsupportedGeneticCode {
                id: actual,
            }) => assert_eq!(actual, id),
            other => panic!("expected unsupported table {id}, got {other:?}"),
        }
    }
}

#[test]
fn duplicate_interpretation_qualifiers_fail_even_when_values_agree() {
    let source = dna("ATGGCTTAA").unwrap();
    for (name, value) in [
        ("codon_start", "1"),
        ("transl_table", "1"),
        ("translation", "MA"),
    ] {
        let cds = complete_feature("ATGGCTTAA", &[(name, Some(value)), (name, Some(value))]);
        match cds.translate_cds(&source).unwrap_err() {
            CdsTranslationError::DuplicateQualifier { name: actual } => assert_eq!(actual, name),
            other => panic!("expected duplicate /{name}, got {other:?}"),
        }
    }
}

#[test]
fn ordinary_metadata_and_case_distinct_qualifier_names_do_not_change_translation() {
    let cds = complete_feature(
        "ATGGCTTAA",
        &[
            ("gene", Some("demo")),
            ("note", Some("first")),
            ("note", Some("second")),
            ("db_xref", Some("DB:1")),
            ("db_xref", Some("DB:2")),
            ("Codon_start", Some("not interpreted")),
            ("TRANSL_TABLE", Some("999")),
            ("Translation", Some("wrong")),
        ],
    );
    let original = cds.clone();
    let sequence = dna("ATGGCTTAA").unwrap();
    assert_eq!(cds.translate_cds(&sequence).unwrap().peptide(), "MA");
    assert_eq!(cds, original);
}

#[test]
fn all_exception_and_pseudogene_qualifiers_fail_under_every_reference_policy() {
    let source = dna("ATGGCTTAA").unwrap();
    for name in [
        "pseudo",
        "pseudogene",
        "transl_except",
        "exception",
        "trans_splicing",
    ] {
        for value in [None, Some(""), Some("retained exception text")] {
            let cds = complete_feature("ATGGCTTAA", &[(name, value)]);
            for policy in [
                CdsReferencePolicy::CheckIfPresent,
                CdsReferencePolicy::RequireMatch,
                CdsReferencePolicy::Ignore,
            ] {
                match cds
                    .translate_cds_with(&source, options(GeneticCode::Standard, policy))
                    .unwrap_err()
                {
                    CdsTranslationError::UnsupportedQualifier { name: actual } => {
                        assert_eq!(actual, name);
                    }
                    other => panic!("expected refusal of /{name} under {policy:?}, got {other:?}"),
                }
            }
        }
    }
}

#[test]
fn wrapped_reference_normalization_is_exact_and_returns_the_checked_reference() {
    let source = dna("ATGGCTTAA").unwrap();
    let cds = complete_feature("ATGGCTTAA", &[("translation", Some(" \tM\n A\r\n"))]);
    for policy in [
        CdsReferencePolicy::CheckIfPresent,
        CdsReferencePolicy::RequireMatch,
    ] {
        let result = cds
            .translate_cds_with(&source, options(GeneticCode::Standard, policy))
            .unwrap();
        assert_eq!(result.peptide(), "MA");
        assert_eq!(result.reference_translation(), Some("MA"));
        assert_eq!(result.reference_status(), CdsReferenceStatus::Matched);
    }
}

#[test]
fn invalid_reference_residues_case_stops_and_non_ascii_space_are_refused() {
    let source = dna("ATGGCTTAA").unwrap();
    for value in [
        None,
        Some(""),
        Some(" \t\r\n"),
        Some("ma"),
        Some("Ma"),
        Some("MA*"),
        Some("M-A"),
        Some("M?"),
        Some("M\u{a0}A"),
        Some("MÄ"),
    ] {
        let error = complete_feature("ATGGCTTAA", &[("translation", value)])
            .translate_cds(&source)
            .unwrap_err();
        assert_invalid_qualifier(error, "translation");
    }
}

#[test]
fn reference_mismatch_reports_exact_strings_and_the_first_differing_position() {
    let source = dna("ATGGCTTAA").unwrap();
    for (reference, difference) in [("M", 1), ("MAC", 2), ("MC", 1), ("AA", 0), ("MX", 1)] {
        let error = complete_feature("ATGGCTTAA", &[("translation", Some(reference))])
            .translate_cds(&source)
            .unwrap_err();
        match error {
            CdsTranslationError::ReferenceMismatch {
                reference: actual,
                computed,
                first_difference,
            } => {
                assert_eq!(actual, reference);
                assert_eq!(computed, "MA");
                assert_eq!(first_difference, difference);
            }
            other => panic!("expected reference mismatch for {reference}, got {other:?}"),
        }
    }
}

#[test]
fn requiring_a_reference_and_ignoring_a_reference_are_distinct_policies() {
    let source = dna("ATGGCTTAA").unwrap();
    assert!(matches!(
        complete_feature("ATGGCTTAA", &[]).translate_cds_with(
            &source,
            options(GeneticCode::Standard, CdsReferencePolicy::RequireMatch)
        ),
        Err(CdsTranslationError::MissingReference)
    ));
    let cds = complete_feature(
        "ATGGCTTAA",
        &[
            ("translation", None),
            ("translation", Some("not a protein*")),
        ],
    );
    let result = cds
        .translate_cds_with(
            &source,
            options(GeneticCode::Standard, CdsReferencePolicy::Ignore),
        )
        .unwrap();
    assert_eq!(result.peptide(), "MA");
    assert_eq!(result.reference_status(), CdsReferenceStatus::Ignored);
    assert_eq!(result.reference_translation(), None);

    let missing = complete_feature("ATGGCTTAA", &[])
        .translate_cds_with(
            &source,
            options(GeneticCode::Standard, CdsReferencePolicy::Ignore),
        )
        .unwrap();
    assert_eq!(missing.reference_status(), CdsReferenceStatus::Ignored);
    assert_eq!(missing.reference_translation(), None);
}

#[test]
fn unsupported_and_out_of_bounds_locations_keep_their_extraction_errors() {
    let source = dna("ATGGCTTAA").unwrap();
    for location in [
        "<1..9",
        "1..>9",
        "OTHER.1:1..9",
        "order(1..3,4..9)",
        "3^4",
        "join(1..3,OTHER.1:4..9)",
    ] {
        let error = feature(location, &[]).translate_cds(&source).unwrap_err();
        assert!(error.source().is_some());
        match error {
            CdsTranslationError::Extraction(inner) => {
                assert_eq!(inner.kind(), IoErrorKind::Unsupported, "{location}");
                assert_eq!(inner.format(), "insdc-location");
            }
            other => panic!("expected unsupported extraction for {location}, got {other:?}"),
        }
    }
    match feature("1..10", &[]).translate_cds(&source).unwrap_err() {
        CdsTranslationError::Extraction(inner) => {
            assert_eq!(inner.kind(), IoErrorKind::LengthMismatch);
        }
        other => panic!("expected a checked bounds error, got {other:?}"),
    }
    assert_eq!(source.to_string(), "ATGGCTTAA");
}

#[test]
fn error_categories_and_underlying_sources_remain_inspectable() {
    let source = dna("ATGGCTTAA").unwrap();
    for (cds, config, code, has_source) in [
        (
            {
                let mut cds = complete_feature("ATGGCTTAA", &[]);
                cds.key = "gene".to_owned();
                cds
            },
            CdsTranslationOptions::default(),
            "invalid_cds_feature",
            false,
        ),
        (
            complete_feature("ATGGCTTAA", &[("pseudo", None)]),
            CdsTranslationOptions::default(),
            "unsupported_cds_qualifier",
            false,
        ),
        (
            complete_feature("ATGGCTTAA", &[("codon_start", Some("0"))]),
            CdsTranslationOptions::default(),
            "invalid_cds_qualifier",
            false,
        ),
        (
            complete_feature("ATGGCTTAA", &[]),
            options(GeneticCode::Standard, CdsReferencePolicy::RequireMatch),
            "missing_cds_translation",
            false,
        ),
        (
            complete_feature("ATGGCTTAA", &[("translation", Some("MC"))]),
            CdsTranslationOptions::default(),
            "cds_translation_mismatch",
            false,
        ),
        (
            complete_feature("ATGGCTTAA", &[("transl_table", Some("7"))]),
            CdsTranslationOptions::default(),
            "unsupported_genetic_code",
            true,
        ),
        (
            feature("<1..9", &[]),
            CdsTranslationOptions::default(),
            "unsupported_format_feature",
            true,
        ),
        (
            feature("1..3", &[]),
            CdsTranslationOptions::default(),
            "cds_translation_error",
            true,
        ),
    ] {
        let error = cds.translate_cds_with(&source, config).unwrap_err();
        assert_eq!(error.code(), code);
        assert_eq!(error.source().is_some(), has_source);
        assert!(!error.to_string().is_empty());
    }
    let short = feature("1", &[("codon_start", Some("3"))]);
    assert!(matches!(
        short.translate_cds(&dna("A").unwrap()),
        Err(CdsTranslationError::Translation(
            TranslationError::FrameOutOfBounds {
                frame: 2,
                length: 1
            }
        ))
    ));
}

#[test]
fn actual_genbank_and_embl_models_translate_without_mutating_source_records() {
    let genbank = parse_genbank::<Dna>(include_str!("fixtures/annotated-example.gb")).unwrap();
    let embl = parse_embl::<Dna>(include_str!("fixtures/annotated-example.embl")).unwrap();
    for record in genbank.iter().chain(&embl) {
        let before = record.clone();
        let cds = &record.features()[1];
        let result = cds
            .translate_cds_with(
                record.record().sequence(),
                options(
                    GeneticCode::VertebrateMitochondrial,
                    CdsReferencePolicy::RequireMatch,
                ),
            )
            .unwrap();
        assert_eq!(result.peptide(), "MA");
        assert_eq!(result.genetic_code(), GeneticCode::Standard);
        assert_eq!(result.reference_translation(), Some("MA"));
        assert_eq!(result.reference_status(), CdsReferenceStatus::Matched);
        assert_eq!(record, &before);
    }
}

#[test]
fn public_partial_p53_annotations_remain_explicitly_untranslatable() {
    let genbank = parse_genbank::<Dna>(include_str!("fixtures/annotated-dq067560.1.gb")).unwrap();
    let embl = parse_embl::<Dna>(include_str!("fixtures/annotated-dq067560.1.embl")).unwrap();
    for record in genbank.iter().chain(&embl) {
        let original = record.clone();
        let cds = record
            .features()
            .iter()
            .find(|feature| feature.key == "CDS")
            .unwrap();
        match cds.translate_cds(record.record().sequence()).unwrap_err() {
            CdsTranslationError::Extraction(inner) => {
                assert_eq!(inner.kind(), IoErrorKind::Unsupported);
            }
            other => panic!("fuzzy source annotation must not be guessed: {other:?}"),
        }
        assert_eq!(record, &original);
    }
}

#[test]
fn public_phix174_cds_products_match_all_eleven_published_proteins() {
    let snapshot = include_str!("fixtures/annotated-nc_001422.1.gb");
    assert_eq!(
        format!("{:x}", Sha256::digest(snapshot.as_bytes())),
        "0a8099e207633f8c3ba998bfe2a198518e1de356c7956691101b2512324b807c"
    );
    let records = parse_genbank::<Dna>(snapshot).unwrap();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    let reference_bytes = include_bytes!("fixtures/cds-biopython-1.88.json");
    assert_eq!(
        format!("{:x}", Sha256::digest(reference_bytes)),
        "0c493b52a3b5a6f809fd982b6df8f30c3cf3b31e22925aa39111c3aeaaae5420"
    );
    let differential: serde_json::Value = serde_json::from_slice(reference_bytes).unwrap();
    assert_eq!(differential["reference"], "Biopython");
    assert_eq!(differential["version"], "1.88");
    assert_eq!(differential["interpreter"], "3.12.3");
    let differential_rows = differential["cds"].as_array().unwrap();
    let original = record.clone();
    assert_eq!(record.record().id(), "NC_001422.1");
    assert_eq!(record.record().sequence().len(), 5_386);
    // Lengths follow the independent one-based inclusive source locations,
    // subtracting the separately encoded final stop. The first three genes
    // explicitly cross the genome's origin through ordered joins.
    let expected = [
        ("NP_040703.1", "join(3981..5386,1..136)", 513),
        ("NP_040704.1", "join(4497..5386,1..136)", 341),
        ("NP_040705.1", "join(5075..5386,1..51)", 120),
        ("NP_040706.1", "51..221", 56),
        ("NP_040707.1", "133..393", 86),
        ("NP_040708.1", "390..848", 152),
        ("NP_040709.1", "568..843", 91),
        ("NP_040710.1", "848..964", 38),
        ("NP_040711.1", "1001..2284", 427),
        ("NP_040712.1", "2395..2922", 175),
        ("NP_040713.1", "2931..3917", 328),
    ];
    let coding: Vec<_> = record
        .features()
        .iter()
        .filter(|cds| cds.key == "CDS")
        .collect();
    assert_eq!(coding.len(), expected.len());
    assert_eq!(differential_rows.len(), expected.len());
    for ((cds, (protein_id, location, length)), reference) in
        coding.iter().zip(expected).zip(differential_rows)
    {
        assert_eq!(cds.location.raw(), location);
        assert_eq!(
            cds.qualifiers
                .iter()
                .find(|qualifier| qualifier.name == "protein_id")
                .unwrap()
                .value
                .as_deref(),
            Some(protein_id)
        );
        let published = cds
            .qualifiers
            .iter()
            .find(|qualifier| qualifier.name == "translation")
            .unwrap()
            .value
            .as_deref()
            .unwrap();
        let published: String = published
            .chars()
            .filter(|symbol| !symbol.is_ascii_whitespace())
            .collect();
        let result = cds
            .translate_cds_with(
                record.record().sequence(),
                options(GeneticCode::Standard, CdsReferencePolicy::RequireMatch),
            )
            .unwrap();
        assert_eq!(result.peptide(), published, "{protein_id}");
        assert_eq!(reference["peptide"], result.peptide(), "{protein_id}");
        assert_eq!(reference["reference_match"], true, "{protein_id}");
        assert_eq!(reference["genetic_code"], result.genetic_code().id());
        assert_eq!(reference["codon_start"], result.codon_start());
        assert_eq!(reference["record_id"], record.record().id());
        let original_index = record
            .features()
            .iter()
            .position(|feature| feature == *cds)
            .unwrap()
            + 1;
        assert_eq!(reference["feature_index"], original_index);
        assert_eq!(result.peptide().len(), length, "{protein_id}");
        assert_eq!(result.reference_status(), CdsReferenceStatus::Matched);
        assert_eq!(result.genetic_code(), GeneticCode::BacterialArchaealPlastid);
        assert_eq!(result.feature_length(), 3 * (length + 1));
        if protein_id == "NP_040709.1" {
            assert_eq!(
                result.peptide(),
                concat!(
                    "MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNL",
                    "RKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE"
                )
            );
        }
    }
    assert_eq!(record, &original);
}

// Each pair is one independently selected, synonymous-free encoding from the
// standard genetic code. Generated peptide expectations come from the selected
// residues, never from another call to the translation engine.
const ENCODINGS: &[(&str, &str)] = &[
    ("GCT", "A"),
    ("CGT", "R"),
    ("AAT", "N"),
    ("GAT", "D"),
    ("TGT", "C"),
    ("CAA", "Q"),
    ("GAA", "E"),
    ("GGT", "G"),
    ("CAT", "H"),
    ("ATT", "I"),
    ("CTG", "L"),
    ("AAA", "K"),
    ("ATG", "M"),
    ("TTT", "F"),
    ("CCT", "P"),
    ("TCT", "S"),
    ("ACT", "T"),
    ("TGG", "W"),
    ("TAT", "Y"),
    ("GTT", "V"),
];

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn generated_cds_survives_segmentation_orientation_and_post_extraction_frame(
        residues in prop::collection::vec(0usize..ENCODINGS.len(), 0..40),
        leading in 0usize..3,
        reverse in any::<bool>(),
        split_seed in any::<usize>(),
    ) {
        let mut coding = "C".repeat(leading);
        coding.push_str("ATG");
        let mut expected = String::from("M");
        for residue in residues {
            coding.push_str(ENCODINGS[residue].0);
            expected.push_str(ENCODINGS[residue].1);
        }
        coding.push_str("TAA");
        let genomic = if reverse { reverse_complement(&coding) } else { coding.clone() };
        let split = 1 + split_seed % (genomic.len() - 1);
        let parent = format!("CC{}GGG{}TT", &genomic[..split], &genomic[split..]);
        let join = format!("join(3..{},{}..{})", 2 + split, 6 + split, 5 + genomic.len());
        let location = if reverse { format!("complement({join})") } else { join };
        let codon_start = (leading + 1).to_string();
        let cds = feature(&location, &[("codon_start", Some(&codon_start)), ("translation", Some(&expected))]);
        let source = dna(&parent).unwrap();
        let before = source.clone();
        let result = cds.translate_cds(&source).unwrap();
        prop_assert_eq!(result.peptide(), expected.as_str());
        prop_assert_eq!(result.feature_length(), coding.len());
        prop_assert_eq!(result.translation().skipped_leading_bases(), leading);
        prop_assert_eq!(result.reference_status(), CdsReferenceStatus::Matched);
        prop_assert_eq!(&source, &before);
        let rna_source = rna(parent.replace('T', "U")).unwrap();
        let rna_result = cds.translate_cds(&rna_source).unwrap();
        prop_assert_eq!(rna_result.peptide(), expected.as_str());
    }

    #[test]
    fn every_generated_reference_substitution_is_reported_at_its_actual_position(
        residues in prop::collection::vec(0usize..ENCODINGS.len(), 0..40),
        position_seed in any::<usize>(),
    ) {
        let mut coding = String::from("ATG");
        let mut expected = String::from("M");
        for residue in residues {
            coding.push_str(ENCODINGS[residue].0);
            expected.push_str(ENCODINGS[residue].1);
        }
        coding.push_str("TAA");
        let position = position_seed % expected.len();
        let replacement = if expected.as_bytes()[position] == b'A' { "C" } else { "A" };
        let mut reference = expected.clone();
        reference.replace_range(position..=position, replacement);
        let cds = complete_feature(&coding, &[("translation", Some(&reference))]);
        match cds.translate_cds(&dna(&coding).unwrap()).unwrap_err() {
            CdsTranslationError::ReferenceMismatch { reference: actual, computed, first_difference } => {
                prop_assert_eq!(actual, reference);
                prop_assert_eq!(computed, expected);
                prop_assert_eq!(first_difference, position);
            }
            other => return Err(TestCaseError::fail(format!("unexpected error: {other:?}"))),
        }
    }
}

#[test]
fn ribosomal_slippage_translates_the_annotated_join() {
    // join(1..5,5..11) reads base 5 twice, a -1 frameshift:
    // ATGGC + CAAATAA = ATG GCC AAA TAA.
    let source = dna("ATGGCAAATAA").unwrap();
    let cds = feature(
        "join(1..5,5..11)",
        &[("ribosomal_slippage", None), ("translation", Some("MAK"))],
    );
    let product = cds.translate_cds(&source).unwrap();
    assert_eq!(product.peptide(), "MAK");
    assert!(product.ribosomal_slippage());
    assert_eq!(product.reference_status(), CdsReferenceStatus::Matched);
    assert_eq!(product.feature_length(), 12);

    // The same location without the qualifier translates identically but does
    // not report it; the qualifier changes reporting, not the arithmetic.
    let plain = feature("join(1..5,5..11)", &[])
        .translate_cds(&source)
        .unwrap();
    assert_eq!(plain.peptide(), "MAK");
    assert!(!plain.ribosomal_slippage());

    // Complemented joins also carry the shift: base 7 is read twice, and the
    // reverse complement of TTATTTG + GCCAT is ATG GCC AAA TAA.
    let reverse = dna("TTATTTGCCAT").unwrap();
    let complemented = feature(
        "complement(join(1..7,7..11))",
        &[("ribosomal_slippage", None)],
    );
    assert_eq!(
        complemented.translate_cds(&reverse).unwrap().peptide(),
        "MAK"
    );
}

#[test]
fn ribosomal_slippage_is_one_flag_and_translates_the_location_as_annotated() {
    let source = dna("ATGGCTTAA").unwrap();
    // Some records flag a CDS whose own location carries no shift; it translates
    // as annotated, and its reference is checked as usual.
    let unjoined = complete_feature(
        "ATGGCTTAA",
        &[("ribosomal_slippage", None), ("translation", Some("MA"))],
    );
    let product = unjoined.translate_cds(&source).unwrap();
    assert_eq!(product.peptide(), "MA");
    assert!(product.ribosomal_slippage());
    // A product read without the shift it needs fails the usual checks.
    let mismatched = complete_feature(
        "ATGGCTTAA",
        &[("ribosomal_slippage", None), ("translation", Some("MAK"))],
    );
    assert!(matches!(
        mismatched.translate_cds(&source).unwrap_err(),
        CdsTranslationError::ReferenceMismatch { .. }
    ));
    for value in ["", "retained text"] {
        let valued = feature("join(1..3,4..9)", &[("ribosomal_slippage", Some(value))]);
        assert_invalid_qualifier(
            valued.translate_cds(&source).unwrap_err(),
            "ribosomal_slippage",
        );
    }
    let repeated = feature(
        "join(1..3,4..9)",
        &[("ribosomal_slippage", None), ("ribosomal_slippage", None)],
    );
    match repeated.translate_cds(&source).unwrap_err() {
        CdsTranslationError::DuplicateQualifier { name } => assert_eq!(name, "ribosomal_slippage"),
        other => panic!("expected a duplicate /ribosomal_slippage, got {other:?}"),
    }
}

#[test]
fn newly_supported_tables_translate_their_reassigned_codons() {
    // Code 4 reads TGA as Trp; code 6 reads TAA and TAG as Gln and stops only at
    // TGA; code 16 reads TAG as Leu.
    for (table, sequence, peptide) in [
        ("4", "ATGTGATAA", "MW"),
        ("6", "ATGTAATGA", "MQ"),
        ("16", "ATGTAGTAA", "ML"),
    ] {
        let source = dna(sequence).unwrap();
        let cds = complete_feature(
            sequence,
            &[
                ("transl_table", Some(table)),
                ("translation", Some(peptide)),
            ],
        );
        let product = cds.translate_cds(&source).unwrap();
        assert_eq!(product.peptide(), peptide, "table {table}");
        assert_eq!(product.genetic_code().id().to_string(), table);
    }
}

#[test]
fn conditional_stops_end_only_a_complete_cds() {
    // Tables 27, 28 and 31: an internal codon that can end translation reads as
    // its amino acid, and the same codon can be the terminal stop.
    for (table, sequence, peptide) in [
        ("27", "ATGTGAGCTTGA", "MWA"),
        ("28", "ATGTAATAGTGA", "MQQ"),
        ("31", "ATGTAATAG", "ME"),
    ] {
        let source = dna(sequence).unwrap();
        let cds = complete_feature(sequence, &[("transl_table", Some(table))]);
        assert_eq!(
            cds.translate_cds(&source).unwrap().peptide(),
            peptide,
            "table {table}"
        );
    }
    // A codon that can never end translation is still not a terminal stop.
    let source = dna("ATGTGAGCT").unwrap();
    let cds = complete_feature("ATGTGAGCT", &[("transl_table", Some("27"))]);
    assert!(matches!(
        cds.translate_cds(&source).unwrap_err(),
        CdsTranslationError::Translation(TranslationError::MissingTerminalStop { .. })
    ));
}
