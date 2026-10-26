//! Independent NCBI-table, ambiguity, frame and compatibility checks.

use std::collections::BTreeSet;
use std::error::Error;
use std::io::Cursor;
use std::path::Path;
use std::process::Command;

use bioferris::SequenceError;
use bioferris::io::{FastaReader, one_record};
use bioferris::prelude::*;
use proptest::prelude::*;

const CODES: [GeneticCode; 3] = [
    GeneticCode::Standard,
    GeneticCode::VertebrateMitochondrial,
    GeneticCode::BacterialArchaealPlastid,
];
const FRAMES: [ReadingFrame; 3] = [ReadingFrame::Zero, ReadingFrame::One, ReadingFrame::Two];
const MODES: [TranslationMode; 3] = [
    TranslationMode::Full,
    TranslationMode::ToStop,
    TranslationMode::Cds,
];

// Biological data from the public NCBI table page, revision 2024-09-23,
// consulted 2026-10-15. Groups and explicit reassignments are independent of
// the production engine's indexed 64-byte tables and bit-mask expansion.
// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
const STANDARD_GROUPS: &[(u8, &str)] = &[
    (b'A', "GCT GCC GCA GCG"),
    (b'R', "CGT CGC CGA CGG AGA AGG"),
    (b'N', "AAT AAC"),
    (b'D', "GAT GAC"),
    (b'C', "TGT TGC"),
    (b'Q', "CAA CAG"),
    (b'E', "GAA GAG"),
    (b'G', "GGT GGC GGA GGG"),
    (b'H', "CAT CAC"),
    (b'I', "ATT ATC ATA"),
    (b'L', "TTA TTG CTT CTC CTA CTG"),
    (b'K', "AAA AAG"),
    (b'M', "ATG"),
    (b'F', "TTT TTC"),
    (b'P', "CCT CCC CCA CCG"),
    (b'S', "TCT TCC TCA TCG AGT AGC"),
    (b'T', "ACT ACC ACA ACG"),
    (b'W', "TGG"),
    (b'Y', "TAT TAC"),
    (b'V', "GTT GTC GTA GTG"),
    (b'*', "TAA TAG TGA"),
];

const BASE_SETS: &[(u8, &str)] = &[
    (b'A', "A"),
    (b'C', "C"),
    (b'G', "G"),
    (b'T', "T"),
    (b'R', "AG"),
    (b'Y', "CT"),
    (b'S', "CG"),
    (b'W', "AT"),
    (b'K', "GT"),
    (b'M', "AC"),
    (b'B', "CGT"),
    (b'D', "AGT"),
    (b'H', "ACT"),
    (b'V', "ACG"),
    (b'N', "ACGT"),
];

fn config(genetic_code: GeneticCode, mode: TranslationMode) -> TranslationConfig {
    TranslationConfig {
        genetic_code,
        mode,
        ..TranslationConfig::default()
    }
}

fn reference_residue(codon: &[u8], code: GeneticCode) -> u8 {
    if code == GeneticCode::VertebrateMitochondrial {
        match codon {
            b"ATA" => return b'M',
            b"TGA" => return b'W',
            b"AGA" | b"AGG" => return b'*',
            _ => {}
        }
    }
    STANDARD_GROUPS
        .iter()
        .find_map(|(residue, group)| {
            group
                .split_whitespace()
                .any(|entry| entry.as_bytes() == codon)
                .then_some(*residue)
        })
        .expect("independent groups cover every concrete codon")
}

fn reference_start(codon: &[u8], code: GeneticCode) -> bool {
    match code {
        GeneticCode::Standard => [*b"TTG", *b"CTG", *b"ATG"].contains(&codon.try_into().unwrap()),
        GeneticCode::VertebrateMitochondrial => {
            [*b"ATT", *b"ATC", *b"ATA", *b"ATG", *b"GTG"].contains(&codon.try_into().unwrap())
        }
        GeneticCode::BacterialArchaealPlastid => [
            *b"TTG", *b"CTG", *b"ATT", *b"ATC", *b"ATA", *b"ATG", *b"GTG",
        ]
        .contains(&codon.try_into().unwrap()),
        // The other tables are checked against NCBI gc.prt in ncbi_genetic_codes.rs.
        other => panic!("no independent initiation set for table {}", other.id()),
    }
}

fn reference_outcome(codon: &[u8], code: GeneticCode) -> (u8, bool) {
    let expansions: Vec<_> = codon
        .iter()
        .map(|base| {
            BASE_SETS
                .iter()
                .find(|(symbol, _)| symbol == base)
                .unwrap()
                .1
        })
        .collect();
    let mut residues = BTreeSet::new();
    let mut all_start = true;
    for first in expansions[0].bytes() {
        for second in expansions[1].bytes() {
            for third in expansions[2].bytes() {
                let concrete = [first, second, third];
                residues.insert(reference_residue(&concrete, code));
                all_start &= reference_start(&concrete, code);
            }
        }
    }
    let outcomes: Vec<_> = residues.into_iter().collect();
    let residue = match outcomes.as_slice() {
        [residue] => *residue,
        [b'D', b'N'] => b'B',
        [b'E', b'Q'] => b'Z',
        [b'I', b'L'] => b'J',
        _ => b'X',
    };
    (residue, all_start)
}

fn all_iupac_codons() -> Vec<[u8; 3]> {
    let mut codons = Vec::new();
    for (first, _) in BASE_SETS {
        for (second, _) in BASE_SETS {
            for (third, _) in BASE_SETS {
                codons.push([*first, *second, *third]);
            }
        }
    }
    assert_eq!(codons.len(), 3375);
    codons
}

fn rna_error(error: TranslationError) -> TranslationError {
    let to_rna = |codon: [u8; 3]| codon.map(|base| if base == b'T' { b'U' } else { base });
    match error {
        TranslationError::InvalidStartCodon { codon } => TranslationError::InvalidStartCodon {
            codon: to_rna(codon),
        },
        TranslationError::InvalidStartCodonForCode {
            codon,
            genetic_code,
        } => TranslationError::InvalidStartCodonForCode {
            codon: to_rna(codon),
            genetic_code,
        },
        TranslationError::MissingTerminalStop { codon } => TranslationError::MissingTerminalStop {
            codon: to_rna(codon),
        },
        other => other,
    }
}

fn assert_pair(input: &str, options: TranslationConfig, expected: Result<&str, TranslationError>) {
    let dna_sequence = dna(input).unwrap();
    let rna_sequence = rna(input.replace('T', "U")).unwrap();
    for (is_rna, actual, expected) in [
        (
            false,
            dna_sequence.translate_with(options),
            expected.clone(),
        ),
        (
            true,
            rna_sequence.translate_with(options),
            expected.map_err(rna_error),
        ),
    ] {
        match expected {
            Ok(peptide) => {
                let translated = actual
                    .unwrap_or_else(|error| panic!("{input}, {options:?}, RNA={is_rna}: {error}"));
                assert_eq!(translated.as_str(), peptide, "{input}, {options:?}");
                assert_eq!(translated.genetic_code(), options.genetic_code);
                assert_eq!(translated.reading_frame(), options.reading_frame);
                assert_eq!(
                    translated.skipped_leading_bases(),
                    options.reading_frame.offset()
                );
                assert_eq!(
                    translated.ignored_trailing_bases(),
                    (input.len() - options.reading_frame.offset()) % 3
                );
            }
            Err(error) => {
                assert_eq!(actual, Err(error), "{input}, {options:?}, RNA={is_rna}");
            }
        }
    }
    assert_eq!(dna_sequence.to_string(), input.to_ascii_uppercase());
}

fn start_error(codon: [u8; 3], genetic_code: GeneticCode) -> TranslationError {
    if genetic_code == GeneticCode::Standard {
        TranslationError::InvalidStartCodon { codon }
    } else {
        TranslationError::InvalidStartCodonForCode {
            codon,
            genetic_code,
        }
    }
}

#[test]
fn defaults_numeric_selectors_and_table_names_are_explicit() {
    assert_eq!(GeneticCode::default(), GeneticCode::Standard);
    assert_eq!(ReadingFrame::default(), ReadingFrame::Zero);
    let defaults = TranslationConfig::default();
    assert_eq!(defaults.genetic_code, GeneticCode::Standard);
    assert_eq!(defaults.reading_frame, ReadingFrame::Zero);
    assert_eq!(defaults.mode, TranslationMode::Full);
    assert_eq!(defaults.partial_codons, PartialCodonPolicy::Reject);
    for (code, id, name) in [
        (GeneticCode::Standard, 1, "Standard"),
        (
            GeneticCode::VertebrateMitochondrial,
            2,
            "Vertebrate Mitochondrial",
        ),
        (
            GeneticCode::BacterialArchaealPlastid,
            11,
            "Bacterial, Archaeal and Plant Plastid",
        ),
    ] {
        assert_eq!(code.id(), id);
        assert_eq!(code.name(), name);
        assert_eq!(GeneticCode::try_from(id), Ok(code));
    }
    let supported = GeneticCode::ALL.map(GeneticCode::id);
    assert_eq!(
        supported,
        [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33
        ]
    );
    for id in 0..=u16::MAX {
        if !supported.contains(&id) {
            assert_eq!(
                GeneticCode::try_from(id),
                Err(TranslationError::UnsupportedGeneticCode { id })
            );
        }
    }
    for (offset, frame) in FRAMES.into_iter().enumerate() {
        assert_eq!(frame.offset(), offset);
        assert_eq!(ReadingFrame::try_from(offset), Ok(frame));
    }
    for frame in [3, 4, 100, usize::MAX] {
        assert_eq!(
            ReadingFrame::try_from(frame),
            Err(TranslationError::InvalidReadingFrame { frame })
        );
    }
}

#[test]
fn all_64_concrete_codons_match_each_published_table_in_dna_and_rna() {
    for code in CODES {
        let mut seen = BTreeSet::new();
        for (_, group) in STANDARD_GROUPS {
            for codon in group.split_whitespace() {
                assert!(seen.insert(codon));
                let expected = char::from(reference_residue(codon.as_bytes(), code)).to_string();
                assert_pair(codon, config(code, TranslationMode::Full), Ok(&expected));
            }
        }
        assert_eq!(seen.len(), 64);
    }
}

#[test]
fn all_3375_iupac_codons_match_independent_expansion_in_every_table_and_molecule() {
    for code in CODES {
        for codon in all_iupac_codons() {
            let text = std::str::from_utf8(&codon).unwrap();
            let (residue, _) = reference_outcome(&codon, code);
            let full = char::from(residue).to_string();
            assert_pair(text, config(code, TranslationMode::Full), Ok(&full));
            let to_stop = if residue == b'*' {
                String::new()
            } else {
                format!("{full}K")
            };
            assert_pair(
                &format!("{text}AAA"),
                config(code, TranslationMode::ToStop),
                Ok(&to_stop),
            );
        }
    }
}

#[test]
fn all_3375_iupac_starts_terminals_and_internal_stops_use_the_selected_table() {
    for code in CODES {
        for codon in all_iupac_codons() {
            let text = std::str::from_utf8(&codon).unwrap();
            let (residue, all_start) = reference_outcome(&codon, code);
            let start = if all_start {
                Ok("MK")
            } else {
                Err(start_error(codon, code))
            };
            assert_pair(
                &format!("{text}AAATAA"),
                config(code, TranslationMode::Cds),
                start,
            );
            let terminal = if residue == b'*' {
                Ok("MK")
            } else {
                Err(TranslationError::MissingTerminalStop { codon })
            };
            assert_pair(
                &format!("ATGAAA{text}"),
                config(code, TranslationMode::Cds),
                terminal,
            );
            let peptide = format!("M{}", char::from(residue));
            let interior = if residue == b'*' {
                Err(TranslationError::InternalStop { offset: 3 })
            } else {
                Ok(peptide.as_str())
            };
            assert_pair(
                &format!("ATG{text}TAA"),
                config(code, TranslationMode::Cds),
                interior,
            );
        }
    }
}

#[test]
fn initiation_recodes_only_the_first_codon_and_all_expansions_must_be_starts() {
    for code in CODES {
        for (codon, allowed) in [
            ("HTG", code != GeneticCode::VertebrateMitochondrial),
            ("ATN", code != GeneticCode::Standard),
            ("NTG", code == GeneticCode::BacterialArchaealPlastid),
        ] {
            let bytes = codon.as_bytes().try_into().unwrap();
            let expected = if allowed {
                Ok("MK")
            } else {
                Err(start_error(bytes, code))
            };
            assert_pair(
                &format!("{codon}AAATAA"),
                config(code, TranslationMode::Cds),
                expected,
            );
        }
    }
    assert_pair(
        "GTGGTGTAA",
        config(GeneticCode::BacterialArchaealPlastid, TranslationMode::Cds),
        Ok("MV"),
    );
    assert_pair(
        "ATGATATAA",
        config(GeneticCode::BacterialArchaealPlastid, TranslationMode::Cds),
        Ok("MI"),
    );
    assert_pair(
        "ATGATATAA",
        config(GeneticCode::VertebrateMitochondrial, TranslationMode::Cds),
        Ok("MM"),
    );
}

#[test]
fn ambiguity_keeps_b_z_j_and_mixed_sense_stop_x_without_false_termination() {
    for code in CODES {
        assert_pair("RAYSARMTY", config(code, TranslationMode::Full), Ok("BZJ"));
        assert_pair("TANAAA", config(code, TranslationMode::ToStop), Ok("XK"));
        assert_pair("ATGTANTAA", config(code, TranslationMode::Cds), Ok("MX"));
        assert_pair("TARAAA", config(code, TranslationMode::ToStop), Ok(""));
    }
    assert_pair(
        "AGRTRA",
        config(GeneticCode::Standard, TranslationMode::Full),
        Ok("R*"),
    );
    assert_pair(
        "AGRTRA",
        config(GeneticCode::VertebrateMitochondrial, TranslationMode::Full),
        Ok("*X"),
    );
}

#[test]
fn published_biopython_188_translation_example_matches_supported_tables() {
    // Public Bio.Seq documentation, not execution or source of its implementation.
    let input = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
    for (code, mode, expected) in [
        (
            GeneticCode::Standard,
            TranslationMode::Full,
            "VAIVMGR*KGAR*",
        ),
        (GeneticCode::Standard, TranslationMode::ToStop, "VAIVMGR"),
        (
            GeneticCode::VertebrateMitochondrial,
            TranslationMode::Full,
            "VAIVMGRWKGAR*",
        ),
        (
            GeneticCode::VertebrateMitochondrial,
            TranslationMode::ToStop,
            "VAIVMGRWKGAR",
        ),
        (
            GeneticCode::VertebrateMitochondrial,
            TranslationMode::Cds,
            "MAIVMGRWKGAR",
        ),
    ] {
        assert_pair(input, config(code, mode), Ok(expected));
    }
    assert_pair(
        input,
        config(GeneticCode::Standard, TranslationMode::Cds),
        Err(TranslationError::InvalidStartCodon { codon: *b"GTG" }),
    );
}

#[test]
fn frames_skip_a_prefix_and_report_trailing_remainder_separately() {
    for code in CODES {
        for frame in FRAMES {
            let prefix = "CC".get(..frame.offset()).unwrap();
            let options = TranslationConfig {
                reading_frame: frame,
                ..config(code, TranslationMode::Cds)
            };
            assert_pair(&format!("{prefix}ATGAAATAA"), options, Ok("MK"));
            for tail in ["", "A", "AC"] {
                let options = TranslationConfig {
                    mode: TranslationMode::Full,
                    partial_codons: PartialCodonPolicy::Ignore,
                    ..options
                };
                assert_pair(&format!("{prefix}ATGAAATAA{tail}"), options, Ok("MK*"));
            }
        }
    }
}

#[test]
fn frame_errors_use_suffix_lengths_but_original_internal_stop_offsets() {
    for code in CODES {
        for frame in FRAMES {
            let prefix = "CC".get(..frame.offset()).unwrap();
            let options = TranslationConfig {
                reading_frame: frame,
                ..config(code, TranslationMode::Cds)
            };
            assert_pair(
                &format!("{prefix}ATGTAGTAA"),
                options,
                Err(TranslationError::InternalStop {
                    offset: frame.offset() + 3,
                }),
            );
            assert_pair(
                &format!("{prefix}ATG"),
                options,
                Err(TranslationError::CdsTooShort { length: 3 }),
            );
            for policy in [PartialCodonPolicy::Reject, PartialCodonPolicy::Ignore] {
                let options = TranslationConfig {
                    partial_codons: policy,
                    ..options
                };
                assert_pair(
                    &format!("{prefix}ATGAAATAAA"),
                    options,
                    Err(TranslationError::PartialCodon {
                        length: 10,
                        remainder: 1,
                    }),
                );
                assert_pair(
                    &format!("{prefix}AA"),
                    options,
                    Err(TranslationError::PartialCodon {
                        length: 2,
                        remainder: 2,
                    }),
                );
            }
        }
    }
}

#[test]
fn equal_length_frame_selects_empty_suffix_and_larger_frames_are_errors() {
    for code in CODES {
        for frame in FRAMES {
            for length in 0..=2 {
                for mode in MODES {
                    for partial_codons in [PartialCodonPolicy::Reject, PartialCodonPolicy::Ignore] {
                        let options = TranslationConfig {
                            reading_frame: frame,
                            partial_codons,
                            ..config(code, mode)
                        };
                        let input = "AC".get(..length).unwrap();
                        match frame.offset().cmp(&length) {
                            std::cmp::Ordering::Greater => assert_pair(
                                input,
                                options,
                                Err(TranslationError::FrameOutOfBounds {
                                    frame: frame.offset(),
                                    length,
                                }),
                            ),
                            std::cmp::Ordering::Equal => {
                                let expected = if mode == TranslationMode::Cds {
                                    Err(TranslationError::CdsTooShort { length: 0 })
                                } else {
                                    Ok("")
                                };
                                assert_pair(input, options, expected);
                            }
                            std::cmp::Ordering::Less => {}
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn to_stop_checks_partial_policy_before_early_exit_and_reports_only_remainder() {
    for code in CODES {
        let options = TranslationConfig {
            reading_frame: ReadingFrame::One,
            ..config(code, TranslationMode::ToStop)
        };
        assert_pair(
            "CATGTAAGGGA",
            options,
            Err(TranslationError::PartialCodon {
                length: 10,
                remainder: 1,
            }),
        );
        assert_pair(
            "CATGTAAGGGA",
            TranslationConfig {
                partial_codons: PartialCodonPolicy::Ignore,
                ..options
            },
            Ok("M"),
        );
    }
}

fn check_legacy_equivalence(input: &str) {
    let sequence = dna(input).unwrap();
    let transcript = sequence.transcribe().unwrap();
    for mode in MODES {
        for partial_codons in [PartialCodonPolicy::Reject, PartialCodonPolicy::Ignore] {
            // Deliberately no struct update: existing two-field literals remain valid.
            let options = TranslationOptions {
                mode,
                partial_codons,
            };
            let configured = TranslationConfig::from(options);
            assert_eq!(configured.genetic_code, GeneticCode::Standard);
            assert_eq!(configured.reading_frame, ReadingFrame::Zero);
            assert_eq!(
                sequence.translate_with(configured),
                sequence.translate_standard(options),
                "{input}"
            );
            assert_eq!(
                transcript.translate_with(configured),
                transcript.translate_standard(options),
                "{input}"
            );
        }
    }
    assert_eq!(
        sequence.translate_with(TranslationConfig::default()),
        sequence.translate()
    );
    assert_eq!(
        sequence.translate_with(config(GeneticCode::Standard, TranslationMode::Cds)),
        sequence.translate_cds()
    );
    assert_eq!(
        sequence.translate_with(config(GeneticCode::Standard, TranslationMode::ToStop)),
        sequence.translate_to_stop()
    );
}

#[test]
fn configured_standard_defaults_preserve_old_results_and_errors() {
    for input in [
        "",
        "A",
        "AT",
        "ATG",
        "ATGTAA",
        "GTGAAATAA",
        "ATGTGATAA",
        "ATGAAA",
        "ATGTAAA",
        "ATGTAAAA",
    ] {
        check_legacy_equivalence(input);
    }
    for codon in all_iupac_codons() {
        let codon = std::str::from_utf8(&codon).unwrap();
        check_legacy_equivalence(&format!("{codon}ATGTAA"));
    }
}

struct LegacySequence(DnaSequence);

// Existing implementors do not need to implement ConfiguredTranslationExt.
impl NucleotideSequenceExt for LegacySequence {
    fn complement(&self) -> Result<Self, SequenceError> {
        self.0.complement().map(Self)
    }

    fn reverse_complement(&self) -> Result<Self, SequenceError> {
        self.0.reverse_complement().map(Self)
    }

    fn translate_standard(
        &self,
        options: TranslationOptions,
    ) -> Result<Translation, TranslationError> {
        self.0.translate_standard(options)
    }
}

#[test]
fn existing_nucleotide_extension_implementors_need_no_new_required_method() {
    let sequence = LegacySequence(dna("ATGAAATAA").unwrap());
    assert_eq!(sequence.translate().unwrap().as_str(), "MK*");
    assert_eq!(sequence.translate_cds().unwrap().as_str(), "MK");
    assert_eq!(sequence.translate_to_stop().unwrap().as_str(), "MK");
    assert_eq!(sequence.complement().unwrap().0.to_string(), "TACTTTATT");
    assert_eq!(
        sequence.reverse_complement().unwrap().0.to_string(),
        "TTATTTCAT"
    );
}

#[test]
fn complete_input_validation_rejects_invalid_skipped_and_ignored_symbols() {
    let options = TranslationConfig {
        reading_frame: ReadingFrame::Two,
        partial_codons: PartialCodonPolicy::Ignore,
        ..TranslationConfig::default()
    };
    for (input, expected_offset) in [("!CATGA", 0), ("C!ATGA", 1), ("CCATG!", 5), ("CCATGU", 5)] {
        let error = dna(input)
            .map(|sequence| sequence.translate_with(options))
            .unwrap_err();
        assert_eq!(error.index(), expected_offset);
    }
    for (input, expected_offset) in [("!CAUGA", 0), ("C!AUGA", 1), ("CCAUG!", 5), ("CCAUGT", 5)] {
        let error = rna(input)
            .map(|sequence| sequence.translate_with(options))
            .unwrap_err();
        assert_eq!(error.index(), expected_offset);
    }
}

#[test]
fn translation_product_retains_metadata_and_requires_stop_free_protein_conversion() {
    let input = dna("CATATGAAGA").unwrap();
    let options = TranslationConfig {
        genetic_code: GeneticCode::VertebrateMitochondrial,
        reading_frame: ReadingFrame::One,
        ..TranslationConfig::default()
    };
    let product = input.translate_with(options).unwrap();
    assert_eq!(product.as_str(), "MW*");
    assert_eq!(product.to_string(), "MW*");
    assert_eq!(product.len(), 3);
    assert!(!product.is_empty());
    assert_eq!(product.genetic_code().id(), 2);
    assert_eq!(product.reading_frame().offset(), 1);
    let error = product.to_protein().unwrap_err();
    assert_eq!((error.index(), error.byte()), (2, b'*'));
    let stopped = input
        .translate_with(TranslationConfig {
            mode: TranslationMode::ToStop,
            ..options
        })
        .unwrap();
    assert_eq!(stopped.to_protein().unwrap().to_string(), "MW");
    assert_eq!(input.to_string(), "CATATGAAGA");
}

#[test]
fn typed_errors_preserve_precise_fields_and_human_messages() {
    for (error, text) in [
        (
            TranslationError::UnsupportedGeneticCode { id: 7 },
            "unsupported genetic code 7; supported NCBI tables are 1-6, 9-16 and 21-33",
        ),
        (
            TranslationError::InvalidReadingFrame { frame: 3 },
            "invalid reading frame 3; expected 0, 1 or 2",
        ),
        (
            TranslationError::FrameOutOfBounds {
                frame: 2,
                length: 1,
            },
            "reading frame 2 exceeds sequence length 1",
        ),
        (
            TranslationError::InvalidStartCodonForCode {
                codon: *b"TTG",
                genetic_code: GeneticCode::VertebrateMitochondrial,
            },
            "invalid table 2 start codon TTG",
        ),
        (
            TranslationError::InvalidStartCodonForCode {
                codon: *b"TTA",
                genetic_code: GeneticCode::BacterialArchaealPlastid,
            },
            "invalid table 11 start codon TTA",
        ),
        (
            TranslationError::InvalidStartCodon { codon: *b"GTG" },
            "invalid table 1 start codon GTG",
        ),
    ] {
        assert_eq!(error.to_string(), text);
        assert!(error.source().is_none());
    }
}

#[test]
fn real_phix174_endolysin_cds_matches_annotated_table_11_protein() {
    let record = one_record(FastaReader::<_, Dna>::new(Cursor::new(include_bytes!(
        "fixtures/phix174-e.fasta"
    ))))
    .unwrap();
    assert_eq!(record.sequence().len(), 276);
    let options = config(GeneticCode::BacterialArchaealPlastid, TranslationMode::Cds);
    let product = record.sequence().translate_with(options).unwrap();
    let expected = concat!(
        "MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNL",
        "RKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE"
    );
    assert_eq!(product.as_str(), expected);
    assert_eq!(product.len(), 91);
    assert_eq!(
        product.genetic_code(),
        GeneticCode::BacterialArchaealPlastid
    );
    assert_eq!(product.to_protein().unwrap().to_string(), expected);
    assert_eq!(
        record
            .sequence()
            .transcribe()
            .unwrap()
            .translate_with(options)
            .unwrap(),
        product
    );
}

fn assert_shifted_result(
    actual: Result<Translation, TranslationError>,
    expected: Result<Translation, TranslationError>,
    code: GeneticCode,
    frame: ReadingFrame,
) {
    match expected {
        Ok(expected) => {
            let actual = actual.unwrap();
            assert_eq!(actual.as_str(), expected.as_str());
            assert_eq!(
                actual.ignored_trailing_bases(),
                expected.ignored_trailing_bases()
            );
            assert_eq!(actual.skipped_leading_bases(), frame.offset());
            assert_eq!(actual.reading_frame(), frame);
            assert_eq!(actual.genetic_code(), code);
        }
        Err(TranslationError::InternalStop { offset }) => assert_eq!(
            actual,
            Err(TranslationError::InternalStop {
                offset: offset + frame.offset()
            })
        ),
        Err(error) => assert_eq!(actual, Err(error)),
    }
}

proptest! {
    #[test]
    fn selected_frame_matches_explicit_suffix_without_losing_original_offsets(
        suffix in "[ACGTRYSWKMBDHVN]{0,75}",
        prefix in "[ACGTRYSWKMBDHVN]{2}",
        frame_index in 0_usize..3,
    ) {
        let frame = FRAMES[frame_index];
        let input = format!("{}{suffix}", &prefix[..frame.offset()]);
        for code in CODES {
            for mode in MODES {
                for partial_codons in [PartialCodonPolicy::Reject, PartialCodonPolicy::Ignore] {
                    let options = TranslationConfig { partial_codons, ..config(code, mode) };
                    let framed = TranslationConfig { reading_frame: frame, ..options };
                    assert_shifted_result(dna(&input).unwrap().translate_with(framed), dna(&suffix).unwrap().translate_with(options), code, frame);
                    assert_shifted_result(rna(input.replace('T', "U")).unwrap().translate_with(framed), rna(suffix.replace('T', "U")).unwrap().translate_with(options), code, frame);
                }
            }
        }
    }

    #[test]
    fn arbitrary_complete_codons_match_the_independent_oracle_and_normalize_case(
        codons in proptest::collection::vec("[ACGTRYSWKMBDHVN]{3}", 0..30),
    ) {
        let input: String = codons.iter().map(String::as_str).collect();
        for code in CODES {
            let full: String = codons.iter().map(|codon| char::from(reference_outcome(codon.as_bytes(), code).0)).collect();
            let to_stop = full.split('*').next().unwrap();
            assert_pair(&input, config(code, TranslationMode::Full), Ok(&full));
            assert_pair(&input, config(code, TranslationMode::ToStop), Ok(to_stop));
            let lowercase = dna(input.to_ascii_lowercase()).unwrap().translate_with(config(code, TranslationMode::Full)).unwrap();
            prop_assert_eq!(lowercase.as_str(), full.as_str());
        }
    }

    #[test]
    fn invalid_input_cannot_hide_in_a_frame_prefix_or_ignored_suffix(
        mut input in proptest::collection::vec(proptest::sample::select(b"ACGT".to_vec()), 3..60),
        position in any::<usize>(),
        frame_index in 0_usize..3,
    ) {
        let invalid_offset = position % input.len();
        input[invalid_offset] = b'!';
        let options = TranslationConfig { reading_frame: FRAMES[frame_index], partial_codons: PartialCodonPolicy::Ignore, ..TranslationConfig::default() };
        let error = dna(&input).map(|sequence| sequence.translate_with(options)).unwrap_err();
        prop_assert_eq!(error.index(), invalid_offset);
        prop_assert_eq!(error.byte(), b'!');
    }
}

fn guide_fence<'a>(guide: &'a str, marker: &str, language: &str) -> &'a str {
    let marker = format!("<!-- {marker} -->");
    assert_eq!(guide.matches(&marker).count(), 1);
    let (_, section) = guide.split_once(&marker).unwrap();
    let (_, block) = section.split_once(&format!("```{language}\n")).unwrap();
    block.split_once("```").unwrap().0
}

#[test]
fn genetic_code_guide_runs_the_actual_cargo_example() -> Result<(), Box<dyn Error>> {
    let guide = include_str!("../../docs/api/genetic-codes-v1.md");
    let program = include_str!("../examples/genetic_codes.rs");
    let (_, program) = program.split_once("\n\n").unwrap();
    assert_eq!(
        guide_fence(guide, "genetic-codes-program", "rust").trim(),
        program.trim()
    );
    assert_eq!(
        guide_fence(guide, "genetic-codes-command", "bash"),
        "cargo bio genetic_codes\n"
    );
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let output = Command::new(env!("CARGO"))
        .args(["bio", "genetic_codes", "--offline", "--locked"])
        .current_dir(root)
        .output()?;
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    assert_eq!(
        output.stdout,
        guide_fence(guide, "genetic-codes-output", "text").as_bytes()
    );
    Ok(())
}
