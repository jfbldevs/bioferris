//! Specification-derived translation, biological invariants and real CDS tests.

use std::collections::BTreeSet;

use bioferris::prelude::*;
use proptest::prelude::*;

fn options(mode: TranslationMode) -> TranslationOptions {
    TranslationOptions {
        mode,
        ..TranslationOptions::default()
    }
}

fn translate(input: &str, mode: TranslationMode) -> Result<Translation, TranslationError> {
    DnaSequence::try_from(input)
        .expect("test input is IUPAC DNA")
        .translate_standard(options(mode))
}

// Independently represented standard-code data: amino-acid groups rather than
// the production engine's 64-byte position-indexed table. Source: NCBI table 1.
const CODON_GROUPS: &[(u8, &str)] = &[
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

// Independent string-set IUPAC expansion, not the production bit masks.
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

fn reference_residue(codon: &[u8]) -> u8 {
    CODON_GROUPS
        .iter()
        .find_map(|(residue, codons)| {
            codons
                .split_whitespace()
                .any(|entry| entry.as_bytes() == codon)
                .then_some(*residue)
        })
        .expect("independent table covers every concrete codon")
}

fn expanded_outcomes(bases: [&str; 3]) -> (u8, bool) {
    let mut residues = BTreeSet::new();
    let mut all_start = true;
    for first in bases[0].bytes() {
        for second in bases[1].bytes() {
            for third in bases[2].bytes() {
                let codon = [first, second, third];
                residues.insert(reference_residue(&codon));
                all_start &= [*b"ATG", *b"TTG", *b"CTG"].contains(&codon);
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

#[test]
fn all_64_published_codons_match_for_dna_and_rna() {
    let mut codons_seen = BTreeSet::new();
    for (expected, codons) in CODON_GROUPS {
        for codon in codons.split_whitespace() {
            assert!(codons_seen.insert(codon));
            let dna = DnaSequence::try_from(codon).unwrap();
            let translation = dna
                .translate_standard(TranslationOptions::default())
                .unwrap();
            assert_eq!(translation.as_str().as_bytes(), &[*expected], "{codon}");
            assert_eq!(translation.ignored_trailing_bases(), 0);
            assert_eq!(
                dna.transcribe()
                    .unwrap()
                    .translate_standard(TranslationOptions::default())
                    .unwrap(),
                translation
            );
        }
    }
    assert_eq!(codons_seen.len(), 64);
}

#[test]
fn all_3375_iupac_codons_match_independent_set_expansion_and_cds_rules() {
    let mut checked = 0;
    for (first, first_set) in BASE_SETS {
        for (second, second_set) in BASE_SETS {
            for (third, third_set) in BASE_SETS {
                let bytes = [*first, *second, *third];
                let codon = std::str::from_utf8(&bytes).unwrap();
                let (expected, all_start) = expanded_outcomes([first_set, second_set, third_set]);
                let full = translate(codon, TranslationMode::Full).unwrap();
                assert_eq!(full.as_str().as_bytes(), &[expected], "{codon}");
                let rna = RnaSequence::try_from(codon.replace('T', "U")).unwrap();
                assert_eq!(
                    rna.translate_standard(TranslationOptions::default())
                        .unwrap(),
                    full,
                    "{codon}"
                );

                let start = translate(&format!("{codon}AAATAA"), TranslationMode::Cds);
                if all_start {
                    assert_eq!(start.unwrap().as_str(), "MK", "{codon}");
                } else {
                    assert!(
                        matches!(start, Err(TranslationError::InvalidStartCodon { .. })),
                        "{codon}"
                    );
                }
                let end = translate(&format!("ATGAAA{codon}"), TranslationMode::Cds);
                let interior = translate(&format!("ATG{codon}TAA"), TranslationMode::Cds);
                let to_stop = translate(&format!("{codon}AAA"), TranslationMode::ToStop).unwrap();
                if expected == b'*' {
                    assert_eq!(end.unwrap().as_str(), "MK", "{codon}");
                    assert_eq!(
                        interior,
                        Err(TranslationError::InternalStop { offset: 3 }),
                        "{codon}"
                    );
                    assert!(to_stop.is_empty(), "{codon}");
                } else {
                    assert!(
                        matches!(end, Err(TranslationError::MissingTerminalStop { .. })),
                        "{codon}"
                    );
                    assert_eq!(
                        interior.unwrap().as_str(),
                        format!("M{}", char::from(expected)),
                        "{codon}"
                    );
                    assert_eq!(
                        to_stop.as_str(),
                        format!("{}K", char::from(expected)),
                        "{codon}"
                    );
                }
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 3375);
}

#[test]
fn complements_cover_every_iupac_symbol_and_preserve_molecule() {
    let dna = DnaSequence::try_from("ACGTRYSWKMBDHVN").unwrap();
    assert_eq!(dna.complement().unwrap().to_string(), "TGCAYRSWMKVHDBN");
    assert_eq!(
        dna.reverse_complement().unwrap().to_string(),
        "NBDHVKMWSRYACGT"
    );
    let rna = dna.transcribe().unwrap();
    assert_eq!(rna.complement().unwrap().to_string(), "UGCAYRSWMKVHDBN");
    assert_eq!(
        rna.reverse_complement().unwrap().to_string(),
        "NBDHVKMWSRYACGU"
    );
    assert_eq!(dna.to_string(), "ACGTRYSWKMBDHVN");
}

#[test]
fn published_biopython_standard_code_example_and_stop_policy() {
    // Bio.Seq 1.88 public documentation; standard table only, not the table 2 example.
    let input = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
    assert_eq!(
        translate(input, TranslationMode::Full).unwrap().as_str(),
        "VAIVMGR*KGAR*"
    );
    assert_eq!(
        translate(input, TranslationMode::ToStop).unwrap().as_str(),
        "VAIVMGR"
    );
    assert_eq!(
        translate("GTGGCCATTGTAATGGGCCGC", TranslationMode::ToStop)
            .unwrap()
            .as_str(),
        "VAIVMGR"
    );
    assert_eq!(
        translate("TANTARGCN", TranslationMode::ToStop)
            .unwrap()
            .as_str(),
        "X"
    );
}

#[test]
fn cds_recodes_only_the_start_and_drops_only_the_terminal_stop() {
    assert_eq!(
        translate("TTGAAATAG", TranslationMode::Full)
            .unwrap()
            .as_str(),
        "LK*"
    );
    assert_eq!(
        translate("TTGAAATAG", TranslationMode::ToStop)
            .unwrap()
            .as_str(),
        "LK"
    );
    assert_eq!(
        translate("TTGAAATAG", TranslationMode::Cds)
            .unwrap()
            .as_str(),
        "MK"
    );
    assert_eq!(
        translate("HTGAAATAR", TranslationMode::Cds)
            .unwrap()
            .as_str(),
        "MK"
    );
    assert_eq!(
        translate("ATGTTGTAA", TranslationMode::Cds)
            .unwrap()
            .as_str(),
        "ML"
    );
    assert_eq!(
        translate("ATGTAA", TranslationMode::Cds).unwrap().as_str(),
        "M"
    );
}

#[test]
fn cds_reports_precise_structured_errors() {
    assert_eq!(
        translate("", TranslationMode::Cds),
        Err(TranslationError::CdsTooShort { length: 0 })
    );
    assert_eq!(
        translate("ATG", TranslationMode::Cds),
        Err(TranslationError::CdsTooShort { length: 3 })
    );
    assert_eq!(
        translate("GTGAAATAA", TranslationMode::Cds),
        Err(TranslationError::InvalidStartCodon { codon: *b"GTG" })
    );
    assert_eq!(
        translate("ATGAAATAN", TranslationMode::Cds),
        Err(TranslationError::MissingTerminalStop { codon: *b"TAN" })
    );
    assert_eq!(
        translate("ATGTAATAG", TranslationMode::Cds),
        Err(TranslationError::InternalStop { offset: 3 })
    );
    assert_eq!(
        translate("ATGTANTAA", TranslationMode::Cds)
            .unwrap()
            .as_str(),
        "MX"
    );
}

#[test]
fn partial_codons_require_explicit_policy_and_cannot_bypass_cds() {
    for input in ["A", "AC", "ATGA", "ATGAC", "TAAA"] {
        let sequence = DnaSequence::try_from(input).unwrap();
        let error = TranslationError::PartialCodon {
            length: input.len(),
            remainder: input.len() % 3,
        };
        for mode in [
            TranslationMode::Full,
            TranslationMode::ToStop,
            TranslationMode::Cds,
        ] {
            assert_eq!(
                sequence.translate_standard(options(mode)),
                Err(error.clone())
            );
            let ignored = sequence.translate_standard(TranslationOptions {
                mode,
                partial_codons: PartialCodonPolicy::Ignore,
            });
            if mode == TranslationMode::Cds {
                assert_eq!(ignored, Err(error.clone()));
            } else {
                assert_eq!(ignored.unwrap().ignored_trailing_bases(), input.len() % 3);
            }
        }
    }
    assert_eq!(
        DnaSequence::try_from("ATGAC")
            .unwrap()
            .translate_standard(TranslationOptions {
                partial_codons: PartialCodonPolicy::Ignore,
                ..TranslationOptions::default()
            })
            .unwrap()
            .as_str(),
        "M"
    );
}

#[test]
fn empty_sequences_are_valid_for_transforms_and_non_cds_translation() {
    let empty = DnaSequence::try_from("").unwrap();
    assert!(empty.complement().unwrap().is_empty());
    assert!(empty.reverse_complement().unwrap().is_empty());
    assert!(
        empty
            .transcribe()
            .unwrap()
            .back_transcribe()
            .unwrap()
            .is_empty()
    );
    for mode in [TranslationMode::Full, TranslationMode::ToStop] {
        let translation = empty.translate_standard(options(mode)).unwrap();
        assert!(translation.is_empty());
        assert_eq!(translation.len(), 0);
        assert_eq!(translation.ignored_trailing_bases(), 0);
        assert!(translation.to_protein().unwrap().is_empty());
    }
}

#[test]
fn invalid_input_is_rejected_before_transform_or_translation() {
    for input in ["ACGU", "AC-G", "AC G", "AC?G", "ACéG", "AC*G"] {
        assert!(DnaSequence::try_from(input).is_err(), "{input}");
    }
    assert!(RnaSequence::try_from("ACTG").is_err());
    assert_eq!(
        DnaSequence::try_from("acgtryswkmbdhvn")
            .unwrap()
            .to_string(),
        "ACGTRYSWKMBDHVN"
    );
    assert!(ProteinSequence::try_from("M*").is_err());
}

#[test]
fn conversion_never_silently_removes_stops() {
    let product = translate("ATGTAAGCT", TranslationMode::Full).unwrap();
    let error = product.to_protein().unwrap_err();
    assert_eq!(error.index(), 1);
    assert_eq!(error.byte(), b'*');
    assert_eq!(product.to_string(), "M*A");
    assert_eq!(product.len(), 3);
    assert_eq!(
        translate("RAYSARWTA", TranslationMode::Full)
            .unwrap()
            .to_protein()
            .unwrap()
            .to_string(),
        "BZJ"
    );
}

#[test]
fn error_messages_include_actionable_context() {
    for (error, expected) in [
        (
            TranslationError::PartialCodon {
                length: 4,
                remainder: 1,
            },
            "length 4",
        ),
        (TranslationError::CdsTooShort { length: 0 }, "length 0"),
        (
            TranslationError::InvalidStartCodon { codon: *b"GTG" },
            "GTG",
        ),
        (
            TranslationError::MissingTerminalStop { codon: *b"AAA" },
            "AAA",
        ),
        (TranslationError::InternalStop { offset: 6 }, "nucleotide 6"),
    ] {
        assert!(error.to_string().contains(expected));
        assert!(std::error::Error::source(&error).is_none());
    }
}

#[test]
fn real_phix174_endolysin_cds_matches_refseq_protein() {
    let records =
        bioferris::io::parse_fasta::<Dna>(include_bytes!("fixtures/phix174-e.fasta").as_slice())
            .unwrap();
    assert_eq!(records.len(), 1);
    let dna = records[0].sequence();
    assert_eq!(dna.len(), 276);
    let product = dna
        .translate_standard(options(TranslationMode::Cds))
        .unwrap();
    let expected = concat!(
        "MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNL",
        "RKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE"
    );
    assert_eq!(product.as_str(), expected);
    assert_eq!(product.len(), 91);
    assert_eq!(
        dna.transcribe()
            .unwrap()
            .translate_standard(options(TranslationMode::Cds))
            .unwrap(),
        product
    );
    let protein = product.to_protein().unwrap();
    let run = Pairwise::gotoh()
        .scores(2, -1)
        .gaps(-5, -1)
        .run(&protein, &protein)
        .unwrap();
    assert_eq!(run.result().first_range(), 0..91);
    assert_eq!(run.result().second_range(), 0..91);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn iupac_transform_invariants(input in "[ACGTRYSWKMBDHVN]{0,300}") {
        let dna = DnaSequence::try_from(input.as_str()).unwrap();
        let rna = dna.transcribe().unwrap();
        prop_assert_eq!(dna.complement().unwrap().complement().unwrap(), dna.clone());
        prop_assert_eq!(dna.reverse_complement().unwrap().reverse_complement().unwrap(), dna.clone());
        prop_assert_eq!(rna.back_transcribe().unwrap(), dna.clone());
        prop_assert_eq!(rna.complement().unwrap().complement().unwrap(), rna.clone());
        prop_assert_eq!(rna.reverse_complement().unwrap().reverse_complement().unwrap(), rna.clone());
        prop_assert_eq!(dna.reverse_complement().unwrap().transcribe().unwrap(), rna.reverse_complement().unwrap());
        prop_assert_eq!(dna.complement().unwrap().len(), input.len());
        prop_assert_eq!(dna.to_string(), input);
    }

    #[test]
    fn translation_is_deterministic_and_molecule_independent(mut input in "[ACGTRYSWKMBDHVN]{0,300}") {
        input.truncate(input.len() / 3 * 3);
        let dna = DnaSequence::try_from(input.as_str()).unwrap();
        let policy = TranslationOptions::default();
        let first = dna.translate_standard(policy).unwrap();
        prop_assert_eq!(first.len(), input.len() / 3);
        prop_assert_eq!(first.ignored_trailing_bases(), 0);
        prop_assert_eq!(&first, &dna.translate_standard(policy).unwrap());
        prop_assert_eq!(first, dna.transcribe().unwrap().translate_standard(policy).unwrap());
    }
}
