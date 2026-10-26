//! Independent NCBI codon data, segmented ORF oracle and coordinate checks.

use std::collections::BTreeSet;
use std::error::Error;
use std::io::Cursor;
use std::ops::Range;

use bioferris::io::{FastaReader, annotated::parse_embl, one_record};
use bioferris::prelude::*;
use proptest::prelude::*;

const CODES: [GeneticCode; 3] = [
    GeneticCode::Standard,
    GeneticCode::VertebrateMitochondrial,
    GeneticCode::BacterialArchaealPlastid,
];
const FRAMES: [ReadingFrame; 3] = [ReadingFrame::Zero, ReadingFrame::One, ReadingFrame::Two];
const STARTS: [OrfStartPolicy; 3] = [
    OrfStartPolicy::Atg,
    OrfStartPolicy::GeneticCode,
    OrfStartPolicy::AnySense,
];
const STOPS: [OrfStopPolicy; 2] = [OrfStopPolicy::Require, OrfStopPolicy::AllowPartial];
const STRANDS: [StrandSelection; 3] = [
    StrandSelection::Both,
    StrandSelection::Forward,
    StrandSelection::Reverse,
];

// Public biological data, not reference implementation source. NCBI genetic
// codes, revision 2024-09-23; provenance is recorded in fixtures/README.md.
// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// Grouped literal codons deliberately differ from the engine's indexed tables.
fn reference_residue(codon: &[u8], code: GeneticCode) -> u8 {
    if code == GeneticCode::VertebrateMitochondrial {
        match codon {
            b"ATA" => return b'M',
            b"TGA" => return b'W',
            b"AGA" | b"AGG" => return b'*',
            _ => {}
        }
    }
    match codon {
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"AAT" | b"AAC" => b'N',
        b"GAT" | b"GAC" => b'D',
        b"TGT" | b"TGC" => b'C',
        b"CAA" | b"CAG" => b'Q',
        b"GAA" | b"GAG" => b'E',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        b"CAT" | b"CAC" => b'H',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"AAA" | b"AAG" => b'K',
        b"ATG" => b'M',
        b"TTT" | b"TTC" => b'F',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"TGG" => b'W',
        b"TAT" | b"TAC" => b'Y',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        _ => panic!("not a concrete DNA codon: {codon:?}"),
    }
}

fn reference_start(codon: &[u8], code: GeneticCode, policy: OrfStartPolicy) -> bool {
    match policy {
        OrfStartPolicy::Atg => codon == b"ATG",
        OrfStartPolicy::AnySense => reference_residue(codon, code) != b'*',
        OrfStartPolicy::GeneticCode => match code {
            GeneticCode::Standard => matches!(codon, b"TTG" | b"CTG" | b"ATG"),
            GeneticCode::VertebrateMitochondrial => {
                matches!(codon, b"ATT" | b"ATC" | b"ATA" | b"ATG" | b"GTG")
            }
            GeneticCode::BacterialArchaealPlastid => matches!(
                codon,
                b"TTG" | b"CTG" | b"ATT" | b"ATC" | b"ATA" | b"ATG" | b"GTG"
            ),
            // The other tables are checked against NCBI gc.prt in ncbi_genetic_codes.rs.
            other => panic!("no independent initiation set for table {}", other.id()),
        },
    }
}

fn base_set(base: u8) -> &'static [u8] {
    match base {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => panic!("not a DNA IUPAC symbol: {base}"),
    }
}

fn reference_iupac_residue(codon: &[u8], code: GeneticCode) -> char {
    let mut outcomes = BTreeSet::new();
    for &a in base_set(codon[0]) {
        for &b in base_set(codon[1]) {
            for &c in base_set(codon[2]) {
                outcomes.insert(reference_residue(&[a, b, c], code));
            }
        }
    }
    let outcomes: Vec<_> = outcomes.into_iter().collect();
    char::from(match outcomes.as_slice() {
        [residue] => *residue,
        [b'D', b'N'] => b'B',
        [b'E', b'Q'] => b'Z',
        [b'I', b'L'] => b'J',
        _ => b'X',
    })
}

fn reference_reverse_complement(input: &[u8]) -> Vec<u8> {
    input
        .iter()
        .rev()
        .map(|base| {
            let position = b"ACGTRYSWKMBDHVN".iter().position(|x| x == base).unwrap();
            b"TGCAYRSWMKVHDBN"[position]
        })
        .collect()
}

fn original_range(range: Range<usize>, length: usize, strand: Strand) -> Range<usize> {
    match strand {
        Strand::Forward => range,
        Strand::Reverse => length - range.end..length - range.start,
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ExpectedOrf {
    strand: Strand,
    frame: ReadingFrame,
    code: GeneticCode,
    range: Range<usize>,
    stop: Option<Range<usize>>,
    peptide: String,
}

impl From<&OpenReadingFrame> for ExpectedOrf {
    fn from(orf: &OpenReadingFrame) -> Self {
        Self {
            strand: orf.strand(),
            frame: orf.reading_frame(),
            code: orf.genetic_code(),
            range: orf.range(),
            stop: orf.stop_range(),
            peptide: orf.peptide().to_owned(),
        }
    }
}

// Construct all complete codons, partition by precomputed stop indices, then
// select a suffix of each segment. This is independent of the production
// one-pass candidate state machine and never calls a translation method.
fn reference_orfs(input: &[u8], options: OrfOptions) -> Vec<ExpectedOrf> {
    let mut result = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse] {
        if matches!(
            (options.strands, strand),
            (StrandSelection::Forward, Strand::Reverse)
                | (StrandSelection::Reverse, Strand::Forward)
        ) {
            continue;
        }
        let oriented = match strand {
            Strand::Forward => input.to_vec(),
            Strand::Reverse => reference_reverse_complement(input),
        };
        for frame in FRAMES {
            let start = frame.offset().min(oriented.len());
            let codons: Vec<_> = oriented[start..].chunks_exact(3).collect();
            let ends: Vec<_> = codons
                .iter()
                .enumerate()
                .filter_map(|(i, codon)| {
                    (reference_residue(codon, options.genetic_code) == b'*').then_some(i)
                })
                .chain(std::iter::once(codons.len()))
                .collect();
            let mut segment_start = 0;
            for end in ends {
                let has_stop = end < codons.len();
                if !has_stop && options.stop_policy == OrfStopPolicy::Require {
                    break;
                }
                let segment = &codons[segment_start..end];
                let first = segment.iter().position(|codon| {
                    reference_start(codon, options.genetic_code, options.start_policy)
                });
                if let Some(first) = first {
                    let selected = &segment[first..];
                    if selected.len() >= options.min_amino_acids {
                        let mut residues: Vec<_> = selected
                            .iter()
                            .map(|codon| reference_residue(codon, options.genetic_code))
                            .collect();
                        if options.start_policy != OrfStartPolicy::AnySense {
                            residues[0] = b'M';
                        }
                        let coding = start + (segment_start + first) * 3..start + end * 3;
                        let stop = has_stop.then(|| start + end * 3..start + (end + 1) * 3);
                        result.push(ExpectedOrf {
                            strand,
                            frame,
                            code: options.genetic_code,
                            range: original_range(coding, input.len(), strand),
                            stop: stop.map(|range| original_range(range, input.len(), strand)),
                            peptide: String::from_utf8(residues).unwrap(),
                        });
                    }
                }
                segment_start = end + 1;
            }
        }
    }
    result
}

fn assert_oracle(input: &[u8], options: OrfOptions) -> Vec<OpenReadingFrame> {
    let expected = reference_orfs(input, options);
    let actual = dna(input).unwrap().find_orfs_with(options).unwrap();
    assert_eq!(
        actual.iter().map(ExpectedOrf::from).collect::<Vec<_>>(),
        expected,
        "{}, {options:?}",
        String::from_utf8_lossy(input)
    );
    actual
}

fn assert_six_frames(input: &[u8], code: GeneticCode) {
    let sequence = dna(input).unwrap();
    let dna_frames = sequence.translate_six_frames_with(code);
    let transcript = rna(input
        .iter()
        .map(|b| if *b == b'T' { b'U' } else { *b })
        .collect::<Vec<_>>())
    .unwrap();
    let rna_frames = transcript.translate_six_frames_with(code);
    for (index, (actual, rna_actual)) in dna_frames.iter().zip(&rna_frames).enumerate() {
        let strand = if index < 3 {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        let frame = FRAMES[index % 3];
        let oriented = if strand == Strand::Forward {
            input.to_vec()
        } else {
            reference_reverse_complement(input)
        };
        let leading = frame.offset().min(input.len());
        let remainder = (input.len() - leading) % 3;
        let end = input.len() - remainder;
        let expected: String = oriented[leading..end]
            .chunks_exact(3)
            .map(|codon| reference_iupac_residue(codon, code))
            .collect();
        for result in [actual, rna_actual] {
            assert_eq!(result.strand(), strand);
            assert_eq!(result.reading_frame(), frame);
            assert_eq!(result.genetic_code(), code);
            assert_eq!(
                result.range(),
                original_range(leading..end, input.len(), strand)
            );
            assert_eq!(result.skipped_leading_bases(), leading);
            assert_eq!(result.ignored_trailing_bases(), remainder);
            assert_eq!(result.peptide(), expected, "{input:?}, {code:?}, {index}");
        }
    }
    assert_eq!(sequence.to_string().as_bytes(), input);
}

#[test]
fn default_options_and_conveniences_select_the_documented_contract() {
    let options = OrfOptions::default();
    assert_eq!(options.genetic_code, GeneticCode::Standard);
    assert_eq!(options.strands, StrandSelection::Both);
    assert_eq!(options.start_policy, OrfStartPolicy::Atg);
    assert_eq!(options.stop_policy, OrfStopPolicy::Require);
    assert_eq!(options.min_amino_acids, 1);
    assert_eq!(options.max_orfs, 100_000);
    let sequence = dna("ATGAAATAA").unwrap();
    let found = sequence.find_orfs().unwrap();
    let configured = sequence.find_orfs_with(options).unwrap();
    assert_eq!(
        found.iter().map(ExpectedOrf::from).collect::<Vec<_>>(),
        configured.iter().map(ExpectedOrf::from).collect::<Vec<_>>()
    );
    for (simple, selected) in sequence.translate_six_frames().iter().zip(
        sequence
            .translate_six_frames_with(GeneticCode::Standard)
            .iter(),
    ) {
        assert_eq!(simple.peptide(), selected.peptide());
        assert_eq!(simple.range(), selected.range());
        assert_eq!(simple.strand(), selected.strand());
        assert_eq!(simple.reading_frame(), selected.reading_frame());
        assert_eq!(simple.genetic_code(), GeneticCode::Standard);
    }
}

#[test]
fn published_coordinate_examples_and_literal_six_frame_output_match() {
    let frames = dna("ATGAAATAA").unwrap().translate_six_frames();
    let expected = [
        (Strand::Forward, 0..9, "MK*", 0, 0),
        (Strand::Forward, 1..7, "*N", 1, 2),
        (Strand::Forward, 2..8, "EI", 2, 1),
        (Strand::Reverse, 0..9, "LFH", 0, 0),
        (Strand::Reverse, 2..8, "YF", 1, 2),
        (Strand::Reverse, 1..7, "IS", 2, 1),
    ];
    for (actual, (strand, range, peptide, leading, trailing)) in frames.iter().zip(expected) {
        assert_eq!(actual.strand(), strand);
        assert_eq!(actual.range(), range);
        assert_eq!(actual.peptide(), peptide);
        assert_eq!(actual.skipped_leading_bases(), leading);
        assert_eq!(actual.ignored_trailing_bases(), trailing);
    }
    for (input, strand, range, stop) in [
        ("ATGAAATAA", Strand::Forward, 0..6, 6..9),
        ("TTATTTCAT", Strand::Reverse, 3..9, 0..3),
    ] {
        let results = dna(input).unwrap().find_orfs().unwrap();
        assert_eq!(results.len(), 1);
        let actual = &results[0];
        assert_eq!(actual.strand(), strand);
        assert_eq!(actual.range(), range);
        assert_eq!(actual.stop_range(), Some(stop));
        assert_eq!(actual.peptide(), "MK");
        assert_eq!(actual.to_protein().unwrap().to_string(), "MK");
    }
}

#[test]
fn short_inputs_keep_six_frames_and_account_for_every_base() {
    for input in ["", "A", "AT", "ATG", "ATGA", "ATGAA", "ATGAAA"] {
        for code in CODES {
            assert_six_frames(input.as_bytes(), code);
        }
    }
    let frames = dna("A").unwrap().translate_six_frames();
    assert_eq!(frames[2].range(), 1..1);
    assert_eq!(frames[5].range(), 0..0);
    assert_eq!(frames[2].reading_frame(), ReadingFrame::Two);
    assert_eq!(frames[2].skipped_leading_bases(), 1);
    assert_eq!(frames[2].ignored_trailing_bases(), 0);
    assert_eq!(frames[0].ignored_trailing_bases(), 1);
}

#[test]
fn every_iupac_codon_translates_on_both_strands_in_all_supported_codes() {
    for &a in b"ACGTRYSWKMBDHVN" {
        for &b in b"ACGTRYSWKMBDHVN" {
            for &c in b"ACGTRYSWKMBDHVN" {
                for code in CODES {
                    assert_six_frames(&[a, b, c], code);
                }
            }
        }
    }
    let sequence = dna("RAYSARMTYTAN").unwrap();
    assert_eq!(sequence.translate_six_frames()[0].peptide(), "BZJX");
}

#[test]
fn every_concrete_start_and_terminal_assignment_matches_independent_tables() {
    for &a in b"ACGT" {
        for &b in b"ACGT" {
            for &c in b"ACGT" {
                for terminal in [b"TAA", b"TAG", b"TGA", b"AGA", b"AGG"] {
                    let input = [&[a, b, c][..], b"ATG", &terminal[..]].concat();
                    let transcript = rna(input
                        .iter()
                        .map(|b| if *b == b'T' { b'U' } else { *b })
                        .collect::<Vec<_>>())
                    .unwrap();
                    for genetic_code in CODES {
                        for start_policy in STARTS {
                            for stop_policy in STOPS {
                                let options = OrfOptions {
                                    genetic_code,
                                    start_policy,
                                    stop_policy,
                                    ..OrfOptions::default()
                                };
                                let expected = assert_oracle(&input, options);
                                let actual = transcript.find_orfs_with(options).unwrap();
                                assert_eq!(
                                    actual.iter().map(ExpectedOrf::from).collect::<Vec<_>>(),
                                    expected.iter().map(ExpectedOrf::from).collect::<Vec<_>>()
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn exhaustive_short_sequences_cover_codes_policies_strands_and_minimums() {
    // 5,461 sequences of lengths 0..=6, including every stop/start pair.
    // All 108 policy/code/strand/minimum combinations are checked per input.
    for length in 0..=6_u32 {
        for mut encoded in 0..4_usize.pow(length) {
            let mut input = vec![b'A'; usize::try_from(length).unwrap()];
            for base in &mut input {
                *base = b"ACGT"[encoded % 4];
                encoded /= 4;
            }
            for genetic_code in CODES {
                for start_policy in STARTS {
                    for stop_policy in STOPS {
                        for strands in STRANDS {
                            for min_amino_acids in [1, 2] {
                                assert_oracle(
                                    &input,
                                    OrfOptions {
                                        genetic_code,
                                        strands,
                                        start_policy,
                                        stop_policy,
                                        min_amino_acids,
                                        ..OrfOptions::default()
                                    },
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn first_start_wins_without_nested_candidates_and_stops_reset_the_interval() {
    let input = b"AAAATGATGAAATAAATGCCCTAG";
    let options = OrfOptions {
        strands: StrandSelection::Forward,
        ..OrfOptions::default()
    };
    let results = assert_oracle(input, options);
    let zero: Vec<_> = results
        .iter()
        .filter(|orf| orf.reading_frame() == ReadingFrame::Zero)
        .collect();
    assert_eq!(zero.len(), 2);
    assert_eq!(
        (zero[0].range(), zero[0].stop_range(), zero[0].peptide()),
        (3..12, Some(12..15), "MMK")
    );
    assert_eq!(
        (zero[1].range(), zero[1].stop_range(), zero[1].peptide()),
        (15..21, Some(21..24), "MP")
    );
    let reverse = reference_reverse_complement(input);
    let reversed = assert_oracle(
        &reverse,
        OrfOptions {
            strands: StrandSelection::Reverse,
            ..options
        },
    );
    let zero: Vec<_> = reversed
        .iter()
        .filter(|orf| orf.reading_frame() == ReadingFrame::Zero)
        .collect();
    assert_eq!(zero[0].range(), 12..21);
    assert_eq!(zero[0].stop_range(), Some(9..12));
    assert_eq!(zero[1].range(), 3..9);
    assert_eq!(zero[1].stop_range(), Some(0..3));
}

#[test]
fn start_modes_recode_only_the_initiator_and_any_sense_never_recodes() {
    let sequence = dna("GTGGTGTAA").unwrap();
    let options = OrfOptions {
        strands: StrandSelection::Forward,
        genetic_code: GeneticCode::BacterialArchaealPlastid,
        ..OrfOptions::default()
    };
    assert!(sequence.find_orfs_with(options).unwrap().is_empty());
    for (start_policy, expected) in [
        (OrfStartPolicy::GeneticCode, "MV"),
        (OrfStartPolicy::AnySense, "VV"),
    ] {
        let results = sequence
            .find_orfs_with(OrfOptions {
                start_policy,
                ..options
            })
            .unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].peptide(), expected);
    }
    assert_eq!(
        sequence.translate_six_frames_with(options.genetic_code)[0].peptide(),
        "VV*"
    );
    for (code, expected) in [
        (GeneticCode::Standard, "M"),
        (GeneticCode::VertebrateMitochondrial, "MW"),
        (GeneticCode::BacterialArchaealPlastid, "M"),
    ] {
        let results = dna("ATGTGATAA")
            .unwrap()
            .find_orfs_with(OrfOptions {
                genetic_code: code,
                ..options
            })
            .unwrap();
        assert_eq!(results[0].peptide(), expected);
    }
}

#[test]
fn partial_policy_keeps_only_complete_codons_and_does_not_infer_a_start() {
    let options = OrfOptions {
        strands: StrandSelection::Forward,
        stop_policy: OrfStopPolicy::AllowPartial,
        ..OrfOptions::default()
    };
    for (input, end, peptide) in [
        ("ATG", 3, "M"),
        ("ATGA", 3, "M"),
        ("ATGAA", 3, "M"),
        ("ATGAAA", 6, "MK"),
        ("ATGAAAT", 6, "MK"),
        ("ATGAAATT", 6, "MK"),
    ] {
        let results = assert_oracle(input.as_bytes(), options);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].range(), 0..end);
        assert_eq!(results[0].stop_range(), None);
        assert_eq!(results[0].peptide(), peptide);
        assert!(
            dna(input)
                .unwrap()
                .find_orfs_with(OrfOptions {
                    stop_policy: OrfStopPolicy::Require,
                    ..options
                })
                .unwrap()
                .is_empty()
        );
    }
    for input in ["", "A", "AT", "AAA", "AAATAA", "TAATAATAA"] {
        assert!(
            dna(input)
                .unwrap()
                .find_orfs_with(options)
                .unwrap()
                .is_empty()
        );
    }
    let results = assert_oracle(b"ATGTAAATG", options);
    assert_eq!(results.len(), 2);
    assert_eq!(results[0].stop_range(), Some(3..6));
    assert_eq!(results[1].stop_range(), None);
    assert_eq!(results[1].range(), 6..9);
}

#[test]
fn ambiguity_is_rejected_at_the_first_original_offset_for_every_policy() {
    for &ambiguous in b"RYSWKMBDHVN" {
        for position in [0, 1, 8, 9, 10] {
            let mut input = b"CATGAAATAAC".to_vec();
            input[position] = ambiguous;
            let transcript: Vec<_> = input
                .iter()
                .map(|b| if *b == b'T' { b'U' } else { *b })
                .collect();
            for strands in STRANDS {
                for start_policy in STARTS {
                    for stop_policy in STOPS {
                        let options = OrfOptions {
                            strands,
                            start_policy,
                            stop_policy,
                            ..OrfOptions::default()
                        };
                        for error in [
                            dna(&input).unwrap().find_orfs_with(options).unwrap_err(),
                            rna(&transcript)
                                .unwrap()
                                .find_orfs_with(options)
                                .unwrap_err(),
                        ] {
                            assert_eq!(
                                error,
                                OrfError::AmbiguousBase {
                                    offset: position,
                                    byte: ambiguous
                                }
                            );
                        }
                    }
                }
            }
        }
    }
    assert_eq!(
        dna("NATGTAAR").unwrap().find_orfs().unwrap_err(),
        OrfError::AmbiguousBase {
            offset: 0,
            byte: b'N'
        }
    );
    assert_eq!(
        dna("aN").unwrap().find_orfs().unwrap_err(),
        OrfError::AmbiguousBase {
            offset: 1,
            byte: b'N'
        }
    );
    assert_eq!(dna("ATG-").unwrap_err().index(), 3);
    assert_eq!(rna("AUGT").unwrap_err().index(), 3);
}

#[test]
fn options_are_validated_before_ambiguity_and_errors_implement_error() {
    let error = dna("N")
        .unwrap()
        .find_orfs_with(OrfOptions {
            min_amino_acids: 0,
            max_orfs: 0,
            ..OrfOptions::default()
        })
        .unwrap_err();
    assert_eq!(error, OrfError::InvalidMinimumLength);
    assert_eq!(
        dna("ATGAAATAAN")
            .unwrap()
            .find_orfs_with(OrfOptions {
                max_orfs: 0,
                ..OrfOptions::default()
            })
            .unwrap_err(),
        OrfError::AmbiguousBase {
            offset: 9,
            byte: b'N'
        }
    );
    for error in [
        error,
        OrfError::AmbiguousBase {
            offset: 1,
            byte: b'R',
        },
        OrfError::ResultLimitExceeded { limit: 3 },
    ] {
        assert!(error.source().is_none());
        assert!(!error.to_string().is_empty());
    }
}

#[test]
fn result_limits_are_inclusive_and_apply_after_minimum_length_filtering() {
    let options = OrfOptions {
        strands: StrandSelection::Forward,
        ..OrfOptions::default()
    };
    let sequence = dna("ATGAAATAAATGTAA").unwrap();
    assert_eq!(
        sequence
            .find_orfs_with(OrfOptions {
                max_orfs: 0,
                ..options
            })
            .unwrap_err(),
        OrfError::ResultLimitExceeded { limit: 0 }
    );
    assert_eq!(
        sequence
            .find_orfs_with(OrfOptions {
                max_orfs: 1,
                ..options
            })
            .unwrap_err(),
        OrfError::ResultLimitExceeded { limit: 1 }
    );
    assert_eq!(
        sequence
            .find_orfs_with(OrfOptions {
                max_orfs: 2,
                ..options
            })
            .unwrap()
            .len(),
        2
    );
    let filtered = sequence
        .find_orfs_with(OrfOptions {
            min_amino_acids: 2,
            max_orfs: 1,
            ..options
        })
        .unwrap();
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].peptide(), "MK");
    for input in ["", "A", "ATG", "AAA", "TAA"] {
        assert!(
            dna(input)
                .unwrap()
                .find_orfs_with(OrfOptions {
                    max_orfs: 0,
                    ..options
                })
                .unwrap()
                .is_empty()
        );
    }
    assert!(
        sequence
            .find_orfs_with(OrfOptions {
                min_amino_acids: usize::MAX,
                max_orfs: 0,
                ..options
            })
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        sequence
            .find_orfs_with(OrfOptions {
                max_orfs: usize::MAX,
                ..options
            })
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn global_result_budget_spans_both_strands_without_counting_nested_starts() {
    let sequence = dna("ATGAAATAATTATTTCAT").unwrap();
    let options = OrfOptions {
        max_orfs: 1,
        ..OrfOptions::default()
    };
    for strands in [StrandSelection::Forward, StrandSelection::Reverse] {
        assert_eq!(
            sequence
                .find_orfs_with(OrfOptions { strands, ..options })
                .unwrap()
                .len(),
            1
        );
    }
    assert_eq!(
        sequence.find_orfs_with(options).unwrap_err(),
        OrfError::ResultLimitExceeded { limit: 1 }
    );

    // A long interval with 10,000 accepted starts still emits only one ORF.
    // This checks semantics without making a timing or memory-performance claim.
    let nested = dna(format!("{}TAA", "ATG".repeat(10_000))).unwrap();
    let results = nested
        .find_orfs_with(OrfOptions {
            strands: StrandSelection::Forward,
            ..options
        })
        .unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].range(), 0..30_000);
    assert_eq!(results[0].stop_range(), Some(30_000..30_003));
    assert_eq!(results[0].peptide(), "M".repeat(10_000));
}

#[test]
fn public_phix174_record_retains_annotated_protein_on_both_orientations() {
    // NC_001422.1:568..843 (one-based inclusive); NP_040709.1 translation.
    // This local public-data fixture predates this increment; no campaign data.
    let record = one_record(FastaReader::<_, Dna>::new(Cursor::new(include_bytes!(
        "fixtures/phix174-e.fasta"
    ))))
    .unwrap();
    let expected = concat!(
        "MVRWTLWDTLAFLLLLSLLLPSLLIMFIPSTFKRPVSSWKALNL",
        "RKTLLMASSVRLKPLNCSRLPCVYAQETLTFLLTQKKTCVKNYVRKE"
    );
    let input = record.sequence().to_string();
    assert_eq!(input.len(), 276);
    for (input, strand, range, stop) in [
        (input.into_bytes(), Strand::Forward, 0..273, 273..276),
        (
            reference_reverse_complement(record.sequence().to_string().as_bytes()),
            Strand::Reverse,
            3..276,
            0..3,
        ),
    ] {
        let options = OrfOptions {
            genetic_code: GeneticCode::BacterialArchaealPlastid,
            min_amino_acids: 91,
            ..OrfOptions::default()
        };
        let candidates = assert_oracle(&input, options);
        let coding = candidates
            .iter()
            .find(|orf| orf.strand() == strand && orf.range() == range)
            .unwrap();
        assert_eq!(coding.peptide(), expected);
        assert_eq!(coding.stop_range(), Some(stop));
        assert_eq!(coding.to_protein().unwrap().to_string(), expected);
        let frames = dna(&input)
            .unwrap()
            .translate_six_frames_with(options.genetic_code);
        let selected = frames
            .iter()
            .find(|frame| frame.strand() == strand && frame.reading_frame() == ReadingFrame::Zero)
            .unwrap();
        assert_eq!(selected.peptide(), format!("{expected}*"));
    }
}

#[test]
fn public_puc19_embl_sequence_uses_the_same_explicit_linear_scan_contract() {
    // ENA M77789.2, 2,686-base pUC19; provenance is in annotated-provenance.md.
    // Circular record metadata does not imply circular ORF wrapping.
    let records = parse_embl::<Dna>(include_str!("fixtures/annotated-m77789.2.embl")).unwrap();
    let record = &records[0];
    assert_eq!(record.record().id(), "M77789.2");
    let sequence = record.record().sequence();
    assert_eq!(sequence.len(), 2686);
    for genetic_code in CODES {
        assert_six_frames(sequence.to_string().as_bytes(), genetic_code);
        for start_policy in STARTS {
            for stop_policy in STOPS {
                let results = assert_oracle(
                    sequence.to_string().as_bytes(),
                    OrfOptions {
                        genetic_code,
                        start_policy,
                        stop_policy,
                        min_amino_acids: 20,
                        ..OrfOptions::default()
                    },
                );
                assert!(!results.is_empty());
                assert!(results.iter().all(|orf| orf.range().end <= sequence.len()));
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn arbitrary_iupac_frames_match_independent_translation_and_existing_explicit_composition(
        input in "[ACGTRYSWKMBDHVN]{0,90}",
        code_index in 0_usize..3,
    ) {
        let code = CODES[code_index];
        assert_six_frames(input.as_bytes(), code);
        let sequence = dna(&input).unwrap();
        let reverse = sequence.reverse_complement().unwrap();
        for frame in sequence.translate_six_frames_with(code) {
            let oriented = if frame.strand() == Strand::Forward { &sequence } else { &reverse };
            if frame.reading_frame().offset() <= input.len() {
                let composed = oriented.translate_with(TranslationConfig {
                    genetic_code: code,
                    reading_frame: frame.reading_frame(),
                    partial_codons: PartialCodonPolicy::Ignore,
                    ..TranslationConfig::default()
                }).unwrap();
                prop_assert_eq!(frame.peptide(), composed.as_str());
                prop_assert_eq!(frame.ignored_trailing_bases(), composed.ignored_trailing_bases());
            }
        }
    }

    #[test]
    fn arbitrary_orfs_match_segmented_oracle_rna_coordinates_and_deterministic_order(
        input in "[ACGT]{0,180}",
        code_index in 0_usize..3,
        start_index in 0_usize..3,
        stop_index in 0_usize..2,
        strand_index in 0_usize..3,
        minimum in 1_usize..12,
    ) {
        let options = OrfOptions {
            genetic_code: CODES[code_index],
            start_policy: STARTS[start_index],
            stop_policy: STOPS[stop_index],
            strands: STRANDS[strand_index],
            min_amino_acids: minimum,
            ..OrfOptions::default()
        };
        let candidates = assert_oracle(input.as_bytes(), options);
        let transcript = rna(input.replace('T', "U")).unwrap().find_orfs_with(options).unwrap();
        let repeated = dna(input.to_ascii_lowercase()).unwrap().find_orfs_with(options).unwrap();
        let actual: Vec<_> = candidates.iter().map(ExpectedOrf::from).collect();
        prop_assert_eq!(&actual, &transcript.iter().map(ExpectedOrf::from).collect::<Vec<_>>());
        prop_assert_eq!(&actual, &repeated.iter().map(ExpectedOrf::from).collect::<Vec<_>>());
        let mut previous_order = None;
        for candidate in &candidates {
            let range = candidate.range();
            prop_assert_eq!(range.len(), candidate.peptide().len() * 3);
            prop_assert!(candidate.peptide().len() >= minimum);
            prop_assert!(!candidate.peptide().contains('*'));
            prop_assert_eq!(candidate.to_protein().unwrap().to_string(), candidate.peptide());
            let (oriented, oriented_start, strand_order) = if candidate.strand() == Strand::Forward {
                (input.as_bytes()[range.clone()].to_vec(), range.start, 0)
            } else {
                (reference_reverse_complement(&input.as_bytes()[range.clone()]), input.len() - range.end, 1)
            };
            prop_assert_eq!(oriented_start % 3, candidate.reading_frame().offset());
            let mut expected: Vec<_> = oriented.chunks_exact(3).map(|codon| reference_residue(codon, options.genetic_code)).collect();
            if options.start_policy != OrfStartPolicy::AnySense { expected[0] = b'M'; }
            prop_assert_eq!(String::from_utf8(expected).unwrap(), candidate.peptide());
            if let Some(stop) = candidate.stop_range() {
                prop_assert_eq!(stop.len(), 3);
                let stop_codon = if candidate.strand() == Strand::Forward {
                    prop_assert_eq!(stop.start, range.end);
                    input.as_bytes()[stop].to_vec()
                } else {
                    prop_assert_eq!(stop.end, range.start);
                    reference_reverse_complement(&input.as_bytes()[stop])
                };
                prop_assert_eq!(reference_residue(&stop_codon, options.genetic_code), b'*');
            } else {
                prop_assert_eq!(options.stop_policy, OrfStopPolicy::AllowPartial);
            }
            let order = (strand_order, candidate.reading_frame().offset(), oriented_start);
            if let Some(previous) = previous_order { prop_assert!(previous < order); }
            previous_order = Some(order);
        }
    }

    #[test]
    fn minimum_length_filters_without_reordering_and_limits_never_truncate(
        input in "[ACGT]{0,150}",
        minimum in 1_usize..10,
        limit in 0_usize..12,
        code_index in 0_usize..3,
    ) {
        let sequence = dna(input).unwrap();
        let options = OrfOptions {
            genetic_code: CODES[code_index],
            start_policy: OrfStartPolicy::AnySense,
            stop_policy: OrfStopPolicy::AllowPartial,
            ..OrfOptions::default()
        };
        let all = sequence.find_orfs_with(options).unwrap();
        let expected: Vec<_> = all.iter().filter(|orf| orf.peptide().len() >= minimum).map(ExpectedOrf::from).collect();
        let found = sequence.find_orfs_with(OrfOptions { min_amino_acids: minimum, ..options }).unwrap();
        prop_assert_eq!(&expected, &found.iter().map(ExpectedOrf::from).collect::<Vec<_>>());
        let limited = sequence.find_orfs_with(OrfOptions { min_amino_acids: minimum, max_orfs: limit, ..options });
        if expected.len() <= limit {
            prop_assert_eq!(expected, limited.unwrap().iter().map(ExpectedOrf::from).collect::<Vec<_>>());
        } else {
            prop_assert_eq!(limited.unwrap_err(), OrfError::ResultLimitExceeded { limit });
        }
    }
}

#[test]
fn orf_discovery_refuses_codes_without_an_unconditional_stop() {
    let sequence = dna("ATGAAATGAAAATAGAAATAA").unwrap();
    for code in GeneticCode::ALL {
        let result = sequence.find_orfs_with(OrfOptions {
            genetic_code: code,
            ..OrfOptions::default()
        });
        if [27, 28, 31].contains(&code.id()) {
            assert_eq!(
                result,
                Err(OrfError::ConditionalStops { genetic_code: code })
            );
        } else {
            assert!(result.is_ok(), "table {}", code.id());
        }
    }
    assert_eq!(
        OrfError::ConditionalStops {
            genetic_code: GeneticCode::KaryorelictNuclear
        }
        .to_string(),
        "genetic code 27 has no codon that always ends translation, so ORF boundaries are undefined"
    );
}
