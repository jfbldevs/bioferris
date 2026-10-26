//! Independent exact-search oracles, range/record invariants and worked examples.

use std::error::Error;
use std::io::Cursor;
use std::iter::FusedIterator;
use std::path::Path;
use std::process::Command;

use bioferris::io::{FastaReader, FastqReader, IoErrorKind, one_record, write_record_fastq};
use bioferris::prelude::*;
use proptest::prelude::*;

// Deliberately naive byte comparison, independent of the production search
// preprocessing and iterator state. All inputs here are validated ASCII.
fn naive_matches(text: &[u8], pattern: &[u8], overlapping: bool) -> Vec<usize> {
    if pattern.is_empty() {
        return (0..=text.len()).collect();
    }
    let Some(last_start) = text.len().checked_sub(pattern.len()) else {
        return Vec::new();
    };
    let mut cursor = 0;
    let mut matches = Vec::new();
    while cursor <= last_start {
        if text[cursor..cursor + pattern.len()] == *pattern {
            matches.push(cursor);
            cursor += if overlapping { 1 } else { pattern.len() };
        } else {
            cursor += 1;
        }
    }
    matches
}

fn require_fused<T: FusedIterator<Item = usize>>(_: &T) {}

fn check_search<A: Alphabet>(input: &str, pattern: &str) {
    let sequence = Sequence::<A>::from_ascii(input).unwrap();
    let normalized = input.to_ascii_uppercase();
    let normalized_pattern = pattern.to_ascii_uppercase();
    let non_overlapping =
        naive_matches(normalized.as_bytes(), normalized_pattern.as_bytes(), false);
    let overlapping = naive_matches(normalized.as_bytes(), normalized_pattern.as_bytes(), true);

    assert_eq!(
        sequence.find(pattern).unwrap(),
        overlapping.first().copied()
    );
    assert_eq!(sequence.contains(pattern).unwrap(), !overlapping.is_empty());
    assert_eq!(sequence.count(pattern).unwrap(), non_overlapping.len());
    assert_eq!(
        sequence.count_overlapping(pattern).unwrap(),
        overlapping.len()
    );
    for (mode, expected) in [
        (MatchOverlap::NonOverlapping, non_overlapping),
        (MatchOverlap::Overlapping, overlapping),
    ] {
        let mut matches = sequence.find_all(pattern, mode).unwrap();
        require_fused(&matches);
        assert_eq!(matches.by_ref().collect::<Vec<_>>(), expected);
        assert_eq!(matches.next(), None);
        assert_eq!(matches.next(), None);
    }
    assert_eq!(sequence.to_string(), normalized);
}

#[test]
fn documented_counts_distinguish_overlapping_from_non_overlapping() {
    let sequence = dna("AAAA").unwrap();
    assert_eq!(sequence.find("AA").unwrap(), Some(0));
    assert_eq!(sequence.count("AA").unwrap(), 2);
    assert_eq!(sequence.count_overlapping("AA").unwrap(), 3);
    for (input, pattern) in [
        ("AAAA", "AA"),
        ("ATATAT", "ATA"),
        ("ATATACATATAC", "ATATAC"),
        ("ACACACAG", "ACACAG"),
        ("ACACACAC", "ACACAG"),
        ("ACGT", "GT"),
        ("ACGT", "ACGTA"),
        ("A", "A"),
        ("A", "C"),
        ("", "A"),
    ] {
        check_search::<Dna>(input, pattern);
    }
}

#[test]
fn ambiguity_is_literal_and_patterns_normalize_case() {
    let sequence = dna("acgtnryn").unwrap();
    assert_eq!(sequence.find("n").unwrap(), Some(4));
    assert_eq!(sequence.find("ry").unwrap(), Some(5));
    assert_eq!(sequence.count_overlapping("N").unwrap(), 2);
    assert_eq!(sequence.find("NN").unwrap(), None);
    assert_eq!(dna("ACGT").unwrap().find("N").unwrap(), None);
    assert_eq!(protein("ILNQ").unwrap().find("J").unwrap(), None);
    assert_eq!(protein("MABZXJUO").unwrap().find("zxj").unwrap(), Some(3));
    check_search::<Rna>("augunnry", "un");
    check_search::<Protein>("MAJXMAJX", "majx");
}

#[test]
fn empty_patterns_match_every_boundary_including_empty_sequences() {
    for input in ["", "A", "ACGTN"] {
        check_search::<Dna>(input, "");
        let sequence = dna(input).unwrap();
        assert_eq!(sequence.find("").unwrap(), Some(0));
        assert!(sequence.contains("").unwrap());
        assert_eq!(sequence.count("").unwrap(), input.len() + 1);
        assert_eq!(sequence.count_overlapping("").unwrap(), input.len() + 1);
    }
    check_search::<Rna>("", "");
    check_search::<Protein>("", "");
}

fn verify_invalid_pattern<A: Alphabet>(pattern: &[u8], invalid_index: usize, byte: u8) {
    // Validation must happen even for an empty source that cannot contain the
    // nonempty pattern: invalid input must not be reported as an ordinary miss.
    let sequence = Sequence::<A>::from_ascii("").unwrap();
    let results = [
        sequence.find(pattern).map(|_| ()),
        sequence.contains(pattern).map(|_| ()),
        sequence.count(pattern).map(|_| ()),
        sequence.count_overlapping(pattern).map(|_| ()),
        sequence
            .find_all(pattern, MatchOverlap::Overlapping)
            .map(|_| ()),
        sequence
            .find_all(pattern, MatchOverlap::NonOverlapping)
            .map(|_| ()),
    ];
    for result in results {
        let error = result.unwrap_err();
        assert_eq!(error.code(), "invalid_sequence");
        assert!(error.source().is_some());
        let SequenceOperationError::InvalidSequence(source) = error else {
            panic!("expected the invalid-pattern sequence error");
        };
        assert_eq!(source.alphabet(), A::NAME);
        assert_eq!(source.index(), invalid_index);
        assert_eq!(source.byte(), byte);
    }
}

#[test]
fn every_search_method_validates_the_entire_pattern_and_reports_its_offset() {
    for byte in [b'U', b'X', b'-', b' ', b'\n', 0, 0xff] {
        verify_invalid_pattern::<Dna>(&[b'a', byte], 1, byte);
    }
    verify_invalid_pattern::<Rna>(b"auT", 2, b'T');
    verify_invalid_pattern::<Protein>(b"MA*", 2, b'*');
    verify_invalid_pattern::<Dna>("Aé".as_bytes(), 1, 0xc3);
}

#[test]
fn all_byte_patterns_follow_the_declared_alphabet_without_wildcards() {
    fn check<A: Alphabet>(symbols: &[u8]) {
        let sequence = Sequence::<A>::from_ascii(symbols).unwrap();
        for byte in 0..=u8::MAX {
            let expected = symbols
                .iter()
                .position(|symbol| *symbol == byte.to_ascii_uppercase());
            match expected {
                Some(index) => assert_eq!(sequence.find([byte]).unwrap(), Some(index)),
                None => verify_invalid_pattern::<A>(&[byte], 0, byte),
            }
        }
    }
    check::<Dna>(b"ACGTRYSWKMBDHVN");
    check::<Rna>(b"ACGURYSWKMBDHVN");
    check::<Protein>(b"ARNDCQEGHILKMFPSTWYVBZXJUO");
}

fn all_words(alphabet: &[u8], maximum_length: usize) -> Vec<String> {
    let mut words = vec![String::new()];
    let mut layer = vec![String::new()];
    for _ in 0..maximum_length {
        layer = layer
            .iter()
            .flat_map(|prefix| {
                alphabet.iter().map(move |byte| {
                    let mut word = prefix.clone();
                    word.push(char::from(*byte));
                    word
                })
            })
            .collect();
        words.extend(layer.iter().cloned());
    }
    words
}

#[test]
fn exhaustive_short_searches_match_an_independent_naive_oracle() {
    let sources = all_words(b"ACN", 5);
    let patterns = all_words(b"ACN", 4);
    assert_eq!((sources.len(), patterns.len()), (364, 121));
    for source in &sources {
        for pattern in &patterns {
            check_search::<Dna>(source, pattern);
        }
    }
}

#[test]
fn match_iterators_own_the_pattern_and_support_independent_partial_iteration() {
    let sequence = dna("ATATAT").unwrap();
    let mut matches = {
        let temporary_pattern = String::from("ata");
        sequence
            .find_all(temporary_pattern, MatchOverlap::Overlapping)
            .unwrap()
    };
    assert_eq!(matches.next(), Some(0));
    assert_eq!(sequence.count("ATA").unwrap(), 1);
    assert_eq!(matches.next(), Some(2));
    assert_eq!(matches.next(), None);
    assert_eq!(matches.next(), None);
}

#[test]
fn repeated_prefixes_and_lazy_partial_results_cover_longer_search_state() {
    let input = format!("{}C{}C", "A".repeat(10_000), "A".repeat(10_000));
    let pattern = format!("{}C", "A".repeat(500));
    let sequence = dna(&input).unwrap();
    assert_eq!(
        sequence
            .find_all(&pattern, MatchOverlap::Overlapping)
            .unwrap()
            .collect::<Vec<_>>(),
        [9_500, 19_501]
    );
    assert_eq!(sequence.count(&pattern).unwrap(), 2);
    assert_eq!(
        sequence
            .find_all("A", MatchOverlap::Overlapping)
            .unwrap()
            .take(3)
            .collect::<Vec<_>>(),
        [0, 1, 2]
    );
}

fn check_every_slice<A: Alphabet>(input: &str) {
    let original = Sequence::<A>::from_ascii(input).unwrap();
    let normalized = input.to_ascii_uppercase();
    for start in 0..=input.len() {
        for end in start..=input.len() {
            let slice = original.slice(start..end).unwrap();
            assert_eq!(slice.to_string(), normalized[start..end]);
            assert_eq!(slice.len(), end - start);
            let prefix = original.slice(0..start).unwrap();
            let suffix = original.slice(end..input.len()).unwrap();
            assert_eq!(
                prefix.concat(&slice).unwrap().concat(&suffix).unwrap(),
                original
            );
        }
    }
    assert_eq!(original.to_string(), normalized);
}

#[test]
fn slices_accept_exact_empty_full_and_terminal_ranges_for_each_alphabet() {
    check_every_slice::<Dna>("aCGtn");
    check_every_slice::<Rna>("acGun");
    check_every_slice::<Protein>("MaJxO");
    check_every_slice::<Dna>("");
    check_every_slice::<Rna>("");
    check_every_slice::<Protein>("");
}

#[test]
fn invalid_ranges_report_the_original_coordinates_without_clipping_or_panicking() {
    let sequence = dna("ACGT").unwrap();
    for (start, end) in [(3, 2), (0, 5), (5, 5), (usize::MAX, 0), (0, usize::MAX)] {
        let error = sequence.slice(start..end).unwrap_err();
        assert_eq!(error.code(), "invalid_range");
        assert!(error.source().is_none());
        assert!(matches!(
            error,
            SequenceOperationError::InvalidRange { start: actual_start, end: actual_end, length: 4 }
                if actual_start == start && actual_end == end
        ));
    }
    assert!(dna("").unwrap().slice(0..1).is_err());
    assert_eq!(sequence.to_string(), "ACGT");
}

#[test]
fn concatenation_has_identity_and_associativity_and_does_not_mutate_sources() {
    fn check<A: Alphabet>(first: &str, second: &str, third: &str) {
        let first = Sequence::<A>::from_ascii(first).unwrap();
        let second = Sequence::<A>::from_ascii(second).unwrap();
        let third = Sequence::<A>::from_ascii(third).unwrap();
        let empty = Sequence::<A>::from_ascii("").unwrap();
        let originals = (first.clone(), second.clone(), third.clone());
        assert_eq!(first.concat(&empty).unwrap(), first);
        assert_eq!(empty.concat(&first).unwrap(), first);
        assert_eq!(empty.concat(&empty).unwrap(), empty);
        assert_eq!(
            first.concat(&second).unwrap().concat(&third).unwrap(),
            first.concat(&second.concat(&third).unwrap()).unwrap()
        );
        assert_eq!((first, second, third), originals);
    }
    check::<Dna>("ac", "n", "t");
    check::<Rna>("ac", "n", "u");
    check::<Protein>("Ma", "JX", "UO");
}

#[test]
fn operation_errors_have_stable_codes_and_standard_error_display() {
    let overflow = SequenceOperationError::LengthOverflow;
    assert_eq!(overflow.code(), "length_overflow");
    assert!(overflow.source().is_none());
    assert!(!overflow.to_string().is_empty());
    let range = dna("AC").unwrap().slice(0..3).unwrap_err();
    assert!(!range.to_string().is_empty());
}

#[test]
fn record_slices_preserve_header_case_quality_and_the_original_record() {
    let original = Record::<Dna>::from_ascii(
        "  read-1\tMixed\tcase  ",
        "aCgTn",
        Some(vec![0, 1, 40, 62, 93]),
    )
    .unwrap();
    let snapshot = original.clone();
    for start in 0..=5 {
        for end in start..=5 {
            let region = original.slice(start..end).unwrap();
            assert_eq!(region.header(), original.header());
            assert_eq!(region.id(), "read-1");
            assert_eq!(region.description(), Some("Mixed\tcase"));
            assert_eq!(
                region.source_sequence(),
                &original.source_sequence()[start..end]
            );
            assert_eq!(region.sequence().to_string(), "ACGTN"[start..end]);
            assert_eq!(region.qualities(), Some(&[0, 1, 40, 62, 93][start..end]));
        }
    }
    assert_eq!(original, snapshot);
}

#[test]
fn records_preserve_absent_qualities_and_empty_present_qualities() {
    let rna_record = Record::<Rna>::from_ascii("rna", "aUg", None).unwrap();
    let rna_slice = rna_record.slice(1..3).unwrap();
    assert_eq!(rna_slice.source_sequence(), "Ug");
    assert_eq!(rna_slice.sequence().to_string(), "UG");
    assert_eq!(rna_slice.description(), None);
    assert_eq!(rna_slice.qualities(), None);

    let protein_record = Record::<Protein>::from_ascii("protein", "MaX", None).unwrap();
    assert_eq!(protein_record.slice(1..3).unwrap().source_sequence(), "aX");

    for qualities in [None, Some(vec![])] {
        let record = Record::<Dna>::from_ascii("empty", "", qualities.clone()).unwrap();
        let region = record.slice(0..0).unwrap();
        assert_eq!(region.qualities(), qualities.as_deref());
        assert_eq!(region, record);
        assert!(record.slice(0..1).is_err());
    }
}

#[test]
fn invalid_record_ranges_are_checked_before_copying_metadata_or_qualities() {
    let record = Record::<Dna>::from_ascii("r1", "aC", Some(vec![10, 20])).unwrap();
    let snapshot = record.clone();
    for (start, end) in [(1, 0), (0, 3), (3, 3), (usize::MAX, usize::MAX)] {
        let error = record.slice(start..end).unwrap_err();
        assert!(matches!(
            error,
            SequenceOperationError::InvalidRange { start: actual_start, end: actual_end, length: 2 }
                if actual_start == start && actual_end == end
        ));
    }
    assert_eq!(record, snapshot);
}

#[test]
fn real_phix174_record_slicing_integrates_with_existing_translation() {
    // Existing public RefSeq fixture; no new network fetch or altered provenance.
    let record = one_record(FastaReader::<_, Dna>::new(Cursor::new(include_bytes!(
        "fixtures/phix174-e.fasta"
    ))))
    .unwrap();
    let prefix = record.slice(0..12).unwrap();
    assert_eq!(prefix.id(), "NC_001422.1:568-843");
    assert_eq!(prefix.header(), record.header());
    assert_eq!(prefix.source_sequence(), "ATGGTACGCTGG");
    assert_eq!(
        prefix
            .sequence()
            .translate_standard(TranslationOptions::default())
            .unwrap()
            .as_str(),
        "MVRW"
    );
    let suffix = record.slice(12..record.sequence().len()).unwrap();
    assert_eq!(
        &prefix.sequence().concat(suffix.sequence()).unwrap(),
        record.sequence()
    );
    assert_eq!(record.sequence().find("ATGGTACGCTGG").unwrap(), Some(0));
    assert_eq!(record.sequence().len(), 276);
}

#[test]
fn fastq_regions_round_trip_with_matching_quality_and_reject_empty_writes() {
    let record = one_record(FastqReader::<_, Dna>::new(Cursor::new(include_bytes!(
        "fixtures/io-example.fastq"
    ))))
    .unwrap();
    let region = record.slice(1..3).unwrap();
    let mut output = Vec::new();
    write_record_fastq(&mut output, &region).unwrap();
    assert_eq!(output, b"@read info\nCg\n+\nI+\n");
    let recovered = one_record(FastqReader::<_, Dna>::new(Cursor::new(output))).unwrap();
    assert_eq!(recovered, region);
    assert_eq!(recovered.qualities(), Some([40, 10].as_slice()));

    let empty = record.slice(4..4).unwrap();
    let mut output = Vec::new();
    assert_eq!(
        write_record_fastq(&mut output, &empty).unwrap_err().kind(),
        IoErrorKind::Unsupported
    );
    assert!(output.is_empty());
    assert_eq!(record.source_sequence(), "aCgT");
    assert_eq!(record.qualities(), Some([0, 40, 10, 93].as_slice()));
}

fn check_generated_operations<A: Alphabet>(
    input: &str,
    pattern: &str,
    first: usize,
    second: usize,
) {
    check_search::<A>(input, pattern);
    let first = first % (input.len() + 1);
    let second = second % (input.len() + 1);
    let start = first.min(second);
    let end = first.max(second);
    check_search::<A>(input, &input[start..end]);

    let sequence = Sequence::<A>::from_ascii(input).unwrap();
    let slice = sequence.slice(start..end).unwrap();
    assert_eq!(slice.to_string(), input[start..end].to_ascii_uppercase());
    let other = Sequence::<A>::from_ascii(pattern).unwrap();
    assert_eq!(
        sequence.concat(&other).unwrap().to_string(),
        format!(
            "{}{}",
            input.to_ascii_uppercase(),
            pattern.to_ascii_uppercase()
        )
    );
    let qualities: Vec<_> = (0..input.len())
        .map(|index| u8::try_from(index % 94).unwrap())
        .collect();
    let record =
        Record::<A>::from_ascii("generated mixed case", input, Some(qualities.clone())).unwrap();
    let region = record.slice(start..end).unwrap();
    assert_eq!(region.source_sequence(), &input[start..end]);
    assert_eq!(region.qualities(), Some(&qualities[start..end]));
    assert_eq!(region.sequence(), &slice);
    assert_eq!(region.header(), record.header());
    assert_eq!(record.source_sequence(), input);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn generated_dna_matches_search_slice_concat_and_record_oracles(
        input in "[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{0,160}",
        pattern in "[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{0,24}",
        first in any::<usize>(),
        second in any::<usize>(),
    ) {
        check_generated_operations::<Dna>(&input, &pattern, first, second);
    }

    #[test]
    fn generated_rna_matches_search_slice_concat_and_record_oracles(
        input in "[ACGURYSWKMBDHVNacguryswkmbdhvn]{0,160}",
        pattern in "[ACGURYSWKMBDHVNacguryswkmbdhvn]{0,24}",
        first in any::<usize>(),
        second in any::<usize>(),
    ) {
        check_generated_operations::<Rna>(&input, &pattern, first, second);
    }

    #[test]
    fn generated_protein_matches_search_slice_concat_and_record_oracles(
        input in "[ARNDCQEGHILKMFPSTWYVBZXJUOarndcqeghilkmfpstwyvbzxjuo]{0,160}",
        pattern in "[ARNDCQEGHILKMFPSTWYVBZXJUOarndcqeghilkmfpstwyvbzxjuo]{0,24}",
        first in any::<usize>(),
        second in any::<usize>(),
    ) {
        check_generated_operations::<Protein>(&input, &pattern, first, second);
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
fn sequence_operations_guide_runs_the_actual_cargo_example() -> Result<(), Box<dyn Error>> {
    let guide = include_str!("../../docs/api/sequence-operations-v1.md");
    let program = include_str!("../examples/sequence_operations.rs");
    let (_, program) = program.split_once("\n\n").unwrap();
    assert_eq!(
        guide_fence(guide, "sequence-operations-program", "rust").trim(),
        program.trim()
    );
    assert_eq!(
        guide_fence(guide, "sequence-operations-command", "bash"),
        "cargo bio sequence_operations\n"
    );
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let output = Command::new(env!("CARGO"))
        .args(["bio", "sequence_operations", "--offline", "--locked"])
        .current_dir(root)
        .output()?;
    assert_eq!(
        output.status.code(),
        Some(0),
        "example failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    assert_eq!(
        output.stdout,
        guide_fence(guide, "sequence-operations-output", "text").as_bytes()
    );
    Ok(())
}
