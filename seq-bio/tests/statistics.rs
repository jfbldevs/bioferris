//! Independent composition, concrete-expansion GC and bounded-window checks.

use std::cell::Cell;
use std::collections::BTreeMap;
use std::error::Error;
use std::io::{Cursor, Write};
use std::path::Path;
use std::process::Command;

use bioferris::io::{BioIoError, FastaReader, FastqReader, IoErrorKind, decode_input};
use bioferris::prelude::*;
use flate2::{Compression, write::GzEncoder};
use proptest::prelude::*;

const GC_TOLERANCE: f64 = 1e-12;

// Concrete sets, rather than the production GC weights or compact symbols.
// Equal weight per listed base is the independent weighted-GC oracle.
fn expansion(symbol: u8) -> &'static [u8] {
    match symbol.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
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
        _ => panic!("test expected an IUPAC nucleotide"),
    }
}

fn small_float(value: usize) -> f64 {
    f64::from(u32::try_from(value).expect("test inputs fit exactly in u32"))
}

fn gc_oracle(input: &str, policy: GcPolicy) -> (f64, u64, f64) {
    // All concrete-set sizes divide 12. Accumulate exact twelfths so the
    // independent oracle does not introduce order-dependent floating drift.
    let mut numerator_twelfths = 0;
    let mut denominator = 0_usize;
    for byte in input.bytes() {
        let bases = expansion(byte);
        let gc_bases = bases.iter().filter(|base| b"GC".contains(base)).count();
        let definite_gc = gc_bases == bases.len();
        match policy {
            GcPolicy::Remove => {
                if definite_gc || gc_bases == 0 {
                    denominator += 1;
                    numerator_twelfths += if definite_gc { 12 } else { 0 };
                }
            }
            GcPolicy::Ignore => {
                denominator += 1;
                numerator_twelfths += if definite_gc { 12 } else { 0 };
            }
            GcPolicy::Weighted => {
                denominator += 1;
                numerator_twelfths += gc_bases * (12 / bases.len());
            }
        }
    }
    let numerator = small_float(numerator_twelfths) / 12.0;
    let fraction = if denominator == 0 {
        0.0
    } else {
        numerator / small_float(denominator)
    };
    (numerator, u64::try_from(denominator).unwrap(), fraction)
}

fn close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= GC_TOLERANCE,
        "actual {actual:?}, expected {expected:?}"
    );
}

fn assert_gc(stats: &SequenceStatistics, input: &str, policy: GcPolicy) {
    let gc = stats.gc.as_ref().expect("nucleotide summaries expose GC");
    let (numerator, denominator, fraction) = gc_oracle(input, policy);
    assert_eq!(gc.policy, policy);
    assert_eq!(gc.denominator, denominator);
    close(gc.gc_equivalent, numerator);
    close(gc.fraction, fraction);
    assert!((0.0..=1.0).contains(&gc.fraction));
}

fn ambiguous(byte: u8, nucleotide: bool) -> bool {
    let byte = byte.to_ascii_uppercase();
    if nucleotide {
        !b"ACGTU".contains(&byte)
    } else {
        b"BZJX".contains(&byte)
    }
}

fn kmer_oracle(
    inputs: &[String],
    k: usize,
    policy: KmerAmbiguity,
    nucleotide: bool,
) -> (u64, u64, BTreeMap<String, u64>) {
    let mut windows = 0;
    let mut skipped = 0;
    let mut counts = BTreeMap::new();
    for input in inputs {
        let normalized = input.to_ascii_uppercase();
        for window in normalized.as_bytes().windows(k) {
            windows += 1;
            if policy == KmerAmbiguity::Skip
                && window.iter().any(|byte| ambiguous(*byte, nucleotide))
            {
                skipped += 1;
            } else {
                let word = std::str::from_utf8(window).unwrap().to_owned();
                *counts.entry(word).or_insert(0) += 1;
            }
        }
    }
    (windows, skipped, counts)
}

fn options(k: usize, ambiguity: KmerAmbiguity) -> StatisticsOptions {
    StatisticsOptions {
        kmers: Some(KmerOptions {
            ambiguity,
            ..KmerOptions::new(k)
        }),
        ..StatisticsOptions::default()
    }
}

#[test]
fn defaults_and_zero_records_are_distinct_from_one_empty_sequence() {
    let defaults = StatisticsOptions::default();
    assert_eq!(defaults.gc_policy, GcPolicy::Remove);
    assert!(defaults.kmers.is_none());
    let kmer_defaults = KmerOptions::new(3);
    assert_eq!(kmer_defaults.k, 3);
    assert_eq!(kmer_defaults.ambiguity, KmerAmbiguity::Skip);
    assert_eq!(kmer_defaults.max_windows, 1_000_000);
    assert_eq!(kmer_defaults.max_distinct, 100_000);
    assert_eq!(kmer_defaults.max_key_bytes, 8 * 1024 * 1024);

    let empty = dna("").unwrap().statistics().unwrap();
    assert_eq!((empty.records, empty.length, empty.ambiguous), (1, 0, 0));
    assert_eq!(
        (empty.minimum_length, empty.maximum_length),
        (Some(0), Some(0))
    );
    assert!(empty.composition.is_empty());
    assert!(empty.kmers.is_none());
    assert_gc(&empty, "", GcPolicy::Remove);

    let none = summarize_sequences(
        std::iter::empty::<&DnaSequence>(),
        options(2, KmerAmbiguity::Skip),
    )
    .unwrap();
    assert_eq!((none.records, none.length, none.ambiguous), (0, 0, 0));
    assert_eq!((none.minimum_length, none.maximum_length), (None, None));
    assert!(none.composition.is_empty());
    assert_gc(&none, "", GcPolicy::Remove);
    let kmers = none.kmers.unwrap();
    assert_eq!((kmers.windows, kmers.counted, kmers.skipped), (0, 0, 0));
    assert!(kmers.counts.is_empty());

    let protein_none = summarize_sequences(
        std::iter::empty::<&ProteinSequence>(),
        StatisticsOptions::default(),
    )
    .unwrap();
    assert!(protein_none.gc.is_none());
    assert_eq!(protein_none.records, 0);
    assert!(protein("").unwrap().statistics().unwrap().gc.is_none());
}

#[test]
fn every_iupac_triplet_matches_independent_concrete_expansion_for_dna_and_rna() {
    let symbols = b"ACGTRYSWKMBDHVN";
    let mut checked = 0;
    for first in symbols {
        for second in symbols {
            for third in symbols {
                let input = String::from_utf8(vec![*first, *second, *third]).unwrap();
                let sequence = dna(&input).unwrap();
                let rna_input = input.replace('T', "U");
                for policy in [GcPolicy::Remove, GcPolicy::Ignore, GcPolicy::Weighted] {
                    let options = StatisticsOptions {
                        gc_policy: policy,
                        ..StatisticsOptions::default()
                    };
                    let stats = sequence.statistics_with(options).unwrap();
                    assert_gc(&stats, &input, policy);
                    let rna_stats = rna(&rna_input).unwrap().statistics_with(options).unwrap();
                    assert_gc(&rna_stats, &rna_input, policy);
                    assert_eq!(stats.ambiguous, rna_stats.ambiguous);
                }
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 3375);
}

#[test]
fn public_biopython_gc_examples_match_the_documented_policies() {
    // Public Bio.SeqUtils 1.88 examples, not execution of its implementation.
    for (input, remove, ignore, weighted) in [
        ("ACTG", 0.5, 0.5, 0.5),
        ("ACTGSSSS", 0.75, 0.75, 0.75),
        ("ACTGN", 0.5, 0.4, 0.5),
        ("GDVV", 1.0, 0.25, 2.0 / 3.0),
        ("", 0.0, 0.0, 0.0),
    ] {
        for (policy, expected) in [
            (GcPolicy::Remove, remove),
            (GcPolicy::Ignore, ignore),
            (GcPolicy::Weighted, weighted),
        ] {
            let summary = dna(input)
                .unwrap()
                .statistics_with(StatisticsOptions {
                    gc_policy: policy,
                    ..StatisticsOptions::default()
                })
                .unwrap();
            close(summary.gc.unwrap().fraction, expected);
        }
    }
    close(
        rna("GGAUCUUCGGAUCU")
            .unwrap()
            .statistics()
            .unwrap()
            .gc
            .unwrap()
            .fraction,
        0.5,
    );
}

#[test]
fn composition_is_normalized_sorted_observed_only_and_ambiguity_is_explicit() {
    let stats = dna("acgtRYSWKMBDHVN").unwrap().statistics().unwrap();
    assert_eq!((stats.length, stats.ambiguous), (15, 11));
    assert_eq!(stats.composition.len(), 15);
    assert!(stats.composition.values().all(|value| *value == 1));
    assert_eq!(
        stats.composition.keys().collect::<String>(),
        "ABCDGHKMNRSTVWY"
    );
    assert_gc(&stats, "ACGTRYSWKMBDHVN", GcPolicy::Remove);
    let gc = stats.gc.unwrap();
    assert_eq!(gc.denominator, 6); // A/C/G/T/S/W, despite S/W being ambiguous.
    close(gc.gc_equivalent, 3.0);
    let only_n = dna("nnn").unwrap().statistics().unwrap();
    assert_eq!(only_n.composition, BTreeMap::from([('N', 3)]));
    assert_eq!(only_n.ambiguous, 3);
    assert_gc(&only_n, "NNN", GcPolicy::Remove);
}

#[test]
fn protein_u_and_o_are_concrete_but_b_z_j_x_are_ambiguous() {
    let sequence = protein("arndcqeghilkmfpstwyvbzxjuo").unwrap();
    for policy in [GcPolicy::Remove, GcPolicy::Ignore, GcPolicy::Weighted] {
        let stats = sequence
            .statistics_with(StatisticsOptions {
                gc_policy: policy,
                ..options(1, KmerAmbiguity::Skip)
            })
            .unwrap();
        assert_eq!((stats.length, stats.ambiguous), (26, 4));
        assert!(stats.gc.is_none());
        let kmers = stats.kmers.unwrap();
        assert_eq!((kmers.windows, kmers.counted, kmers.skipped), (26, 22, 4));
        assert_eq!(kmers.counts.get("U"), Some(&1));
        assert_eq!(kmers.counts.get("O"), Some(&1));
        for word in ["B", "Z", "J", "X"] {
            assert!(!kmers.counts.contains_key(word));
        }
    }
}

#[test]
fn aggregate_gc_uses_counts_and_kmers_never_cross_record_boundaries() {
    let sequences = [dna("G").unwrap(), dna("AAAA").unwrap(), dna("").unwrap()];
    let stats = summarize_sequences(&sequences, options(2, KmerAmbiguity::Skip)).unwrap();
    assert_eq!((stats.records, stats.length), (3, 5));
    assert_eq!(
        (stats.minimum_length, stats.maximum_length),
        (Some(0), Some(4))
    );
    close(stats.gc.unwrap().fraction, 0.2);
    let kmers = stats.kmers.unwrap();
    assert_eq!((kmers.windows, kmers.counted, kmers.skipped), (3, 3, 0));
    assert_eq!(kmers.counts, BTreeMap::from([("AA".to_owned(), 3)]));
    let short = [dna("AA").unwrap(), dna("TT").unwrap()];
    let stats = summarize_sequences(&short, options(3, KmerAmbiguity::Literal)).unwrap();
    assert_eq!(stats.kmers.unwrap().windows, 0);
}

#[test]
fn overlapping_kmers_skip_or_retain_literal_ambiguity_without_canonicalization() {
    for (input, policy, windows, counted, skipped, expected) in [
        ("AAAA", KmerAmbiguity::Skip, 3, 3, 0, vec![("AA", 3)]),
        (
            "ACGN",
            KmerAmbiguity::Skip,
            3,
            2,
            1,
            vec![("AC", 1), ("CG", 1)],
        ),
        (
            "ACGN",
            KmerAmbiguity::Literal,
            3,
            3,
            0,
            vec![("AC", 1), ("CG", 1), ("GN", 1)],
        ),
        ("SW", KmerAmbiguity::Skip, 1, 0, 1, vec![]),
        (
            "ATGCAT",
            KmerAmbiguity::Skip,
            5,
            5,
            0,
            vec![("AT", 2), ("CA", 1), ("GC", 1), ("TG", 1)],
        ),
    ] {
        let stats = dna(input)
            .unwrap()
            .statistics_with(options(2, policy))
            .unwrap();
        let kmers = stats.kmers.unwrap();
        assert_eq!(kmers.k, 2);
        assert_eq!(kmers.ambiguity, policy);
        assert_eq!(
            (kmers.windows, kmers.counted, kmers.skipped),
            (windows, counted, skipped)
        );
        let expected: BTreeMap<_, _> = expected
            .into_iter()
            .map(|(word, count)| (word.to_owned(), count))
            .collect();
        assert_eq!(kmers.counts, expected);
    }
}

fn assert_limit(error: &StatisticsError, resource: &str, limit: u64, observed: u64) {
    assert_eq!(error.code(), "limit_exceeded");
    assert!(!error.to_string().is_empty());
    assert!(error.source().is_none());
    assert!(matches!(error,
        StatisticsError::LimitExceeded { resource: actual_resource, limit: actual_limit, observed: actual_observed }
            if *actual_resource == resource && *actual_limit == limit && *actual_observed == observed
    ));
}

#[test]
fn k_is_checked_before_any_input_is_consumed() {
    for k in [0, 65, usize::MAX] {
        assert!(matches!(
            options(k, KmerAmbiguity::Skip).validate(),
            Err(StatisticsError::InvalidOptions(_))
        ));
        let consumed = Cell::new(0);
        let sequence = dna("ACGT").unwrap();
        let input = std::iter::once(&sequence).inspect(|_| consumed.set(consumed.get() + 1));
        let error = summarize_sequences(input, options(k, KmerAmbiguity::Skip)).unwrap_err();
        assert_eq!(error.code(), "invalid_statistics_options");
        assert!(matches!(error, StatisticsError::InvalidOptions(_)));
        assert_eq!(consumed.get(), 0);

        let input = std::iter::once(Record::<Dna>::from_ascii("read", "AC", None))
            .inspect(|_| consumed.set(consumed.get() + 1));
        assert!(summarize_records(input, options(k, KmerAmbiguity::Skip)).is_err());
        assert_eq!(consumed.get(), 0);
    }
    let stats = dna("A".repeat(65))
        .unwrap()
        .statistics_with(options(64, KmerAmbiguity::Skip))
        .unwrap();
    let kmers = stats.kmers.unwrap();
    assert_eq!((kmers.windows, kmers.counted), (2, 2));
    assert_eq!(kmers.counts, BTreeMap::from([("A".repeat(64), 2)]));
}

#[test]
fn windows_budget_includes_skips_and_preflights_the_entire_next_record() {
    let limited = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_windows: 1,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    for input in ["AAAA", "NNNN"] {
        assert_limit(
            &dna(input).unwrap().statistics_with(limited).unwrap_err(),
            "kmer_windows",
            1,
            3,
        );
    }
    let sequences = [dna("AA").unwrap(), dna("AAA").unwrap(), dna("T").unwrap()];
    let consumed = Cell::new(0);
    let input = sequences
        .iter()
        .inspect(|_| consumed.set(consumed.get() + 1));
    assert_limit(
        &summarize_sequences(input, limited).unwrap_err(),
        "kmer_windows",
        1,
        3,
    );
    assert_eq!(consumed.get(), 2);
    let skipped = [dna("NN").unwrap(), dna("NN").unwrap()];
    assert_limit(
        &summarize_sequences(&skipped, limited).unwrap_err(),
        "kmer_windows",
        1,
        2,
    );
}

#[test]
fn distinct_key_budget_is_aggregate_and_repeated_keys_do_not_consume_it_twice() {
    let limited = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_distinct: 1,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    let repeated = [dna("AA").unwrap(), dna("AA").unwrap()];
    assert_eq!(
        summarize_sequences(&repeated, limited)
            .unwrap()
            .kmers
            .unwrap()
            .counts
            .get("AA"),
        Some(&2)
    );
    let distinct = [dna("AA").unwrap(), dna("AC").unwrap(), dna("AG").unwrap()];
    let consumed = Cell::new(0);
    let input = distinct
        .iter()
        .inspect(|_| consumed.set(consumed.get() + 1));
    assert_limit(
        &summarize_sequences(input, limited).unwrap_err(),
        "kmer_distinct",
        1,
        2,
    );
    assert_eq!(consumed.get(), 2);
}

#[test]
fn stored_key_byte_budget_is_aggregate_and_checks_before_a_new_key() {
    let limited = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_key_bytes: 3,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    let repeated = [dna("AA").unwrap(), dna("AA").unwrap()];
    assert!(summarize_sequences(&repeated, limited).is_ok());
    let distinct = [dna("AA").unwrap(), dna("AC").unwrap()];
    assert_limit(
        &summarize_sequences(&distinct, limited).unwrap_err(),
        "kmer_key_bytes",
        3,
        4,
    );
    let exact = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_key_bytes: 64,
            ..KmerOptions::new(64)
        }),
        ..StatisticsOptions::default()
    };
    assert!(dna("A".repeat(65)).unwrap().statistics_with(exact).is_ok());
}

#[test]
fn zero_budgets_allow_no_windows_or_no_keys_but_never_waive_a_limit() {
    let no_windows = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_windows: 0,
            max_distinct: 0,
            max_key_bytes: 0,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    for input in ["", "A", "N"] {
        let stats = dna(input).unwrap().statistics_with(no_windows).unwrap();
        assert_eq!(stats.kmers.unwrap().windows, 0);
    }
    assert_limit(
        &dna("NN").unwrap().statistics_with(no_windows).unwrap_err(),
        "kmer_windows",
        0,
        1,
    );
    let no_keys = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_distinct: 0,
            max_key_bytes: 0,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    let skipped = dna("NNN")
        .unwrap()
        .statistics_with(no_keys)
        .unwrap()
        .kmers
        .unwrap();
    assert_eq!(
        (skipped.windows, skipped.skipped, skipped.counted),
        (2, 2, 0)
    );
    assert_limit(
        &dna("AC").unwrap().statistics_with(no_keys).unwrap_err(),
        "kmer_distinct",
        0,
        1,
    );
    let no_bytes = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_key_bytes: 0,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    assert_limit(
        &dna("AC").unwrap().statistics_with(no_bytes).unwrap_err(),
        "kmer_key_bytes",
        0,
        2,
    );
}

#[test]
fn input_failure_retains_context_and_stops_before_later_records() {
    let consumed = Cell::new(0);
    let records = [
        Record::<Dna>::from_ascii("good", "ACGT", None),
        Err(BioIoError::new(IoErrorKind::InvalidQuality, "fastq-sanger", "bad quality").at(2, 8)),
        Record::<Dna>::from_ascii("unvisited", "AAAA", None),
    ];
    let input = records
        .into_iter()
        .inspect(|_| consumed.set(consumed.get() + 1));
    let error = summarize_records(input, StatisticsOptions::default()).unwrap_err();
    assert_eq!(consumed.get(), 2);
    assert_eq!(error.code(), "invalid_quality");
    assert!(error.source().is_some());
    let StatisticsError::Input(source) = error else {
        panic!("input error must be retained")
    };
    assert_eq!(source.kind(), IoErrorKind::InvalidQuality);
    assert_eq!((source.record(), source.line()), (Some(2), Some(8)));
    assert_eq!(source.message(), "bad quality");
}

#[test]
fn actual_late_fasta_error_never_yields_a_successful_partial_summary() {
    let input = ">first\nACGT\n>invalid\nU\n>later\nAAAA\n";
    let consumed = Cell::new(0);
    let records = FastaReader::<_, Dna>::new(Cursor::new(input))
        .inspect(|_| consumed.set(consumed.get() + 1));
    let error = summarize_records(records, options(2, KmerAmbiguity::Skip)).unwrap_err();
    assert_eq!(consumed.get(), 2);
    assert_eq!(error.code(), "invalid_sequence");
    let StatisticsError::Input(source) = error else {
        panic!("reader error must be retained")
    };
    assert_eq!(source.record(), Some(2));
    assert_eq!(source.format(), "fasta");
}

#[test]
fn record_budget_errors_stop_input_iteration_at_the_first_excess_record() {
    for (limits, resource, limit, observed) in [
        (
            KmerOptions {
                max_windows: 1,
                ..KmerOptions::new(2)
            },
            "kmer_windows",
            1,
            2,
        ),
        (
            KmerOptions {
                max_distinct: 1,
                ..KmerOptions::new(2)
            },
            "kmer_distinct",
            1,
            2,
        ),
        (
            KmerOptions {
                max_key_bytes: 2,
                ..KmerOptions::new(2)
            },
            "kmer_key_bytes",
            2,
            4,
        ),
    ] {
        let consumed = Cell::new(0);
        let input = ["AA", "AC", "AG"]
            .into_iter()
            .map(|sequence| Record::<Dna>::from_ascii("r", sequence, None))
            .inspect(|_| consumed.set(consumed.get() + 1));
        let error = summarize_records(
            input,
            StatisticsOptions {
                kmers: Some(limits),
                ..StatisticsOptions::default()
            },
        )
        .unwrap_err();
        assert_limit(&error, resource, limit, observed);
        assert_eq!(consumed.get(), 2);
    }
}

#[test]
fn plain_and_gzip_fastq_summaries_agree_and_qualities_do_not_weight_counts() {
    let input = include_bytes!("fixtures/io-example.fastq");
    let plain = summarize_records(
        FastqReader::<_, Dna>::new(Cursor::new(input)),
        options(2, KmerAmbiguity::Skip),
    )
    .unwrap();
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(input).unwrap();
    let compressed = encoder.finish().unwrap();
    let reader = decode_input(Cursor::new(compressed.clone())).unwrap();
    let zipped = summarize_records(
        FastqReader::<_, Dna>::new(reader),
        options(2, KmerAmbiguity::Skip),
    )
    .unwrap();
    assert_eq!(zipped, plain);
    assert_eq!(
        plain.composition,
        BTreeMap::from([('A', 1), ('C', 1), ('G', 1), ('T', 1)])
    );
    assert_gc(&plain, "ACGT", GcPolicy::Remove);
    assert_eq!(
        plain.kmers.unwrap().counts,
        BTreeMap::from([
            ("AC".to_owned(), 1),
            ("CG".to_owned(), 1),
            ("GT".to_owned(), 1)
        ])
    );

    let truncated = &compressed[..compressed.len() - 2];
    let reader = decode_input(Cursor::new(truncated.to_vec())).unwrap();
    assert!(matches!(
        summarize_records(
            FastqReader::<_, Dna>::new(reader),
            StatisticsOptions::default()
        ),
        Err(StatisticsError::Input(_))
    ));
}

#[test]
fn real_phix174_fasta_summary_matches_independent_source_symbol_counts() {
    let input = include_str!("fixtures/phix174-e.fasta");
    let sequence: String = input
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let stats = summarize_records(
        FastaReader::<_, Dna>::new(Cursor::new(input)),
        StatisticsOptions::default(),
    )
    .unwrap();
    assert_eq!((stats.records, stats.length, stats.ambiguous), (1, 276, 0));
    assert_eq!(
        (stats.minimum_length, stats.maximum_length),
        (Some(276), Some(276))
    );
    let mut expected = BTreeMap::new();
    for symbol in sequence.chars() {
        *expected.entry(symbol).or_insert(0) += 1;
    }
    assert_eq!(stats.composition, expected);
    assert_gc(&stats, &sequence, GcPolicy::Remove);
}

#[test]
fn count_overflow_is_an_explicit_error_category() {
    let error = StatisticsError::CountOverflow;
    assert_eq!(error.code(), "count_overflow");
    assert!(error.source().is_none());
    assert!(!error.to_string().is_empty());
}

fn verify_generated_summary<A: StatisticsAlphabet>(
    inputs: &[String],
    k: usize,
    ambiguity: KmerAmbiguity,
    nucleotide: bool,
) {
    let sequences: Vec<_> = inputs
        .iter()
        .map(|input| Sequence::<A>::from_ascii(input).unwrap())
        .collect();
    let joined: String = inputs
        .iter()
        .map(|input| input.to_ascii_uppercase())
        .collect();
    let mut composition = BTreeMap::new();
    for symbol in joined.chars() {
        *composition.entry(symbol).or_insert(0) += 1;
    }
    let expected_ambiguous = u64::try_from(
        joined
            .bytes()
            .filter(|byte| ambiguous(*byte, nucleotide))
            .count(),
    )
    .unwrap();
    let (windows, skipped, expected_kmers) = kmer_oracle(inputs, k, ambiguity, nucleotide);
    for policy in [GcPolicy::Remove, GcPolicy::Ignore, GcPolicy::Weighted] {
        let settings = StatisticsOptions {
            gc_policy: policy,
            ..options(k, ambiguity)
        };
        let stats = summarize_sequences(&sequences, settings).unwrap();
        assert_eq!(stats.records, u64::try_from(inputs.len()).unwrap());
        assert_eq!(stats.length, u64::try_from(joined.len()).unwrap());
        assert_eq!(stats.ambiguous, expected_ambiguous);
        assert_eq!(stats.minimum_length, inputs.iter().map(String::len).min());
        assert_eq!(stats.maximum_length, inputs.iter().map(String::len).max());
        assert_eq!(stats.composition, composition);
        if nucleotide {
            assert_gc(&stats, &joined, policy);
        } else {
            assert!(stats.gc.is_none());
        }
        let kmers = stats.kmers.as_ref().unwrap();
        assert_eq!(kmers.k, k);
        assert_eq!(kmers.ambiguity, ambiguity);
        assert_eq!(
            (kmers.windows, kmers.counted, kmers.skipped),
            (windows, windows - skipped, skipped)
        );
        assert_eq!(kmers.counts, expected_kmers);
        assert_eq!(kmers.counts.values().sum::<u64>(), kmers.counted);
        let reversed = summarize_sequences(sequences.iter().rev(), settings).unwrap();
        assert_eq!(stats, reversed);
    }
    for (original, sequence) in inputs.iter().zip(sequences) {
        assert_eq!(sequence.to_string(), original.to_ascii_uppercase());
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn generated_dna_summaries_match_independent_gc_and_window_oracles(
        inputs in prop::collection::vec("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{0,80}", 0..7),
        k in 1_usize..=12,
        literal in any::<bool>(),
    ) {
        verify_generated_summary::<Dna>(&inputs, k, if literal { KmerAmbiguity::Literal } else { KmerAmbiguity::Skip }, true);
    }

    #[test]
    fn generated_rna_summaries_match_independent_gc_and_window_oracles(
        inputs in prop::collection::vec("[ACGURYSWKMBDHVNacguryswkmbdhvn]{0,80}", 0..7),
        k in 1_usize..=12,
        literal in any::<bool>(),
    ) {
        verify_generated_summary::<Rna>(&inputs, k, if literal { KmerAmbiguity::Literal } else { KmerAmbiguity::Skip }, true);
    }

    #[test]
    fn generated_protein_summaries_match_independent_symbol_and_window_oracles(
        inputs in prop::collection::vec("[ARNDCQEGHILKMFPSTWYVBZXJUOarndcqeghilkmfpstwyvbzxjuo]{0,80}", 0..7),
        k in 1_usize..=12,
        literal in any::<bool>(),
    ) {
        verify_generated_summary::<Protein>(&inputs, k, if literal { KmerAmbiguity::Literal } else { KmerAmbiguity::Skip }, false);
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
fn statistics_guide_runs_the_actual_cargo_example() -> Result<(), Box<dyn Error>> {
    let guide = include_str!("../../docs/api/sequence-statistics-v1.md");
    let program = include_str!("../examples/sequence_statistics.rs");
    let (_, program) = program.split_once("\n\n").unwrap();
    assert_eq!(
        guide_fence(guide, "statistics-program", "rust").trim(),
        program.trim()
    );
    assert_eq!(
        guide_fence(guide, "statistics-command", "bash"),
        "cargo bio sequence_statistics\n"
    );
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let output = Command::new(env!("CARGO"))
        .args(["bio", "sequence_statistics", "--offline", "--locked"])
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
        guide_fence(guide, "statistics-output", "text").as_bytes()
    );
    Ok(())
}
