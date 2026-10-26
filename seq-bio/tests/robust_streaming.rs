//! Bounded generated-corpus checks of public import and compression interfaces.
//!
//! These are deterministic regression/property checks, not scientific qualification
//! or an assertion that every possible byte stream is safe.

mod robust_support;

use std::fmt::Write as _;
use std::io::{Cursor, Read, Write};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use bioferris::io::alignment::{AlignmentFormat, parse_alignment};
use bioferris::io::annotated::{parse_embl, parse_genbank};
use bioferris::io::bgzf::{BGZF_EOF_MARKER, BgzfError, BgzfReader, BgzfWriteOptions, BgzfWriter};
use bioferris::io::blast::{
    BlastDocument, BlastQuery, BlastReadLimits, BlastTabularFormat, BlastTabularOptions,
    BlastTabularReader, BlastXmlFormat, BlastXmlOptions, BlastXmlReader, parse_blast_tabular,
    parse_blast_xml,
};
use bioferris::io::graph::{parse_dot_reader, parse_gfa_reader};
use bioferris::io::gzip::GzipWriter;
use bioferris::io::result_json::{ImportedResult, parse_result_json_reader};
use bioferris::io::{
    FastaReader, FastqReader, IoErrorKind, ReadLimits, decode_input, write_record_fasta,
    write_record_fastq,
};
use bioferris::prelude::{PhylipLayout, PhylipNameStyle};
use bioferris::{Dna, Record};
use robust_support::{
    AdversarialRead, CHUNKS, Generator, Schedule, StreamOutcome, bgzf_member, buffered, crc32,
    drain, gzip_member, mutations, outcome, random_corpus,
};
use serde_json::json;

const FASTA: &[u8] = b">sample mixed case\r\nAc\r\ngTN\r\n>second\nRY\n";
const FASTQ: &[u8] =
    b"@sample mixed case\r\nAc\r\ngTN\r\n+sample mixed case\r\n!I\r\n~+@\r\n@second\nRY\n+\n!~\n";
const MAX_DECODED: usize = 16 * 1024;

fn sequence_stream(
    input: &[u8],
    fastq: bool,
    schedule: Schedule,
    limits: ReadLimits,
) -> StreamOutcome<Record<Dna>> {
    let (reader, calls) = buffered(input, schedule);
    if fastq {
        drain(
            FastqReader::<_, Dna>::with_limits(reader, limits),
            input.len(),
            &calls,
        )
    } else {
        drain(
            FastaReader::<_, Dna>::with_limits(reader, limits),
            input.len(),
            &calls,
        )
    }
}

fn sequence(input: &[u8], fastq: bool, schedule: Schedule) -> StreamOutcome<Record<Dna>> {
    sequence_stream(input, fastq, schedule, ReadLimits::default())
}

fn assert_record_invariants(record: &Record<Dna>, fastq: bool) {
    assert!(!record.id().is_empty());
    assert_eq!(record.source_sequence().len(), record.sequence().len());
    assert_eq!(
        record.source_sequence().to_ascii_uppercase(),
        record.sequence().to_string()
    );
    let mut written = Vec::new();
    if fastq {
        let qualities = record.qualities().unwrap();
        assert_eq!(qualities.len(), record.sequence().len());
        assert!(qualities.iter().all(|quality| *quality <= 93));
        write_record_fastq(&mut written, record).unwrap();
    } else {
        assert!(record.qualities().is_none());
        write_record_fasta(&mut written, record, bioferris::io::LossPolicy::Reject).unwrap();
    }
    assert_eq!(
        sequence(&written, fastq, Schedule::chunks(1)).values,
        [record.clone()]
    );
}

#[test]
fn generated_record_grammars_match_independent_records_at_every_chunk_schedule() {
    let mut generator = Generator::new(0x925a_10f3);
    for case in 0..16 {
        let mut sequence_text = String::new();
        let mut quality_text = String::new();
        let mut plain_records = Vec::new();
        let mut scored_records = Vec::new();
        for record in 0..3 {
            let header = format!("read_{case}_{record}_β  generated");
            let length = 1 + generator.index(31);
            let source = generator.dna(length);
            let scores: Vec<_> = (0..source.len()).map(|_| generator.byte() % 94).collect();
            let wrap = 1 + generator.index(9);
            writeln!(sequence_text, ">{header}\r").unwrap();
            writeln!(quality_text, "@{header}\r").unwrap();
            for bytes in source.as_bytes().chunks(wrap) {
                let line = std::str::from_utf8(bytes).unwrap();
                writeln!(sequence_text, "{line}\r").unwrap();
                writeln!(quality_text, "{line}\r").unwrap();
            }
            writeln!(quality_text, "+{header}\r").unwrap();
            for values in scores.chunks(wrap) {
                for value in values {
                    quality_text.push(char::from(value + 33));
                }
                quality_text.push_str("\r\n");
            }
            plain_records.push(Record::from_ascii(&header, &source, None).unwrap());
            scored_records.push(Record::from_ascii(&header, &source, Some(scores)).unwrap());
        }
        for chunk in CHUNKS {
            for (input, fastq, expected) in [
                (sequence_text.as_bytes(), false, &plain_records),
                (quality_text.as_bytes(), true, &scored_records),
            ] {
                let observed = sequence(input, fastq, Schedule::chunks(chunk));
                assert_eq!(&observed.values, expected, "case {case}, chunk {chunk}");
                assert!(observed.failure.is_none(), "{observed:?}");
            }
        }
    }
}

#[test]
fn seeded_sequence_byte_edits_and_arbitrary_bytes_preserve_semantics_or_fuse() {
    for (input, fastq) in [(FASTA, false), (FASTQ, true)] {
        let corpus = mutations(input, 0x81e4_23a7)
            .into_iter()
            .chain(random_corpus());
        for (case, bytes) in corpus.enumerate() {
            let expected = sequence(&bytes, fastq, Schedule::chunks(64));
            for record in &expected.values {
                assert_record_invariants(record, fastq);
            }
            for chunk in CHUNKS {
                assert_eq!(
                    sequence(&bytes, fastq, Schedule::chunks(chunk)),
                    expected,
                    "case {case}, fastq {fastq}, chunk {chunk}, bytes {bytes:?}"
                );
            }
        }
    }
}

#[test]
fn line_record_and_source_limits_do_not_depend_on_transport_chunks() {
    for (input, fastq) in [(FASTA, false), (FASTQ, true)] {
        for maximum in 0..=input.len() + 1 {
            for limits in [
                ReadLimits {
                    max_line_bytes: maximum,
                    ..ReadLimits::default()
                },
                ReadLimits {
                    max_record_symbols: maximum,
                    ..ReadLimits::default()
                },
                ReadLimits {
                    max_total_bytes: maximum,
                    ..ReadLimits::default()
                },
            ] {
                let expected = sequence_stream(input, fastq, Schedule::chunks(64), limits);
                if let Some(error) = &expected.failure {
                    assert_eq!(error.kind, IoErrorKind::LimitExceeded);
                }
                for chunk in CHUNKS {
                    assert_eq!(
                        sequence_stream(input, fastq, Schedule::chunks(chunk), limits),
                        expected,
                        "maximum {maximum}, fastq {fastq}, chunk {chunk}"
                    );
                }
            }
        }
    }
}

#[test]
fn interrupted_sequence_reads_retry_without_losing_a_partial_line() {
    for (input, fastq) in [(FASTA, false), (FASTQ, true)] {
        let expected = sequence(input, fastq, Schedule::chunks(64));
        for position in 0..=input.len() {
            assert_eq!(
                sequence(input, fastq, Schedule::interrupted(position)),
                expected,
                "transient interruption at {position}, fastq {fastq}"
            );
        }
    }
}

#[test]
fn permanent_sequence_transport_errors_are_typed_and_fused_at_every_offset() {
    for (input, fastq) in [(FASTA, false), (FASTQ, true)] {
        for position in 0..=input.len() {
            let observed = sequence(input, fastq, Schedule::failing(position));
            assert_eq!(observed.failure.unwrap().kind, IoErrorKind::Io);
        }
    }
}

fn alignment_inputs(first: &str, second: &str) -> Vec<(AlignmentFormat, String)> {
    let width = first.len();
    let mut inputs = vec![
        (
            AlignmentFormat::Fasta,
            format!(">first\n{first}\n>second\n{second}\n"),
        ),
        (
            AlignmentFormat::Clustal,
            format!("CLUSTAL W\n\nfirst {first}\nsecond {second}\n"),
        ),
        (
            AlignmentFormat::Stockholm,
            format!("# STOCKHOLM 1.0\nfirst {first}\nsecond {second}\n//\n"),
        ),
        (
            AlignmentFormat::Nexus,
            format!(
                "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR={width};\nFORMAT DATATYPE=DNA GAP=- MISSING=?;\nMATRIX\nfirst {first}\nsecond {second}\n;\nEND;\n"
            ),
        ),
    ];
    for names in [PhylipNameStyle::Strict, PhylipNameStyle::Relaxed] {
        for layout in [PhylipLayout::Sequential, PhylipLayout::Interleaved] {
            let rows = if names == PhylipNameStyle::Strict {
                format!("first     {first}\nsecond    {second}\n")
            } else {
                format!("first {first}\nsecond {second}\n")
            };
            inputs.push((
                AlignmentFormat::Phylip { names, layout },
                format!("2 {width}\n{rows}"),
            ));
        }
    }
    inputs
}

#[test]
fn generated_alignment_profiles_keep_exact_rectangles_across_chunks_and_mutations() {
    let mut generator = Generator::new(0x549b_c621);
    for case in 0..6 {
        let width = 4 + case * 3;
        let first = generator.dna(width).to_ascii_uppercase();
        let mut second = generator.dna(width).to_ascii_uppercase().into_bytes();
        second[case % width] = b'-';
        let second = String::from_utf8(second).unwrap();
        for (format, input) in alignment_inputs(&first, &second) {
            let expected = parse_alignment::<Dna>(input.as_bytes(), format).unwrap();
            assert_eq!(expected.width(), width);
            assert_eq!(expected.rows()[0].aligned(), first);
            assert_eq!(expected.rows()[1].aligned(), second);
            for bytes in mutations(input.as_bytes(), 0x65b8_a2d9) {
                let baseline = outcome(parse_alignment::<Dna>(Cursor::new(&bytes), format));
                for chunk in CHUNKS {
                    let (reader, _) = buffered(&bytes, Schedule::chunks(chunk));
                    assert_eq!(outcome(parse_alignment::<Dna>(reader, format)), baseline);
                }
            }
        }
    }
}

fn graph_inputs(sequence: &str) -> (String, String) {
    let mut gfa = String::from("H\tVN:Z:1.0\n");
    let mut dot = String::from("digraph POASTA {\n  rankdir=LR;\n");
    for (index, symbol) in sequence.chars().enumerate() {
        writeln!(gfa, "S\tn{index}\t{symbol}\tCL:i:{index}").unwrap();
        writeln!(
            dot,
            "  n{index} [label=\"n{index}\\n{symbol}\\nclass {index}\"];"
        )
        .unwrap();
    }
    for index in 1..sequence.len() {
        let previous = index - 1;
        writeln!(gfa, "L\tn{previous}\t+\tn{index}\t+\t0M\tRC:i:1").unwrap();
        writeln!(dot, "  n{previous} -> n{index} [label=\"1\"];").unwrap();
    }
    let path = (0..sequence.len())
        .map(|index| format!("n{index}+"))
        .collect::<Vec<_>>()
        .join(",");
    writeln!(gfa, "P\tgenerated\t{path}\t*").unwrap();
    dot.push_str("}\n");
    (gfa, dot)
}

#[test]
fn generated_graph_and_result_imports_preserve_models_under_mutated_transport() {
    let mut generator = Generator::new(0x415f_23cd);
    for length in 1..=8 {
        let source = generator.dna(length).to_ascii_uppercase();
        let (gfa, dot) = graph_inputs(&source);
        let graph = parse_gfa_reader::<Dna>(gfa.as_bytes()).unwrap();
        let topology = parse_dot_reader::<Dna>(dot.as_bytes()).unwrap();
        assert_eq!(graph.path_sequence(0).as_deref(), Some(source.as_str()));
        assert_eq!(graph.nodes(), topology.nodes());
        assert_eq!(graph.edges(), topology.edges());
        let text = &source[1..];
        let value = json!({"schema":"align-output-v1", "kind":"myers_distance",
            "algorithm":"Myers bit-vector", "mode":"global unit edit distance", "alphabet":"DNA",
            "pattern":{"length":length,"sequence":source}, "text":{"length":length-1,"sequence":text}, "distance":1});
        let result = value.to_string();
        let ImportedResult::Alignment(imported) =
            parse_result_json_reader(result.as_bytes()).unwrap()
        else {
            panic!("expected native result");
        };
        assert_eq!(imported.value(), &value);
        for (profile, input) in [&gfa, &dot, &result].into_iter().enumerate() {
            for bytes in mutations(input.as_bytes(), 0x931c_62ab) {
                for chunk in CHUNKS {
                    let (reader, _) = buffered(&bytes, Schedule::chunks(chunk));
                    if profile == 0 {
                        assert_eq!(
                            outcome(parse_gfa_reader::<Dna>(reader)),
                            outcome(parse_gfa_reader::<Dna>(Cursor::new(&bytes)))
                        );
                    } else if profile == 1 {
                        assert_eq!(
                            outcome(parse_dot_reader::<Dna>(reader)),
                            outcome(parse_dot_reader::<Dna>(Cursor::new(&bytes)))
                        );
                    } else {
                        assert_eq!(
                            outcome(parse_result_json_reader(reader)),
                            outcome(parse_result_json_reader(Cursor::new(&bytes)))
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn eager_document_transports_retry_interruptions_and_reject_permanent_failures() {
    let alignment = ">first\nAC-GT\n>second\nACTGT\n";
    let (gfa, dot) = graph_inputs("ACGT");
    let result = json!({"schema":"align-output-v1", "kind":"myers_distance", "algorithm":"Myers bit-vector", "mode":"global unit edit distance", "alphabet":"DNA", "pattern":{"length":4,"sequence":"ACGT"},"text":{"length":3,"sequence":"CGT"},"distance":1}).to_string();
    for (profile, text) in [alignment, &gfa, &dot, &result].into_iter().enumerate() {
        for position in 0..=text.len() {
            for failing in [false, true] {
                let schedule = if failing {
                    Schedule::failing(position)
                } else {
                    Schedule::interrupted(position)
                };
                let (reader, _) = buffered(text.as_bytes(), schedule);
                let failure = match profile {
                    0 => parse_alignment::<Dna>(reader, AlignmentFormat::Fasta).err(),
                    1 => parse_gfa_reader::<Dna>(reader).err(),
                    2 => parse_dot_reader::<Dna>(reader).err(),
                    _ => parse_result_json_reader(reader).err(),
                };
                if failing {
                    assert_eq!(
                        failure.unwrap().kind(),
                        IoErrorKind::Io,
                        "{profile}, {position}"
                    );
                } else {
                    assert!(failure.is_none(), "{profile}, {position}: {failure:?}");
                }
            }
        }
    }
}

fn reverse_complement(input: &str) -> String {
    input
        .bytes()
        .rev()
        .map(|byte| {
            char::from(match byte.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'N' => b'N',
                b'R' => b'Y',
                b'Y' => b'R',
                _ => panic!("generator alphabet"),
            })
        })
        .collect()
}

#[test]
fn generated_annotated_grammars_keep_case_qualifiers_and_join_orientation() {
    let mut generator = Generator::new(0x3f92_a8e1);
    for length in 2..=24 {
        let source = generator.dna(length);
        let split = length / 2;
        let location = format!("complement(join(1..{split},{}..{length}))", split + 1);
        let genbank = format!(
            "LOCUS       GENERATED {length} bp DNA linear SYN 18-OCT-2026\nDEFINITION  Generated robustness input.\nFEATURES             Location/Qualifiers\n     misc_feature    {location}\n                     /note=\"seeded {length}\"\nORIGIN\n        1 {source}\n//\n"
        );
        let embl = format!(
            "ID   GENERATED; SV 1; linear; genomic DNA; STD; SYN; {length} BP.\nDE   Generated robustness input.\nFH   Key             Location/Qualifiers\nFT   misc_feature    {location}\nFT                   /note=\"seeded {length}\"\nSQ   Sequence {length} BP;\n     {source} {length}\n//\n"
        );
        for (is_embl, input) in [(false, genbank), (true, embl)] {
            let parse = |text: &str| {
                if is_embl {
                    parse_embl::<Dna>(text)
                } else {
                    parse_genbank::<Dna>(text)
                }
            };
            let records = parse(&input).unwrap();
            assert_eq!(records.len(), 1);
            let record = &records[0];
            assert_eq!(record.record().source_sequence(), source);
            assert_eq!(
                record.features()[0].qualifiers[0].value,
                Some(format!("seeded {length}"))
            );
            assert_eq!(
                record.features()[0]
                    .extract(record.record().sequence())
                    .unwrap()
                    .to_string(),
                reverse_complement(&source)
            );
            let wrong_length = input
                .replacen(&format!("{length} bp"), &format!("{} bp", length + 1), 1)
                .replacen(&format!("{length} BP."), &format!("{} BP.", length + 1), 1);
            assert_eq!(
                parse(&wrong_length).unwrap_err().kind(),
                IoErrorKind::LengthMismatch
            );
            for bytes in mutations(input.as_bytes(), 0xe153_9c42) {
                // These public string APIs cannot receive invalid UTF-8; no fake
                // transport/UTF-8 coverage is claimed for this subtest.
                if let Ok(text) = std::str::from_utf8(&bytes) {
                    if let Ok(records) = parse(text) {
                        for record in records {
                            assert_record_invariants(record.record(), false);
                        }
                    }
                }
            }
        }
    }
}

fn blast_limits() -> BlastReadLimits {
    BlastReadLimits {
        read: ReadLimits {
            max_line_bytes: 4096,
            max_record_symbols: 256,
            max_total_bytes: 8192,
        },
        max_queries: 8,
        max_hits_per_query: 8,
        max_hsps_per_query: 16,
        max_query_bytes: 4096,
        max_query_id_bytes: 256,
    }
}

fn tabular_options(commented: bool) -> BlastTabularOptions {
    BlastTabularOptions {
        format: if commented {
            BlastTabularFormat::Commented
        } else {
            BlastTabularFormat::Uncommented
        },
        limits: blast_limits(),
        ..BlastTabularOptions::default()
    }
}

fn tabular_input(commented: bool) -> String {
    let rows = "q_β\tsubject\t100\t4\t0\t0\t1\t4\t9\t6\t1.00e-4\t8.0\nq_β\tsubject\t100\t4\t0\t0\t5\t8\t1\t4\t2.00e-4\t7.0\n";
    if commented {
        format!(
            "# BLASTN synthetic\n# Query: q_β source description\n# Database: generated\n# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n# 2 hits found\n{rows}# BLASTN synthetic\n# Query: empty no matches\n# 0 hits found\n# BLAST processed 2 queries\n"
        )
    } else {
        rows.to_owned()
    }
}

fn tabular_stream(input: &[u8], commented: bool, schedule: Schedule) -> StreamOutcome<BlastQuery> {
    let (reader, calls) = buffered(input, schedule);
    drain(
        BlastTabularReader::new(reader, tabular_options(commented)).unwrap(),
        input.len(),
        &calls,
    )
}

fn xml_input(format: BlastXmlFormat) -> String {
    let hsp = [
        ("num", "1"),
        ("bit-score", "8.0"),
        ("score", "4"),
        ("evalue", "1.00e-4"),
        ("query-from", "1"),
        ("query-to", "4"),
        ("hit-from", "9"),
        ("hit-to", "6"),
        ("query-frame", "1"),
        ("hit-frame", "-1"),
        ("identity", "4"),
        ("positive", "4"),
        ("gaps", "0"),
        ("align-len", "4"),
        ("qseq", "ACGT"),
        ("hseq", "ACGT"),
        ("midline", "||||"),
    ];
    let prefix = if format == BlastXmlFormat::Legacy {
        "Hsp_"
    } else {
        ""
    };
    let mut fields = String::new();
    for (name, value) in hsp {
        write!(fields, "<{prefix}{name}>{value}</{prefix}{name}>").unwrap();
    }
    if format == BlastXmlFormat::Legacy {
        format!(
            "<?xml version=\"1.0\"?><BlastOutput><BlastOutput_program>blastn</BlastOutput_program><BlastOutput_version>synthetic</BlastOutput_version><BlastOutput_reference>generated &amp; local</BlastOutput_reference><BlastOutput_db>db</BlastOutput_db><BlastOutput_query-ID>q_β</BlastOutput_query-ID><BlastOutput_query-def>source &amp; description</BlastOutput_query-def><BlastOutput_query-len>8</BlastOutput_query-len><BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect><Parameters_gap-open>0</Parameters_gap-open><Parameters_gap-extend>0</Parameters_gap-extend></Parameters></BlastOutput_param><BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_query-ID>q_β</Iteration_query-ID><Iteration_query-def>source &amp; description</Iteration_query-def><Iteration_query-len>8</Iteration_query-len><Iteration_hits><Hit><Hit_num>1</Hit_num><Hit_id>subject</Hit_id><Hit_def>generated subject</Hit_def><Hit_accession>subject</Hit_accession><Hit_len>9</Hit_len><Hit_hsps><Hsp>{fields}</Hsp></Hit_hsps></Hit></Iteration_hits></Iteration><Iteration><Iteration_iter-num>2</Iteration_iter-num><Iteration_query-ID>empty</Iteration_query-ID><Iteration_query-len>8</Iteration_query-len><Iteration_hits/><Iteration_message>No hits found</Iteration_message></Iteration></BlastOutput_iterations></BlastOutput>"
        )
    } else {
        format!(
            "<?xml version=\"1.0\"?><BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><BlastOutput2><report><Report><program>blastn</program><version>synthetic</version><reference>generated &amp; local</reference><search-target><Target><db>db</db></Target></search-target><params><Parameters><expect>10</expect></Parameters></params><results><Results><search><Search><query-id>q_β</query-id><query-title>source &amp; description</query-title><query-len>8</query-len><hits><Hit><num>1</num><description><HitDescr><id>subject</id><accession>subject</accession><title>generated subject</title></HitDescr></description><len>9</len><hsps><Hsp>{fields}</Hsp></hsps></Hit></hits></Search></search></Results></results></Report></report></BlastOutput2></BlastXML2>"
        )
    }
}

fn xml_options(format: BlastXmlFormat) -> BlastXmlOptions {
    BlastXmlOptions {
        format,
        limits: blast_limits(),
        max_depth: 32,
        max_event_bytes: 4096,
        max_metadata_bytes: 4096,
    }
}

fn xml_stream(
    input: &[u8],
    format: BlastXmlFormat,
    schedule: Schedule,
) -> StreamOutcome<BlastQuery> {
    let (reader, calls) = buffered(input, schedule);
    drain(
        BlastXmlReader::new(reader, xml_options(format)).unwrap(),
        input.len(),
        &calls,
    )
}

fn assert_blast_model(document: &BlastDocument, queries: usize) {
    assert_eq!(document.queries().len(), queries);
    let query = &document.queries()[0];
    assert_eq!(query.id(), "q_β");
    assert_eq!(query.hits()[0].id(), "subject");
    let hsp = &query.hits()[0].hsps()[0];
    assert_eq!(hsp.query_coordinates().unwrap().range(), 0..4);
    let hit = hsp.hit_coordinates().unwrap();
    assert_eq!(
        (hit.reported_start(), hit.reported_end(), hit.range()),
        (9, 6, 5..9)
    );
    assert_eq!(hsp.evalue(), Some(0.0001));
}

#[test]
fn blast_tabular_mutations_keep_source_lexemes_queries_coordinates_and_fused_errors() {
    for commented in [false, true] {
        let input = tabular_input(commented);
        let document = parse_blast_tabular(input.as_bytes(), tabular_options(commented)).unwrap();
        assert_blast_model(&document, if commented { 2 } else { 1 });
        assert_eq!(
            document.queries()[0].hits()[0].hsps()[0].field("evalue"),
            Some("1.00e-4")
        );
        if commented {
            assert!(document.queries()[1].hits().is_empty());
            assert_eq!(document.processed_queries(), Some(2));
        }
        for bytes in mutations(input.as_bytes(), 0xa319_827d)
            .into_iter()
            .chain(random_corpus())
        {
            let expected = tabular_stream(&bytes, commented, Schedule::chunks(64));
            for chunk in CHUNKS {
                assert_eq!(
                    tabular_stream(&bytes, commented, Schedule::chunks(chunk)),
                    expected
                );
            }
        }
    }
}

#[test]
fn blast_xml_mutations_keep_namespace_entities_and_whole_query_models_chunk_invariant() {
    for format in [BlastXmlFormat::Legacy, BlastXmlFormat::Xml2] {
        let input = xml_input(format);
        let document = parse_blast_xml(input.as_bytes(), xml_options(format)).unwrap();
        assert_blast_model(
            &document,
            if format == BlastXmlFormat::Legacy {
                2
            } else {
                1
            },
        );
        assert_eq!(
            document.queries()[0].description(),
            Some("source & description")
        );
        for bytes in mutations(input.as_bytes(), 0xf713_c4a9)
            .into_iter()
            .chain(random_corpus())
        {
            let expected = xml_stream(&bytes, format, Schedule::chunks(64));
            for chunk in CHUNKS {
                assert_eq!(
                    xml_stream(&bytes, format, Schedule::chunks(chunk)),
                    expected,
                    "format {format:?}, chunk {chunk}"
                );
            }
        }
        if format == BlastXmlFormat::Xml2 {
            let attributes = concat!(
                " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"",
                " xsi:schemaLocation=\"schema\tline\nnext&#x9;part&#xA;&amp;\""
            );
            let input = xml_root_attributes(&input, format, attributes);
            let expected = parse_blast_xml(input.as_bytes(), xml_options(format)).unwrap();
            let source = expected
                .source_fields()
                .iter()
                .find(|field| field.name().ends_with("@xsi:schemaLocation"))
                .unwrap();
            assert_eq!(source.value(), "schema\tline\nnext\tpart\n&");
            assert_blast_model(&expected, 1);
            for chunk in CHUNKS {
                let (reader, _) = buffered(input.as_bytes(), Schedule::chunks(chunk));
                assert_eq!(
                    parse_blast_xml(reader, xml_options(format)).unwrap(),
                    expected
                );
            }
        }
    }
}

#[test]
fn blast_transient_and_permanent_transport_failures_have_explicit_outcomes() {
    let tabular = tabular_input(false);
    let expected = tabular_stream(tabular.as_bytes(), false, Schedule::chunks(64));
    for position in 0..=tabular.len() {
        assert_eq!(
            tabular_stream(tabular.as_bytes(), false, Schedule::interrupted(position)),
            expected,
            "tabular interrupt {position}"
        );
        assert_eq!(
            tabular_stream(tabular.as_bytes(), false, Schedule::failing(position))
                .failure
                .unwrap()
                .kind,
            IoErrorKind::Io
        );
    }
    for format in [BlastXmlFormat::Legacy, BlastXmlFormat::Xml2] {
        let input = xml_input(format);
        let expected = xml_stream(input.as_bytes(), format, Schedule::chunks(64));
        // Uniform one-byte chunks cover every byte seam. Transport faults use
        // 17 evenly spread positions plus EOF, not a quadratic large campaign.
        for position in (0..=16).map(|step| input.len() * step / 16) {
            assert_eq!(
                xml_stream(input.as_bytes(), format, Schedule::interrupted(position)),
                expected,
                "XML {format:?}, interrupt {position}"
            );
            assert_eq!(
                xml_stream(input.as_bytes(), format, Schedule::failing(position))
                    .failure
                    .unwrap()
                    .kind,
                IoErrorKind::Io
            );
        }
    }
}

fn xml_root_attributes(input: &str, format: BlastXmlFormat, attributes: &str) -> String {
    let start = if format == BlastXmlFormat::Legacy {
        "<BlastOutput>"
    } else {
        "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\">"
    };
    let root = start.strip_suffix('>').unwrap();
    input.replacen(start, &format!("{root}{attributes}>"), 1)
}

fn namespace_attribute_probe() {
    for format in [BlastXmlFormat::Legacy, BlastXmlFormat::Xml2] {
        let baseline = xml_input(format);
        let expected = parse_blast_xml(baseline.as_bytes(), xml_options(format)).unwrap();
        // XML2 already has one default namespace declaration on its root.
        let allowance = if format == BlastXmlFormat::Legacy {
            256
        } else {
            255
        };
        let mut options = xml_options(format);
        options.limits.read.max_line_bytes = 8192;
        options.max_event_bytes = 8192;
        options.max_metadata_bytes = 8192;
        for count in [allowance, allowance + 1] {
            let mut attributes = String::new();
            for index in 0..count {
                write!(attributes, " xmlns:n{index}=\"u\"").unwrap();
            }
            let input = xml_root_attributes(&baseline, format, &attributes);
            assert!(input.len() < robust_support::MAX_INPUT_BYTES);
            let (reader, _) = buffered(input.as_bytes(), Schedule::chunks(1));
            let result = parse_blast_xml(reader, options);
            if count == allowance {
                assert_eq!(result.unwrap().queries(), expected.queries());
            } else {
                assert_eq!(result.unwrap_err().kind(), IoErrorKind::InvalidStructure);
            }
        }
        let duplicate = xml_root_attributes(&baseline, format, " xmlns:n=\"u\" xmlns:n=\"v\"");
        assert_eq!(
            parse_blast_xml(duplicate.as_bytes(), options)
                .unwrap_err()
                .kind(),
            IoErrorKind::InvalidStructure
        );
        let mut attributes = String::new();
        for index in 0..512 {
            write!(attributes, " a{index}=\"x\"").unwrap();
        }
        let input = xml_root_attributes(&baseline, format, &attributes);
        assert!(input.len() < robust_support::MAX_INPUT_BYTES);
        assert_eq!(
            parse_blast_xml(input.as_bytes(), options)
                .unwrap_err()
                .kind(),
            IoErrorKind::Unsupported
        );
    }
}

#[test]
fn blast_xml_namespace_and_attribute_limits_complete_under_a_subprocess_watchdog() {
    const PROBE_ENV: &str = "BIOFERRIS_ROBUST_XML_NAMESPACE_PROBE";
    if std::env::var(PROBE_ENV).as_deref() == Ok("bounded-v1") {
        namespace_attribute_probe();
        return;
    }
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "blast_xml_namespace_and_attribute_limits_complete_under_a_subprocess_watchdog",
            "--nocapture",
        ])
        .env(PROBE_ENV, "bounded-v1")
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let started = Instant::now();
    loop {
        if child.try_wait().unwrap().is_some() {
            let result = child.wait_with_output().unwrap();
            assert!(
                result.status.success(),
                "namespace probe failed: stdout={} stderr={}",
                String::from_utf8_lossy(&result.stdout),
                String::from_utf8_lossy(&result.stderr)
            );
            break;
        }
        if started.elapsed() >= Duration::from_secs(10) {
            child.kill().unwrap();
            let result = child.wait_with_output().unwrap();
            panic!("bounded namespace/attribute probe exceeded 10-second watchdog: {result:?}");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

fn decoded_gzip(bytes: &[u8], schedule: Schedule) -> std::io::Result<Vec<u8>> {
    let source = AdversarialRead::new(bytes, schedule);
    let mut reader = decode_input(source).map_err(std::io::Error::other)?;
    let mut output = Vec::new();
    reader
        .by_ref()
        .take(u64::try_from(MAX_DECODED + 1).unwrap())
        .read_to_end(&mut output)?;
    assert!(
        output.len() <= MAX_DECODED,
        "decoded corpus budget exceeded"
    );
    Ok(output)
}

fn decoded_bgzf(bytes: &[u8], schedule: Schedule) -> (Vec<u8>, Option<BgzfError>) {
    let mut reader = BgzfReader::new(AdversarialRead::new(bytes, schedule));
    let mut output = Vec::new();
    let mut buffer = [0xa5; 19];
    for _ in 0..=MAX_DECODED {
        buffer.fill(0xa5);
        match reader.read_bytes(&mut buffer) {
            Ok(0) => return (output, None),
            Ok(count) => {
                output.extend_from_slice(&buffer[..count]);
                assert!(output.len() <= MAX_DECODED);
            }
            Err(error) => {
                assert_eq!(buffer, [0xa5; 19], "failed member modified caller buffer");
                assert!(matches!(
                    reader.read_bytes(&mut buffer),
                    Err(BgzfError::FailedState)
                ));
                return (output, Some(error));
            }
        }
    }
    panic!("BGZF read-call budget exceeded");
}

#[test]
fn independent_gzip_and_bgzf_frames_survive_all_short_reads_and_interruptions() {
    assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
    let mut generator = Generator::new(0xc621_7b3a);
    let first = generator.bytes(37);
    let second = generator.bytes(53);
    let expected = [first.as_slice(), second.as_slice()].concat();
    let gzip = [gzip_member(&first), gzip_member(b""), gzip_member(&second)].concat();
    let bgzf = [
        bgzf_member(&first),
        BGZF_EOF_MARKER.to_vec(),
        bgzf_member(&second),
        BGZF_EOF_MARKER.to_vec(),
    ]
    .concat();
    for chunk in CHUNKS {
        assert_eq!(
            decoded_gzip(&gzip, Schedule::chunks(chunk)).unwrap(),
            expected
        );
        let (observed, failure) = decoded_bgzf(&bgzf, Schedule::chunks(chunk));
        assert!(failure.is_none(), "{failure:?}");
        assert_eq!(observed, expected);
    }
    for position in 0..=gzip.len() {
        assert_eq!(
            decoded_gzip(&gzip, Schedule::interrupted(position)).unwrap(),
            expected
        );
    }
    for position in 0..=bgzf.len() {
        let (observed, failure) = decoded_bgzf(&bgzf, Schedule::interrupted(position));
        assert!(failure.is_none(), "position {position}: {failure:?}");
        assert_eq!(observed, expected);
        let (_, failure) = decoded_bgzf(&bgzf, Schedule::failing(position));
        assert!(
            matches!(failure, Some(BgzfError::Io { .. })),
            "position {position}: {failure:?}"
        );
    }
}

#[test]
fn corrupt_or_truncated_bgzf_members_never_release_their_unverified_payload() {
    let safe = bgzf_member(b"safe prefix:");
    let secret = b"unverified payload must stay hidden";
    let member = bgzf_member(secret);
    for end in 1..member.len() {
        let input = [safe.as_slice(), &member[..end]].concat();
        let (observed, failure) = decoded_bgzf(&input, Schedule::chunks(1));
        assert!(failure.is_some(), "truncated member end {end}");
        assert_eq!(observed, b"safe prefix:");
    }
    let mut rejected = 0;
    for position in 0..member.len() {
        let mut mutated = member.clone();
        mutated[position] ^= 0x80;
        let input = [safe.as_slice(), mutated.as_slice(), &BGZF_EOF_MARKER].concat();
        for chunk in CHUNKS {
            let (observed, failure) = decoded_bgzf(&input, Schedule::chunks(chunk));
            if failure.is_some() {
                rejected += 1;
                assert_eq!(
                    observed, b"safe prefix:",
                    "position {position}, chunk {chunk}"
                );
            } else {
                // Some header metadata/padding edits remain valid; data must not change.
                assert_eq!(observed, [b"safe prefix:".as_slice(), secret].concat());
            }
        }
    }
    assert!(
        rejected > member.len(),
        "mutation corpus did not exercise rejection"
    );
}

#[test]
fn compression_writers_round_trip_seeded_bytes_for_distinct_write_schedules() {
    let mut generator = Generator::new(0x78ba_301f);
    for length in [0, 1, 7, 31, 255, 1024, 4096] {
        let bytes = generator.bytes(length);
        for chunk in CHUNKS {
            let mut gzip = GzipWriter::new(Vec::new());
            let mut bgzf = BgzfWriter::with_options(
                Vec::new(),
                BgzfWriteOptions {
                    block_size: 113,
                    ..BgzfWriteOptions::default()
                },
            )
            .unwrap();
            for piece in bytes.chunks(chunk) {
                gzip.write_all(piece).unwrap();
                bgzf.write_all(piece).unwrap();
            }
            let gzip = gzip.finish().unwrap();
            let bgzf = bgzf.finish().unwrap();
            assert_eq!(decoded_gzip(&gzip, Schedule::chunks(3)).unwrap(), bytes);
            let (decoded, failure) = decoded_bgzf(&bgzf, Schedule::chunks(3));
            assert!(failure.is_none(), "{failure:?}");
            assert_eq!(decoded, bytes);
        }
    }
}

#[test]
fn truncated_gzip_and_bad_trailers_are_not_successful_complete_imports() {
    let member = gzip_member(FASTA);
    for end in 2..member.len() {
        assert!(
            decoded_gzip(&member[..end], Schedule::chunks(1)).is_err(),
            "truncation {end}"
        );
    }
    for position in member.len() - 8..member.len() {
        let mut corrupted = member.clone();
        corrupted[position] ^= 0x40;
        let source = decode_input(AdversarialRead::new(&corrupted, Schedule::chunks(1))).unwrap();
        let mut iterator = FastaReader::<_, Dna>::new(source);
        let mut failure = None;
        for record in iterator.by_ref() {
            if let Err(error) = record {
                failure = Some(error.kind());
                break;
            }
        }
        assert_eq!(failure, Some(IoErrorKind::Io));
        assert!(iterator.next().is_none());
    }
}

fn documentation_fence<'a>(guide: &'a str, name: &str, language: &str) -> &'a str {
    guide
        .split_once(&format!("<!-- robustness-{name} -->"))
        .unwrap()
        .1
        .split_once(&format!("```{language}\n"))
        .unwrap()
        .1
        .split_once("```")
        .unwrap()
        .0
}

#[test]
fn documented_robustness_recipe() {
    let guide = include_str!("../../docs/testing/robustness.md");
    let input = b"@sample\nAcGT\n+\n!I~+\n";
    let malformed = b"@sample\nAcGT\n+\n!I ~\n";
    assert_eq!(
        documentation_fence(guide, "input", "fastq").as_bytes(),
        input
    );
    assert_eq!(
        documentation_fence(guide, "invalid", "fastq").as_bytes(),
        malformed
    );
    assert_eq!(
        documentation_fence(guide, "command", "bash").trim(),
        "cargo +1.85.0 test -p bioferris --test robust_streaming --offline --locked --jobs 2 documented_robustness_recipe -- --exact --nocapture"
    );
    let expected = Record::<Dna>::from_ascii("sample", "AcGT", Some(vec![0, 40, 93, 10])).unwrap();
    for chunk in CHUNKS {
        let parsed = sequence(input, true, Schedule::chunks(chunk));
        assert_eq!(parsed.values, [expected.clone()]);
        assert!(parsed.failure.is_none());
    }
    let invalid = sequence(malformed, true, Schedule::chunks(1));
    assert!(invalid.values.is_empty());
    let kind = invalid.failure.unwrap().kind;
    assert_eq!(kind, IoErrorKind::InvalidQuality);
    let output = format!(
        "records: 1\nid: {}\nsource: {}\nnormalized: {}\nphred: {:?}\nchunk schedules: {}\nfused after invalid quality: true\ninvalid quality code: {}\n",
        expected.id(),
        expected.source_sequence(),
        expected.sequence(),
        expected.qualities().unwrap(),
        CHUNKS.len(),
        kind.code()
    );
    assert_eq!(output, documentation_fence(guide, "output", "text"));
    print!("{output}");
}
