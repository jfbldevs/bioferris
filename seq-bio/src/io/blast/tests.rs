use std::io::{self, BufRead, Cursor, Read};

use super::*;
use crate::io::IoErrorKind;

const ROW: &str = "q\ts\t100\t4\t0\t0\t1\t4\t8\t5\t1e-20\t42\n";
const HEADER: &str = "# BLASTN 2.15.0+\n# Query: q source query\n# Database: test\n# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n";

fn options(names: &[&str]) -> BlastTabularOptions {
    BlastTabularOptions {
        fields: Some(names.iter().map(|name| (*name).to_owned()).collect()),
        ..BlastTabularOptions::default()
    }
}

fn commented() -> BlastTabularOptions {
    BlastTabularOptions {
        format: BlastTabularFormat::Commented,
        ..BlastTabularOptions::default()
    }
}

fn parse(text: &str, options: BlastTabularOptions) -> Result<BlastDocument, BioIoError> {
    parse_blast_tabular(Cursor::new(text.as_bytes()), options)
}

#[test]
fn ncbi_cookbook_rows_preserve_scores_ids_and_reverse_endpoints() {
    let input = include_str!("../../../tests/fixtures/blast-ncbi-cookbook.tsv");
    let document = parse(input, BlastTabularOptions::default()).unwrap();
    assert_eq!(document.queries().len(), 1);
    let query = &document.queries()[0];
    assert_eq!(query.id(), "AE000111.1");
    assert_eq!(query.identity_field(), "qseqid");
    assert_eq!(query.hits().len(), 3);
    let hsp = &query.hits()[1].hsps()[0];
    assert_eq!(hsp.evalue(), Some(7.69e-30));
    assert_eq!(hsp.field("evalue"), Some("7.69e-30"));
    assert_eq!(hsp.bit_score(), Some(132.0));
    assert_eq!(hsp.query_coordinates().unwrap().range(), 5564..5671);
    let coordinates = hsp.hit_coordinates().unwrap();
    assert_eq!(coordinates.range(), 6820..6928);
    assert_eq!(coordinates.reported_start(), 6928);
    assert_eq!(coordinates.reported_end(), 6821);
    assert_eq!(coordinates.direction(), BlastDirection::Decreasing);
    assert_eq!(hsp.hit_strand(), None);
    assert_eq!(hsp.source_fields().len(), 12);
}

#[test]
fn custom_two_column_results_never_fabricate_scores_or_coordinates() {
    let document = parse("q\ts\n", options(&["qseqid", "sseqid"])).unwrap();
    let hsp = &document.queries()[0].hits()[0].hsps()[0];
    assert_eq!(hsp.query_coordinates(), None);
    assert_eq!(hsp.hit_coordinates(), None);
    assert_eq!(hsp.evalue(), None);
    assert_eq!(hsp.bit_score(), None);
    assert_eq!(hsp.alignment_length(), None);
    assert_eq!(hsp.query_strand(), None);
    assert_eq!(hsp.hit_strand(), None);
}

#[test]
fn unknown_custom_fields_and_numeric_lexemes_survive_verbatim() {
    let document = parse(
        "q\ts\t01.2300e-400\tUnicode β and spaces\n",
        options(&["qseqid", "sseqid", "evalue", "custom_note"]),
    )
    .unwrap();
    let hsp = &document.queries()[0].hits()[0].hsps()[0];
    assert_eq!(hsp.evalue(), Some(0.0));
    assert_eq!(hsp.field("evalue"), Some("01.2300e-400"));
    assert_eq!(hsp.field("custom_note"), Some("Unicode β and spaces"));
    assert_eq!(hsp.source_fields()[3].name(), "custom_note");
}

#[test]
fn query_contiguity_and_subject_first_occurrence_order_are_explicit() {
    let document = parse(
        "q\ts2\nq\ts1\nq\ts2\nr\ts3\n",
        options(&["qseqid", "sseqid"]),
    )
    .unwrap();
    assert_eq!(
        document
            .queries()
            .iter()
            .map(BlastQuery::id)
            .collect::<Vec<_>>(),
        ["q", "r"]
    );
    let hits = document.queries()[0].hits();
    assert_eq!(
        hits.iter().map(BlastHit::id).collect::<Vec<_>>(),
        ["s2", "s1"]
    );
    assert_eq!(hits[0].hsps().len(), 2);
    assert_eq!(hits[0].hsps()[0].line(), 1);
    assert_eq!(hits[0].hsps()[1].line(), 3);
    let mut reader = BlastTabularReader::new(
        Cursor::new("q\ts\nr\ts\nq\ts\n"),
        options(&["qseqid", "sseqid"]),
    )
    .unwrap();
    assert_eq!(reader.next().unwrap().unwrap().id(), "q");
    assert_eq!(reader.next().unwrap().unwrap().id(), "r");
    assert_eq!(
        reader.next().unwrap().unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
    assert!(reader.next().is_none());
    assert!(reader.next().is_none());
}

#[test]
fn commented_blocks_retain_unknown_comments_metadata_and_footer() {
    let text = format!(
        "{HEADER}# RID: public-id\n# an unrecognized comment β\n# 1 hits found\n{ROW}# BLAST processed 1 queries\n# trailing note\n"
    );
    let document = parse(&text, commented()).unwrap();
    let query = &document.queries()[0];
    assert_eq!(query.reported_id(), Some("q"));
    assert_eq!(query.description(), Some("source query"));
    assert_eq!(query.program(), Some("BLASTN"));
    assert_eq!(query.version(), Some("2.15.0+"));
    assert_eq!(query.target(), Some("test"));
    assert_eq!(query.declared_rows(), Some(1));
    assert!(
        query
            .comments()
            .iter()
            .any(|comment| comment.text() == "# an unrecognized comment β")
    );
    assert!(
        query
            .comments()
            .iter()
            .any(|comment| comment.text() == "# RID: public-id")
    );
    assert_eq!(document.processed_queries(), Some(1));
    assert_eq!(
        document
            .comments()
            .iter()
            .map(BlastComment::text)
            .collect::<Vec<_>>(),
        ["# BLAST processed 1 queries", "# trailing note"]
    );
}

#[test]
fn no_hit_and_zero_query_documents_preserve_their_source_information() {
    let text = "# BLASTP 2.15.0+\n# Query: absent description\n# Database: proteins\n# 0 hits found\n# BLAST processed 1 queries\n";
    let document = parse(text, commented()).unwrap();
    let query = &document.queries()[0];
    assert_eq!(query.id(), "absent");
    assert_eq!(query.identity_field(), "query-header");
    assert!(query.hits().is_empty());
    assert_eq!(query.declared_rows(), Some(0));
    assert_eq!(query.fields(), &[] as &[String]);
    let zero = parse("# BLAST processed 0 queries\n", commented()).unwrap();
    assert!(zero.queries().is_empty());
    assert_eq!(zero.processed_queries(), Some(0));
    assert_eq!(zero.comments().len(), 1);
    assert!(
        parse("", BlastTabularOptions::default())
            .unwrap()
            .queries()
            .is_empty()
    );
}

#[test]
fn commented_count_is_hsp_rows_instead_of_distinct_subjects() {
    let document = parse(&format!("{HEADER}# 2 hits found\n{ROW}{ROW}"), commented()).unwrap();
    assert_eq!(document.queries()[0].hits().len(), 1);
    assert_eq!(document.queries()[0].hits()[0].hsps().len(), 2);
}

#[test]
fn accession_identity_can_differ_from_the_comment_query_id() {
    let input = "# BLASTN 2.2.18+\n# Query: gi|1786181|gb|AE000111.1|AE000111\n# Database: ecoli\n# Fields: query acc., subject acc., evalue, q. start, q. end, s. start, s. end\n# 1 hits found\nAE000111\tAE000174\t8e-30\t5565\t5671\t6928\t6821\n";
    let document = parse(input, commented()).unwrap();
    let query = &document.queries()[0];
    assert_eq!(query.id(), "AE000111");
    assert_eq!(query.identity_field(), "qacc");
    assert_eq!(
        query.reported_id(),
        Some("gi|1786181|gb|AE000111.1|AE000111")
    );
    assert_eq!(query.hits()[0].identity_field(), "sacc");
}

#[test]
fn custom_comment_labels_need_explicit_mapping_and_known_order_must_match() {
    let text = "# Query: q\n# Fields: query id, subject id, custom annotation\n# 1 hits found\nq\ts\tvalue\n";
    assert_eq!(
        parse(text, commented()).unwrap_err().kind(),
        IoErrorKind::Unsupported
    );
    let mut configured = options(&["qseqid", "sseqid", "note"]);
    configured.format = BlastTabularFormat::Commented;
    let document = parse(text, configured.clone()).unwrap();
    assert_eq!(
        document.queries()[0].hits()[0].hsps()[0].field("note"),
        Some("value")
    );
    configured.fields = Some(vec![
        "sseqid".to_owned(),
        "qseqid".to_owned(),
        "note".to_owned(),
    ]);
    assert_eq!(
        parse(text, configured).unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
}

#[test]
fn source_frames_strands_and_translated_units_are_preserved() {
    let names = [
        "qseqid", "sseqid", "qlen", "slen", "qstart", "qend", "sstart", "send", "frames", "qframe",
        "sframe", "sstrand", "length", "qseq", "sseq",
    ];
    let text = "q\ts\t1000\t500\t129\t100\t1\t10\t-2/0\t-2\t0\tN/A\t10\tACDEFGHIKL\tAcDEFGHIKL\n";
    let document = parse(text, options(&names)).unwrap();
    let hsp = &document.queries()[0].hits()[0].hsps()[0];
    assert_eq!(hsp.query_coordinates().unwrap().range(), 99..129);
    assert_eq!(hsp.hit_coordinates().unwrap().range(), 0..10);
    assert_eq!(hsp.query_frame(), Some(-2));
    assert_eq!(hsp.hit_frame(), Some(0));
    assert_eq!(hsp.query_strand(), Some(BlastStrand::Minus));
    assert_eq!(hsp.hit_strand(), None);
    assert_eq!(hsp.hit_alignment(), Some("AcDEFGHIKL"));
    assert_eq!(document.queries()[0].length(), Some(1000));
    assert_eq!(document.queries()[0].hits()[0].length(), Some(500));
}

#[test]
fn all_typed_scores_and_counts_remain_optional_and_exactly_reported() {
    let names = [
        "qseqid", "sseqid", "score", "bitscore", "evalue", "length", "nident", "pident",
        "mismatch", "positive", "ppos", "gapopen", "gaps", "qcovs", "qcovhsp", "qcovus",
    ];
    let text = "q\ts\t-2.5\t3.25\t1.000e-100\t10\t5\t50.000\t3\t7\t70\t1\t2\t50\t25\t40\n";
    let document = parse(text, options(&names)).unwrap();
    let hsp = &document.queries()[0].hits()[0].hsps()[0];
    assert_eq!(hsp.raw_score(), Some(-2.5));
    assert_eq!(hsp.bit_score(), Some(3.25));
    assert_eq!(hsp.identical(), Some(5));
    assert_eq!(hsp.identity_percent(), Some(50.0));
    assert_eq!(hsp.mismatches(), Some(3));
    assert_eq!(hsp.positives(), Some(7));
    assert_eq!(hsp.positive_percent(), Some(70.0));
    assert_eq!(hsp.gap_openings(), Some(1));
    assert_eq!(hsp.gaps(), Some(2));
    assert_eq!(hsp.field("evalue"), Some("1.000e-100"));
}

#[test]
fn explicit_delimiter_std_expansion_crlf_and_final_unterminated_line_work() {
    let mut configured = options(&["std", "custom"]);
    configured.delimiter = b'@';
    let text = format!("{}@text", ROW.trim_end().replace('\t', "@"));
    let document = parse(&text, configured).unwrap();
    assert_eq!(
        document.queries()[0].hits()[0].hsps()[0].field("custom"),
        Some("text")
    );
    let document = parse(&ROW.replace('\n', "\r\n"), BlastTabularOptions::default()).unwrap();
    assert_eq!(
        document.queries()[0].hits()[0].hsps()[0].alignment_length(),
        Some(4)
    );
}

#[test]
fn options_reject_duplicate_missing_identity_and_partial_coordinate_fields() {
    for names in [
        vec!["std", "qseqid"],
        vec!["qseqid"],
        vec!["sseqid"],
        vec!["qseqid", "sallseqid"],
        vec!["qseqid", "sseqid", "qstart"],
        vec!["qseqid", "sseqid", "send"],
        vec!["qseqid", "sseqid", "invalid-name"],
    ] {
        assert!(BlastTabularReader::new(Cursor::new(""), options(&names)).is_err());
    }
    for delimiter in [0, b'\n', b'\r', b'#', 255] {
        let configured = BlastTabularOptions {
            delimiter,
            ..BlastTabularOptions::default()
        };
        assert_eq!(
            BlastTabularReader::new(Cursor::new(""), configured)
                .err()
                .unwrap()
                .kind(),
            IoErrorKind::Unsupported
        );
    }
}

#[test]
fn checked_integer_coordinate_extremes_do_not_overflow() {
    let text = format!("q\ts\t{}\t{}\n", u64::MAX, u64::MAX);
    let document = parse(&text, options(&["qseqid", "sseqid", "qstart", "qend"])).unwrap();
    let coordinates = document.queries()[0].hits()[0].hsps()[0]
        .query_coordinates()
        .unwrap();
    assert_eq!(coordinates.range(), (u64::MAX - 1)..u64::MAX);
    assert_eq!(coordinates.direction(), BlastDirection::Point);
    for value in ["0", "-1", "1.5", "18446744073709551616"] {
        assert!(
            parse(
                &format!("q\ts\t{value}\t1\n"),
                options(&["qseqid", "sseqid", "qstart", "qend"])
            )
            .is_err()
        );
    }
}

#[test]
fn invalid_nonfinite_negative_and_out_of_range_numbers_fail() {
    for (name, value) in [
        ("evalue", "NaN"),
        ("evalue", "inf"),
        ("evalue", "-1"),
        ("bitscore", "-1"),
        ("score", "1e400"),
        ("pident", "100.1"),
        ("ppos", "-1"),
        ("qcovs", "101"),
        ("length", "0"),
        ("qlen", "0"),
        ("slen", "0"),
        ("mismatch", "1.5"),
        ("gaps", "-1"),
    ] {
        let failure = parse(
            &format!("q\ts\t{value}\n"),
            options(&["qseqid", "sseqid", name]),
        )
        .unwrap_err();
        assert_eq!(
            failure.kind(),
            IoErrorKind::InvalidStructure,
            "{name} {value}"
        );
        assert_eq!(failure.line(), Some(1));
    }
}

#[test]
fn contradictory_strands_frames_counts_and_lengths_fail() {
    for (names, row) in [
        (
            vec!["qseqid", "sseqid", "qstart", "qend", "qframe"],
            "q\ts\t1\t4\t-1\n",
        ),
        (
            vec!["qseqid", "sseqid", "sstart", "send", "sstrand"],
            "q\ts\t8\t5\tplus\n",
        ),
        (
            vec!["qseqid", "sseqid", "frames", "qframe"],
            "q\ts\t1/-1\t2\n",
        ),
        (
            vec!["qseqid", "sseqid", "sframe", "sstrand"],
            "q\ts\t-2\tplus\n",
        ),
        (
            vec!["qseqid", "sseqid", "qlen", "qstart", "qend"],
            "q\ts\t5\t1\t6\n",
        ),
        (vec!["qseqid", "sseqid", "length", "nident"], "q\ts\t5\t6\n"),
        (
            vec!["qseqid", "sseqid", "positive", "nident"],
            "q\ts\t4\t5\n",
        ),
        (vec!["qseqid", "sseqid", "gaps", "gapopen"], "q\ts\t0\t1\n"),
        (vec!["qseqid", "sseqid", "qseq", "sseq"], "q\ts\tACG\tAC\n"),
        (vec!["qseqid", "sseqid", "length", "qseq"], "q\ts\t3\tAC\n"),
        (vec!["qseqid", "sseqid", "qseq"], "q\ts\tA.C\n"),
        (vec!["qseqid", "sseqid", "qframe"], "q\ts\t4\n"),
        (vec!["qseqid", "sseqid", "frames"], "q\ts\t1/2/3\n"),
    ] {
        assert!(parse(row, options(&names)).is_err(), "{row}");
    }
}

#[test]
fn query_and_subject_metadata_cannot_change_between_rows() {
    for (names, text) in [
        (vec!["qseqid", "sseqid", "qlen"], "q\ts\t10\nq\tt\t11\n"),
        (vec!["qseqid", "sseqid", "slen"], "q\ts\t10\nq\ts\t11\n"),
        (vec!["qseqid", "sseqid", "stitle"], "q\ts\tone\nq\ts\ttwo\n"),
        (vec!["qseqid", "sseqid", "qacc"], "q\ts\tA\nq\ts\tB\n"),
    ] {
        assert_eq!(
            parse(text, options(&names)).unwrap_err().kind(),
            IoErrorKind::InvalidStructure
        );
    }
}

#[test]
fn malformed_comment_order_counts_headers_and_footers_fail() {
    for text in [
        format!("{HEADER}{ROW}"),
        format!("{HEADER}# 0 hits found\n{ROW}"),
        format!("{HEADER}# 2 hits found\n{ROW}"),
        format!("{HEADER}# 1 hits found\n{ROW}# BLAST processed 2 queries\n"),
        format!("{HEADER}# 1 hits found\n{ROW}# 1 hits found\n"),
        format!("{HEADER}# 1 hits found\n{ROW}# Fields: query id, subject id\n"),
        format!("{HEADER}# 1 hits found\n{ROW}# BLAST processed 1 queries\n{ROW}"),
        "# Query: \n# 0 hits found\n".to_owned(),
        "# BLASTN 2.15.0+\n".to_owned(),
        "# Query: q\n# 0 hits found\n# BLAST processed 1 queries\n# BLAST processed 1 queries\n"
            .to_owned(),
        ROW.to_owned(),
    ] {
        assert!(parse(&text, commented()).is_err(), "{text}");
    }
    assert_eq!(
        parse(
            &format!(
                "# Query: another\n{}# 1 hits found\n{ROW}",
                HEADER.lines().skip(2).collect::<Vec<_>>().join("\n") + "\n"
            ),
            commented()
        )
        .unwrap_err()
        .kind(),
        IoErrorKind::InvalidStructure
    );
    assert_eq!(
        parse("# comment\n", BlastTabularOptions::default())
            .unwrap_err()
            .kind(),
        IoErrorKind::Unsupported
    );
}

#[test]
fn malformed_columns_controls_and_utf8_fail_once_with_context() {
    for text in [
        "q\ts\textra\n",
        "\ts\n",
        "query with spaces\ts\n",
        "q\ts\0\n",
    ] {
        let mut reader =
            BlastTabularReader::new(Cursor::new(text), options(&["qseqid", "sseqid"])).unwrap();
        let failure = reader.next().unwrap().unwrap_err();
        assert_eq!(failure.record(), Some(1));
        assert_eq!(failure.line(), Some(1));
        assert!(reader.next().is_none());
    }
    let mut reader =
        BlastTabularReader::new(Cursor::new(b"q\t\xff\n"), options(&["qseqid", "sseqid"])).unwrap();
    assert_eq!(
        reader.next().unwrap().unwrap_err().kind(),
        IoErrorKind::InvalidStructure
    );
    assert!(reader.next().is_none());
}

#[test]
fn every_resource_limit_has_an_observable_boundary() {
    for resource in 0..5 {
        let mut configured = options(&["qseqid", "sseqid"]);
        match resource {
            0 => configured.limits.max_queries = 0,
            1 => configured.limits.max_hits_per_query = 0,
            2 => configured.limits.max_hsps_per_query = 0,
            3 => configured.limits.max_query_bytes = 3,
            4 => configured.limits.max_query_id_bytes = 0,
            _ => unreachable!(),
        }
        assert_eq!(
            parse("q\ts\n", configured).unwrap_err().kind(),
            IoErrorKind::LimitExceeded
        );
    }
    let mut exact = options(&["qseqid", "sseqid"]);
    exact.limits.max_query_bytes = 4;
    exact.limits.max_query_id_bytes = 1;
    exact.limits.max_queries = 1;
    exact.limits.max_hits_per_query = 1;
    exact.limits.max_hsps_per_query = 1;
    assert!(parse("q\ts\n", exact.clone()).is_ok());
    assert_eq!(
        parse("q\ts\nq\ts\n", exact.clone()).unwrap_err().kind(),
        IoErrorKind::LimitExceeded
    );
    exact.limits.read.max_total_bytes = 3;
    assert_eq!(
        parse("q\ts\n", exact).unwrap_err().kind(),
        IoErrorKind::LimitExceeded
    );
}

#[test]
fn declared_and_header_resource_limits_fail_before_rows() {
    let mut configured = commented();
    configured.limits.max_hsps_per_query = 0;
    assert_eq!(
        parse("# Query: q\n# 1 hits found\n", configured)
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
    let mut configured = commented();
    configured.limits.max_query_bytes = 8;
    assert_eq!(
        parse("# very long unknown preamble\n", configured)
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
    let mut configured = options(&["qseqid", "sseqid"]);
    configured.limits.read.max_line_bytes = 32;
    assert_eq!(
        parse(&format!("{}\ts\n", "q".repeat(64)), configured)
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
}

struct FailedInput;
impl Read for FailedInput {
    fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::other("injected transport failure"))
    }
}
impl BufRead for FailedInput {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Err(io::Error::other("injected transport failure"))
    }
    fn consume(&mut self, _: usize) {}
}

#[test]
fn transport_error_is_contextual_and_fuses_iteration() {
    let mut reader = BlastTabularReader::new(FailedInput, BlastTabularOptions::default()).unwrap();
    let failure = reader.next().unwrap().unwrap_err();
    assert_eq!(failure.kind(), IoErrorKind::Io);
    assert_eq!(failure.record(), Some(1));
    assert_eq!(failure.line(), Some(1));
    assert!(reader.next().is_none());
    assert!(reader.next().is_none());
}
