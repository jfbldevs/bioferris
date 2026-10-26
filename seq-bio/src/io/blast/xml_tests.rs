use std::io::{self, BufRead, Cursor, Read};

use super::*;
use crate::io::IoErrorKind;
use serde_json::{Value, json};

const NUCLEOTIDE: &str = include_str!("../../../tests/fixtures/blast-ncbi-2.9.0-nucleotide.xml");
const NUCLEOTIDE_XML2: &str =
    include_str!("../../../tests/fixtures/blast-ncbi-2.9.0-nucleotide-xml2.xml");
const TRANSLATED: &str = include_str!("../../../tests/fixtures/blast-ncbi-2.9.0-translated.xml");
const TRANSLATED_XML2: &str =
    include_str!("../../../tests/fixtures/blast-ncbi-2.9.0-translated-xml2.xml");

fn parse(text: &str, format: BlastXmlFormat) -> Result<BlastDocument, BioIoError> {
    parse_blast_xml(
        Cursor::new(text),
        BlastXmlOptions {
            format,
            ..Default::default()
        },
    )
}

fn legacy_query(id: &str) -> String {
    format!(
        "<Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_query-ID>{id}</Iteration_query-ID><Iteration_query-def>none</Iteration_query-def><Iteration_query-len>4</Iteration_query-len><Iteration_hits/><Iteration_message>No hits found</Iteration_message></Iteration>"
    )
}

fn legacy(queries: &str) -> String {
    format!(
        "<?xml version=\"1.0\"?><BlastOutput><BlastOutput_program>blastn</BlastOutput_program><BlastOutput_version>BLASTN synthetic</BlastOutput_version><BlastOutput_reference>A &amp; B</BlastOutput_reference><BlastOutput_db>db</BlastOutput_db><BlastOutput_query-ID>root</BlastOutput_query-ID><BlastOutput_query-def>root title</BlastOutput_query-def><BlastOutput_query-len>4</BlastOutput_query-len><BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect><Parameters_gap-open>0</Parameters_gap-open><Parameters_gap-extend>0</Parameters_gap-extend></Parameters></BlastOutput_param><BlastOutput_iterations>{queries}</BlastOutput_iterations></BlastOutput>"
    )
}

fn xml2(search: &str) -> String {
    format!(
        "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><BlastOutput2><report><Report><program>blastn</program><version>synthetic</version><reference>public format</reference><search-target><Target><db>db</db></Target></search-target><params><Parameters><expect>10</expect></Parameters></params><results><Results><search>{search}</search></Results></results></Report></report></BlastOutput2></BlastXML2>"
    )
}

#[test]
fn public_nucleotide_reports_preserve_original_and_normalized_coordinates() {
    for (input, format) in [
        (NUCLEOTIDE, BlastXmlFormat::Legacy),
        (NUCLEOTIDE_XML2, BlastXmlFormat::Xml2),
    ] {
        let document = parse(input, format).unwrap();
        let query = &document.queries()[0];
        assert_eq!(document.queries().len(), 1);
        assert_eq!(query.id(), "G26684.1");
        assert_eq!(query.length(), Some(285));
        assert_eq!(query.program(), Some("blastn"));
        assert_eq!(query.version(), Some("BLASTN 2.9.0+"));
        assert_eq!(query.hits().len(), 10);
        let first = &query.hits()[0].hsps()[0];
        assert_eq!(first.query_coordinates().unwrap().range(), 133..166);
        assert_eq!(
            first.hit_coordinates().unwrap().range(),
            101_449_143..101_449_177
        );
        assert_eq!(
            first.hit_coordinates().unwrap().reported_start(),
            101_449_177
        );
        assert_eq!(first.hit_strand(), Some(BlastStrand::Minus));
        assert_eq!(first.bit_score(), Some(40.9604));
        assert_eq!(first.evalue(), Some(0.375_311));
        assert_eq!(first.alignment_length(), Some(34));
        assert_eq!(first.midline(), Some("|||||||||||||| || |||||| || ||||||"));
        assert!(
            query
                .report_fields()
                .iter()
                .any(|field| field.value().contains("Sch&auml;ffer"))
        );
        assert_eq!(first.mismatches(), None);
        assert_eq!(first.identity_percent(), None);
    }
}

#[test]
fn public_translated_reports_keep_source_units_and_explicit_frames() {
    let one = parse(TRANSLATED, BlastXmlFormat::Legacy).unwrap();
    let two = parse(TRANSLATED_XML2, BlastXmlFormat::Xml2).unwrap();
    for document in [&one, &two] {
        let query = &document.queries()[0];
        assert_eq!(query.id(), "AI021773.1");
        assert_eq!(query.program(), Some("blastx"));
        assert_eq!(query.hits().len(), 10);
        for hit in query.hits() {
            for hsp in hit.hsps() {
                assert!(hsp.query_frame().is_some());
                assert!(hsp.query_strand().is_some());
                assert_eq!(
                    hsp.query_alignment().unwrap().len(),
                    hsp.hit_alignment().unwrap().len()
                );
            }
        }
    }
    for (a, b) in one.queries()[0].hits().iter().zip(two.queries()[0].hits()) {
        assert_eq!(a.id(), b.id());
        assert_eq!(a.length(), b.length());
        assert_eq!(a.hsps().len(), b.hsps().len());
        for (a, b) in a.hsps().iter().zip(b.hsps()) {
            assert_eq!(a.query_coordinates(), b.query_coordinates());
            assert_eq!(a.hit_coordinates(), b.hit_coordinates());
            assert_eq!(a.query_frame(), b.query_frame());
            assert_eq!(a.evalue(), b.evalue());
            assert_eq!(a.bit_score(), b.bit_score());
            assert_eq!(a.query_alignment(), b.query_alignment());
        }
    }
}

#[test]
fn multiple_and_no_hit_queries_share_report_metadata_without_identity_merging() {
    let input = legacy(&(legacy_query("one") + &legacy_query("two")));
    let document = parse(&input, BlastXmlFormat::Legacy).unwrap();
    assert_eq!(document.queries().len(), 2);
    assert!(
        document
            .queries()
            .iter()
            .all(|query| query.hits().is_empty())
    );
    assert_eq!(document.queries()[0].reported_id(), Some("one"));
    assert_eq!(document.queries()[1].reported_id(), Some("two"));
    assert!(std::ptr::eq(
        document.queries()[0].report_fields(),
        document.queries()[1].report_fields()
    ));
    assert!(
        document.queries()[0]
            .source_fields()
            .iter()
            .any(|field| field.name() == "Iteration_message" && field.value() == "No hits found")
    );
    assert!(
        parse(&legacy(""), BlastXmlFormat::Legacy)
            .unwrap()
            .queries()
            .is_empty()
    );
}

#[test]
fn xml2_optional_fields_no_hits_and_multiple_reports_remain_distinct() {
    let first =
        xml2("<Search><query-id>one</query-id><hits/><message>No hits found</message></Search>");
    let second = first.replace("<query-id>one", "<query-id>two");
    let combined = first.replace(
        "</BlastXML2>",
        &second.replace("<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\">", ""),
    );
    let document = parse(&combined, BlastXmlFormat::Xml2).unwrap();
    assert_eq!(document.queries().len(), 2);
    for query in document.queries() {
        assert_eq!(query.length(), None);
        assert_eq!(query.description(), None);
        assert!(query.hits().is_empty());
    }
    assert!(!std::ptr::eq(
        document.queries()[0].report_fields(),
        document.queries()[1].report_fields()
    ));
}

#[test]
fn xml2_multiple_hit_descriptions_preserve_order_and_taxonomy() {
    let extra = "</HitDescr><HitDescr><id>second</id><title>Another &amp; title</title><taxid>42</taxid><sciname>species</sciname></HitDescr>";
    let input = NUCLEOTIDE_XML2.replacen("</HitDescr>", extra, 1);
    let document = parse(&input, BlastXmlFormat::Xml2).unwrap();
    let hit = &document.queries()[0].hits()[0];
    assert_eq!(hit.id(), "gi|372099107|ref|NC_000069.6|");
    assert!(
        hit.source_fields()
            .iter()
            .any(|field| field.name() == "description/HitDescr[1]/title"
                && field.value() == "Another & title")
    );
    assert!(
        hit.source_fields()
            .iter()
            .any(|field| field.name() == "description/HitDescr[1]/taxid" && field.value() == "42")
    );
}

#[test]
fn encodings_namespaces_unknown_semantics_and_iterative_modes_are_explicit_errors() {
    let source = legacy(&legacy_query("one"));
    let unsupported = [
        source.replace("version=\"1.0\"", "version=\"1.1\""),
        source.replace("version=\"1.0\"", "version=\"1.0\" encoding=\"UTF-16\""),
        source.replace("<BlastOutput>", "<BlastOutput xmlns=\"urn:wrong\">"),
        source.replace(
            "<BlastOutput_program>blastn",
            "<BlastOutput_program>psiblast",
        ),
        source.replace("<Iteration_hits/>", "<unknown/>"),
        source.replace("A &amp; B", "A &invented; B"),
        source.replace(
            "<BlastOutput>",
            "<!DOCTYPE BlastOutput [<!ENTITY x 'bad'>]><BlastOutput>",
        ),
        source.replace("<BlastOutput>", "<?process value?><BlastOutput>"),
    ];
    for input in unsupported {
        assert_eq!(
            parse(&input, BlastXmlFormat::Legacy).unwrap_err().kind(),
            IoErrorKind::Unsupported,
            "{input}"
        );
    }
    for container in ["iterations", "bl2seq"] {
        let input = xml2("<Search><query-id>one</query-id></Search>")
            .replace("<search>", &format!("<{container}>"))
            .replace("</search>", &format!("</{container}>"));
        assert_eq!(
            parse(&input, BlastXmlFormat::Xml2).unwrap_err().kind(),
            IoErrorKind::Unsupported
        );
    }
    assert!(parse(NUCLEOTIDE, BlastXmlFormat::Xml2).is_err());
    assert!(parse(NUCLEOTIDE_XML2, BlastXmlFormat::Legacy).is_err());
}

#[test]
fn malformed_nested_attributes_entities_counts_and_coordinates_fail() {
    let source = legacy(&legacy_query("one"));
    for input in [
        source.replace("<BlastOutput>", "<BlastOutput x=\"1\">"),
        source.replace("<BlastOutput>", "<BlastOutput xmlns=\"\" xmlns=\"\">"),
        source.replace("</Iteration>", "</broken>"),
        source.replace(
            "<Iteration_hits/>",
            "<Iteration_hits><broken></Iteration_hits>",
        ),
        source.replace("A &amp; B", "A &#0; B"),
        source.replace("A &amp; B", "A \u{0} B"),
        source.replace("<Iteration_query-len>4", "<Iteration_query-len>0"),
        source.replace("<Iteration_hits/>", "<Iteration_hits/><Iteration_hits/>"),
        NUCLEOTIDE.replace("<Hsp_evalue>0.375311", "<Hsp_evalue>NaN"),
        NUCLEOTIDE.replace("<Hsp_hit-from>101449177", "<Hsp_hit-from>0"),
        NUCLEOTIDE.replace("<Hsp_identity>30", "<Hsp_identity>999"),
        NUCLEOTIDE.replace("<Hsp_hit-frame>-1", "<Hsp_hit-frame>1"),
    ] {
        assert!(
            parse(&input, BlastXmlFormat::Legacy).is_err(),
            "accepted malformed input: {input}"
        );
    }
}

#[test]
fn stream_returns_complete_query_then_reports_invalid_suffix_and_fuses() {
    let source = legacy(&(legacy_query("one") + &legacy_query("two")));
    let truncated = &source[..source.len() - 5];
    let mut reader =
        BlastXmlReader::new(Cursor::new(truncated), BlastXmlOptions::default()).unwrap();
    assert_eq!(reader.next().unwrap().unwrap().id(), "one");
    assert_eq!(reader.next().unwrap().unwrap().id(), "two");
    let error = reader.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::Truncated);
    assert!(error.line().is_some());
    assert!(error.record().is_some());
    assert!(reader.next().is_none());
    assert!(reader.next().is_none());
    assert!(parse(truncated, BlastXmlFormat::Legacy).is_err());
    assert!(parse(&(source + "<extra/>"), BlastXmlFormat::Legacy).is_err());
}

#[test]
fn every_xml_resource_limit_has_a_checked_boundary() {
    let source = legacy(&legacy_query("one"));
    for resource in 0..8 {
        let mut options = BlastXmlOptions::default();
        match resource {
            0 => options.limits.max_queries = 0,
            1 => options.limits.max_query_bytes = 1,
            2 => options.limits.max_query_id_bytes = 0,
            3 => options.limits.read.max_total_bytes = source.len() - 1,
            4 => options.limits.read.max_line_bytes = source.len() - 1,
            5 => options.max_depth = 1,
            6 => options.max_event_bytes = 1,
            7 => options.max_metadata_bytes = 1,
            _ => unreachable!(),
        }
        assert_eq!(
            parse_blast_xml(Cursor::new(&source), options)
                .unwrap_err()
                .kind(),
            IoErrorKind::LimitExceeded,
            "resource {resource}"
        );
    }
    for resource in 0..2 {
        let mut options = BlastXmlOptions::default();
        if resource == 0 {
            options.limits.max_hits_per_query = 0;
        } else {
            options.limits.max_hsps_per_query = 0;
        }
        assert_eq!(
            parse_blast_xml(Cursor::new(NUCLEOTIDE), options)
                .unwrap_err()
                .kind(),
            IoErrorKind::LimitExceeded
        );
    }
    let mut exact = BlastXmlOptions::default();
    exact.limits.read.max_total_bytes = source.len();
    exact.limits.read.max_line_bytes = source.len();
    assert!(parse_blast_xml(Cursor::new(&source), exact).is_ok());
    let empty = "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"/>";
    let mut options = BlastXmlOptions {
        format: BlastXmlFormat::Xml2,
        ..BlastXmlOptions::default()
    };
    options.limits.max_queries = 0;
    assert!(
        parse_blast_xml(Cursor::new(empty), options)
            .unwrap()
            .queries()
            .is_empty()
    );
    options.max_depth = 0;
    assert_eq!(
        parse_blast_xml(Cursor::new(empty), options)
            .unwrap_err()
            .kind(),
        IoErrorKind::LimitExceeded
    );
}

#[test]
fn document_comments_cdata_and_only_predefined_entities_are_preserved() {
    let source = legacy(&legacy_query("one"))
        .replace("A &amp; B", "<![CDATA[A < B]]> &#x1F9EC;")
        .replace("<Iteration_hits/>", "<!--query note--><Iteration_hits/>");
    let document = parse(&(source + "<!--after-->"), BlastXmlFormat::Legacy).unwrap();
    assert!(
        document
            .comments()
            .iter()
            .any(|comment| comment.text() == "after")
    );
    assert_eq!(document.queries()[0].comments()[0].text(), "query note");
    assert!(
        document.queries()[0]
            .report_fields()
            .iter()
            .any(|field| field.value() == "A < B 🧬")
    );
}

#[test]
fn source_error_reports_are_not_empty_successes() {
    let input = "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><BlastOutput2><error><Err><code>12</code><message>bad query</message></Err></error></BlastOutput2></BlastXML2>";
    let failure = parse(input, BlastXmlFormat::Xml2).unwrap_err();
    assert!(failure.message().contains("12: bad query"));
}

#[test]
fn xml_declaration_attributes_and_ascii_source_encoding_are_checked() {
    let source = legacy(&legacy_query("one"));
    for attributes in [
        "version=\"1.0\" extra=\"bad\"",
        "version=\"1.0\" version=\"1.0\"",
        "version=\"1.0\" standalone=\"maybe\"",
        "version=\"1.0\" standalone=\"yes\" encoding=\"UTF-8\"",
    ] {
        let input = source.replace("version=\"1.0\"", attributes);
        assert!(
            parse(&input, BlastXmlFormat::Legacy).is_err(),
            "{attributes}"
        );
    }
    let ascii = source.replace(
        "version=\"1.0\"",
        "version=\"1.0\" encoding=\"US-ASCII\" standalone=\"yes\"",
    );
    assert!(parse(&ascii, BlastXmlFormat::Legacy).is_ok());
    assert!(parse(&ascii.replace("A &amp; B", "A β B"), BlastXmlFormat::Legacy).is_err());
    assert!(
        parse(
            &ascii.replace("A &amp; B", "A &#x3B2; B"),
            BlastXmlFormat::Legacy
        )
        .is_ok()
    );
}

#[test]
fn query_and_subject_limits_and_explicit_coordinate_bounds_are_checked() {
    for input in [
        NUCLEOTIDE.replace("<Iteration_query-len>285", "<Iteration_query-len>1"),
        NUCLEOTIDE.replace("<Hit_len>160039680", "<Hit_len>10"),
        NUCLEOTIDE.replace(
            "<Hsp_query-from>134",
            "<Hsp_query-from>18446744073709551616",
        ),
        NUCLEOTIDE_XML2.replace("<query-strand>Plus", "<query-strand>Other"),
    ] {
        let format = if input.contains("<BlastXML2") {
            BlastXmlFormat::Xml2
        } else {
            BlastXmlFormat::Legacy
        };
        assert!(parse(&input, format).is_err());
    }
    let repeated = legacy(&(legacy_query("same") + &legacy_query("same")));
    assert_eq!(
        parse(&repeated, BlastXmlFormat::Legacy).unwrap_err().kind(),
        IoErrorKind::Unsupported
    );
}

struct FailedInput;
impl Read for FailedInput {
    fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::other("injected XML transport failure"))
    }
}
impl BufRead for FailedInput {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Err(io::Error::other("injected XML transport failure"))
    }
    fn consume(&mut self, _: usize) {}
}

#[test]
fn xml_transport_error_retains_context_and_fuses() {
    let mut reader = BlastXmlReader::new(FailedInput, BlastXmlOptions::default()).unwrap();
    let error = reader.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::Io);
    assert_eq!(error.record(), Some(1));
    assert_eq!(error.line(), Some(1));
    assert!(reader.next().is_none());
}

#[test]
fn retained_biopython_observations_match_all_declared_comparable_fields() {
    let observations: Value = serde_json::from_str(include_str!(
        "../../../tests/fixtures/blast-biopython-1.88.json"
    ))
    .unwrap();
    assert_eq!(observations["version"], "1.88");
    assert_eq!(observations["python"], "3.12.3");
    let mut queries = 0;
    let mut hsps = 0;
    for case in observations["cases"].as_array().unwrap() {
        let file = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(case["file"].as_str().unwrap());
        let format = case["format"].as_str().unwrap();
        let document = match format {
            "blast-xml" | "blast-xml2" => read_blast_xml(
                file,
                BlastXmlOptions {
                    format: if format == "blast-xml" {
                        BlastXmlFormat::Legacy
                    } else {
                        BlastXmlFormat::Xml2
                    },
                    ..BlastXmlOptions::default()
                },
            )
            .unwrap(),
            "blast-tab6" | "blast-tab7" => read_blast_tabular(
                file,
                BlastTabularOptions {
                    format: if format == "blast-tab6" {
                        BlastTabularFormat::Uncommented
                    } else {
                        BlastTabularFormat::Commented
                    },
                    ..BlastTabularOptions::default()
                },
            )
            .unwrap(),
            _ => panic!("unreviewed fixture format"),
        };
        let expected = case["queries"].as_array().unwrap();
        assert_eq!(document.queries().len(), expected.len());
        for (query, expected) in document.queries().iter().zip(expected) {
            queries += 1;
            assert_eq!(query.id(), expected["id"].as_str().unwrap());
            if let Some(length) = expected.get("length") {
                assert_eq!(json!(query.length()), *length);
            }
            let expected_hits = expected["hits"].as_array().unwrap();
            assert_eq!(query.hits().len(), expected_hits.len());
            for (hit, expected) in query.hits().iter().zip(expected_hits) {
                assert_eq!(hit.id(), expected["id"].as_str().unwrap());
                if let Some(length) = expected.get("length") {
                    assert_eq!(json!(hit.length()), *length);
                }
                let expected_hsps = expected["hsps"].as_array().unwrap();
                assert_eq!(hit.hsps().len(), expected_hsps.len());
                for (hsp, expected) in hit.hsps().iter().zip(expected_hsps) {
                    hsps += 1;
                    compare_reference_hsp(hsp, expected);
                }
            }
        }
    }
    assert_eq!(queries, 11);
    assert_eq!(hsps, 52);
}

fn compare_reference_hsp(hsp: &BlastHsp, expected: &Value) {
    let query = hsp.query_coordinates().unwrap();
    let hit = hsp.hit_coordinates().unwrap();
    assert_eq!(
        json!([query.range().start, query.range().end]),
        expected["query_range"]
    );
    assert_eq!(
        json!([hit.range().start, hit.range().end]),
        expected["hit_range"]
    );
    if let Some(reported) = expected.get("reported_query") {
        assert_eq!(
            json!([query.reported_start(), query.reported_end()]),
            *reported
        );
    }
    if let Some(reported) = expected.get("reported_hit") {
        assert_eq!(json!([hit.reported_start(), hit.reported_end()]), *reported);
    }
    for (name, value) in [
        ("evalue", json!(hsp.evalue())),
        ("bit_score", json!(hsp.bit_score())),
        ("raw_score", json!(hsp.raw_score())),
        ("alignment_length", json!(hsp.alignment_length())),
        ("identical", json!(hsp.identical())),
        ("positives", json!(hsp.positives())),
        ("gaps", json!(hsp.gaps())),
        ("query_alignment", json!(hsp.query_alignment())),
        ("hit_alignment", json!(hsp.hit_alignment())),
        ("midline", json!(hsp.midline())),
        ("mismatches", json!(hsp.mismatches())),
        ("gap_openings", json!(hsp.gap_openings())),
        ("identity_percent", json!(hsp.identity_percent())),
    ] {
        if let Some(expected) = expected.get(name) {
            if name == "positives" && hsp.positives().is_none() {
                // NCBIXML fills an omitted XML2 positive count from identity.
                // The importer deliberately keeps absent source values absent.
                assert!(
                    !hsp.source_fields()
                        .iter()
                        .any(|field| field.name() == "positive" || field.name() == "Hsp_positive")
                );
                assert_eq!(*expected, json!(hsp.identical()));
                continue;
            }
            assert_eq!(value, *expected, "{name}");
        }
    }
    if let Some(frames) = expected.get("observed_frame") {
        // Biopython may synthesize absent frames. Only explicit source frames
        // participate in equality; absent fields remain absent in BioFerris.
        if let Some(frame) = hsp.query_frame() {
            assert_eq!(json!(frame), frames[0]);
        }
        if let Some(frame) = hsp.hit_frame() {
            assert_eq!(json!(frame), frames[1]);
        }
    }
}
