//! Structured JSON and facade file-output contract tests.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use align::output::JSON_SCHEMA_V1;
use align::prelude::*;
use serde_json::Value;

fn dna_records() -> Vec<SequenceRecord<Dna>> {
    [
        ("reference", "ACGT"),
        ("sample_2", "AGT"),
        ("deletion", "ACT"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, DnaSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

fn json_from(write: impl FnOnce(&mut Vec<u8>) -> align::Result<()>) -> (String, Value) {
    let mut output = Vec::new();
    write(&mut output).unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(text.ends_with("}\n"));
    let value = serde_json::from_str(&text).unwrap();
    (text, value)
}

#[test]
fn pairwise_alignment_and_value_json_preserve_objective_semantics() {
    let first = DnaSequence::try_from("GATTACA").unwrap();
    let second = DnaSequence::try_from("GCATGCT").unwrap();
    let run = Pairwise::gotoh().run(&first, &second).unwrap();
    let (text, value) = json_from(|writer| run.write_json_to(writer));
    assert_eq!(value["schema"], JSON_SCHEMA_V1);
    assert_eq!(value["kind"], "pairwise_alignment");
    assert_eq!(value["objective"]["name"], "score");
    assert_eq!(value["ranges"]["first"]["start"], 0);
    assert_eq!(value["ranges"]["second"]["end"], second.len());
    assert_eq!(value["cigar"], run.cigar());
    assert_eq!(value["alignment"]["first"].as_str().unwrap().len(), 7);

    let (again, _) = json_from(|writer| run.write_json_to(writer));
    assert_eq!(text, again);

    let value_run = Pairwise::wfa().run_penalty(&first, &second).unwrap();
    let (_, value) = json_from(|writer| value_run.write_json_to(writer));
    assert_eq!(value["kind"], "pairwise_value");
    assert_eq!(value["objective"]["name"], "penalty");
    assert_eq!(value["objective"]["direction"], "minimize");
    assert_eq!(value["has_traceback"], false);
    assert!(value.get("alignment").is_none());
}

#[test]
fn multiple_json_is_rectangular_and_escapes_identifiers() {
    let identifier = "sample \"one\\two\u{2028}x";
    let records = [
        SequenceRecord::new(identifier, DnaSequence::try_from("ACGT").unwrap()).unwrap(),
        SequenceRecord::new("sample_2", DnaSequence::try_from("AGT").unwrap()).unwrap(),
    ];
    let run = Multiple::mafft(MafftStrategy::FftNs1)
        .run(&records)
        .unwrap();
    let (text, value) = json_from(|writer| run.write_json_to(writer));
    assert_eq!(value["kind"], "multiple_alignment");
    assert_eq!(value["profile"], MafftStrategy::FftNs1.profile_identifier());
    assert_eq!(value["alignment"][0]["identifier"], identifier);
    assert_eq!(value["statistics"]["sequences"], 2);
    assert!(text.contains("\\\"one\\\\two\\u2028x"));
    let width = usize::try_from(value["statistics"]["columns"].as_u64().unwrap()).unwrap();
    assert!(
        value["alignment"]
            .as_array()
            .unwrap()
            .iter()
            .all(|row| row["sequence"].as_str().unwrap().len() == width)
    );
}

#[test]
fn myers_distance_and_search_have_distinct_json_shapes() {
    let pattern = DnaSequence::try_from("GATTACA").unwrap();
    let text = DnaSequence::try_from("TTTGATTACATTT").unwrap();
    let distance = Pairwise::myers().run_distance(&pattern, &text).unwrap();
    let (_, value) = json_from(|writer| distance.write_json_to(writer));
    assert_eq!(value["kind"], "myers_distance");
    assert_eq!(value["distance"], distance.distance());

    let search = Pairwise::myers().run_find_all(&pattern, &text, 2).unwrap();
    let (_, value) = json_from(|writer| search.write_json_to(writer));
    assert_eq!(value["kind"], "myers_search");
    assert_eq!(value["maximum_distance"], 2);
    assert_eq!(
        value["matches"].as_array().unwrap().len(),
        search.matches().len()
    );
    for (index, matched) in search.matches().iter().copied().enumerate() {
        assert_eq!(value["matches"][index]["text_end"], matched.end());
        assert_eq!(value["matches"][index]["distance"], matched.distance());
    }
}

#[test]
fn poasta_query_json_retains_graph_trace_and_search_evidence() {
    let records = dna_records();
    let graph_run = Graph::poasta().run(&records).unwrap();
    let query = DnaSequence::try_from("ACGTT").unwrap();
    let run = Graph::poasta()
        .run_query(graph_run.result().graph(), &query)
        .unwrap();
    let (_, value) = json_from(|writer| run.write_json_to(writer));
    assert_eq!(value["kind"], "poasta_graph_query");
    assert_eq!(
        value["graph"]["nodes"],
        u64::try_from(run.graph().node_count()).unwrap()
    );
    assert_eq!(
        value["operations"].as_array().unwrap().len(),
        run.result().operations().len()
    );
    assert_eq!(
        value["aligned_pairs"].as_array().unwrap().len(),
        run.result().aligned_pairs().len()
    );
    assert!(value["search"]["possible_states"].as_u64().is_some());
}

#[test]
fn run_file_methods_cover_human_machine_phylogenetic_and_graph_outputs() {
    let directory = TemporaryDirectory::new();
    let records = dna_records();
    let msa = Multiple::mafft(MafftStrategy::FftNs1)
        .run(&records)
        .unwrap();

    let paths = [
        directory.path().join("alignment.txt"),
        directory.path().join("alignment.json"),
        directory.path().join("alignment.fasta"),
        directory.path().join("alignment.aln"),
        directory.path().join("alignment.phy"),
        directory.path().join("alignment.nex"),
        directory.path().join("alignment.sto"),
    ];
    msa.write_text(&paths[0]).unwrap();
    msa.write_json(&paths[1]).unwrap();
    msa.write_fasta(&paths[2]).unwrap();
    msa.write_clustal(&paths[3]).unwrap();
    msa.write_phylip(&paths[4]).unwrap();
    msa.write_nexus(&paths[5]).unwrap();
    msa.write_stockholm(&paths[6]).unwrap();

    let strict = PhylipOptions::new().with_name_style(PhylipNameStyle::Strict);
    msa.write_phylip_with(directory.path().join("strict.phy"), strict)
        .unwrap();

    let first = DnaSequence::try_from("ACGT").unwrap();
    let second = DnaSequence::try_from("AGT").unwrap();
    Pairwise::gotoh()
        .run(&first, &second)
        .unwrap()
        .write_fasta(directory.path().join("pairwise.fasta"), "target", "query")
        .unwrap();
    Pairwise::myers()
        .run_find_all(&second, &first, 1)
        .unwrap()
        .write_tsv(directory.path().join("matches.tsv"))
        .unwrap();

    let graph = Graph::poasta().run(&records).unwrap();
    graph.write_gfa(directory.path().join("graph.gfa")).unwrap();
    graph.write_dot(directory.path().join("graph.dot")).unwrap();

    for path in fs::read_dir(directory.path()).unwrap() {
        let path = path.unwrap().path();
        assert!(fs::metadata(&path).unwrap().len() > 0, "{}", path.display());
    }
}

static TEMPORARY_DIRECTORY_ID: AtomicUsize = AtomicUsize::new(0);

struct TemporaryDirectory {
    path: PathBuf,
}

impl TemporaryDirectory {
    fn new() -> Self {
        let id = TEMPORARY_DIRECTORY_ID.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("align-output-tests-{}-{id}", std::process::id()));
        fs::create_dir(&path).unwrap();
        Self { path }
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}
