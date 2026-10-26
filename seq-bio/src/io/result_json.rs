//! Strict import of versioned result documents, not execution requests.
//!
//! Validates schema, shapes, alphabets, coordinates, trace consistency and
//! derived statistics. Imported scores, algorithm labels and execution evidence
//! are assertions in a file, not proof that an algorithm was run or is optimal.
//! CLI envelopes are responses written by the `bioferris` command-line program,
//! which has since been removed; this reader keeps archived responses readable.
//!
//! ```
//! use bioferris::io::result_json::{parse_result_json, ImportedResult, AlignmentResultKind};
//! let input = r#"{"schema":"align-output-v1","kind":"myers_distance",
//! "algorithm":"Myers bit-vector","mode":"global unit edit distance",
//! "alphabet":"DNA","pattern":{"length":4,"sequence":"ACGT"},
//! "text":{"length":3,"sequence":"AGT"},"distance":1}"#;
//! let ImportedResult::Alignment(document) = parse_result_json(input)? else {
//!     panic!("expected a native alignment result");
//! };
//! assert_eq!(document.kind(), AlignmentResultKind::MyersDistance);
//! assert_eq!(document.value()["distance"], 1);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::collections::BTreeSet;
use std::io::BufRead;
use std::path::Path;

use align_core::{Alphabet, Dna, Protein, Rna};
use serde_json::{Map, Value};

use crate::{
    DnaSequenceExt, NucleotideSequenceExt, PartialCodonPolicy, RnaSequenceExt, TranslationMode,
    TranslationOptions,
};

use super::{BioIoError, IoErrorKind, open_input, read_text};

const LIMIT: usize = 64 * 1024 * 1024;

/// Supported result discriminators in the native `align-output-v1` schema.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentResultKind {
    /// Pairwise traceback, ranges, aligned rows and statistics.
    PairwiseAlignment,
    /// A score/distance/penalty without a fabricated traceback.
    PairwiseValue,
    /// Rectangular multiple alignment and metadata.
    MultipleAlignment,
    /// Global unit edit distance, not a traceback.
    MyersDistance,
    /// Ordered approximate-search endpoints, not start coordinates.
    MyersSearch,
    /// Sequence-to-graph operations without the underlying graph topology.
    PoastaGraphQuery,
}

/// Supported biological/error responses in a CLI schema-version 1 envelope.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CliResultKind {
    /// Complement, reverse complement or transcription response.
    SequenceTransform,
    /// Translation with explicit code/frame/mode/partial policy.
    Translation,
    /// Pairwise response containing a native alignment document.
    Pairwise,
    /// Structured error response; not a successful biological result.
    Error,
}

/// An immutable native result whose document-level invariants were checked.
#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentResultDocument {
    kind: AlignmentResultKind,
    value: Value,
}

impl AlignmentResultDocument {
    /// Returns the validated discriminator.
    #[must_use]
    pub const fn kind(&self) -> AlignmentResultKind {
        self.kind
    }

    /// Borrows the complete immutable payload, including uninterpreted metadata.
    #[must_use]
    pub fn value(&self) -> &Value {
        &self.value
    }
}

/// An immutable, validated biological/error CLI response envelope.
#[derive(Clone, Debug, PartialEq)]
pub struct CliResultDocument {
    kind: CliResultKind,
    value: Value,
}

impl CliResultDocument {
    /// Returns the supported response category.
    #[must_use]
    pub const fn kind(&self) -> CliResultKind {
        self.kind
    }

    /// Borrows the complete immutable envelope.
    #[must_use]
    pub fn value(&self) -> &Value {
        &self.value
    }
}

/// Distinguishes a native result from a CLI response; requests are not accepted.
#[derive(Clone, Debug, PartialEq)]
pub enum ImportedResult {
    /// A validated native alignment-result document.
    Alignment(AlignmentResultDocument),
    /// A validated CLI response envelope.
    Cli(CliResultDocument),
}

/// Reads a native/CLI result from a file, including gzip input.
///
/// # Errors
/// Returns I/O, size, schema, unsupported-kind or content-consistency errors.
pub fn read_result_json(path: impl AsRef<Path>) -> Result<ImportedResult, BioIoError> {
    parse_result_json_reader(open_input(path)?)
}

/// Reads one complete result document from a buffered stream, capped at 64 MiB.
///
/// # Errors
/// Rejects trailing JSON documents, duplicate object keys, bad schemas and
/// inconsistent supported payloads. Does not execute a requested operation.
pub fn parse_result_json_reader(reader: impl BufRead) -> Result<ImportedResult, BioIoError> {
    parse_result_json(&read_text(reader.take(LIMIT as u64 + 1), "JSON result")?)
}

/// Parses a supported result document without treating it as an execution proof.
///
/// # Errors
/// Returns [`IoErrorKind::InvalidSchema`] for malformed shapes, versions, keys or
/// inconsistent values; [`IoErrorKind::Unsupported`] for unimplemented kinds or
/// operations; and [`IoErrorKind::LimitExceeded`] above 64 MiB.
pub fn parse_result_json(input: &str) -> Result<ImportedResult, BioIoError> {
    if input.len() > LIMIT {
        return Err(BioIoError::new(
            IoErrorKind::LimitExceeded,
            "JSON result",
            "document exceeds 64 MiB",
        ));
    }
    let value: Value = serde_json::from_str(input)
        .map_err(|error| invalid(error.to_string()).at(1, error.line()))?;
    reject_duplicate_keys(input)?;
    let root = object(&value)?;
    match (
        root.contains_key("schema"),
        root.contains_key("schema_version"),
    ) {
        (true, false) => {
            let kind = validate_alignment(&value)?;
            Ok(ImportedResult::Alignment(AlignmentResultDocument {
                kind,
                value,
            }))
        }
        (false, true) => {
            let kind = validate_cli(&value)?;
            Ok(ImportedResult::Cli(CliResultDocument { kind, value }))
        }
        _ => Err(invalid(
            "expected exactly one discriminator: schema or schema_version",
        )),
    }
}

fn validate_alignment(value: &Value) -> Result<AlignmentResultKind, BioIoError> {
    if string(field(value, "schema")?)? != "align-output-v1" {
        return Err(invalid("unsupported native schema version"));
    }
    nonempty(field(value, "algorithm")?)?;
    let alphabet = string(field(value, "alphabet")?)?;
    alphabet_symbols(alphabet)?;
    match string(field(value, "kind")?)? {
        "pairwise_alignment" => {
            validate_pairwise(value, alphabet, true)?;
            Ok(AlignmentResultKind::PairwiseAlignment)
        }
        "pairwise_value" => {
            validate_pairwise(value, alphabet, false)?;
            Ok(AlignmentResultKind::PairwiseValue)
        }
        "multiple_alignment" => {
            validate_multiple(value, alphabet)?;
            Ok(AlignmentResultKind::MultipleAlignment)
        }
        "myers_distance" => {
            validate_myers(value, alphabet, false)?;
            Ok(AlignmentResultKind::MyersDistance)
        }
        "myers_search" => {
            validate_myers(value, alphabet, true)?;
            Ok(AlignmentResultKind::MyersSearch)
        }
        "poasta_graph_query" => {
            validate_poasta(value, alphabet)?;
            Ok(AlignmentResultKind::PoastaGraphQuery)
        }
        kind => Err(unsupported(format!(
            "unsupported native result kind {kind:?}"
        ))),
    }
}

fn validate_pairwise(value: &Value, alphabet: &str, traceback: bool) -> Result<(), BioIoError> {
    let mut keys = vec![
        "schema",
        "kind",
        "algorithm",
        "mode",
        "alphabet",
        "scoring",
        "objective",
        "inputs",
        "details",
    ];
    if traceback {
        keys.extend([
            "ranges",
            "native_alignment_score",
            "statistics",
            "cigar",
            "alignment",
        ]);
    } else {
        keys.push("has_traceback");
    }
    exact_keys(value, &keys)?;
    nonempty(field(value, "mode")?)?;
    nonempty(field(value, "scoring")?)?;
    let objective = objective(field(value, "objective")?, false)?;
    report_fields(field(value, "details")?)?;
    let inputs = field(value, "inputs")?;
    exact_keys(inputs, &["first", "second"])?;
    let first = input_sequence(field(inputs, "first")?, alphabet)?;
    let second = input_sequence(field(inputs, "second")?, alphabet)?;
    if !traceback {
        if field(value, "has_traceback")?.as_bool() != Some(false) {
            return Err(invalid(
                "value-only results must declare has_traceback=false",
            ));
        }
        return Ok(());
    }
    let score = integer(field(value, "native_alignment_score")?)?;
    if objective.0 == "score" && score != objective.1 {
        return Err(invalid(
            "objective score differs from native_alignment_score",
        ));
    }
    let ranges = field(value, "ranges")?;
    exact_keys(ranges, &["first", "second"])?;
    let first_range = range(field(ranges, "first")?, first.len())?;
    let second_range = range(field(ranges, "second")?, second.len())?;
    let aligned = field(value, "alignment")?;
    exact_keys(aligned, &["first", "second"])?;
    let aligned_first = string(field(aligned, "first")?)?;
    let aligned_second = string(field(aligned, "second")?)?;
    validate_symbols(aligned_first, alphabet, true)?;
    validate_symbols(aligned_second, alphabet, true)?;
    if aligned_first.len() != aligned_second.len() {
        return Err(invalid("pairwise rows have different widths"));
    }
    if aligned_first
        .bytes()
        .filter(|byte| *byte != b'-')
        .ne(first[first_range].bytes())
        || aligned_second
            .bytes()
            .filter(|byte| *byte != b'-')
            .ne(second[second_range].bytes())
    {
        return Err(invalid(
            "ungapped rows do not match the declared input ranges",
        ));
    }
    validate_pairwise_columns(value, aligned_first, aligned_second)
}

#[allow(clippy::cast_precision_loss)]
fn validate_pairwise_columns(value: &Value, first: &str, second: &str) -> Result<(), BioIoError> {
    let mut counts = [0_usize; 7];
    counts[0] = first.len();
    let mut operations = Vec::with_capacity(first.len());
    let mut previous = None;
    for (left, right) in first.bytes().zip(second.bytes()) {
        let (operation, slot) = match (left, right) {
            (b'-', b'-') => return Err(invalid("pairwise alignment contains an all-gap column")),
            (b'-', _) => ('I', 3),
            (_, b'-') => ('D', 5),
            _ if left == right => ('=', 1),
            _ => ('X', 2),
        };
        counts[slot] += 1;
        if matches!(operation, 'I' | 'D') && previous != Some(operation) {
            counts[slot + 1] += 1;
        }
        previous = Some(operation);
        operations.push(operation);
    }
    if string(field(value, "cigar")?)? != run_length(&operations) {
        return Err(invalid("extended CIGAR disagrees with aligned rows"));
    }
    let names = [
        "columns",
        "matches",
        "substitutions",
        "insertion_residues",
        "insertion_runs",
        "deletion_residues",
        "deletion_runs",
        "identity",
    ];
    let stats = field(value, "statistics")?;
    exact_keys(stats, &names)?;
    for (name, expected) in names.iter().zip(counts) {
        if natural(field(stats, name)?)? != expected {
            return Err(invalid(format!(
                "pairwise statistic {name} disagrees with aligned rows"
            )));
        }
    }
    let identity = field(stats, "identity")?;
    if first.is_empty() {
        if !identity.is_null() {
            return Err(invalid("empty alignment identity must be null"));
        }
    } else {
        let actual = identity
            .as_f64()
            .ok_or_else(|| invalid("identity must be a number"))?;
        let expected = counts[1] as f64 / first.len() as f64;
        if (actual - expected).abs() > 1e-12 || !(0.0..=1.0).contains(&actual) {
            return Err(invalid("identity differs from matches divided by columns"));
        }
    }
    Ok(())
}

fn validate_multiple(value: &Value, alphabet: &str) -> Result<(), BioIoError> {
    exact_keys(
        value,
        &[
            "schema",
            "kind",
            "algorithm",
            "profile",
            "mode",
            "alphabet",
            "statistics",
            "alignment",
            "summary",
            "details",
        ],
    )?;
    nonempty(field(value, "profile")?)?;
    if !field(value, "mode")?.is_null() {
        nonempty(field(value, "mode")?)?;
    }
    report_fields(field(value, "summary")?)?;
    report_fields(field(value, "details")?)?;
    let stats = field(value, "statistics")?;
    exact_keys(
        stats,
        &[
            "sequences",
            "columns",
            "conserved_columns",
            "gap_characters",
        ],
    )?;
    let rows = array(field(value, "alignment")?)?;
    let width = natural(field(stats, "columns")?)?;
    if natural(field(stats, "sequences")?)? != rows.len() {
        return Err(invalid(
            "multiple-alignment row count disagrees with statistics",
        ));
    }
    let mut names = BTreeSet::new();
    let mut sequences = Vec::new();
    let mut gaps = 0;
    for row in rows {
        exact_keys(row, &["identifier", "sequence"])?;
        let name = nonempty(field(row, "identifier")?)?;
        if name.bytes().any(|byte| byte.is_ascii_control()) || !names.insert(name) {
            return Err(invalid("invalid or duplicate alignment identifier"));
        }
        let sequence = string(field(row, "sequence")?)?;
        validate_symbols(sequence, alphabet, true)?;
        if sequence.len() != width {
            return Err(invalid("multiple-alignment rows are not rectangular"));
        }
        gaps += sequence.bytes().filter(|byte| *byte == b'-').count();
        sequences.push(sequence.as_bytes());
    }
    if rows.is_empty() {
        return Err(invalid(
            "native multiple alignments require at least one row",
        ));
    }
    let mut conserved = 0;
    for column in 0..width {
        if sequences.iter().all(|row| row[column] == b'-') {
            return Err(invalid(
                "native multiple alignment contains an all-gap column",
            ));
        }
        let first = sequences[0][column];
        conserved += usize::from(first != b'-' && sequences.iter().all(|row| row[column] == first));
    }
    if natural(field(stats, "gap_characters")?)? != gaps
        || natural(field(stats, "conserved_columns")?)? != conserved
    {
        return Err(invalid(
            "multiple-alignment statistics disagree with the rows",
        ));
    }
    Ok(())
}

fn validate_myers(value: &Value, alphabet: &str, search: bool) -> Result<(), BioIoError> {
    let mut keys = vec![
        "schema",
        "kind",
        "algorithm",
        "mode",
        "alphabet",
        "pattern",
        "text",
    ];
    if search {
        keys.extend(["maximum_distance", "matches"]);
    } else {
        keys.push("distance");
    }
    exact_keys(value, &keys)?;
    if string(field(value, "algorithm")?)? != "Myers bit-vector"
        || string(field(value, "mode")?)?
            != if search {
                "approximate suffix search"
            } else {
                "global unit edit distance"
            }
    {
        return Err(invalid(
            "Myers result has inconsistent algorithm/mode metadata",
        ));
    }
    let pattern = input_sequence(field(value, "pattern")?, alphabet)?;
    let text = input_sequence(field(value, "text")?, alphabet)?;
    if search {
        let maximum = natural(field(value, "maximum_distance")?)?;
        let mut previous = None;
        for (index, matched) in array(field(value, "matches")?)?.iter().enumerate() {
            exact_keys(matched, &["rank", "text_end", "distance"])?;
            let end = natural(field(matched, "text_end")?)?;
            let distance = natural(field(matched, "distance")?)?;
            if natural(field(matched, "rank")?)? != index + 1
                || end > text.len()
                || previous.is_some_and(|previous| previous >= end)
                || distance > maximum
                || distance > pattern.len()
            {
                return Err(invalid(
                    "Myers matches violate ranks, endpoint ordering/bounds or distance bounds",
                ));
            }
            previous = Some(end);
        }
    } else {
        let distance = natural(field(value, "distance")?)?;
        if distance < pattern.len().abs_diff(text.len())
            || distance > pattern.len().max(text.len())
            || (distance == 0) != (pattern == text)
        {
            return Err(invalid(
                "global edit distance violates sequence-length or equality bounds",
            ));
        }
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn validate_poasta(value: &Value, alphabet: &str) -> Result<(), BioIoError> {
    exact_keys(
        value,
        &[
            "schema",
            "kind",
            "algorithm",
            "profile",
            "mode",
            "alphabet",
            "costs",
            "objective",
            "backend",
            "graph",
            "query",
            "search",
            "trace",
            "operations",
            "aligned_pairs",
        ],
    )?;
    if string(field(value, "algorithm")?)? != "POASTA"
        || string(field(value, "mode")?)? != "global affine sequence-to-graph"
    {
        return Err(invalid("inconsistent POASTA algorithm/mode metadata"));
    }
    nonempty(field(value, "profile")?)?;
    if !matches!(
        string(field(value, "backend")?)?,
        "sparse" | "dense" | "auto"
    ) {
        return Err(invalid("unknown POASTA backend"));
    }
    let costs = field(value, "costs")?;
    exact_keys(costs, &["mismatch", "gap_opening", "gap_extension"])?;
    align::advanced::poa::PoastaAffineCosts::new(
        integer(field(costs, "mismatch")?)?,
        integer(field(costs, "gap_opening")?)?,
        integer(field(costs, "gap_extension")?)?,
    )
    .map_err(|error| invalid(format!("invalid POASTA costs: {error}")))?;
    let metric = objective(field(value, "objective")?, true)?;
    if metric.0 != "score"
        || metric.1 > 0
        || string(field(field(value, "objective")?, "meaning")?)? != "negated affine path cost"
    {
        return Err(invalid("invalid POASTA objective"));
    }
    let graph = field(value, "graph")?;
    exact_keys(graph, &["nodes", "edges", "classes", "paths"])?;
    let nodes = natural(field(graph, "nodes")?)?;
    let classes = natural(field(graph, "classes")?)?;
    let paths = natural(field(graph, "paths")?)?;
    let edges = natural(field(graph, "edges")?)?;
    if classes > nodes
        || (nodes == 0) != (classes == 0)
        || (nodes > 0 && paths == 0)
        || (nodes == 0 && edges != 0)
    {
        return Err(invalid("inconsistent graph dimensions"));
    }
    let query = input_sequence(field(value, "query")?, alphabet)?;
    let search = field(value, "search")?;
    let counters = [
        "possible_states",
        "discovered_states",
        "expanded_states",
        "astar_expansions",
        "dfs_expansions",
        "greedy_matches",
        "bubble_prunes",
        "peak_queue_entries",
        "modeled_peak_bytes",
    ];
    let mut search_keys = counters.to_vec();
    search_keys.extend(["exact_match_dfs", "superbubble_pruning"]);
    exact_keys(search, &search_keys)?;
    for name in counters {
        natural(field(search, name)?)?;
    }
    for name in ["exact_match_dfs", "superbubble_pruning"] {
        boolean(field(search, name)?)?;
    }
    if natural(field(search, "discovered_states")?)? > natural(field(search, "possible_states")?)? {
        return Err(invalid("discovered states exceed possible states"));
    }
    let mut expected_query = 0;
    let mut seen_nodes = BTreeSet::new();
    let mut aligned_pairs = Vec::new();
    let mut codes = Vec::new();
    for operation in array(field(value, "operations")?)? {
        let name = string(field(operation, "operation")?)?;
        let (uses_node, uses_query, code) = match name {
            "aligned" => (true, true, 'A'),
            "graph_only" => (true, false, 'G'),
            "query_only" => (false, true, 'Q'),
            _ => return Err(invalid("unsupported graph operation")),
        };
        let mut keys = vec!["operation"];
        if uses_node {
            keys.push("node");
        }
        if uses_query {
            keys.push("query_index");
        }
        exact_keys(operation, &keys)?;
        let node = if uses_node {
            let node = natural(field(operation, "node")?)?;
            if node >= nodes || !seen_nodes.insert(node) {
                return Err(invalid("graph operation has an absent or repeated node"));
            }
            Some(node)
        } else {
            None
        };
        if uses_query {
            let index = natural(field(operation, "query_index")?)?;
            if index != expected_query || index >= query.len() {
                return Err(invalid(
                    "graph operations must consume query indices exactly once in order",
                ));
            }
            expected_query += 1;
            if let Some(node) = node {
                aligned_pairs.push((node, index));
            }
        }
        codes.push(code);
    }
    if expected_query != query.len() || string(field(value, "trace")?)? != run_length(&codes) {
        return Err(invalid("graph trace or query consumption is inconsistent"));
    }
    let pairs = array(field(value, "aligned_pairs")?)?;
    if pairs.len() != aligned_pairs.len() {
        return Err(invalid("aligned pair count differs from operations"));
    }
    for (pair, expected) in pairs.iter().zip(aligned_pairs) {
        exact_keys(pair, &["node", "query_index"])?;
        if (
            natural(field(pair, "node")?)?,
            natural(field(pair, "query_index")?)?,
        ) != expected
        {
            return Err(invalid("aligned pairs disagree with operations"));
        }
    }
    Ok(())
}

fn validate_cli(value: &Value) -> Result<CliResultKind, BioIoError> {
    if natural(field(value, "schema_version")?)? != 1 {
        return Err(invalid("unsupported CLI schema version"));
    }
    if !boolean(field(value, "ok")?)? {
        exact_keys(value, &["schema_version", "ok", "error"])?;
        let error = field(value, "error")?;
        exact_keys(error, &["code", "message", "details"])?;
        let code = nonempty(field(error, "code")?)?;
        if !code
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_')
        {
            return Err(invalid(
                "error code must use lowercase ASCII, digits or underscores",
            ));
        }
        nonempty(field(error, "message")?)?;
        object(field(error, "details")?)?;
        return Ok(CliResultKind::Error);
    }
    exact_keys(value, &["schema_version", "ok", "result"])?;
    let result = field(value, "result")?;
    match string(field(result, "operation")?)? {
        "pairwise" => {
            validate_cli_pairwise(result)?;
            Ok(CliResultKind::Pairwise)
        }
        "translate" => {
            validate_cli_translation(result)?;
            Ok(CliResultKind::Translation)
        }
        "complement" | "reverse-complement" | "transcribe" | "back-transcribe" => {
            validate_cli_transform(result)?;
            Ok(CliResultKind::SequenceTransform)
        }
        operation => Err(unsupported(format!(
            "CLI response operation {operation:?} is not supported by this result reader"
        ))),
    }
}

fn validate_cli_pairwise(result: &Value) -> Result<(), BioIoError> {
    exact_keys(
        result,
        &[
            "operation",
            "molecule",
            "target",
            "query",
            "algorithm",
            "scoring",
            "alignment",
        ],
    )?;
    let alphabet = cli_alphabet(string(field(result, "molecule")?)?)?;
    let target = string(field(result, "target")?)?;
    let query = string(field(result, "query")?)?;
    validate_symbols(target, alphabet, false)?;
    validate_symbols(query, alphabet, false)?;
    let nested = field(result, "alignment")?;
    if validate_alignment(nested)? != AlignmentResultKind::PairwiseAlignment
        || string(field(nested, "alphabet")?)? != alphabet
    {
        return Err(invalid(
            "CLI pairwise result has a different nested kind/alphabet",
        ));
    }
    let inputs = field(nested, "inputs")?;
    if input_sequence(field(inputs, "first")?, alphabet)? != target
        || input_sequence(field(inputs, "second")?, alphabet)? != query
    {
        return Err(invalid(
            "CLI target/query differ from the nested alignment inputs",
        ));
    }
    let scoring = field(result, "scoring")?;
    let algorithm = string(field(result, "algorithm")?)?;
    let expected_algorithm = match algorithm {
        "gotoh" => "Gotoh",
        "needleman-wunsch" => "Needleman-Wunsch",
        "smith-waterman" => "Smith-Waterman",
        _ => return Err(unsupported("unsupported CLI pairwise algorithm")),
    };
    if string(field(nested, "algorithm")?)? != expected_algorithm {
        return Err(invalid("CLI algorithm differs from the nested alignment"));
    }
    let (expected_mode, expected_scoring) = match algorithm {
        "gotoh" => (
            "global affine",
            "match=2 mismatch=-1 gap_open=-5 gap_extend=-1",
        ),
        "needleman-wunsch" => ("global linear", "match=2 mismatch=-1 gap=-2"),
        _ => ("local linear", "match=2 mismatch=-1 gap=-2"),
    };
    if string(field(nested, "mode")?)? != expected_mode
        || string(field(nested, "scoring")?)? != expected_scoring
        || string(field(field(nested, "objective")?, "name")?)? != "score"
    {
        return Err(invalid(
            "nested CLI mode/scoring/objective differs from its contract",
        ));
    }
    if algorithm != "smith-waterman" {
        let ranges = field(nested, "ranges")?;
        if range(field(ranges, "first")?, target.len())? != (0..target.len())
            || range(field(ranges, "second")?, query.len())? != (0..query.len())
        {
            return Err(invalid(
                "global CLI alignment must span both complete inputs",
            ));
        }
    }
    let keys = if algorithm == "gotoh" {
        vec![
            "match",
            "mismatch",
            "gap_open",
            "gap_extend",
            "gap_convention",
        ]
    } else {
        vec!["match", "mismatch", "gap"]
    };
    exact_keys(scoring, &keys)?;
    if integer(field(scoring, "match")?)? != 2 || integer(field(scoring, "mismatch")?)? != -1 {
        return Err(invalid("CLI scoring differs from the version 1 contract"));
    }
    if algorithm == "gotoh" {
        if integer(field(scoring, "gap_open")?)? != -5
            || integer(field(scoring, "gap_extend")?)? != -1
            || string(field(scoring, "gap_convention")?)? != "open + (length - 1) * extension"
        {
            return Err(invalid(
                "CLI affine gaps differ from the version 1 contract",
            ));
        }
    } else if integer(field(scoring, "gap")?)? != -2 {
        return Err(invalid(
            "CLI linear gap differs from the version 1 contract",
        ));
    }
    Ok(())
}

fn validate_cli_transform(result: &Value) -> Result<(), BioIoError> {
    exact_keys(
        result,
        &[
            "operation",
            "molecule",
            "input",
            "input_length",
            "output_molecule",
            "sequence",
            "length",
        ],
    )?;
    let molecule = string(field(result, "molecule")?)?;
    let input = string(field(result, "input")?)?;
    validate_symbols(input, cli_alphabet(molecule)?, false)?;
    let operation = string(field(result, "operation")?)?;
    let output_molecule = match operation {
        "transcribe" => "rna",
        "back-transcribe" => "dna",
        _ => molecule,
    };
    if string(field(result, "output_molecule")?)? != output_molecule {
        return Err(invalid("incorrect output molecule"));
    }
    let expected = match (molecule, operation) {
        ("dna", "complement") => crate::dna(input).and_then(|sequence| sequence.complement()),
        ("dna", "reverse-complement") => {
            crate::dna(input).and_then(|sequence| sequence.reverse_complement())
        }
        ("rna", "back-transcribe") => {
            crate::rna(input).and_then(|sequence| sequence.back_transcribe())
        }
        ("rna", "complement" | "reverse-complement") | ("dna", "transcribe") => {
            let expected = match operation {
                "complement" => crate::rna(input).and_then(|sequence| sequence.complement()),
                "reverse-complement" => {
                    crate::rna(input).and_then(|sequence| sequence.reverse_complement())
                }
                _ => crate::dna(input).and_then(|sequence| sequence.transcribe()),
            }
            .map_err(|error| invalid(error.to_string()))?
            .to_string();
            return check_transform_output(result, input, &expected);
        }
        _ => return Err(invalid("operation is incompatible with input molecule")),
    }
    .map_err(|error| invalid(error.to_string()))?
    .to_string();
    check_transform_output(result, input, &expected)
}

fn check_transform_output(result: &Value, input: &str, expected: &str) -> Result<(), BioIoError> {
    if natural(field(result, "input_length")?)? != input.len()
        || natural(field(result, "length")?)? != expected.len()
        || string(field(result, "sequence")?)? != expected
    {
        return Err(invalid(
            "sequence transform output or lengths disagree with the input",
        ));
    }
    Ok(())
}

fn validate_cli_translation(result: &Value) -> Result<(), BioIoError> {
    exact_keys(
        result,
        &[
            "operation",
            "molecule",
            "input",
            "input_length",
            "genetic_code",
            "reading_frame",
            "mode",
            "partial",
            "peptide",
            "length",
            "ignored_trailing_bases",
        ],
    )?;
    if natural(field(result, "genetic_code")?)? != 1
        || natural(field(result, "reading_frame")?)? != 0
    {
        return Err(unsupported(
            "only genetic code 1 and frame zero are supported",
        ));
    }
    let mode = match string(field(result, "mode")?)? {
        "full" => TranslationMode::Full,
        "to-stop" => TranslationMode::ToStop,
        "cds" => TranslationMode::Cds,
        _ => return Err(invalid("unknown translation mode")),
    };
    let partial_codon = match string(field(result, "partial")?)? {
        "reject" => PartialCodonPolicy::Reject,
        "ignore" => PartialCodonPolicy::Ignore,
        _ => return Err(invalid("unknown partial-codon policy")),
    };
    let input = string(field(result, "input")?)?;
    let molecule = string(field(result, "molecule")?)?;
    validate_symbols(input, cli_alphabet(molecule)?, false)?;
    let options = TranslationOptions {
        mode,
        partial_codons: partial_codon,
    };
    let translated = match molecule {
        "dna" => crate::dna(input)
            .map_err(|error| invalid(error.to_string()))?
            .translate_standard(options),
        "rna" => crate::rna(input)
            .map_err(|error| invalid(error.to_string()))?
            .translate_standard(options),
        _ => return Err(invalid("translation requires DNA or RNA")),
    }
    .map_err(|error| invalid(error.to_string()))?;
    if natural(field(result, "input_length")?)? != input.len()
        || string(field(result, "peptide")?)? != translated.as_str()
        || natural(field(result, "length")?)? != translated.len()
        || natural(field(result, "ignored_trailing_bases")?)? != translated.ignored_trailing_bases()
    {
        return Err(invalid(
            "translation output or lengths disagree with its explicit policies",
        ));
    }
    Ok(())
}

fn objective(value: &Value, meaning: bool) -> Result<(&str, i64), BioIoError> {
    let mut keys = vec!["name", "direction", "value"];
    if meaning {
        keys.push("meaning");
    }
    exact_keys(value, &keys)?;
    let name = string(field(value, "name")?)?;
    let direction = string(field(value, "direction")?)?;
    let number = integer(field(value, "value")?)?;
    if !matches!(
        (name, direction),
        ("score", "maximize") | ("distance" | "penalty", "minimize")
    ) || (name != "score" && number < 0)
    {
        return Err(invalid("invalid objective name, direction or sign"));
    }
    Ok((name, number))
}

fn report_fields(value: &Value) -> Result<(), BioIoError> {
    for item in array(value)? {
        exact_keys(item, &["label", "value"])?;
        string(field(item, "label")?)?;
        string(field(item, "value")?)?;
    }
    Ok(())
}

fn input_sequence<'a>(value: &'a Value, alphabet: &str) -> Result<&'a str, BioIoError> {
    exact_keys(value, &["length", "sequence"])?;
    let sequence = string(field(value, "sequence")?)?;
    validate_symbols(sequence, alphabet, false)?;
    if natural(field(value, "length")?)? != sequence.len() {
        return Err(invalid(
            "declared sequence length differs from actual length",
        ));
    }
    Ok(sequence)
}

fn range(value: &Value, length: usize) -> Result<std::ops::Range<usize>, BioIoError> {
    exact_keys(value, &["start", "end"])?;
    let start = natural(field(value, "start")?)?;
    let end = natural(field(value, "end")?)?;
    if start > end || end > length {
        return Err(invalid("invalid zero-based half-open sequence range"));
    }
    Ok(start..end)
}

fn validate_symbols(sequence: &str, alphabet: &str, gaps: bool) -> Result<(), BioIoError> {
    let symbols = alphabet_symbols(alphabet)?;
    if !sequence
        .bytes()
        .all(|byte| symbols.contains(&byte) || (gaps && byte == b'-'))
    {
        return Err(invalid(format!(
            "sequence contains noncanonical {alphabet} symbols"
        )));
    }
    Ok(())
}

fn alphabet_symbols(alphabet: &str) -> Result<&'static [u8], BioIoError> {
    match alphabet {
        "DNA" => Ok(Dna::SYMBOLS),
        "RNA" => Ok(Rna::SYMBOLS),
        "protein" => Ok(Protein::SYMBOLS),
        _ => Err(unsupported(format!("unsupported alphabet {alphabet:?}"))),
    }
}

fn cli_alphabet(molecule: &str) -> Result<&'static str, BioIoError> {
    match molecule {
        "dna" => Ok("DNA"),
        "rna" => Ok("RNA"),
        "protein" => Ok("protein"),
        _ => Err(unsupported("unsupported molecule")),
    }
}

fn run_length(operations: &[char]) -> String {
    let mut result = String::new();
    let mut offset = 0;
    while let Some(operation) = operations.get(offset) {
        let count = operations[offset..]
            .iter()
            .take_while(|next| *next == operation)
            .count();
        result.push_str(&count.to_string());
        result.push(*operation);
        offset += count;
    }
    result
}

fn object(value: &Value) -> Result<&Map<String, Value>, BioIoError> {
    value
        .as_object()
        .ok_or_else(|| invalid("expected a JSON object"))
}
fn field<'a>(value: &'a Value, name: &str) -> Result<&'a Value, BioIoError> {
    object(value)?
        .get(name)
        .ok_or_else(|| invalid(format!("missing field {name:?}")))
}
fn string(value: &Value) -> Result<&str, BioIoError> {
    value.as_str().ok_or_else(|| invalid("expected a string"))
}
fn nonempty(value: &Value) -> Result<&str, BioIoError> {
    let text = string(value)?;
    if text.trim().is_empty() {
        Err(invalid("expected a nonempty string"))
    } else {
        Ok(text)
    }
}
fn natural(value: &Value) -> Result<usize, BioIoError> {
    value
        .as_u64()
        .and_then(|number| usize::try_from(number).ok())
        .ok_or_else(|| invalid("expected a platform-sized nonnegative integer"))
}
fn integer(value: &Value) -> Result<i64, BioIoError> {
    value
        .as_i64()
        .ok_or_else(|| invalid("expected a signed 64-bit integer"))
}
fn boolean(value: &Value) -> Result<bool, BioIoError> {
    value.as_bool().ok_or_else(|| invalid("expected a boolean"))
}
fn array(value: &Value) -> Result<&[Value], BioIoError> {
    value
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| invalid("expected an array"))
}

fn exact_keys(value: &Value, keys: &[&str]) -> Result<(), BioIoError> {
    let object = object(value)?;
    if object.len() != keys.len() || object.keys().any(|key| !keys.contains(&key.as_str())) {
        return Err(invalid(format!(
            "expected exactly these fields: {}",
            keys.join(", ")
        )));
    }
    Ok(())
}

// The JSON parser has already verified grammar. This second pass rejects duplicate
// keys (including escaped aliases) that Value would otherwise silently overwrite.
fn reject_duplicate_keys(input: &str) -> Result<(), BioIoError> {
    let bytes = input.as_bytes();
    let mut stack: Vec<Option<BTreeSet<String>>> = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        match bytes[index] {
            b'{' => stack.push(Some(BTreeSet::new())),
            b'[' => stack.push(None),
            b'}' | b']' => {
                stack.pop();
            }
            b'"' => {
                let start = index;
                index += 1;
                while index < bytes.len() && bytes[index] != b'"' {
                    if bytes[index] == b'\\' {
                        index += 1;
                    }
                    index += 1;
                }
                let end = index + 1;
                let mut next = end;
                while next < bytes.len() && bytes[next].is_ascii_whitespace() {
                    next += 1;
                }
                if bytes.get(next) == Some(&b':') {
                    let key: String = serde_json::from_str(&input[start..end])
                        .map_err(|error| invalid(error.to_string()))?;
                    if let Some(Some(keys)) = stack.last_mut() {
                        if !keys.insert(key.clone()) {
                            return Err(invalid(format!("duplicate JSON key {key:?}")));
                        }
                    }
                }
            }
            _ => {}
        }
        index += 1;
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> BioIoError {
    BioIoError::new(IoErrorKind::InvalidSchema, "JSON result", message)
}
fn unsupported(message: impl Into<String>) -> BioIoError {
    BioIoError::new(IoErrorKind::Unsupported, "JSON result", message)
}
