//! Direct file and stream output for ergonomic facade runs.
//!
//! Human text, structured JSON, alignment interchange formats, and POASTA graph
//! formats are kept separate. Every writer consumes the native typed result; no
//! writer parses a human report to recover scientific data.

use core::fmt;
use core::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::{Alphabet, Sequence};
use align_formats::{FormatError, NexusAlphabet, PhylipOptions};
use align_poa::{PoaOperation, PoastaBackend, PoastaPoaResult};

use crate::Result;
use crate::report::{
    MultipleRun, MyersDistanceRun, MyersSearchRun, PairwiseMetric, PairwiseRun, PairwiseValueRun,
    PoastaQueryRun, ReportField,
};
use crate::result::{MultipleAlignmentView, PairwiseAlignmentExt, PairwiseAlignmentView};

/// Stable schema identifier emitted by every facade JSON writer in this version.
pub const JSON_SCHEMA_V1: &str = "align-output-v1";

type FormatResult = core::result::Result<(), FormatError>;

impl<A, R> PairwiseRun<'_, A, R>
where
    A: Alphabet,
    R: PairwiseAlignmentView,
{
    /// Writes the normal human-readable report as UTF-8 text.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the normal human-readable report as UTF-8 text.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes a self-contained versioned JSON alignment result.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when the traceback cannot be rendered or output
    /// fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_pairwise_json(writer, self))
    }

    /// Streams a self-contained versioned JSON alignment result.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when the traceback cannot be rendered or output
    /// fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_pairwise_json(writer, self)?;
        Ok(())
    }

    /// Writes the traceback as two aligned FASTA records.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] for invalid identifiers, an inconsistent
    /// traceback, or failed output.
    pub fn write_fasta(
        &self,
        path: impl AsRef<Path>,
        first_identifier: &str,
        second_identifier: &str,
    ) -> Result<()> {
        align_formats::write_pairwise_fasta(
            path,
            first_identifier,
            self.first(),
            second_identifier,
            self.second(),
            self.alignment(),
        )?;
        Ok(())
    }
}

impl<A: Alphabet, R> PairwiseValueRun<'_, A, R> {
    /// Writes the human-readable score, distance, or penalty summary.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the human-readable score, distance, or penalty summary.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes a versioned JSON value result without fabricating a traceback.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_pairwise_value_json(writer, self))
    }

    /// Streams a versioned JSON value result without fabricating a traceback.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when output fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_pairwise_value_json(writer, self)?;
        Ok(())
    }
}

impl<A, R> MultipleRun<A, R>
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    /// Writes the normal human-readable multiple-alignment report.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the normal human-readable multiple-alignment report.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes a self-contained versioned JSON multiple alignment.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_multiple_json(writer, self))
    }

    /// Streams a self-contained versioned JSON multiple alignment.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when output fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_multiple_json(writer, self)?;
        Ok(())
    }

    /// Writes the rectangular alignment as aligned FASTA.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_fasta(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_alignment_fasta(path, self.alignment())?;
        Ok(())
    }

    /// Writes the rectangular alignment in CLUSTAL format.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_clustal(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_clustal(path, self.alignment())?;
        Ok(())
    }

    /// Writes relaxed sequential PHYLIP without truncating taxon names.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] for a nonportable taxon identifier or failed
    /// output.
    pub fn write_phylip(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_phylip(path, self.alignment())?;
        Ok(())
    }

    /// Writes PHYLIP with explicit strict/relaxed and layout options.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] for an invalid taxon identifier or failed output.
    pub fn write_phylip_with(&self, path: impl AsRef<Path>, options: PhylipOptions) -> Result<()> {
        align_formats::write_phylip_with_options(path, self.alignment(), options)?;
        Ok(())
    }

    /// Writes the rectangular alignment as unwrapped Stockholm 1.0.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] for a whitespace-containing sequence name or
    /// failed output.
    pub fn write_stockholm(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_stockholm(path, self.alignment())?;
        Ok(())
    }
}

impl<A, R> MultipleRun<A, R>
where
    A: NexusAlphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    /// Writes a NEXUS `DATA` block with the exact typed molecule datatype.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] for an invalid taxon label or failed output.
    pub fn write_nexus(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_nexus(path, self.alignment())?;
        Ok(())
    }
}

impl<A: Alphabet> MultipleRun<A, PoastaPoaResult<A>> {
    /// Writes the authoritative progressive POASTA graph as GFA 1.0.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_gfa(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_poa_gfa(path, self.result().graph())?;
        Ok(())
    }

    /// Writes the authoritative progressive POASTA graph as Graphviz DOT.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_dot(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_poa_dot(path, self.result().graph())?;
        Ok(())
    }
}

impl<A: Alphabet> MyersDistanceRun<'_, A> {
    /// Writes the human-readable global edit-distance summary.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the human-readable global edit-distance summary.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes the global edit distance as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_myers_distance_json(writer, self))
    }

    /// Streams the global edit distance as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when output fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_myers_distance_json(writer, self)?;
        Ok(())
    }
}

impl<A: Alphabet> MyersSearchRun<'_, A> {
    /// Writes the human-readable approximate-search table.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the human-readable approximate-search table.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes the ordered approximate-search matches as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_myers_search_json(writer, self))
    }

    /// Streams the ordered approximate-search matches as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when output fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_myers_search_json(writer, self)?;
        Ok(())
    }

    /// Writes the ordered approximate-search matches as tab-separated values.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_tsv(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_myers_tsv(path, self.matches())?;
        Ok(())
    }
}

impl<A: Alphabet> PoastaQueryRun<'_, A> {
    /// Writes the normal human-readable graph-query report.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or file output fails.
    pub fn write_text(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_text_report(path, &self.report())?;
        Ok(())
    }

    /// Streams the normal human-readable graph-query report.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when rendering or output fails.
    pub fn write_text_to(&self, writer: impl Write) -> Result<()> {
        align_formats::write_text_report_to(writer, &self.report())?;
        Ok(())
    }

    /// Writes the exact graph-query traceback and evidence as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_json(&self, path: impl AsRef<Path>) -> Result<()> {
        write_json_path(path, |writer| write_poasta_query_json(writer, self))
    }

    /// Streams the exact graph-query traceback and evidence as versioned JSON.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when output fails.
    pub fn write_json_to(&self, writer: impl Write) -> Result<()> {
        write_poasta_query_json(writer, self)?;
        Ok(())
    }

    /// Writes the authoritative input graph as GFA 1.0.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_gfa(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_poa_gfa(path, self.graph())?;
        Ok(())
    }

    /// Writes the authoritative input graph as Graphviz DOT.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error`] when file output fails.
    pub fn write_dot(&self, path: impl AsRef<Path>) -> Result<()> {
        align_formats::write_poa_dot(path, self.graph())?;
        Ok(())
    }
}

fn write_json_path(
    path: impl AsRef<Path>,
    render: impl FnOnce(&mut BufWriter<File>) -> FormatResult,
) -> Result<()> {
    let file = File::create(path).map_err(FormatError::from)?;
    let mut writer = BufWriter::new(file);
    render(&mut writer)?;
    writer.flush().map_err(FormatError::from)?;
    Ok(())
}

fn write_pairwise_json<A, R>(writer: impl Write, run: &PairwiseRun<'_, A, R>) -> FormatResult
where
    A: Alphabet,
    R: PairwiseAlignmentView,
{
    let alignment = run.alignment();
    let rendered = alignment.render(run.first(), run.second())?;
    let statistics = run.statistics();
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(",\n  \"kind\": \"pairwise_alignment\",\n  \"algorithm\": ")?;
    json.string(run.info().algorithm())?;
    json.raw(",\n  \"mode\": ")?;
    json.string(run.info().mode())?;
    json.raw(",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"scoring\": ")?;
    json.string(run.info().scoring())?;
    json.raw(",\n  \"objective\": ")?;
    write_metric(&mut json, run.info().metric())?;
    json.raw(",\n  \"inputs\": {\n    \"first\": {\"length\": ")?;
    json.display(run.first().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.first())?;
    json.raw("},\n    \"second\": {\"length\": ")?;
    json.display(run.second().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.second())?;
    json.raw("}\n  },\n  \"ranges\": {\n    \"first\": ")?;
    write_range(&mut json, alignment.first_range())?;
    json.raw(",\n    \"second\": ")?;
    write_range(&mut json, alignment.second_range())?;
    json.raw("\n  },\n  \"native_alignment_score\": ")?;
    json.display(alignment.score())?;
    json.raw(",\n  \"statistics\": {\"columns\": ")?;
    json.display(statistics.columns())?;
    json.raw(", \"matches\": ")?;
    json.display(statistics.matches())?;
    json.raw(", \"substitutions\": ")?;
    json.display(statistics.substitutions())?;
    json.raw(", \"insertion_residues\": ")?;
    json.display(statistics.insertion_residues())?;
    json.raw(", \"insertion_runs\": ")?;
    json.display(statistics.insertion_runs())?;
    json.raw(", \"deletion_residues\": ")?;
    json.display(statistics.deletion_residues())?;
    json.raw(", \"deletion_runs\": ")?;
    json.display(statistics.deletion_runs())?;
    json.raw(", \"identity\": ")?;
    if let Some(identity) = statistics.identity() {
        json.display(identity)?;
    } else {
        json.raw("null")?;
    }
    json.raw("},\n  \"cigar\": ")?;
    json.string(&run.cigar())?;
    json.raw(",\n  \"alignment\": {\"first\": ")?;
    json.string(rendered.first())?;
    json.raw(", \"second\": ")?;
    json.string(rendered.second())?;
    json.raw("},\n  \"details\": ")?;
    write_fields(&mut json, run.info().details())?;
    json.raw("\n}\n")
}

fn write_pairwise_value_json<A: Alphabet, R>(
    writer: impl Write,
    run: &PairwiseValueRun<'_, A, R>,
) -> FormatResult {
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(",\n  \"kind\": \"pairwise_value\",\n  \"algorithm\": ")?;
    json.string(run.info().algorithm())?;
    json.raw(",\n  \"mode\": ")?;
    json.string(run.info().mode())?;
    json.raw(",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"scoring\": ")?;
    json.string(run.info().scoring())?;
    json.raw(",\n  \"objective\": ")?;
    write_metric(&mut json, run.info().metric())?;
    json.raw(",\n  \"has_traceback\": false,\n  \"inputs\": {\n    \"first\": {\"length\": ")?;
    json.display(run.first().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.first())?;
    json.raw("},\n    \"second\": {\"length\": ")?;
    json.display(run.second().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.second())?;
    json.raw("}\n  },\n  \"details\": ")?;
    write_fields(&mut json, run.info().details())?;
    json.raw("\n}\n")
}

fn write_multiple_json<A, R>(writer: impl Write, run: &MultipleRun<A, R>) -> FormatResult
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    let statistics = run.statistics();
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(",\n  \"kind\": \"multiple_alignment\",\n  \"algorithm\": ")?;
    json.string(run.info().algorithm())?;
    json.raw(",\n  \"profile\": ")?;
    json.string(run.info().profile())?;
    json.raw(",\n  \"mode\": ")?;
    if let Some(mode) = run.info().mode() {
        json.string(mode)?;
    } else {
        json.raw("null")?;
    }
    json.raw(",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"statistics\": {\"sequences\": ")?;
    json.display(statistics.sequences())?;
    json.raw(", \"columns\": ")?;
    json.display(statistics.columns())?;
    json.raw(", \"conserved_columns\": ")?;
    json.display(statistics.conserved_columns())?;
    json.raw(", \"gap_characters\": ")?;
    json.display(statistics.gap_characters())?;
    json.raw("},\n  \"alignment\": [")?;
    for (index, row) in run.alignment().rows().iter().enumerate() {
        if index != 0 {
            json.raw(",")?;
        }
        json.raw("\n    {\"identifier\": ")?;
        json.string(row.identifier())?;
        json.raw(", \"sequence\": ")?;
        json.aligned_row(row.symbols())?;
        json.raw("}")?;
    }
    if !run.alignment().rows().is_empty() {
        json.raw("\n  ")?;
    }
    json.raw("],\n  \"summary\": ")?;
    write_fields(&mut json, run.info().summary())?;
    json.raw(",\n  \"details\": ")?;
    write_fields(&mut json, run.info().details())?;
    json.raw("\n}\n")
}

fn write_myers_distance_json<A: Alphabet>(
    writer: impl Write,
    run: &MyersDistanceRun<'_, A>,
) -> FormatResult {
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(",\n  \"kind\": \"myers_distance\",\n  \"algorithm\": \"Myers bit-vector\",\n  \"mode\": \"global unit edit distance\",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"pattern\": {\"length\": ")?;
    json.display(run.pattern().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.pattern())?;
    json.raw("},\n  \"text\": {\"length\": ")?;
    json.display(run.text().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.text())?;
    json.raw("},\n  \"distance\": ")?;
    json.display(run.distance())?;
    json.raw("\n}\n")
}

fn write_myers_search_json<A: Alphabet>(
    writer: impl Write,
    run: &MyersSearchRun<'_, A>,
) -> FormatResult {
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(",\n  \"kind\": \"myers_search\",\n  \"algorithm\": \"Myers bit-vector\",\n  \"mode\": \"approximate suffix search\",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"pattern\": {\"length\": ")?;
    json.display(run.pattern().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.pattern())?;
    json.raw("},\n  \"text\": {\"length\": ")?;
    json.display(run.text().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.text())?;
    json.raw("},\n  \"maximum_distance\": ")?;
    json.display(run.maximum_distance())?;
    json.raw(",\n  \"matches\": [")?;
    for (index, matched) in run.matches().iter().copied().enumerate() {
        if index != 0 {
            json.raw(",")?;
        }
        json.raw("\n    {\"rank\": ")?;
        json.display(index + 1)?;
        json.raw(", \"text_end\": ")?;
        json.display(matched.end())?;
        json.raw(", \"distance\": ")?;
        json.display(matched.distance())?;
        json.raw("}")?;
    }
    if !run.matches().is_empty() {
        json.raw("\n  ")?;
    }
    json.raw("]\n}\n")
}

fn write_poasta_query_json<A: Alphabet>(
    writer: impl Write,
    run: &PoastaQueryRun<'_, A>,
) -> FormatResult {
    let result = run.result();
    let statistics = result.statistics();
    let memory = result.memory_estimate();
    let options = result.options();
    let costs = run.costs();
    let mut json = JsonWriter::new(writer);
    json.raw("{\n  \"schema\": ")?;
    json.string(JSON_SCHEMA_V1)?;
    json.raw(
        ",\n  \"kind\": \"poasta_graph_query\",\n  \"algorithm\": \"POASTA\",\n  \"profile\": ",
    )?;
    json.string(run.profile())?;
    json.raw(",\n  \"mode\": \"global affine sequence-to-graph\",\n  \"alphabet\": ")?;
    json.string(A::NAME)?;
    json.raw(",\n  \"costs\": {\"mismatch\": ")?;
    json.display(costs.mismatch())?;
    json.raw(", \"gap_opening\": ")?;
    json.display(costs.gap_opening())?;
    json.raw(", \"gap_extension\": ")?;
    json.display(costs.gap_extension())?;
    json.raw(
        "},\n  \"objective\": {\"name\": \"score\", \"direction\": \"maximize\", \"value\": ",
    )?;
    json.display(result.score())?;
    json.raw(", \"meaning\": \"negated affine path cost\"},\n  \"backend\": ")?;
    json.string(poasta_backend(result.backend()))?;
    json.raw(",\n  \"graph\": {\"nodes\": ")?;
    json.display(run.graph().node_count())?;
    json.raw(", \"edges\": ")?;
    json.display(run.graph().edge_count())?;
    json.raw(", \"classes\": ")?;
    json.display(run.graph().class_count())?;
    json.raw(", \"paths\": ")?;
    json.display(run.graph().path_count())?;
    json.raw("},\n  \"query\": {\"length\": ")?;
    json.display(run.query().len())?;
    json.raw(", \"sequence\": ")?;
    json.sequence(run.query())?;
    json.raw("},\n  \"search\": {\"possible_states\": ")?;
    json.display(statistics.total_possible_states())?;
    json.raw(", \"discovered_states\": ")?;
    json.display(statistics.discovered_states())?;
    json.raw(", \"expanded_states\": ")?;
    json.display(statistics.expanded_states())?;
    json.raw(", \"astar_expansions\": ")?;
    json.display(statistics.astar_expanded_states())?;
    json.raw(", \"dfs_expansions\": ")?;
    json.display(statistics.dfs_expanded_states())?;
    json.raw(", \"greedy_matches\": ")?;
    json.display(statistics.greedy_exact_matches())?;
    json.raw(", \"bubble_prunes\": ")?;
    json.display(statistics.bubble_pruned_states())?;
    json.raw(", \"peak_queue_entries\": ")?;
    json.display(statistics.peak_queue_entries())?;
    json.raw(", \"modeled_peak_bytes\": ")?;
    json.display(memory.modeled_peak_bytes())?;
    json.raw(", \"exact_match_dfs\": ")?;
    json.boolean(options.exact_match_dfs())?;
    json.raw(", \"superbubble_pruning\": ")?;
    json.boolean(options.superbubble_pruning())?;
    json.raw("},\n  \"trace\": ")?;
    json.string(&poasta_trace(result.operations()))?;
    json.raw(",\n  \"operations\": [")?;
    for (index, operation) in result.operations().iter().copied().enumerate() {
        if index != 0 {
            json.raw(",")?;
        }
        json.raw("\n    ")?;
        write_poa_operation(&mut json, operation)?;
    }
    if !result.operations().is_empty() {
        json.raw("\n  ")?;
    }
    json.raw("],\n  \"aligned_pairs\": [")?;
    for (index, pair) in result.aligned_pairs().iter().copied().enumerate() {
        if index != 0 {
            json.raw(",")?;
        }
        json.raw("\n    {\"node\": ")?;
        json.display(pair.node().index())?;
        json.raw(", \"query_index\": ")?;
        json.display(pair.query_index())?;
        json.raw("}")?;
    }
    if !result.aligned_pairs().is_empty() {
        json.raw("\n  ")?;
    }
    json.raw("]\n}\n")
}

fn write_metric<W: Write>(writer: &mut JsonWriter<W>, metric: PairwiseMetric) -> FormatResult {
    match metric {
        PairwiseMetric::Score(value) => {
            writer.raw("{\"name\": \"score\", \"direction\": \"maximize\", \"value\": ")?;
            writer.display(value)?;
        }
        PairwiseMetric::Distance(value) => {
            writer.raw("{\"name\": \"distance\", \"direction\": \"minimize\", \"value\": ")?;
            writer.display(value)?;
        }
        PairwiseMetric::Penalty(value) => {
            writer.raw("{\"name\": \"penalty\", \"direction\": \"minimize\", \"value\": ")?;
            writer.display(value)?;
        }
    }
    writer.raw("}")
}

fn write_range<W: Write>(
    writer: &mut JsonWriter<W>,
    range: core::ops::Range<usize>,
) -> FormatResult {
    writer.raw("{\"start\": ")?;
    writer.display(range.start)?;
    writer.raw(", \"end\": ")?;
    writer.display(range.end)?;
    writer.raw("}")
}

fn write_fields<W: Write>(writer: &mut JsonWriter<W>, fields: &[ReportField]) -> FormatResult {
    writer.raw("[")?;
    for (index, field) in fields.iter().enumerate() {
        if index != 0 {
            writer.raw(", ")?;
        }
        writer.raw("{\"label\": ")?;
        writer.string(field.label())?;
        writer.raw(", \"value\": ")?;
        writer.string(field.value())?;
        writer.raw("}")?;
    }
    writer.raw("]")
}

fn write_poa_operation<W: Write>(
    writer: &mut JsonWriter<W>,
    operation: PoaOperation,
) -> FormatResult {
    match operation {
        PoaOperation::Aligned { node, query_index } => {
            writer.raw("{\"operation\": \"aligned\", \"node\": ")?;
            writer.display(node.index())?;
            writer.raw(", \"query_index\": ")?;
            writer.display(query_index)?;
        }
        PoaOperation::GraphOnly { node } => {
            writer.raw("{\"operation\": \"graph_only\", \"node\": ")?;
            writer.display(node.index())?;
        }
        PoaOperation::QueryOnly { query_index } => {
            writer.raw("{\"operation\": \"query_only\", \"query_index\": ")?;
            writer.display(query_index)?;
        }
    }
    writer.raw("}")
}

const fn poasta_backend(backend: PoastaBackend) -> &'static str {
    match backend {
        PoastaBackend::Sparse => "sparse",
        PoastaBackend::Dense => "dense",
        PoastaBackend::Auto => "auto",
    }
}

fn poasta_trace(operations: &[PoaOperation]) -> String {
    let Some(first) = operations.first().copied() else {
        return String::new();
    };
    let mut output = String::new();
    let mut selected = poa_operation_code(first);
    let mut count = 0_usize;
    for operation in operations.iter().copied() {
        let code = poa_operation_code(operation);
        if code == selected {
            count += 1;
        } else {
            write!(output, "{count}{selected}").expect("writing to a String cannot fail");
            selected = code;
            count = 1;
        }
    }
    write!(output, "{count}{selected}").expect("writing to a String cannot fail");
    output
}

const fn poa_operation_code(operation: PoaOperation) -> char {
    match operation {
        PoaOperation::Aligned { .. } => 'A',
        PoaOperation::GraphOnly { .. } => 'G',
        PoaOperation::QueryOnly { .. } => 'Q',
    }
}

struct JsonWriter<W> {
    writer: W,
}

impl<W: Write> JsonWriter<W> {
    const fn new(writer: W) -> Self {
        Self { writer }
    }

    fn raw(&mut self, value: &str) -> FormatResult {
        self.writer.write_all(value.as_bytes())?;
        Ok(())
    }

    fn display(&mut self, value: impl fmt::Display) -> FormatResult {
        write!(self.writer, "{value}")?;
        Ok(())
    }

    fn boolean(&mut self, value: bool) -> FormatResult {
        self.raw(if value { "true" } else { "false" })
    }

    fn string(&mut self, value: &str) -> FormatResult {
        self.quoted_chars(value.chars())
    }

    fn sequence<A: Alphabet>(&mut self, sequence: &Sequence<A>) -> FormatResult {
        self.quoted_chars(sequence.iter().map(|symbol| char::from(symbol.as_ascii())))
    }

    fn aligned_row<A: Alphabet>(
        &mut self,
        symbols: &[Option<align_core::Symbol<A>>],
    ) -> FormatResult {
        self.quoted_chars(
            symbols
                .iter()
                .map(|cell| cell.map_or('-', |symbol| char::from(symbol.as_ascii()))),
        )
    }

    fn quoted_chars(&mut self, characters: impl Iterator<Item = char>) -> FormatResult {
        self.raw("\"")?;
        for character in characters {
            match character {
                '"' => self.raw("\\\"")?,
                '\\' => self.raw("\\\\")?,
                '\u{08}' => self.raw("\\b")?,
                '\u{0c}' => self.raw("\\f")?,
                '\n' => self.raw("\\n")?,
                '\r' => self.raw("\\r")?,
                '\t' => self.raw("\\t")?,
                '\u{2028}' | '\u{2029}' => {
                    write!(self.writer, "\\u{:04x}", u32::from(character))?;
                }
                value if value.is_control() => {
                    write!(self.writer, "\\u{:04x}", u32::from(value))?;
                }
                value => {
                    let mut encoded = [0_u8; 4];
                    self.writer
                        .write_all(value.encode_utf8(&mut encoded).as_bytes())?;
                }
            }
        }
        self.raw("\"")
    }
}
