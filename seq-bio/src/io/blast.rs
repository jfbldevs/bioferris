//! Bounded BLAST result import with source fields and original coordinates.
//!
//! Import preserves reported search results; it does not run BLAST. Tabular
//! profiles are explicit, omitted values stay absent, and source strings remain
//! available alongside their checked numerical interpretation.

use std::io::BufRead;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;

use super::{BioIoError, ReadLimits, open_input};

#[path = "blast/fields.rs"]
mod fields;
#[path = "blast/tabular.rs"]
mod tabular;
pub use tabular::BlastTabularReader;
#[path = "blast/xml.rs"]
mod xml;
pub use xml::BlastXmlReader;

#[cfg(test)]
#[path = "blast/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "blast/xml_tests.rs"]
mod xml_tests;

/// An explicitly selected BLAST tabular input profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlastTabularFormat {
    /// Outfmt 6: rows only, with caller-selected columns or the standard twelve.
    Uncommented,
    /// Outfmt 7: comment-delimited query blocks with declared columns/counts.
    Commented,
}

/// Explicit source, grouping and metadata budgets for BLAST input.
///
/// These cap source bytes and retained entry counts, not exact process memory.
#[derive(Clone, Copy, Debug)]
pub struct BlastReadLimits {
    /// Physical-line and total decoded-byte limits; `max_record_symbols` is unused.
    pub read: ReadLimits,
    /// Maximum query blocks, including zero-hit queries; default 100,000.
    pub max_queries: usize,
    /// Maximum distinct subject hits in one query; default 100,000.
    pub max_hits_per_query: usize,
    /// Maximum HSP rows in one query; default 1,000,000.
    pub max_hsps_per_query: usize,
    /// Maximum source bytes retained/processed for one query; default 64 MiB.
    pub max_query_bytes: usize,
    /// Maximum cumulative bytes in retained completed-query IDs; default 16 MiB.
    pub max_query_id_bytes: usize,
}

impl Default for BlastReadLimits {
    fn default() -> Self {
        Self {
            read: ReadLimits::default(),
            max_queries: 100_000,
            max_hits_per_query: 100_000,
            max_hsps_per_query: 1_000_000,
            max_query_bytes: 64 * 1024 * 1024,
            max_query_id_bytes: 16 * 1024 * 1024,
        }
    }
}

/// Tabular dialect, optional column mapping and explicit parser limits.
#[derive(Clone, Debug)]
pub struct BlastTabularOptions {
    /// Explicit uncommented/commented format; default uncommented outfmt 6.
    pub format: BlastTabularFormat,
    /// Explicit source column names; `None` uses default twelve or a comment header.
    ///
    /// `std` expands to the standard twelve. Names must be unique after expansion;
    /// at most 128 columns are accepted. Unknown custom columns remain raw strings.
    pub fields: Option<Vec<String>>,
    /// Literal one-byte separator; default tab, or printable ASCII other than `#`.
    pub delimiter: u8,
    /// Source and retained-entry budgets.
    pub limits: BlastReadLimits,
}

impl Default for BlastTabularOptions {
    fn default() -> Self {
        Self {
            format: BlastTabularFormat::Uncommented,
            fields: None,
            delimiter: b'\t',
            limits: BlastReadLimits::default(),
        }
    }
}

/// An explicitly selected NCBI XML result profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlastXmlFormat {
    /// Outfmt 5, following the legacy `NCBI_BlastOutput.mod.dtd` names.
    Legacy,
    /// Single-file outfmt 16, following the namespace-aware `schema_alt` names.
    Xml2,
}

/// XML profile and explicit transport/event/tree/metadata budgets.
#[derive(Clone, Copy, Debug)]
pub struct BlastXmlOptions {
    /// Explicit XML profile; the default is legacy outfmt 5.
    pub format: BlastXmlFormat,
    /// Shared decoded-source and query/hit/HSP limits.
    pub limits: BlastReadLimits,
    /// Maximum total XML nesting depth, including root containers; default 64.
    pub max_depth: usize,
    /// Maximum source bytes consumed in one XML event; default 8 MiB.
    pub max_event_bytes: usize,
    /// Maximum retained document/report metadata bytes; default 16 MiB.
    pub max_metadata_bytes: usize,
}

impl Default for BlastXmlOptions {
    fn default() -> Self {
        Self {
            format: BlastXmlFormat::Legacy,
            limits: BlastReadLimits::default(),
            max_depth: 64,
            max_event_bytes: 8 * 1024 * 1024,
            max_metadata_bytes: 16 * 1024 * 1024,
        }
    }
}

/// Original coordinate ordering, independent of biological molecule/strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlastDirection {
    /// Original start is below original end.
    Increasing,
    /// Original start is above original end.
    Decreasing,
    /// Original start and end are equal; no direction follows from their order.
    Point,
}

/// A biological strand explicitly reported or implied by a nonzero source frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlastStrand {
    /// Plus strand.
    Plus,
    /// Minus strand.
    Minus,
}

/// Checked original one-based inclusive endpoints and their normalized interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BlastCoordinates {
    pub(crate) reported_start: u64,
    pub(crate) reported_end: u64,
}

impl BlastCoordinates {
    /// Returns ascending zero-based half-open coordinates in source sequence units.
    #[must_use]
    pub fn range(&self) -> Range<u64> {
        self.reported_start.min(self.reported_end) - 1..self.reported_start.max(self.reported_end)
    }

    /// Returns the original positive one-based inclusive starting coordinate.
    #[must_use]
    pub const fn reported_start(&self) -> u64 {
        self.reported_start
    }

    /// Returns the original positive one-based inclusive ending coordinate.
    #[must_use]
    pub const fn reported_end(&self) -> u64 {
        self.reported_end
    }

    /// Returns source ordering without claiming a strand for a protein axis.
    #[must_use]
    pub const fn direction(&self) -> BlastDirection {
        if self.reported_start < self.reported_end {
            BlastDirection::Increasing
        } else if self.reported_start > self.reported_end {
            BlastDirection::Decreasing
        } else {
            BlastDirection::Point
        }
    }
}

/// One source field, retained as original text with its canonical name or XML path.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BlastSourceField {
    pub(crate) name: String,
    pub(crate) value: String,
}

impl BlastSourceField {
    /// Returns the selected column name or explicit source metadata path.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the unmodified source value, including numeric lexical precision.
    #[must_use]
    pub fn value(&self) -> &str {
        &self.value
    }
}

/// A preserved source comment without its physical line ending.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BlastComment {
    pub(crate) text: String,
    pub(crate) line: usize,
}

impl BlastComment {
    /// Returns the original comment text, including the tabular `#` prefix.
    #[must_use]
    pub fn text(&self) -> &str {
        &self.text
    }

    /// Returns the one-based physical source line.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

macro_rules! optional_getters {
    ($($name:ident : $ty:ty => $documentation:literal),* $(,)?) => {
        $(#[doc = $documentation]
        #[must_use]
        pub const fn $name(&self) -> Option<$ty> {
            self.$name
        })*
    };
}

/// One source HSP, with optional typed statistics and every original source field.
#[derive(Clone, Debug, PartialEq)]
pub struct BlastHsp {
    pub(crate) line: usize,
    pub(crate) source_fields: Vec<BlastSourceField>,
    pub(crate) query_coordinates: Option<BlastCoordinates>,
    pub(crate) hit_coordinates: Option<BlastCoordinates>,
    pub(crate) query_frame: Option<i8>,
    pub(crate) hit_frame: Option<i8>,
    pub(crate) query_strand: Option<BlastStrand>,
    pub(crate) hit_strand: Option<BlastStrand>,
    pub(crate) evalue: Option<f64>,
    pub(crate) bit_score: Option<f64>,
    pub(crate) raw_score: Option<f64>,
    pub(crate) alignment_length: Option<u64>,
    pub(crate) identical: Option<u64>,
    pub(crate) identity_percent: Option<f64>,
    pub(crate) mismatches: Option<u64>,
    pub(crate) positives: Option<u64>,
    pub(crate) positive_percent: Option<f64>,
    pub(crate) gap_openings: Option<u64>,
    pub(crate) gaps: Option<u64>,
    pub(crate) query_alignment: Option<String>,
    pub(crate) hit_alignment: Option<String>,
    pub(crate) midline: Option<String>,
}

impl BlastHsp {
    optional_getters! {
        query_coordinates: BlastCoordinates => "Returns original and normalized query coordinates when present.",
        hit_coordinates: BlastCoordinates => "Returns original and normalized subject coordinates when present.",
        query_frame: i8 => "Returns the explicitly supplied query frame, including zero.",
        hit_frame: i8 => "Returns the explicitly supplied subject frame, including zero.",
        query_strand: BlastStrand => "Returns a reported query strand or one implied by a nonzero frame.",
        hit_strand: BlastStrand => "Returns a reported subject strand or one implied by a nonzero frame.",
        evalue: f64 => "Returns the finite nonnegative e-value, or None when omitted; source text is retained.",
        bit_score: f64 => "Returns the finite nonnegative bit score, or None when omitted.",
        raw_score: f64 => "Returns the finite raw score, or None when omitted.",
        alignment_length: u64 => "Returns the explicitly reported alignment column count, including gaps.",
        identical: u64 => "Returns the explicitly reported number of identical aligned residues.",
        identity_percent: f64 => "Returns the reported percentage identity in 0..100 without recomputation.",
        mismatches: u64 => "Returns the explicitly reported mismatch count.",
        positives: u64 => "Returns the explicitly reported positive-scoring residue count.",
        positive_percent: f64 => "Returns the reported percentage positives in 0..100.",
        gap_openings: u64 => "Returns the explicitly reported number of gap openings.",
        gaps: u64 => "Returns the explicitly reported number of gap columns."
    }

    /// Returns the one-based source line of the HSP row or opening element.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }

    /// Returns source values in their original column/metadata order.
    #[must_use]
    pub fn source_fields(&self) -> &[BlastSourceField] {
        &self.source_fields
    }

    /// Looks up an exact source column/path name without interpreting its value.
    #[must_use]
    pub fn field(&self, name: &str) -> Option<&str> {
        self.source_fields
            .iter()
            .find(|field| field.name == name)
            .map(|field| field.value.as_str())
    }

    /// Returns the source query alignment string with original case and gaps.
    #[must_use]
    pub fn query_alignment(&self) -> Option<&str> {
        self.query_alignment.as_deref()
    }

    /// Returns the source subject alignment string with original case and gaps.
    #[must_use]
    pub fn hit_alignment(&self) -> Option<&str> {
        self.hit_alignment.as_deref()
    }

    /// Returns a supplied XML alignment midline, preserving spaces.
    #[must_use]
    pub fn midline(&self) -> Option<&str> {
        self.midline.as_deref()
    }
}

/// One subject hit, containing HSPs in source order.
#[derive(Clone, Debug, PartialEq)]
pub struct BlastHit {
    pub(crate) id: String,
    pub(crate) identity_field: String,
    pub(crate) description: Option<String>,
    pub(crate) length: Option<u64>,
    pub(crate) source_fields: Vec<BlastSourceField>,
    pub(crate) hsps: Vec<BlastHsp>,
}

impl BlastHit {
    /// Returns the source subject identity; no generated-ID heuristic is applied.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Returns the source column/XML path selected for this identity.
    #[must_use]
    pub fn identity_field(&self) -> &str {
        &self.identity_field
    }
    /// Returns a source description when explicitly present.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Returns the explicitly supplied complete subject length.
    #[must_use]
    pub const fn length(&self) -> Option<u64> {
        self.length
    }
    /// Returns source metadata, retaining multiple XML descriptions separately.
    #[must_use]
    pub fn source_fields(&self) -> &[BlastSourceField] {
        &self.source_fields
    }
    /// Returns all source HSPs, in source occurrence order within this subject.
    #[must_use]
    pub fn hsps(&self) -> &[BlastHsp] {
        &self.hsps
    }
}

/// One complete query group with source metadata and ordered subject hits.
#[derive(Clone, Debug, PartialEq)]
pub struct BlastQuery {
    pub(crate) id: String,
    pub(crate) identity_field: String,
    pub(crate) reported_id: Option<String>,
    pub(crate) description: Option<String>,
    pub(crate) length: Option<u64>,
    pub(crate) program: Option<Arc<str>>,
    pub(crate) version: Option<Arc<str>>,
    pub(crate) target: Option<Arc<str>>,
    pub(crate) fields: Vec<String>,
    pub(crate) source_fields: Vec<BlastSourceField>,
    pub(crate) report_fields: Arc<[BlastSourceField]>,
    pub(crate) comments: Vec<BlastComment>,
    pub(crate) declared_rows: Option<usize>,
    pub(crate) hits: Vec<BlastHit>,
    pub(crate) line: usize,
}

impl BlastQuery {
    /// Returns the chosen source query identity, without accession reconstruction.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Returns the column/XML path selected for identity, or `query-header` for no hits.
    #[must_use]
    pub fn identity_field(&self) -> &str {
        &self.identity_field
    }
    /// Returns a separately supplied raw query-header ID, including generated IDs.
    #[must_use]
    pub fn reported_id(&self) -> Option<&str> {
        self.reported_id.as_deref()
    }
    /// Returns the source query description, if present.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Returns the explicitly supplied complete query length.
    #[must_use]
    pub const fn length(&self) -> Option<u64> {
        self.length
    }
    /// Returns a source program name, if supplied; it is not inferred from rows.
    #[must_use]
    pub fn program(&self) -> Option<&str> {
        self.program.as_deref()
    }
    /// Returns the original program version string, if supplied.
    #[must_use]
    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }
    /// Returns the source target database text, if supplied.
    #[must_use]
    pub fn target(&self) -> Option<&str> {
        self.target.as_deref()
    }
    /// Returns canonical tabular column names in source order; XML uses source paths.
    #[must_use]
    pub fn fields(&self) -> &[String] {
        &self.fields
    }
    /// Returns source metadata in source order.
    #[must_use]
    pub fn source_fields(&self) -> &[BlastSourceField] {
        &self.source_fields
    }
    /// Returns shared XML report metadata such as parameters and source references.
    ///
    /// Sharing keeps a large report header from being cloned for every query.
    /// Tabular queries have no separate XML report fields.
    #[must_use]
    pub fn report_fields(&self) -> &[BlastSourceField] {
        &self.report_fields
    }
    /// Returns all comments assigned to this query, including unrecognized comments.
    #[must_use]
    pub fn comments(&self) -> &[BlastComment] {
        &self.comments
    }
    /// Returns the source-declared HSP row count, if supplied.
    #[must_use]
    pub const fn declared_rows(&self) -> Option<usize> {
        self.declared_rows
    }
    /// Returns hits in their first-occurrence order.
    #[must_use]
    pub fn hits(&self) -> &[BlastHit] {
        &self.hits
    }
    /// Returns the one-based physical line where this query group began.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

/// An eager collected result retaining document-level comments and metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct BlastDocument {
    pub(crate) queries: Vec<BlastQuery>,
    pub(crate) comments: Vec<BlastComment>,
    pub(crate) source_fields: Vec<BlastSourceField>,
    pub(crate) processed_queries: Option<usize>,
}

impl BlastDocument {
    /// Returns all query groups in source order, including no-hit queries.
    #[must_use]
    pub fn queries(&self) -> &[BlastQuery] {
        &self.queries
    }
    /// Returns comments outside query groups, including a processed-query footer.
    #[must_use]
    pub fn comments(&self) -> &[BlastComment] {
        &self.comments
    }
    /// Returns document-level source metadata, such as XML declarations/statistics.
    #[must_use]
    pub fn source_fields(&self) -> &[BlastSourceField] {
        &self.source_fields
    }
    /// Returns the explicitly reported total processed-query count, if present.
    #[must_use]
    pub const fn processed_queries(&self) -> Option<usize> {
        self.processed_queries
    }
}

/// Collects a bounded tabular input, retaining all queries and document comments.
///
/// # Errors
/// Returns the first configuration, transport, syntax, consistency or limit error.
/// On error no partial document is returned. Use [`BlastTabularReader`] to consume
/// complete earlier queries explicitly while retaining later failure information.
pub fn parse_blast_tabular(
    reader: impl BufRead,
    options: BlastTabularOptions,
) -> Result<BlastDocument, BioIoError> {
    let mut reader = BlastTabularReader::new(reader, options)?;
    let queries = reader.by_ref().collect::<Result<Vec<_>, _>>()?;
    Ok(BlastDocument {
        queries,
        comments: reader.comments().to_vec(),
        source_fields: Vec::new(),
        processed_queries: reader.processed_queries(),
    })
}

/// Opens a named plain/gzip source and collects its explicit tabular profile.
///
/// # Errors
/// Returns open/decompression errors or the first parser/configuration/limit error.
/// The whole input is consumed, including gzip trailers and query-count footers.
pub fn read_blast_tabular(
    path: impl AsRef<Path>,
    options: BlastTabularOptions,
) -> Result<BlastDocument, BioIoError> {
    parse_blast_tabular(open_input(path)?, options)
}

/// Collects one complete, explicitly selected NCBI XML result document.
///
/// # Errors
/// Returns malformed/unsupported XML, failed model checks, transport errors or
/// exhausted budgets. No partial document is returned. XML parsing never fetches
/// DTDs, expands declared entities or performs `XInclude`.
pub fn parse_blast_xml(
    reader: impl BufRead,
    options: BlastXmlOptions,
) -> Result<BlastDocument, BioIoError> {
    let mut reader = BlastXmlReader::new(reader, options)?;
    let queries = reader.by_ref().collect::<Result<Vec<_>, _>>()?;
    Ok(BlastDocument {
        queries,
        comments: reader.comments().to_vec(),
        source_fields: reader.source_fields().to_vec(),
        processed_queries: None,
    })
}

/// Opens a named plain/gzip source and collects its explicit NCBI XML profile.
///
/// # Errors
/// Returns open/decompression errors or the first parser/model/limit error.
/// Consume succeeds only after the complete XML document and transport are read.
pub fn read_blast_xml(
    path: impl AsRef<Path>,
    options: BlastXmlOptions,
) -> Result<BlastDocument, BioIoError> {
    parse_blast_xml(open_input(path)?, options)
}
