//! Complete-query streaming for explicit BLAST tabular dialects.

use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::iter::FusedIterator;
use std::sync::Arc;

use super::fields::{self, ParsedRow, error};
use super::{BlastComment, BlastHit, BlastQuery, BlastTabularFormat, BlastTabularOptions};
use crate::io::{BioIoError, IoErrorKind, Lines};

struct SourceLine {
    text: String,
    line: usize,
    bytes: usize,
}

struct BuildingQuery {
    query: BlastQuery,
    rows: usize,
    bytes: usize,
    fields_declared: bool,
    hit_indices: HashMap<String, usize>,
}

/// A bounded streaming iterator over complete BLAST tabular query groups.
///
/// One query and a boundary lookahead line are buffered. Completed query IDs are
/// retained under an explicit aggregate byte budget to reject noncontiguous
/// reappearances. Errors are returned once and terminate iteration. Consume to
/// EOF before trusting a footer, later groups or compression integrity.
pub struct BlastTabularReader<R: BufRead> {
    lines: Lines<R>,
    options: BlastTabularOptions,
    explicit_fields: Option<Vec<String>>,
    row_fields: Vec<String>,
    pending: Option<SourceLine>,
    current: Option<BuildingQuery>,
    preamble: Vec<BlastComment>,
    preamble_bytes: usize,
    preamble_program: Option<String>,
    preamble_version: Option<String>,
    comments: Vec<BlastComment>,
    comment_bytes: usize,
    completed: HashSet<String>,
    completed_id_bytes: usize,
    queries: usize,
    processed_queries: Option<usize>,
    done: bool,
}

impl<R: BufRead> BlastTabularReader<R> {
    /// Validates parser options without consuming input and creates an iterator.
    ///
    /// # Errors
    /// Returns invalid delimiter/field errors before reading any source bytes.
    pub fn new(reader: R, options: BlastTabularOptions) -> Result<Self, BioIoError> {
        if options.delimiter != b'\t'
            && (!(b' '..=b'~').contains(&options.delimiter) || options.delimiter == b'#')
        {
            return Err(error(
                IoErrorKind::Unsupported,
                "delimiter must be tab or printable ASCII other than '#'",
            ));
        }
        let explicit_fields = options.fields.as_deref().map(fields::fields).transpose()?;
        let row_fields = explicit_fields
            .clone()
            .unwrap_or_else(|| fields::DEFAULT_FIELDS.map(str::to_owned).to_vec());
        Ok(Self {
            lines: Lines::new(reader, options.limits.read, fields::FORMAT),
            options,
            explicit_fields,
            row_fields,
            pending: None,
            current: None,
            preamble: Vec::new(),
            preamble_bytes: 0,
            preamble_program: None,
            preamble_version: None,
            comments: Vec::new(),
            comment_bytes: 0,
            completed: HashSet::new(),
            completed_id_bytes: 0,
            queries: 0,
            processed_queries: None,
            done: false,
        })
    }

    /// Returns parsed comments outside query blocks, including a final footer.
    ///
    /// This grows only as the iterator is consumed; read to EOF for completeness.
    #[must_use]
    pub fn comments(&self) -> &[BlastComment] {
        &self.comments
    }

    /// Returns the parsed source footer count; it is absent before the footer.
    #[must_use]
    pub const fn processed_queries(&self) -> Option<usize> {
        self.processed_queries
    }

    fn next_line(&mut self) -> Result<Option<SourceLine>, BioIoError> {
        if self.pending.is_some() {
            return Ok(self.pending.take());
        }
        let text = self
            .lines
            .next_line()
            .map_err(|failure| failure.at(self.queries.max(1), self.lines.line + 1))?;
        Ok(text.map(|text| SourceLine {
            text,
            line: self.lines.line,
            bytes: usize::try_from(self.lines.byte_position() - self.lines.line_start())
                .unwrap_or(usize::MAX),
        }))
    }

    fn charge_query_bytes(&mut self, bytes: usize) -> Result<(), BioIoError> {
        if let Some(current) = &mut self.current {
            current.bytes = current
                .bytes
                .checked_add(bytes)
                .filter(|bytes| *bytes <= self.options.limits.max_query_bytes)
                .ok_or_else(|| {
                    error(
                        IoErrorKind::LimitExceeded,
                        "per-query source-byte limit exceeded",
                    )
                })?;
        }
        Ok(())
    }

    fn begin(
        &mut self,
        id: String,
        identity_field: String,
        reported_id: Option<String>,
        description: Option<String>,
        line: usize,
    ) -> Result<(), BioIoError> {
        if self.queries == self.options.limits.max_queries {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "query-count limit exceeded",
            ));
        }
        self.queries += 1;
        let commented = self.options.format == BlastTabularFormat::Commented;
        let fields = if commented {
            Vec::new()
        } else {
            self.row_fields.clone()
        };
        self.current = Some(BuildingQuery {
            query: BlastQuery {
                id,
                identity_field,
                reported_id,
                description,
                length: None,
                program: self.preamble_program.take().map(Into::into),
                version: self.preamble_version.take().map(Into::into),
                target: None,
                fields,
                source_fields: Vec::new(),
                report_fields: Arc::default(),
                comments: std::mem::take(&mut self.preamble),
                declared_rows: None,
                hits: Vec::new(),
                line,
            },
            rows: 0,
            bytes: self.preamble_bytes,
            fields_declared: !commented,
            hit_indices: HashMap::new(),
        });
        self.preamble_bytes = 0;
        Ok(())
    }

    fn finish_current(&mut self) -> Result<BlastQuery, BioIoError> {
        let current = self
            .current
            .take()
            .ok_or_else(|| error(IoErrorKind::InvalidStructure, "missing current query"))?;
        if self.options.format == BlastTabularFormat::Commented {
            let declared = current.query.declared_rows.ok_or_else(|| {
                error(
                    IoErrorKind::Truncated,
                    "query lacks a '# N hits found' declaration",
                )
            })?;
            if declared != current.rows {
                return Err(error(
                    IoErrorKind::LengthMismatch,
                    format!("declared {declared} HSP rows, observed {}", current.rows),
                ));
            }
        }
        if self.completed.contains(&current.query.id) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!(
                    "query {:?} reappears after a completed query group",
                    current.query.id
                ),
            ));
        }
        self.completed_id_bytes = self
            .completed_id_bytes
            .checked_add(current.query.id.len())
            .filter(|bytes| *bytes <= self.options.limits.max_query_id_bytes)
            .ok_or_else(|| {
                error(
                    IoErrorKind::LimitExceeded,
                    "completed-query identifier-byte limit exceeded",
                )
            })?;
        self.completed.insert(current.query.id.clone());
        Ok(current.query)
    }

    fn add_row(&mut self, row: ParsedRow) -> Result<(), BioIoError> {
        let current = self.current.as_mut().ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                "row appears before a query header",
            )
        })?;
        if current.rows == self.options.limits.max_hsps_per_query {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "per-query HSP limit exceeded",
            ));
        }
        if current.rows == 0 {
            if row.query_identity == "qseqid"
                && current
                    .query
                    .reported_id
                    .as_ref()
                    .is_some_and(|id| id != &row.query_id)
            {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "qseqid differs from query-header ID",
                ));
            }
            current.query.id = row.query_id.clone();
            current.query.identity_field = row.query_identity;
            current.query.length = row.query_length;
            current.query.source_fields = row.query_fields;
        } else if current.query.id != row.query_id
            || current.query.length != row.query_length
            || current.query.source_fields != row.query_fields
        {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "query identity/metadata changes within a query block",
            ));
        }
        if self.completed.contains(&current.query.id) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "query ID reappears after a completed group",
            ));
        }
        if current
            .query
            .declared_rows
            .is_some_and(|count| current.rows >= count)
        {
            return Err(error(
                IoErrorKind::LengthMismatch,
                "more rows than the declared HSP count",
            ));
        }
        let index = if let Some(index) = current.hit_indices.get(&row.hit_id) {
            *index
        } else {
            if current.query.hits.len() == self.options.limits.max_hits_per_query {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "per-query hit limit exceeded",
                ));
            }
            let index = current.query.hits.len();
            current.hit_indices.insert(row.hit_id.clone(), index);
            current.query.hits.push(BlastHit {
                id: row.hit_id,
                identity_field: row.hit_identity,
                description: row.hit_description.clone(),
                length: row.hit_length,
                source_fields: row.hit_fields.clone(),
                hsps: Vec::new(),
            });
            index
        };
        let hit = &mut current.query.hits[index];
        if hit.length != row.hit_length
            || hit.description != row.hit_description
            || hit.source_fields != row.hit_fields
        {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "subject metadata changes between HSP rows",
            ));
        }
        hit.hsps.push(row.hsp);
        current.rows += 1;
        Ok(())
    }

    fn preserve_global(&mut self, line: SourceLine) -> Result<(), BioIoError> {
        self.comment_bytes = self
            .comment_bytes
            .checked_add(line.bytes)
            .filter(|bytes| *bytes <= self.options.limits.max_query_bytes)
            .ok_or_else(|| {
                error(
                    IoErrorKind::LimitExceeded,
                    "document-comment byte limit exceeded",
                )
            })?;
        self.comments.push(BlastComment {
            text: line.text,
            line: line.line,
        });
        Ok(())
    }

    // One dispatch point enforces ordering for all source header declarations.
    #[allow(clippy::too_many_lines)]
    fn comment(&mut self, line: SourceLine) -> Result<Option<BlastQuery>, BioIoError> {
        let value = line
            .text
            .strip_prefix('#')
            .unwrap_or(&line.text)
            .trim_start();
        if let Some(summary) = value.strip_prefix("BLAST processed ") {
            if self.processed_queries.is_some() {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "repeated processed-query footer",
                ));
            }
            let count = summary
                .strip_suffix(" queries")
                .or_else(|| summary.strip_suffix(" query"))
                .ok_or_else(|| {
                    error(
                        IoErrorKind::InvalidStructure,
                        "malformed processed-query footer",
                    )
                })?;
            let count = usize::try_from(fields::unsigned(count, "processed query count", false)?)
                .map_err(|_| {
                error(
                    IoErrorKind::LimitExceeded,
                    "processed-query count is not addressable",
                )
            })?;
            if count != self.queries {
                return Err(error(
                    IoErrorKind::LengthMismatch,
                    "processed-query footer disagrees with query block count",
                ));
            }
            self.processed_queries = Some(count);
            return self.preserve_global(line).map(|()| None);
        }
        let program = program_header(value);
        let query_header = value.strip_prefix("Query:");
        if self.processed_queries.is_some() {
            if program.is_some() || query_header.is_some() {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "query header after processed-query footer",
                ));
            }
            return self.preserve_global(line).map(|()| None);
        }
        if self.current.is_some() && (program.is_some() || query_header.is_some()) {
            self.pending = Some(line);
            return self.finish_current().map(Some);
        }
        if let Some(header) = query_header {
            let header = header.trim();
            let (id, description) = header
                .split_once(char::is_whitespace)
                .map_or((header, None), |(id, rest)| {
                    (id, Some(rest.trim_start().to_owned()))
                });
            if id.is_empty() {
                return Err(error(IoErrorKind::InvalidHeader, "empty query-header ID"));
            }
            self.begin(
                id.to_owned(),
                "query-header".to_owned(),
                Some(id.to_owned()),
                description,
                line.line,
            )?;
        }
        if self.current.is_none() {
            if let Some((name, version)) = program {
                self.preamble_program = Some(name.to_owned());
                self.preamble_version = version.map(str::to_owned);
            }
            self.preamble_bytes = self
                .preamble_bytes
                .checked_add(line.bytes)
                .filter(|bytes| *bytes <= self.options.limits.max_query_bytes)
                .ok_or_else(|| {
                    error(
                        IoErrorKind::LimitExceeded,
                        "query-preamble byte limit exceeded",
                    )
                })?;
            self.preamble.push(BlastComment {
                text: line.text,
                line: line.line,
            });
            return Ok(None);
        }
        self.charge_query_bytes(line.bytes)?;
        let current = self
            .current
            .as_mut()
            .ok_or_else(|| error(IoErrorKind::InvalidStructure, "missing query block"))?;
        if let Some(database) = value.strip_prefix("Database:") {
            if current.rows != 0 || current.query.target.is_some() {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "repeated/late database declaration",
                ));
            }
            current.query.target = Some(database.trim_start().into());
        } else if let Some(labels) = value.strip_prefix("Fields:") {
            if current.rows != 0 || current.fields_declared {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "repeated/late field declaration",
                ));
            }
            current.query.fields = fields::header_fields(labels, self.explicit_fields.as_deref())?;
            current.fields_declared = true;
        } else if let Some(count) = value
            .strip_suffix(" hits found")
            .or_else(|| value.strip_suffix(" hit found"))
        {
            if current.rows != 0 || current.query.declared_rows.is_some() {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "repeated/late row-count declaration",
                ));
            }
            let count = usize::try_from(fields::unsigned(count, "declared HSP count", false)?)
                .map_err(|_| {
                    error(
                        IoErrorKind::LimitExceeded,
                        "declared HSP count is not addressable",
                    )
                })?;
            if count > self.options.limits.max_hsps_per_query {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "declared HSP count exceeds per-query limit",
                ));
            }
            current.query.declared_rows = Some(count);
        }
        current.query.comments.push(BlastComment {
            text: line.text,
            line: line.line,
        });
        Ok(None)
    }

    fn read_query(&mut self) -> Result<Option<BlastQuery>, BioIoError> {
        loop {
            let Some(line) = self.next_line()? else {
                if self.current.is_some() {
                    return self.finish_current().map(Some);
                }
                if self.preamble_program.is_some() {
                    return Err(error(
                        IoErrorKind::Truncated,
                        "program header has no query block",
                    ));
                }
                self.comment_bytes = self
                    .comment_bytes
                    .checked_add(self.preamble_bytes)
                    .filter(|bytes| *bytes <= self.options.limits.max_query_bytes)
                    .ok_or_else(|| {
                        error(
                            IoErrorKind::LimitExceeded,
                            "document-comment byte limit exceeded",
                        )
                    })?;
                for comment in std::mem::take(&mut self.preamble) {
                    self.comments.push(comment);
                }
                return Ok(None);
            };
            let line_number = line.line;
            let outcome = self.consume_line(line);
            match outcome {
                Ok(Some(query)) => return Ok(Some(query)),
                Ok(None) => {}
                Err(failure) => return Err(failure.at(self.queries.max(1), line_number)),
            }
        }
    }

    fn consume_line(&mut self, line: SourceLine) -> Result<Option<BlastQuery>, BioIoError> {
        if line.text.trim().is_empty() {
            self.charge_query_bytes(line.bytes)?;
            return Ok(None);
        }
        if line.text.starts_with('#') {
            if self.options.format != BlastTabularFormat::Commented {
                return Err(error(
                    IoErrorKind::Unsupported,
                    "comment line in explicit outfmt 6 input",
                ));
            }
            return self.comment(line);
        }
        if self.processed_queries.is_some() {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "data row after processed-query footer",
            ));
        }
        if self.options.format == BlastTabularFormat::Commented {
            let current = self.current.as_ref().ok_or_else(|| {
                error(
                    IoErrorKind::InvalidStructure,
                    "outfmt 7 data requires a query header",
                )
            })?;
            if !current.fields_declared || current.query.declared_rows.is_none() {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "outfmt 7 data requires Fields and hit-count declarations",
                ));
            }
            let row = fields::parse_row(
                &line.text,
                &current.query.fields,
                self.options.delimiter,
                line.line,
            )?;
            self.charge_query_bytes(line.bytes)?;
            self.add_row(row)?;
        } else {
            let row = fields::parse_row(
                &line.text,
                &self.row_fields,
                self.options.delimiter,
                line.line,
            )?;
            if self
                .current
                .as_ref()
                .is_some_and(|current| current.query.id != row.query_id)
            {
                self.pending = Some(line);
                return self.finish_current().map(Some);
            }
            if self.current.is_none() {
                self.begin(
                    row.query_id.clone(),
                    row.query_identity.clone(),
                    None,
                    None,
                    line.line,
                )?;
            }
            self.charge_query_bytes(line.bytes)?;
            self.add_row(row)?;
        }
        Ok(None)
    }
}

impl<R: BufRead> Iterator for BlastTabularReader<R> {
    type Item = Result<BlastQuery, BioIoError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.read_query() {
            Ok(Some(query)) => Some(Ok(query)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(failure) => {
                self.done = true;
                self.current = None;
                self.pending = None;
                let failure = if failure.line().is_none() {
                    failure.at(self.queries.max(1), self.lines.line.max(1))
                } else {
                    failure
                };
                Some(Err(failure))
            }
        }
    }
}

impl<R: BufRead> FusedIterator for BlastTabularReader<R> {}

fn program_header(value: &str) -> Option<(&str, Option<&str>)> {
    let (name, version) = value
        .split_once(char::is_whitespace)
        .map_or((value, None), |(name, version)| {
            (name, Some(version.trim_start()))
        });
    matches!(
        name,
        "BLASTN"
            | "BLASTP"
            | "BLASTX"
            | "TBLASTN"
            | "TBLASTX"
            | "PSIBLAST"
            | "RPSBLAST"
            | "DELTABLAST"
    )
    .then_some((name, version))
}
