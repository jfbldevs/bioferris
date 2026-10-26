//! Tabular field names, source-preserving conversion and shared HSP checks.

use std::collections::HashSet;

use super::{BlastCoordinates, BlastDirection, BlastHsp, BlastSourceField, BlastStrand};
use crate::io::{BioIoError, IoErrorKind};

pub(super) const FORMAT: &str = "BLAST tabular";
pub(super) const DEFAULT_FIELDS: [&str; 12] = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore",
];
const QUERY_IDENTITIES: [&str; 4] = ["qseqid", "qaccver", "qacc", "qgi"];
const HIT_IDENTITIES: [&str; 4] = ["sseqid", "saccver", "sacc", "sgi"];

pub(super) fn error(kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, FORMAT, message)
}

pub(super) fn fields(input: &[String]) -> Result<Vec<String>, BioIoError> {
    let mut names = Vec::new();
    let mut seen = HashSet::new();
    for field in input {
        let expansion: Vec<&str> = if field == "std" {
            DEFAULT_FIELDS.to_vec()
        } else {
            vec![field]
        };
        for name in expansion {
            if name.len() > 256 {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "tabular field names are limited to 256 bytes",
                ));
            }
            if name.is_empty()
                || !name.as_bytes()[0].is_ascii_alphabetic()
                || !name
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
            {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "invalid tabular field name",
                ));
            }
            if names.len() == 128 {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "at most 128 tabular fields are supported",
                ));
            }
            if !seen.insert(name.to_owned()) {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    format!("duplicate field {name}"),
                ));
            }
            names.push(name.to_owned());
        }
    }
    if !QUERY_IDENTITIES.iter().any(|name| seen.contains(*name))
        || !HIT_IDENTITIES.iter().any(|name| seen.contains(*name))
    {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "fields require a query and subject identity (seqid, accver, acc or gi)",
        ));
    }
    for (start, end) in [("qstart", "qend"), ("sstart", "send")] {
        if seen.contains(start) != seen.contains(end) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("{start} and {end} must occur together"),
            ));
        }
    }
    Ok(names)
}

pub(super) fn label_name(label: &str) -> Option<&'static str> {
    match label.trim().to_ascii_lowercase().as_str() {
        "query id" => Some("qseqid"),
        "query gi" => Some("qgi"),
        "query acc." => Some("qacc"),
        "query acc.ver" => Some("qaccver"),
        "query length" => Some("qlen"),
        "subject id" => Some("sseqid"),
        "subject ids" => Some("sallseqid"),
        "subject gi" => Some("sgi"),
        "subject gis" => Some("sallgi"),
        "subject acc." => Some("sacc"),
        "subject acc.ver" => Some("saccver"),
        "subject accs." => Some("sallacc"),
        "subject length" => Some("slen"),
        "q. start" => Some("qstart"),
        "q. end" => Some("qend"),
        "s. start" => Some("sstart"),
        "s. end" => Some("send"),
        "query seq" => Some("qseq"),
        "subject seq" => Some("sseq"),
        "evalue" | "expect value" => Some("evalue"),
        "bit score" => Some("bitscore"),
        "score" => Some("score"),
        "alignment length" => Some("length"),
        "% identity" => Some("pident"),
        "identical" => Some("nident"),
        "mismatches" => Some("mismatch"),
        "positives" => Some("positive"),
        "gap opens" => Some("gapopen"),
        "gaps" => Some("gaps"),
        "% positives" => Some("ppos"),
        "query/sbjct frames" => Some("frames"),
        "query frame" => Some("qframe"),
        "sbjct frame" | "subject frame" => Some("sframe"),
        "btop" => Some("btop"),
        "subject tax ids" => Some("staxids"),
        "subject sci names" => Some("sscinames"),
        "subject com names" => Some("scomnames"),
        "subject blast names" => Some("sblastnames"),
        "subject super kingdoms" => Some("sskingdoms"),
        "subject title" => Some("stitle"),
        "subject titles" => Some("salltitles"),
        "subject strand" => Some("sstrand"),
        "% query coverage per subject" | "query coverage per subject" => Some("qcovs"),
        "% query coverage per hsp" | "query coverage per hsp" => Some("qcovhsp"),
        "% query coverage per uniq subject" | "query coverage per uniq subject" => Some("qcovus"),
        _ => None,
    }
}

pub(super) fn header_fields(
    labels: &str,
    explicit: Option<&[String]>,
) -> Result<Vec<String>, BioIoError> {
    let labels: Vec<_> = labels.split(',').map(str::trim).collect();
    if let Some(explicit) = explicit {
        if labels.len() != explicit.len() {
            return Err(error(
                IoErrorKind::LengthMismatch,
                "# Fields count differs from explicit column mapping",
            ));
        }
        for (label, name) in labels.iter().zip(explicit) {
            if label_name(label).is_some_and(|canonical| canonical != name) {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    format!("# Fields label {label:?} conflicts with {name}"),
                ));
            }
        }
        return Ok(explicit.to_vec());
    }
    let names = labels
        .iter()
        .map(|label| {
            label_name(label).map(str::to_owned).ok_or_else(|| {
                error(
                    IoErrorKind::Unsupported,
                    format!("unknown # Fields label {label:?}; provide an explicit column mapping"),
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    fields(&names)
}

pub(super) struct ParsedRow {
    pub query_id: String,
    pub query_identity: String,
    pub query_length: Option<u64>,
    pub query_fields: Vec<BlastSourceField>,
    pub hit_id: String,
    pub hit_identity: String,
    pub hit_length: Option<u64>,
    pub hit_description: Option<String>,
    pub hit_fields: Vec<BlastSourceField>,
    pub hsp: BlastHsp,
}

// The flat conversion keeps each supported source column visible in one place.
#[allow(clippy::too_many_lines)]
pub(super) fn parse_row(
    line: &str,
    names: &[String],
    delimiter: u8,
    source_line: usize,
) -> Result<ParsedRow, BioIoError> {
    if line
        .bytes()
        .any(|byte| byte.is_ascii_control() && byte != delimiter)
    {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "control byte in tabular row",
        ));
    }
    let values: Vec<_> = line.split(char::from(delimiter)).collect();
    if values.len() != names.len() {
        return Err(error(
            IoErrorKind::LengthMismatch,
            format!("expected {} columns, found {}", names.len(), values.len()),
        ));
    }
    let get = |name: &str| {
        names
            .iter()
            .position(|field| field == name)
            .map(|index| values[index])
    };
    let identity = |alternatives: &[&str]| -> Result<(String, String), BioIoError> {
        let (name, value) = alternatives
            .iter()
            .find_map(|name| get(name).map(|value| (*name, value)))
            .ok_or_else(|| error(IoErrorKind::InvalidStructure, "missing identity column"))?;
        if value.is_empty() || value.chars().any(char::is_whitespace) {
            return Err(error(
                IoErrorKind::InvalidHeader,
                format!("empty or whitespace-containing {name}"),
            ));
        }
        Ok((name.to_owned(), value.to_owned()))
    };
    let (query_identity, query_id) = identity(&QUERY_IDENTITIES)?;
    let (hit_identity, hit_id) = identity(&HIT_IDENTITIES)?;
    let query_length = get("qlen")
        .map(|value| unsigned(value, "qlen", true))
        .transpose()?;
    let hit_length = get("slen")
        .map(|value| unsigned(value, "slen", true))
        .transpose()?;
    let query_coordinates = coordinates(get("qstart"), get("qend"), query_length, "query")?;
    let hit_coordinates = coordinates(get("sstart"), get("send"), hit_length, "subject")?;
    let mut query_frame = get("qframe")
        .map(|value| frame(value, "qframe"))
        .transpose()?;
    let mut hit_frame = get("sframe")
        .map(|value| frame(value, "sframe"))
        .transpose()?;
    if let Some(frames) = get("frames") {
        let (query, subject) = frames.split_once('/').ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                "frames requires query/subject frames",
            )
        })?;
        let query = frame(query, "frames query")?;
        let subject = frame(subject, "frames subject")?;
        if query_frame.is_some_and(|existing| existing != query)
            || hit_frame.is_some_and(|existing| existing != subject)
        {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "frames disagrees with qframe/sframe",
            ));
        }
        query_frame = Some(query);
        hit_frame = Some(subject);
    }
    let query_strand = frame_strand(query_frame);
    let frame_hit_strand = frame_strand(hit_frame);
    let explicit_hit_strand = get("sstrand").map(strand).transpose()?.flatten();
    if explicit_hit_strand.is_some()
        && frame_hit_strand.is_some()
        && explicit_hit_strand != frame_hit_strand
    {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "sstrand disagrees with subject frame",
        ));
    }
    let hit_strand = explicit_hit_strand.or(frame_hit_strand);
    validate_direction(query_coordinates, query_strand, "query")?;
    validate_direction(hit_coordinates, hit_strand, "subject")?;
    let alignment_length = get("length")
        .map(|value| unsigned(value, "length", true))
        .transpose()?;
    let integer = |name| {
        get(name)
            .map(|value| unsigned(value, name, false))
            .transpose()
    };
    let percentage = |name| get(name).map(|value| percent(value, name)).transpose();
    for name in ["qcovs", "qcovhsp", "qcovus"] {
        percentage(name)?;
    }
    let hsp = BlastHsp {
        line: source_line,
        source_fields: names
            .iter()
            .zip(&values)
            .map(|(name, value)| BlastSourceField {
                name: name.clone(),
                value: (*value).to_owned(),
            })
            .collect(),
        query_coordinates,
        hit_coordinates,
        query_frame,
        hit_frame,
        query_strand,
        hit_strand,
        evalue: get("evalue")
            .map(|value| finite(value, "evalue", true))
            .transpose()?,
        bit_score: get("bitscore")
            .map(|value| finite(value, "bitscore", true))
            .transpose()?,
        raw_score: get("score")
            .map(|value| finite(value, "score", false))
            .transpose()?,
        alignment_length,
        identical: integer("nident")?,
        identity_percent: percentage("pident")?,
        mismatches: integer("mismatch")?,
        positives: integer("positive")?,
        positive_percent: percentage("ppos")?,
        gap_openings: integer("gapopen")?,
        gaps: integer("gaps")?,
        query_alignment: get("qseq").map(str::to_owned),
        hit_alignment: get("sseq").map(str::to_owned),
        midline: None,
    };
    validate_hsp(&hsp)?;
    let source = |selected: &[&str]| {
        names
            .iter()
            .zip(&values)
            .filter(|(name, _)| selected.contains(&name.as_str()))
            .map(|(name, value)| BlastSourceField {
                name: name.clone(),
                value: (*value).to_owned(),
            })
            .collect()
    };
    Ok(ParsedRow {
        query_id,
        query_identity,
        query_length,
        query_fields: source(&["qseqid", "qgi", "qacc", "qaccver", "qlen"]),
        hit_id,
        hit_identity,
        hit_length,
        hit_description: get("stitle").map(str::to_owned),
        hit_fields: source(&[
            "sseqid",
            "sallseqid",
            "sgi",
            "sallgi",
            "sacc",
            "saccver",
            "sallacc",
            "slen",
            "stitle",
            "salltitles",
            "staxids",
            "sscinames",
            "scomnames",
            "sblastnames",
            "sskingdoms",
            "qcovs",
            "qcovus",
        ]),
        hsp,
    })
}

pub(super) fn unsigned(value: &str, name: &str, positive: bool) -> Result<u64, BioIoError> {
    let parsed = if !value.is_empty() && value.bytes().all(|byte| byte.is_ascii_digit()) {
        value.parse::<u64>().ok()
    } else {
        None
    };
    parsed
        .filter(|value| !positive || *value > 0)
        .ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                format!(
                    "{name} must be a {}unsigned 64-bit integer",
                    if positive {
                        "positive "
                    } else {
                        "nonnegative "
                    }
                ),
            )
        })
}

pub(super) fn finite(value: &str, name: &str, nonnegative: bool) -> Result<f64, BioIoError> {
    value
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && (!nonnegative || *value >= 0.0))
        .ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                format!(
                    "{name} must be a finite {}number",
                    if nonnegative { "nonnegative " } else { "" }
                ),
            )
        })
}

fn percent(value: &str, name: &str) -> Result<f64, BioIoError> {
    let value = finite(value, name, true)?;
    if value > 100.0 {
        return Err(error(
            IoErrorKind::InvalidStructure,
            format!("{name} must lie in 0..100"),
        ));
    }
    Ok(value)
}

pub(super) fn frame(value: &str, name: &str) -> Result<i8, BioIoError> {
    value
        .parse::<i8>()
        .ok()
        .filter(|value| (-3..=3).contains(value))
        .ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                format!("{name} must be an integer in -3..3"),
            )
        })
}

pub(super) fn strand(value: &str) -> Result<Option<BlastStrand>, BioIoError> {
    match value.to_ascii_lowercase().as_str() {
        "plus" => Ok(Some(BlastStrand::Plus)),
        "minus" => Ok(Some(BlastStrand::Minus)),
        "n/a" => Ok(None),
        _ => Err(error(
            IoErrorKind::InvalidStructure,
            "strand must be plus, minus or N/A",
        )),
    }
}

pub(super) fn frame_strand(frame: Option<i8>) -> Option<BlastStrand> {
    match frame {
        Some(1..=3) => Some(BlastStrand::Plus),
        Some(-3..=-1) => Some(BlastStrand::Minus),
        _ => None,
    }
}

pub(super) fn coordinates(
    start: Option<&str>,
    end: Option<&str>,
    length: Option<u64>,
    name: &str,
) -> Result<Option<BlastCoordinates>, BioIoError> {
    match (start, end) {
        (None, None) => Ok(None),
        (Some(start), Some(end)) => {
            let reported_start = unsigned(start, name, true)?;
            let reported_end = unsigned(end, name, true)?;
            if length.is_some_and(|length| reported_start > length || reported_end > length) {
                return Err(error(
                    IoErrorKind::LengthMismatch,
                    format!("{name} coordinates exceed declared sequence length"),
                ));
            }
            Ok(Some(BlastCoordinates {
                reported_start,
                reported_end,
            }))
        }
        _ => Err(error(
            IoErrorKind::InvalidStructure,
            format!("incomplete {name} coordinate pair"),
        )),
    }
}

pub(super) fn validate_direction(
    coordinates: Option<BlastCoordinates>,
    strand: Option<BlastStrand>,
    name: &str,
) -> Result<(), BioIoError> {
    if let (Some(coordinates), Some(strand)) = (coordinates, strand) {
        if matches!(
            (coordinates.direction(), strand),
            (BlastDirection::Increasing, BlastStrand::Minus)
                | (BlastDirection::Decreasing, BlastStrand::Plus)
        ) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("{name} coordinate direction contradicts strand/frame"),
            ));
        }
    }
    Ok(())
}

pub(super) fn validate_hsp(hsp: &BlastHsp) -> Result<(), BioIoError> {
    if let Some(length) = hsp.alignment_length {
        for count in [
            hsp.identical,
            hsp.mismatches,
            hsp.positives,
            hsp.gap_openings,
            hsp.gaps,
        ]
        .into_iter()
        .flatten()
        {
            if count > length {
                return Err(error(
                    IoErrorKind::LengthMismatch,
                    "HSP count exceeds alignment length",
                ));
            }
        }
    }
    if let (Some(length), Some(identical), Some(mismatches), Some(gaps)) = (
        hsp.alignment_length,
        hsp.identical,
        hsp.mismatches,
        hsp.gaps,
    ) {
        if u128::from(identical) + u128::from(mismatches) + u128::from(gaps) != u128::from(length) {
            return Err(error(
                IoErrorKind::LengthMismatch,
                "identical, mismatch and gap counts do not sum to alignment length",
            ));
        }
    }
    if matches!((hsp.identical, hsp.positives), (Some(identical), Some(positive)) if identical > positive)
    {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "identities exceed positives",
        ));
    }
    if matches!((hsp.gap_openings, hsp.gaps), (Some(openings), Some(gaps)) if openings > gaps) {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "gap openings exceed gap count",
        ));
    }
    for alignment in [hsp.query_alignment.as_deref(), hsp.hit_alignment.as_deref()]
        .into_iter()
        .flatten()
    {
        if alignment.is_empty()
            || !alignment
                .bytes()
                .all(|byte| byte.is_ascii_alphabetic() || matches!(byte, b'*' | b'-'))
        {
            return Err(error(
                IoErrorKind::InvalidSequence,
                "alignment strings require ASCII letters, '*' or '-'",
            ));
        }
        if hsp
            .alignment_length
            .is_some_and(|length| u64::try_from(alignment.len()).ok() != Some(length))
        {
            return Err(error(
                IoErrorKind::LengthMismatch,
                "alignment string length disagrees with reported length",
            ));
        }
    }
    if matches!((&hsp.query_alignment, &hsp.hit_alignment), (Some(query), Some(hit)) if query.len() != hit.len())
    {
        return Err(error(
            IoErrorKind::LengthMismatch,
            "query and subject alignment strings differ in length",
        ));
    }
    Ok(())
}
