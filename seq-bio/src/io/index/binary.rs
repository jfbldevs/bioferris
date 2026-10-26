//! Schema-2 binary sidecars: fixed-width record locations that stay on disk.
//!
//! The layout is specified in `docs/algorithms/binary-sequence-index-v1.md`. A
//! binary sidecar is read whole once when opened, to validate its structure and
//! checksum and to load the identifiers; afterwards only the selected record's
//! 64-byte location is read for each lookup.

use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use sha2::{Digest, Sha256};

use super::{
    BUFFER_BYTES, HashingReader, IndexAlphabet, IndexError, IndexOptions, IndexedFormat, Keys,
    Location, RecordReader, SourceState, SpanCheck, absolute_path, allocation_error, as_u64,
    at_eof, check_limit, check_unchanged, hash_section, open_regular, profile_error, read_entry,
    valid_id, verify_source,
};

/// The first eight bytes of every binary sidecar.
pub(super) const MAGIC: [u8; 8] = *b"\x89BFIDX\r\n";
const SCHEMA_VERSION: u32 = 2;
const KEY_POLICY_CODE: u8 = 1;
const HEADER_BYTES: usize = 88;
const ENTRY_BYTES: usize = 64;
const DIGEST_BYTES: usize = 32;
/// The most records a quick check hashes and parses.
pub(super) const QUICK_CHECK_RECORDS: usize = 32;

/// A source modification time as seconds and nanoseconds from the Unix epoch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Stamp {
    seconds: i64,
    nanoseconds: u32,
}

impl Stamp {
    pub(super) fn of(time: SystemTime) -> Result<Self, IndexError> {
        match time.duration_since(SystemTime::UNIX_EPOCH) {
            Ok(after) => Ok(Self {
                seconds: i64::try_from(after.as_secs()).map_err(|_| IndexError::CountOverflow)?,
                nanoseconds: after.subsec_nanos(),
            }),
            Err(error) => {
                let before = error.duration();
                let seconds =
                    i64::try_from(before.as_secs()).map_err(|_| IndexError::CountOverflow)?;
                // Nanoseconds count forward from the second, as after the epoch.
                Ok(if before.subsec_nanos() == 0 {
                    Self {
                        seconds: -seconds,
                        nanoseconds: 0,
                    }
                } else {
                    Self {
                        seconds: (-seconds).checked_sub(1).ok_or(IndexError::CountOverflow)?,
                        nanoseconds: 1_000_000_000 - before.subsec_nanos(),
                    }
                })
            }
        }
    }
}

/// The fixed 88-byte header that starts a binary sidecar.
#[derive(Debug)]
pub(super) struct Header {
    pub(super) format: IndexedFormat,
    pub(super) records: usize,
    id_bytes: u64,
    pub(super) source_len: u64,
    pub(super) source_modified: Stamp,
    pub(super) source_sha256: [u8; DIGEST_BYTES],
}

impl Header {
    fn encode(&self, molecule: &str) -> Result<[u8; HEADER_BYTES], IndexError> {
        let mut bytes = [0; HEADER_BYTES];
        bytes[..8].copy_from_slice(&MAGIC);
        bytes[8..12].copy_from_slice(&SCHEMA_VERSION.to_le_bytes());
        bytes[12] = molecule_code(molecule)?;
        bytes[13] = format_code(self.format);
        bytes[14] = KEY_POLICY_CODE;
        bytes[16..24].copy_from_slice(&as_u64(self.records)?.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.id_bytes.to_le_bytes());
        bytes[32..40].copy_from_slice(&self.source_len.to_le_bytes());
        bytes[40..48].copy_from_slice(&self.source_modified.seconds.to_le_bytes());
        bytes[48..52].copy_from_slice(&self.source_modified.nanoseconds.to_le_bytes());
        bytes[56..88].copy_from_slice(&self.source_sha256);
        Ok(bytes)
    }

    fn decode<A: IndexAlphabet>(bytes: &[u8; HEADER_BYTES]) -> Result<Self, IndexError> {
        if bytes[..8] != MAGIC || u32_at(bytes, 8) != SCHEMA_VERSION || bytes[14] != KEY_POLICY_CODE
        {
            return Err(invalid(
                "unsupported format ID, schema version or key policy",
            ));
        }
        let molecule = match bytes[12] {
            1 => "dna",
            2 => "rna",
            3 => "protein",
            _ => return Err(invalid("unknown molecule code in binary sidecar")),
        };
        if molecule != A::MOLECULE {
            return Err(IndexError::InvalidIndex(format!(
                "index molecule {molecule:?} does not match {}",
                A::MOLECULE
            )));
        }
        let format = match bytes[13] {
            1 => IndexedFormat::Fasta,
            2 => IndexedFormat::FastqSanger,
            _ => return Err(invalid("unknown source format code in binary sidecar")),
        };
        let nanoseconds = u32_at(bytes, 48);
        if bytes[15] != 0 || bytes[52..56] != [0; 4] || nanoseconds >= 1_000_000_000 {
            return Err(invalid(
                "reserved header bytes or the modification time are invalid",
            ));
        }
        let mut source_sha256 = [0; DIGEST_BYTES];
        source_sha256.copy_from_slice(&bytes[56..88]);
        Ok(Self {
            format,
            records: usize::try_from(u64_at(bytes, 16)).map_err(|_| IndexError::CountOverflow)?,
            id_bytes: u64_at(bytes, 24),
            source_len: u64_at(bytes, 32),
            source_modified: Stamp {
                seconds: i64::from_le_bytes(word(bytes, 40)),
                nanoseconds,
            },
            source_sha256,
        })
    }

    /// Header, entries, newline-terminated IDs, lookup positions and checksum.
    fn file_len(&self) -> Option<u64> {
        let records = u64::try_from(self.records).ok()?;
        let fixed = u64::try_from(HEADER_BYTES + DIGEST_BYTES).ok()?;
        let per_record = u64::try_from(ENTRY_BYTES + 1 + 8).ok()?;
        records
            .checked_mul(per_record)?
            .checked_add(self.id_bytes)?
            .checked_add(fixed)
    }
}

impl Location {
    fn encode(&self) -> [u8; ENTRY_BYTES] {
        let mut bytes = [0; ENTRY_BYTES];
        bytes[..8].copy_from_slice(&self.offset.to_le_bytes());
        bytes[8..16].copy_from_slice(&self.byte_length.to_le_bytes());
        bytes[16..24].copy_from_slice(&self.sequence_length.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.header_line.to_le_bytes());
        bytes[32..].copy_from_slice(&self.sha256);
        bytes
    }

    fn decode(bytes: &[u8; ENTRY_BYTES]) -> Self {
        let mut sha256 = [0; DIGEST_BYTES];
        sha256.copy_from_slice(&bytes[32..]);
        Self {
            offset: u64_at(bytes, 0),
            byte_length: u64_at(bytes, 8),
            sequence_length: u64_at(bytes, 16),
            header_line: u64_at(bytes, 24),
            sha256,
        }
    }
}

/// An opened binary sidecar whose locations are read on demand.
#[derive(Debug)]
pub(super) struct Sidecar {
    file: File,
    path: PathBuf,
    state: SourceState,
}

impl Sidecar {
    /// Reads the location of the record at a zero-based source position.
    pub(super) fn location(&self, position: usize) -> Result<Location, IndexError> {
        let offset = as_u64(position)?
            .checked_mul(as_u64(ENTRY_BYTES)?)
            .and_then(|offset| offset.checked_add(u64::try_from(HEADER_BYTES).ok()?))
            .ok_or(IndexError::CountOverflow)?;
        let mut file = &self.file;
        file.seek(SeekFrom::Start(offset))?;
        let mut bytes = [0; ENTRY_BYTES];
        file.read_exact(&mut bytes)?;
        Ok(Location::decode(&bytes))
    }

    /// Streams every location in source order.
    pub(super) fn locations(&self, records: usize) -> Result<Locations<'_>, IndexError> {
        let mut file = &self.file;
        file.seek(SeekFrom::Start(as_u64(HEADER_BYTES)?))?;
        Ok(Locations {
            reader: BufReader::with_capacity(BUFFER_BYTES, file),
            remaining: records,
        })
    }

    /// Rejects a sidecar whose identity, size or timestamps changed.
    pub(super) fn verify(&self) -> Result<(), IndexError> {
        check_unchanged(
            &self.file,
            &self.path,
            &self.state,
            "sidecar",
            IndexError::InvalidIndex,
        )
    }
}

/// Source-ordered locations streamed from a binary sidecar.
pub(super) struct Locations<'a> {
    reader: BufReader<&'a File>,
    remaining: usize,
}

impl Locations<'_> {
    pub(super) fn next_location(&mut self) -> Result<Option<Location>, IndexError> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let mut bytes = [0; ENTRY_BYTES];
        self.reader.read_exact(&mut bytes)?;
        self.remaining -= 1;
        Ok(Some(Location::decode(&bytes)))
    }
}

/// A validated binary sidecar: its header, identifiers and retained file.
pub(super) struct Loaded {
    pub(super) header: Header,
    pub(super) keys: Keys,
    pub(super) sidecar: Sidecar,
}

/// Returns whether a sidecar path starts with the binary magic bytes.
pub(super) fn is_binary(path: &Path) -> Result<bool, IndexError> {
    let (file, _) = open_regular(path)?;
    let mut prefix = Vec::with_capacity(MAGIC.len());
    file.take(as_u64(MAGIC.len())?).read_to_end(&mut prefix)?;
    Ok(prefix == MAGIC)
}

/// Reads and validates a whole binary sidecar, keeping only its identifiers.
///
/// Every location is checked against the declared source length and format, the
/// identifiers and their lookup order are validated, and the trailing SHA-256
/// must cover every preceding byte.
pub(super) fn load<A: IndexAlphabet>(
    path: &Path,
    options: IndexOptions,
) -> Result<Loaded, IndexError> {
    let path = absolute_path(path)?;
    let (file, state) = open_regular(&path)?;
    check_limit("index_bytes", options.max_index_bytes, state.len)?;
    if state.len < as_u64(HEADER_BYTES + DIGEST_BYTES)? {
        return Err(invalid(
            "binary sidecar is shorter than its header and checksum",
        ));
    }
    // The checksum covers every byte before the final digest.
    let body_len = state.len - as_u64(DIGEST_BYTES)?;
    let (header, keys) = {
        let mut input = HashingReader {
            reader: BufReader::with_capacity(BUFFER_BYTES, &file).take(body_len),
            digest: Sha256::new(),
        };
        let mut bytes = [0; HEADER_BYTES];
        input.read_exact(&mut bytes)?;
        let header = Header::decode::<A>(&bytes)?;
        check_limit(
            "index_records",
            options.max_records,
            as_u64(header.records)?,
        )?;
        check_limit("index_id_bytes", options.max_id_bytes, header.id_bytes)?;
        if header.file_len() != Some(state.len) {
            return Err(invalid(
                "binary sidecar length does not match its record and ID counts",
            ));
        }
        check_locations(&mut input, &header)?;
        let keys = read_keys(&mut input, &header)?;
        let HashingReader { reader, digest } = input;
        let mut reader = reader.into_inner();
        let mut checksum = [0; DIGEST_BYTES];
        reader.read_exact(&mut checksum)?;
        if digest.finalize().as_slice() != checksum || !at_eof(&mut reader)? {
            return Err(invalid(
                "binary sidecar checksum does not match its contents",
            ));
        }
        (header, keys)
    };
    check_unchanged(&file, &path, &state, "sidecar", IndexError::InvalidIndex)?;
    Ok(Loaded {
        header,
        keys,
        sidecar: Sidecar { file, path, state },
    })
}

fn check_locations(input: &mut impl Read, header: &Header) -> Result<(), IndexError> {
    let mut spans = SpanCheck::new(header.source_len, header.format);
    let mut bytes = [0; ENTRY_BYTES];
    for _ in 0..header.records {
        input.read_exact(&mut bytes)?;
        let location = Location::decode(&bytes);
        spans.check(
            location.offset,
            location.byte_length,
            location.sequence_length,
            location.header_line,
        )?;
    }
    Ok(())
}

fn read_keys(input: &mut impl Read, header: &Header) -> Result<Keys, IndexError> {
    let section = header
        .id_bytes
        .checked_add(as_u64(header.records)?)
        .ok_or(IndexError::CountOverflow)?;
    let mut text = Vec::new();
    text.try_reserve_exact(usize::try_from(section).map_err(|_| IndexError::CountOverflow)?)
        .map_err(allocation_error)?;
    input.by_ref().take(section).read_to_end(&mut text)?;
    if as_u64(text.len())? != section {
        return Err(IndexError::Io(io::ErrorKind::UnexpectedEof.into()));
    }
    let text = String::from_utf8(text).map_err(|_| invalid("record IDs must be UTF-8"))?;
    let mut ends = Vec::new();
    ends.try_reserve_exact(header.records)
        .map_err(allocation_error)?;
    let mut start = 0;
    for (end, _) in text.match_indices('\n') {
        if ends.len() == header.records {
            return Err(id_count_error());
        }
        if !valid_id(&text[start..end]) {
            return Err(invalid(
                "entry ID must be a nonempty first-token identifier",
            ));
        }
        ends.push(end);
        start = end + 1;
    }
    if ends.len() != header.records || start != text.len() {
        return Err(id_count_error());
    }
    let mut sorted = Vec::new();
    sorted
        .try_reserve_exact(header.records)
        .map_err(allocation_error)?;
    let mut word = [0; 8];
    for _ in 0..header.records {
        input.read_exact(&mut word)?;
        let position = usize::try_from(u64::from_le_bytes(word))
            .ok()
            .filter(|position| *position < header.records)
            .ok_or_else(|| invalid("lookup position is outside the entry table"))?;
        sorted.push(position);
    }
    let keys = Keys { text, ends, sorted };
    // Strictly ascending IDs over in-range positions also make the lookup a
    // permutation, so every record is reachable exactly once.
    for pair in keys.sorted.windows(2) {
        match keys.id(pair[0]).cmp(keys.id(pair[1])) {
            Ordering::Less => {}
            Ordering::Equal => return Err(invalid("sidecar contains duplicate record IDs")),
            Ordering::Greater => {
                return Err(invalid("lookup positions are not in ascending ID order"));
            }
        }
    }
    Ok(keys)
}

/// Hashes the whole source and every saved span, then reparses every record.
///
/// The source must already have the saved length. A different source digest is
/// a stale source; matching bytes with different spans, IDs or metadata mean the
/// sidecar does not describe them.
pub(super) fn verify_full<A: IndexAlphabet>(
    source: &mut File,
    source_path: &Path,
    source_state: &SourceState,
    loaded: &Loaded,
    options: IndexOptions,
) -> Result<(), IndexError> {
    let Loaded {
        header,
        keys,
        sidecar,
    } = loaded;
    let (digest, spans_match) = fingerprint_saved(
        source,
        header.source_len,
        &mut sidecar.locations(header.records)?,
    )?;
    verify_source(source, source_path, source_state)?;
    if digest.as_slice() != header.source_sha256 {
        return Err(IndexError::StaleSource(
            "source SHA-256 differs from the index".into(),
        ));
    }
    if !spans_match {
        return Err(saved_entries_differ());
    }
    source.seek(SeekFrom::Start(0))?;
    let compared = compare_records::<A>(
        BufReader::new(&mut *source),
        header,
        keys,
        &mut sidecar.locations(header.records)?,
        options,
    );
    verify_source(source, source_path, source_state)?;
    compared
}

fn fingerprint_saved(
    file: &mut File,
    source_len: u64,
    saved: &mut Locations<'_>,
) -> Result<([u8; DIGEST_BYTES], bool), IndexError> {
    file.seek(SeekFrom::Start(0))?;
    let mut source_digest = Sha256::new();
    let mut offset = 0;
    let mut buffer = [0; BUFFER_BYTES];
    let mut spans_match = true;
    while let Some(location) = saved.next_location()? {
        // Loading checked that spans are ordered, nonoverlapping and in bounds.
        let gap = location
            .offset
            .checked_sub(offset)
            .ok_or(IndexError::CountOverflow)?;
        hash_section(file, gap, &mut source_digest, None, &mut buffer)?;
        let mut span_digest = Sha256::new();
        hash_section(
            file,
            location.byte_length,
            &mut source_digest,
            Some(&mut span_digest),
            &mut buffer,
        )?;
        spans_match &= span_digest.finalize().as_slice() == location.sha256;
        offset = location
            .offset
            .checked_add(location.byte_length)
            .ok_or(IndexError::CountOverflow)?;
    }
    let tail = source_len
        .checked_sub(offset)
        .ok_or(IndexError::CountOverflow)?;
    hash_section(file, tail, &mut source_digest, None, &mut buffer)?;
    if !at_eof(file)? {
        return Err(IndexError::StaleSource(
            "source grew during fingerprinting".into(),
        ));
    }
    Ok((source_digest.finalize().into(), spans_match))
}

fn compare_records<A: IndexAlphabet>(
    source: impl io::BufRead,
    header: &Header,
    keys: &Keys,
    saved: &mut Locations<'_>,
    options: IndexOptions,
) -> Result<(), IndexError> {
    let mut reader = RecordReader::<_, A>::new(source, header.format, options.read_limits);
    let mut position = 0;
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| profile_error(error.into()))?;
        let span = reader.record_span().ok_or_else(|| {
            IndexError::InvalidIndex("reader returned a record without a source span".into())
        })?;
        let Some(location) = saved.next_location()? else {
            return Err(saved_entries_differ());
        };
        if record.id() != keys.id(position)
            || span.start != location.offset
            || span.end.checked_sub(span.start) != Some(location.byte_length)
            || as_u64(record.sequence().len())? != location.sequence_length
            || as_u64(span.header_line)? != location.header_line
        {
            return Err(saved_entries_differ());
        }
        position += 1;
    }
    if position != header.records {
        return Err(saved_entries_differ());
    }
    Ok(())
}

/// Compares the source length and modification time with the saved ones and
/// verifies the sampled records exactly as a lookup would.
pub(super) fn verify_quick<A: IndexAlphabet>(
    source: &mut File,
    source_state: &SourceState,
    loaded: &Loaded,
    options: IndexOptions,
) -> Result<(), IndexError> {
    let Loaded {
        header,
        keys,
        sidecar,
    } = loaded;
    if Stamp::of(source_state.modified)? != header.source_modified {
        return Err(IndexError::StaleSource(
            "source modification time differs from the index".into(),
        ));
    }
    for position in sampled_positions(header.records) {
        let location = sidecar.location(position)?;
        read_entry::<A>(
            source,
            keys.id(position),
            &location,
            header.format,
            options.read_limits,
        )?;
    }
    Ok(())
}

/// The first and last records and others evenly spaced by record number.
pub(super) fn sampled_positions(records: usize) -> impl Iterator<Item = usize> {
    let samples = records.min(QUICK_CHECK_RECORDS);
    let last = records.saturating_sub(1);
    let steps = samples.saturating_sub(1).max(1);
    // floor(sample * last / steps), split so the product cannot overflow.
    (0..samples).map(move |sample| last / steps * sample + last % steps * sample / steps)
}

/// The source facts a binary sidecar records besides its entries.
pub(super) struct Source<'a> {
    pub(super) format: IndexedFormat,
    pub(super) len: u64,
    pub(super) modified: SystemTime,
    pub(super) sha256: &'a [u8; DIGEST_BYTES],
}

/// Writes a complete binary sidecar from source-ordered locations.
pub(super) fn write<A: IndexAlphabet>(
    writer: &mut dyn Write,
    source: &Source<'_>,
    keys: &Keys,
    mut next_location: impl FnMut() -> Result<Option<Location>, IndexError>,
) -> Result<(), IndexError> {
    let header = Header {
        format: source.format,
        records: keys.len(),
        id_bytes: as_u64(keys.text.len() - keys.len())?,
        source_len: source.len,
        source_modified: Stamp::of(source.modified)?,
        source_sha256: *source.sha256,
    };
    let mut output = HashingWriter {
        writer,
        digest: Sha256::new(),
    };
    output.write_all(&header.encode(A::MOLECULE)?)?;
    let mut written = 0;
    while let Some(location) = next_location()? {
        output.write_all(&location.encode())?;
        written += 1;
    }
    if written != keys.len() {
        return Err(IndexError::InvalidIndex(
            "location count differs from the ID count".into(),
        ));
    }
    output.write_all(keys.text.as_bytes())?;
    for position in &keys.sorted {
        output.write_all(&as_u64(*position)?.to_le_bytes())?;
    }
    let HashingWriter { writer, digest } = output;
    writer.write_all(&digest.finalize())?;
    Ok(())
}

struct HashingWriter<W> {
    writer: W,
    digest: Sha256,
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let count = self.writer.write(bytes)?;
        self.digest.update(&bytes[..count]);
        Ok(count)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

fn molecule_code(molecule: &str) -> Result<u8, IndexError> {
    match molecule {
        "dna" => Ok(1),
        "rna" => Ok(2),
        "protein" => Ok(3),
        _ => Err(invalid("unknown molecule identifier")),
    }
}

const fn format_code(format: IndexedFormat) -> u8 {
    match format {
        IndexedFormat::Fasta => 1,
        IndexedFormat::FastqSanger => 2,
    }
}

fn u32_at(bytes: &[u8], offset: usize) -> u32 {
    let mut word = [0; 4];
    word.copy_from_slice(&bytes[offset..offset + 4]);
    u32::from_le_bytes(word)
}

fn u64_at(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(word(bytes, offset))
}

fn word(bytes: &[u8], offset: usize) -> [u8; 8] {
    let mut word = [0; 8];
    word.copy_from_slice(&bytes[offset..offset + 8]);
    word
}

fn invalid(message: &str) -> IndexError {
    IndexError::InvalidIndex(message.into())
}

fn id_count_error() -> IndexError {
    invalid("the ID section must hold one line-feed-terminated ID per record")
}

fn saved_entries_differ() -> IndexError {
    invalid("saved entries differ from the fully rederived source records")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn samples_cover_the_first_and_last_records_evenly() {
        assert_eq!(sampled_positions(0).count(), 0);
        assert_eq!(sampled_positions(1).collect::<Vec<_>>(), [0]);
        assert_eq!(sampled_positions(2).collect::<Vec<_>>(), [0, 1]);
        assert_eq!(
            sampled_positions(32).collect::<Vec<_>>(),
            (0..32).collect::<Vec<_>>()
        );
        let spread = sampled_positions(1_000).collect::<Vec<_>>();
        assert_eq!(spread.len(), 32);
        assert_eq!((spread[0], spread[1], spread[31]), (0, 32, 999));
        assert!(spread.windows(2).all(|pair| pair[0] < pair[1]));
        let huge = sampled_positions(usize::MAX).collect::<Vec<_>>();
        assert_eq!((huge[0], huge[31]), (0, usize::MAX - 1));
        assert!(huge.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn stamps_count_nanoseconds_forward_on_both_sides_of_the_epoch() {
        use std::time::Duration;

        let epoch = SystemTime::UNIX_EPOCH;
        let after = Stamp::of(epoch + Duration::new(5, 250)).unwrap();
        assert_eq!((after.seconds, after.nanoseconds), (5, 250));
        let whole = Stamp::of(epoch - Duration::new(5, 0)).unwrap();
        assert_eq!((whole.seconds, whole.nanoseconds), (-5, 0));
        let before = Stamp::of(epoch - Duration::new(5, 250)).unwrap();
        assert_eq!((before.seconds, before.nanoseconds), (-6, 999_999_750));
    }
}
