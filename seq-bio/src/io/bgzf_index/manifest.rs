use std::fmt;
use std::io::{Read, Write};

use super::{
    BTreeMap, BgzfIndexError, BgzfIndexOptions, BgzfVirtualOffset, Deserialize, FORMAT_ID,
    IndexAlphabet, IndexError, IndexedFormat, KEY_POLICY, Path, PhantomData, Serialize,
    allocation_error, as_u64, check_limit, invalid, io, open_regular, verify_source,
};

#[derive(Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(deny_unknown_fields)]
pub(super) struct Manifest {
    pub format_id: String,
    pub schema_version: u32,
    pub molecule: String,
    pub source_format: IndexedFormat,
    pub key_policy: String,
    pub source_len: u64,
    pub source_sha256: String,
    pub uncompressed_len: u64,
    #[serde(deserialize_with = "objects")]
    pub blocks: Vec<BlockEntry>,
    #[serde(deserialize_with = "objects")]
    pub entries: Vec<Entry>,
}

#[derive(Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(deny_unknown_fields)]
pub(super) struct BlockEntry {
    pub compressed_offset: u64,
    pub compressed_size: u32,
    pub uncompressed_offset: u64,
    pub uncompressed_size: u32,
    pub canonical_eof_marker: bool,
}

#[derive(Debug, Deserialize, Serialize, Eq, PartialEq)]
#[serde(deny_unknown_fields)]
pub(super) struct Entry {
    pub id: String,
    pub offset: u64,
    pub byte_length: u64,
    pub virtual_start: u64,
    pub virtual_end: u64,
    pub sequence_length: u64,
    pub record: usize,
    pub header_line: usize,
    pub sha256: String,
}

impl Manifest {
    pub fn virtual_offset(&self, offset: u64) -> Result<BgzfVirtualOffset, BgzfIndexError> {
        if offset > self.uncompressed_len {
            return Err(invalid("decompressed position exceeds source size"));
        }
        let next = self
            .blocks
            .partition_point(|block| block.uncompressed_offset < offset);
        if let Some(block) = self
            .blocks
            .get(next)
            .filter(|block| block.uncompressed_offset == offset)
        {
            return BgzfVirtualOffset::new(block.compressed_offset, 0).map_err(Into::into);
        }
        let Some(block) = next.checked_sub(1).and_then(|index| self.blocks.get(index)) else {
            return Err(invalid("no verified block contains decompressed position"));
        };
        let within = offset
            .checked_sub(block.uncompressed_offset)
            .ok_or(IndexError::CountOverflow)?;
        if within >= u64::from(block.uncompressed_size) {
            return Err(invalid(
                "decompressed position falls outside a verified block",
            ));
        }
        BgzfVirtualOffset::new(
            block.compressed_offset,
            u32::try_from(within).map_err(|_| IndexError::CountOverflow)?,
        )
        .map_err(Into::into)
    }

    pub fn uncompressed_offset(&self, offset: BgzfVirtualOffset) -> Result<u64, BgzfIndexError> {
        if offset.compressed_offset() == self.source_len && offset.uncompressed_offset() == 0 {
            return Ok(self.uncompressed_len);
        }
        let found = self
            .blocks
            .binary_search_by_key(&offset.compressed_offset(), |block| block.compressed_offset)
            .map_err(|_| invalid("virtual compressed address is not a verified block boundary"))?;
        let block = &self.blocks[found];
        if u32::from(offset.uncompressed_offset()) > block.uncompressed_size {
            return Err(invalid("virtual data offset exceeds its verified block"));
        }
        block
            .uncompressed_offset
            .checked_add(u64::from(offset.uncompressed_offset()))
            .ok_or_else(|| IndexError::CountOverflow.into())
    }

    pub fn validate<A: IndexAlphabet>(
        &self,
        options: BgzfIndexOptions,
    ) -> Result<(), BgzfIndexError> {
        if self.format_id != FORMAT_ID
            || self.schema_version != 1
            || self.key_policy != KEY_POLICY
            || self.molecule != A::MOLECULE
        {
            return Err(invalid(
                "unsupported format, version, key policy or molecule",
            ));
        }
        digest(&self.source_sha256)?;
        check_limit(
            "compressed_source_bytes",
            options.max_compressed_bytes,
            self.source_len,
        )?;
        check_limit(
            "source_bytes",
            as_u64(options.index.read_limits.max_total_bytes)?,
            self.uncompressed_len,
        )?;
        check_limit(
            "bgzf_blocks",
            as_u64(options.max_blocks)?,
            as_u64(self.blocks.len())?,
        )?;
        check_limit(
            "index_records",
            as_u64(options.index.max_records)?,
            as_u64(self.entries.len())?,
        )?;
        self.validate_blocks()?;
        self.validate_entries(options)?;
        check_limit(
            "index_bytes",
            as_u64(options.index.max_index_bytes)?,
            serialized_size(self)?,
        )
    }

    fn validate_blocks(&self) -> Result<(), BgzfIndexError> {
        let mut compressed = 0_u64;
        let mut uncompressed = 0_u64;
        for block in &self.blocks {
            if block.compressed_offset != compressed
                || block.uncompressed_offset != uncompressed
                || !(26..=65_536).contains(&block.compressed_size)
                || block.uncompressed_size > 65_536
                || (block.canonical_eof_marker
                    && (block.compressed_size != 28 || block.uncompressed_size != 0))
            {
                return Err(invalid("block table coordinates or sizes are inconsistent"));
            }
            BgzfVirtualOffset::new(block.compressed_offset, 0)?;
            compressed = compressed
                .checked_add(u64::from(block.compressed_size))
                .ok_or(IndexError::CountOverflow)?;
            uncompressed = uncompressed
                .checked_add(u64::from(block.uncompressed_size))
                .ok_or(IndexError::CountOverflow)?;
        }
        BgzfVirtualOffset::new(compressed, 0)?;
        if compressed != self.source_len
            || uncompressed != self.uncompressed_len
            || self
                .blocks
                .last()
                .is_none_or(|block| block.uncompressed_size != 0)
        {
            return Err(invalid(
                "block table does not cover the source and final empty member",
            ));
        }
        Ok(())
    }

    fn validate_entries(&self, options: BgzfIndexOptions) -> Result<(), BgzfIndexError> {
        let mut end = 0_u64;
        let mut line = 0;
        let mut id_bytes = 0_u64;
        let mut ids = BTreeMap::new();
        for (index, entry) in self.entries.iter().enumerate() {
            let next = entry
                .offset
                .checked_add(entry.byte_length)
                .ok_or_else(|| invalid("record span overflows"))?;
            if entry.id.is_empty()
                || entry
                    .id
                    .chars()
                    .any(|character| character.is_whitespace() || character.is_control())
                || entry.record != index.checked_add(1).ok_or(IndexError::CountOverflow)?
                || entry.header_line <= line
                || entry.offset < end
                || entry.byte_length == 0
                || next > self.uncompressed_len
                || entry.sequence_length > entry.byte_length
                || (self.source_format == IndexedFormat::FastqSanger && entry.sequence_length == 0)
            {
                return Err(invalid(
                    "record coordinates, identity or order are inconsistent",
                ));
            }
            if self.virtual_offset(entry.offset)?.raw() != entry.virtual_start
                || self.virtual_offset(next)?.raw() != entry.virtual_end
            {
                return Err(invalid(
                    "record virtual position is not its canonical verified mapping",
                ));
            }
            digest(&entry.sha256)?;
            id_bytes = id_bytes
                .checked_add(as_u64(entry.id.len())?)
                .ok_or(IndexError::CountOverflow)?;
            check_limit(
                "index_id_bytes",
                as_u64(options.index.max_id_bytes)?,
                id_bytes,
            )?;
            if ids.insert(&entry.id, entry.record).is_some() {
                return Err(invalid("sidecar contains duplicate record keys"));
            }
            end = next;
            line = entry.header_line;
        }
        Ok(())
    }
}

fn digest(value: &str) -> Result<(), BgzfIndexError> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(invalid(
            "SHA-256 must contain exactly 64 lowercase hexadecimal characters",
        ));
    }
    Ok(())
}

struct ObjectOnly<T>(T);
impl<'de, T: Deserialize<'de>> Deserialize<'de> for ObjectOnly<T> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Visitor<T>(PhantomData<T>);
        impl<'de, T: Deserialize<'de>> serde::de::Visitor<'de> for Visitor<T> {
            type Value = ObjectOnly<T>;
            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a BGZF-index JSON object")
            }
            fn visit_map<M: serde::de::MapAccess<'de>>(
                self,
                map: M,
            ) -> Result<Self::Value, M::Error> {
                T::deserialize(serde::de::value::MapAccessDeserializer::new(map)).map(ObjectOnly)
            }
        }
        deserializer.deserialize_map(Visitor(PhantomData))
    }
}

fn objects<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Vec<T>, D::Error> {
    Vec::<ObjectOnly<T>>::deserialize(deserializer)
        .map(|objects| objects.into_iter().map(|object| object.0).collect())
}

pub(super) fn read_manifest(
    path: &Path,
    options: BgzfIndexOptions,
) -> Result<Manifest, BgzfIndexError> {
    let (mut file, state) = open_regular(path)?;
    let limit = as_u64(options.index.max_index_bytes)?;
    check_limit("index_bytes", limit, state.len)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(usize::try_from(state.len).map_err(|_| IndexError::CountOverflow)?)
        .map_err(allocation_error)?;
    (&mut file).take(limit).read_to_end(&mut bytes)?;
    let mut extra = [0];
    loop {
        match file.read(&mut extra) {
            Ok(0) => break,
            Ok(_) => {
                return Err(IndexError::LimitExceeded {
                    resource: "index_bytes",
                    limit,
                    observed: limit.checked_add(1).ok_or(IndexError::CountOverflow)?,
                }
                .into());
            }
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error.into()),
        }
    }
    verify_source(&file, path, &state)
        .map_err(|error| invalid(format!("sidecar changed while reading: {error}")))?;
    serde_json::from_slice::<ObjectOnly<Manifest>>(&bytes)
        .map(|object| object.0)
        .map_err(json_error)
}

#[derive(Default)]
struct Counter(u64);
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(u64::try_from(bytes.len()).map_err(io::Error::other)?)
            .ok_or_else(|| io::Error::other("index size overflow"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
pub(super) fn serialized_size(value: &impl Serialize) -> Result<u64, BgzfIndexError> {
    let mut counter = Counter::default();
    serde_json::to_writer(&mut counter, value).map_err(json_error)?;
    Ok(counter.0)
}
pub(super) fn json_error(error: serde_json::Error) -> BgzfIndexError {
    if let Some(kind) = error.io_error_kind() {
        io::Error::new(kind, error).into()
    } else {
        invalid(error.to_string())
    }
}
