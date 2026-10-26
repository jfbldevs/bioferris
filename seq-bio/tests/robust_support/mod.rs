//! Deterministic, finite adversarial transports and independent test framing.

use std::cell::Cell;
use std::io::{self, BufReader, Read};
use std::rc::Rc;

use bioferris::io::{BioIoError, IoErrorKind};

pub const CHUNKS: [usize; 7] = [1, 2, 3, 5, 8, 13, 64];
pub const MUTATIONS: usize = 32;
pub const RANDOM_CASES: usize = 48;
pub const MAX_RANDOM_BYTES: usize = 96;
pub const MAX_INPUT_BYTES: usize = 8192;

/// Xorshift32 with an explicit, nonzero seed; no ambient entropy or clock.
pub struct Generator(u32);

impl Generator {
    pub const fn new(seed: u32) -> Self {
        assert!(seed != 0);
        Self(seed)
    }

    fn word(&mut self) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 17;
        self.0 ^= self.0 << 5;
        self.0
    }

    pub fn byte(&mut self) -> u8 {
        self.word().to_le_bytes()[0]
    }

    pub fn index(&mut self, bound: usize) -> usize {
        usize::try_from(self.word()).unwrap() % bound
    }

    pub fn dna(&mut self, length: usize) -> String {
        (0..length)
            .map(|_| char::from(b"ACGTacgtNRY"[self.index(11)]))
            .collect()
    }

    pub fn bytes(&mut self, length: usize) -> Vec<u8> {
        (0..length).map(|_| self.byte()).collect()
    }
}

/// Full short-input corpus, not an unbounded fuzz campaign.
pub fn random_corpus() -> Vec<Vec<u8>> {
    let mut generator = Generator::new(0x6d2b_79f5);
    (0..RANDOM_CASES)
        .map(|index| generator.bytes(index * 17 % (MAX_RANDOM_BYTES + 1)))
        .collect()
}

/// Includes the original, four truncations and 32 seeded single-site edits.
pub fn mutations(input: &[u8], seed: u32) -> Vec<Vec<u8>> {
    assert!(!input.is_empty());
    assert!(input.len() < MAX_INPUT_BYTES);
    let mut corpus = vec![input.to_vec()];
    for length in [0, 1, input.len() / 2, input.len() - 1] {
        corpus.push(input[..length].to_vec());
    }
    let mut generator = Generator::new(seed);
    for index in 0..MUTATIONS {
        let mut changed = input.to_vec();
        let position = generator.index(input.len());
        match index % 4 {
            0 => {
                changed.remove(position);
            }
            1 => changed.insert(position, generator.byte()),
            2 => changed[position] ^= 1 << (generator.byte() % 8),
            _ => {
                changed[position] =
                    [0, b'\n', b'\r', b'\t', 0xff, b'"', b'<', b'>'][generator.index(8)];
            }
        }
        corpus.push(changed);
    }
    corpus
}

#[derive(Clone, Copy, Debug)]
pub struct Schedule {
    pub chunk: usize,
    pub interrupt_at: Option<usize>,
    pub fail_at: Option<usize>,
}

impl Schedule {
    pub const fn chunks(chunk: usize) -> Self {
        Self {
            chunk,
            interrupt_at: None,
            fail_at: None,
        }
    }

    pub const fn interrupted(position: usize) -> Self {
        Self {
            chunk: 5,
            interrupt_at: Some(position),
            fail_at: None,
        }
    }

    pub const fn failing(position: usize) -> Self {
        Self {
            chunk: 3,
            interrupt_at: None,
            fail_at: Some(position),
        }
    }
}

/// A legal short reader, with one optional transient interruption or hard error.
pub struct AdversarialRead {
    bytes: Vec<u8>,
    position: usize,
    schedule: Schedule,
    interrupted: bool,
    pub calls: Rc<Cell<usize>>,
}

impl AdversarialRead {
    pub fn new(bytes: &[u8], schedule: Schedule) -> Self {
        assert!(schedule.chunk > 0);
        assert!(bytes.len() <= MAX_INPUT_BYTES);
        Self {
            bytes: bytes.to_vec(),
            position: 0,
            schedule,
            interrupted: false,
            calls: Rc::new(Cell::new(0)),
        }
    }
}

impl Read for AdversarialRead {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if destination.is_empty() {
            return Ok(0);
        }
        self.calls.set(self.calls.get() + 1);
        if self.schedule.fail_at == Some(self.position) {
            return Err(io::Error::other("injected finite-transport failure"));
        }
        if !self.interrupted && self.schedule.interrupt_at == Some(self.position) {
            self.interrupted = true;
            return Err(io::ErrorKind::Interrupted.into());
        }
        let mut end = self
            .position
            .saturating_add(destination.len().min(self.schedule.chunk))
            .min(self.bytes.len());
        for boundary in [
            self.schedule.fail_at,
            self.schedule.interrupt_at.filter(|_| !self.interrupted),
        ]
        .into_iter()
        .flatten()
        {
            if boundary > self.position {
                end = end.min(boundary);
            }
        }
        let count = end - self.position;
        destination[..count].copy_from_slice(&self.bytes[self.position..end]);
        self.position = end;
        Ok(count)
    }
}

pub fn buffered(bytes: &[u8], schedule: Schedule) -> (BufReader<AdversarialRead>, Rc<Cell<usize>>) {
    let source = AdversarialRead::new(bytes, schedule);
    let calls = source.calls.clone();
    (BufReader::with_capacity(7, source), calls)
}

#[derive(Debug, PartialEq, Eq)]
pub struct Failure {
    pub kind: IoErrorKind,
    format: &'static str,
    record: Option<usize>,
    line: Option<usize>,
}

impl From<BioIoError> for Failure {
    fn from(error: BioIoError) -> Self {
        Self {
            kind: error.kind(),
            format: error.format(),
            record: error.record(),
            line: error.line(),
        }
    }
}

pub fn outcome<T>(result: Result<T, BioIoError>) -> Result<T, Failure> {
    result.map_err(Failure::from)
}

#[derive(Debug, PartialEq)]
pub struct StreamOutcome<T> {
    pub values: Vec<T>,
    pub failure: Option<Failure>,
}

/// Bound item production and check both logical fusion and no post-error reads.
pub fn drain<T>(
    mut iterator: impl Iterator<Item = Result<T, BioIoError>>,
    input_bytes: usize,
    calls: &Cell<usize>,
) -> StreamOutcome<T> {
    let mut result = StreamOutcome {
        values: Vec::new(),
        failure: None,
    };
    for _ in 0..=input_bytes + 1 {
        match iterator.next() {
            Some(Ok(value)) => result.values.push(value),
            Some(Err(error)) => {
                result.failure = Some(error.into());
                break;
            }
            None => break,
        }
    }
    assert!(
        result.values.len() <= input_bytes + 1,
        "item budget exceeded"
    );
    let before = calls.get();
    assert!(iterator.next().is_none(), "iterator did not fuse");
    assert!(
        iterator.next().is_none(),
        "iterator resumed after termination"
    );
    assert_eq!(before, calls.get(), "fused iterator touched its source");
    result
}

/// Bitwise RFC 1952 CRC32, deliberately not the product's CRC implementation.
pub fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0_u32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ if crc & 1 == 1 { 0xedb8_8320 } else { 0 };
        }
    }
    !crc
}

fn stored_deflate(input: &[u8]) -> Vec<u8> {
    let length = u16::try_from(input.len()).unwrap();
    let mut bytes = vec![1]; // RFC 1951 BFINAL=1, BTYPE=00 with byte padding.
    bytes.extend_from_slice(&length.to_le_bytes());
    bytes.extend_from_slice(&(!length).to_le_bytes());
    bytes.extend_from_slice(input);
    bytes
}

pub fn gzip_member(input: &[u8]) -> Vec<u8> {
    let mut bytes = vec![31, 139, 8, 0, 0, 0, 0, 0, 0, 255];
    bytes.extend_from_slice(&stored_deflate(input));
    bytes.extend_from_slice(&crc32(input).to_le_bytes());
    bytes.extend_from_slice(&u32::try_from(input.len()).unwrap().to_le_bytes());
    bytes
}

pub fn bgzf_member(input: &[u8]) -> Vec<u8> {
    let payload = stored_deflate(input);
    let size = 18 + payload.len() + 8;
    let mut bytes = vec![31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0];
    bytes.extend_from_slice(&u16::try_from(size - 1).unwrap().to_le_bytes());
    bytes.extend_from_slice(&payload);
    bytes.extend_from_slice(&crc32(input).to_le_bytes());
    bytes.extend_from_slice(&u32::try_from(input.len()).unwrap().to_le_bytes());
    bytes
}
