use align_core::{Alphabet, Sequence};

use crate::MsaError;

/// One uniquely nameable, ungapped biological input sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceRecord<A: Alphabet> {
    identifier: String,
    sequence: Sequence<A>,
}

impl<A: Alphabet> SequenceRecord<A> {
    /// Creates a named sequence after validating its identifier.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] when the identifier is empty or contains an ASCII
    /// control character. Collection-level duplicate detection occurs in an
    /// alignment algorithm.
    pub fn new(identifier: impl Into<String>, sequence: Sequence<A>) -> Result<Self, MsaError> {
        let identifier = identifier.into();
        validate_identifier(&identifier)?;
        Ok(Self {
            identifier,
            sequence,
        })
    }

    /// Returns the stable record identifier.
    #[must_use]
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    /// Returns the validated ungapped sequence.
    #[must_use]
    pub const fn sequence(&self) -> &Sequence<A> {
        &self.sequence
    }

    /// Separates the identifier and sequence.
    #[must_use]
    pub fn into_parts(self) -> (String, Sequence<A>) {
        (self.identifier, self.sequence)
    }
}

pub(crate) fn validate_identifier(identifier: &str) -> Result<(), MsaError> {
    if identifier.trim().is_empty() {
        return Err(MsaError::EmptyIdentifier);
    }
    if let Some((index, _)) = identifier
        .as_bytes()
        .iter()
        .enumerate()
        .find(|(_, byte)| byte.is_ascii_control())
    {
        return Err(MsaError::ControlCharacterInIdentifier {
            identifier: identifier.to_owned(),
            index,
        });
    }
    Ok(())
}
