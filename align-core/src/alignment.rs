use core::{fmt, ops::Range};

use crate::{AlignmentStructureError, Alphabet, Score, Sequence};

/// One column-producing operation in a pairwise alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentOperation {
    /// Equal symbols; consumes one symbol from each sequence.
    Match,
    /// Unequal symbols; consumes one symbol from each sequence.
    Substitution,
    /// A symbol from the first sequence aligned to a gap in the second.
    Deletion,
    /// A gap in the first sequence aligned to a symbol from the second.
    Insertion,
}

impl AlignmentOperation {
    /// Returns whether this operation consumes the first sequence.
    #[must_use]
    pub const fn consumes_first(self) -> bool {
        matches!(self, Self::Match | Self::Substitution | Self::Deletion)
    }

    /// Returns whether this operation consumes the second sequence.
    #[must_use]
    pub const fn consumes_second(self) -> bool {
        matches!(self, Self::Match | Self::Substitution | Self::Insertion)
    }
}

/// A scored, operation-based pairwise alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairwiseAlignment {
    score: Score,
    first_range: Range<usize>,
    second_range: Range<usize>,
    operations: Vec<AlignmentOperation>,
}

impl PairwiseAlignment {
    /// Constructs an alignment after checking operation consumption.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentStructureError`] for inverted ranges or when the
    /// operations do not consume exactly the declared ranges.
    pub fn try_new(
        score: Score,
        first_range: Range<usize>,
        second_range: Range<usize>,
        operations: Vec<AlignmentOperation>,
    ) -> Result<Self, AlignmentStructureError> {
        if first_range.start > first_range.end {
            return Err(AlignmentStructureError::InvalidRange {
                sequence: "first",
                start: first_range.start,
                end: first_range.end,
            });
        }
        if second_range.start > second_range.end {
            return Err(AlignmentStructureError::InvalidRange {
                sequence: "second",
                start: second_range.start,
                end: second_range.end,
            });
        }
        let (consumed_first, consumed_second) =
            operations
                .iter()
                .fold((0_usize, 0_usize), |(first, second), operation| {
                    (
                        first + usize::from(operation.consumes_first()),
                        second + usize::from(operation.consumes_second()),
                    )
                });
        let expected_first = first_range.len();
        let expected_second = second_range.len();
        if consumed_first != expected_first {
            return Err(AlignmentStructureError::ConsumptionMismatch {
                sequence: "first",
                expected: expected_first,
                actual: consumed_first,
            });
        }
        if consumed_second != expected_second {
            return Err(AlignmentStructureError::ConsumptionMismatch {
                sequence: "second",
                expected: expected_second,
                actual: consumed_second,
            });
        }
        Ok(Self {
            score,
            first_range,
            second_range,
            operations,
        })
    }

    /// Returns the optimal score.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Returns the consumed range in the first sequence.
    #[must_use]
    pub fn first_range(&self) -> Range<usize> {
        self.first_range.clone()
    }

    /// Returns the consumed range in the second sequence.
    #[must_use]
    pub fn second_range(&self) -> Range<usize> {
        self.second_range.clone()
    }

    /// Returns traceback operations in forward alignment order.
    #[must_use]
    pub fn operations(&self) -> &[AlignmentOperation] {
        &self.operations
    }

    /// Returns the number of alignment columns.
    #[must_use]
    pub fn columns(&self) -> usize {
        self.operations.len()
    }

    /// Renders both aligned sequences with `-` for gaps.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentStructureError`] if a stored range is outside a supplied
    /// sequence or if output memory cannot be reserved.
    pub fn render<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<RenderedAlignment, AlignmentStructureError> {
        if self.first_range.end > first.len() {
            return Err(AlignmentStructureError::RangeOutOfBounds {
                sequence: "first",
                end: self.first_range.end,
                length: first.len(),
            });
        }
        if self.second_range.end > second.len() {
            return Err(AlignmentStructureError::RangeOutOfBounds {
                sequence: "second",
                end: self.second_range.end,
                length: second.len(),
            });
        }

        let mut rendered_first = String::new();
        let mut rendered_second = String::new();
        rendered_first
            .try_reserve_exact(self.operations.len())
            .map_err(|_| AlignmentStructureError::AllocationFailed {
                columns: self.operations.len(),
            })?;
        rendered_second
            .try_reserve_exact(self.operations.len())
            .map_err(|_| AlignmentStructureError::AllocationFailed {
                columns: self.operations.len(),
            })?;

        let mut first_index = self.first_range.start;
        let mut second_index = self.second_range.start;
        for operation in &self.operations {
            if operation.consumes_first() {
                let symbol =
                    first
                        .get(first_index)
                        .ok_or(AlignmentStructureError::RangeOutOfBounds {
                            sequence: "first",
                            end: first_index.saturating_add(1),
                            length: first.len(),
                        })?;
                rendered_first.push(char::from(symbol.as_ascii()));
                first_index += 1;
            } else {
                rendered_first.push('-');
            }

            if operation.consumes_second() {
                let symbol =
                    second
                        .get(second_index)
                        .ok_or(AlignmentStructureError::RangeOutOfBounds {
                            sequence: "second",
                            end: second_index.saturating_add(1),
                            length: second.len(),
                        })?;
                rendered_second.push(char::from(symbol.as_ascii()));
                second_index += 1;
            } else {
                rendered_second.push('-');
            }
        }

        Ok(RenderedAlignment {
            first: rendered_first,
            second: rendered_second,
        })
    }
}

/// Human-readable gapped sequences produced from an operation traceback.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RenderedAlignment {
    first: String,
    second: String,
}

impl RenderedAlignment {
    /// Returns the gapped first sequence.
    #[must_use]
    pub fn first(&self) -> &str {
        &self.first
    }

    /// Returns the gapped second sequence.
    #[must_use]
    pub fn second(&self) -> &str {
        &self.second
    }
}

impl fmt::Display for RenderedAlignment {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}\n{}", self.first, self.second)
    }
}

#[cfg(test)]
mod tests {
    use crate::{AlignmentOperation, DnaSequence, PairwiseAlignment};

    #[test]
    fn constructor_rejects_wrong_consumption() {
        let error = PairwiseAlignment::try_new(0, 0..2, 0..1, vec![AlignmentOperation::Match])
            .expect_err("one match cannot consume two first-sequence symbols");
        assert!(error.to_string().contains("expected 2"));
    }

    #[test]
    fn constructor_rejects_inverted_ranges() {
        let start = 2;
        let end = 1;
        let error = PairwiseAlignment::try_new(0, start..end, 0..0, vec![])
            .expect_err("inverted ranges are invalid");
        assert!(error.to_string().contains("starts at 2"));
    }

    #[test]
    fn operations_render_without_embedding_gaps_in_sequences() {
        let first = DnaSequence::from_ascii(b"AC").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"AG").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(
            -1,
            0..2,
            0..2,
            vec![
                AlignmentOperation::Match,
                AlignmentOperation::Deletion,
                AlignmentOperation::Insertion,
            ],
        )
        .expect("operation consumption is valid");
        let rendered = alignment.render(&first, &second).expect("render is valid");
        assert_eq!(rendered.first(), "AC-");
        assert_eq!(rendered.second(), "A-G");
        assert_eq!(rendered.to_string(), "AC-\nA-G");
    }
}
