use align_core::{AlignmentOperation, PairwiseAlignment};
use align_msa::MultipleAlignment;

/// Exact statistics derived from one canonical pairwise traceback.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PairwiseStatistics {
    columns: usize,
    matches: usize,
    substitutions: usize,
    insertion_residues: usize,
    insertion_runs: usize,
    deletion_residues: usize,
    deletion_runs: usize,
}

impl PairwiseStatistics {
    pub(crate) fn from_alignment(alignment: &PairwiseAlignment) -> Self {
        let mut statistics = Self {
            columns: alignment.columns(),
            matches: 0,
            substitutions: 0,
            insertion_residues: 0,
            insertion_runs: 0,
            deletion_residues: 0,
            deletion_runs: 0,
        };
        let mut previous = None;
        for operation in alignment.operations().iter().copied() {
            match operation {
                AlignmentOperation::Match => statistics.matches += 1,
                AlignmentOperation::Substitution => statistics.substitutions += 1,
                AlignmentOperation::Insertion => {
                    statistics.insertion_residues += 1;
                    if previous != Some(AlignmentOperation::Insertion) {
                        statistics.insertion_runs += 1;
                    }
                }
                AlignmentOperation::Deletion => {
                    statistics.deletion_residues += 1;
                    if previous != Some(AlignmentOperation::Deletion) {
                        statistics.deletion_runs += 1;
                    }
                }
            }
            previous = Some(operation);
        }
        statistics
    }

    /// Returns the number of alignment columns.
    #[must_use]
    pub const fn columns(self) -> usize {
        self.columns
    }

    /// Returns exact-match columns.
    #[must_use]
    pub const fn matches(self) -> usize {
        self.matches
    }

    /// Returns literal-substitution columns.
    #[must_use]
    pub const fn substitutions(self) -> usize {
        self.substitutions
    }

    /// Returns residues occurring only in the second sequence.
    #[must_use]
    pub const fn insertion_residues(self) -> usize {
        self.insertion_residues
    }

    /// Returns contiguous insertion runs.
    #[must_use]
    pub const fn insertion_runs(self) -> usize {
        self.insertion_runs
    }

    /// Returns residues occurring only in the first sequence.
    #[must_use]
    pub const fn deletion_residues(self) -> usize {
        self.deletion_residues
    }

    /// Returns contiguous deletion runs.
    #[must_use]
    pub const fn deletion_runs(self) -> usize {
        self.deletion_runs
    }

    /// Returns exact matches divided by every alignment column.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn identity(self) -> Option<f64> {
        (self.columns != 0).then(|| self.matches as f64 / self.columns as f64)
    }
}

/// Exact structural statistics for a rectangular multiple alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MultipleStatistics {
    sequences: usize,
    columns: usize,
    conserved_columns: usize,
    gap_characters: usize,
}

impl MultipleStatistics {
    pub(crate) fn from_alignment<A: align_core::Alphabet>(
        alignment: &MultipleAlignment<A>,
    ) -> Self {
        let mut conserved_columns = 0;
        let mut gap_characters = 0;
        for column in 0..alignment.width() {
            let mut selected = None;
            let mut conserved = true;
            for row in alignment.rows() {
                if let Some(symbol) = row.symbols()[column] {
                    if let Some(expected) = selected {
                        conserved &= symbol == expected;
                    } else {
                        selected = Some(symbol);
                    }
                } else {
                    gap_characters += 1;
                    conserved = false;
                }
            }
            conserved_columns += usize::from(conserved && selected.is_some());
        }
        Self {
            sequences: alignment.sequence_count(),
            columns: alignment.width(),
            conserved_columns,
            gap_characters,
        }
    }

    /// Returns the number of rows.
    #[must_use]
    pub const fn sequences(self) -> usize {
        self.sequences
    }

    /// Returns the rectangular alignment width.
    #[must_use]
    pub const fn columns(self) -> usize {
        self.columns
    }

    /// Returns literal all-row conserved columns without gaps.
    #[must_use]
    pub const fn conserved_columns(self) -> usize {
        self.conserved_columns
    }

    /// Returns the total number of gap cells.
    #[must_use]
    pub const fn gap_characters(self) -> usize {
        self.gap_characters
    }
}
