//! Common inspection helpers across algorithm-specific result types.

use std::path::Path;

use align_core::{AlignmentOperation, Alphabet, PairwiseAlignment};
use align_fast::{WavefrontAlignment, XDropAlignment, XDropSeedAlignment};
use align_msa::MultipleAlignment;

use crate::pairwise::AdaptiveAffineAlignment;

/// Borrows the canonical operation-based pairwise alignment from a result.
pub trait PairwiseAlignmentView {
    /// Returns the canonical scored traceback.
    fn alignment(&self) -> &PairwiseAlignment;
}

impl PairwiseAlignmentView for PairwiseAlignment {
    fn alignment(&self) -> &PairwiseAlignment {
        self
    }
}

impl PairwiseAlignmentView for WavefrontAlignment {
    fn alignment(&self) -> &PairwiseAlignment {
        self.alignment()
    }
}

impl PairwiseAlignmentView for AdaptiveAffineAlignment {
    fn alignment(&self) -> &PairwiseAlignment {
        self.alignment()
    }
}

impl PairwiseAlignmentView for XDropAlignment {
    fn alignment(&self) -> &PairwiseAlignment {
        self.alignment()
    }
}

impl PairwiseAlignmentView for XDropSeedAlignment {
    fn alignment(&self) -> &PairwiseAlignment {
        self.alignment()
    }
}

/// Convenient statistics and extended CIGAR output for any pairwise result.
pub trait PairwiseAlignmentExt: PairwiseAlignmentView {
    /// Returns an extended CIGAR string using `=`, `X`, `D`, and `I` operations.
    #[must_use]
    fn cigar(&self) -> String {
        let operations = self.alignment().operations();
        let mut output = String::with_capacity(operations.len().saturating_mul(2));
        let Some(first) = operations.first().copied() else {
            return output;
        };
        let mut selected = first;
        let mut count = 0_usize;
        for operation in operations.iter().copied() {
            if operation == selected {
                count += 1;
            } else {
                append_cigar_run(&mut output, count, selected);
                selected = operation;
                count = 1;
            }
        }
        append_cigar_run(&mut output, count, selected);
        output
    }

    /// Returns the number of exact-match columns.
    #[must_use]
    fn matches(&self) -> usize {
        self.alignment()
            .operations()
            .iter()
            .filter(|operation| matches!(operation, AlignmentOperation::Match))
            .count()
    }

    /// Returns exact matches divided by alignment columns, or `None` when empty.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    fn identity(&self) -> Option<f64> {
        let columns = self.alignment().columns();
        (columns != 0).then(|| self.matches() as f64 / columns as f64)
    }
}

impl<T: PairwiseAlignmentView + ?Sized> PairwiseAlignmentExt for T {}

fn append_cigar_run(output: &mut String, count: usize, operation: AlignmentOperation) {
    use core::fmt::Write as _;

    let code = match operation {
        AlignmentOperation::Match => '=',
        AlignmentOperation::Substitution => 'X',
        AlignmentOperation::Deletion => 'D',
        AlignmentOperation::Insertion => 'I',
    };
    write!(output, "{count}{code}").expect("writing to a String cannot fail");
}

/// Borrows the rectangular alignment from an algorithm-specific MSA result.
pub trait MultipleAlignmentView {
    /// Alphabet shared by every alignment row.
    type Alphabet: Alphabet;

    /// Returns the canonical rectangular alignment.
    fn alignment(&self) -> &MultipleAlignment<Self::Alphabet>;
}

/// Common output methods for algorithm-specific MSA results.
pub trait MultipleAlignmentExt: MultipleAlignmentView {
    /// Writes the canonical alignment as aligned FASTA.
    ///
    /// # Errors
    ///
    /// Returns [`align_formats::FastaError`] when the file cannot be created or
    /// written.
    fn write_fasta(&self, path: impl AsRef<Path>) -> Result<(), align_formats::FastaError> {
        align_formats::write_alignment_fasta(path, self.alignment())
    }
}

impl<T: MultipleAlignmentView + ?Sized> MultipleAlignmentExt for T {}
