//! Human-readable reports derived from native scientific results.

mod graph;
mod multiple;
mod options;
mod pairwise;
mod search;
mod statistics;
mod value;

pub use graph::{PoastaQueryReport, PoastaQueryRun};
pub use multiple::{MultipleReport, MultipleRun, MultipleRunInfo};
pub use options::{
    ColorMode, CoordinateStyle, DEFAULT_BLOCK_WIDTH, ReportError, ReportOptions, Verbosity,
};
pub use pairwise::{PairwiseMetric, PairwiseReport, PairwiseRun, PairwiseRunInfo, ReportField};
pub use search::{MyersDistanceReport, MyersDistanceRun, MyersSearchReport, MyersSearchRun};
pub use statistics::{MultipleStatistics, PairwiseStatistics};
pub use value::{PairwiseValueReport, PairwiseValueRun};

fn sanitized_label(label: impl Into<String>) -> String {
    let label = label.into();
    let mut output = String::with_capacity(label.len());
    for character in label.chars() {
        if character.is_control() {
            output.extend(character.escape_default());
        } else {
            output.push(character);
        }
    }
    output
}
