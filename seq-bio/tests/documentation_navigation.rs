//! Local link/anchor regression coverage for the toolkit's maintained guides.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

const GUIDES: &[&str] = &[
    "docs/testing/security.md",
    "docs/testing/README.md",
    "docs/testing/robustness.md",
    "docs/testing/contracts-and-science.md",
    "docs/testing/faults-and-resources.md",
    "docs/testing/coverage.md",
    "docs/testing/mutation.md",
    "docs/testing/portability.md",
    "docs/testing/performance.md",
    "docs/algorithms/linear-score-only.md",
    "docs/algorithms/binary-sequence-index-v1.md",
    "docs/algorithms/muscle5.md",
    "docs/examples/alignment-recipes.md",
    "seq-bio/drafts/README.md",
    "README.md",
    "docs/README.md",
    "docs/quick-start.md",
    "docs/algorithms/cds-translation-v1.md",
    "docs/api/cds-translation-v1.md",
    "seq-bio/tests/fixtures/cds-provenance.md",
    "docs/algorithms/protein-properties-v1.md",
    "docs/api/protein-properties-v1.md",
    "seq-bio/tests/fixtures/protein-properties-provenance.md",
    "docs/algorithms/motifs-v1.md",
    "docs/algorithms/alignment-analysis-v1.md",
    "docs/algorithms/blast-input-v1.md",
    "docs/algorithms/compressed-output-v1.md",
    "docs/algorithms/bgzf-v1.md",
    "docs/algorithms/bgzf-index-v1.md",
    "docs/api/motifs-v1.md",
    "docs/api/alignment-analysis-v1.md",
    "docs/api/blast-input-v1.md",
    "docs/api/compressed-io-v1.md",
    "seq-bio/tests/fixtures/motifs-provenance.md",
    "seq-bio/tests/fixtures/blast-provenance.md",
    "seq-bio/tests/fixtures/bgzf-provenance.md",
    "docs/examples/README.md",
    "docs/examples/sequence-index-recipes.md",
    "docs/examples/annotated-output-recipes.md",
    "docs/api/seq-bio-v1.md",
    "docs/api/sequence-operations-v1.md",
    "docs/api/sequence-statistics-v1.md",
    "docs/api/genetic-codes-v1.md",
    "docs/api/reading-frames-v1.md",
    "docs/algorithms/reading-frames-v1.md",
    "docs/api/sequence-index-v1.md",
    "docs/api/bioferris-identity.md",
    "docs/api/rust-api-v1.md",
    "docs/api/phylogenetic-output-v1.md",
    "docs/api/file-io-v1.md",
    "docs/api/alignment-input-v1.md",
    "docs/api/annotated-input-v1.md",
    "docs/api/annotated-output-v1.md",
    "docs/api/graph-result-input-v1.md",
    "seq-bio/tests/fixtures/README.md",
    "seq-bio/tests/fixtures/annotated-provenance.md",
];

/// CLI test files deleted with the CLI (decision 029) that historical validation
/// records still cite. Those records stay unchanged, so only links from
/// `docs/validation/` to exactly these files may be missing.
const REMOVED_TARGETS: &[&str] = &[];

fn prose_lines(text: &str) -> impl Iterator<Item = &str> {
    let mut fence = false;
    text.lines().filter(move |line| {
        if line.trim_start().starts_with("```") {
            fence = !fence;
            return false;
        }
        !fence
    })
}

fn anchors(text: &str) -> HashSet<String> {
    let mut counts = HashMap::<String, usize>::new();
    let mut output = HashSet::new();
    for line in prose_lines(text) {
        let line = line.trim_start();
        let Some(heading) = line.strip_prefix('#') else {
            continue;
        };
        let heading = heading.trim_start_matches('#').trim();
        let base: String = heading
            .chars()
            .flat_map(char::to_lowercase)
            .filter_map(|character| {
                if character.is_whitespace() {
                    Some('-')
                } else if character.is_alphanumeric() || matches!(character, '-' | '_') {
                    Some(character)
                } else {
                    None
                }
            })
            .collect();
        let count = counts.entry(base.clone()).or_default();
        output.insert(if *count == 0 {
            base.clone()
        } else {
            format!("{base}-{count}")
        });
        *count += 1;
    }
    output
}

fn check_guide(
    guide: &Path,
    targets: &mut HashMap<PathBuf, HashSet<String>>,
    links: &mut usize,
    fragments: &mut usize,
    removed: &mut HashSet<&'static str>,
) {
    let historical = guide
        .parent()
        .unwrap()
        .iter()
        .any(|part| part == "validation");
    let content = fs::read_to_string(guide).unwrap();
    for line in prose_lines(&content) {
        let mut rest = line;
        while let Some((_, after)) = rest.split_once("](") {
            let Some((destination, tail)) = after.split_once(')') else {
                break;
            };
            rest = tail;
            if destination.contains("://") || destination.starts_with("mailto:") {
                continue;
            }
            let destination = destination.trim_matches(['<', '>']);
            let (file, fragment) = destination
                .split_once('#')
                .map_or((destination, None), |(file, fragment)| {
                    (file, Some(fragment))
                });
            let path = if file.is_empty() {
                guide.to_path_buf()
            } else {
                guide.parent().unwrap().join(file)
            };
            if !path.exists() && historical {
                if let Some(target) = REMOVED_TARGETS
                    .iter()
                    .find(|target| **target == destination)
                {
                    removed.insert(target);
                    continue;
                }
            }
            assert!(
                path.exists(),
                "{} has missing target {destination}",
                guide.display()
            );
            *links += 1;
            if let Some(fragment) = fragment.filter(|fragment| !fragment.is_empty()) {
                let target = path.canonicalize().unwrap();
                let found = targets
                    .entry(target.clone())
                    .or_insert_with(|| anchors(&fs::read_to_string(&target).unwrap()));
                assert!(
                    found.contains(fragment),
                    "{} has missing anchor {destination}; available: {found:?}",
                    guide.display()
                );
                *fragments += 1;
            }
        }
    }
}

#[test]
fn maintained_guides_have_existing_local_links_and_heading_anchors() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let mut targets = HashMap::new();
    let mut links = 0usize;
    let mut fragments = 0usize;
    let mut removed = HashSet::new();
    for guide in GUIDES {
        check_guide(
            &root.join(guide),
            &mut targets,
            &mut links,
            &mut fragments,
            &mut removed,
        );
    }
    println!(
        "{} guides; {links} local links; {fragments} heading references validated; \
         {} removed CLI tests cited by historical records",
        GUIDES.len(),
        removed.len()
    );
    assert!(links > 200 && fragments > 40);
    // Every allowance must still be needed: none may be stale.
    assert_eq!(removed.len(), REMOVED_TARGETS.len());
}
