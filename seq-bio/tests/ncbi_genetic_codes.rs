//! Every NCBI genetic code against NCBI's own `gc.prt`, parsed independently.

use bioferris::prelude::*;

/// One table of the `gc.prt` fixture, read without the production tables.
struct NcbiTable {
    id: u16,
    name: String,
    amino_acids: Vec<u8>,
    starts: Vec<u8>,
}

/// The quoted value after `key`, where `key` starts a whitespace-delimited token.
fn quoted<'a>(block: &'a str, key: &str) -> Option<&'a str> {
    let (index, _) = block
        .match_indices(key)
        .find(|(index, _)| block[..*index].ends_with(char::is_whitespace))?;
    let rest = &block[index + key.len()..];
    let open = rest.find('"')?;
    let close = rest[open + 1..].find('"')?;
    Some(&rest[open + 1..open + 1 + close])
}

fn ncbi_tables() -> Vec<NcbiTable> {
    // Comments are whole lines starting with `--`; the strings contain `--`.
    let text = include_str!("fixtures/ncbi-gc.prt")
        .lines()
        .filter(|line| !line.trim_start().starts_with("--"))
        .collect::<Vec<_>>()
        .join("\n");
    text.split('}')
        .filter_map(|block| {
            let (index, _) = block
                .match_indices("id ")
                .find(|(index, _)| block[..*index].ends_with(char::is_whitespace))?;
            let digits = block[index + 3..]
                .chars()
                .skip_while(|c| c.is_whitespace())
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            Some(NcbiTable {
                id: digits.parse().unwrap(),
                name: quoted(block, "name")?
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" "),
                amino_acids: quoted(block, "ncbieaa")?.as_bytes().to_vec(),
                starts: quoted(block, "sncbieaa")?.as_bytes().to_vec(),
            })
        })
        .collect()
}

fn codon(index: usize) -> String {
    let bases = b"TCAG";
    [bases[index / 16], bases[index / 4 % 4], bases[index % 4]]
        .iter()
        .map(|base| char::from(*base))
        .collect()
}

fn translate(sequence: &str, code: GeneticCode, mode: TranslationMode) -> Result<String, String> {
    dna(sequence)
        .unwrap()
        .translate_with(TranslationConfig {
            genetic_code: code,
            mode,
            ..TranslationConfig::default()
        })
        .map(|product| product.to_string())
        .map_err(|error| error.to_string())
}

#[test]
fn every_table_matches_ncbi_gc_prt() {
    let tables = ncbi_tables();
    assert_eq!(
        tables.iter().map(|table| table.id).collect::<Vec<_>>(),
        GeneticCode::ALL.map(GeneticCode::id).to_vec()
    );
    for table in &tables {
        assert_eq!(table.amino_acids.len(), 64);
        assert_eq!(table.starts.len(), 64);
        let code = GeneticCode::try_from(table.id).unwrap();
        assert_eq!(code.name(), table.name, "table {}", table.id);
        let conditional = table
            .starts
            .iter()
            .zip(&table.amino_acids)
            .any(|(start, residue)| *start == b'*' && *residue != b'*');
        assert_eq!(code.has_conditional_stops(), conditional);
        let terminator = codon(table.starts.iter().position(|mark| *mark == b'*').unwrap());
        for index in 0..64 {
            let triplet = codon(index);
            let context = format!("table {} codon {triplet}", table.id);
            // A lone codon translates to gc.prt's amino acid, or `*`.
            assert_eq!(
                translate(&triplet, code, TranslationMode::Full).unwrap(),
                char::from(table.amino_acids[index]).to_string(),
                "{context}"
            );
            // A complete CDS starts only at gc.prt's initiators ...
            let started = translate(
                &format!("{triplet}GCT{terminator}"),
                code,
                TranslationMode::Cds,
            );
            assert_eq!(started.is_ok(), table.starts[index] == b'M', "{context}");
            // ... and ends at every codon gc.prt marks as able to end translation.
            let ended = translate(&format!("ATGGCT{triplet}"), code, TranslationMode::Cds);
            assert_eq!(ended.is_ok(), table.starts[index] == b'*', "{context}");
        }
    }
}

#[test]
fn conditional_stops_read_as_amino_acids_except_at_the_end_of_a_cds() {
    for (code, sequence, full, cds) in [
        (GeneticCode::KaryorelictNuclear, "ATGTGATGA", "MWW", "MW"),
        (GeneticCode::CondylostomaNuclear, "ATGTAGTAA", "MQQ", "MQ"),
        (
            GeneticCode::BlastocrithidiaNuclear,
            "ATGTAATAG",
            "MEE",
            "ME",
        ),
    ] {
        assert_eq!(
            translate(sequence, code, TranslationMode::Full).unwrap(),
            full
        );
        assert_eq!(
            translate(sequence, code, TranslationMode::ToStop).unwrap(),
            full
        );
        assert_eq!(
            translate(sequence, code, TranslationMode::Cds).unwrap(),
            cds
        );
    }
}
