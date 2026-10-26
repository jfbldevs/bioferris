# MUSCLE 5.1 black-box observations

`nature-default-v1.json` is inert data produced in the isolated sibling laboratory
`../rust-sequence-analysis-reference-validation/muscle5`. Cargo never invokes Python,
WSL, or MUSCLE. The official executable was used only through its documented CLI;
its source was not downloaded or inspected.

- Profile: `MUSCLE-5.1.0-Nature-2022-default-v1`
- Executable identity: `muscle 5.1.linux64 [12f0e2]`
- Executable SHA-256: `2bd82928a5cc72b5086b60928a3a820a0ce066ca9fbc0fecd1543eebca981e76`
- Generator SHA-256: `1124b017be24cb63405a3a088abea46fd7254bf904a14528227f6c4b04fd8277`
- Case specification SHA-256: `b2e8cda7286208191c11c30021d2c544a75639bb0520e479a25d45e90849f400`
- Fixture SHA-256: `dba37876d62cde08faa5a7d04052c3a92f468ace0e2da40caa45befb1938e2e8`

The 15 cases cover PPP `Align`, `Super5`, all guide-tree permutations, perturbed
seeds, protein, DNA, RNA, terminal/internal gaps, repeats, redundancy, and divergent
families. Tests compare rows by identifier and report exact aligned strings plus SP
and TC; the external alignment, never the Rust implementation, is the metric oracle.
