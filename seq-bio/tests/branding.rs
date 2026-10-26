//! Public identity of the Cargo package.

#[test]
fn cargo_package_exposes_the_bioferris_identity() {
    assert_eq!(env!("CARGO_PKG_NAME"), "bioferris");
    assert_eq!(
        env!("CARGO_PKG_DESCRIPTION"),
        "Rust-native bioinformatics, from sequences to workflows."
    );
}
