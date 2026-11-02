//! Regression on real chewBBACA hashed calls (Listeria, Pasteur cgMLST, 1748
//! loci): CRC32 hashes above i32::MAX used to be read as missing, collapsing
//! the distances (282 and 11 instead of 1196 and 35). Expected values come
//! from an independent per-locus comparison: `*N` novel alleles missing by
//! default, compared as alleles with --novel-as-allele.
use std::process::Command;

#[test]
fn listeria_hashed_profiles() {
    let out = Command::new(env!("CARGO_BIN_EXE_cgmlst-dists"))
        .args(["-s", "-i", concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/listeria_hashed_3.tsv")])
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows[0], ["cgmlst-dists", "ERR10109288", "ERR10109289", "ERR10109291"]);
    assert_eq!(rows[1][1..], ["0", "1196", "1199"]);
    assert_eq!(rows[2][1..], ["1196", "0", "35"]);
    assert_eq!(rows[3][1..], ["1199", "35", "0"]);
}

#[test]
fn listeria_hashed_profiles_novel_as_allele() {
    let out = Command::new(env!("CARGO_BIN_EXE_cgmlst-dists"))
        .args(["-s", "-N", "-i", concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/listeria_hashed_3.tsv")])
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows[1][1..], ["0", "1234", "1233"]);
    assert_eq!(rows[2][1..], ["1234", "0", "59"]);
    assert_eq!(rows[3][1..], ["1233", "59", "0"]);
}
