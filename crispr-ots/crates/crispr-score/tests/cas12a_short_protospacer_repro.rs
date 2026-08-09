//! Repro for a Cas12a CFD scoring defect observed with a 20-nt protospacer.
//!
//! Found while scoring the Liu et al. 2019 mini-human AsCpf1 library (20-nt
//! spacers) against an hg38 TTTN/20-nt index. 107,718 of 7.8M reported
//! off-targets came back with CFD exactly 1.000000, including sites the scanner
//! called mm=3 and mm=4. For a 20-nt guide that is arithmetically impossible:
//! every off-diagonal entry of the penalty table is populated from the CSV, and
//! the largest single-mismatch penalty at positions 1..20 of the 2xNLS matrix is
//! 0.8936, so three or four of them cannot multiply to 1.0.
//!
//! The two real cases are hg38 sites reconstructed from the reference and
//! verified base-by-base by hand.

use crispr_encoding::{Base, Site};
use crispr_score::{Cas12aCfd, Cas12aMatrix};

const PROTO_LEN: u32 = 20;
/// kmers_csv encodes non-ACGT PAM bases as `A`; PAM bits are not scored.
const PAM: &str = "TTTA";

fn encode(seq: &str) -> Site {
    let bases: Vec<Base> = seq
        .bytes()
        .map(|c| match c {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            other => panic!("bad base {}", other as char),
        })
        .collect();
    Site::encode_slice(&bases)
}

fn score(cfd: &Cas12aCfd, guide_proto: &str, ot_proto: &str) -> f64 {
    cfd.score_pair_with_len(
        encode(&format!("{PAM}{guide_proto}")),
        encode(&format!("{PAM}{ot_proto}")),
        PROTO_LEN,
    )
}

fn hand_mismatches(a: &str, b: &str) -> Vec<usize> {
    a.bytes()
        .zip(b.bytes())
        .enumerate()
        .filter(|(_, (g, o))| g != o)
        .map(|(i, _)| i + 1)
        .collect()
}

fn matrix() -> Cas12aCfd {
    Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).expect("bundled 2xNLS matrix parses")
}

/// Real hg38 off-targets that crispr-ots scored at CFD = 1.000000.
/// A guide/off-target pair with any mismatch inside the protospacer must score
/// strictly below 1.0, whatever the exact penalties are.
#[test]
fn real_offtargets_with_mismatches_do_not_score_one() {
    let cfd = matrix();
    let guide = "AAACAAACACAATCAGAAGC";
    let cases = [
        ("chr2:24579095 (+), scanner mm=4", "AAACATACACATACAGAACC", 4),
        ("chr11:30221314 (-), scanner mm=3", "AAACAAACAGTAACAGAAGC", 3),
    ];

    let mut failures = Vec::new();
    for (name, ot, scanner_mm) in cases {
        let positions = hand_mismatches(guide, ot);
        assert_eq!(
            positions.len(),
            scanner_mm,
            "{name}: hand-counted mismatches disagree with the scanner"
        );
        let got = score(&cfd, guide, ot);
        eprintln!("{name}: mismatches at {positions:?} -> CFD {got:.6}");
        if got >= 1.0 {
            failures.push(format!("{name}: {} mismatches scored {got:.6}", positions.len()));
        }
    }
    assert!(
        failures.is_empty(),
        "off-targets with real mismatches scored 1.0:\n  {}",
        failures.join("\n  ")
    );
}

/// Look up the raw matrix value for a mismatch, applying the complementing key
/// convention: the CSV's `d` base is on the strand the crRNA pairs with, so the
/// key for an off-target protospacer base `ot` is `d<complement(ot)>`.
fn matrix_value(guide: u8, ot: u8, pos_1indexed: usize) -> f64 {
    let rna = if guide == b'T' { 'U' } else { guide as char };
    let d = match ot {
        b'A' => 'T',
        b'C' => 'G',
        b'G' => 'C',
        b'T' => 'A',
        _ => panic!("bad base"),
    };
    let key = format!("r{rna}:d{d}");
    for line in Cas12aMatrix::TwoXNls.csv().lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() == 4 && f[1] == pos_1indexed.to_string() && f[2] == key {
            return f[3].parse().unwrap();
        }
    }
    panic!("no matrix entry for {key} at position {pos_1indexed}");
}

/// Every protospacer position of a 20-nt guide must be scored with the correct
/// matrix entry. Comparing against the matrix rather than asserting "< 1.0" is
/// deliberate: 1.0 is a legal measured value (2xNLS pos 21 `rA:dC`, pos 23
/// `rC:dC`), so a bare inequality would both miss real defects and flag
/// legitimate entries.
#[test]
fn every_position_of_a_20nt_guide_is_scored() {
    let cfd = matrix();
    let proto = "ACACACACACACACACACAC";

    let mut wrong = Vec::new();
    for i in 0..PROTO_LEN as usize {
        let mut mutated: Vec<u8> = proto.bytes().collect();
        mutated[i] = if mutated[i] == b'A' { b'G' } else { b'T' };
        let got = score(&cfd, proto, &String::from_utf8(mutated.clone()).unwrap());
        let want = matrix_value(proto.as_bytes()[i], mutated[i], i + 1);
        if (got - want).abs() > 1e-12 {
            wrong.push(format!("pos {}: expected {want:.6}, got {got:.6}", i + 1));
        }
    }
    assert!(
        wrong.is_empty(),
        "20-nt CFD disagrees with the matrix:\n  {}",
        wrong.join("\n  ")
    );
}

/// Same check at the matrix-native 23-nt length.
#[test]
fn every_position_of_a_23nt_guide_is_scored() {
    let cfd = matrix();
    let proto = "ACACACACACACACACACACACA";

    let mut wrong = Vec::new();
    for i in 0..23usize {
        let mut mutated: Vec<u8> = proto.bytes().collect();
        mutated[i] = if mutated[i] == b'A' { b'G' } else { b'T' };
        let got = cfd.score_pair_with_len(
            encode(&format!("{PAM}{proto}")),
            encode(&format!("{PAM}{}", String::from_utf8(mutated.clone()).unwrap())),
            23,
        );
        let want = matrix_value(proto.as_bytes()[i], mutated[i], i + 1);
        if (got - want).abs() > 1e-12 {
            wrong.push(format!("pos {}: expected {want:.6}, got {got:.6}", i + 1));
        }
    }
    assert!(
        wrong.is_empty(),
        "23-nt CFD disagrees with the matrix:\n  {}",
        wrong.join("\n  ")
    );
}
