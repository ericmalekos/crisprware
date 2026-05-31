use crispr_encoding::{Base, Site};
use crispr_score::{Cas12aCfd, Cas12aMatrix};

fn encode(s: &str) -> Site {
    let bases: Vec<Base> = s.bytes().map(|c| match c.to_ascii_uppercase() {
        b'A' => Base::A, b'C' => Base::C, b'G' => Base::G, b'T' => Base::T,
        _ => panic!("non-ACGT"),
    }).collect();
    Site::encode(bases)
}

#[test]
fn validate_against_parasol() {
    let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
    #[allow(clippy::unreadable_literal)] // these match the parasol output verbatim
    let cases = [
        ("perfect", "TTTAAAAAAAAAAAAAAAAAAAAA", "TTTAAAAAAAAAAAAAAAAAAAAA", 1.0),
        ("1mm",     "TTTAAAAAAAAAAAAAAAAAAAAA", "TTTAACAAAAAAAAAAAAAAAAAA", 0.1029411765),
        ("2mm",     "TTTAAAAAAAAAAAAAAAAAAAAA", "TTTACCAAAAAAAAAAAAAAAAAA", 0.0209099265),
        ("3mm",     "TTTACGTACGTACGTACGTACGTA", "TTTACGAACGTACGTACGTACTTA", 0.8505747126),
    ];
    for (label, t, o, exp) in cases {
        let s = scorer.score_pair(encode(t), encode(o));
        let d = (s - exp).abs();
        println!("  {label}: ours = {s:.10}, parasol = {exp:.10}, |Δ| = {d:.3e}");
        assert!(d < 1e-6, "mismatch on {label}: {s} vs {exp}");
    }
}
