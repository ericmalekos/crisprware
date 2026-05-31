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
        // 27-mer Cas12a sites: 4-bp TTTN PAM + 23-bp protospacer.
        ("perfect",    "TTTAAAAAAAAAAAAAAAAAAAAAAAA", "TTTAAAAAAAAAAAAAAAAAAAAAAAA", 1.0),
        ("1mm pos2",   "TTTAAAAAAAAAAAAAAAAAAAAAAAA", "TTTAACAAAAAAAAAAAAAAAAAAAAA", 0.10294117647058823),
        ("2mm pos1,2", "TTTAAAAAAAAAAAAAAAAAAAAAAAA", "TTTACCAAAAAAAAAAAAAAAAAAAAA", 0.020909926470588234),
        // 4 mismatches in PAM-distal positions 18-21 (1-indexed in the
        // protospacer). Exercises matrix positions 21..23 that the
        // legacy 20-nt enzyme couldn't reach.
        ("distal 4mm", "TTTAAAAAAAAAAAAAAAAAAAAAAAA", "TTTAAAAAAAAAAAAAAAAAACCCCAA", 0.10851063829787233),
        ("realistic",  "TTTACGTACGTACGTACGTACGTACGT", "TTTACGAACGTACGTACGTACTTACGT", 0.8505747126436781),
    ];
    for (label, t, o, exp) in cases {
        let s = scorer.score_pair(encode(t), encode(o));
        let d = (s - exp).abs();
        println!("  {label}: ours = {s:.10}, parasol = {exp:.10}, |Δ| = {d:.3e}");
        assert!(d < 1e-6, "mismatch on {label}: {s} vs {exp}");
    }
}
