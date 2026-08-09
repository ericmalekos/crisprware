//! Golden values for the Cas12a CFD scorer.
//!
//! **These values intentionally NO LONGER match `parasol_scripts/score_flashfry_cfd.py`.**
//! That reference implementation builds its lookup key without complementing the
//! off-target base:
//!
//! ```python
//! mm = f"r{dna_to_rna_base(tb)}:d{ob}"   # ob not complemented
//! score *= lookup.get((pos, mm), 1.0)    # a missing key silently scores 1.0
//! ```
//!
//! The matrices use the complementary-strand convention: the `d` base is on the
//! strand the crRNA pairs with, so the four Watson-Crick pairs (`rA:dT`,
//! `rC:dG`, `rG:dC`, `rU:dA`) are absent because they are matches, while the
//! same-letter pairs are present as genuine mismatches. Under the reference's
//! uncomplemented key, one of every three mismatch types misses the table and,
//! thanks to the `.get(..., 1.0)` default, is scored with **no penalty at all**.
//! See `crispr-score/src/cas12a.rs::parse_mm_key` for the full evidence, and
//! `crispr-score/src/cfd.rs::parse_mm_key`, which already handled this correctly
//! for the SpCas9 matrix.
//!
//! Expected values below were derived directly from the bundled 2xNLS CSV under
//! the corrected convention. The `parasol` column records what the buggy
//! reference produced, to document the magnitude of the correction.

use crispr_encoding::{Base, Site};
use crispr_score::{Cas12aCfd, Cas12aMatrix};

fn encode(s: &str) -> Site {
    let bases: Vec<Base> = s
        .bytes()
        .map(|c| match c.to_ascii_uppercase() {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'T' => Base::T,
            _ => panic!("non-ACGT"),
        })
        .collect();
    Site::encode(bases)
}

#[test]
fn cas12a_cfd_golden_values() {
    let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
    #[allow(clippy::unreadable_literal)]
    // (label, guide site, off-target site, corrected CFD, what parasol produced)
    let cases = [
        // 27-mer Cas12a sites: 4-bp TTTN PAM + 23-bp protospacer.
        (
            "perfect",
            "TTTAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTAAAAAAAAAAAAAAAAAAAAAAAA",
            1.0,
            1.0,
        ),
        // pos 2, guide A vs OT C -> key rA:dG = 0.0800 (parasol used rA:dC).
        (
            "1mm pos2",
            "TTTAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTAACAAAAAAAAAAAAAAAAAAAAA",
            0.08,
            0.10294117647058823,
        ),
        // rA:dG at pos 1 (0.1351) x pos 2 (0.0800).
        (
            "2mm pos1,2",
            "TTTAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTACCAAAAAAAAAAAAAAAAAAAAA",
            0.010810810810810811,
            0.020909926470588234,
        ),
        // 4 mismatches at PAM-distal positions 18-21, exercising matrix
        // positions a 20-nt enzyme cannot reach.
        (
            "distal 4mm",
            "TTTAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTAAAAAAAAAAAAAAAAAACCCCAA",
            0.05624844930982941,
            0.10851063829787233,
        ),
        // The case that shows why this matters. Two mismatches:
        //   pos 3  guide T vs OT A -> rU:dT = 0.1096
        //          parasol keyed rU:dA, which is ABSENT -> scored 1.0, free pass
        //   pos 18 guide G vs OT T -> rG:dA = 0.3333
        //          parasol keyed rG:dT = 0.8506, the most-tolerated entry
        // Net: parasol overestimated this off-target's activity by ~23x,
        // reporting a near-perfect-match 0.851 for a site whose real predicted
        // activity is 0.037.
        (
            "realistic",
            "TTTACGTACGTACGTACGTACGTACGT",
            "TTTACGAACGTACGTACGTACTTACGT",
            0.0365296803652968,
            0.8505747126436781,
        ),
    ];

    for (label, t, o, expected, parasol) in cases {
        let s = scorer.score_pair(encode(t), encode(o));
        let d = (s - expected).abs();
        println!(
            "  {label}: ours = {s:.10}, expected = {expected:.10}, \
             (buggy parasol reference = {parasol:.10})"
        );
        assert!(d < 1e-9, "mismatch on {label}: {s} vs {expected}");
    }
}
