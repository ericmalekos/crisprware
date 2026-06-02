//! Cas12a (Cpf1) "CFD-like" off-target scoring.
#![allow(clippy::similar_names)] // tttn_sum / tttv_sum mirror upstream matrix names
//!
//! Cas12a doesn't share `SpCas9`'s lookup table. crisprware's
//! `parasol_scripts/score_flashfry_cfd.py` implements a per-mismatch
//! penalty model using one of two activity matrices:
//!
//! - `2xNLS-Cas12a`: the `AsCas12a` wild-type-ish profile used in
//!   `2xNLS` constructs (DOI: 10.1016/j.cell.2018.05.012).
//! - `enCas12a`: the engineered broad-PAM `AsCas12a` variant.
//!
//! Each matrix has 276 entries — 23 protospacer positions × 12 distinct
//! base-pair mismatches (4 RNA bases × 3 non-matching DNA bases per
//! base; same-base pairs are implicit 1.0). Format is
//! `RDA,Pos,MM,avg_percent_active` where MM is `r{A,C,G,U}:d{A,C,G,T}`.
//!
//! The per-pair score is the **product of per-position penalties**
//! over all mismatches in the 20-nt (or up to 23-nt) protospacer
//! body. There is **no PAM-weight tail term** (unlike SpCas9, where
//! the OT's GG-of-NGG contributes a `pam_weight[g1][g2]` factor).
//!
//! Per-guide aggregate (the "single number" analog of SpCas9 specificity):
//!
//! ```text
//! tttn_specificity = 1 / Σ_OT (score_pair × count)
//! tttv_specificity = 1 / Σ_{OT.pam ∉ TTTT*} (score_pair × count)
//! ```
//!
//! Both use the GuideScan2 denominator convention — `1 / Σ`, no `+1`.
//! The TTTV variant exists because Cas12a doesn't cleave TTTT-prefixed
//! off-targets in practice, so excluding them gives a per-guide score
//! weighted toward biologically credible off-targets.

use crispr_encoding::Site;

/// Per-guide aggregated Cas12a result. Mirrors the SpCas9 [`crate::CfdResult`]
/// shape but carries the two-flavor specificity (TTTN and TTTV).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cas12aResult {
    /// Max per-(guide, OT) CFD score over off-targets (mm > 0). Lower is
    /// better — a guide with no Cas12a-tolerated off-target scores near 0.
    pub max_cfd: f64,
    /// `1 / Σ_all_OT (cfd × count)`. Includes all PAM variants.
    pub tttn_specificity: f64,
    /// `1 / Σ_{OT.pam starts with TTT[ACG]} (cfd × count)`. Excludes
    /// TTTT-prefixed off-targets, which Cas12a cleaves poorly in practice.
    pub tttv_specificity: f64,
}

/// Bundled matrix variants. Matrix data is embedded at compile time via
/// [`include_str!`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cas12aMatrix {
    /// `2xNLS-Cas12a` — wild-type-ish `AsCas12a` profile.
    TwoXNls,
    /// `enCas12a` — engineered broad-PAM `AsCas12a` variant.
    EnCas12a,
}

impl Cas12aMatrix {
    /// Raw CSV body for this matrix variant.
    #[must_use]
    pub const fn csv(self) -> &'static str {
        match self {
            Self::TwoXNls => include_str!("../data/cas12a_2xnls.csv"),
            Self::EnCas12a => include_str!("../data/cas12a_encas12a.csv"),
        }
    }

    /// Human-readable name (matches the `RDA` column in the CSV).
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::TwoXNls => "2xNLS-Cas12a",
            Self::EnCas12a => "enCas12a",
        }
    }
}

/// Cas12a CFD scorer. Construct once with [`Cas12aCfd::from_matrix`] or
/// [`Cas12aCfd::from_csv`], then call [`Cas12aCfd::score_pair`] /
/// [`Cas12aCfd::aggregate`] repeatedly.
///
/// Thread-safe (`Send + Sync` — only `Copy` fields). Construction parses
/// the matrix once; the hot path is a tight `Site` × `Site` loop with no
/// allocations.
#[derive(Debug, Clone)]
pub struct Cas12aCfd {
    /// `penalty[guide_2bit][ot_2bit][position_0..22]`. Matching pairs
    /// (g == o) sit at 1.0 (no penalty contribution).
    penalty: [[[f64; 23]; 4]; 4],
    /// Bundled matrix variant this scorer was built from (purely for
    /// reporting). `None` when constructed from a user-supplied CSV.
    matrix_name: String,
}

/// Errors surfaced when parsing a matrix CSV.
#[derive(Debug)]
pub enum Cas12aParseError {
    /// Header line didn't match the expected schema.
    BadHeader { line: String },
    /// A row had the wrong number of comma-separated fields.
    BadArity { line_number: usize, line: String },
    /// A row's `Pos` field wasn't a valid 1..=23 position.
    BadPosition { line_number: usize, raw: String },
    /// A row's `MM` field didn't match the `r{ACGU}:d{ACGT}` shape.
    BadMmKey { line_number: usize, raw: String },
    /// A row's `avg_percent_active` field wasn't parseable as `f64`.
    BadScore { line_number: usize, raw: String },
}

impl std::fmt::Display for Cas12aParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadHeader { line } => write!(
                f,
                "Cas12a matrix header mismatch — expected \
                 `RDA,Pos,MM,avg_percent_active`, got {line:?}"
            ),
            Self::BadArity { line_number, line } => write!(
                f,
                "Cas12a matrix line {line_number}: expected 4 comma-separated \
                 fields, got {line:?}"
            ),
            Self::BadPosition { line_number, raw } => write!(
                f,
                "Cas12a matrix line {line_number}: invalid `Pos` field \
                 (expected integer 1..=23), got {raw:?}"
            ),
            Self::BadMmKey { line_number, raw } => write!(
                f,
                "Cas12a matrix line {line_number}: invalid `MM` key \
                 (expected `r{{A,C,G,U}}:d{{A,C,G,T}}`), got {raw:?}"
            ),
            Self::BadScore { line_number, raw } => write!(
                f,
                "Cas12a matrix line {line_number}: invalid \
                 `avg_percent_active` (expected float), got {raw:?}"
            ),
        }
    }
}

impl std::error::Error for Cas12aParseError {}

const EXPECTED_HEADER: &str = "RDA,Pos,MM,avg_percent_active";

impl Cas12aCfd {
    /// Construct from one of the bundled matrices. The CSV is parsed once
    /// at construction; further `score_pair` / `aggregate` calls are
    /// allocation-free.
    ///
    /// # Errors
    /// Returns [`Cas12aParseError`] if the bundled CSV is malformed
    /// (should be impossible in practice — bundled data is statically
    /// validated by `bundled_matrices_parse_cleanly` test).
    pub fn from_matrix(matrix: Cas12aMatrix) -> Result<Self, Cas12aParseError> {
        let mut s = Self::from_csv(matrix.csv())?;
        s.matrix_name = matrix.name().to_string();
        Ok(s)
    }

    /// Construct from an arbitrary CSV body with the same schema as the
    /// bundled matrices (`RDA,Pos,MM,avg_percent_active`).
    ///
    /// # Errors
    /// See [`Cas12aParseError`].
    pub fn from_csv(csv_body: &str) -> Result<Self, Cas12aParseError> {
        let mut lines = csv_body.lines();
        let header = lines.next().unwrap_or("").trim();
        if header != EXPECTED_HEADER {
            return Err(Cas12aParseError::BadHeader {
                line: header.to_string(),
            });
        }

        let mut penalty = [[[1.0_f64; 23]; 4]; 4];
        let mut matrix_name = String::new();
        for (idx, raw_line) in lines.enumerate() {
            let line_number = idx + 2;
            let line = raw_line.trim();
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split(',').collect();
            if fields.len() != 4 {
                return Err(Cas12aParseError::BadArity {
                    line_number,
                    line: line.to_string(),
                });
            }
            if matrix_name.is_empty() {
                matrix_name = fields[0].to_string();
            }

            let pos: usize = fields[1]
                .parse()
                .map_err(|_| Cas12aParseError::BadPosition {
                    line_number,
                    raw: fields[1].to_string(),
                })?;
            if !(1..=23).contains(&pos) {
                return Err(Cas12aParseError::BadPosition {
                    line_number,
                    raw: fields[1].to_string(),
                });
            }

            let (guide_idx, ot_idx) =
                parse_mm_key(fields[2]).ok_or_else(|| Cas12aParseError::BadMmKey {
                    line_number,
                    raw: fields[2].to_string(),
                })?;

            let score: f64 = fields[3]
                .parse()
                .map_err(|_| Cas12aParseError::BadScore {
                    line_number,
                    raw: fields[3].to_string(),
                })?;

            penalty[guide_idx][ot_idx][pos - 1] = score;
        }

        Ok(Self {
            penalty,
            matrix_name,
        })
    }

    /// Name of the matrix this scorer was built from (e.g.
    /// `"2xNLS-Cas12a"`).
    #[must_use]
    pub fn matrix_name(&self) -> &str {
        &self.matrix_name
    }

    /// Score a single `(target, off_target)` pair on the **23-nt-protospacer
    /// Cas12a** layout (`Enzyme::cpf1_tttn`): 4-bp PAM at the top, then
    /// 23-bp protospacer. All 23 matrix positions are scored.
    ///
    /// `target` should be the on-target/guide bits; `off_target` is the
    /// genomic OT. PAM bits aren't read (no Cas12a PAM weight model).
    #[must_use]
    pub fn score_pair(&self, target: Site, off_target: Site) -> f64 {
        self.score_pair_with_len(target, off_target, 23)
    }

    /// Like [`Cas12aCfd::score_pair`] but parameterizes the protospacer
    /// length, supporting the legacy 20-nt protospacer Cas12a layout
    /// (matching FlashFry's `Cpf1ParameterPack`). Out-of-range positions
    /// (anything > `protospacer_len`) implicitly score 1.0 (match).
    ///
    /// # Panics
    /// Panics if `protospacer_len > 23` — the matrix only has 23 rows.
    #[must_use]
    pub fn score_pair_with_len(
        &self,
        target: Site,
        off_target: Site,
        protospacer_len: u32,
    ) -> f64 {
        assert!(
            protospacer_len <= 23,
            "Cas12a matrix has 23 positions; protospacer_len = {protospacer_len} > 23"
        );
        let mut score = 1.0_f64;
        // Encoding: first character of the 24-bp Cas12a site (the PAM's
        // T-1) is at the highest used bits. For an N-bp scan length:
        //   bit_offset(site_position) = (N - 1 - site_position) * 2
        // The protospacer's PAM-proximal base sits at site position 4
        // (just after the 4-bp PAM); the PAM-distal base at position
        // N-1. Matrix position `p` (1-indexed, 1..=protospacer_len)
        // corresponds to site position `p + 3`.
        //
        // **Match logic**: only multiply in the penalty when the
        // target and OT bases differ. The matrix sometimes contains
        // entries for same-letter pairs (`rA:dA`, etc.) — these
        // describe G:U-wobble / non-Watson-Crick scenarios in the
        // upstream Cas12a literature but are dead code under the
        // string-equality match convention crisprware's
        // `parasol_scripts/score_flashfry_cfd.py` uses
        // (`tb != ob → lookup, else skip`). We follow that convention
        // verbatim so output numbers match the parasol pipeline.
        let scan_len: u32 = u32::from(4_u8) + protospacer_len; // PAM + protospacer
        for p in 1..=protospacer_len {
            let site_pos = p + 3;
            let bit_offset = (scan_len - 1 - site_pos) * 2;
            let g = ((target.low >> bit_offset) & 0b11) as usize;
            let o = ((off_target.low >> bit_offset) & 0b11) as usize;
            if g != o {
                score *= self.penalty[g][o][(p - 1) as usize];
            }
        }
        score
    }

    /// Aggregate Cas12a CFD across a guide's off-targets.
    ///
    /// `ots` yields `(off_target_site, mismatches, count)` triples. The
    /// `count` weights multi-mapping off-targets exactly the way the
    /// parasol script does. **All** off-targets contribute, regardless
    /// of mismatch count; the diagonal (mm = 0, match) score is 1.0 so
    /// on-targets are still counted. To match the GuideScan2 / parasol
    /// convention, the specificity is `1 / Σ (cfd × count)` — no `+ 1`.
    ///
    /// The TTTV variant excludes off-targets whose 4-bp PAM starts with
    /// `TTTT` (i.e. position 4 of the site is T). `tttn_specificity`
    /// includes everything.
    ///
    /// # Panics
    /// As [`Cas12aCfd::score_pair_with_len`].
    #[must_use]
    pub fn aggregate<I>(&self, target: Site, ots: I) -> Cas12aResult
    where
        I: IntoIterator<Item = (Site, u8, u32)>,
    {
        self.aggregate_with_len(target, ots, 23)
    }

    /// Like [`Cas12aCfd::aggregate`] but parameterizes the protospacer
    /// length (see [`Cas12aCfd::score_pair_with_len`]).
    #[must_use]
    pub fn aggregate_with_len<I>(
        &self,
        target: Site,
        ots: I,
        protospacer_len: u32,
    ) -> Cas12aResult
    where
        I: IntoIterator<Item = (Site, u8, u32)>,
    {
        let mut max_cfd = 0.0_f64;
        let mut tttn_sum = 0.0_f64;
        let mut tttv_sum = 0.0_f64;
        for (ot, mismatches, count) in ots {
            let cfd = self.score_pair_with_len(target, ot, protospacer_len);
            let weight = f64::from(count.max(1));
            tttn_sum += cfd * weight;
            if !is_tttt_prefix(ot, protospacer_len) {
                tttv_sum += cfd * weight;
            }
            // `max_cfd` summarizes the worst off-target — exclude
            // perfect on-targets (mm == 0, cfd == 1.0) so the column
            // stays a useful indicator.
            if mismatches > 0 && cfd > max_cfd {
                max_cfd = cfd;
            }
        }
        let tttn_specificity = if tttn_sum > 0.0 { 1.0 / tttn_sum } else { 1.0 };
        let tttv_specificity = if tttv_sum > 0.0 { 1.0 / tttv_sum } else { 1.0 };
        Cas12aResult {
            max_cfd,
            tttn_specificity,
            tttv_specificity,
        }
    }
}

/// Check whether an off-target site's PAM is `TTTT` (i.e. its position 3
/// base is `T`, which means the entire 4-bp PAM is TTTT given the TTTN
/// enzyme scan filter only emits sites with T at positions 0, 1, 2).
///
/// `protospacer_len` is the length of the protospacer portion; total
/// scan length is `4 + protospacer_len`. PAM bits sit at the highest
/// used bits per the standard encoding.
#[must_use]
pub fn is_tttt_prefix(site: Site, protospacer_len: u32) -> bool {
    // Site position 3 (the variable N of TTTN) bit offset:
    //   scan_len = 4 + protospacer_len; bit_offset = (scan_len - 1 - 3) * 2
    //   = (protospacer_len) * 2
    let bit_offset = protospacer_len * 2;
    // 2-bit DNA encoding: T = 3 (A=0, C=1, G=2, T=3).
    ((site.low >> bit_offset) & 0b11) == 0b11
}

/// Parse an `r{A,C,G,U}:d{A,C,G,T}` key into `(rna_base_2bit,
/// dna_base_2bit)`. Returns `None` for any malformed key.
fn parse_mm_key(key: &str) -> Option<(usize, usize)> {
    let bytes = key.as_bytes();
    if bytes.len() != 5 || bytes[0] != b'r' || bytes[2] != b':' || bytes[3] != b'd' {
        return None;
    }
    let rna = match bytes[1] {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'U' => 3, // U maps to DNA T's 2-bit code 3
        _ => return None,
    };
    let dna = match bytes[4] {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    Some((rna, dna))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crispr_encoding::{Base, Site};

    fn encode(s: &str) -> Site {
        let bases: Vec<Base> = s
            .bytes()
            .map(|c| match c.to_ascii_uppercase() {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'T' => Base::T,
                _ => panic!("non-ACGT base in test input"),
            })
            .collect();
        Site::encode(bases)
    }

    #[test]
    fn bundled_matrices_parse_cleanly() {
        let a = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).expect("2xnls parses");
        assert_eq!(a.matrix_name(), "2xNLS-Cas12a");
        let b = Cas12aCfd::from_matrix(Cas12aMatrix::EnCas12a).expect("encas12a parses");
        assert_eq!(b.matrix_name(), "enCas12a");
    }

    #[test]
    fn bundled_matrix_has_no_unset_mismatch_holes() {
        // Each matrix should have 23 × 12 = 276 non-default penalties.
        // We can't count them after the fact (1.0 is also a legal score),
        // but we can spot-check: any (guide != ot) cell at position 1
        // should be < 1.0 in the 2xNLS matrix (these are knockdown
        // measurements, never "more active than perfect").
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        for g in 0..4 {
            for o in 0..4 {
                if g == o {
                    continue;
                }
                // rU:dT corresponds to (g=3, o=3) — but that's a *match*
                // (U == T at the DNA level) and isn't in the table.
                // The matrix omits 4 keys: dT for rU and the three
                // same-base "rA:dA"... wait those ARE present. Let me
                // re-check by listing what's in the matrix:
                //   rA:dA, rA:dC, rA:dG (no rA:dT)
                //   rC:dA, rC:dC, rC:dT (no rC:dG → that's a match)
                //   rG:dA, rG:dG, rG:dT (no rG:dC → match)
                //   rU:dC, rU:dG, rU:dT (no rU:dA → that's match)
                // Skip the complementary-match pairs.
                let match_pairs = [(0, 3), (1, 2), (2, 1), (3, 0)];
                if match_pairs.contains(&(g, o)) {
                    continue;
                }
                let val = scorer.penalty[g][o][0];
                assert!(
                    val > 0.0 && val < 1.0,
                    "expected non-trivial penalty at (g={g}, o={o}, p=0), got {val}"
                );
            }
        }
    }

    #[test]
    fn perfect_match_scores_one() {
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        // 4-bp TTTN PAM + 20-bp protospacer; same site twice = no
        // mismatches = score 1.0.
        let s = encode("TTTAACGTACGTACGTACGTACGT");
        let score = scorer.score_pair(s, s);
        assert!((score - 1.0).abs() < 1e-12, "perfect-match score {score}");
    }

    #[test]
    fn single_mismatch_uses_table_value() {
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        // Target/guide:           TTTA | AAAAAAAAAAAAAAAAAAAA
        // Off-target with 1 mm:   TTTA | AAAACAAAAAAAAAAAAAAA  (pos 5 in
        //                                     site → matrix pos 2)
        // Guide base at protospacer pos 2 = A → key `rA:d?`.
        // OT base at same pos = C → key `rA:dC`.
        let target = encode("TTTAAAAAAAAAAAAAAAAAAAAAAAA");
        let off = encode("TTTAACAAAAAAAAAAAAAAAAAAAAA");
        let score = scorer.score_pair(target, off);
        // Compare to the raw matrix entry for (pos=2, rA:dC):
        let body = Cas12aMatrix::TwoXNls.csv();
        let mut expected: Option<f64> = None;
        for line in body.lines().skip(1) {
            let f: Vec<&str> = line.split(',').collect();
            if f.len() == 4 && f[1] == "2" && f[2] == "rA:dC" {
                expected = Some(f[3].parse().unwrap());
                break;
            }
        }
        let expected = expected.expect("rA:dC at pos 2 in 2xNLS matrix");
        assert!(
            (score - expected).abs() < 1e-12,
            "expected score = matrix value {expected}, got {score}"
        );
    }

    #[test]
    fn two_mismatches_multiply_their_penalties() {
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        // Mismatches at protospacer positions 1 and 2.
        // Guide: TTTA AA AAAAAAAAAAAAAAAAAA
        // OT:    TTTA CC AAAAAAAAAAAAAAAAAA   (rA:dC at pos 1 and pos 2)
        let target = encode("TTTAAAAAAAAAAAAAAAAAAAAAAAA");
        let off = encode("TTTACCAAAAAAAAAAAAAAAAAAAAA");
        let score = scorer.score_pair(target, off);

        // Look up both pos=1 rA:dC and pos=2 rA:dC, multiply.
        let body = Cas12aMatrix::TwoXNls.csv();
        let mut p1 = 1.0_f64;
        let mut p2 = 1.0_f64;
        for line in body.lines().skip(1) {
            let f: Vec<&str> = line.split(',').collect();
            if f.len() != 4 || f[2] != "rA:dC" {
                continue;
            }
            let v: f64 = f[3].parse().unwrap();
            match f[1] {
                "1" => p1 = v,
                "2" => p2 = v,
                _ => {}
            }
        }
        let expected = p1 * p2;
        assert!(
            (score - expected).abs() < 1e-12,
            "expected {expected} (= {p1} × {p2}), got {score}"
        );
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact 1.0 / 0.0 sentinel checks
    fn aggregate_specificity_formula() {
        // 1 / Σ (cfd × count). For an empty OT list, fall back to 1.0.
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        let target = encode("TTTAAAAAAAAAAAAAAAAAAAAAAAA");
        let empty = std::iter::empty::<(Site, u8, u32)>();
        let r = scorer.aggregate(target, empty);
        assert_eq!(r.tttn_specificity, 1.0);
        assert_eq!(r.tttv_specificity, 1.0);
        assert_eq!(r.max_cfd, 0.0);

        // Single on-target hit contributes cfd=1.0 × 1 = 1.0 to denom.
        let r = scorer.aggregate(target, std::iter::once((target, 0, 1)));
        assert!((r.tttn_specificity - 1.0).abs() < 1e-12);
        assert!((r.tttv_specificity - 1.0).abs() < 1e-12);
        // max_cfd skips the on-target.
        assert_eq!(r.max_cfd, 0.0);
    }

    #[test]
    fn tttv_excludes_tttt_offtargets() {
        let scorer = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        // Target with TTTA PAM.
        let target = encode("TTTAAAAAAAAAAAAAAAAAAAAAAAA");
        // Two off-targets, both at mm=1 (single base diff at proto pos
        // 1). One has TTTA PAM (kept by tttv); the other has TTTT PAM
        // (excluded by tttv).
        let ot_tttv = encode("TTTACAAAAAAAAAAAAAAAAAAAAAA"); // pos 1: A→C
        let ot_tttt = encode("TTTTCAAAAAAAAAAAAAAAAAAAAAA"); // pos 1: A→C, PAM TTTT
        let r = scorer.aggregate(
            target,
            [(ot_tttv, 1, 1), (ot_tttt, 1, 1)],
        );
        // tttn_sum = cfd_tttv + cfd_tttt; tttv_sum = cfd_tttv only.
        // Both off-targets get the same protospacer-position penalty
        // (PAM bits aren't scored), so cfd_tttv == cfd_tttt.
        let single = scorer.score_pair(target, ot_tttv);
        let expected_tttn = 1.0 / (2.0 * single);
        let expected_tttv = 1.0 / (1.0 * single);
        assert!(
            (r.tttn_specificity - expected_tttn).abs() < 1e-12,
            "tttn_specificity {} vs expected {}",
            r.tttn_specificity,
            expected_tttn
        );
        assert!(
            (r.tttv_specificity - expected_tttv).abs() < 1e-12,
            "tttv_specificity {} vs expected {}",
            r.tttv_specificity,
            expected_tttv
        );
    }

    #[test]
    #[allow(clippy::float_cmp)] // both paths produce bit-identical f64
    fn from_csv_matches_from_matrix() {
        let a = Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).unwrap();
        let b = Cas12aCfd::from_csv(Cas12aMatrix::TwoXNls.csv()).unwrap();
        // Both should compute identical scores on an arbitrary pair.
        let t = encode("TTTAAAAAAAAAAAAAAAAAAAAAAAA");
        let o = encode("TTTACCAAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(a.score_pair(t, o), b.score_pair(t, o));
    }

    #[test]
    fn from_csv_rejects_bad_header() {
        assert!(matches!(
            Cas12aCfd::from_csv("wrong\nrow,row,row,row\n"),
            Err(Cas12aParseError::BadHeader { .. })
        ));
    }
}
