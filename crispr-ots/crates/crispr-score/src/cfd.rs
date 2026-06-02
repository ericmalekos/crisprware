//! Doench 2016 CFD (Cutting Frequency Determination) scoring.
//!
//! CFD is a per-mismatch penalty model trained on SpCas9 cleavage data. It
//! is **SpCas9-specific** and **assumes a 20-nt protospacer + 3-nt PAM** —
//! see `NOTES.md` for the science caveats around applying it to anything
//! else (Cas12a, longer guides, alternative PAMs).
//!
//! ## Why this looks different from FlashFry's `Doench2016CFDScore.scala`
//!
//! FlashFry's CFD inner loop decodes each off-target's 64-bit encoding to
//! a `String`, walks it character-by-character against the guide string,
//! and looks up `HashMap[String, Double]` keyed on `"rC:dG,9"`-style
//! formatted strings. That's three of our four "Tier 1" hot-path hazards
//! at once (string allocation per OT, `HashMap` lookup per mismatch,
//! redundant string formatting per lookup — see optimization #3 in
//! `snappy-skipping-pillow.md`).
//!
//! Our `Cfd` operates directly on bit fields of the encoded `Site`. The
//! mismatch table is a flat `[[[f64; 20]; 4]; 4]` indexed by
//! `(guide_2bit, ot_2bit, position)`, populated once at construction from
//! FlashFry's published values. No strings, no maps, no formatting in the
//! hot path — just `(low >> bit_offset) & 0b11`, two array indexes, and a
//! multiply.

use crispr_encoding::Site;

/// Per-guide aggregated CFD result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CfdResult {
    /// Maximum CFD score across the guide's off-targets (excluding the
    /// on-target). Range `[0.0, 1.0]`. **Lower is better** (good guide:
    /// no off-target is well-tolerated by Cas9).
    pub max_cfd: f64,
    /// Aggregate specificity. Range `(0.0, 1.0]`. **Higher is better.**
    /// The exact formula depends on the [`SpecConvention`] passed to
    /// [`Cfd::aggregate`] — see that enum's docs.
    pub specificity: f64,
}

/// Which paper's specificity-denominator convention to use.
///
/// CFD itself (the per-(guide, OT) score) is well-defined and identical
/// across implementations. What differs is how a per-guide *aggregate*
/// specificity is rolled up from those per-pair scores:
///
/// - **FlashFry** publishes `1 / (1 + Σ_OT cfd × count)`. The `+ 1`
///   accounts for the on-target as a "perfect" CFD contribution; the
///   sum is over off-targets only. A guide with no off-targets scores
///   exactly 1.0.
/// - **GuideScan2** publishes `1 / Σ_all_rows cfd`. The sum includes
///   every row their enumerator emits — the on-target site and any
///   multi-mapping copies of it (each contributing CFD = 1.0), plus
///   off-targets. A guide with `n` on-target copies and no off-targets
///   scores `1/n`.
///
/// The two formulas differ only by which terms appear in the
/// denominator; the per-pair CFD math is the same.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecConvention {
    /// `1 / (1 + Σ_offtarget cfd × count)`. The FlashFry convention.
    Flashfry,
    /// `1 / Σ_all_rows cfd`, i.e. on-target multi-mappers and off-targets
    /// each contribute their (count × cfd) to the denominator. Matches
    /// the GuideScan2 paper.
    Guidescan,
}

/// CFD scorer. Construct once with `Cfd::new()`, then call `score_pair`
/// or `aggregate` repeatedly. Thread-safe (`Send + Sync` — only `Copy`
/// fields).
#[derive(Debug, Clone)]
pub struct Cfd {
    /// `penalty[guide_2bit][ot_2bit][position_0indexed]`. For matching
    /// pairs (guide == ot) the value is `1.0` (no penalty contribution).
    penalty: [[[f64; 20]; 4]; 4],
    /// `pam_weight[pam_pos1_2bit][pam_pos2_2bit]` for the last 2 bases of
    /// the 3-bp PAM (the GG of NGG). The first base of the PAM (the N) is
    /// not part of the lookup.
    pam_weight: [[f64; 4]; 4],
}

impl Default for Cfd {
    fn default() -> Self {
        Self::new()
    }
}

impl Cfd {
    /// Build the scorer. Populates the lookup tables from FlashFry's
    /// published values (see [`RAW_MM_PENALTIES`] / [`RAW_PAM_WEIGHTS`]).
    #[must_use]
    pub fn new() -> Self {
        let mut penalty = [[[1.0_f64; 20]; 4]; 4];
        for &(key, score) in RAW_MM_PENALTIES {
            let (g, o, p) = parse_mm_key(key);
            penalty[g][o][p] = score;
        }
        let mut pam_weight = [[0.0_f64; 4]; 4];
        for &(key, score) in RAW_PAM_WEIGHTS {
            let (b1, b2) = parse_pam_key(key);
            pam_weight[b1][b2] = score;
        }
        Self {
            penalty,
            pam_weight,
        }
    }

    /// The mismatch penalty table, indexed `[guide_2bit][ot_2bit][position]`.
    /// Exposed so out-of-crate scorers (e.g. the GPU kernel) can upload the
    /// exact same FlashFry-derived values rather than re-transcribing them.
    #[must_use]
    pub fn penalty_table(&self) -> &[[[f64; 20]; 4]; 4] {
        &self.penalty
    }

    /// The PAM-tail weight table, indexed `[pam_pos1_2bit][pam_pos2_2bit]`
    /// (the last two bases of the 3-bp PAM). Companion to
    /// [`penalty_table`](Self::penalty_table) for the GPU table upload.
    #[must_use]
    pub fn pam_table(&self) -> &[[f64; 4]; 4] {
        &self.pam_weight
    }

    /// Score a single `(guide, off-target)` pair. **Both** sites must be
    /// encoded as 23-bp (20-nt protospacer + 3-nt PAM) — the layout
    /// produced by `Enzyme::spcas9_ngg()` and friends.
    ///
    /// Returns the multiplied per-position mismatch penalty × PAM weight.
    /// In the limit where every protospacer base matches and the OT PAM
    /// is NGG, this is `1.0` (full cleavage); for a poor off-target it's
    /// near zero.
    #[must_use]
    pub fn score_pair(&self, guide: Site, off_target: Site) -> f64 {
        let mut score = 1.0_f64;
        // Iterate the 20 protospacer positions. In a 23-bp Cas9 site:
        // - Position 0 (5' protospacer, "PAM-distal") sits at bits 44-45.
        // - Position 19 (3' protospacer, "PAM-proximal") sits at bits 6-7.
        // The PAM occupies bits 0-5. (See encoding docs in `Site`.)
        for p in 0..20_u32 {
            let bit_offset = (22 - p) * 2;
            let g = ((guide.low >> bit_offset) & 0b11) as usize;
            let o = ((off_target.low >> bit_offset) & 0b11) as usize;
            score *= self.penalty[g][o][p as usize];
        }
        // PAM tail: last two bases of the OT PAM (the GG of NGG), at
        // bits 0-3 of `low`. FlashFry's lookup is keyed on the string-order
        // pair, i.e. "{base_at_pos1}{base_at_pos2}" → bits 2-3 then bits 0-1.
        let pam_pos1 = ((off_target.low >> 2) & 0b11) as usize;
        let pam_pos2 = (off_target.low & 0b11) as usize;
        score * self.pam_weight[pam_pos1][pam_pos2]
    }

    /// Aggregate CFD across a guide's off-targets using the FlashFry
    /// convention. Equivalent to `aggregate_with(guide, ots,
    /// SpecConvention::Flashfry)`; kept as the original signature so
    /// existing callers don't break.
    #[must_use]
    pub fn aggregate<I>(&self, guide: Site, ots: I) -> CfdResult
    where
        I: IntoIterator<Item = (Site, u8, u32)>,
    {
        self.aggregate_with(guide, ots, SpecConvention::Flashfry)
    }

    /// Aggregate CFD across a guide's off-targets under the requested
    /// specificity convention.
    ///
    /// `ots` yields `(site, mismatches, count)` triples — typically built
    /// from the scanner's `Hit` results after collapsing by off-target
    /// sequence so that `count` is the number of *genomic positions*
    /// sharing that sequence.
    ///
    /// Behavior under each convention:
    /// - [`SpecConvention::Flashfry`]: off-target rows (`mismatches > 0`)
    ///   contribute `cfd × count` to the denominator; on-target rows
    ///   are skipped because the `+ 1` term in the closed form already
    ///   accounts for one perfect on-target. `max_cfd` is computed over
    ///   off-targets only.
    /// - [`SpecConvention::Guidescan`]: every row contributes
    ///   `cfd × count` to the denominator — on-target rows contribute
    ///   `1.0 × count` each (their `score_pair` is identically 1.0).
    ///   `max_cfd` is still computed over off-targets only so the column
    ///   stays meaningful as a worst-off-target indicator.
    ///
    /// The explicit `count` argument matters: a single multi-mapping
    /// off-target with 6 genomic positions contributes 6 × cfd to the
    /// denominator. Using `Site::count()` here would silently default to
    /// 1 in pipelines (like ours) that don't pre-merge identical
    /// sequences during DB build.
    #[must_use]
    pub fn aggregate_with<I>(
        &self,
        guide: Site,
        ots: I,
        convention: SpecConvention,
    ) -> CfdResult
    where
        I: IntoIterator<Item = (Site, u8, u32)>,
    {
        let mut max_cfd = 0.0_f64;
        let mut off_sum = 0.0_f64;
        let mut on_count: u32 = 0;
        let mut any = false;
        for (ot, mismatches, count) in ots {
            let weight = f64::from(count.max(1));
            if mismatches == 0 {
                on_count = on_count.saturating_add(count.max(1));
                continue;
            }
            let cfd = self.score_pair(guide, ot);
            off_sum += cfd * weight;
            if cfd > max_cfd {
                max_cfd = cfd;
            }
            any = true;
        }
        let specificity = match convention {
            SpecConvention::Flashfry => {
                if any && off_sum > 0.0 {
                    1.0 / (1.0 + off_sum)
                } else {
                    1.0
                }
            }
            SpecConvention::Guidescan => {
                // GuideScan2 includes one CFD-1.0 contribution per
                // on-target row. A guide with zero on-target rows
                // (shouldn't happen in our pipeline, but handle it)
                // falls back to the FlashFry +1 so specificity stays
                // well-defined.
                let on_term = if on_count == 0 {
                    1.0
                } else {
                    f64::from(on_count)
                };
                let denom = on_term + off_sum;
                if denom > 0.0 {
                    1.0 / denom
                } else {
                    1.0
                }
            }
        };
        CfdResult {
            max_cfd,
            specificity,
        }
    }
}

/// Parse a FlashFry-style mismatch key `"rX:dY,N"` into `(guide_idx,
/// ot_idx, position_0indexed)`.
///
/// Key layout:
/// - `rX` where X is the RNA guide base (A, C, G, U); maps to DNA 2-bit
///   code (U → T).
/// - `dY` where Y is `specialReverseCompBase(ot_base_after_T_to_U)`. The
///   table key uses the *complement of the OT DNA base after T→U folding*
///   so e.g. an OT with `T` at position 5 produces key `dA` (U → A);
///   `A` produces `dT`; `C` produces `dG`; `G` produces `dC`. We invert
///   that mapping here to recover the OT's 2-bit code.
/// - `,N` where N is 1-indexed position 1..=20.
const fn parse_mm_key(key: &str) -> (usize, usize, usize) {
    let bytes = key.as_bytes();
    // Layout: r X : d Y , [digits]
    debug_assert_byte(bytes[0], b'r');
    let r_char = bytes[1];
    debug_assert_byte(bytes[2], b':');
    debug_assert_byte(bytes[3], b'd');
    let d_char = bytes[4];
    debug_assert_byte(bytes[5], b',');

    let guide_idx = match r_char {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'U' => 3,
        _ => panic!("invalid r-base in mismatch key"),
    };
    // Invert specialReverseCompBase (A→T, C→G, G→C, U→A) to recover the
    // OT's 2-bit code.
    let ot_idx = match d_char {
        b'T' => 0, // dT comes from OT base A
        b'G' => 1, // dG ← C
        b'C' => 2, // dC ← G
        b'A' => 3, // dA ← T (T→U→A)
        _ => panic!("invalid d-base in mismatch key"),
    };

    let pos = parse_position(bytes, 6);
    (guide_idx, ot_idx, pos - 1)
}

/// Parse the 1-indexed position number from `bytes[start..]`.
const fn parse_position(bytes: &[u8], start: usize) -> usize {
    let mut acc = 0_usize;
    let mut i = start;
    while i < bytes.len() {
        let b = bytes[i];
        if b < b'0' || b > b'9' {
            break;
        }
        acc = acc * 10 + (b - b'0') as usize;
        i += 1;
    }
    acc
}

/// Parse a 2-character PAM-tail key like `"GG"` into the 2-bit indices.
const fn parse_pam_key(key: &str) -> (usize, usize) {
    let bytes = key.as_bytes();
    (ascii_base_to_idx(bytes[0]), ascii_base_to_idx(bytes[1]))
}

const fn ascii_base_to_idx(c: u8) -> usize {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => panic!("invalid base in PAM key"),
    }
}

const fn debug_assert_byte(actual: u8, expected: u8) {
    assert!(actual == expected, "malformed CFD lookup key");
}

/// Raw mismatch penalties transcribed verbatim from FlashFry's
/// `Doench2016CFDScore.scala:173-210`. Format: `(key, multiplier)`.
///
/// The float literals carry many significant digits from the Python pickle
/// FlashFry originally read from; rather than re-grouping each one we
/// silence the readability lint at the item level.
#[allow(clippy::unreadable_literal)]
const RAW_MM_PENALTIES: &[(&str, f64)] = &[
    ("rC:dC,9", 0.619047619),
    ("rC:dC,8", 0.642857143),
    ("rG:dA,8", 0.625),
    ("rG:dG,19", 0.448275862),
    ("rG:dG,18", 0.476190476),
    ("rG:dG,15", 0.272727273),
    ("rG:dG,14", 0.428571429),
    ("rG:dG,17", 0.235294118),
    ("rG:dG,16", 0.0),
    ("rC:dC,20", 0.058823529),
    ("rG:dT,20", 0.9375),
    ("rG:dG,13", 0.421052632),
    ("rG:dG,12", 0.529411765),
    ("rC:dT,13", 0.384615385),
    ("rC:dT,18", 0.538461538),
    ("rC:dC,3", 0.5),
    ("rU:dG,12", 0.947368421),
    ("rA:dG,13", 0.210526316),
    ("rA:dG,12", 0.263157895),
    ("rA:dG,11", 0.4),
    ("rA:dG,10", 0.333333333),
    ("rA:dA,19", 0.538461538),
    ("rA:dA,18", 0.5),
    ("rA:dG,15", 0.272727273),
    ("rA:dG,14", 0.214285714),
    ("rA:dA,15", 0.2),
    ("rA:dA,14", 0.533333333),
    ("rA:dA,17", 0.133333333),
    ("rA:dA,16", 0.0),
    ("rA:dA,11", 0.307692308),
    ("rA:dA,10", 0.882352941),
    ("rA:dA,13", 0.3),
    ("rA:dA,12", 0.333333333),
    ("rG:dA,13", 0.3),
    ("rG:dA,12", 0.384615385),
    ("rG:dA,11", 0.384615385),
    ("rG:dA,10", 0.8125),
    ("rG:dA,17", 0.25),
    ("rG:dA,16", 0.0),
    ("rG:dA,15", 0.142857143),
    ("rG:dA,14", 0.266666667),
    ("rG:dT,10", 0.933333333),
    ("rG:dA,19", 0.666666667),
    ("rG:dA,18", 0.666666667),
    ("rC:dT,15", 0.222222222),
    ("rA:dC,4", 0.625),
    ("rA:dG,9", 0.571428571),
    ("rA:dG,8", 0.428571429),
    ("rC:dA,3", 0.6875),
    ("rC:dA,2", 0.909090909),
    ("rC:dA,1", 1.0),
    ("rC:dA,7", 0.8125),
    ("rC:dA,6", 0.928571429),
    ("rC:dA,5", 0.636363636),
    ("rC:dA,4", 0.8),
    ("rC:dA,9", 0.875),
    ("rC:dA,8", 0.875),
    ("rA:dG,20", 0.227272727),
    ("rA:dG,3", 0.428571429),
    ("rA:dG,2", 0.785714286),
    ("rG:dA,20", 0.7),
    ("rC:dT,20", 0.5),
    ("rC:dT,12", 0.714285714),
    ("rG:dT,17", 0.933333333),
    ("rC:dA,17", 0.466666667),
    ("rC:dA,16", 0.307692308),
    ("rC:dA,15", 0.066666667),
    ("rC:dA,14", 0.733333333),
    ("rC:dA,13", 0.7),
    ("rC:dA,12", 0.538461538),
    ("rC:dA,11", 0.307692308),
    ("rC:dA,10", 0.941176471),
    ("rG:dG,11", 0.428571429),
    ("rA:dC,20", 0.764705882),
    ("rC:dA,19", 0.461538462),
    ("rG:dG,10", 0.4),
    ("rU:dG,17", 0.705882353),
    ("rU:dG,16", 0.666666667),
    ("rU:dG,15", 0.272727273),
    ("rU:dG,14", 0.285714286),
    ("rU:dG,13", 0.789473684),
    ("rU:dC,20", 0.176470588),
    ("rU:dG,11", 0.666666667),
    ("rU:dG,10", 0.533333333),
    ("rG:dA,7", 0.571428571),
    ("rG:dA,6", 0.666666667),
    ("rG:dA,5", 0.3),
    ("rG:dA,4", 0.363636364),
    ("rG:dA,3", 0.5),
    ("rG:dA,2", 0.636363636),
    ("rG:dA,1", 1.0),
    ("rG:dT,7", 1.0),
    ("rG:dT,4", 0.9),
    ("rG:dG,6", 0.681818182),
    ("rU:dT,20", 0.5625),
    ("rC:dC,15", 0.05),
    ("rC:dC,14", 0.0),
    ("rC:dC,17", 0.058823529),
    ("rC:dC,16", 0.153846154),
    ("rC:dC,11", 0.25),
    ("rC:dC,10", 0.388888889),
    ("rC:dC,13", 0.136363636),
    ("rC:dC,12", 0.444444444),
    ("rC:dA,20", 0.3),
    ("rC:dC,19", 0.125),
    ("rC:dC,18", 0.133333333),
    ("rA:dA,1", 1.0),
    ("rA:dA,3", 0.705882353),
    ("rA:dA,2", 0.727272727),
    ("rA:dA,5", 0.363636364),
    ("rA:dA,4", 0.636363636),
    ("rA:dA,7", 0.4375),
    ("rA:dA,6", 0.714285714),
    ("rA:dA,9", 0.6),
    ("rA:dA,8", 0.428571429),
    ("rU:dG,20", 0.090909091),
    ("rU:dT,12", 0.8),
    ("rU:dT,13", 0.692307692),
    ("rU:dT,10", 0.857142857),
    ("rU:dT,11", 0.75),
    ("rU:dT,16", 0.909090909),
    ("rU:dT,17", 0.533333333),
    ("rU:dT,14", 0.619047619),
    ("rU:dT,15", 0.578947368),
    ("rC:dC,1", 0.913043478),
    ("rU:dT,18", 0.666666667),
    ("rC:dC,2", 0.695652174),
    ("rC:dC,5", 0.6),
    ("rC:dC,4", 0.5),
    ("rC:dC,7", 0.470588235),
    ("rC:dC,6", 0.5),
    ("rA:dC,2", 0.8),
    ("rU:dT,8", 0.8),
    ("rU:dT,9", 0.928571429),
    ("rA:dC,11", 0.65),
    ("rA:dC,19", 0.375),
    ("rA:dC,18", 0.4),
    ("rA:dC,17", 0.176470588),
    ("rA:dC,16", 0.192307692),
    ("rA:dC,15", 0.65),
    ("rU:dT,3", 0.714285714),
    ("rU:dT,4", 0.476190476),
    ("rC:dA,18", 0.642857143),
    ("rU:dT,6", 0.866666667),
    ("rA:dC,10", 0.555555556),
    ("rC:dT,10", 0.866666667),
    ("rU:dT,5", 0.5),
    ("rC:dT,8", 0.65),
    ("rC:dT,9", 0.857142857),
    ("rC:dT,6", 0.928571429),
    ("rC:dT,7", 0.75),
    ("rC:dT,4", 0.842105263),
    ("rC:dT,5", 0.571428571),
    ("rC:dT,2", 0.727272727),
    ("rA:dC,9", 0.666666667),
    ("rC:dT,1", 1.0),
    ("rA:dC,8", 0.733333333),
    ("rU:dT,1", 1.0),
    ("rC:dT,14", 0.35),
    ("rU:dT,2", 0.846153846),
    ("rU:dG,19", 0.275862069),
    ("rG:dT,14", 0.75),
    ("rG:dT,15", 0.941176471),
    ("rG:dT,16", 1.0),
    ("rA:dC,14", 0.466666667),
    ("rG:dG,20", 0.428571429),
    ("rG:dT,11", 1.0),
    ("rG:dT,12", 0.933333333),
    ("rG:dT,13", 0.923076923),
    ("rA:dG,7", 0.4375),
    ("rA:dC,13", 0.652173913),
    ("rA:dG,5", 0.5),
    ("rA:dG,4", 0.352941176),
    ("rG:dT,18", 0.692307692),
    ("rG:dT,19", 0.714285714),
    ("rA:dG,1", 0.857142857),
    ("rA:dC,12", 0.722222222),
    ("rG:dG,1", 0.714285714),
    ("rG:dT,3", 0.75),
    ("rG:dG,3", 0.384615385),
    ("rG:dG,2", 0.692307692),
    ("rG:dG,5", 0.785714286),
    ("rG:dG,4", 0.529411765),
    ("rG:dG,7", 0.6875),
    ("rG:dT,5", 0.866666667),
    ("rG:dG,9", 0.538461538),
    ("rG:dG,8", 0.615384615),
    ("rG:dT,8", 1.0),
    ("rG:dT,9", 0.642857143),
    ("rU:dG,18", 0.428571429),
    ("rU:dT,7", 0.875),
    ("rA:dG,6", 0.454545455),
    ("rG:dT,6", 1.0),
    ("rA:dA,20", 0.6),
    ("rA:dC,5", 0.72),
    ("rA:dG,17", 0.176470588),
    ("rU:dC,8", 0.733333333),
    ("rA:dG,16", 0.0),
    ("rG:dT,2", 0.846153846),
    ("rA:dG,19", 0.206896552),
    ("rU:dG,3", 0.428571429),
    ("rU:dG,2", 0.857142857),
    ("rU:dG,1", 0.857142857),
    ("rA:dG,18", 0.19047619),
    ("rU:dG,7", 0.6875),
    ("rU:dG,6", 0.909090909),
    ("rU:dG,5", 1.0),
    ("rU:dG,4", 0.647058824),
    ("rU:dG,9", 0.923076923),
    ("rU:dG,8", 1.0),
    ("rC:dT,11", 0.75),
    ("rC:dT,3", 0.866666667),
    ("rU:dC,19", 0.25),
    ("rU:dC,18", 0.333333333),
    ("rU:dC,13", 0.260869565),
    ("rU:dC,12", 0.5),
    ("rU:dC,11", 0.4),
    ("rU:dC,10", 0.5),
    ("rU:dC,17", 0.117647059),
    ("rU:dC,16", 0.346153846),
    ("rU:dC,15", 0.05),
    ("rU:dC,14", 0.0),
    ("rU:dC,7", 0.588235294),
    ("rU:dC,6", 0.571428571),
    ("rU:dC,5", 0.64),
    ("rU:dC,4", 0.625),
    ("rU:dC,3", 0.5),
    ("rU:dC,2", 0.84),
    ("rU:dC,1", 0.956521739),
    ("rC:dT,17", 0.466666667),
    ("rA:dC,3", 0.611111111),
    ("rC:dT,19", 0.428571429),
    ("rA:dC,1", 1.0),
    ("rA:dC,7", 0.705882353),
    ("rA:dC,6", 0.714285714),
    ("rU:dC,9", 0.619047619),
    ("rG:dA,9", 0.533333333),
    ("rU:dT,19", 0.285714286),
    ("rC:dT,16", 1.0),
    ("rG:dT,1", 0.9),
];

/// Raw PAM-tail weights transcribed verbatim from FlashFry's
/// `Doench2016CFDScore.scala:212-214`.
#[allow(clippy::unreadable_literal)]
const RAW_PAM_WEIGHTS: &[(&str, f64)] = &[
    ("AA", 0.0),
    ("AC", 0.0),
    ("GT", 0.016129032),
    ("AG", 0.259259259),
    ("CC", 0.0),
    ("CA", 0.0),
    ("CG", 0.107142857),
    ("TC", 0.0),
    ("GG", 1.0),
    ("GC", 0.022222222),
    ("AT", 0.0),
    ("GA", 0.069444444),
    ("TG", 0.038961039),
    ("CT", 0.0),
    ("TT", 0.0),
    ("TA", 0.0),
];

/// CRISPOR's recommended threshold below which an off-target's CFD is
/// considered negligible. FlashFry uses the same value at
/// `Doench2016CFDScore.scala:171`.
pub const CFD_MINIMUM_THRESHOLD: f64 = 0.023;

#[cfg(test)]
mod tests {
    use super::*;
    use crispr_encoding::{Base, Site};

    /// Build a 23-bp Cas9 NGG site from an ASCII string.
    fn site(ascii: &[u8]) -> Site {
        assert_eq!(ascii.len(), 23, "Cas9 NGG site must be 23 bp");
        Site::encode_ascii(ascii)
    }

    #[test]
    fn cfd_is_constructible() {
        let _ = Cfd::new();
    }

    #[test]
    fn lookup_table_has_no_zero_holes_for_matches() {
        let cfd = Cfd::new();
        // Self-comparison (guide == ot) implies all penalty[g][g][p] = 1.0
        // so the protospacer product is 1.0. Verified by `self_score_is_pam_weight`.
        for g in 0..4 {
            for p in 0..20 {
                assert!(
                    (cfd.penalty[g][g][p] - 1.0).abs() < f64::EPSILON,
                    "match penalty at g={g} p={p} should be 1.0, got {}",
                    cfd.penalty[g][g][p]
                );
            }
        }
    }

    #[test]
    fn self_score_with_ngg_pam_is_one() {
        // A guide compared to itself with NGG PAM → CFD = 1.0 × pam_weight["GG"]
        // = 1.0 × 1.0 = 1.0.
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let score = cfd.score_pair(guide, guide);
        assert!((score - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pam_weight_zero_for_non_canonical_pam() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        // OT has identical protospacer but PAM AAA → pam_weight["AA"] = 0.0
        // → total CFD = 0.0.
        let ot = site(b"AAAAAAAAAAAAAAAAAAAAAAA");
        let score = cfd.score_pair(guide, ot);
        assert!(score.abs() < 1e-12, "AA PAM scores 0.0, got {score}");
    }

    #[test]
    fn pam_tg_weight_matches_flashfry() {
        // PAM TGG → pam_weight["GG"] = 1.0 (the TG is the *first two*; we
        // use the *last two*, which for TGG is "GG"). Confusingly,
        // pam_weight["TG"] = 0.039 is what you'd get for a PAM ending in
        // "...TG" — not relevant for NGG-flanked OTs. Verify GG path:
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let ot = site(b"AAAAAAAAAAAAAAAAAAAATGG");
        let score = cfd.score_pair(guide, ot);
        assert!((score - 1.0).abs() < 1e-12);
    }

    #[test]
    fn single_mismatch_uses_table_value() {
        // Build a guide ending in "...AAAAAGG" (positions 14-19 are A's,
        // PAM AGG). OT differs at protospacer position 19 (the PAM-
        // proximal base) — guide T, OT A. Expected lookup: key "rU:dT,20"
        // → mm_lookup table value.
        //
        // Wait, position 19 is the *20th* protospacer base (1-indexed 20).
        // Guide: ...AAAAAT  AGG → 20-mer ends in T at 0-indexed pos 19.
        // OT:    ...AAAAAA  AGG → ends in A at 0-indexed pos 19.
        // Key: rU (guide T → rU), dT (OT A → specialReverseCompBase(A) = T → dT), pos 20.
        // Expected: mm_lookup["rU:dT,20"] = (let's check the table).
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAATAGG"); // pos 19 = T
        let ot = site(b"AAAAAAAAAAAAAAAAAAAAAGG"); // pos 19 = A
        let score = cfd.score_pair(guide, ot);

        // The expected penalty is the value at key "rU:dT,20".
        let expected = RAW_MM_PENALTIES
            .iter()
            .find(|&&(k, _)| k == "rU:dT,20")
            .expect("table contains rU:dT,20")
            .1;
        assert!(
            (score - expected).abs() < 1e-12,
            "expected {expected}, got {score}",
        );
    }

    #[test]
    fn two_mismatches_multiply_their_penalties() {
        let cfd = Cfd::new();
        // Two mismatches at positions 19 and 20:
        // - Pos 19 (0-indexed): guide T, OT A → rU:dT,20
        // - Pos 18 (0-indexed): guide T, OT A → rU:dT,19
        let guide = site(b"AAAAAAAAAAAAAAAAAATTAGG"); // pos 18,19 = TT
        let ot = site(b"AAAAAAAAAAAAAAAAAAAAAGG"); // pos 18,19 = AA
        let score = cfd.score_pair(guide, ot);

        let p1 = RAW_MM_PENALTIES
            .iter()
            .find(|&&(k, _)| k == "rU:dT,20")
            .unwrap()
            .1;
        let p2 = RAW_MM_PENALTIES
            .iter()
            .find(|&&(k, _)| k == "rU:dT,19")
            .unwrap()
            .1;
        let expected = p1 * p2;
        assert!(
            (score - expected).abs() < 1e-12,
            "expected {expected}, got {score}",
        );
    }

    #[test]
    fn parse_mm_key_handles_all_corners() {
        // First/last position, every (r, d) combination is exercised by
        // the construction. Sanity-check a handful by hand.
        assert_eq!(parse_mm_key("rA:dT,1"), (0, 0, 0)); // guide A, ot A
        assert_eq!(parse_mm_key("rU:dA,20"), (3, 3, 19)); // guide T, ot T
        assert_eq!(parse_mm_key("rC:dG,5"), (1, 1, 4)); // guide C, ot C
        assert_eq!(parse_mm_key("rG:dC,10"), (2, 2, 9)); // guide G, ot G
    }

    #[test]
    fn parse_pam_key_corners() {
        assert_eq!(parse_pam_key("GG"), (2, 2));
        assert_eq!(parse_pam_key("AA"), (0, 0));
        assert_eq!(parse_pam_key("TT"), (3, 3));
        assert_eq!(parse_pam_key("CG"), (1, 2));
    }

    /// Strict equality is appropriate here even with clippy's float-cmp
    /// caveat: these tests check *exact* table values or *exact*
    /// arithmetic identities (e.g. self-score = 1.0 × 1.0), not the
    /// outputs of any computation that could introduce rounding.
    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-12,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn aggregate_empty_yields_spec_1() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let result = cfd.aggregate(guide, std::iter::empty());
        assert_close(result.max_cfd, 0.0);
        assert_close(result.specificity, 1.0);
    }

    #[test]
    fn aggregate_skips_on_target() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        // One "off-target" hit that is actually the on-target (mm = 0).
        let result = cfd.aggregate(guide, [(guide, 0_u8, 1)]);
        assert_close(result.max_cfd, 0.0);
        assert_close(result.specificity, 1.0);
    }

    #[test]
    fn aggregate_specificity_formula() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        // One off-target with a single mismatch at the PAM-proximal base.
        let ot1 = site(b"AAAAAAAAAAAAAAAAAAATAGG"); // 1 mm vs guide
        let score1 = cfd.score_pair(guide, ot1);
        let result = cfd.aggregate(guide, [(ot1, 1_u8, 1)]);
        assert!((result.max_cfd - score1).abs() < 1e-12);
        let expected_spec = 1.0 / (1.0 + score1 * 1.0);
        assert!(
            (result.specificity - expected_spec).abs() < 1e-12,
            "expected {expected_spec}, got {}",
            result.specificity
        );
    }

    #[test]
    fn aggregate_weights_by_multi_mapping_count() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let ot = site(b"AAAAAAAAAAAAAAAAAAATAGG"); // 1 mm at pos 19
        let cfd_one = cfd.score_pair(guide, ot);
        // Same OT counted 5×: specificity sum should be 5 × cfd_one.
        let single = cfd.aggregate(guide, [(ot, 1_u8, 1)]);
        let multi = cfd.aggregate(guide, [(ot, 1_u8, 5)]);
        assert!((single.max_cfd - multi.max_cfd).abs() < 1e-12);
        let expected_multi_spec = 1.0 / (1.0 + cfd_one * 5.0);
        assert!(
            (multi.specificity - expected_multi_spec).abs() < 1e-12,
            "multi-map spec drift: expected {expected_multi_spec}, got {}",
            multi.specificity
        );
    }

    #[test]
    fn aggregate_max_picks_largest() {
        let cfd = Cfd::new();
        let guide = site(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let ot_bad = site(b"AAAAAAAAAAAAAAAAAAATAGG"); // 1 mm at pos 19 → high CFD
        let ot_worse = site(b"TTTTTTTTTTTTTTTTTTTTAGG"); // many mismatches → low CFD
        let s_bad = cfd.score_pair(guide, ot_bad);
        let s_worse = cfd.score_pair(guide, ot_worse);
        assert!(s_bad > s_worse, "single-mm should outscore all-mm");
        let result = cfd.aggregate(guide, [(ot_bad, 1_u8, 1), (ot_worse, 20_u8, 1)]);
        assert!((result.max_cfd - s_bad).abs() < 1e-12);
    }

    #[test]
    fn raw_mm_penalties_table_has_expected_size() {
        // FlashFry's table has 240 mismatch entries: 4 sgRNA bases × 3
        // non-matching DNA bases × 20 positions.
        assert_eq!(RAW_MM_PENALTIES.len(), 240);
    }

    #[test]
    fn raw_pam_weights_covers_all_16_pairs() {
        assert_eq!(RAW_PAM_WEIGHTS.len(), 16);
    }

    #[test]
    fn protospacer_iteration_order_matches_flashfry() {
        // FlashFry uses position 1 for the 5'-most protospacer base
        // (PAM-distal). Verify our bit-offset arithmetic matches: a
        // single mismatch at position 1 should use the key "...:...,1".
        let cfd = Cfd::new();
        // Guide first base = T; OT first base = A; rest identical.
        let guide = site(b"TAAAAAAAAAAAAAAAAAAAAGG"); // pos 0 = T
        let ot = site(b"AAAAAAAAAAAAAAAAAAAAAGG"); // pos 0 = A
        let score = cfd.score_pair(guide, ot);
        let expected = RAW_MM_PENALTIES
            .iter()
            .find(|&&(k, _)| k == "rU:dT,1")
            .unwrap()
            .1;
        assert!(
            (score - expected).abs() < 1e-12,
            "expected {expected}, got {score}",
        );
    }

    // Tying the construction to base codes — a sanity check the enum
    // numbering hasn't drifted out from under us.
    #[test]
    fn base_codes_match_table_indexing() {
        assert_eq!(Base::A.bits(), 0);
        assert_eq!(Base::C.bits(), 1);
        assert_eq!(Base::G.bits(), 2);
        assert_eq!(Base::T.bits(), 3);
    }
}
