//! Enzyme configuration: PAM pattern, protospacer length, comparison mask.
//!
//! Data-driven by design (unlike FlashFry's hard-coded `case object`
//! parameter packs). IUPAC PAM codes are expanded to literal `Base`
//! sequences at construction time, so the FASTA scanner and the bit-level
//! mismatch loop only ever deal with concrete A/C/G/T.
//!
//! Compare-masks are derived from `pam_side`, `pam_len`, and
//! `protospacer_len` — adding a new enzyme is a one-line constructor call,
//! not a code change to the encoding layer.

use serde::{Deserialize, Serialize};

use crate::base::{Base, IupacCode};
use crate::site::{SiteMask, LOW_SEQ_MASK, MAX_SITE_LEN};

/// Which end of the protospacer the PAM sits on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PamSide {
    /// PAM follows the protospacer (3' end). SpCas9-style. The PAM ends up
    /// in the *lowest* bits of the encoded `Site` (since the string's last
    /// characters are at the LSB).
    ThreePrime,
    /// PAM precedes the protospacer (5' end). Cpf1/Cas12a-style. The PAM
    /// ends up in the *highest* used bits of the encoded `Site`.
    FivePrime,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Enzyme {
    pub name: String,
    pub pam_side: PamSide,
    pub pam_len: u8,
    pub protospacer_len: u8,
    /// Bits set on protospacer positions; PAM positions cleared. Used as
    /// the `mask` argument to [`crate::mismatches_masked`].
    pub compare_mask: SiteMask,
    /// All literal PAM sequences the enzyme accepts. Built by expanding the
    /// configured IUPAC pattern (or by direct construction for multi-pattern
    /// enzymes like SpCas9 NGG+NAG).
    pub canonical_pams: Vec<Vec<Base>>,
}

impl Enzyme {
    /// Construct from an IUPAC PAM pattern.
    ///
    /// # Panics
    /// Panics if `pam.len() + protospacer_len > MAX_SITE_LEN` or if the
    /// total exceeds `u8::MAX`.
    #[must_use]
    pub fn from_iupac(
        name: impl Into<String>,
        pam: &[IupacCode],
        pam_side: PamSide,
        protospacer_len: u8,
    ) -> Self {
        let pam_len = u8::try_from(pam.len()).expect("PAM length must fit in u8");
        let canonical = expand_iupac(pam);
        Self::from_canonical(name, canonical, pam_side, pam_len, protospacer_len)
    }

    /// Construct directly from pre-expanded literal PAMs. Use this for
    /// enzymes whose accepted PAM set isn't a single IUPAC pattern (e.g.
    /// SpCas9 NGG+NAG = `{AGG, CGG, GGG, TGG, AAG, CAG, GAG, TAG}`).
    ///
    /// # Panics
    /// Panics if `canonical_pams` is empty, if any entry's length differs
    /// from `pam_len`, or if `pam_len + protospacer_len > MAX_SITE_LEN`.
    #[must_use]
    pub fn from_canonical(
        name: impl Into<String>,
        canonical_pams: Vec<Vec<Base>>,
        pam_side: PamSide,
        pam_len: u8,
        protospacer_len: u8,
    ) -> Self {
        assert!(
            !canonical_pams.is_empty(),
            "enzyme must have at least one canonical PAM"
        );
        let total = usize::from(pam_len) + usize::from(protospacer_len);
        assert!(
            total <= MAX_SITE_LEN,
            "PAM ({pam_len}) + protospacer ({protospacer_len}) = {total} \
             exceeds the {MAX_SITE_LEN}-base site limit"
        );
        for pam in &canonical_pams {
            assert_eq!(
                pam.len(),
                usize::from(pam_len),
                "canonical PAM length mismatch (expected {pam_len}, got {})",
                pam.len()
            );
        }
        let compare_mask = build_compare_mask(pam_side, pam_len, protospacer_len);
        Self {
            name: name.into(),
            pam_side,
            pam_len,
            protospacer_len,
            compare_mask,
            canonical_pams,
        }
    }

    /// Total length (PAM + protospacer) of a site for this enzyme.
    #[inline]
    #[must_use]
    pub const fn total_scan_len(&self) -> u8 {
        self.pam_len + self.protospacer_len
    }

    // ---- FlashFry-equivalent presets ----

    /// SpCas9 with strict NGG PAM, 20-nt protospacer.
    /// Equivalent to FlashFry's `Cas9NGGParameterPack`.
    #[must_use]
    pub fn spcas9_ngg() -> Self {
        Self::from_iupac(
            "SpCas9-NGG",
            &[IupacCode::N, IupacCode::G, IupacCode::G],
            PamSide::ThreePrime,
            20,
        )
    }

    /// SpCas9 with strict NAG PAM, 20-nt protospacer.
    /// Equivalent to FlashFry's `Cas9NAGParameterPack`.
    #[must_use]
    pub fn spcas9_nag() -> Self {
        Self::from_iupac(
            "SpCas9-NAG",
            &[IupacCode::N, IupacCode::A, IupacCode::G],
            PamSide::ThreePrime,
            20,
        )
    }

    /// SpCas9 accepting either NGG or NAG, 20-nt protospacer.
    /// Equivalent to FlashFry's `Cas9ParameterPack` (`pam = ["GG", "AG"]`).
    #[must_use]
    pub fn spcas9_ngg_or_nag() -> Self {
        let mut canonical = expand_iupac(&[IupacCode::N, IupacCode::G, IupacCode::G]);
        canonical.extend(expand_iupac(&[IupacCode::N, IupacCode::A, IupacCode::G]));
        Self::from_canonical("SpCas9-NGG-or-NAG", canonical, PamSide::ThreePrime, 3, 20)
    }

    /// Cpf1 (Cas12a) with TTTN PAM, **23-nt protospacer** (total scan
    /// length 27 bp).
    ///
    /// We diverge from FlashFry's `Cpf1ParameterPack` (which caps at 20-nt
    /// because its 24-bp `Long` encoding can't fit anything longer) and
    /// match the protospacer length the published Cas12a activity matrices
    /// were measured at (`crispr-score::Cas12aMatrix`'s 23 positions).
    /// Our `MAX_SITE_LEN = 28` supports the 27-bp encoded site directly,
    /// so all 23 matrix positions enter the CFD-like score.
    ///
    /// Note FlashFry treats the variable `N` as part of a 4-bp PAM region,
    /// *not* as the first base of the protospacer — we follow that
    /// convention.
    #[must_use]
    pub fn cpf1_tttn() -> Self {
        Self::from_iupac(
            "Cpf1-TTTN",
            &[IupacCode::T, IupacCode::T, IupacCode::T, IupacCode::N],
            PamSide::FivePrime,
            23,
        )
    }
}

/// Expand an IUPAC pattern to every literal `Vec<Base>` it matches.
///
/// `&[N]` → 4 entries, `&[N, G, G]` → 4 entries, `&[T, T, T, V]` → 3 entries,
/// `&[]` → one empty `Vec` (identity element of the cartesian product).
#[must_use]
pub fn expand_iupac(pattern: &[IupacCode]) -> Vec<Vec<Base>> {
    pattern.iter().fold(vec![Vec::new()], |acc, code| {
        let bases = code.bases();
        let mut next = Vec::with_capacity(acc.len() * bases.len());
        for prefix in &acc {
            for &b in bases {
                let mut extended = prefix.clone();
                extended.push(b);
                next.push(extended);
            }
        }
        next
    })
}

/// Build the comparison mask: protospacer bits set, PAM bits cleared.
///
/// Encoding convention recap (see `site.rs` docs):
/// - The *last* character of the input string is at bits 0-1 of `low`.
/// - For 3' PAM (Cas9-style), the PAM occupies the lowest `pam_len * 2` bits.
/// - For 5' PAM (Cpf1-style), the PAM occupies the highest bits of the used
///   range (positions `0..pam_len` of the string → bits at the top of `low`,
///   and into `high` for total length > 24).
fn build_compare_mask(side: PamSide, pam_len: u8, protospacer_len: u8) -> SiteMask {
    let total = usize::from(pam_len) + usize::from(protospacer_len);
    let full = SiteMask::for_len(total);
    let pam_bits = u32::from(pam_len) * 2;

    match side {
        PamSide::ThreePrime => {
            // PAM at low end. Clear bits 0..pam_bits.
            let pam_mask = (1u64 << pam_bits) - 1;
            SiteMask {
                low: full.low & !pam_mask,
                high: full.high,
            }
        }
        PamSide::FivePrime => {
            let proto_bits = u32::from(protospacer_len) * 2;
            if total <= 24 {
                // PAM occupies bits proto_bits..(total*2) of `low`. `high`
                // is zero in `full` so leave it alone.
                let proto_mask = (1u64 << proto_bits) - 1;
                SiteMask {
                    low: full.low & proto_mask,
                    high: 0,
                }
            } else if usize::from(protospacer_len) <= 24 {
                // Protospacer fits entirely in `low`; PAM occupies the
                // remaining high bits of `low` plus all of `high`.
                let proto_mask = (1u64 << proto_bits) - 1;
                SiteMask {
                    low: full.low & proto_mask,
                    high: 0,
                }
            } else {
                // Protospacer spans low + part of high; PAM only in high.
                let proto_in_high = u32::from(protospacer_len) - 24;
                let high_keep = (1u64 << (proto_in_high * 2)) - 1;
                SiteMask {
                    low: full.low & LOW_SEQ_MASK,
                    high: full.high & high_keep,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expand_n_yields_four() {
        let pams = expand_iupac(&[IupacCode::N]);
        assert_eq!(pams.len(), 4);
        for b in Base::ALL {
            assert!(pams.contains(&vec![b]));
        }
    }

    #[test]
    fn expand_ngg() {
        let pams = expand_iupac(&[IupacCode::N, IupacCode::G, IupacCode::G]);
        assert_eq!(pams.len(), 4);
        for pam in &pams {
            assert_eq!(pam.len(), 3);
            assert_eq!(pam[1], Base::G);
            assert_eq!(pam[2], Base::G);
        }
    }

    #[test]
    fn expand_tttv() {
        let pams = expand_iupac(&[IupacCode::T, IupacCode::T, IupacCode::T, IupacCode::V]);
        assert_eq!(pams.len(), 3, "V expands to A, C, G (not T)");
        for pam in &pams {
            assert_eq!(pam[0], Base::T);
            assert_eq!(pam[1], Base::T);
            assert_eq!(pam[2], Base::T);
            assert_ne!(pam[3], Base::T);
        }
    }

    #[test]
    fn expand_empty_is_identity() {
        let pams = expand_iupac(&[]);
        assert_eq!(pams, vec![Vec::<Base>::new()]);
    }

    #[test]
    fn expand_two_ns() {
        let pams = expand_iupac(&[IupacCode::N, IupacCode::N]);
        assert_eq!(pams.len(), 16);
    }

    // ---- FlashFry mask-constant equivalence ----

    #[test]
    fn spcas9_ngg_mask_matches_flashfry() {
        let e = Enzyme::spcas9_ngg();
        // FlashFry: Cas9NGGParameterPack.comparisonBitEncoding = 0x3FFFFFFFFFC0L
        assert_eq!(e.compare_mask.low, 0x3FFF_FFFF_FFC0);
        assert_eq!(e.compare_mask.high, 0);
        assert_eq!(e.total_scan_len(), 23);
        assert_eq!(e.canonical_pams.len(), 4);
    }

    #[test]
    fn spcas9_nag_mask_matches_flashfry() {
        let e = Enzyme::spcas9_nag();
        assert_eq!(e.compare_mask.low, 0x3FFF_FFFF_FFC0);
        assert_eq!(e.total_scan_len(), 23);
        assert_eq!(e.canonical_pams.len(), 4);
        for pam in &e.canonical_pams {
            assert_eq!(pam[1], Base::A);
            assert_eq!(pam[2], Base::G);
        }
    }

    #[test]
    fn spcas9_ngg_or_nag_mask_matches_flashfry() {
        let e = Enzyme::spcas9_ngg_or_nag();
        assert_eq!(e.compare_mask.low, 0x3FFF_FFFF_FFC0);
        assert_eq!(e.canonical_pams.len(), 8);
    }

    #[test]
    fn cpf1_tttn_23nt_protospacer_layout() {
        // We diverge from FlashFry's 20-nt cap (chosen so the encoded
        // site fits in a single Long): we run 23-nt protospacer for full
        // matrix coverage. PAM occupies bits 46-47 of `low` plus bits 0-5
        // of `high`; the protospacer fills bits 0-45 of `low`.
        let e = Enzyme::cpf1_tttn();
        assert_eq!(e.total_scan_len(), 27);
        assert_eq!(e.compare_mask.low, (1u64 << 46) - 1);
        assert_eq!(e.compare_mask.high, 0);
        assert_eq!(e.canonical_pams.len(), 4);
    }

    // ---- Forward-looking flexibility ----

    #[test]
    fn cas12a_tttv_with_23nt_protospacer() {
        // The headline motivating case for 2 × u64 encoding: 4 + 23 = 27 bp.
        let e = Enzyme::from_iupac(
            "Cas12a-TTTV-23",
            &[IupacCode::T, IupacCode::T, IupacCode::T, IupacCode::V],
            PamSide::FivePrime,
            23,
        );
        assert_eq!(e.total_scan_len(), 27);
        // Protospacer (23 bases) occupies bits 0..46 of `low`.
        assert_eq!(e.compare_mask.low, (1u64 << 46) - 1);
        // PAM (4 bases) occupies the rest: bits 46-47 of `low` plus bits
        // 0-5 of `high`. The compare-mask zeros all of those.
        assert_eq!(e.compare_mask.high, 0);
        assert_eq!(e.canonical_pams.len(), 3);
    }

    #[test]
    fn cas9_19mer_mask_matches_flashfry() {
        let e = Enzyme::from_iupac(
            "SpCas9-NGG-19",
            &[IupacCode::N, IupacCode::G, IupacCode::G],
            PamSide::ThreePrime,
            19,
        );
        // FlashFry: Cas9NGG19ParameterPack.comparisonBitEncoding = 0x0FFFFFFFFFC0L
        assert_eq!(e.compare_mask.low, 0x0FFF_FFFF_FFC0);
    }

    #[test]
    #[should_panic(expected = "exceeds the 28-base site limit")]
    fn rejects_oversized_enzyme() {
        let _ = Enzyme::from_iupac(
            "Too-Long",
            &[
                IupacCode::N,
                IupacCode::N,
                IupacCode::N,
                IupacCode::N,
                IupacCode::N,
            ],
            PamSide::ThreePrime,
            24,
        );
    }

    #[test]
    #[should_panic(expected = "at least one canonical PAM")]
    fn rejects_empty_canonical_set() {
        let _ = Enzyme::from_canonical("Empty", Vec::new(), PamSide::ThreePrime, 3, 20);
    }
}
