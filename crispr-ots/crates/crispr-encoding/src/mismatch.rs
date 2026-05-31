//! Bit-parallel Hamming-mismatch counting between two `Site`s under a mask.
//!
//! The core algorithm is FlashFry's POPCNT trick
//! (`bitcoding/BitEncoding.scala:127-132`), extended over two `u64`s so it
//! works for sites longer than 24 bases:
//!
//! 1. `xor` the two encodings → bit pairs differ iff at least one of the
//!    two bits in the pair is set.
//! 2. AND with the comparison mask to zero out PAM / unused positions.
//! 3. For each 2-bit pair, set one "any-bit-set" indicator by OR-ing the
//!    pair's low and high bits into the same high position.
//! 4. `popcount` those indicators → number of differing bases.
//!
//! On any platform with hardware POPCNT (x86-64, aarch64, RISC-V with Zbb),
//! `u64::count_ones` compiles to a single instruction; the whole primitive
//! is a handful of ALU ops per `u64`.

use crate::site::{Site, SiteMask};

/// Bits at the "high" position of every 2-bit pair within a `u64`.
const UPPER_BITS: u64 = 0xAAAA_AAAA_AAAA_AAAA;

/// Count base-level mismatches between `a` and `b`, considering only the
/// 2-bit positions where `mask` is set.
///
/// PAM bits should already be cleared in `mask` by the caller — typically
/// via [`crate::enzyme::Enzyme::compare_mask`].
#[inline]
#[must_use]
pub fn mismatches_masked(a: Site, b: Site, mask: SiteMask) -> u32 {
    let cl = (a.low ^ b.low) & mask.low;
    let ch = (a.high ^ b.high) & mask.high;
    pop_pairs(cl) + pop_pairs(ch)
}

/// Population count of differing 2-bit pairs within a single `u64`.
///
/// For each pair `(b1, b0)`, sets a single indicator bit at the high
/// position iff either of `b1` or `b0` is set.
#[inline]
const fn pop_pairs(x: u64) -> u32 {
    ((x & UPPER_BITS) | ((x << 1) & UPPER_BITS)).count_ones()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::site::{Site, SiteMask, MAX_SITE_LEN};

    #[test]
    fn identical_sites_zero_mismatches() {
        let s = Site::encode_ascii(b"ACGTACGT");
        assert_eq!(mismatches_masked(s, s, SiteMask::for_len(8)), 0);
    }

    #[test]
    fn one_base_diff_one_mismatch() {
        let a = Site::encode_ascii(b"ACGTACGT");
        let b = Site::encode_ascii(b"ACGTACGA"); // T→A at last position
        assert_eq!(mismatches_masked(a, b, SiteMask::for_len(8)), 1);
    }

    #[test]
    fn opposite_strings_n_mismatches() {
        // A↔T differ in both bit positions of the pair (00 vs 11), which
        // exercises the OR-of-both-bits path.
        let a = Site::encode_ascii(b"AAAAAA");
        let b = Site::encode_ascii(b"TTTTTT");
        assert_eq!(mismatches_masked(a, b, SiteMask::for_len(6)), 6);
    }

    #[test]
    fn ag_vs_ct_two_mismatches_per_pair() {
        // A=00, G=10, C=01, T=11. A↔C differs in low bit (00 vs 01); G↔T
        // differs in low bit (10 vs 11). Each differs by exactly one bit
        // but should still count as one *base* mismatch.
        let a = Site::encode_ascii(b"AG");
        let b = Site::encode_ascii(b"CT");
        assert_eq!(mismatches_masked(a, b, SiteMask::for_len(2)), 2);
    }

    #[test]
    fn mask_can_isolate_a_single_position() {
        let a = Site::encode_ascii(b"ACGT");
        let b = Site::encode_ascii(b"AAAA"); // differs at positions 1, 2, 3

        // Mask covering only the last position (bits 0-1).
        let last_only = SiteMask { low: 0b11, high: 0 };
        assert_eq!(mismatches_masked(a, b, last_only), 1);

        // Mask covering only the first position (bits 6-7 of a 4-base site).
        let first_only = SiteMask {
            low: 0b11 << 6,
            high: 0,
        };
        assert_eq!(mismatches_masked(a, b, first_only), 0); // both start with A
    }

    #[test]
    fn flashfry_spcas9_mask_excludes_pam() {
        // 20 A's protospacer, two different PAMs (AGG vs TGG).
        let on = Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG");
        let off = Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAATGG");

        // FlashFry's Cas9ParameterPack.comparisonBitEncoding.
        let pam_excluded = SiteMask {
            low: 0x3FFF_FFFF_FFC0,
            high: 0,
        };
        assert_eq!(
            mismatches_masked(on, off, pam_excluded),
            0,
            "PAM-only difference must be masked out"
        );

        // Without the mask, the N-base difference is visible.
        assert_eq!(mismatches_masked(on, off, SiteMask::for_len(23)), 1);
    }

    #[test]
    fn flashfry_cpf1_mask_excludes_5prime_pam() {
        // 4-bp PAM (TTTN) + 20 A's protospacer. PAMs differ only at N.
        let on = Site::encode_ascii(b"TTTAAAAAAAAAAAAAAAAAAAAA");
        let off = Site::encode_ascii(b"TTTCAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(on.decode_ascii(24).len(), 24);

        // FlashFry's Cpf1ParameterPack.comparisonBitEncoding.
        let pam_excluded = SiteMask {
            low: 0x00FF_FFFF_FFFF,
            high: 0,
        };
        assert_eq!(
            mismatches_masked(on, off, pam_excluded),
            0,
            "5' PAM difference must be masked out"
        );
    }

    #[test]
    fn mismatch_count_uses_high_for_long_sites() {
        // 27-character sites (e.g., Cas12a TTTV + 23-nt protospacer).
        // Differ at position 0 only, which lives in `high`.
        let a = Site::encode_ascii(b"TTTAACGTACGTACGTACGTACGTACG"); // 27 chars, first = T
        let b = Site::encode_ascii(b"ATTAACGTACGTACGTACGTACGTACG"); // first = A
        assert_eq!(a.decode_ascii(27).len(), 27);

        let mask = SiteMask::for_len(27);
        assert_eq!(
            mismatches_masked(a, b, mask),
            1,
            "diff in `high` must be counted"
        );
    }

    #[test]
    fn empty_mask_yields_zero() {
        let a = Site::encode_ascii(b"AAAA");
        let b = Site::encode_ascii(b"TTTT");
        assert_eq!(mismatches_masked(a, b, SiteMask::EMPTY), 0);
    }

    #[test]
    fn maximum_length_28_all_different() {
        let a = Site::encode_ascii(&[b'A'; MAX_SITE_LEN]);
        let b = Site::encode_ascii(&[b'T'; MAX_SITE_LEN]);
        let mask = SiteMask::for_len(MAX_SITE_LEN);
        assert_eq!(mismatches_masked(a, b, mask), 28);
    }
}
