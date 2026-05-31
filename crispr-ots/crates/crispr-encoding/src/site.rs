//! Packed CRISPR site: 2 × `u64` holding up to 28 bases plus a `u16` count.
//!
//! # Layout
//!
//! - `low` (u64): bases at positions `(n-24)..(n-1)` of the input string,
//!   two bits per base. The *last* character of the input string sits at
//!   bits 0-1 (LSB); the 24th-from-last at bits 46-47. Bits 48-63 are
//!   unused.
//! - `high` (u64):
//!   - bits 0-7: bases at positions `(n-28)..(n-25)` — zero unless `n > 24`
//!   - bits 8-15: reserved
//!   - bits 16-31: `u16` count (how many genomic positions share this seq)
//!   - bits 32-63: reserved for future metadata
//!
//! The "last char at LSB" convention is inherited from FlashFry so that the
//! PAM, which is always at one end of the string, ends up in a contiguous
//! bit range that can be masked with a single literal — and so that any
//! FlashFry-derived constants (lookup keys, mask values) transcribe
//! unchanged.

use serde::{Deserialize, Serialize};

use crate::base::Base;

/// Bit-packed CRISPR site. See module docs for layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[repr(C)]
pub struct Site {
    pub low: u64,
    pub high: u64,
}

/// Bitmask over a `Site`. Set bits identify 2-bit positions that should be
/// *included* when computing mismatches. PAM bits are typically cleared so
/// PAM differences (which the FASTA scanner already verified literally)
/// don't count as off-target mismatches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[repr(C)]
pub struct SiteMask {
    pub low: u64,
    pub high: u64,
}

/// Maximum bases a `Site` can hold for comparison. 24 in `low` + 4 in `high`.
pub const MAX_SITE_LEN: usize = 28;

/// 48 bits in `low` are available for sequence (24 bases × 2).
pub const LOW_SEQ_BITS: u32 = 48;
/// 8 bits in `high` are available for sequence (4 bases × 2).
pub const HIGH_SEQ_BITS: u32 = 8;
/// Bit offset of the `u16` count within `high`.
pub const COUNT_SHIFT: u32 = 16;
/// Mask isolating the count bits in `high`.
pub const COUNT_MASK: u64 = 0xFFFF_0000;
/// Mask covering all sequence-bearing bits in `low` (24 bases).
pub const LOW_SEQ_MASK: u64 = 0x0000_FFFF_FFFF_FFFF;
/// Mask covering all sequence-bearing bits in `high` (4 bases).
pub const HIGH_SEQ_MASK: u64 = 0x0000_0000_0000_00FF;

impl Site {
    /// Encode an iterator of `Base`s.
    ///
    /// # Panics
    /// Panics if the input yields more than [`MAX_SITE_LEN`] bases.
    pub fn encode<I: IntoIterator<Item = Base>>(bases: I) -> Self {
        let collected: Vec<Base> = bases.into_iter().collect();
        Self::encode_slice(&collected)
    }

    /// Encode a slice of `Base`s.
    ///
    /// # Panics
    /// Panics if `bases.len() > MAX_SITE_LEN`.
    #[must_use]
    pub fn encode_slice(bases: &[Base]) -> Self {
        assert!(
            bases.len() <= MAX_SITE_LEN,
            "site sequence of {} bases exceeds the {MAX_SITE_LEN}-base encoding limit",
            bases.len()
        );
        let n = bases.len();
        if n <= 24 {
            let mut low: u64 = 0;
            for b in bases {
                low = (low << 2) | u64::from(b.bits());
            }
            Self { low, high: 0 }
        } else {
            let n_high = n - 24;
            let mut high: u64 = 0;
            for b in &bases[..n_high] {
                high = (high << 2) | u64::from(b.bits());
            }
            let mut low: u64 = 0;
            for b in &bases[n_high..] {
                low = (low << 2) | u64::from(b.bits());
            }
            Self { low, high }
        }
    }

    /// Encode an ASCII sequence.
    ///
    /// # Panics
    /// Panics if any byte is not in `{A, C, G, T}` (case-insensitive) or if
    /// the sequence exceeds [`MAX_SITE_LEN`].
    #[must_use]
    pub fn encode_ascii(seq: &[u8]) -> Self {
        let bases: Vec<Base> = seq
            .iter()
            .map(|&c| {
                Base::from_ascii(c)
                    .unwrap_or_else(|| panic!("non-ACGT byte {c:#x} in encode_ascii"))
            })
            .collect();
        Self::encode_slice(&bases)
    }

    /// Read the count field (bits 16-31 of `high`).
    #[inline]
    #[must_use]
    pub const fn count(self) -> u16 {
        ((self.high & COUNT_MASK) >> COUNT_SHIFT) as u16
    }

    /// Return a copy of this site with the count field replaced.
    #[inline]
    #[must_use]
    pub const fn with_count(mut self, count: u16) -> Self {
        self.high = (self.high & !COUNT_MASK) | ((count as u64) << COUNT_SHIFT);
        self
    }

    /// Decode `n` bases. Bits beyond position `n` are ignored.
    ///
    /// # Panics
    /// Panics if `n > MAX_SITE_LEN`.
    #[must_use]
    pub fn decode(self, n: usize) -> Vec<Base> {
        assert!(
            n <= MAX_SITE_LEN,
            "decode length {n} exceeds {MAX_SITE_LEN}"
        );
        let mut out = Vec::with_capacity(n);
        let n_high = n.saturating_sub(24);
        // Iterate from highest bit position (first char) to lowest (last
        // char). Within `high`, the first char of the input is at bit
        // position `(n_high - 1) * 2`.
        for i in (0..n_high).rev() {
            let bits = ((self.high >> (i * 2)) & 0b11) as u8;
            out.push(Base::from_bits(bits).expect("valid 2-bit base"));
        }
        let low_count = n - n_high;
        for i in (0..low_count).rev() {
            let bits = ((self.low >> (i * 2)) & 0b11) as u8;
            out.push(Base::from_bits(bits).expect("valid 2-bit base"));
        }
        out
    }

    /// Decode to an ASCII string of length `n`.
    #[must_use]
    pub fn decode_ascii(self, n: usize) -> String {
        self.decode(n)
            .into_iter()
            .map(|b| b.to_ascii() as char)
            .collect()
    }
}

impl SiteMask {
    /// Empty mask (no bits set).
    pub const EMPTY: Self = Self { low: 0, high: 0 };

    /// A mask that covers exactly the first `n` bases of an encoded site.
    /// Use as a base mask; clear PAM bits separately if desired.
    ///
    /// # Panics
    /// Panics if `n > MAX_SITE_LEN`.
    #[must_use]
    pub fn for_len(n: usize) -> Self {
        assert!(n <= MAX_SITE_LEN, "mask length {n} exceeds {MAX_SITE_LEN}");
        let low_bits = n.min(24) * 2;
        // Safe: `n.min(24) * 2` is at most 48, so the shift is in range.
        let low = if low_bits == 0 {
            0
        } else {
            (1u64 << low_bits) - 1
        };
        let high = if n > 24 {
            (1u64 << ((n - 24) * 2)) - 1
        } else {
            0
        };
        Self { low, high }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_site_is_default() {
        let s = Site::encode_slice(&[]);
        assert_eq!(s.low, 0);
        assert_eq!(s.high, 0);
        assert_eq!(s, Site::default());
    }

    #[test]
    fn last_char_at_lsb() {
        let site = Site::encode_ascii(b"ACGT");
        // String "ACGT": last char T → bits 0-1 of low.
        assert_eq!(site.low & 0b11, u64::from(Base::T.bits()));
        assert_eq!((site.low >> 2) & 0b11, u64::from(Base::G.bits()));
        assert_eq!((site.low >> 4) & 0b11, u64::from(Base::C.bits()));
        assert_eq!((site.low >> 6) & 0b11, u64::from(Base::A.bits()));
        assert_eq!(site.high, 0);
    }

    #[test]
    fn flashfry_aagct_layout() {
        // FlashFry's BitEncoding encodes "AAGCT" as:
        //   start at 0
        //   A: (0 << 2) | 0 = 0
        //   A: (0 << 2) | 0 = 0
        //   G: (0 << 2) | 2 = 0b10
        //   C: (0b10 << 2) | 1 = 0b1001
        //   T: (0b1001 << 2) | 3 = 0b100111 = 0x27
        let site = Site::encode_ascii(b"AAGCT");
        assert_eq!(site.low, 0x27);
        assert_eq!(site.high, 0);
    }

    #[test]
    fn roundtrip_short() {
        let s = "ACGTAC";
        let site = Site::encode_ascii(s.as_bytes());
        assert_eq!(site.decode_ascii(s.len()), s);
    }

    #[test]
    fn roundtrip_at_capacity_24() {
        let s = "ACGTACGTACGTACGTACGTACGT"; // 24 chars
        assert_eq!(s.len(), 24);
        let site = Site::encode_ascii(s.as_bytes());
        assert_eq!(site.high & HIGH_SEQ_MASK, 0, "no overflow into high");
        assert_eq!(site.decode_ascii(s.len()), s);
    }

    #[test]
    fn roundtrip_beyond_24_uses_high() {
        let s = "GTACACGTACGTACGTACGTACGTACGT"; // 28 chars
        assert_eq!(s.len(), 28);
        let site = Site::encode_ascii(s.as_bytes());
        assert_ne!(
            site.high & HIGH_SEQ_MASK,
            0,
            "first 4 chars spill into high"
        );
        assert_eq!(site.decode_ascii(s.len()), s);
    }

    #[test]
    fn roundtrip_27_chars() {
        // Cas12a TTTV (4 bp) + 23-nt protospacer = 27 chars — our motivating
        // case for the 2 × u64 layout.
        let s = "TTTAACGTACGTACGTACGTACGTACG"; // 27 chars
        assert_eq!(s.len(), 27);
        let site = Site::encode_ascii(s.as_bytes());
        assert_eq!(site.decode_ascii(s.len()), s);
    }

    #[test]
    #[should_panic(expected = "exceeds the 28-base encoding limit")]
    fn encoding_29_panics() {
        let s = vec![b'A'; 29];
        let _ = Site::encode_ascii(&s);
    }

    #[test]
    #[should_panic(expected = "non-ACGT byte")]
    fn encoding_n_panics() {
        let _ = Site::encode_ascii(b"ACGTN");
    }

    #[test]
    fn count_roundtrip() {
        let s = Site::default().with_count(12_345);
        assert_eq!(s.count(), 12_345);

        // u16::MAX fits.
        let s = Site::default().with_count(u16::MAX);
        assert_eq!(s.count(), u16::MAX);

        // Setting the count does not perturb the sequence bits.
        let site = Site::encode_ascii(b"ACGT").with_count(99);
        assert_eq!(site.decode_ascii(4), "ACGT");
        assert_eq!(site.count(), 99);
    }

    #[test]
    fn count_overwrites_cleanly() {
        let s = Site::default().with_count(0xABCD).with_count(0x1234);
        assert_eq!(s.count(), 0x1234);
        // No bleed into sequence bits.
        assert_eq!(s.low, 0);
        assert_eq!(s.high & !COUNT_MASK, 0);
    }

    #[test]
    fn mask_for_len_20() {
        let m = SiteMask::for_len(20);
        assert_eq!(m.low, (1u64 << 40) - 1);
        assert_eq!(m.high, 0);
    }

    #[test]
    fn mask_for_len_23() {
        // 23 bases × 2 bits = 46 bits.
        let m = SiteMask::for_len(23);
        assert_eq!(m.low, 0x3FFF_FFFF_FFFF);
        assert_eq!(m.high, 0);
    }

    #[test]
    fn mask_for_len_24_is_full_low() {
        let m = SiteMask::for_len(24);
        assert_eq!(m.low, LOW_SEQ_MASK);
        assert_eq!(m.high, 0);
    }

    #[test]
    fn mask_for_len_27_uses_high() {
        let m = SiteMask::for_len(27);
        assert_eq!(m.low, LOW_SEQ_MASK);
        // 27 - 24 = 3 bases in high → 6 bits.
        assert_eq!(m.high, (1u64 << 6) - 1);
    }

    #[test]
    fn mask_for_len_28_full() {
        let m = SiteMask::for_len(28);
        assert_eq!(m.low, LOW_SEQ_MASK);
        assert_eq!(m.high, HIGH_SEQ_MASK);
    }

    #[test]
    fn mask_for_len_0() {
        let m = SiteMask::for_len(0);
        assert_eq!(m, SiteMask::EMPTY);
    }
}
