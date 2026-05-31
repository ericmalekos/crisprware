//! `PackedEntry`: 24-byte on-wire and in-memory record for one off-target
//! site.
//!
//! The internal `BinTable` storage uses this directly, and the new v2
//! on-disk format is a contiguous sequence of these per bin. The packed
//! layout lets us:
//!
//! - `mmap` a database file and treat the payload as `&[PackedEntry]`
//!   without any per-entry deserialization (via `bytemuck::cast_slice`).
//! - Shrink the in-memory `BinTable` modestly (24 bytes per site vs ~28
//!   for `(Site, Position)` with `Strand` enum padding under bincode).
//! - Keep `Site` byte-for-byte intact so the SIMD inner loop still loads
//!   `site_low` / `site_high` contiguously.
//!
//! The `Position` field is packed into the third `u64`:
//!
//! ```text
//! bit  0..30  offset            u31, 0..2^31      (human chr1 ≈ 2^28 < 2^31)
//! bit  31     strand            0 = +, 1 = -
//! bits 32..47 contig_id         u16, 0..65535     (genomes have ≪ 65K contigs)
//! bits 48..63 reserved          must be zero in v2
//! ```

use bytemuck::{Pod, Zeroable};
use crispr_encoding::Site;
use serde::{Deserialize, Serialize};

use crate::position::{Position, Strand};

/// One stored off-target site: bit-packed `Site` + bit-packed `Position`.
#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Pod, Zeroable, Serialize, Deserialize,
)]
#[repr(C)]
pub struct PackedEntry {
    /// Same encoding as `crispr_encoding::Site::low`.
    pub site_low: u64,
    /// Same encoding as `crispr_encoding::Site::high`.
    pub site_high: u64,
    /// Bit-packed `Position`; see module docs for the layout.
    pub packed_position: u64,
}

/// 31-bit mask isolating the `offset` field of `packed_position`.
const OFFSET_MASK: u64 = 0x7FFF_FFFF;
/// Bit position of the strand flag within `packed_position`.
const STRAND_BIT: u32 = 31;
/// Bit position of the `contig_id` field within `packed_position`.
const CONTIG_SHIFT: u32 = 32;
/// 16-bit mask for `contig_id` after shifting.
const CONTIG_MASK: u64 = 0xFFFF;
/// Largest `contig_id` that can be stored. Any genome with more than
/// this many contigs is unsupported by the v2 on-disk format.
pub const MAX_CONTIG_ID: u32 = 0xFFFF;
/// Largest offset that can be stored (u31 = ~2.1 G). Human chr1 is the
/// largest known mammalian chromosome at ~248 Mbp, comfortably below.
pub const MAX_OFFSET: u32 = 0x7FFF_FFFF;

impl PackedEntry {
    /// Build from in-memory `Site` and `Position`.
    ///
    /// # Panics
    /// Panics if `position.contig_id > MAX_CONTIG_ID` or
    /// `position.offset > MAX_OFFSET`. Both limits are well above any
    /// realistic genome; the CLI layer pre-validates inputs.
    #[must_use]
    pub fn from_site_position(site: Site, position: Position) -> Self {
        assert!(
            position.contig_id <= MAX_CONTIG_ID,
            "contig_id {} exceeds v2 limit {}",
            position.contig_id,
            MAX_CONTIG_ID
        );
        assert!(
            position.offset <= MAX_OFFSET,
            "offset {} exceeds v2 limit {} (~2.1 G)",
            position.offset,
            MAX_OFFSET
        );
        let strand_bit = match position.strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        };
        let packed_position = u64::from(position.offset)
            | (strand_bit << STRAND_BIT)
            | (u64::from(position.contig_id) << CONTIG_SHIFT);
        Self {
            site_low: site.low,
            site_high: site.high,
            packed_position,
        }
    }

    /// Recover the `Site` half. Free — just reads two `u64` fields.
    #[inline]
    #[must_use]
    pub const fn site(self) -> Site {
        Site {
            low: self.site_low,
            high: self.site_high,
        }
    }

    /// Recover the `Position` half. ~3 ALU ops; no allocations.
    #[inline]
    #[must_use]
    pub const fn position(self) -> Position {
        let p = self.packed_position;
        let offset = (p & OFFSET_MASK) as u32;
        let strand = if (p >> STRAND_BIT) & 1 == 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        // Cast to u16 first to clamp, then widen back to u32.
        let contig_id = ((p >> CONTIG_SHIFT) & CONTIG_MASK) as u32;
        Position {
            contig_id,
            offset,
            strand,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(low: u64, high: u64) -> Site {
        Site { low, high }
    }

    #[test]
    fn roundtrip_forward() {
        let s = site(0xDEAD_BEEF_FEED_FACE, 0x1234_5678);
        let p = Position {
            contig_id: 17,
            offset: 12_345_678,
            strand: Strand::Forward,
        };
        let e = PackedEntry::from_site_position(s, p);
        assert_eq!(e.site(), s);
        assert_eq!(e.position(), p);
    }

    #[test]
    fn roundtrip_reverse() {
        let s = site(0x42, 0);
        let p = Position {
            contig_id: 0,
            offset: 0,
            strand: Strand::Reverse,
        };
        let e = PackedEntry::from_site_position(s, p);
        assert_eq!(e.position(), p);
    }

    #[test]
    fn max_offset_roundtrips() {
        let p = Position {
            contig_id: MAX_CONTIG_ID,
            offset: MAX_OFFSET,
            strand: Strand::Reverse,
        };
        let e = PackedEntry::from_site_position(Site::default(), p);
        assert_eq!(e.position(), p);
    }

    #[test]
    #[should_panic(expected = "offset")]
    fn offset_overflow_panics() {
        let p = Position {
            contig_id: 0,
            offset: MAX_OFFSET + 1, // would overflow into the strand bit
            strand: Strand::Forward,
        };
        let _ = PackedEntry::from_site_position(Site::default(), p);
    }

    #[test]
    fn size_is_24_bytes() {
        assert_eq!(std::mem::size_of::<PackedEntry>(), 24);
        assert_eq!(std::mem::align_of::<PackedEntry>(), 8);
    }

    #[test]
    fn bytemuck_castable() {
        // Sanity: the bytemuck::Pod derive lets us safely cast a slice of
        // PackedEntry to bytes and back — this is the property the v2
        // on-disk format relies on for zero-copy mmap reads.
        let entries = [
            PackedEntry::from_site_position(
                Site { low: 1, high: 2 },
                Position {
                    contig_id: 3,
                    offset: 4,
                    strand: Strand::Forward,
                },
            ),
            PackedEntry::from_site_position(
                Site { low: 5, high: 6 },
                Position {
                    contig_id: 7,
                    offset: 8,
                    strand: Strand::Reverse,
                },
            ),
        ];
        let bytes: &[u8] = bytemuck::cast_slice(&entries);
        assert_eq!(bytes.len(), 2 * 24);
        let back: &[PackedEntry] = bytemuck::cast_slice(bytes);
        assert_eq!(back, &entries);
    }
}
