//! Genomic coordinate types.
//!
//! Convention: `offset` is always the 0-indexed position of the *leftmost*
//! base of the site on the *forward* strand, regardless of which strand the
//! protospacer lives on. This matches SAM/BED conventions and lets the
//! aggregator stay strand-agnostic.

use serde::{Deserialize, Serialize};

/// Which DNA strand a CRISPR site sits on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// SAM/BED-style single-character representation: `+` or `-`.
    #[inline]
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::Forward => '+',
            Self::Reverse => '-',
        }
    }
}

/// Location of a single CRISPR site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Position {
    pub contig_id: u32,
    pub offset: u32,
    pub strand: Strand,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strand_chars() {
        assert_eq!(Strand::Forward.as_char(), '+');
        assert_eq!(Strand::Reverse.as_char(), '-');
    }

    #[test]
    fn position_is_pod() {
        let p = Position {
            contig_id: 7,
            offset: 12345,
            strand: Strand::Reverse,
        };
        assert_eq!(p.contig_id, 7);
        assert_eq!(p.offset, 12345);
        assert_eq!(p.strand, Strand::Reverse);
    }
}
