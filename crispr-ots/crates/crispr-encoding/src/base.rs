//! 2-bit DNA base encoding and IUPAC ambiguity codes.
//!
//! The 2-bit assignments (A=0, C=1, G=2, T=3) match FlashFry's convention
//! so lookup tables (CFD penalty maps, Doench coefficients) transcribe
//! directly to flat arrays indexed by the bit pattern.

use std::fmt;

use serde::{Deserialize, Serialize};

/// A single DNA base, encoded as 2 bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    /// All four bases in canonical order.
    pub const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Raw 2-bit value used in `Site` encoding.
    #[inline]
    #[must_use]
    pub const fn bits(self) -> u8 {
        self as u8
    }

    /// Decode from a 2-bit value. Returns `None` for inputs `> 3`.
    #[inline]
    #[must_use]
    pub const fn from_bits(bits: u8) -> Option<Self> {
        match bits {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            _ => None,
        }
    }

    /// Decode from ASCII (case-insensitive). Returns `None` for anything
    /// other than `A`, `C`, `G`, `T`.
    #[inline]
    #[must_use]
    pub const fn from_ascii(c: u8) -> Option<Self> {
        match c {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            _ => None,
        }
    }

    /// Uppercase ASCII byte representation.
    #[inline]
    #[must_use]
    pub const fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }

    /// Watson-Crick complement.
    #[inline]
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_ascii() as char)
    }
}

/// IUPAC nucleotide ambiguity codes. Used to describe PAM patterns at the
/// configuration layer; expanded to literal `Base` sequences before any
/// runtime matching so the hot path only ever sees A/C/G/T.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[allow(clippy::upper_case_acronyms)]
pub enum IupacCode {
    A,
    C,
    G,
    T,
    /// A or G
    R,
    /// C or T
    Y,
    /// G or C
    S,
    /// A or T
    W,
    /// G or T
    K,
    /// A or C
    M,
    /// C or G or T (not A)
    B,
    /// A or G or T (not C)
    D,
    /// A or C or T (not G)
    H,
    /// A or C or G (not T)
    V,
    /// any of A, C, G, T
    N,
}

impl IupacCode {
    /// Bases this code matches.
    #[must_use]
    pub const fn bases(self) -> &'static [Base] {
        match self {
            Self::A => &[Base::A],
            Self::C => &[Base::C],
            Self::G => &[Base::G],
            Self::T => &[Base::T],
            Self::R => &[Base::A, Base::G],
            Self::Y => &[Base::C, Base::T],
            Self::S => &[Base::C, Base::G],
            Self::W => &[Base::A, Base::T],
            Self::K => &[Base::G, Base::T],
            Self::M => &[Base::A, Base::C],
            Self::B => &[Base::C, Base::G, Base::T],
            Self::D => &[Base::A, Base::G, Base::T],
            Self::H => &[Base::A, Base::C, Base::T],
            Self::V => &[Base::A, Base::C, Base::G],
            Self::N => &Base::ALL,
        }
    }

    /// Parse from an ASCII character (case-insensitive). Returns `None` for
    /// any character outside the IUPAC nucleotide set.
    #[inline]
    #[must_use]
    pub const fn from_ascii(c: u8) -> Option<Self> {
        match c {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            b'R' | b'r' => Some(Self::R),
            b'Y' | b'y' => Some(Self::Y),
            b'S' | b's' => Some(Self::S),
            b'W' | b'w' => Some(Self::W),
            b'K' | b'k' => Some(Self::K),
            b'M' | b'm' => Some(Self::M),
            b'B' | b'b' => Some(Self::B),
            b'D' | b'd' => Some(Self::D),
            b'H' | b'h' => Some(Self::H),
            b'V' | b'v' => Some(Self::V),
            b'N' | b'n' => Some(Self::N),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bits_roundtrip() {
        for b in Base::ALL {
            assert_eq!(Base::from_bits(b.bits()), Some(b));
        }
    }

    #[test]
    fn ascii_roundtrip() {
        for b in Base::ALL {
            assert_eq!(Base::from_ascii(b.to_ascii()), Some(b));
        }
        // Case-insensitive parse.
        assert_eq!(Base::from_ascii(b'a'), Some(Base::A));
        assert_eq!(Base::from_ascii(b't'), Some(Base::T));
        // Out of charset.
        assert_eq!(Base::from_ascii(b'N'), None);
        assert_eq!(Base::from_ascii(b'X'), None);
    }

    #[test]
    fn complement_is_involution() {
        for b in Base::ALL {
            assert_eq!(b.complement().complement(), b);
        }
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::C.complement(), Base::G);
    }

    #[test]
    fn iupac_n_expands_to_four() {
        assert_eq!(IupacCode::N.bases().len(), 4);
        for b in Base::ALL {
            assert!(IupacCode::N.bases().contains(&b));
        }
    }

    #[test]
    fn iupac_v_excludes_t() {
        let v = IupacCode::V.bases();
        assert_eq!(v.len(), 3);
        assert!(!v.contains(&Base::T));
    }

    #[test]
    fn iupac_single_bases_are_singletons() {
        for (code, base) in [
            (IupacCode::A, Base::A),
            (IupacCode::C, Base::C),
            (IupacCode::G, Base::G),
            (IupacCode::T, Base::T),
        ] {
            assert_eq!(code.bases(), &[base]);
        }
    }

    #[test]
    fn iupac_two_base_codes() {
        assert_eq!(IupacCode::R.bases(), &[Base::A, Base::G]);
        assert_eq!(IupacCode::Y.bases(), &[Base::C, Base::T]);
        assert_eq!(IupacCode::W.bases(), &[Base::A, Base::T]);
    }
}
