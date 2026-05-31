//! 2-bit DNA encoding, IUPAC PAM expansion, and enzyme configuration.
//!
//! Foundation crate for the workspace. Hosts:
//! - [`Base`] (2-bit DNA enum) and [`IupacCode`] (PAM ambiguity codes)
//! - [`Site`] (2 × `u64`, up to 28 bp comparable) and [`SiteMask`]
//! - [`Enzyme`] configuration with IUPAC PAM expansion
//! - [`mismatches_masked`] — the core bit-parallel mismatch primitive
//!
//! See the workspace `NOTES.md` for the caveats inherited from FlashFry's
//! `BitEncoding` (24-base hard cap, last-char-at-LSB convention) and
//! GuideScan2's PAM handling (run-global 5'/3' flag).

pub mod base;
pub mod enzyme;
pub mod mismatch;
pub mod site;

pub use base::{Base, IupacCode};
pub use enzyme::{expand_iupac, Enzyme, PamSide};
pub use mismatch::mismatches_masked;
pub use site::{Site, SiteMask, MAX_SITE_LEN};
