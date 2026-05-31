//! Off-target scan engine.
//!
//! Phase 1c provides the scalar reference implementation:
//! - [`Guide`] and [`Hit`] result types
//! - the [`Scanner`] trait
//! - [`BinScanner`], the FlashFry-style bin-scan implementation against an
//!   in-memory [`crispr_db::BinTable`]
//!
//! Later phases will add SIMD batching, rayon-parallel bin iteration, and
//! (optionally, behind the same trait) an FM-index `Scanner` for users who
//! want GuideScan2's flexibility profile.

pub mod bin_scanner;

pub use bin_scanner::{BinScanner, Guide, Hit, Scanner};
