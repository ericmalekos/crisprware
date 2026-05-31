//! Scoring metrics for CRISPR off-targets.
//!
//! Day 1: Doench 2016 [`Cfd`], computed directly on bit-encoded `Site`s
//! (no string decoding — see the `NOTES.md` caveat about FlashFry and
//! GuideScan2 both paying that cost in their CFD inner loops).
//!
//! Future scorers (MIT/Hsu 2013, ClosestHit, JostAndSantos, …) will live
//! behind a `ScoreModel` trait once we have more than one to abstract
//! over.

pub mod cfd;

pub use cfd::{Cfd, CfdResult, SpecConvention, CFD_MINIMUM_THRESHOLD};
