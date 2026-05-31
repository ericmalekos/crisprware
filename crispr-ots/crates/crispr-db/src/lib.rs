//! On-disk database format and the build-time site discovery pipeline.
//!
//! - [`fasta`] — streaming reader (uncompressed + gzipped)
//! - [`position`] — `Position` + `Strand`
//! - [`packed`] — `PackedEntry` (24-byte on-disk + in-memory record)
//! - [`finder`] — `SiteFinder` (enzyme + sequence → sites on both strands)
//! - [`bin`] — `BinTable` (bucket sites by protospacer 5'-end prefix)
//! - [`mmap_db`] — zero-copy reader for the v2 on-disk format

pub mod bin;
pub mod fasta;
pub mod finder;
pub mod mmap_db;
pub mod packed;
pub mod position;
pub mod source;

pub use bin::BinTable;
pub use fasta::{read_fasta, read_fasta_from, Contig};
pub use finder::{SiteFinder, SiteSink};
pub use mmap_db::{write_binary_db, MmapDb, MmapDbError};
pub use packed::PackedEntry;
pub use position::{Position, Strand};
pub use source::BinSource;
