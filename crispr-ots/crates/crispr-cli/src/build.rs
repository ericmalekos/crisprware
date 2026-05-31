//! `crispr-ots build`: turn a reference FASTA into a persisted database.
//!
//! Pure function from `(reference path, enzyme, output path)` to a file on
//! disk in the crispr-ots format. The file is intended to be reused across
//! many subsequent `discover` invocations against the same reference and
//! enzyme — the build is the expensive step (~3 s for chr22 NGG, would
//! be tens of seconds for a full human genome) while a discover against a
//! pre-built database is fast.

use std::path::PathBuf;

use crispr_db::{read_fasta, BinTable, SiteFinder};
use crispr_encoding::Enzyme;

/// Default bin-prefix width. Empirically near-optimal for a chr22-sized
/// reference (~5 M sites). For the full human genome (~250 M sites),
/// width 10 is the sweet spot — see `BENCHMARKS.md` for the cost model
/// and the empirical sweep.
pub const DEFAULT_BIN_WIDTH: u8 = 7;

/// Arguments for one `build` invocation.
#[derive(Debug, Clone)]
pub struct BuildConfig {
    /// Reference FASTA path. Gzip auto-detected by `.gz` extension.
    pub reference: PathBuf,
    /// Enzyme to scan for. Determines PAM, protospacer length, and the
    /// `compare_mask` baked into the persisted database.
    pub enzyme: Enzyme,
    /// Output path. The conventional extension is `.crot` (CRISPR-OTS).
    pub output: PathBuf,
    /// Width of the bin-prefix key, in bases. `None` uses
    /// [`DEFAULT_BIN_WIDTH`]. The width is recorded in the persisted
    /// database so `discover` automatically uses the matching prefix.
    ///
    /// Choosing this is a space-vs-speed trade-off:
    /// - Bin metadata size scales as `4^width × 16 B`. At width 7 the
    ///   offset table is 256 KB; at width 12 it is 256 MB.
    /// - Scan time drops as the prefilter sharpens (fewer false
    ///   positives per bin) until the prefilter outer loop dominates.
    /// - The empirical sweet spot on chr22 (~5 M sites) is around 8-9;
    ///   on full human (~250 M sites), around 10.
    pub bin_width: Option<u8>,
}

#[derive(Debug)]
pub enum BuildError {
    Io(std::io::Error),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl std::error::Error for BuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
        }
    }
}

impl From<std::io::Error> for BuildError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

/// Build a database from a reference FASTA and persist it to disk.
///
/// # Errors
/// Surfaces FASTA-read errors, FASTA-encoding errors (non-ACGT bases on
/// otherwise-valid PAM-flanked windows are filtered, but unrecoverable I/O
/// problems propagate), and serialization failures.
pub fn run_build(config: &BuildConfig) -> Result<(), BuildError> {
    let table = build_table_in_memory(config)?;
    table.save_to_path(&config.output)?;
    Ok(())
}

/// Build a `BinTable` without persisting it. Used by integration tests and
/// any future caller that wants to keep the table in memory rather than
/// reload it from disk.
///
/// # Errors
/// Same as [`run_build`] minus the serialization step.
///
/// # Panics
/// Panics (via `BinTable::with_width`) if `config.bin_width` is set to a
/// value outside `1..=15` or greater than the enzyme's protospacer length.
/// The CLI layer pre-validates the range; library callers must do the same.
pub fn build_table_in_memory(config: &BuildConfig) -> Result<BinTable, BuildError> {
    let contigs = read_fasta(&config.reference)?;
    let finder = SiteFinder::new(config.enzyme.clone());
    let mut table = BinTable::with_width(
        config.enzyme.clone(),
        config.bin_width.unwrap_or(DEFAULT_BIN_WIDTH),
    );
    for contig in &contigs {
        let id = table.register_contig(contig.name.clone());
        finder.scan(&contig.sequence, id, &mut table);
    }
    Ok(table)
}
