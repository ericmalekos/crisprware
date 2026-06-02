//! Library half of the `crispr-ots` CLI.
//!
//! The thin binary in `main.rs` parses arguments via `clap` and dispatches
//! to functions exposed here. Splitting it this way lets integration tests
//! exercise the full pipeline in-process without spawning child processes.

pub mod build;
pub mod csv_out;
pub mod discover;
pub mod kmers_csv;
pub mod ot_stream;
pub mod tsv;

pub use build::{build_table_in_memory, run_build, BuildConfig, BuildError};
pub use crispr_score::SpecConvention;
pub use discover::{
    run_discover, run_discover_with, DiscoverConfig, DiscoverError, DiscoverInput, OutputFormat,
    OtFormat, OutputMode, ScannerKind, ScoreMetric,
};
