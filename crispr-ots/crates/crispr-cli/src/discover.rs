//! `crispr-ots enumerate` (a.k.a. `discover`): build per-query off-target
//! reports against an already-built `BinTable` / `MmapDb`.
//!
//! Two input shapes are supported:
//! - [`DiscoverInput::QueryFasta`] — the original native form. Sites are
//!   discovered inside the FASTA via `SiteFinder` and one report row is
//!   emitted per discovered site.
//! - [`DiscoverInput::KmersCsv`] — the GuideScan2 / crisprware form. Each
//!   row of the CSV is one explicit guide with its own id, chromosome,
//!   and position.
//!
//! Three output shapes are supported:
//! - [`OutputFormat::Csv`] — GuideScan2-style: one row per (guide, hit)
//!   with columns `id,sequence,match_chrm,match_position,match_strand,
//!   match_distance,specificity`. The format crisprware's
//!   `score_guides.py` reads back.
//! - [`OutputFormat::Tsv`] — FlashFry-style: one row per guide with
//!   `otCount`, `otSequences` (full off-target listing within
//!   `--mismatches`), and optional CFD columns.
//! - [`OutputFormat::Both`] — write the CSV at `output` and the TSV at
//!   `<output>.detail.tsv`.
//!
//! The specificity column's denominator convention is selectable via
//! [`SpecConvention`] (see the crispr-score docs).

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crispr_db::{read_fasta, BinSource, Position, SiteFinder, Strand};
use crispr_encoding::Site;
use crispr_scan::{BinScanner, Guide, Hit, Scanner};
use crispr_score::{Cfd, CfdResult, SpecConvention};

use crate::csv_out::{self, CsvRow};
use crate::kmers_csv::{self, KmerEntry};
use crate::tsv::{write_rows as write_tsv_rows, Cas12aColumns, CfdColumns, DiscoverRow};
use crispr_score::{Cas12aCfd, Cas12aMatrix, Cas12aResult};

/// Scoring metrics that can be applied during discover.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreMetric {
    /// Doench 2016 CFD. SpCas9-only; the model has no biological
    /// interpretation for other enzymes.
    Cfd,
    /// Cas12a "CFD-like" multiplicative score backed by the
    /// `2xNLS-Cas12a` activity matrix (the `AsCas12a` wild-type-ish
    /// profile from `parasol_scripts/off_targ_2xNLS_Cas12a.csv`).
    Cas12aTwoXNls,
    /// Cas12a CFD-like with the `enCas12a` engineered broad-PAM matrix.
    Cas12aEnCas12a,
}

/// Output-file shape for `enumerate`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// One row per (guide, off-target hit). GuideScan2-compatible.
    Csv,
    /// One row per guide. FlashFry-style; lists every off-target up to
    /// `--mismatches` in a single comma-separated column.
    Tsv,
    /// Both. CSV at `--output`; TSV at `<output>.detail.tsv`.
    Both,
}

/// Scalable streaming scored-output modes (`--output-mode`). Unlike the legacy
/// [`OutputFormat`] CSV/TSV (which materialize every hit in memory before
/// writing), these stream guide batches to disk with bounded memory.
/// [`OutputMode::None`] (default) keeps the legacy `--format` path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputMode {
    /// Use the legacy `--format` path (no streaming).
    #[default]
    None,
    /// Mode 1: one row per guide with the aggregated specificity, computed
    /// over *all* off-targets (no CFD floor).
    Aggregated,
    /// Mode 2: one record per (guide, off-target) — chrom/start/strand/CFD —
    /// floored at `--cfd-threshold` for disk control.
    PerOffTarget,
    /// Both Mode 1 and Mode 2, from a single streaming pass.
    Both,
}

/// On-disk encoding for the Mode-2 per-off-target file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OtFormat {
    /// Compact 12-byte binary records (default) — smallest on disk; read with
    /// the bundled aggregator.
    #[default]
    Binary,
    /// Plain TSV rows (`guide_id\tchrom\tstart\tstrand\tcfd`) — greppable,
    /// ~2-3× larger.
    Tsv,
    /// Both binary and TSV.
    Both,
}

/// Which scan backend to run the off-target search on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScannerKind {
    /// Multithreaded SIMD CPU bin-scanner (AVX-512 VPOPCNTDQ where available,
    /// otherwise AVX2/scalar via `wide`). The default.
    #[default]
    Cpu,
    /// CUDA GPU kernel (`crispr-scan-gpu`). Requires a binary built with
    /// `--features gpu` and an available NVIDIA device.
    Gpu,
}

/// Input-file shape for `enumerate`.
#[derive(Debug, Clone)]
pub enum DiscoverInput {
    /// FASTA from which the scanner discovers candidate sites.
    QueryFasta(PathBuf),
    /// Per-guide CSV in crisprware/GuideScan2 format. See [`kmers_csv`].
    KmersCsv(PathBuf),
}

/// Configuration for one `enumerate` invocation. The bin source is passed
/// separately to [`run_discover`] so callers choose whether to load from
/// disk or build in memory.
#[derive(Debug, Clone)]
pub struct DiscoverConfig {
    pub input: DiscoverInput,
    pub max_mismatches: u8,
    pub scores: Vec<ScoreMetric>,
    pub output: PathBuf,
    pub format: OutputFormat,
    pub spec_convention: SpecConvention,
    /// If `Some(t)`, drop any guide that has an off-target within `t`
    /// mismatches (a 0-mm duplicate, or any 1..=t-mm match). Matches
    /// `guidescan enumerate --threshold`. In the streaming output path this
    /// drives a cheap mm≤t pre-screen: such guides are dropped *before* the
    /// full mm≤`max_mismatches` CFD scan, which then runs only on the
    /// survivors. `None` disables filtering.
    pub threshold: Option<u8>,
    /// If `Some(n)`, keep at most `n` distinct off-target sequences per
    /// `(guide, mismatch-count)` bin. The cap is applied while the raw
    /// hit stream is being grouped — once a bin reaches `n` distinct
    /// sequences, additional sequences at that mismatch level for the
    /// same guide are dropped. `None` disables the cap.
    ///
    /// Affects `otCount`, the `otSequences` listing in TSV output, and
    /// the CFD denominator (since aggregation sums over kept sequences
    /// only). Useful for taming low-complexity guides that explode to
    /// 100k+ off-targets at high mismatch counts.
    pub max_per_bin: Option<u32>,
    /// Streaming scored-output mode (`--output-mode`). [`OutputMode::None`]
    /// (default) routes to the legacy `--format` path; any other value runs
    /// the batched [`run_discover_streaming`] driver (GPU backend).
    pub output_mode: OutputMode,
    /// Encoding for the Mode-2 per-off-target file (`--ot-format`).
    pub ot_format: OtFormat,
    /// CFD floor applied to the Mode-2 file *only* (Mode 1 is never floored):
    /// off-targets with CFD below this are dropped from disk. `0.0` keeps all.
    pub cfd_threshold: f64,
    /// Guides per streaming batch ceiling (`None`/`Some(0)` = a built-in
    /// default).
    pub batch_size: Option<usize>,
    /// Max predicted hits per streaming batch before flushing — the real
    /// memory bound (host buffer ≈ `hit_budget × 12 B`).
    pub hit_budget: u64,
    /// Per-guide off-target cap for streamed scored output (`--max-off-targets`),
    /// applied on every backend (CPU, GPU-emit, GPU-CFD-kernel). Once a guide
    /// accumulates this many on+off-target hits, scoring stops and the guide is
    /// flagged `saturated` (its specificity is ≈ 0 regardless). Bounds per-guide
    /// work so a few homopolymer guides can't gate genome-scale wall time.
    /// `0` = unlimited (exact, but a handful of guides can be very slow).
    pub max_off_targets: u32,
    /// Specificity floor for streamed scored output (`--min-specificity`),
    /// applied on every backend. When `> 0`, a guide is flagged `saturated` and
    /// scanning stops the moment its off-target CFD sum proves specificity has
    /// dropped below this value (`off_sum > 1/min_specificity - 1`). The
    /// principled cap: cuts on the score itself, and exits the worst guides
    /// almost immediately. `0` = disabled. Acts together with `max_off_targets`
    /// (stop on either).
    pub min_specificity: f64,
    /// When [`threshold`](Self::threshold) drives a pre-screen, keep the dropped
    /// (non-specific) guides in the Mode-1 output with a `dropped` marker rather
    /// than omitting them. `false` (default) = omit; output is survivors only.
    pub keep_dropped: bool,
}

impl Default for DiscoverConfig {
    fn default() -> Self {
        Self {
            input: DiscoverInput::KmersCsv(PathBuf::new()),
            max_mismatches: 4,
            scores: Vec::new(),
            output: PathBuf::new(),
            format: OutputFormat::Csv,
            spec_convention: SpecConvention::Guidescan,
            threshold: None,
            max_per_bin: Some(500),
            output_mode: OutputMode::None,
            ot_format: OtFormat::Binary,
            cfd_threshold: 0.023,
            batch_size: None,
            hit_budget: 64_000_000,
            max_off_targets: 0,
            min_specificity: 0.0,
            keep_dropped: false,
        }
    }
}

/// Errors surfaced by the discover pipeline. Downstream layers wrap with
/// `anyhow` for friendly CLI messages.
#[derive(Debug)]
pub enum DiscoverError {
    Io(std::io::Error),
    QueryEmpty,
    KmersCsv(kmers_csv::KmerCsvError),
    /// GPU backend requested but unavailable (binary built without
    /// `--features gpu`), or a CUDA error occurred during setup or scanning.
    Gpu(String),
}

impl std::fmt::Display for DiscoverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::QueryEmpty => write!(f, "query FASTA contained no contigs"),
            Self::KmersCsv(e) => write!(f, "{e}"),
            Self::Gpu(m) => write!(f, "GPU scan error: {m}"),
        }
    }
}

impl std::error::Error for DiscoverError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::QueryEmpty | Self::Gpu(_) => None,
            Self::KmersCsv(e) => Some(e),
        }
    }
}

impl From<std::io::Error> for DiscoverError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<kmers_csv::KmerCsvError> for DiscoverError {
    fn from(e: kmers_csv::KmerCsvError) -> Self {
        Self::KmersCsv(e)
    }
}

/// Run the discover pipeline against a pre-built bin source.
///
/// # Errors
/// Returns [`DiscoverError`] for I/O failures, empty query inputs, or
/// malformed kmers-CSV rows.
///
/// # Panics
/// Panics if any `usize` index doesn't fit in `u32` — only possible with
/// pathologically large inputs (billions of contigs or guides). Treated
/// as unreachable for realistic genome-scale workloads.
pub fn run_discover(source: &dyn BinSource, config: &DiscoverConfig) -> Result<(), DiscoverError> {
    run_discover_with(source, config, ScannerKind::Cpu)
}

/// As [`run_discover`], but with an explicit scan backend ([`ScannerKind`]).
///
/// # Errors
/// As [`run_discover`], plus [`DiscoverError::Gpu`] when the GPU backend is
/// requested on a binary built without `--features gpu`, or a CUDA error
/// occurs during GPU setup or scanning.
///
/// # Panics
/// As [`run_discover`].
#[allow(clippy::too_many_lines)]
pub fn run_discover_with(
    source: &dyn BinSource,
    config: &DiscoverConfig,
    scanner: ScannerKind,
) -> Result<(), DiscoverError> {
    let trace = std::env::var("CRISPR_OTS_TRACE").is_ok();
    let t_start = std::time::Instant::now();

    let enzyme = source.enzyme().clone();
    let total_len = usize::from(enzyme.total_scan_len());

    // ---- Build the guide vector and side-table of metadata. ----
    let guide_records: Vec<GuideRecord> = match &config.input {
        DiscoverInput::QueryFasta(path) => collect_from_fasta(path, source)?,
        DiscoverInput::KmersCsv(path) => collect_from_kmers_csv(path, source)?,
    };

    if guide_records.is_empty() {
        return Err(DiscoverError::QueryEmpty);
    }

    // ---- Scan the reference for off-targets. ----
    let guides: Vec<Guide> = guide_records
        .iter()
        .enumerate()
        .map(|(i, g)| Guide {
            id: format!("g{i}"),
            site: g.site,
        })
        .collect();
    // ---- Streaming scored-output modes bypass the in-memory hit path. ----
    if config.output_mode != OutputMode::None {
        return run_discover_streaming(source, config, &guide_records, &guides, scanner);
    }

    let t_scan = std::time::Instant::now();
    let hits = match scanner {
        ScannerKind::Cpu => BinScanner::new(source).scan(&guides, config.max_mismatches),
        ScannerKind::Gpu => scan_on_gpu(source, &guides, config.max_mismatches)?,
    };
    if trace {
        eprintln!(
            "[trace] scan {} guides: {:.2}s ({} raw hits)",
            guides.len(),
            t_scan.elapsed().as_secs_f64(),
            hits.len()
        );
    }

    // ---- Group hits by (guide, off-target sequence) for aggregation. ----
    // We collapse over genomic positions to match FlashFry's reporting,
    // where multi-mapping OT sequences appear once with a count.
    //
    // `max_per_bin` (the `--max-off-targets-per-bin` flag) short-circuits
    // here: once a `(guide, mismatch-count)` bin has accumulated `n`
    // distinct off-target sequences, any *new* sequences at that mismatch
    // level for that guide are dropped. Sequences already in the bin
    // still have their position count incremented for genuine
    // multi-mappers — only *new* sequences are gated.
    let mut grouped: HashMap<usize, HashMap<Site, GroupAgg>> = HashMap::new();
    let mut bin_seq_count: HashMap<(usize, u8), u32> = HashMap::new();
    let mut dropped_seqs: u64 = 0;
    for hit in &hits {
        let gi = hit.guide_index as usize;
        let inner = grouped.entry(gi).or_default();
        if let Some(existing) = inner.get_mut(&hit.off_target) {
            existing.count += 1;
            debug_assert_eq!(existing.mismatches, hit.mismatches);
            continue;
        }
        if let Some(cap) = config.max_per_bin {
            let n = bin_seq_count.entry((gi, hit.mismatches)).or_insert(0);
            if *n >= cap {
                dropped_seqs += 1;
                continue;
            }
            *n += 1;
        }
        inner.insert(
            hit.off_target,
            GroupAgg {
                mismatches: hit.mismatches,
                count: 1,
            },
        );
    }
    if trace && dropped_seqs > 0 {
        eprintln!(
            "[trace] dropped {dropped_seqs} off-target sequences past --max-off-targets-per-bin cap"
        );
    }

    // ---- Index raw hits by guide for per-position CSV output. ----
    // Sorted by (contig_id, offset, strand) for stable output. We only
    // build this for CSV/Both modes since TSV doesn't need positions.
    let need_per_hit = matches!(config.format, OutputFormat::Csv | OutputFormat::Both);
    let hits_by_guide = if need_per_hit {
        Some(index_hits_by_guide(&hits))
    } else {
        None
    };

    let cfd = if config.scores.contains(&ScoreMetric::Cfd) {
        Some(Cfd::new())
    } else {
        None
    };
    let cas12a_scorer = if config.scores.contains(&ScoreMetric::Cas12aTwoXNls) {
        Some(Cas12aCfd::from_matrix(Cas12aMatrix::TwoXNls).expect("bundled 2xNLS matrix parses"))
    } else if config.scores.contains(&ScoreMetric::Cas12aEnCas12a) {
        Some(
            Cas12aCfd::from_matrix(Cas12aMatrix::EnCas12a)
                .expect("bundled enCas12a matrix parses"),
        )
    } else {
        None
    };
    let cas12a_protospacer_len: u32 = u32::from(enzyme.protospacer_len);

    // ---- Compute the per-guide CFD result + apply --threshold filter. ----
    let mut keep_guide: Vec<bool> = vec![true; guide_records.len()];
    let mut cfd_results: Vec<Option<CfdResult>> = vec![None; guide_records.len()];
    let mut cas12a_results: Vec<Option<Cas12aResult>> = vec![None; guide_records.len()];
    for (idx, rec) in guide_records.iter().enumerate() {
        let group = grouped.get(&idx);
        let all_offs: Vec<(Site, u8, u32)> = group
            .map(|m| {
                m.iter()
                    .map(|(s, agg)| (*s, agg.mismatches, agg.count))
                    .collect()
            })
            .unwrap_or_default();

        if let Some(t) = config.threshold {
            // Drop the guide if either:
            //   (a) some off-target sequence (mm > 0) has mm <= t, OR
            //   (b) the guide's own sequence appears at >= 2 genomic
            //       positions at mm = 0 (multi-mapping on-target — at
            //       least one of those copies is at an unintended
            //       location, so for the threshold's purposes it acts
            //       as a 0-mm off-target).
            // Matches guidescan2's --threshold semantics: per its docs,
            // a guide is dropped if any non-on-target hit is within
            // `threshold` mismatches, including 0-mm multi-mappers.
            let close_off = all_offs.iter().any(|(_, mm, count)| {
                (*mm > 0 && *mm <= t) || (*mm == 0 && *count >= 2)
            });
            if close_off {
                keep_guide[idx] = false;
                continue;
            }
        }

        if let Some(c) = cas12a_scorer.as_ref() {
            cas12a_results[idx] = Some(c.aggregate_with_len(
                rec.site,
                all_offs.iter().map(|(s, mm, count)| (*s, *mm, *count)),
                cas12a_protospacer_len,
            ));
        }
        if let Some(c) = cfd.as_ref() {
            cfd_results[idx] = Some(c.aggregate_with(
                rec.site,
                all_offs.iter().map(|(s, mm, count)| (*s, *mm, *count)),
                config.spec_convention,
            ));
        }
    }

    // ---- Emit outputs. ----
    match config.format {
        OutputFormat::Csv => {
            write_csv_output(
                &config.output,
                source,
                &guide_records,
                hits_by_guide.as_ref().unwrap(),
                &cfd_results,
                &cas12a_results,
                &keep_guide,
            )?;
        }
        OutputFormat::Tsv => {
            write_tsv_output(
                &config.output,
                &guide_records,
                &grouped,
                &cfd_results,
                &cas12a_results,
                &keep_guide,
                total_len,
            )?;
        }
        OutputFormat::Both => {
            write_csv_output(
                &config.output,
                source,
                &guide_records,
                hits_by_guide.as_ref().unwrap(),
                &cfd_results,
                &cas12a_results,
                &keep_guide,
            )?;
            let tsv_path = sidecar_tsv_path(&config.output);
            write_tsv_output(
                &tsv_path,
                &guide_records,
                &grouped,
                &cfd_results,
                &cas12a_results,
                &keep_guide,
                total_len,
            )?;
        }
    }

    if trace {
        eprintln!("[trace] total: {:.2}s", t_start.elapsed().as_secs_f64());
    }
    Ok(())
}

/// Run the scan on the CUDA backend. Present only when built with
/// `--features gpu`.
#[cfg(feature = "gpu")]
fn scan_on_gpu(
    source: &dyn BinSource,
    guides: &[Guide],
    max_mismatches: u8,
) -> Result<Vec<Hit>, DiscoverError> {
    let scanner = crispr_scan_gpu::GpuScanner::new(source)
        .map_err(|e| DiscoverError::Gpu(e.to_string()))?;
    Ok(scanner.scan(guides, max_mismatches))
}

/// Stub for builds without GPU support: always errors with a build hint.
#[cfg(not(feature = "gpu"))]
fn scan_on_gpu(
    _source: &dyn BinSource,
    _guides: &[Guide],
    _max_mismatches: u8,
) -> Result<Vec<Hit>, DiscoverError> {
    Err(DiscoverError::Gpu(
        "this crispr-ots binary was built without GPU support; rebuild with \
         `cargo build --features gpu` to use `--scanner gpu`"
            .to_string(),
    ))
}

// ======================================================================
// Streaming scored-output driver (`--output-mode`)
// ======================================================================

/// Per-guide bucket array length: counts off-targets with exactly `i`
/// mismatches for `i` in `0..MM_BUCKETS` (index 0 is unused — the mm-0 bucket
/// is derived from `on_count`). 8 covers every realistic `--mismatches`.
const MM_BUCKETS: usize = 8;

/// Primary + optional secondary Cas12a CFD scorer for one dual-matrix scan.
/// The PRIMARY drives `off_sum`/`tttv_sum`, `max_cfd`, the saturation cap, the
/// Mode-2 `cfd` column, and the `--cfd-threshold` floor; the SECONDARY (when
/// present) drives only the `*_2` specificity pair and the Mode-2 `cfd2` column.
struct Cas12aScorers {
    primary: Cas12aCfd,
    secondary: Option<Cas12aCfd>,
}

/// Per-guide CFD accumulators for the streaming driver (Mode 1).
struct Accum {
    /// Σ cfd over off-targets (Cas9), or Σ TTTN (Cas12a primary matrix).
    off_sum: Vec<f64>,
    /// Σ over off-targets whose PAM ≠ TTTT (Cas12a primary TTTV only).
    tttv_sum: Vec<f64>,
    /// Σ TTTN over off-targets for the secondary Cas12a matrix (dual runs only).
    off_sum2: Vec<f64>,
    /// Σ TTTV over off-targets for the secondary Cas12a matrix (dual runs only).
    tttv_sum2: Vec<f64>,
    /// Largest off-target CFD seen per guide (primary matrix).
    max_cfd: Vec<f64>,
    /// On-target (mm == 0) multiplicity per guide.
    on_count: Vec<u32>,
    /// Off-target (mm > 0) position count per guide.
    off_count: Vec<u32>,
    /// Off-targets bucketed by exact mismatch count; index `m` = number of
    /// off-targets at `m` mismatches (`1..=max_mm` used; mm-0 derived from
    /// `on_count`). Feeds the track's `_mismatchCounts`.
    mm_counts: Vec<[u32; MM_BUCKETS]>,
    /// Set when a per-guide cap was hit (any backend) — the guide's sums are
    /// partial (specificity is an upper bound, ≈ 0).
    saturated: Vec<bool>,
}

/// Mode-2 streaming writers (binary and/or TSV).
struct OtWriters {
    bin: Option<crate::ot_stream::OtBinaryWriter>,
    tsv: Option<crate::ot_stream::OtTsvWriter>,
}

/// Apply one off-target hit: accumulate into Mode-1 sums (no floor) and, if
/// per-off-target output is on and the CFD clears the floor, write its
/// Mode-2 record(s).
#[allow(clippy::too_many_arguments)]
fn accumulate_hit(
    abs_guide: usize,
    guide_site: Site,
    ot_site: Site,
    ot_pos: Position,
    mm: u8,
    cfd: Option<&Cfd>,
    cas12a: Option<&Cas12aScorers>,
    proto_len: u32,
    cfd_threshold: f64,
    ot_cap: u32,
    off_sum_cap: f64,
    acc: &mut Accum,
    writers: &mut OtWriters,
    source: &dyn BinSource,
) -> Result<(), DiscoverError> {
    // Once capped, stop accumulating + writing Mode-2 records for this guide —
    // mirrors the GPU CFD kernel's early-exit so every backend agrees.
    if acc.saturated[abs_guide] {
        return Ok(());
    }
    if mm == 0 {
        acc.on_count[abs_guide] += 1;
        flag_if_saturated(acc, abs_guide, ot_cap, off_sum_cap);
        return Ok(());
    }
    // `is_tttt` is matrix-independent (the PAM isn't scored) — compute it once.
    // The primary CFD drives all aggregates + the Mode-2 `cfd`; the secondary
    // (dual Cas12a only) drives only `off_sum2`/`tttv_sum2` + the `cfd2` column.
    let (cfd_val, cfd2_opt, is_tttt) = if let Some(c) = cas12a {
        let is_tttt = crispr_score::is_tttt_prefix(ot_site, proto_len);
        let primary = c.primary.score_pair_with_len(guide_site, ot_site, proto_len);
        let secondary = c
            .secondary
            .as_ref()
            .map(|s| s.score_pair_with_len(guide_site, ot_site, proto_len));
        (primary, secondary, is_tttt)
    } else if let Some(c) = cfd {
        (c.score_pair(guide_site, ot_site), None, false)
    } else {
        (0.0, None, false)
    };
    acc.off_sum[abs_guide] += cfd_val;
    acc.off_count[abs_guide] += 1;
    if (mm as usize) < MM_BUCKETS {
        acc.mm_counts[abs_guide][mm as usize] += 1;
    }
    if cas12a.is_some() && !is_tttt {
        acc.tttv_sum[abs_guide] += cfd_val;
    }
    if let Some(cfd2_val) = cfd2_opt {
        acc.off_sum2[abs_guide] += cfd2_val;
        if !is_tttt {
            acc.tttv_sum2[abs_guide] += cfd2_val;
        }
    }
    if cfd_val > acc.max_cfd[abs_guide] {
        acc.max_cfd[abs_guide] = cfd_val;
    }
    // Floor the off-target LIST on EITHER matrix: list the off-target if the
    // primary OR the secondary CFD clears the threshold (counts/specificity are
    // unfloored regardless). Single-matrix runs reduce to the primary check.
    let pass_floor = cfd_val >= cfd_threshold || cfd2_opt.is_some_and(|c2| c2 >= cfd_threshold);
    if (writers.bin.is_some() || writers.tsv.is_some()) && pass_floor {
        let gid = u32::try_from(abs_guide).unwrap_or(u32::MAX);
        let strand_rev = matches!(ot_pos.strand, Strand::Reverse);
        if let Some(w) = writers.bin.as_mut() {
            w.write_record(crate::ot_stream::OtRecord {
                guide_id: gid,
                contig: u16::try_from(ot_pos.contig_id).unwrap_or(u16::MAX),
                offset: ot_pos.offset,
                strand_reverse: strand_rev,
                cfd_q: crate::ot_stream::OtRecord::quantize_cfd(cfd_val),
                cfd2_q: cfd2_opt.map_or(0, crate::ot_stream::OtRecord::quantize_cfd),
                mm,
            })?;
        }
        if let Some(w) = writers.tsv.as_mut() {
            let chrom = source.contig_name(ot_pos.contig_id).unwrap_or("");
            w.write_row(gid, chrom, ot_pos.offset, strand_rev, mm, cfd_val, cfd2_opt)?;
        }
    }
    flag_if_saturated(acc, abs_guide, ot_cap, off_sum_cap);
    Ok(())
}

/// Flag a guide `saturated` once either per-guide cap is hit: count
/// (`on + off >= ot_cap`) or off-target CFD sum (`off_sum > off_sum_cap`,
/// where `off_sum_cap = 1/min_specificity - 1`). `0` / `0.0` disable the
/// respective cap. Shared by the CPU and GPU-emit accumulation paths so they
/// agree with the GPU CFD kernel's in-kernel early-exit.
fn flag_if_saturated(acc: &mut Accum, g: usize, ot_cap: u32, off_sum_cap: f64) {
    if (ot_cap != 0 && acc.on_count[g] + acc.off_count[g] >= ot_cap)
        || (off_sum_cap > 0.0 && acc.off_sum[g] > off_sum_cap)
    {
        acc.saturated[g] = true;
    }
}

/// SpCas9 specificity from the off-target CFD sum + on-target multiplicity.
fn specificity(convention: SpecConvention, off_sum: f64, on_count: u32) -> f64 {
    match convention {
        SpecConvention::Flashfry => {
            if off_sum > 0.0 {
                1.0 / (1.0 + off_sum)
            } else {
                1.0
            }
        }
        SpecConvention::Guidescan => {
            let on = if on_count == 0 { 1.0 } else { f64::from(on_count) };
            let denom = on + off_sum;
            if denom > 0.0 {
                1.0 / denom
            } else {
                1.0
            }
        }
    }
}

/// `<output>.<suffix>` (e.g. `scored.csv` → `scored.csv.ot.bin`).
fn suffix_path(output: &Path, suffix: &str) -> PathBuf {
    let mut s = output.as_os_str().to_owned();
    s.push(".");
    s.push(suffix);
    PathBuf::from(s)
}

/// Streaming/batched scored-output driver. Scans `guides` in batches against
/// the resident (GPU) or mmap'd (CPU) index and streams the aggregated
/// (Mode 1) and/or per-off-target (Mode 2) outputs to disk; peak memory is
/// `O(#guides + batch hits)`, independent of total hit volume.
fn run_discover_streaming(
    source: &dyn BinSource,
    config: &DiscoverConfig,
    guide_records: &[GuideRecord],
    guides: &[Guide],
    scanner: ScannerKind,
) -> Result<(), DiscoverError> {
    let enzyme = source.enzyme();
    let proto_len = u32::from(enzyme.protospacer_len);

    let cfd = if config.scores.contains(&ScoreMetric::Cfd) {
        Some(Cfd::new())
    } else {
        None
    };
    // Collect the requested Cas12a matrices in listed order (dedup). The first is
    // the primary; the second distinct one (if any) is the secondary "WT" matrix,
    // scored in the same scan. Any 3rd distinct matrix is ignored.
    let cas12a = {
        let mut mats: Vec<Cas12aMatrix> = Vec::new();
        for s in &config.scores {
            let m = match s {
                ScoreMetric::Cas12aEnCas12a => Some(Cas12aMatrix::EnCas12a),
                ScoreMetric::Cas12aTwoXNls => Some(Cas12aMatrix::TwoXNls),
                _ => None,
            };
            if let Some(m) = m {
                if !mats.contains(&m) {
                    mats.push(m);
                }
            }
        }
        mats.first().map(|&m0| Cas12aScorers {
            primary: Cas12aCfd::from_matrix(m0).expect("bundled Cas12a matrix parses"),
            secondary: mats
                .get(1)
                .map(|&m1| Cas12aCfd::from_matrix(m1).expect("bundled Cas12a matrix parses")),
        })
    };
    let is_cas12a = cas12a.is_some();
    let cas12a_dual = cas12a.as_ref().is_some_and(|c| c.secondary.is_some());

    let n = guide_records.len();
    let mut acc = Accum {
        off_sum: vec![0.0; n],
        tttv_sum: vec![0.0; n],
        off_sum2: vec![0.0; n],
        tttv_sum2: vec![0.0; n],
        max_cfd: vec![0.0; n],
        on_count: vec![0; n],
        off_count: vec![0; n],
        mm_counts: vec![[0u32; MM_BUCKETS]; n],
        saturated: vec![false; n],
    };

    let want_agg = matches!(config.output_mode, OutputMode::Aggregated | OutputMode::Both);
    let want_ot = matches!(config.output_mode, OutputMode::PerOffTarget | OutputMode::Both);
    let mut writers = OtWriters {
        bin: if want_ot && matches!(config.ot_format, OtFormat::Binary | OtFormat::Both) {
            Some(crate::ot_stream::OtBinaryWriter::create(
                &suffix_path(&config.output, "ot.bin"),
                crate::ot_stream::OtHeader {
                    max_mm: config.max_mismatches,
                    cfd_threshold: config.cfd_threshold,
                    n_guides: u64::try_from(n).unwrap_or(u64::MAX),
                },
            )?)
        } else {
            None
        },
        tsv: if want_ot && matches!(config.ot_format, OtFormat::Tsv | OtFormat::Both) {
            Some(crate::ot_stream::OtTsvWriter::create(
                &suffix_path(&config.output, "ot.tsv"),
                cas12a_dual,
            )?)
        } else {
            None
        },
    };

    let default_batch = 50_000usize;
    let batch_cap = config.batch_size.filter(|&b| b > 0).unwrap_or(default_batch);
    // Per-guide saturation caps, applied uniformly across the CPU, GPU-emit, and
    // GPU-CFD-kernel paths. `0` / `0.0` = disabled (the default). off_sum_cap is
    // the off-target CFD sum above which specificity is below --min-specificity.
    let ot_cap = config.max_off_targets;
    let off_sum_cap = if config.min_specificity > 0.0 {
        1.0 / config.min_specificity - 1.0
    } else {
        0.0
    };

    // Optional pre-screen + scan. --threshold t drops guides that have an
    // off-target within t mismatches; the full CFD scan then runs only on the
    // survivors. On GPU one resident index serves both the screen (early-exit
    // scan_screen) and the full pass — no second upload. `abs_map` maps a
    // survivor's scan-slice index back to its absolute guide index.
    let mut dropped = vec![false; n];
    let screened = config.threshold.is_some();
    match scanner {
        ScannerKind::Cpu => {
            if let Some(t) = config.threshold {
                dropped = screen_dropped_cpu(source, guides, t);
            }
            let survivor_abs = survivor_indices(&dropped, screened);
            let survivor_guides: Vec<Guide> =
                survivor_abs.iter().map(|&i| guides[i].clone()).collect();
            let (scan_guides, abs_map): (&[Guide], Option<&[usize]>) = if screened {
                (&survivor_guides, Some(survivor_abs.as_slice()))
            } else {
                (guides, None)
            };
            let mut base = 0usize;
            while base < scan_guides.len() {
                let end = (base + batch_cap).min(scan_guides.len());
                let hits =
                    BinScanner::new(source).scan(&scan_guides[base..end], config.max_mismatches);
                for h in &hits {
                    let rel = base + h.guide_index as usize;
                    let abs = abs_map.map_or(rel, |m| m[rel]);
                    accumulate_hit(
                        abs,
                        scan_guides[rel].site,
                        h.off_target,
                        h.position,
                        h.mismatches,
                        cfd.as_ref(),
                        cas12a.as_ref(),
                        proto_len,
                        config.cfd_threshold,
                        ot_cap,
                        off_sum_cap,
                        &mut acc,
                        &mut writers,
                        source,
                    )?;
                }
                base = end;
            }
        }
        ScannerKind::Gpu => {
            streaming_gpu(
                source,
                guides,
                config.threshold,
                want_ot,
                is_cas12a,
                batch_cap,
                config.max_mismatches,
                config.cfd_threshold,
                ot_cap,
                off_sum_cap,
                cfd.as_ref(),
                cas12a.as_ref(),
                proto_len,
                &mut acc,
                &mut writers,
                &mut dropped,
            )?;
        }
    }

    // ---- write Mode 1 + sidecar ----
    if want_agg {
        let mut agg = crate::ot_stream::AggWriter::create(&config.output, is_cas12a, cas12a_dual)?;
        for (i, rec) in guide_records.iter().enumerate() {
            if dropped[i] {
                // Pre-screened out (off-target within --threshold): omit, unless
                // --keep-dropped asked for a complete record of every guide.
                if config.keep_dropped {
                    agg.write_dropped(&rec.id)?;
                }
                continue;
            }
            // Per-mismatch-bucket counts `[mm0, mm1, …, mm_maxmm]`: mm-0 is the
            // number of *other* perfect-match loci (`on_count - 1`, excluding the
            // guide's own site); the rest come from the scan accumulator.
            let max_mm = usize::from(config.max_mismatches).min(MM_BUCKETS - 1);
            let mut buckets = Vec::with_capacity(max_mm + 1);
            buckets.push(acc.on_count[i].saturating_sub(1));
            buckets.extend_from_slice(&acc.mm_counts[i][1..=max_mm]);

            if is_cas12a {
                // Same on-target-regularized convention as SpCas9: 1/(on_count + Σ)
                // for GuideScan, 1/(1 + Σ) for FlashFry. The streaming accumulator
                // keeps the on-target out of `off_sum` (it lives in `on_count`), so
                // the on-target term must be added back here — exactly what
                // `specificity()` does. (Previously `1/off_sum`, which dropped the
                // on-target term entirely: unbounded, and off-target-free guides
                // landed mid-distribution instead of at the maximum.)
                let tttn = specificity(config.spec_convention, acc.off_sum[i], acc.on_count[i]);
                let tttv = specificity(config.spec_convention, acc.tttv_sum[i], acc.on_count[i]);
                let (tttn2, tttv2) = if cas12a_dual {
                    (
                        Some(specificity(config.spec_convention, acc.off_sum2[i], acc.on_count[i])),
                        Some(specificity(config.spec_convention, acc.tttv_sum2[i], acc.on_count[i])),
                    )
                } else {
                    (None, None)
                };
                agg.write_cas12a(
                    &rec.id,
                    tttn,
                    tttv,
                    tttn2,
                    tttv2,
                    acc.max_cfd[i],
                    acc.off_count[i],
                    &buckets,
                    acc.saturated[i],
                )?;
            } else {
                let spec = specificity(config.spec_convention, acc.off_sum[i], acc.on_count[i]);
                agg.write_cas9(
                    &rec.id,
                    spec,
                    acc.max_cfd[i],
                    acc.off_count[i],
                    &buckets,
                    acc.saturated[i],
                )?;
            }
        }
        agg.finish()?;
    }
    if let Some(w) = writers.bin.take() {
        w.finish()?;
    }
    if let Some(w) = writers.tsv.take() {
        w.finish()?;
    }
    if want_ot {
        let mut side =
            crate::ot_stream::GuidesSidecarWriter::create(&suffix_path(&config.output, "guides.tsv"))?;
        for (i, rec) in guide_records.iter().enumerate() {
            if dropped[i] {
                continue;
            }
            side.write_row(
                u32::try_from(i).unwrap_or(u32::MAX),
                &rec.id,
                &rec.contig_name,
                rec.offset,
                matches!(rec.strand, Strand::Reverse),
                &rec.sequence_for_output,
                acc.on_count[i],
            )?;
        }
        side.finish()?;
    }
    Ok(())
}

/// Pre-screen: flag every guide that has an off-target within `t` mismatches —
/// a 0-mm duplicate (`mm == 0` with count > 1) or any `1..=t`-mm match — using a
/// cheap mm≤t count scan. The candidate-bin set grows steeply with the mismatch
/// budget, so a mm≤2 screen is far cheaper than the full mm≤`max_mismatches`
/// scan; dropping these guides up front means the full CFD scan runs only on the
/// survivors. Returns one flag per guide (`true` = drop, i.e. non-specific).
fn screen_dropped_cpu(source: &dyn BinSource, guides: &[Guide], t: u8) -> Vec<bool> {
    let n = guides.len();
    let mut dropped = vec![false; n];
    let batch = 50_000usize;
    let mut base = 0usize;
    while base < n {
        let end = (base + batch).min(n);
        // BinScanner at max_mm = t returns only mm ≤ t hits.
        let hits = BinScanner::new(source).scan(&guides[base..end], t);
        let mut on = vec![0u32; end - base];
        let mut close = vec![false; end - base];
        for h in &hits {
            let r = h.guide_index as usize;
            if h.mismatches == 0 {
                on[r] += 1;
            } else {
                close[r] = true;
            }
        }
        for r in 0..(end - base) {
            if on[r] > 1 || close[r] {
                dropped[base + r] = true;
            }
        }
        base = end;
    }
    dropped
}

/// Survivor (non-dropped) absolute guide indices, or empty when no screen ran.
fn survivor_indices(dropped: &[bool], screened: bool) -> Vec<usize> {
    if screened {
        (0..dropped.len()).filter(|&i| !dropped[i]).collect()
    } else {
        Vec::new()
    }
}

/// GPU streaming loop: build the resident scanner once, then two-pass each
/// batch — `scan_counts_prefilter` to size the output buffer exactly, then
/// `scan_raw` to emit — and accumulate via [`accumulate_hit`].
#[cfg(feature = "gpu")]
#[allow(clippy::too_many_arguments)]
fn gpu_emit_pass(
    scanner: &crispr_scan_gpu::GpuScanner,
    source: &dyn BinSource,
    guides: &[Guide],
    abs_map: Option<&[usize]>,
    batch_cap: usize,
    max_mm: u8,
    cfd_threshold: f64,
    ot_cap: u32,
    off_sum_cap: f64,
    cfd: Option<&Cfd>,
    cas12a: Option<&Cas12aScorers>,
    proto_len: u32,
    acc: &mut Accum,
    writers: &mut OtWriters,
) -> Result<(), DiscoverError> {
    let mut base = 0usize;
    while base < guides.len() {
        let end = (base + batch_cap).min(guides.len());
        let batch = &guides[base..end];
        let counts = scanner.scan_counts_prefilter(batch, max_mm);
        let predicted: u64 = counts.iter().map(|&c| u64::from(c)).sum();
        let out_cap = u32::try_from(predicted).unwrap_or(u32::MAX);
        // Fast genome-scale hit emission (candidate bins only); identical hit
        // set to brute-force scan_raw, sized exactly by the count pass above.
        let raw = scanner
            .scan_prefilter_raw(batch, max_mm, out_cap)
            .map_err(|e| DiscoverError::Gpu(e.to_string()))?;
        for k in 0..raw.guide.len() {
            let rel = base + raw.guide[k] as usize;
            let abs = abs_map.map_or(rel, |m| m[rel]);
            let e = scanner.entry(raw.entry[k]);
            let mm = u8::try_from(raw.mm[k]).unwrap_or(u8::MAX);
            accumulate_hit(
                abs,
                guides[rel].site,
                e.site(),
                e.position(),
                mm,
                cfd,
                cas12a,
                proto_len,
                cfd_threshold,
                ot_cap,
                off_sum_cap,
                acc,
                writers,
                source,
            )?;
        }
        base = end;
    }
    Ok(())
}

/// Mode-1 fast path (SpCas9 only): fill the per-guide CFD accumulators straight
/// from the GPU CFD-accumulation kernel using the shared resident `scanner`. No
/// hits ever leave the device — output is `O(#guides)` — so this is the path
/// that scales aggregated specificity to a full genome. Guides are processed in
/// large chunks (the chunk only caps the per-launch guide upload).
#[cfg(feature = "gpu")]
fn gpu_cfd_pass(
    scanner: &crispr_scan_gpu::GpuScanner,
    guides: &[Guide],
    abs_map: Option<&[usize]>,
    max_mm: u8,
    ot_cap: u32,
    off_sum_cap: f64,
    acc: &mut Accum,
) -> Result<(), DiscoverError> {
    const CFD_CHUNK: usize = 1 << 20; // ~1M guides/launch
    let mut base = 0usize;
    while base < guides.len() {
        let end = (base + CFD_CHUNK).min(guides.len());
        let aggs = scanner
            .scan_cfd_specificity(&guides[base..end], max_mm, ot_cap, off_sum_cap)
            .map_err(|e| DiscoverError::Gpu(e.to_string()))?;
        for (j, a) in aggs.iter().enumerate() {
            let abs = abs_map.map_or(base + j, |m| m[base + j]);
            acc.off_sum[abs] = a.off_sum;
            acc.max_cfd[abs] = a.max_cfd;
            acc.on_count[abs] = a.on_count;
            acc.off_count[abs] = a.off_count;
            acc.saturated[abs] = a.saturated;
        }
        base = end;
    }
    Ok(())
}

/// GPU streaming driver: build ONE resident index, run the optional early-exit
/// pre-screen (`scan_screen`) and the full scan against it (no second upload),
/// and fill `dropped`. Dispatches to the no-download CFD kernel for Mode-1-only
/// SpCas9, else the hit-emit path.
#[cfg(feature = "gpu")]
#[allow(clippy::too_many_arguments)]
fn streaming_gpu(
    source: &dyn BinSource,
    guides: &[Guide],
    threshold: Option<u8>,
    want_ot: bool,
    is_cas12a: bool,
    batch_cap: usize,
    max_mm: u8,
    cfd_threshold: f64,
    ot_cap: u32,
    off_sum_cap: f64,
    cfd: Option<&Cfd>,
    cas12a: Option<&Cas12aScorers>,
    proto_len: u32,
    acc: &mut Accum,
    writers: &mut OtWriters,
    dropped: &mut [bool],
) -> Result<(), DiscoverError> {
    let scanner =
        crispr_scan_gpu::GpuScanner::new(source).map_err(|e| DiscoverError::Gpu(e.to_string()))?;
    let screened = threshold.is_some();
    if let Some(t) = threshold {
        for (i, f) in scanner.scan_screen(guides, t).into_iter().enumerate() {
            dropped[i] = f != 0;
        }
    }
    let survivor_abs = survivor_indices(dropped, screened);
    let survivor_guides: Vec<Guide> = survivor_abs.iter().map(|&i| guides[i].clone()).collect();
    let (sg, am): (&[Guide], Option<&[usize]>) = if screened {
        (&survivor_guides, Some(survivor_abs.as_slice()))
    } else {
        (guides, None)
    };
    if !want_ot && cfd.is_some() && !is_cas12a {
        gpu_cfd_pass(&scanner, sg, am, max_mm, ot_cap, off_sum_cap, acc)
    } else {
        gpu_emit_pass(
            &scanner, source, sg, am, batch_cap, max_mm, cfd_threshold, ot_cap, off_sum_cap, cfd,
            cas12a, proto_len, acc, writers,
        )
    }
}

/// Stub for non-GPU builds: streaming with `--scanner gpu` needs the feature.
#[cfg(not(feature = "gpu"))]
#[allow(clippy::too_many_arguments)]
fn streaming_gpu(
    _source: &dyn BinSource,
    _guides: &[Guide],
    _threshold: Option<u8>,
    _want_ot: bool,
    _is_cas12a: bool,
    _batch_cap: usize,
    _max_mm: u8,
    _cfd_threshold: f64,
    _ot_cap: u32,
    _off_sum_cap: f64,
    _cfd: Option<&Cfd>,
    _cas12a: Option<&Cas12aScorers>,
    _proto_len: u32,
    _acc: &mut Accum,
    _writers: &mut OtWriters,
    _dropped: &mut [bool],
) -> Result<(), DiscoverError> {
    Err(DiscoverError::Gpu(
        "streaming --output-mode with --scanner gpu requires a binary built \
         with `cargo build --features gpu`"
            .to_string(),
    ))
}

/// Side-channel metadata for one input guide. The variant carries the
/// extra fields the GuideScan2 CSV writer needs to echo back.
#[derive(Debug, Clone)]
struct GuideRecord {
    /// Output identifier. For FASTA input: `"<contig>:<offset>:<strand>"`
    /// synthesized from the discovered site. For kmers CSV input:
    /// the `id` column verbatim.
    id: String,
    /// Bit-encoded scan target (protospacer + PAM).
    site: Site,
    /// Sequence column as it should appear in CSV output. For FASTA
    /// input: the decoded scan target. For kmers CSV input: the literal
    /// `sequence + pam` string from the input (preserves the PAM's `N`
    /// characters etc.).
    sequence_for_output: String,
    /// On-target contig name.
    contig_name: String,
    /// On-target offset. For FASTA input this is 0-indexed (Position
    /// convention); for kmers CSV input it's whatever the user passed.
    offset: u32,
    /// On-target strand.
    strand: Strand,
}

#[derive(Debug, Clone, Copy)]
struct GroupAgg {
    mismatches: u8,
    count: u32,
}

fn collect_from_fasta(
    path: &Path,
    source: &dyn BinSource,
) -> Result<Vec<GuideRecord>, DiscoverError> {
    let enzyme = source.enzyme().clone();
    let finder = SiteFinder::new(enzyme.clone());
    let query_contigs = read_fasta(path)?;
    if query_contigs.is_empty() {
        return Err(DiscoverError::QueryEmpty);
    }
    let total_len = usize::from(enzyme.total_scan_len());

    let mut records: Vec<GuideRecord> = Vec::new();
    for contig in &query_contigs {
        let mut sink: Vec<(Site, Position)> = Vec::new();
        finder.scan(&contig.sequence, 0, &mut sink);
        for (site, pos) in sink {
            let strand_char = pos.strand.as_char();
            records.push(GuideRecord {
                id: format!(
                    "{name}:{off}:{strand_char}",
                    name = contig.name,
                    off = pos.offset
                ),
                site,
                sequence_for_output: site.decode_ascii(total_len),
                contig_name: contig.name.clone(),
                offset: pos.offset,
                strand: pos.strand,
            });
        }
    }
    Ok(records)
}

fn collect_from_kmers_csv(
    path: &Path,
    source: &dyn BinSource,
) -> Result<Vec<GuideRecord>, DiscoverError> {
    let enzyme = source.enzyme().clone();
    let entries = kmers_csv::read_from_path(path, &enzyme)?;
    let records = entries
        .into_iter()
        .map(
            |KmerEntry {
                 id,
                 sequence,
                 pam,
                 chromosome,
                 position,
                 sense,
                 site,
             }| {
                let strand = match sense {
                    '+' => Strand::Forward,
                    _ => Strand::Reverse,
                };
                let mut sequence_for_output = sequence;
                sequence_for_output.push_str(&pam);
                GuideRecord {
                    id,
                    site,
                    sequence_for_output,
                    contig_name: chromosome,
                    offset: u32::try_from(position).unwrap_or(u32::MAX),
                    strand,
                }
            },
        )
        .collect();
    Ok(records)
}

/// Group raw hits by `guide_index`, sorted within each group by
/// `(contig_id, offset, strand, mismatches)` for stable output. Returns
/// a vec indexed by guide; an empty inner vec means that guide had no
/// hits at all (not even the on-target — possible in pathological test
/// fixtures).
fn index_hits_by_guide(hits: &[Hit]) -> Vec<Vec<Hit>> {
    let mut by_guide: Vec<Vec<Hit>> = Vec::new();
    for hit in hits {
        let idx = hit.guide_index as usize;
        if by_guide.len() <= idx {
            by_guide.resize_with(idx + 1, Vec::new);
        }
        by_guide[idx].push(*hit);
    }
    for v in &mut by_guide {
        v.sort_by(|a, b| {
            a.position
                .contig_id
                .cmp(&b.position.contig_id)
                .then(a.position.offset.cmp(&b.position.offset))
                .then(a.position.strand.as_char().cmp(&b.position.strand.as_char()))
                .then(a.mismatches.cmp(&b.mismatches))
        });
    }
    by_guide
}

#[allow(clippy::too_many_arguments)]
fn write_csv_output(
    path: &Path,
    source: &dyn BinSource,
    guide_records: &[GuideRecord],
    hits_by_guide: &[Vec<Hit>],
    cfd_results: &[Option<CfdResult>],
    cas12a_results: &[Option<Cas12aResult>],
    keep_guide: &[bool],
) -> Result<(), DiscoverError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let mut rows: Vec<CsvRow<'_>> = Vec::new();
    let empty_hits: Vec<Hit> = Vec::new();
    for (idx, rec) in guide_records.iter().enumerate() {
        if !keep_guide[idx] {
            continue;
        }
        // Specificity for the CSV's `specificity` column. SpCas9 takes
        // precedence when both scorers were requested in the same run
        // (uncommon, since they target different enzymes); when only
        // Cas12a is populated we surface its TTTN specificity. Falls
        // back to 1.0 when no score metric was requested.
        let spec = match (cfd_results[idx], cas12a_results[idx]) {
            (Some(r), _) => r.specificity,
            (None, Some(r)) => r.tttn_specificity,
            (None, None) => 1.0,
        };
        let hits_for_guide: &[Hit] = hits_by_guide
            .get(idx)
            .map_or(empty_hits.as_slice(), Vec::as_slice);
        for hit in hits_for_guide {
            let match_chrm = source
                .contig_name(hit.position.contig_id)
                .unwrap_or("");
            rows.push(CsvRow {
                id: &rec.id,
                sequence: &rec.sequence_for_output,
                match_chrm,
                match_position: hit.position.offset,
                match_strand: hit.position.strand.as_char(),
                match_distance: hit.mismatches,
                specificity: spec,
            });
        }
    }

    csv_out::write_rows(&mut writer, &rows)?;
    writer.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_tsv_output(
    path: &Path,
    guide_records: &[GuideRecord],
    grouped: &HashMap<usize, HashMap<Site, GroupAgg>>,
    cfd_results: &[Option<CfdResult>],
    cas12a_results: &[Option<Cas12aResult>],
    keep_guide: &[bool],
    total_len: usize,
) -> Result<(), DiscoverError> {
    let mut rows: Vec<DiscoverRow> = Vec::with_capacity(guide_records.len());
    for (idx, rec) in guide_records.iter().enumerate() {
        if !keep_guide[idx] {
            continue;
        }
        let group = grouped.get(&idx);
        let all_offs: Vec<(Site, u8, u32)> = group
            .map(|m| {
                m.iter()
                    .map(|(s, agg)| (*s, agg.mismatches, agg.count))
                    .collect()
            })
            .unwrap_or_default();
        let ot_count: usize = all_offs.iter().map(|(_, _, c)| *c as usize).sum();
        let ot_sequences = all_offs
            .iter()
            .map(|(site, mm, count)| format!("{}_{}_{}", site.decode_ascii(total_len), count, mm))
            .collect::<Vec<_>>()
            .join(",");
        let cfd = cfd_results[idx].map(|r| CfdColumns {
            cfd_max: r.max_cfd,
            cfd_specificity: r.specificity,
        });
        let cas12a = cas12a_results[idx].map(|r| Cas12aColumns {
            cas12a_max: r.max_cfd,
            cas12a_spec_tttn: r.tttn_specificity,
            cas12a_spec_tttv: r.tttv_specificity,
        });
        rows.push(DiscoverRow {
            contig: rec.contig_name.clone(),
            start: rec.offset,
            stop: rec.offset + u32::try_from(total_len).expect("len fits in u32"),
            target: rec.site.decode_ascii(total_len),
            strand: rec.strand,
            ot_count,
            ot_sequences,
            cfd,
            cas12a,
        });
    }
    rows.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.start.cmp(&b.start))
            .then(a.strand.as_char().cmp(&b.strand.as_char()))
    });
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_tsv_rows(&mut writer, &rows)?;
    writer.flush()?;
    Ok(())
}

/// Compute the sidecar `<output>.detail.tsv` path used by
/// `--format both`. Always appends the suffix verbatim regardless of the
/// original extension — the goal is a predictable name, not a clever one.
fn sidecar_tsv_path(csv_path: &Path) -> PathBuf {
    let mut s: OsString = csv_path.as_os_str().to_owned();
    s.push(".detail.tsv");
    PathBuf::from(s)
}

/// Map a CLI enzyme name to an `Enzyme` preset. Names are case-insensitive
/// and match FlashFry's `--enzyme` argument shape.
#[must_use]
pub fn enzyme_from_name(name: &str) -> Option<crispr_encoding::Enzyme> {
    use crispr_encoding::Enzyme;
    match name.to_ascii_lowercase().as_str() {
        "spcas9ngg" | "spcas9-ngg" => Some(Enzyme::spcas9_ngg()),
        "spcas9nag" | "spcas9-nag" => Some(Enzyme::spcas9_nag()),
        "spcas9" => Some(Enzyme::spcas9_ngg_or_nag()),
        "cpf1" | "cas12a" => Some(Enzyme::cpf1_tttn()),
        _ => None,
    }
}

/// Build an `Enzyme` directly from a PAM string + protospacer length + PAM
/// orientation — the flexible, preset-free index interface. `pam` is an IUPAC
/// motif (e.g. `NGG`, `TTTV`); `five_prime` selects a Cas12a-style 5' PAM, else
/// a 3' SpCas9-style PAM. Off-target sites and mismatch geometry derive from
/// these three values alone, so any PAM/length/orientation works without a
/// named preset (the `.crot` index stores them, so `enumerate` needs no enzyme
/// argument).
///
/// # Errors
/// Returns a message if `pam` is empty, contains a non-IUPAC base, or if
/// `pam.len() + protospacer_len` exceeds the 28-base site limit.
pub fn enzyme_from_pam(pam: &str, protospacer_len: u8, five_prime: bool) -> Result<crispr_encoding::Enzyme, String> {
    use crispr_encoding::{Enzyme, IupacCode, PamSide, MAX_SITE_LEN};
    if pam.is_empty() {
        return Err("--pam must be a non-empty IUPAC string (e.g. NGG, TTTV)".to_string());
    }
    let codes: Vec<IupacCode> = pam
        .bytes()
        .map(|b| {
            IupacCode::from_ascii(b.to_ascii_uppercase())
                .ok_or_else(|| format!("invalid PAM base '{}' (expected IUPAC: A C G T R Y S W K M B D H V N)", b as char))
        })
        .collect::<Result<_, _>>()?;
    let total = codes.len() + usize::from(protospacer_len);
    if total > MAX_SITE_LEN {
        return Err(format!(
            "PAM ({}) + protospacer ({protospacer_len}) = {total} exceeds the {MAX_SITE_LEN}-base site limit",
            codes.len()
        ));
    }
    let side = if five_prime { PamSide::FivePrime } else { PamSide::ThreePrime };
    let name = format!(
        "PAM-{}-{}-{}nt",
        pam.to_ascii_uppercase(),
        if five_prime { "5p" } else { "3p" },
        protospacer_len
    );
    Ok(Enzyme::from_iupac(&name, &codes, side, protospacer_len))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_dir() -> PathBuf {
        let base = std::env::temp_dir();
        let path = base.join(format!(
            "crispr-ots-discover-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        ));
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    fn write_file(path: &Path, body: &str) {
        let mut f = File::create(path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
    }

    #[test]
    fn enzyme_lookup_table() {
        assert!(enzyme_from_name("spcas9ngg").is_some());
        assert!(enzyme_from_name("SPCAS9-NGG").is_some());
        assert!(enzyme_from_name("cpf1").is_some());
        assert!(enzyme_from_name("CAS12a").is_some());
        assert!(enzyme_from_name("unknown").is_none());
    }

    #[test]
    fn enzyme_from_pam_builds_arbitrary_enzymes() {
        use crispr_encoding::PamSide;
        // SpCas9-style: NGG, 20-nt, 3' PAM.
        let cas9 = enzyme_from_pam("NGG", 20, false).expect("NGG parses");
        assert_eq!(cas9.pam_len, 3);
        assert_eq!(cas9.protospacer_len, 20);
        assert!(matches!(cas9.pam_side, PamSide::ThreePrime));
        // Cas12a-style: TTTV, 23-nt, 5' PAM (case-insensitive).
        let cas12a = enzyme_from_pam("tttv", 23, true).expect("TTTV parses");
        assert_eq!(cas12a.pam_len, 4);
        assert_eq!(cas12a.protospacer_len, 23);
        assert!(matches!(cas12a.pam_side, PamSide::FivePrime));
        // Errors: empty, bad base, over the 28-base site limit.
        assert!(enzyme_from_pam("", 20, false).is_err());
        assert!(enzyme_from_pam("NXG", 20, false).is_err());
        assert!(enzyme_from_pam("NGG", 30, false).is_err());
    }

    #[test]
    fn sidecar_path_appends_suffix() {
        assert_eq!(
            sidecar_tsv_path(Path::new("/tmp/out.csv")),
            PathBuf::from("/tmp/out.csv.detail.tsv")
        );
    }

    #[test]
    fn single_guide_end_to_end_tsv() {
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let queries = dir.join("query.fa");
        let output = dir.join("out.tsv");
        write_file(&reference, ">chrA\nAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTT\n");
        write_file(&queries, ">q\nAAAAAAAAAAAAAAAAAAAAAGG\n");

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries),
                max_mismatches: 0,
                scores: vec![ScoreMetric::Cfd],
                output: output.clone(),
                format: OutputFormat::Tsv,
                spec_convention: SpecConvention::Flashfry,
                threshold: None,
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();

        let tsv = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 2);
        let fields: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(fields[0], "q");
        assert_eq!(fields[3], "AAAAAAAAAAAAAAAAAAAAAGG");
        assert_eq!(fields[4], "+");
        // FlashFry-style otCount: 1 for the on-target hit.
        assert_eq!(fields[5], "1");
        assert_eq!(fields[7], "0");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn single_guide_end_to_end_csv() {
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let queries = dir.join("query.fa");
        let output = dir.join("out.csv");
        write_file(&reference, ">chrA\nAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTT\n");
        write_file(&queries, ">q\nAAAAAAAAAAAAAAAAAAAAAGG\n");

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries),
                max_mismatches: 0,
                scores: vec![ScoreMetric::Cfd],
                output: output.clone(),
                format: OutputFormat::Csv,
                spec_convention: SpecConvention::Guidescan,
                threshold: None,
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();

        let csv = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        // Header + 1 on-target row (the only hit at 0 mm).
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            "id,sequence,match_chrm,match_position,match_strand,match_distance,specificity"
        );
        let fields: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(fields[0], "q:0:+");
        assert_eq!(fields[1], "AAAAAAAAAAAAAAAAAAAAAGG");
        assert_eq!(fields[2], "chrA");
        assert_eq!(fields[3], "0");
        assert_eq!(fields[4], "+");
        assert_eq!(fields[5], "0");
        // GuideScan convention with one on-target: 1 / 1.0 = 1.0
        assert!(fields[6].starts_with("1.0"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn kmers_csv_input_drives_csv_output() {
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let kmers = dir.join("kmers.csv");
        let output = dir.join("out.csv");
        write_file(&reference, ">chrA\nAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTT\n");
        write_file(
            &kmers,
            "id,sequence,pam,chromosome,position,sense\n\
             test:0:+,AAAAAAAAAAAAAAAAAAAA,NGG,chrA,0,+\n",
        );

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::KmersCsv(kmers),
                max_mismatches: 0,
                scores: vec![ScoreMetric::Cfd],
                output: output.clone(),
                format: OutputFormat::Csv,
                spec_convention: SpecConvention::Guidescan,
                threshold: None,
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();

        let csv = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 2);
        // id echoes the input verbatim, and sequence preserves NGG.
        let fields: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(fields[0], "test:0:+");
        assert_eq!(fields[1], "AAAAAAAAAAAAAAAAAAAANGG");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn threshold_drops_guides_with_close_offtargets() {
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        // Reference has two near-identical NGG sites separated by 30 bp;
        // they're 1 mismatch apart. With --threshold 1, the guide should
        // be dropped from the output.
        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let queries = dir.join("q.fa");
        let output = dir.join("out.tsv");
        write_file(
            &reference,
            ">chrA\nAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAACAGGTTTTTTTT\n",
        );
        write_file(&queries, ">q\nAAAAAAAAAAAAAAAAAAAAAGG\n");

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries),
                max_mismatches: 4,
                scores: vec![ScoreMetric::Cfd],
                output: output.clone(),
                format: OutputFormat::Tsv,
                spec_convention: SpecConvention::Flashfry,
                threshold: Some(1),
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();

        let tsv = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        // header only — both query sites are dropped because each one's
        // closest off-target is the other (1 mm).
        assert_eq!(lines.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn threshold_drops_multimapping_on_target() {
        // The guide's exact sequence appears at TWO genomic positions
        // (both mm = 0). guidescan2 treats this as a hit at distance 0
        // from another PAM-adjacent site and drops the guide under any
        // non-negative `--threshold`. Verify we match that semantic.
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let queries = dir.join("q.fa");
        let output = dir.join("out.tsv");
        // Same NGG-flanked 20-mer twice, 30 bp apart, no mismatches.
        write_file(
            &reference,
            ">chrA\nAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAGGTTTTTTTT\n",
        );
        write_file(&queries, ">q\nAAAAAAAAAAAAAAAAAAAAAGG\n");

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries),
                max_mismatches: 4,
                scores: vec![ScoreMetric::Cfd],
                output: output.clone(),
                format: OutputFormat::Tsv,
                spec_convention: SpecConvention::Flashfry,
                threshold: Some(0),
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();

        let tsv = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        // header only — both discoveries of the guide are dropped because
        // the same sequence is found at two genomic positions.
        assert_eq!(lines.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn max_per_bin_caps_distinct_off_target_sequences() {
        // Reference contains the query at offset 0 (the on-target) plus
        // three distinct 1-mismatch off-targets at positions 30, 60, 90.
        // With `max_per_bin = Some(2)` the third should be dropped from
        // the per-(guide, mm=1) bin. The on-target (mm=0) stays.
        use crate::build::{build_table_in_memory, BuildConfig};
        use crispr_encoding::Enzyme;

        let dir = fixture_dir();
        let reference = dir.join("ref.fa");
        let queries = dir.join("q.fa");
        let output = dir.join("out.tsv");
        // Each protospacer is 20 As + AGG = 23 chars. Variations at
        // position 1 keep the PAM (`AGG`) intact and create a single
        // mismatch each.
        let mut body = String::from(">chrA\n");
        body.push_str("AAAAAAAAAAAAAAAAAAAAAGG"); // on-target
        body.push_str(&"T".repeat(7)); // spacer to next NGG
        body.push_str("CAAAAAAAAAAAAAAAAAAAAGG"); // 1mm
        body.push_str(&"T".repeat(7));
        body.push_str("GAAAAAAAAAAAAAAAAAAAAGG"); // 1mm
        body.push_str(&"T".repeat(7));
        body.push_str("TAAAAAAAAAAAAAAAAAAAAGG"); // 1mm (dropped at cap)
        body.push('\n');
        write_file(&reference, &body);
        write_file(&queries, ">q\nAAAAAAAAAAAAAAAAAAAAAGG\n");

        let table = build_table_in_memory(&BuildConfig {
            reference,
            enzyme: Enzyme::spcas9_ngg(),
            output: dir.join("ignored.crot"),
            bin_width: None,
        })
        .unwrap();

        // Reference run: cap disabled → otCount counts all four hits
        // (1 on-target + 3 off-targets at mm=1).
        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries.clone()),
                max_mismatches: 4,
                scores: vec![],
                output: output.clone(),
                format: OutputFormat::Tsv,
                spec_convention: SpecConvention::Flashfry,
                threshold: None,
                max_per_bin: None,
                ..Default::default()
            },
        )
        .unwrap();
        let unbounded = std::fs::read_to_string(&output).unwrap();
        let unbounded_otcount: usize = unbounded
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(5)
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(unbounded_otcount, 4);

        // Capped run: keep at most 2 distinct sequences per (guide, mm)
        // bin. Mm=0 keeps the 1 on-target; mm=1 keeps 2 of 3 → total 3.
        run_discover(
            &table,
            &DiscoverConfig {
                input: DiscoverInput::QueryFasta(queries),
                max_mismatches: 4,
                scores: vec![],
                output: output.clone(),
                format: OutputFormat::Tsv,
                spec_convention: SpecConvention::Flashfry,
                threshold: None,
                max_per_bin: Some(2),
                ..Default::default()
            },
        )
        .unwrap();
        let capped = std::fs::read_to_string(&output).unwrap();
        let capped_otcount: usize = capped
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .nth(5)
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(capped_otcount, 3);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
