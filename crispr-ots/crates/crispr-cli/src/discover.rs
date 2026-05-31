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
use crate::tsv::{write_rows as write_tsv_rows, CfdColumns, DiscoverRow};

/// Scoring metrics that can be applied during discover.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreMetric {
    Cfd,
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
    /// If `Some(t)`, drop any guide whose nearest *off*-target (mm > 0)
    /// has mismatch ≤ t. Matches `guidescan enumerate --threshold`.
    /// `None` disables filtering.
    pub threshold: Option<u8>,
}

/// Errors surfaced by the discover pipeline. Downstream layers wrap with
/// `anyhow` for friendly CLI messages.
#[derive(Debug)]
pub enum DiscoverError {
    Io(std::io::Error),
    QueryEmpty,
    KmersCsv(kmers_csv::KmerCsvError),
}

impl std::fmt::Display for DiscoverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::QueryEmpty => write!(f, "query FASTA contained no contigs"),
            Self::KmersCsv(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for DiscoverError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::QueryEmpty => None,
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
#[allow(clippy::too_many_lines)]
pub fn run_discover(source: &dyn BinSource, config: &DiscoverConfig) -> Result<(), DiscoverError> {
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
    let scanner = BinScanner::new(source);
    let t_scan = std::time::Instant::now();
    let hits = scanner.scan(&guides, config.max_mismatches);
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
    let mut grouped: HashMap<usize, HashMap<Site, GroupAgg>> = HashMap::new();
    for hit in &hits {
        let entry = grouped
            .entry(hit.guide_index as usize)
            .or_default()
            .entry(hit.off_target)
            .or_insert(GroupAgg {
                mismatches: hit.mismatches,
                count: 0,
            });
        entry.count += 1;
        debug_assert_eq!(entry.mismatches, hit.mismatches);
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

    // ---- Compute the per-guide CFD result + apply --threshold filter. ----
    let mut keep_guide: Vec<bool> = vec![true; guide_records.len()];
    let mut cfd_results: Vec<Option<CfdResult>> = vec![None; guide_records.len()];
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
            // Drop the guide if any off-target row (mm > 0) has mm <= t.
            // The on-target rows (mm == 0) never trigger the filter.
            let close_off = all_offs
                .iter()
                .any(|(_, mm, _)| *mm > 0 && *mm <= t);
            if close_off {
                keep_guide[idx] = false;
                continue;
            }
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
                &keep_guide,
            )?;
        }
        OutputFormat::Tsv => {
            write_tsv_output(
                &config.output,
                &guide_records,
                &grouped,
                &cfd_results,
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
                &keep_guide,
            )?;
            let tsv_path = sidecar_tsv_path(&config.output);
            write_tsv_output(
                &tsv_path,
                &guide_records,
                &grouped,
                &cfd_results,
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
        let spec = cfd_results[idx].map_or(1.0, |r| r.specificity);
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

fn write_tsv_output(
    path: &Path,
    guide_records: &[GuideRecord],
    grouped: &HashMap<usize, HashMap<Site, GroupAgg>>,
    cfd_results: &[Option<CfdResult>],
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
        rows.push(DiscoverRow {
            contig: rec.contig_name.clone(),
            start: rec.offset,
            stop: rec.offset + u32::try_from(total_len).expect("len fits in u32"),
            target: rec.site.decode_ascii(total_len),
            strand: rec.strand,
            ot_count,
            ot_sequences,
            cfd,
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
}
