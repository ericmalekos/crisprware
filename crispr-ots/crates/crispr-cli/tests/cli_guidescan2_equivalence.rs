//! Cross-tool equivalence test against GuideScan2 on the 1000-gRNA chr22
//! fixture. Auto-skips unless both fixtures are present.
//!
//! What this validates and why
//! ===========================
//!
//! GuideScan2 uses an FM-index over the raw genome and recursive
//! backward-search to enumerate off-targets, which is a totally different
//! algorithm from our (and FlashFry's) k-mer-bin + bit-encoded scan.
//! Matching it cross-validates that our search finds the *same set of
//! genomic positions* a fundamentally different algorithm finds. That's a
//! strong correctness signal independent of FlashFry.
//!
//! Two semantic differences from FlashFry/ours that the test handles:
//!   1. GuideScan2 enumerates one row per genomic position. We collapse
//!      identical off-target *sequences* into a single row carrying a
//!      `count`. So we compare *total position count per guide*, not
//!      *distinct sequence count*.
//!   2. CFD specificity convention. FlashFry: `1 / (1 + Σ cfd × count)`
//!      with the on-target as a singular +1. GuideScan2: `1 / Σ_all_rows`
//!      with one CFD-1.0 entry per on-target position. We support both
//!      via [`SpecConvention`]; this test passes
//!      [`SpecConvention::Guidescan`] and expects bit-identical
//!      specificity (to GS2's print precision).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crispr_cli::{
    build_table_in_memory, run_discover, BuildConfig, DiscoverConfig, DiscoverInput, OutputFormat,
    ScoreMetric, SpecConvention,
};
use crispr_encoding::Enzyme;

const CHR22_PATH: &str = "/tmp/ffry-quickstart/chr22.fa.gz";
const QUERY_FASTA: &str = "/tmp/ffry-quickstart/random_1000.fasta";
const GS2_CSV_OUT: &str = "/tmp/ffry-quickstart/random_1000.gs2.ngg.csv.out";
const MAX_MISMATCHES: u8 = 4;

/// GuideScan2 prints specificity to 6 significant digits; ~1e-6 absolute
/// is the formatting floor. We give ourselves an order of magnitude of
/// headroom for floating-point rounding noise from multiplication order.
const SPEC_TOL: f64 = 1e-5;

fn fixture_present() -> bool {
    Path::new(CHR22_PATH).exists()
        && Path::new(QUERY_FASTA).exists()
        && Path::new(GS2_CSV_OUT).exists()
}

#[derive(Debug, Default)]
struct Gs2Guide {
    total_rows: u64,
    specificity: f64,
}

/// Parse GuideScan2's CSV: one row per (guide, off-target). All rows for
/// a given guide share the same specificity, so we record it once.
fn parse_guidescan2(path: &Path) -> HashMap<String, Gs2Guide> {
    let text = std::fs::read_to_string(path).expect("read GS2 output");
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let cols: Vec<&str> = header.split(',').collect();
    let id_idx = cols.iter().position(|c| *c == "id").expect("id col");
    let spec_idx = cols
        .iter()
        .position(|c| *c == "specificity")
        .expect("specificity col");
    let mut out: HashMap<String, Gs2Guide> = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        let gid = fields[id_idx].to_string();
        let spec: f64 = fields[spec_idx].parse().expect("spec float");
        let entry = out.entry(gid).or_default();
        entry.total_rows += 1;
        entry.specificity = spec; // overwrite — same on every row of a guide
    }
    out
}

/// Parse just the `otCount` column out of our TSV output (per-guide
/// FlashFry-style format). Used for the position-count cross-check;
/// specificity comes from the CSV path.
fn parse_ours_tsv_otcount(path: &Path) -> HashMap<String, u64> {
    let text = std::fs::read_to_string(path).expect("read ours TSV");
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let cols: Vec<&str> = header.split('\t').collect();
    let tgt_idx = cols.iter().position(|c| *c == "target").expect("target");
    let ot_idx = cols.iter().position(|c| *c == "otCount").expect("otCount");
    let mut out: HashMap<String, u64> = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let target = fields[tgt_idx].to_string();
        out.insert(target, fields[ot_idx].parse().unwrap_or(0));
    }
    out
}

/// Read the query FASTA to build `gs2_id -> 23-bp target` map.
fn read_id_to_target(path: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(path).expect("read query FASTA");
    let mut out = HashMap::new();
    let mut cur: Option<String> = None;
    for line in text.lines() {
        if let Some(stripped) = line.strip_prefix('>') {
            cur = Some(stripped.split_whitespace().next().unwrap_or("").to_string());
        } else if let Some(id) = cur.take() {
            out.insert(id, line.trim().to_string());
        }
    }
    out
}

#[test]
#[allow(clippy::too_many_lines)] // integration test reads linearly
fn match_set_and_specificity_match_guidescan2_under_gs2_convention() {
    if !fixture_present() {
        eprintln!(
            "Skipping GuideScan2 equivalence test: fixture missing.\n\
             Setup:\n\
             \n  1. Build a GuideScan2 chr22 index (one-time, ~13 s):\n     \
                cd /tmp/ffry-quickstart\n     \
                gunzip -k chr22.fa.gz\n     \
                /home/eric/miniconda3/envs/crisprlib/bin/guidescan index --index chr22.gs2 chr22.fa\n\
             \n  2. Convert random_1000.fasta to GuideScan2's CSV input format:\n     \
                # see python helper in NOTES.md\n\
             \n  3. Run enumerate with NGG alt-PAMs to match our search semantics:\n     \
                /home/eric/miniconda3/envs/crisprlib/bin/guidescan enumerate \\\n       \
                    --kmers-file random_1000.gs2.csv --mismatches 4 \\\n       \
                    --alt-pam AGG --alt-pam CGG --alt-pam GGG --alt-pam TGG \\\n       \
                    --format csv --mode complete \\\n       \
                    --output random_1000.gs2.ngg.csv.out chr22.gs2"
        );
        return;
    }

    // Drive our pipeline twice: once for the TSV (carries otCount for the
    // position-count check) and once for the CSV (carries the per-row
    // specificity under the GS2 convention). Both runs share one
    // in-memory bin table.
    let tmp = std::env::temp_dir().join(format!(
        "crispr-ots-gs2-test-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos())
    ));
    std::fs::create_dir_all(&tmp).expect("temp dir");
    let tsv_path = tmp.join("ours.tsv");
    let csv_path = tmp.join("ours.csv");

    let table = build_table_in_memory(&BuildConfig {
        reference: PathBuf::from(CHR22_PATH),
        enzyme: Enzyme::spcas9_ngg(),
        output: PathBuf::from("/dev/null"),
        bin_width: None,
    })
    .expect("build chr22 table");

    run_discover(
        &table,
        &DiscoverConfig {
            input: DiscoverInput::QueryFasta(PathBuf::from(QUERY_FASTA)),
            max_mismatches: MAX_MISMATCHES,
            scores: vec![ScoreMetric::Cfd],
            output: tsv_path.clone(),
            format: OutputFormat::Tsv,
            spec_convention: SpecConvention::Flashfry,
            threshold: None,
            max_per_bin: None,
            ..Default::default()
        },
    )
    .expect("discover (TSV) succeeds");

    run_discover(
        &table,
        &DiscoverConfig {
            input: DiscoverInput::QueryFasta(PathBuf::from(QUERY_FASTA)),
            max_mismatches: MAX_MISMATCHES,
            scores: vec![ScoreMetric::Cfd],
            output: csv_path.clone(),
            format: OutputFormat::Csv,
            spec_convention: SpecConvention::Guidescan,
            threshold: None,
            max_per_bin: None,
            ..Default::default()
        },
    )
    .expect("discover (CSV) succeeds");

    let id_to_target = read_id_to_target(Path::new(QUERY_FASTA));
    let gs2 = parse_guidescan2(Path::new(GS2_CSV_OUT));
    let ours_otcount = parse_ours_tsv_otcount(&tsv_path);
    let ours_spec = parse_ours_csv(&csv_path);

    let mut common = 0_usize;
    let mut pos_mismatches = 0_usize;
    let mut worst_spec_drift = 0.0_f64;
    let mut spec_drift_count = 0_usize;
    let mut examples: Vec<String> = Vec::new();

    for (gid, target) in &id_to_target {
        let Some(g) = gs2.get(gid) else { continue };
        let Some(o_ot) = ours_otcount.get(target) else {
            continue;
        };
        let Some(o_spec) = ours_spec.get(gid) else {
            // The GS2-style CSV is keyed by id, and our id matches the
            // input FASTA contig name (e.g. `chr22_site_42:offset:+`).
            // Skip if for some reason the row didn't make it.
            continue;
        };
        common += 1;

        // ---- Equivalence check 1: total off-target positions per guide. ----
        if g.total_rows != *o_ot {
            pos_mismatches += 1;
            if examples.len() < 3 {
                examples.push(format!(
                    "  position-count drift {gid}: GS2={}, ours={}",
                    g.total_rows, o_ot
                ));
            }
        }

        // ---- Equivalence check 2: bit-identical specificity. ----
        // Both tools now compute `1 / Σ_all_rows(cfd)`. Differences are
        // limited to floating-point rounding plus GS2's 6-sigfig print
        // precision.
        let diff = (o_spec - g.specificity).abs();
        if diff > worst_spec_drift {
            worst_spec_drift = diff;
        }
        if diff > SPEC_TOL {
            spec_drift_count += 1;
            if examples.len() < 6 {
                examples.push(format!(
                    "  spec drift {gid}: gs2={:.6}, ours={:.6}, Δ={diff:e}",
                    g.specificity, o_spec
                ));
            }
        }
    }

    eprintln!("GS2 equivalence on {common} common guides:");
    eprintln!("  position-count mismatches: {pos_mismatches}");
    eprintln!("  max |Δ spec|: {worst_spec_drift:e}");
    eprintln!("  spec mismatches > {SPEC_TOL:e}: {spec_drift_count}");
    for line in &examples {
        eprintln!("{line}");
    }

    assert_eq!(
        pos_mismatches, 0,
        "GuideScan2 and crispr-ots must find the same off-target positions per guide"
    );
    assert_eq!(
        spec_drift_count, 0,
        "Under --spec-convention guidescan, our specificity must match GS2's within {SPEC_TOL:e}",
    );

    let _ = std::fs::remove_dir_all(&tmp);
}

/// Parse the per-guide specificity out of our CSV output. CSV has one
/// row per (guide, off-target hit); specificity is repeated on every row
/// so taking the first occurrence is enough.
fn parse_ours_csv(path: &Path) -> HashMap<String, f64> {
    let text = std::fs::read_to_string(path).expect("read ours CSV");
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let cols: Vec<&str> = header.split(',').collect();
    let id_idx = cols.iter().position(|c| *c == "id").expect("id col");
    let spec_idx = cols
        .iter()
        .position(|c| *c == "specificity")
        .expect("specificity col");
    let mut out: HashMap<String, f64> = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() <= spec_idx {
            continue;
        }
        let gid = fields[id_idx].to_string();
        let spec: f64 = fields[spec_idx].parse().unwrap_or(0.0);
        out.entry(gid).or_insert(spec);
    }
    out
}
