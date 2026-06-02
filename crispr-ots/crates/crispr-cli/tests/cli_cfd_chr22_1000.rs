//! Larger-scale FlashFry-equivalence test on 1000 random chr22 NGG sites.
//!
//! The 3-guide EMX1 fixture covers the happy path but missed two
//! semantic differences that only surface at scale:
//!   1. CFD specificity weights off-targets by the number of *genomic
//!      positions* each unique OT sequence appears at. Our pipeline used
//!      to weight by 1.
//!   2. FlashFry's `otCount` column is the total position count
//!      including on-target hits and multi-mappers; we used to report a
//!      simpler "distinct OT sequences with mm > 0" count.
//!
//! Both have been corrected. This test guards against regression.
//!
//! Auto-skips unless the fixture is set up at `/tmp/ffry-quickstart`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crispr_cli::{
    build_table_in_memory, run_discover, BuildConfig, DiscoverConfig, DiscoverInput, OutputFormat,
    ScoreMetric, SpecConvention,
};
use crispr_encoding::Enzyme;

const CHR22_PATH: &str = "/tmp/ffry-quickstart/chr22.fa.gz";
const QUERY_FASTA: &str = "/tmp/ffry-quickstart/random_1000.fasta";
const FLASHFRY_SCORED: &str = "/tmp/ffry-quickstart/random_1000.scored";
const MAX_MISMATCHES: u8 = 4;

/// Float tolerance for CFD specificity. We've validated that the math
/// reproduces FlashFry exactly, so we only need to absorb floating-point
/// rounding noise from multiplication order — a single ULP is ~2.2e-16.
const FLOAT_TOL: f64 = 1e-9;

fn fixture_present() -> bool {
    Path::new(CHR22_PATH).exists()
        && Path::new(QUERY_FASTA).exists()
        && Path::new(FLASHFRY_SCORED).exists()
}

struct ScoredRow {
    ot_count: u64,
    specificity: f64,
}

fn parse_scored(
    path: &Path,
    target_col: &str,
    ot_col: &str,
    spec_col: &str,
) -> HashMap<String, ScoredRow> {
    let text = std::fs::read_to_string(path).expect("read scored TSV");
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let columns: Vec<&str> = header.split('\t').collect();
    let ti = columns
        .iter()
        .position(|c| *c == target_col)
        .unwrap_or_else(|| panic!("missing '{target_col}' column"));
    let oi = columns
        .iter()
        .position(|c| *c == ot_col)
        .unwrap_or_else(|| panic!("missing '{ot_col}' column"));
    let si = columns
        .iter()
        .position(|c| *c == spec_col)
        .unwrap_or_else(|| panic!("missing '{spec_col}' column"));
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let target = fields[ti].to_string();
            let row = ScoredRow {
                ot_count: fields[oi].parse().expect("otCount integer"),
                specificity: fields[si].parse().expect("specificity float"),
            };
            (target, row)
        })
        .collect()
}

#[test]
#[allow(clippy::too_many_lines)] // integration test reads linearly
fn cfd_and_otcount_match_flashfry_on_1000_random_chr22_guides() {
    if !fixture_present() {
        eprintln!(
            "Skipping 1000-gRNA equivalence test: fixture missing.\n\
             Setup (~10 s):\n\
             \n  1. Build the chr22 NGG database (one-time, ~30 s):\n     \
                cd /tmp/ffry-quickstart\n     \
                java -Xmx4g -jar FlashFry-assembly-1.15.jar index \\\n       \
                    --tmpLocation ./tmp --database chr22_cas9ngg_database \\\n       \
                    --reference chr22.fa.gz --enzyme spcas9ngg\n\
             \n  2. Generate the 1000-gRNA query FASTA from chr22:\n     \
                python3 -c 'import gzip,re,random;random.seed(42);\\\n       \
                s=\"\".join(l.strip() for l in gzip.open(\"chr22.fa.gz\",\"rt\") if not l.startswith(\">\"));\\\n       \
                starts=[m.start() for m in re.finditer(r\"(?=([ACGT]{{20}}[ACGT]GG))\", s)];\\\n       \
                f=open(\"random_1000.fasta\",\"w\");\\\n       \
                _=[f.write(f\">chr22_site_{{{{i}}}}\\n{{{{s[start:start+23]}}}}\\n\") for i,start in enumerate(sorted(random.sample(starts,1000)))]'\n\
             \n  3. Run FlashFry to generate the ground-truth scored output:\n     \
                java -Xmx4g -jar FlashFry-assembly-1.15.jar discover \\\n       \
                    --database chr22_cas9ngg_database \\\n       \
                    --fasta random_1000.fasta --output random_1000.output\n     \
                java -Xmx4g -jar FlashFry-assembly-1.15.jar score \\\n       \
                    --input random_1000.output --output random_1000.scored \\\n       \
                    --scoringMetrics doench2016cfd \\\n       \
                    --database chr22_cas9ngg_database"
        );
        return;
    }

    // Run our pipeline.
    let tmp = std::env::temp_dir().join(format!(
        "crispr-ots-1000-test-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos())
    ));
    std::fs::create_dir_all(&tmp).expect("temp dir");
    let our_output = tmp.join("ours.tsv");

    let t_build = std::time::Instant::now();
    let table = build_table_in_memory(&BuildConfig {
        reference: PathBuf::from(CHR22_PATH),
        enzyme: Enzyme::spcas9_ngg(),
        output: PathBuf::from("/dev/null"),
        bin_width: None,
    })
    .expect("build chr22 table");
    let build_wall = t_build.elapsed();

    let t_scan = std::time::Instant::now();
    run_discover(
        &table,
        &DiscoverConfig {
            input: DiscoverInput::QueryFasta(PathBuf::from(QUERY_FASTA)),
            max_mismatches: MAX_MISMATCHES,
            scores: vec![ScoreMetric::Cfd],
            output: our_output.clone(),
            format: OutputFormat::Tsv,
            spec_convention: SpecConvention::Flashfry,
            threshold: None,
            max_per_bin: None,
            ..Default::default()
        },
    )
    .expect("discover succeeds");
    let scan_wall = t_scan.elapsed();
    eprintln!(
        "crispr-ots build (in-memory): {:.2}s; discover --score cfd: {:.2}s on 1000 chr22 NGG guides",
        build_wall.as_secs_f64(),
        scan_wall.as_secs_f64()
    );

    // Compare.
    let ours = parse_scored(&our_output, "target", "otCount", "cfd_specificity");
    let ff = parse_scored(
        Path::new(FLASHFRY_SCORED),
        "target",
        "otCount",
        "DoenchCFD_specificityscore",
    );

    eprintln!(
        "Loaded: ours={} rows, FlashFry={} rows (FlashFry drops OVERFLOW guides during scoring)",
        ours.len(),
        ff.len()
    );

    let common: Vec<&String> = ff.keys().filter(|t| ours.contains_key(*t)).collect();
    assert_eq!(
        common.len(),
        ff.len(),
        "every guide FlashFry scored should also be in our output"
    );

    let mut ot_mismatches = 0_usize;
    let mut worst_spec_drift = 0.0_f64;
    let mut spec_mismatches = 0_usize;
    for target in &common {
        let f = &ff[*target];
        let o = &ours[*target];
        if f.ot_count != o.ot_count {
            ot_mismatches += 1;
            eprintln!(
                "  ot_count drift for {target}: FlashFry={}, ours={}",
                f.ot_count, o.ot_count
            );
        }
        let diff = (f.specificity - o.specificity).abs();
        if diff > worst_spec_drift {
            worst_spec_drift = diff;
        }
        if diff > FLOAT_TOL {
            spec_mismatches += 1;
        }
    }

    eprintln!(
        "Comparison: otCount mismatches = {}/{}, spec mismatches > {FLOAT_TOL:e} = {}/{}, max |Δ spec| = {:e}",
        ot_mismatches,
        common.len(),
        spec_mismatches,
        common.len(),
        worst_spec_drift
    );

    assert_eq!(
        ot_mismatches, 0,
        "otCount must match FlashFry on every scoreable guide"
    );
    assert_eq!(
        spec_mismatches, 0,
        "CFD specificity must match FlashFry within {FLOAT_TOL:e}",
    );

    let _ = std::fs::remove_dir_all(&tmp);
}
