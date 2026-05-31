//! End-to-end CFD equivalence test against FlashFry.
//!
//! Runs the `discover --score cfd` pipeline in-process via `run_discover`
//! and diffs the produced TSV's `cfd_specificity` column against the
//! output of FlashFry's `score --scoringMetrics doench2016cfd` on the
//! same fixture. Auto-skips if either input is missing.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crispr_cli::{
    build_table_in_memory, run_discover, BuildConfig, DiscoverConfig, DiscoverInput, OutputFormat,
    ScoreMetric, SpecConvention,
};
use crispr_encoding::Enzyme;

const CHR22_PATH: &str = "/tmp/ffry-quickstart/chr22.fa.gz";
const EMX1_PATH: &str = "/tmp/ffry-quickstart/EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta";
/// Output of `FlashFry score --scoringMetrics doench2016cfd …`.
const FLASHFRY_SCORED_PATH: &str = "/tmp/ffry-quickstart/EMX1.scored";
const MAX_MISMATCHES: u8 = 4;

/// Numerical tolerance for floating-point equality. We've validated the
/// math reproduces FlashFry exactly, so 1e-12 is loose enough to absorb
/// any sane formatting wobble while still catching real drift.
const FLOAT_TOL: f64 = 1e-12;

fn fixture_present() -> bool {
    Path::new(CHR22_PATH).exists()
        && Path::new(EMX1_PATH).exists()
        && Path::new(FLASHFRY_SCORED_PATH).exists()
}

/// Map FlashFry's `target` column → its `DoenchCFD_specificityscore`.
fn parse_flashfry_scored(path: &str) -> HashMap<String, f64> {
    let text = std::fs::read_to_string(path).expect("read FlashFry scored output");
    let mut lines = text.lines();
    let header = lines.next().expect("header line").to_string();
    let columns: Vec<&str> = header.split('\t').collect();
    let target_idx = columns
        .iter()
        .position(|c| *c == "target")
        .expect("target column present");
    let spec_idx = columns
        .iter()
        .position(|c| *c == "DoenchCFD_specificityscore")
        .expect("DoenchCFD_specificityscore column present");
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let target = fields[target_idx].to_string();
            let spec: f64 = fields[spec_idx].parse().expect("spec is a float");
            (target, spec)
        })
        .collect()
}

/// Map our `target` column → our `cfd_specificity` column.
fn parse_our_tsv(path: &Path) -> HashMap<String, f64> {
    let text = std::fs::read_to_string(path).expect("read our TSV");
    let mut lines = text.lines();
    let header = lines.next().expect("header line").to_string();
    let columns: Vec<&str> = header.split('\t').collect();
    let target_idx = columns
        .iter()
        .position(|c| *c == "target")
        .expect("target column");
    let spec_idx = columns
        .iter()
        .position(|c| *c == "cfd_specificity")
        .expect("cfd_specificity column");
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let target = fields[target_idx].to_string();
            let spec: f64 = fields[spec_idx].parse().expect("our spec is a float");
            (target, spec)
        })
        .collect()
}

#[test]
fn cfd_specificity_matches_flashfry_on_emx1() {
    if !fixture_present() {
        eprintln!(
            "Skipping CFD equivalence test: fixture missing.\n\
             Run FlashFry's score subcommand to generate EMX1.scored:\n  \
             java -Xmx4g -jar FlashFry-assembly-1.15.jar score \\\n    \
                 --input EMX1.output --output EMX1.scored \\\n    \
                 --scoringMetrics doench2016cfd --database chr22_cas9ngg_database"
        );
        return;
    }

    let tmp = std::env::temp_dir().join(format!(
        "crispr-ots-cli-cfd-test-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos())
    ));
    std::fs::create_dir_all(&tmp).expect("temp dir");
    let our_output = tmp.join("emx1.crispr-ots.tsv");

    let table = build_table_in_memory(&BuildConfig {
        reference: PathBuf::from(CHR22_PATH),
        enzyme: Enzyme::spcas9_ngg(),
        output: PathBuf::from("/dev/null"), // unused — we don't call run_build
        bin_width: None,
    })
    .expect("build chr22 table in memory");

    let config = DiscoverConfig {
        input: DiscoverInput::QueryFasta(PathBuf::from(EMX1_PATH)),
        max_mismatches: MAX_MISMATCHES,
        scores: vec![ScoreMetric::Cfd],
        output: our_output.clone(),
        format: OutputFormat::Tsv,
        spec_convention: SpecConvention::Flashfry,
        threshold: None,
                max_per_bin: None,
    };
    run_discover(&table, &config).expect("discover succeeds");

    let ours = parse_our_tsv(&our_output);
    let flashfry = parse_flashfry_scored(FLASHFRY_SCORED_PATH);

    assert_eq!(
        ours.len(),
        flashfry.len(),
        "row count mismatch: ours={}, flashfry={}",
        ours.len(),
        flashfry.len()
    );

    for (target, ff_spec) in &flashfry {
        let our_spec = ours
            .get(target)
            .copied()
            .unwrap_or_else(|| panic!("missing guide in our output: {target}"));
        let diff = (our_spec - ff_spec).abs();
        assert!(
            diff < FLOAT_TOL,
            "CFD specificity drift for {target}: ours={our_spec}, flashfry={ff_spec}, diff={diff:e}",
        );
        eprintln!("  {target}: ours={our_spec}, flashfry={ff_spec} (Δ={diff:e})");
    }
    eprintln!(
        "✓ all {} CFD specificities match FlashFry within {FLOAT_TOL:e}",
        ours.len()
    );

    let _ = std::fs::remove_dir_all(&tmp);
}
