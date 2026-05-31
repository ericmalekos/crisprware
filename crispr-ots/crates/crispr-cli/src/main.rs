//! `crispr-ots` command-line entry point.
//!
//! Thin wrapper over the `crispr_cli` library. Real work lives there so
//! integration tests can call it in-process.
//!
//! The CLI mirrors GuideScan2's `index` / `enumerate` subcommands by
//! design — that lets `crisprware`'s existing `subprocess.run(["guidescan",
//! ...])` calls flip to `crispr-ots` with a single-word swap. The native
//! `build` / `discover` names are kept as aliases for users who already
//! follow the existing CLI.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};

use crispr_cli::discover::enzyme_from_name;
use crispr_cli::{
    run_build, run_discover_with, BuildConfig, DiscoverConfig, DiscoverInput, OutputFormat,
    ScannerKind, ScoreMetric, SpecConvention,
};
use crispr_db::MmapDb;

#[derive(Parser, Debug)]
#[command(
    name = "crispr-ots",
    version,
    about = "Fast CRISPR off-target scanning"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build a persistent off-target database from a reference FASTA.
    /// Native alias: `build`.
    ///
    /// Mirrors `guidescan index --index <prefix> <fasta>`: writes
    /// `<prefix>.crot` (the actual database) plus zero-byte stub files at
    /// `<prefix>.forward`, `<prefix>.reverse`, and `<prefix>.gs` so that
    /// downstream tooling (e.g. crisprware's `check_files_exist`) keeps
    /// working without modification.
    #[command(alias = "build")]
    Index(IndexArgs),

    /// Enumerate off-targets for a set of guides against a previously-built
    /// database. Native alias: `discover`.
    ///
    /// Mirrors `guidescan enumerate ... <index>`: takes the index prefix as
    /// a positional argument and accepts the full GuideScan2 flag surface
    /// (`--kmers-file`, `--mismatches`, `--threads`, `--format`, etc.).
    /// Additional native flags (`--queries`, `--score`, `--spec-convention`)
    /// are available for users who don't need GS2 compatibility.
    #[command(alias = "discover")]
    Enumerate(EnumerateArgs),
}

#[derive(Parser, Debug)]
struct IndexArgs {
    /// Reference genome FASTA (gzip auto-detected by `.gz` extension).
    /// Positional argument; matches GuideScan2's `guidescan index ...
    /// <fasta>` call style.
    #[arg(value_name = "FASTA")]
    fasta: PathBuf,

    /// Output index prefix. The database is written to `<prefix>.crot`
    /// and zero-byte stubs are emitted at `<prefix>.forward`,
    /// `<prefix>.reverse`, and `<prefix>.gs` for downstream
    /// compatibility. If `<prefix>` already ends in `.crot`, no second
    /// extension is appended.
    #[arg(long, value_name = "PREFIX", alias = "output", short = 'o')]
    index: PathBuf,

    /// Enzyme preset. Recognised: spcas9ngg, spcas9nag, spcas9 (NGG+NAG),
    /// cpf1 (a.k.a. cas12a).
    #[arg(short, long, default_value = "spcas9ngg")]
    enzyme: String,

    /// Bin-prefix width in bases (1-15). Higher widths sharpen the
    /// per-bin prefilter at the cost of a larger bin-offset table.
    /// Empirical sweet spots: 7 (default) for chr22-scale genomes,
    /// 9-10 for the full human genome, 11 for full mouse. See
    /// `BENCHMARKS.md`.
    #[arg(short = 'w', long, value_parser = parse_bin_width)]
    bin_width: Option<u8>,
}

#[derive(Parser, Debug)]
struct EnumerateArgs {
    /// Index prefix produced by `crispr-ots index --index <prefix>`.
    /// The database is loaded from `<prefix>.crot`; if `<prefix>` ends
    /// in `.crot` it is used directly.
    #[arg(value_name = "INDEX")]
    index: PathBuf,

    /// Per-guide CSV input in the crisprware / GuideScan2 format: a
    /// single-column TSV whose header is the literal string
    /// `id,sequence,pam,chromosome,position,sense` and whose rows are
    /// comma-separated values of that shape.
    #[arg(long, value_name = "PATH", conflicts_with = "queries")]
    kmers_file: Option<PathBuf>,

    /// Query FASTA (native form). Discovers protospacer sites inside
    /// the FASTA via the same scanner used to build the index, then
    /// enumerates off-targets for each. Conflicts with `--kmers-file`.
    #[arg(short, long, value_name = "PATH")]
    queries: Option<PathBuf>,

    /// Maximum mismatches in reported off-targets. Alias:
    /// `--max-mismatches` (the native CLI's name).
    #[arg(short = 'm', long, alias = "max-mismatches", default_value_t = 4)]
    mismatches: u8,

    /// Number of rayon worker threads. 0 = use rayon's default (all
    /// logical CPUs).
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Scan backend. `cpu` (default) is the multithreaded SIMD bin-scanner
    /// (AVX-512 VPOPCNTDQ where the CPU supports it, otherwise AVX2). `gpu`
    /// runs the CUDA kernel and requires a binary built with `--features gpu`
    /// plus an available NVIDIA device.
    #[arg(long, value_enum, default_value_t = ScannerArg::Cpu)]
    scanner: ScannerArg,

    /// Output format. `csv` emits one row per off-target hit, matching
    /// GuideScan2's enumerate CSV. `tsv` emits one row per guide with
    /// every off-target listed in an `otSequences` column (FlashFry
    /// style). `both` emits the CSV at `--output` and the TSV at
    /// `<output>.detail.tsv`.
    #[arg(short, long, value_enum, default_value_t = OutputFormatArg::Csv)]
    format: OutputFormatArg,

    /// Specificity-formula convention. `guidescan` uses
    /// `1 / Σ_all_rows(cfd)` (the GuideScan2 paper). `flashfry` uses
    /// `1 / (1 + Σ_offtarget cfd × count)` (the FlashFry paper).
    /// Defaults: csv→guidescan, tsv→flashfry, both→guidescan.
    #[arg(long, value_enum)]
    spec_convention: Option<SpecConventionArg>,

    /// Score metrics to apply. Currently supports `cfd`. May be repeated.
    /// Defaults to `cfd` (a specificity column is the whole point of
    /// the GS2-compat output).
    #[arg(long, value_parser = parse_score_metric)]
    score: Vec<ScoreMetric>,

    /// Output path. For `--format both` this is the CSV; the TSV is
    /// written to `<output>.detail.tsv`.
    #[arg(short, long)]
    output: PathBuf,

    /// Per-(guide, mismatch-count) off-target cap. Once a bin reaches
    /// N distinct off-target sequences, additional sequences at that
    /// mismatch level for the same guide are dropped. Affects
    /// `otCount`, the `otSequences` listing, and the CFD specificity
    /// denominator (kept sequences only). Tames low-complexity guides
    /// with 100 k+ off-targets at high mismatch counts. Set to -1 to
    /// disable.
    #[arg(long, default_value_t = 500, allow_hyphen_values = true)]
    max_off_targets_per_bin: i64,

    // ── GuideScan2 compatibility surface ──────────────────────────────
    /// Cap on per-guide off-target rows (-1 = uncapped). Accepted for
    /// GuideScan2 CLI compatibility; today we always emit every
    /// off-target up to `--mismatches`. See also
    /// `--max-off-targets-per-bin` for our finer-grained per-mismatch
    /// cap.
    #[arg(long, default_value_t = -1, allow_hyphen_values = true)]
    max_off_targets: i64,

    /// Discard guides whose nearest *off*-target (excluding the on-target
    /// site) has mismatch ≤ this value. -1 disables filtering. Matches
    /// `guidescan enumerate --threshold`.
    #[arg(long, default_value_t = -1, allow_hyphen_values = true)]
    threshold: i32,

    /// GuideScan2's `--mode succinct|complete`. Accepted as a no-op: we
    /// always enumerate every off-target up to `--mismatches`, which
    /// is what `complete` does.
    #[arg(long, default_value = "complete")]
    mode: String,

    /// Alternate PAM sequence(s). May be repeated. The literal string
    /// `None` is accepted as a no-op for crisprware compatibility.
    #[arg(long, action = clap::ArgAction::Append)]
    alt_pam: Vec<String>,

    /// RNA bulges. Accepted only when 0; non-zero values are rejected
    /// at runtime (bulge support is out of scope for the bin-scan
    /// engine).
    #[arg(long, default_value_t = 0)]
    rna_bulges: u32,

    /// DNA bulges. Same constraint as `--rna-bulges`.
    #[arg(long, default_value_t = 0)]
    dna_bulges: u32,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormatArg {
    Csv,
    Tsv,
    Both,
}

impl From<OutputFormatArg> for OutputFormat {
    fn from(arg: OutputFormatArg) -> Self {
        match arg {
            OutputFormatArg::Csv => Self::Csv,
            OutputFormatArg::Tsv => Self::Tsv,
            OutputFormatArg::Both => Self::Both,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum SpecConventionArg {
    Guidescan,
    Flashfry,
}

impl From<SpecConventionArg> for SpecConvention {
    fn from(arg: SpecConventionArg) -> Self {
        match arg {
            SpecConventionArg::Guidescan => Self::Guidescan,
            SpecConventionArg::Flashfry => Self::Flashfry,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ScannerArg {
    Cpu,
    Gpu,
}

impl From<ScannerArg> for ScannerKind {
    fn from(arg: ScannerArg) -> Self {
        match arg {
            ScannerArg::Cpu => Self::Cpu,
            ScannerArg::Gpu => Self::Gpu,
        }
    }
}

fn parse_score_metric(s: &str) -> Result<ScoreMetric, String> {
    match s.to_ascii_lowercase().as_str() {
        "cfd" | "doench2016cfd" => Ok(ScoreMetric::Cfd),
        "cfd-cas12a" | "cfd-cas12a:2xnls" | "cas12a-2xnls" => Ok(ScoreMetric::Cas12aTwoXNls),
        "cfd-cas12a:encas12a" | "cas12a-encas12a" | "encas12a" => {
            Ok(ScoreMetric::Cas12aEnCas12a)
        }
        other => Err(format!(
            "unknown score metric '{other}' (supported: cfd, cfd-cas12a, \
             cfd-cas12a:encas12a)"
        )),
    }
}

fn parse_bin_width(s: &str) -> Result<u8, String> {
    let w: u8 = s
        .parse()
        .map_err(|_| format!("bin-width must be a small integer, got '{s}'"))?;
    if !(1..=15).contains(&w) {
        return Err(format!(
            "bin-width {w} out of range; must be in 1..=15 (4^15 = ~1 G bins is the u32 key limit)"
        ));
    }
    Ok(w)
}

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn real_main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Index(args) => run_index_cmd(args),
        Command::Enumerate(args) => run_enumerate_cmd(args),
    }
}

fn run_index_cmd(args: IndexArgs) -> Result<()> {
    let enzyme = enzyme_from_name(&args.enzyme)
        .ok_or_else(|| anyhow!("unknown enzyme '{}'", args.enzyme))?;
    let crot_path = resolve_index_crot_path(&args.index);
    let config = BuildConfig {
        reference: args.fasta,
        enzyme,
        output: crot_path.clone(),
        bin_width: args.bin_width,
    };
    run_build(&config).context("index build failed")?;
    emit_stub_index_files(&args.index).context("stub-file emission failed")?;
    Ok(())
}

fn run_enumerate_cmd(args: EnumerateArgs) -> Result<()> {
    if args.rna_bulges != 0 || args.dna_bulges != 0 {
        bail!(
            "non-zero --rna-bulges/--dna-bulges not supported by crispr-ots (out of scope \
             for the bin-scan engine); pass 0 to silence this error"
        );
    }

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .context("could not initialise the rayon thread pool")?;
    }

    let input = match (args.kmers_file, args.queries) {
        (Some(_), Some(_)) => unreachable!("clap conflicts_with prevents this"),
        (Some(kmers), None) => DiscoverInput::KmersCsv(kmers),
        (None, Some(queries)) => DiscoverInput::QueryFasta(queries),
        (None, None) => bail!(
            "one of --kmers-file or --queries is required to enumerate off-targets"
        ),
    };

    let format: OutputFormat = args.format.into();
    let spec_convention = args
        .spec_convention
        .map_or_else(|| spec_default_for(format), Into::into);

    let scores = if args.score.is_empty() {
        vec![ScoreMetric::Cfd]
    } else {
        args.score
    };

    let crot_path = resolve_index_crot_path(&args.index);
    let mmap = MmapDb::open(&crot_path)
        .with_context(|| format!("could not mmap database at {}", crot_path.display()))?;

    let config = DiscoverConfig {
        input,
        max_mismatches: args.mismatches,
        scores,
        output: args.output,
        format,
        spec_convention,
        threshold: if args.threshold < 0 {
            None
        } else {
            Some(u8::try_from(args.threshold).expect("threshold fits in u8"))
        },
        max_per_bin: if args.max_off_targets_per_bin < 0 {
            None
        } else {
            Some(
                u32::try_from(args.max_off_targets_per_bin)
                    .expect("--max-off-targets-per-bin fits in u32"),
            )
        },
    };
    run_discover_with(&mmap, &config, args.scanner.into()).context("enumerate failed")
}

fn spec_default_for(format: OutputFormat) -> SpecConvention {
    match format {
        OutputFormat::Tsv => SpecConvention::Flashfry,
        OutputFormat::Csv | OutputFormat::Both => SpecConvention::Guidescan,
    }
}

/// Resolve the on-disk path for the `.crot` database from a user-supplied
/// index prefix. If the prefix already ends in `.crot` it is used as-is;
/// otherwise the suffix is appended.
fn resolve_index_crot_path(prefix: &Path) -> PathBuf {
    if prefix
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("crot"))
    {
        prefix.to_path_buf()
    } else {
        let mut s = prefix.as_os_str().to_owned();
        s.push(".crot");
        PathBuf::from(s)
    }
}

/// Emit three zero-byte files at `<prefix>.forward`, `<prefix>.reverse`,
/// and `<prefix>.gs` so that downstream tooling (e.g. crisprware's
/// `check_files_exist`) can probe for them without modification.
fn emit_stub_index_files(prefix: &Path) -> Result<()> {
    for suffix in ["forward", "reverse", "gs"] {
        let mut p = prefix.as_os_str().to_owned();
        p.push(".");
        p.push(suffix);
        let path = PathBuf::from(p);
        std::fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(&path)
            .with_context(|| format!("could not create stub file {}", path.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crot_path_appends_when_missing() {
        let p = resolve_index_crot_path(Path::new("/tmp/mm39"));
        assert_eq!(p, PathBuf::from("/tmp/mm39.crot"));
    }

    #[test]
    fn crot_path_left_alone_when_already_present() {
        let p = resolve_index_crot_path(Path::new("/tmp/mm39.crot"));
        assert_eq!(p, PathBuf::from("/tmp/mm39.crot"));
    }

    #[test]
    fn crot_path_handles_dotted_prefix() {
        // crisprware does `output_directory/foo_gscan2/foo_gscan2` —
        // dotted directory components shouldn't break the suffix check.
        let p = resolve_index_crot_path(Path::new("/out/foo.gscan2/foo_gscan2"));
        assert_eq!(p, PathBuf::from("/out/foo.gscan2/foo_gscan2.crot"));
    }

    #[test]
    fn cli_parses_guidescan_enumerate_arg_shape() {
        // The exact invocation crisprware emits today (see
        // score_guides.py:186-211). Validates that we accept it without
        // surfacing any `--rna-bulges`/`--mode`/`--alt-pam` flag errors.
        let parsed = Cli::try_parse_from([
            "crispr-ots",
            "enumerate",
            "--max-off-targets",
            "-1",
            "--threads",
            "8",
            "--mismatches",
            "3",
            "--format",
            "csv",
            "--rna-bulges",
            "0",
            "--dna-bulges",
            "0",
            "--threshold",
            "-1",
            "--mode",
            "complete",
            "--alt-pam",
            "None",
            "--kmers-file",
            "/tmp/q.csv",
            "--output",
            "/tmp/o.csv",
            "/tmp/db",
        ])
        .expect("crisprware-style enumerate args should parse");
        match parsed.command {
            Command::Enumerate(args) => {
                assert_eq!(args.mismatches, 3);
                assert_eq!(args.threads, 8);
                assert_eq!(args.threshold, -1);
                assert!(matches!(args.format, OutputFormatArg::Csv));
                assert_eq!(args.kmers_file.as_deref(), Some(Path::new("/tmp/q.csv")));
                assert_eq!(args.output, PathBuf::from("/tmp/o.csv"));
                assert_eq!(args.index, PathBuf::from("/tmp/db"));
            }
            Command::Index(_) => panic!("expected Enumerate command, got Index"),
        }
    }

    #[test]
    fn cli_parses_guidescan_index_arg_shape() {
        let parsed = Cli::try_parse_from([
            "crispr-ots",
            "index",
            "--index",
            "/tmp/db",
            "/tmp/ref.fa.gz",
        ])
        .expect("crisprware-style index args should parse");
        match parsed.command {
            Command::Index(args) => {
                assert_eq!(args.fasta, PathBuf::from("/tmp/ref.fa.gz"));
                assert_eq!(args.index, PathBuf::from("/tmp/db"));
                assert_eq!(args.enzyme, "spcas9ngg");
            }
            Command::Enumerate(_) => panic!("expected Index command, got Enumerate"),
        }
    }

    #[test]
    fn cli_accepts_build_and_discover_aliases() {
        let parsed = Cli::try_parse_from([
            "crispr-ots",
            "build",
            "--index",
            "/tmp/db",
            "/tmp/ref.fa",
        ])
        .expect("`build` alias should resolve to Index");
        assert!(matches!(parsed.command, Command::Index(_)));

        let parsed = Cli::try_parse_from([
            "crispr-ots",
            "discover",
            "--queries",
            "/tmp/q.fa",
            "--output",
            "/tmp/o.tsv",
            "/tmp/db",
        ])
        .expect("`discover` alias should resolve to Enumerate");
        assert!(matches!(parsed.command, Command::Enumerate(_)));
    }
}
