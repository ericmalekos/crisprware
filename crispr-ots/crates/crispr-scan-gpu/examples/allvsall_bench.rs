//! Genome-wide all-vs-all GPU scaling probe.
//!
//! Uses the index's own sites as guides (the real all-vs-all input) and times
//! `GpuScanner::scan_counts` at increasing guide counts, projecting the
//! per-guide rate to the full site count.
//!
//! Usage: `allvsall_bench <index_prefix|.crot> [max_mm=4] [n1,n2,...]`

#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::doc_markdown
)]

use std::path::Path;
use std::time::Instant;

use crispr_db::MmapDb;
use crispr_scan::Guide;
use crispr_scan_gpu::GpuScanner;

fn main() {
    let mut args = std::env::args().skip(1);
    let idx = args
        .next()
        .expect("usage: allvsall_bench <index> [max_mm] [n1,n2,..]");
    let max_mm: u8 = args.next().map_or(4, |s| s.parse().expect("max_mm integer"));
    let ns: Vec<usize> = args.next().map_or_else(
        || vec![10_000, 100_000, 1_000_000],
        |s| {
            s.split(',')
                .map(|x| x.parse().expect("guide count integer"))
                .collect()
        },
    );
    let path = if Path::new(&idx)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("crot"))
    {
        idx.clone()
    } else {
        format!("{idx}.crot")
    };

    let db = MmapDb::open(Path::new(&path)).expect("open index");
    let total_sites = db.entries().len();
    eprintln!("index: {total_sites} entries; max_mm={max_mm}");

    let t = Instant::now();
    let scanner = GpuScanner::new(&db).expect("build GPU scanner");
    eprintln!(
        "GpuScanner::new (flatten + upload): {:.2}s",
        t.elapsed().as_secs_f64()
    );

    let entries = db.entries();
    let mut validated = false;
    for &n in &ns {
        let n = n.min(total_sites);
        // Sample guides uniformly across the whole index (by stride), NOT the
        // first n: entries are sorted by bin key, so the first n are all
        // low-complexity poly-A-prefix sites with pathological off-target
        // counts. A strided sample is representative of the real all-vs-all mix.
        let step = (total_sites / n).max(1);
        let guides: Vec<Guide> = entries
            .iter()
            .step_by(step)
            .take(n)
            .map(|e| Guide {
                id: String::new(),
                site: e.site(),
            })
            .collect();
        let t = Instant::now();
        let counts = scanner.scan_counts_prefilter(&guides, max_mm);
        let secs = t.elapsed().as_secs_f64();

        // One-time correctness gate: the prefilter must produce identical counts
        // to the brute-force kernel (which visits every entry).
        if !validated {
            let brute = scanner.scan_counts(&guides, max_mm);
            assert!(
                counts == brute,
                "prefilter counts differ from brute force at n={n}"
            );
            eprintln!("validation OK: prefilter == brute-force counts ({n} guides)");
            validated = true;
        }

        let total_ot: u64 = counts.iter().map(|&c| u64::from(c)).sum();
        let per_guide_us = secs * 1e6 / n as f64;
        let proj_full_s = per_guide_us * total_sites as f64 / 1e6;
        eprintln!(
            "n={n:>9}  prefilter={secs:8.3}s  hits_total={total_ot:>13}  \
             per-guide={per_guide_us:7.2}us  -> all-vs-all({total_sites}) ~ {proj_full_s:8.0}s ({:.2} h)",
            proj_full_s / 3600.0
        );
    }
}
