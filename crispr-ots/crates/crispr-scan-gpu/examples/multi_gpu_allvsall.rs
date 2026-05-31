//! Multi-GPU sharded all-vs-all (genome self-scan).
//!
//! Splits the index's entries into N contiguous shards (one per GPU), and on
//! each GPU builds a resident `GpuScanner` and runs the per-guide bin-prefilter
//! over its shard's guides in parallel. Contiguous shards keep each guide's
//! candidate bins L2-resident (the fast path). Reports the multi-GPU wall time
//! and projects to the full self-scan.
//!
//! Usage: `multi_gpu_allvsall <index> [n_gpus=1] [max_mm=4] [frac_percent=100]`
//! With `frac_percent < 100`, a centred slice of that size is scanned (a quick,
//! representative benchmark); 100 scans every site against the genome.

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
        .expect("usage: multi_gpu_allvsall <index> [n_gpus] [max_mm] [frac_percent]");
    let n_gpus: usize = args.next().map_or(1, |s| s.parse().expect("n_gpus"));
    let max_mm: u8 = args.next().map_or(4, |s| s.parse().expect("max_mm"));
    let frac_pct: usize = args.next().map_or(100, |s| s.parse().expect("frac_percent"));
    let path = if Path::new(&idx)
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("crot"))
    {
        idx.clone()
    } else {
        format!("{idx}.crot")
    };

    let db = MmapDb::open(Path::new(&path)).expect("open index");
    let n_entries = db.entries().len();
    // Entries are stored in bin order, so a contiguous slice is bin-adjacent
    // (the L2-friendly access pattern). For a partial benchmark take a centred
    // slice; for the real thing, all entries.
    let (e_lo, e_hi) = if frac_pct >= 100 {
        (0, n_entries)
    } else {
        let range = n_entries * frac_pct / 100;
        let lo = (n_entries - range) / 2;
        (lo, lo + range)
    };
    let range_entries = e_hi - e_lo;
    let per = range_entries / n_gpus;
    eprintln!(
        "index {n_entries} entries; all-vs-all over [{e_lo},{e_hi}) = {range_entries} guides \
         on {n_gpus} GPU(s), mm<={max_mm}"
    );

    let t = Instant::now();
    let results: Vec<(usize, usize, u64, f64)> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..n_gpus)
            .map(|g| {
                let lo = e_lo + g * per;
                let hi = if g == n_gpus - 1 { e_hi } else { e_lo + (g + 1) * per };
                let db_ref = &db;
                s.spawn(move || {
                    let scanner = GpuScanner::new_on_device(db_ref, g).expect("gpu scanner");
                    let guides: Vec<Guide> = db_ref.entries()[lo..hi]
                        .iter()
                        .map(|e| Guide {
                            id: String::new(),
                            site: e.site(),
                        })
                        .collect();
                    let t0 = Instant::now();
                    let counts = scanner.scan_counts_prefilter(&guides, max_mm);
                    let scan_s = t0.elapsed().as_secs_f64();
                    let total: u64 = counts.iter().map(|&c| u64::from(c)).sum();
                    (g, hi - lo, total, scan_s)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let wall = t.elapsed().as_secs_f64();

    let grand: u64 = results.iter().map(|r| r.2).sum();
    let max_scan = results.iter().map(|r| r.3).fold(0.0f64, f64::max);
    for &(g, ng, total, scan_s) in &results {
        eprintln!("  gpu{g}: {ng} guides -> {total} hits, scan {scan_s:.1}s");
    }
    let proj = wall * n_entries as f64 / range_entries as f64;
    eprintln!(
        "multi-GPU: wall {wall:.1}s (slowest per-GPU scan {max_scan:.1}s); grand total {grand} hits"
    );
    eprintln!(
        "  -> full all-vs-all ({n_entries} guides) ~ {proj:.0}s ({:.1} min)",
        proj / 60.0
    );
}
