//! Bin-pair-tiled all-vs-all self-scan: validate + benchmark.
//!
//! Validates `self_scan_counts` against the per-guide `scan_counts_prefilter`
//! on a small bin range, then times a representative middle slice and projects
//! to the full genome self-scan. Pass a 4th arg `full` to also run the whole
//! self-scan.
//!
//! Usage: `selfscan_bench <index> [max_mm=4] [slice_bins=400000] [full]`

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
        .expect("usage: selfscan_bench <index> [max_mm] [slice_bins] [full]");
    let max_mm: u8 = args.next().map_or(4, |s| s.parse().expect("max_mm"));
    let slice_bins: u32 = args.next().map_or(400_000, |s| s.parse().expect("slice_bins"));
    let full = args.next().is_some_and(|s| s == "full" || s == "1");
    let path = if Path::new(&idx)
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("crot"))
    {
        idx.clone()
    } else {
        format!("{idx}.crot")
    };

    let db = MmapDb::open(Path::new(&path)).expect("open index");
    let t = Instant::now();
    let scanner = GpuScanner::new(&db).expect("build GPU scanner");
    eprintln!("GpuScanner::new: {:.2}s", t.elapsed().as_secs_f64());
    let bin_count = scanner.bin_count();
    let stride = max_mm as usize + 1;
    let bos = db.bin_offsets();
    let entries = db.entries();

    // --- correctness: tiled self-scan must match the per-guide prefilter ---
    let v_lo = bin_count / 2;
    let v_hi = (v_lo + 4000).min(bin_count);
    let self_v = scanner.self_scan_counts(v_lo, v_hi, max_mm);
    let mut guides = Vec::new();
    let mut gidx = Vec::new(); // global entry index per guide
    for b in v_lo..v_hi {
        let bo = bos[b as usize];
        for j in bo.start_index..bo.start_index + bo.count {
            guides.push(Guide {
                id: String::new(),
                site: entries[j as usize].site(),
            });
            gidx.push(j as usize);
        }
    }
    let pf = scanner.scan_counts_prefilter(&guides, max_mm);
    let ok = gidx.iter().enumerate().all(|(i, &g)| {
        (0..stride).all(|m| self_v[g * stride + m] == pf[i * stride + m])
    });
    eprintln!(
        "validation ({} guides in bins [{v_lo},{v_hi})): tiled self-scan {} per-guide prefilter",
        guides.len(),
        if ok { "==" } else { "!=" }
    );
    assert!(ok, "tiled self-scan disagrees with per-guide prefilter");

    // --- benchmark a representative middle slice ---
    let s_lo = bin_count * 3 / 8;
    let s_hi = (s_lo + slice_bins).min(bin_count);
    let slice_entries: u64 = (s_lo..s_hi).map(|b| u64::from(bos[b as usize].count)).sum();
    let t = Instant::now();
    let counts = scanner.self_scan_counts(s_lo, s_hi, max_mm);
    let secs = t.elapsed().as_secs_f64();
    let total_hits: u64 = counts.iter().map(|&c| u64::from(c)).sum();
    let proj_by_entries = secs * f64::from(scanner.entry_count()) / (slice_entries.max(1) as f64);
    eprintln!(
        "slice bins [{s_lo},{s_hi}) = {slice_entries} entries: {secs:.3}s total (incl ~5.5GB dl), \
         {total_hits} hits -> full ~ {proj_by_entries:.0}s ({:.2} h) [by entries]",
        proj_by_entries / 3600.0
    );

    if full {
        let t = Instant::now();
        let counts = scanner.self_scan_counts(0, bin_count, max_mm);
        let hits: u64 = counts.iter().map(|&c| u64::from(c)).sum();
        eprintln!(
            "FULL self-scan [0,{bin_count}): {:.1}s, {hits} total hits",
            t.elapsed().as_secs_f64()
        );
    }
}
