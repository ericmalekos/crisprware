//! FlashFry-style bin-scan engine.
//!
//! Per-bin inner loop with a bin-prefix prefilter (FlashFry's
//! `OrderedBinTraversalFactory` idea — optimization #6 in the plan): a
//! guide whose own 7-mer prefix differs from a bin's key by more than
//! the mismatch budget cannot match any site in that bin (all sites in
//! the bin share that 7-mer by construction), so the entire bin is
//! skipped for that guide. Drops the comparison count by ~10–20× on
//! genome-scale workloads.
//!
//! Phase 5: the outer bin loop runs in parallel via rayon. Each thread
//! keeps its own scratch buffers and a local `Vec<Hit>`; results are
//! concatenated at the end. Bins are independent after the prefilter
//! (no shared mutable state in the inner loop), so the parallel speedup
//! is close to the core count on machines where the workload isn't
//! memory-bandwidth-bound.
//!
//! Phase 6: the per-bin inner loop now batches four
//! `(active-guide × site)` mismatch comparisons per iteration via the
//! `wide` crate's `u64x4` (AVX2 on x86-64, NEON on aarch64, scalar
//! everywhere else). Active guides are gathered into struct-of-arrays
//! scratch buffers (`active_guide_lows`, `active_guide_highs`) at the
//! top of each bin so the SIMD chunk loads are contiguous; the broadcast
//! goes the other direction (one site against four guides). The 24-base
//! pair-popcount trick from `crispr_encoding::mismatches_masked` is
//! re-implemented lane-wise via SWAR (Hacker's Delight § 5-1).

use crispr_db::{BinSource, PackedEntry};
use crispr_encoding::{mismatches_masked, Site, SiteMask};
use rayon::prelude::*;
use wide::u64x4;

/// A single guide submitted for off-target search.
#[derive(Debug, Clone)]
pub struct Guide {
    /// User-provided identifier (FASTA name, CSV id, etc.).
    pub id: String,
    /// Encoded protospacer + on-target PAM, using the same `Enzyme` as the
    /// `BinTable` being scanned.
    pub site: Site,
}

/// One discovered off-target hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Index into the input `&[Guide]` slice.
    pub guide_index: u32,
    /// The off-target site, in canonical (5'-protospacer-first) orientation.
    pub off_target: Site,
    /// Genomic location of the off-target.
    pub position: crispr_db::Position,
    /// Hamming distance between guide and off-target, PAM bits excluded.
    pub mismatches: u8,
}

/// Engine that searches a database for off-targets.
pub trait Scanner {
    /// Find every off-target with at most `max_mismatches` mismatches.
    /// Results are returned in **engine-defined order** — callers that need
    /// a stable order should sort afterwards.
    fn scan(&self, guides: &[Guide], max_mismatches: u8) -> Vec<Hit>;
}

/// Bin-scan scanner. Works against any [`BinSource`] — currently
/// implemented for the in-memory `BinTable` (used by integration tests
/// and the build path) and the mmap'd `MmapDb` (used by the CLI's
/// discover path). The scanner pulls all of the bin metadata it needs
/// at the top of each `scan()` invocation, so the trait's virtual
/// dispatch costs one call per scan — no hot-loop indirection.
pub struct BinScanner<'a> {
    source: &'a (dyn BinSource + 'a),
}

impl std::fmt::Debug for BinScanner<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BinScanner")
            .field("enzyme", &self.source.enzyme().name)
            .field("bin_width", &self.source.bin_width())
            .finish()
    }
}

impl<'a> BinScanner<'a> {
    #[must_use]
    pub fn new(source: &'a (dyn BinSource + 'a)) -> Self {
        Self { source }
    }
}

/// Per-worker accumulator threaded through the parallel fold: reusable
/// scratch buffers (guide indices + SoA gather of guide.low / guide.high),
/// local hit collector, and (when tracing) per-thread prefilter counters.
type WorkerState = (
    Vec<u32>, // active guide indices
    Vec<u64>, // active guide.low values, indexed parallel to indices
    Vec<u64>, // active guide.high values
    Vec<Hit>, // per-thread hit accumulator
    u64,      // total (guide, bin) pairs visited
    u64,      // pairs that passed prefilter
);

/// Bits at the "high" position of every 2-bit pair within a `u64`. Same
/// constant as the one in `crispr_encoding::mismatch`; duplicated here so
/// we don't need a public re-export on the hot path.
const UPPER_PAIR_BITS_U64: u64 = 0xAAAA_AAAA_AAAA_AAAA;

/// Lane-wise count of "differing 2-bit pairs": for each 2-bit pair within
/// each `u64` lane, set one indicator bit if either bit of the pair is
/// set, then popcount across the lane. SWAR popcount (Hacker's Delight
/// § 5-1) — no hardware lane-wise popcount needed, so works on any
/// platform `wide` supports.
#[inline]
fn pop_pairs_x4(x: u64x4) -> u64x4 {
    let upper = u64x4::splat(UPPER_PAIR_BITS_U64);
    let indicator = (x & upper) | ((x << 1) & upper);
    swar_popcount_x4(indicator)
}

#[inline]
fn swar_popcount_x4(mut v: u64x4) -> u64x4 {
    let m1 = u64x4::splat(0x5555_5555_5555_5555);
    let m2 = u64x4::splat(0x3333_3333_3333_3333);
    let m4 = u64x4::splat(0x0F0F_0F0F_0F0F_0F0F);
    let h01 = u64x4::splat(0x0101_0101_0101_0101);
    v -= (v >> 1) & m1;
    v = (v & m2) + ((v >> 2) & m2);
    v = (v + (v >> 4)) & m4;
    (v * h01) >> 56
}

/// SIMD mismatch primitive: four `(guide, site)` pairs against the same
/// compare-mask in one call. `guide_lows`/`guide_highs` hold the
/// `Site.low`/`Site.high` fields of four active guides; the single
/// `site_low`/`site_high` are broadcast to all four lanes.
///
/// Mirrors `crispr_encoding::mismatches_masked` exactly, lane-wise.
#[inline]
fn mismatches_x4(
    guide_lows: u64x4,
    guide_highs: u64x4,
    site_low: u64,
    site_high: u64,
    mask_low: u64,
    mask_high: u64,
) -> u64x4 {
    let cl = (guide_lows ^ u64x4::splat(site_low)) & u64x4::splat(mask_low);
    let ch = (guide_highs ^ u64x4::splat(site_high)) & u64x4::splat(mask_high);
    pop_pairs_x4(cl) + pop_pairs_x4(ch)
}

/// Portable per-bin kernel: compare every site in `sites` against the
/// gathered active guides four-at-a-time via `wide::u64x4` (AVX2 on x86-64,
/// NEON on aarch64, scalar elsewhere) with a SWAR pair-popcount, plus a
/// scalar tail for the 0–3 leftover guides. Pushes one [`Hit`] per
/// `(guide, site)` pair within `max_mismatches`. This is the default path on
/// any CPU without AVX-512 VPOPCNTDQ.
fn scan_bin_x4(
    active_lows: &[u64],
    active_highs: &[u64],
    active_idx: &[u32],
    sites: &[PackedEntry],
    compare_mask: SiteMask,
    max_mismatches: u8,
    out: &mut Vec<Hit>,
) {
    let mask_low = compare_mask.low;
    let mask_high = compare_mask.high;
    let max_mm_u32 = u32::from(max_mismatches);
    let max_mm_u64 = u64::from(max_mm_u32);
    let n = active_idx.len();
    let n_chunks = n / 4;
    let tail_start = n_chunks * 4;
    for entry in sites {
        let site = entry.site();
        let position = entry.position();
        let site_low = site.low;
        let site_high = site.high;
        // ---- SIMD body: 4 guides per iteration. ----
        for c in 0..n_chunks {
            let base = c * 4;
            let g_lows = u64x4::new([
                active_lows[base],
                active_lows[base + 1],
                active_lows[base + 2],
                active_lows[base + 3],
            ]);
            let g_highs = u64x4::new([
                active_highs[base],
                active_highs[base + 1],
                active_highs[base + 2],
                active_highs[base + 3],
            ]);
            let mm = mismatches_x4(g_lows, g_highs, site_low, site_high, mask_low, mask_high);
            let mm_arr = mm.to_array();
            for lane in 0..4 {
                if mm_arr[lane] <= max_mm_u64 {
                    out.push(Hit {
                        guide_index: active_idx[base + lane],
                        off_target: site,
                        position,
                        mismatches: u8::try_from(mm_arr[lane]).expect("mismatch fits in u8"),
                    });
                }
            }
        }
        // ---- Scalar tail (0–3 leftover guides). ----
        for k in tail_start..n {
            let g_site = Site {
                low: active_lows[k],
                high: active_highs[k],
            };
            let mm = mismatches_masked(g_site, site, compare_mask);
            if mm <= max_mm_u32 {
                out.push(Hit {
                    guide_index: active_idx[k],
                    off_target: site,
                    position,
                    mismatches: u8::try_from(mm).expect("mismatch fits in u8"),
                });
            }
        }
    }
}

/// AVX-512 VPOPCNTDQ per-bin kernel: same arithmetic as [`scan_bin_x4`] but
/// eight guides per iteration via `__m512i`, replacing the 10-op SWAR
/// pair-popcount with a single hardware `_mm512_popcnt_epi64` per 64-bit
/// lane. Selected only when the running CPU reports `avx512f` +
/// `avx512vpopcntdq` (Genoa / Genoa-X, Ice Lake-SP, Sapphire Rapids, …).
///
/// # Safety
/// The caller must ensure the CPU supports `avx512f` and `avx512vpopcntdq`
/// (the `enable`d target features); calling on any other CPU is UB.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512vpopcntdq")]
#[allow(clippy::cast_possible_wrap, clippy::similar_names)]
unsafe fn scan_bin_avx512(
    active_lows: &[u64],
    active_highs: &[u64],
    active_idx: &[u32],
    sites: &[PackedEntry],
    compare_mask: SiteMask,
    max_mismatches: u8,
    out: &mut Vec<Hit>,
) {
    use std::arch::x86_64::{
        _mm512_add_epi64, _mm512_and_si512, _mm512_loadu_si512, _mm512_or_si512,
        _mm512_popcnt_epi64, _mm512_set1_epi64, _mm512_slli_epi64, _mm512_storeu_si512,
        _mm512_xor_si512,
    };
    let max_mm_u32 = u32::from(max_mismatches);
    let max_mm_u64 = u64::from(max_mm_u32);
    let n = active_idx.len();
    let n_chunks = n / 8;
    let tail_start = n_chunks * 8;
    // Constants broadcast to all eight 64-bit lanes once per bin. These pure
    // compute intrinsics need no `unsafe` block: inside a `#[target_feature]`
    // fn the enabled features are statically guaranteed, so only the pointer
    // loads/stores below (which dereference raw pointers) require `unsafe`.
    let upper = _mm512_set1_epi64(UPPER_PAIR_BITS_U64 as i64);
    let mask_low_v = _mm512_set1_epi64(compare_mask.low as i64);
    let mask_high_v = _mm512_set1_epi64(compare_mask.high as i64);
    for entry in sites {
        let site = entry.site();
        let position = entry.position();
        let site_low_v = _mm512_set1_epi64(site.low as i64);
        let site_high_v = _mm512_set1_epi64(site.high as i64);
        for c in 0..n_chunks {
            let base = c * 8;
            // SAFETY: `base + 8 <= n`, so the unaligned 64-byte loads stay in
            // bounds; all intrinsics require avx512f/avx512vpopcntdq, which the
            // caller guarantees.
            let mm_arr: [u64; 8] = unsafe {
                let g_lows = _mm512_loadu_si512(active_lows.as_ptr().add(base).cast());
                let g_highs = _mm512_loadu_si512(active_highs.as_ptr().add(base).cast());
                let cl = _mm512_and_si512(_mm512_xor_si512(g_lows, site_low_v), mask_low_v);
                let ch = _mm512_and_si512(_mm512_xor_si512(g_highs, site_high_v), mask_high_v);
                let ind_l = _mm512_or_si512(
                    _mm512_and_si512(cl, upper),
                    _mm512_and_si512(_mm512_slli_epi64::<1>(cl), upper),
                );
                let ind_h = _mm512_or_si512(
                    _mm512_and_si512(ch, upper),
                    _mm512_and_si512(_mm512_slli_epi64::<1>(ch), upper),
                );
                let mm = _mm512_add_epi64(_mm512_popcnt_epi64(ind_l), _mm512_popcnt_epi64(ind_h));
                let mut arr = [0u64; 8];
                _mm512_storeu_si512(arr.as_mut_ptr().cast(), mm);
                arr
            };
            for lane in 0..8 {
                if mm_arr[lane] <= max_mm_u64 {
                    out.push(Hit {
                        guide_index: active_idx[base + lane],
                        off_target: site,
                        position,
                        mismatches: u8::try_from(mm_arr[lane]).expect("mismatch fits in u8"),
                    });
                }
            }
        }
        // ---- Scalar tail (0–7 leftover guides). ----
        for k in tail_start..n {
            let g_site = Site {
                low: active_lows[k],
                high: active_highs[k],
            };
            let mm = mismatches_masked(g_site, site, compare_mask);
            if mm <= max_mm_u32 {
                out.push(Hit {
                    guide_index: active_idx[k],
                    off_target: site,
                    position,
                    mismatches: u8::try_from(mm).expect("mismatch fits in u8"),
                });
            }
        }
    }
}

impl Scanner for BinScanner<'_> {
    #[allow(clippy::too_many_lines)] // SIMD inner loop pulls the hot path inline
    fn scan(&self, guides: &[Guide], max_mismatches: u8) -> Vec<Hit> {
        let compare_mask = self.source.enzyme().compare_mask;
        let bin_width = self.source.bin_width();
        let bin_width_bits = u32::from(bin_width) * 2;
        // Saturating shift so an in-principle bin_width of 16 (32 bits)
        // still yields u32::MAX; production bins are capped at 15 by
        // `BinTable::with_width`.
        let bin_mask: u32 = if bin_width_bits >= 32 {
            u32::MAX
        } else {
            (1u32 << bin_width_bits) - 1
        };
        let max_mm_u32 = u32::from(max_mismatches);

        // Precompute each guide's bin-prefix once. Reused on every bin
        // by every worker thread.
        let guide_prefixes: Vec<u32> = guides.iter().map(|g| self.source.bin_key(g.site)).collect();

        // Materialize bins into a Vec so rayon can split the work. The
        // `HashMap::iter` we'd otherwise use isn't a `ParallelIterator`,
        // and `par_bridge` introduces channel overhead per item. For
        // 16 384 bins (width 7) the vec is ~256 KB; trivially worth it.
        let bins: Vec<(u32, &[PackedEntry])> = self.source.iter_bins().collect();

        let trace = std::env::var_os("CRISPR_OTS_SCAN_TRACE").is_some();

        // Per-thread fold: each worker keeps its own scratch buffers and
        // local hit collector. SIMD inner loop processes 4 (guide × site)
        // pairs per iteration; tail handled scalar.
        // Select the SIMD kernel once per scan. AVX-512 VPOPCNTDQ (Genoa /
        // Genoa-X, Ice Lake-SP, Sapphire Rapids, Zen 4+) does the 2-bit-pair
        // popcount in one hardware instruction over eight guides per
        // iteration; every other CPU falls back to the portable `wide::u64x4`
        // SWAR path. `CRISPR_OTS_NO_AVX512` forces the fallback so the two
        // kernels can be A/B-benchmarked on the same binary and node.
        #[cfg(target_arch = "x86_64")]
        let use_avx512 = std::env::var_os("CRISPR_OTS_NO_AVX512").is_none()
            && std::arch::is_x86_feature_detected!("avx512f")
            && std::arch::is_x86_feature_detected!("avx512vpopcntdq");
        let (hits, total_pairs, prefilter_pass) = bins
            .par_iter()
            .fold(
                || {
                    let cap = guides.len() / 16;
                    (
                        Vec::with_capacity(cap),
                        Vec::with_capacity(cap),
                        Vec::with_capacity(cap),
                        Vec::new(),
                        0_u64,
                        0_u64,
                    )
                },
                |(
                    mut active_idx,
                    mut active_lows,
                    mut active_highs,
                    mut local_hits,
                    mut total,
                    mut passed,
                ): WorkerState,
                 (bin_key, sites)| {
                    // ---- Bin-prefix prefilter (gather active guides into SoA). ----
                    active_idx.clear();
                    active_lows.clear();
                    active_highs.clear();
                    for (gi, prefix) in guide_prefixes.iter().enumerate() {
                        total += 1;
                        if bin_prefix_mismatches(*prefix, *bin_key, bin_mask) <= max_mm_u32 {
                            active_idx.push(u32::try_from(gi).expect("guide index fits in u32"));
                            let g_site = guides[gi].site;
                            active_lows.push(g_site.low);
                            active_highs.push(g_site.high);
                            passed += 1;
                        }
                    }
                    if !active_idx.is_empty() {
                        #[cfg(target_arch = "x86_64")]
                        if use_avx512 {
                            // SAFETY: `use_avx512` is only true when this CPU
                            // reported avx512f + avx512vpopcntdq.
                            unsafe {
                                scan_bin_avx512(
                                    &active_lows,
                                    &active_highs,
                                    &active_idx,
                                    sites,
                                    compare_mask,
                                    max_mismatches,
                                    &mut local_hits,
                                );
                            }
                        } else {
                            scan_bin_x4(
                                &active_lows,
                                &active_highs,
                                &active_idx,
                                sites,
                                compare_mask,
                                max_mismatches,
                                &mut local_hits,
                            );
                        }
                        #[cfg(not(target_arch = "x86_64"))]
                        scan_bin_x4(
                            &active_lows,
                            &active_highs,
                            &active_idx,
                            sites,
                            compare_mask,
                            max_mismatches,
                            &mut local_hits,
                        );
                    }
                    (
                        active_idx,
                        active_lows,
                        active_highs,
                        local_hits,
                        total,
                        passed,
                    )
                },
            )
            .map(|(_idx, _lows, _highs, hits, total, passed)| (hits, total, passed))
            .reduce(
                || (Vec::new(), 0_u64, 0_u64),
                |(mut a_hits, a_total, a_pass), (b_hits, b_total, b_pass)| {
                    a_hits.extend(b_hits);
                    (a_hits, a_total + b_total, a_pass + b_pass)
                },
            );

        if trace {
            // Counts can't realistically exceed 2^53 in any genome we'd
            // scan, so the lossy cast is fine for the percentage display.
            #[allow(clippy::cast_precision_loss)]
            let pct = 100.0 * prefilter_pass as f64 / total_pairs.max(1) as f64;
            eprintln!(
                "[scan-trace] {total_pairs} (guide, bin) pairs; \
                 {prefilter_pass} passed prefilter ({pct:.1}%); \
                 {} hits; {} rayon threads",
                hits.len(),
                rayon::current_num_threads(),
            );
        }
        hits
    }
}

impl BinScanner<'_> {
    /// Streaming-aggregate scan: instead of collecting every [`Hit`] into a
    /// `Vec` (which is O(total hits) — hundreds of GB on T2T satellite guides),
    /// fold each hit into a per-thread accumulator `A` inline. `init` makes a
    /// fresh per-thread `A`, `fold_hit` applies one hit, `merge` combines two
    /// `A`s for the reduce. Peak memory is O(#A + one site-chunk of hits per
    /// thread), independent of total hit volume — the CPU analog of the GPU's
    /// in-kernel aggregation. Reuses the exact prefilter + SIMD kernels as
    /// [`Scanner::scan`].
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub fn scan_fold<A, INIT, FOLD, SKIP, MERGE>(
        &self,
        guides: &[Guide],
        max_mismatches: u8,
        init: INIT,
        fold_hit: FOLD,
        skip: SKIP,
        merge: MERGE,
    ) -> A
    where
        A: Send,
        INIT: Fn() -> A + Sync,
        FOLD: Fn(&mut A, &Hit) + Sync,
        // Per-guide early-exit predicate: when it returns true for `(acc, gi)`,
        // guide `gi` is dropped from this bin's prefilter (and every later bin
        // on this thread), so a saturated guide stops being scanned. Pass
        // `|_, _| false` to disable.
        SKIP: Fn(&A, usize) -> bool + Sync,
        MERGE: Fn(A, A) -> A + Sync,
    {
        let compare_mask = self.source.enzyme().compare_mask;
        let bin_width = self.source.bin_width();
        let bin_width_bits = u32::from(bin_width) * 2;
        let bin_mask: u32 = if bin_width_bits >= 32 {
            u32::MAX
        } else {
            (1u32 << bin_width_bits) - 1
        };
        let max_mm_u32 = u32::from(max_mismatches);
        let guide_prefixes: Vec<u32> = guides.iter().map(|g| self.source.bin_key(g.site)).collect();
        let bins: Vec<(u32, &[PackedEntry])> = self.source.iter_bins().collect();
        #[cfg(target_arch = "x86_64")]
        let use_avx512 = std::env::var_os("CRISPR_OTS_NO_AVX512").is_none()
            && std::arch::is_x86_feature_detected!("avx512f")
            && std::arch::is_x86_feature_detected!("avx512vpopcntdq");
        // Drain the transient per-bin hit buffer every SITE_CHUNK sites so a
        // single low-complexity (satellite) bin can't materialize millions of
        // hits at once — keeps peak memory bounded by the chunk, not the bin.
        const SITE_CHUNK: usize = 16_384;

        bins.par_iter()
            .fold(
                || {
                    (
                        Vec::<u32>::new(),
                        Vec::<u64>::new(),
                        Vec::<u64>::new(),
                        Vec::<Hit>::new(),
                        init(),
                    )
                },
                |(mut active_idx, mut active_lows, mut active_highs, mut local_hits, mut acc),
                 (bin_key, sites)| {
                    active_idx.clear();
                    active_lows.clear();
                    active_highs.clear();
                    for (gi, prefix) in guide_prefixes.iter().enumerate() {
                        // Skip guides that have already saturated on this thread —
                        // the cap-aware early-exit that keeps satellite guides from
                        // dominating the scan.
                        if skip(&acc, gi) {
                            continue;
                        }
                        if bin_prefix_mismatches(*prefix, *bin_key, bin_mask) <= max_mm_u32 {
                            active_idx.push(u32::try_from(gi).expect("guide index fits in u32"));
                            let g_site = guides[gi].site;
                            active_lows.push(g_site.low);
                            active_highs.push(g_site.high);
                        }
                    }
                    if !active_idx.is_empty() {
                        for chunk in sites.chunks(SITE_CHUNK) {
                            #[cfg(target_arch = "x86_64")]
                            if use_avx512 {
                                // SAFETY: `use_avx512` is only true when the CPU
                                // reported avx512f + avx512vpopcntdq.
                                unsafe {
                                    scan_bin_avx512(
                                        &active_lows,
                                        &active_highs,
                                        &active_idx,
                                        chunk,
                                        compare_mask,
                                        max_mismatches,
                                        &mut local_hits,
                                    );
                                }
                            } else {
                                scan_bin_x4(
                                    &active_lows,
                                    &active_highs,
                                    &active_idx,
                                    chunk,
                                    compare_mask,
                                    max_mismatches,
                                    &mut local_hits,
                                );
                            }
                            #[cfg(not(target_arch = "x86_64"))]
                            scan_bin_x4(
                                &active_lows,
                                &active_highs,
                                &active_idx,
                                chunk,
                                compare_mask,
                                max_mismatches,
                                &mut local_hits,
                            );
                            for h in local_hits.drain(..) {
                                fold_hit(&mut acc, &h);
                            }
                        }
                    }
                    (active_idx, active_lows, active_highs, local_hits, acc)
                },
            )
            .map(|(_, _, _, _, acc)| acc)
            .reduce(|| init(), |a, b| merge(a, b))
    }

    /// Serial (single-thread) sibling of [`scan_fold`] for ONE guide chunk:
    /// walks every bin **in sequence**, folding each hit into a single `acc` via
    /// `fold_hit`. Because bins are processed serially into one accumulator, a
    /// guide's running count is complete as it grows, so `skip` is an **exact**
    /// per-guide early-exit — a saturated guide is dropped from every later bin,
    /// and `fold_hit` (which writes the off-target) is never called past the cap.
    /// The caller fans `par_iter` over guide chunks, each chunk owning its own
    /// `acc` + output shard: partition-by-guides instead of by-bins, giving an
    /// exact cap and streamed (non-materialized) Mode-2 output. `fold_hit` is
    /// `FnMut` (it writes), so this method itself is single-threaded.
    pub fn scan_chunk_fold<A, FOLD, SKIP>(
        &self,
        guides: &[Guide],
        max_mismatches: u8,
        acc: &mut A,
        mut fold_hit: FOLD,
        skip: SKIP,
    ) where
        FOLD: FnMut(&mut A, &Hit),
        SKIP: Fn(&A, usize) -> bool,
    {
        let compare_mask = self.source.enzyme().compare_mask;
        let bin_width = self.source.bin_width();
        let bin_width_bits = u32::from(bin_width) * 2;
        let bin_mask: u32 = if bin_width_bits >= 32 {
            u32::MAX
        } else {
            (1u32 << bin_width_bits) - 1
        };
        let max_mm_u32 = u32::from(max_mismatches);
        let guide_prefixes: Vec<u32> = guides.iter().map(|g| self.source.bin_key(g.site)).collect();
        #[cfg(target_arch = "x86_64")]
        let use_avx512 = std::env::var_os("CRISPR_OTS_NO_AVX512").is_none()
            && std::arch::is_x86_feature_detected!("avx512f")
            && std::arch::is_x86_feature_detected!("avx512vpopcntdq");
        const SITE_CHUNK: usize = 16_384;
        let mut active_idx: Vec<u32> = Vec::new();
        let mut active_lows: Vec<u64> = Vec::new();
        let mut active_highs: Vec<u64> = Vec::new();
        let mut local_hits: Vec<Hit> = Vec::new();
        for (bin_key, sites) in self.source.iter_bins() {
            active_idx.clear();
            active_lows.clear();
            active_highs.clear();
            for (gi, prefix) in guide_prefixes.iter().enumerate() {
                // Exact early-exit: a guide saturated on this (single) accumulator
                // is dropped from this bin and every later one.
                if skip(acc, gi) {
                    continue;
                }
                if bin_prefix_mismatches(*prefix, bin_key, bin_mask) <= max_mm_u32 {
                    active_idx.push(u32::try_from(gi).expect("guide index fits in u32"));
                    let g_site = guides[gi].site;
                    active_lows.push(g_site.low);
                    active_highs.push(g_site.high);
                }
            }
            if active_idx.is_empty() {
                continue;
            }
            for chunk in sites.chunks(SITE_CHUNK) {
                #[cfg(target_arch = "x86_64")]
                if use_avx512 {
                    // SAFETY: `use_avx512` is only true when the CPU reported
                    // avx512f + avx512vpopcntdq.
                    unsafe {
                        scan_bin_avx512(
                            &active_lows,
                            &active_highs,
                            &active_idx,
                            chunk,
                            compare_mask,
                            max_mismatches,
                            &mut local_hits,
                        );
                    }
                } else {
                    scan_bin_x4(
                        &active_lows,
                        &active_highs,
                        &active_idx,
                        chunk,
                        compare_mask,
                        max_mismatches,
                        &mut local_hits,
                    );
                }
                #[cfg(not(target_arch = "x86_64"))]
                scan_bin_x4(
                    &active_lows,
                    &active_highs,
                    &active_idx,
                    chunk,
                    compare_mask,
                    max_mismatches,
                    &mut local_hits,
                );
                for h in local_hits.drain(..) {
                    fold_hit(acc, &h);
                }
            }
        }
    }
}

/// Count base-level mismatches between two bin-prefix u32 keys (`mask`
/// covers only the lower `bin_width × 2` bits). Same POPCNT-pair trick as
/// `mismatches_masked` for the 2 × `u64` site encoding, specialized to
/// `u32`.
#[inline]
fn bin_prefix_mismatches(guide_prefix: u32, bin_key: u32, mask: u32) -> u32 {
    const UPPER_PAIR_BITS: u32 = 0xAAAA_AAAA;
    let xor = (guide_prefix ^ bin_key) & mask;
    ((xor & UPPER_PAIR_BITS) | ((xor << 1) & UPPER_PAIR_BITS)).count_ones()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crispr_db::{read_fasta_from, BinTable, SiteFinder, Strand};
    use crispr_encoding::Enzyme;

    /// Helper: build a `BinTable` from raw FASTA bytes for a given enzyme.
    fn build_table(fasta: &[u8], enzyme: Enzyme) -> BinTable {
        let contigs = read_fasta_from(fasta).expect("parses FASTA");
        let finder = SiteFinder::new(enzyme.clone());
        let mut table = BinTable::for_enzyme(enzyme);
        for (i, contig) in contigs.iter().enumerate() {
            finder.scan(
                &contig.sequence,
                u32::try_from(i).expect("contig fits in u32"),
                &mut table,
            );
        }
        table
    }

    #[test]
    fn empty_table_yields_no_hits() {
        let table = BinTable::for_enzyme(Enzyme::spcas9_ngg());
        let scanner = BinScanner::new(&table);
        let guide = Guide {
            id: "g0".into(),
            site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG"),
        };
        let hits = scanner.scan(&[guide], 4);
        assert!(hits.is_empty());
    }

    #[test]
    fn empty_guide_list_yields_no_hits() {
        // Even a non-empty table yields no hits with no guides.
        let table = build_table(b">c\nAAAAAAAAAAAAAAAAAAAAAGG\n", Enzyme::spcas9_ngg());
        assert!(table.total() > 0);
        let scanner = BinScanner::new(&table);
        let hits = scanner.scan(&[], 4);
        assert!(hits.is_empty());
    }

    #[test]
    fn guide_finds_its_own_on_target_at_zero_mismatches() {
        let table = build_table(b">c\nAAAAAAAAAAAAAAAAAAAAAGG\n", Enzyme::spcas9_ngg());
        let scanner = BinScanner::new(&table);
        let guide = Guide {
            id: "self".into(),
            site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG"),
        };
        let hits = scanner.scan(&[guide], 0);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].mismatches, 0);
        assert_eq!(hits[0].guide_index, 0);
        assert_eq!(hits[0].position.strand, Strand::Forward);
        assert_eq!(hits[0].position.offset, 0);
    }

    #[test]
    fn guide_finds_off_target_within_budget() {
        // Two sites in the genome:
        //   - on-target at offset 0:        "AAAAAAAAAAAAAAAAAAAA AGG" (matches guide)
        //   - off-target at offset 30:      "AAAAAAAAAAAAAAAAAAAT AGG" (1 mm at pos 19)
        // Plus padding to keep them distinct.
        let fasta = b">c\n\
                      AAAAAAAAAAAAAAAAAAAAAGG\
                      CCCCCCC\
                      AAAAAAAAAAAAAAAAAAATAGG\n";
        let table = build_table(fasta, Enzyme::spcas9_ngg());
        let scanner = BinScanner::new(&table);
        let guide = Guide {
            id: "g".into(),
            site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG"),
        };

        // With budget 0: only the exact match.
        let hits_strict = scanner.scan(std::slice::from_ref(&guide), 0);
        assert_eq!(hits_strict.len(), 1);
        assert_eq!(hits_strict[0].mismatches, 0);

        // With budget 1: both sites — the off-target has one mismatch.
        let hits_lax = scanner.scan(std::slice::from_ref(&guide), 1);
        assert_eq!(hits_lax.len(), 2);
        let mut mm: Vec<u8> = hits_lax.iter().map(|h| h.mismatches).collect();
        mm.sort_unstable();
        assert_eq!(mm, vec![0, 1]);
    }

    #[test]
    fn multiple_guides_get_distinct_hit_indices() {
        let table = build_table(b">c\nAAAAAAAAAAAAAAAAAAAAAGG\n", Enzyme::spcas9_ngg());
        let scanner = BinScanner::new(&table);
        let g0 = Guide {
            id: "g0".into(),
            site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG"),
        };
        let g1 = Guide {
            id: "g1".into(),
            site: Site::encode_ascii(b"TTTTTTTTTTTTTTTTTTTTAGG"), // very different
        };
        let hits = scanner.scan(&[g0, g1], 0);
        // Only g0 should match.
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].guide_index, 0);
    }

    #[test]
    fn pam_difference_alone_does_not_count() {
        // Same protospacer, different PAM N (AGG vs TGG). The compare-mask
        // excludes PAM bits, so this is a 0-mismatch hit.
        let table = build_table(b">c\nAAAAAAAAAAAAAAAAAAAAAGG\n", Enzyme::spcas9_ngg());
        let scanner = BinScanner::new(&table);
        let guide = Guide {
            id: "g".into(),
            site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAATGG"),
        };
        let hits = scanner.scan(&[guide], 0);
        assert_eq!(hits.len(), 1, "PAM N differs but is masked out");
        assert_eq!(hits[0].mismatches, 0);
    }

    #[test]
    fn reverse_strand_site_is_found() {
        // A reverse-strand Cas9 NGG site: forward strand has CCT + 20 ACGT.
        let fasta = b">c\nCCTACGTACGTACGTACGTACGT\n";
        let table = build_table(fasta, Enzyme::spcas9_ngg());

        // Build a guide whose protospacer matches the reverse-strand
        // protospacer (i.e. the rev-comp of "ACGTACGTACGTACGTACGT").
        // rc("ACGTACGTACGTACGTACGT") = "ACGTACGTACGTACGTACGT" — this
        // sequence happens to be its own rev-comp, so we can just use it.
        let guide = Guide {
            id: "g".into(),
            site: Site::encode_ascii(b"ACGTACGTACGTACGTACGTAGG"),
        };

        let scanner = BinScanner::new(&table);
        let hits = scanner.scan(&[guide], 0);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].position.strand, Strand::Reverse);
        assert_eq!(hits[0].mismatches, 0);
    }

    #[test]
    fn cpf1_end_to_end() {
        // TTT + N=A + 23-mer ACGTACGTACGTACGTACGTACG (27 chars total).
        let table = build_table(b">c\nTTTAACGTACGTACGTACGTACGTACG\n", Enzyme::cpf1_tttn());
        let scanner = BinScanner::new(&table);
        let guide = Guide {
            id: "g".into(),
            site: Site::encode_ascii(b"TTTAACGTACGTACGTACGTACGTACG"),
        };
        let hits = scanner.scan(&[guide], 0);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].mismatches, 0);
        assert_eq!(hits[0].position.strand, Strand::Forward);
    }

    /// The AVX-512 8-wide kernel must produce byte-identical hits to the
    /// portable x4 kernel on the same inputs. The small scan tests above only
    /// exercise the scalar tail (1–2 guides → `n_chunks == 0`); this one drives
    /// full 8-lane chunks plus a remainder. Skips on CPUs/builds without
    /// VPOPCNTDQ (e.g. the Zen 2 GPU hosts), where the dispatch never picks it.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx512_kernel_matches_x4_kernel() {
        use crispr_db::{PackedEntry, Position};

        if std::env::var_os("CRISPR_OTS_NO_AVX512").is_some()
            || !std::arch::is_x86_feature_detected!("avx512f")
            || !std::arch::is_x86_feature_detected!("avx512vpopcntdq")
        {
            return;
        }

        // xorshift64 — deterministic, no external rng dependency.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        // 50 random sites; full-width random Site words (both kernels apply the
        // same compare-mask internally, so unrealistic count/reserved bits are
        // harmless and still exercise every lane).
        let sites: Vec<PackedEntry> = (0..50u32)
            .map(|i| {
                PackedEntry::from_site_position(
                    Site {
                        low: next(),
                        high: next(),
                    },
                    Position {
                        contig_id: 0,
                        offset: i,
                        strand: Strand::Forward,
                    },
                )
            })
            .collect();

        // 21 guides → two full 8-lane chunks + a 5-guide scalar tail (and, for
        // the x4 kernel, five 4-lane chunks + a 1-guide tail).
        let active_idx: Vec<u32> = (0..21).collect();
        let active_lows: Vec<u64> = (0..21).map(|_| next()).collect();
        let active_highs: Vec<u64> = (0..21).map(|_| next()).collect();

        let sort_key = |h: &Hit| (h.guide_index, h.position.offset, h.mismatches);
        for enzyme in [Enzyme::spcas9_ngg(), Enzyme::cpf1_tttn()] {
            let mask = enzyme.compare_mask;
            for max_mm in [0u8, 1, 2, 4, 8] {
                let mut hits_x4 = Vec::new();
                scan_bin_x4(
                    &active_lows,
                    &active_highs,
                    &active_idx,
                    &sites,
                    mask,
                    max_mm,
                    &mut hits_x4,
                );
                let mut hits_x8 = Vec::new();
                // SAFETY: VPOPCNTDQ confirmed by the guard at the top of the test.
                unsafe {
                    scan_bin_avx512(
                        &active_lows,
                        &active_highs,
                        &active_idx,
                        &sites,
                        mask,
                        max_mm,
                        &mut hits_x8,
                    );
                }
                let mut keys_x4: Vec<_> = hits_x4.iter().map(sort_key).collect();
                let mut keys_x8: Vec<_> = hits_x8.iter().map(sort_key).collect();
                keys_x4.sort_unstable();
                keys_x8.sort_unstable();
                assert_eq!(
                    keys_x4, keys_x8,
                    "AVX-512 and x4 kernels disagree (enzyme={}, max_mm={max_mm})",
                    enzyme.name
                );
            }
        }
    }
}
