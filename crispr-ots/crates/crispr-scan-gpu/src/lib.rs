//! GPU (CUDA via `cudarc`) implementation of the [`crispr_scan::Scanner`]
//! trait.
//!
//! The off-target database is a flat array of 24-byte `PackedEntry` records
//! (`site_low: u64`, `site_high: u64`, `packed_position: u64`). We upload that
//! array once to VRAM — it stays resident across queries (a mouse `.crot` is
//! ~6.7 GB, comfortably inside a 24 GB A5500 or a 40/80 GB A100) — then for
//! each scan upload the (small) guide set, launch one thread per database
//! entry that inner-loops over every guide computing the same
//! XOR / compare-mask / 2-bit-pair popcount as the CPU `BinScanner`, and
//! atomic-append the matches to a device buffer.
//!
//! The bin-prefix prefilter is intentionally skipped: brute force is correct
//! (it visits exactly the entries the prefilter would have kept plus ones that
//! cannot match, so the hit *set* is identical), each entry is read from VRAM
//! once (the small guide arrays stay L2-resident), and the GPU's throughput
//! makes the extra full-compares cheap. Everything else — `crispr-encoding`,
//! `crispr-db`, `crispr-score`, the threshold/cap logic, and the CSV/TSV
//! writers — is reused unchanged; only the mismatch kernel is re-expressed in
//! CUDA C.

// u32<->usize length/index casts are pervasive in the host glue and lossless
// on the 64-bit targets we support; hardware names (A5500, VRAM, NVRTC, …)
// trip the doc linter.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::doc_markdown
)]

use std::sync::Arc;

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::{compile_ptx_with_opts, CompileOptions};

use crispr_db::{BinSource, PackedEntry};
use crispr_encoding::SiteMask;
use crispr_scan::{Guide, Hit, Scanner};

/// CUDA C source for the brute-force mismatch kernel. One thread per database
/// entry; each thread scans every guide. Mirrors
/// `crispr_encoding::mismatches_masked` exactly: the differing-2-bit-pair
/// indicator `(c & 0xAAAA…) | ((c << 1) & 0xAAAA…)` is population-counted with
/// the hardware `__popcll`, summed over the `low` and `high` site words.
const KERNEL_SRC: &str = r#"
extern "C" __global__ void scan_kernel(
    const unsigned long long* __restrict__ entries, // 3 u64/entry: low, high, pos
    const unsigned int n_entries,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high,
    const unsigned int n_guides,
    const unsigned long long mask_low,
    const unsigned long long mask_high,
    const unsigned int max_mm,
    unsigned int* __restrict__ out_count,
    const unsigned int out_cap,
    unsigned int* __restrict__ out_guide,
    unsigned int* __restrict__ out_entry,
    unsigned int* __restrict__ out_mm)
{
    const unsigned long long ei =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (ei >= n_entries) return;
    const unsigned long long e_low  = entries[3 * ei];
    const unsigned long long e_high = entries[3 * ei + 1];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int gi = 0; gi < n_guides; ++gi) {
        const unsigned long long cl = (guide_low[gi]  ^ e_low)  & mask_low;
        const unsigned long long ch = (guide_high[gi] ^ e_high) & mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= max_mm) {
            const unsigned int slot = atomicAdd(out_count, 1u);
            if (slot < out_cap) {
                out_guide[slot] = gi;
                out_entry[slot] = (unsigned int)ei;
                out_mm[slot] = mm;
            }
        }
    }
}

// Genome-scale path: accumulate per-guide off-target COUNTS by mismatch level
// (atomicAdd into out_counts[gi*(max_mm+1)+mm]) instead of emitting hits. The
// output is n_guides*(max_mm+1) u32 regardless of hit volume — the right shape
// for all-vs-all specificity over millions of guides (a full hit list would be
// petabytes). Note the mm==0 bucket includes each guide's own on-target site.
extern "C" __global__ void scan_count_kernel(
    const unsigned long long* __restrict__ entries, const unsigned int n_entries,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high, const unsigned int n_guides,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, unsigned int* __restrict__ out_counts)
{
    const unsigned long long ei =
        (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (ei >= n_entries) return;
    const unsigned long long e_low = entries[3 * ei];
    const unsigned long long e_high = entries[3 * ei + 1];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    const unsigned int stride = max_mm + 1;
    for (unsigned int gi = 0; gi < n_guides; ++gi) {
        const unsigned long long cl = (guide_low[gi]  ^ e_low)  & mask_low;
        const unsigned long long ch = (guide_high[gi] ^ e_high) & mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= max_mm) {
            atomicAdd(&out_counts[gi * stride + mm], 1u);
        }
    }
}

// ---- Bin-prefilter path ----
// Each guide visits only its candidate bins (bins whose bin-key is within
// max_mm base-mismatches of the guide's bin-key), accumulating per-guide counts
// in registers (no atomics). Candidate bins are enumerated by perturbing 0..k
// of the bin-key's base positions. Same counts as scan_count_kernel, but it
// skips the ~99% of bins that cannot contain an off-target.

// Extract `num_bits` bits starting at `start_bit` from the 2x u64 site — mirror
// of the Rust `extract_bits` (low word holds the low 48 bits of sequence).
__device__ __forceinline__ unsigned int bin_key_of(
    unsigned long long low, unsigned long long high,
    unsigned int start_bit, unsigned int num_bits)
{
    const unsigned int LOW_WIDTH = 48u;
    const unsigned int end_bit = start_bit + num_bits;
    const unsigned long long m =
        (num_bits >= 64u) ? ~0ULL : ((1ULL << num_bits) - 1ULL);
    unsigned long long raw;
    if (end_bit <= LOW_WIDTH) {
        raw = (low >> start_bit) & m;
    } else if (start_bit >= LOW_WIDTH) {
        raw = (high >> (start_bit - LOW_WIDTH)) & m;
    } else {
        const unsigned int low_count = LOW_WIDTH - start_bit;
        const unsigned long long low_part = low >> start_bit;
        const unsigned long long high_part =
            (high & ((1ULL << (num_bits - low_count)) - 1ULL)) << low_count;
        raw = low_part | high_part;
    }
    return (unsigned int)raw;
}

struct ScanCtx {
    const unsigned long long* entries;
    const unsigned int* bin_off; // [2*key]=start, [2*key+1]=count
    unsigned long long e_low, e_high, mask_low, mask_high;
    unsigned int max_mm;
    unsigned int* local; // [max_mm+1] per-thread counters
};

__device__ void process_bin(const ScanCtx* c, unsigned int key) {
    const unsigned int start = c->bin_off[2u * key];
    const unsigned int count = c->bin_off[2u * key + 1u];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int j = 0; j < count; ++j) {
        const unsigned long long idx = (unsigned long long)(start + j);
        const unsigned long long o_low = c->entries[3 * idx];
        const unsigned long long o_high = c->entries[3 * idx + 1];
        const unsigned long long cl = (o_low ^ c->e_low) & c->mask_low;
        const unsigned long long ch = (o_high ^ c->e_high) & c->mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= c->max_mm) {
            c->local[mm] += 1u;
        }
    }
}

// Visit `key` and every bin reachable by mutating <= k_rem base positions at
// index >= start_pos. Strictly increasing positions => each candidate visited
// exactly once.
__device__ void visit_bins(
    const ScanCtx* c, unsigned int key, int w, int k_rem, int start_pos)
{
    process_bin(c, key);
    if (k_rem == 0) return;
    for (int p = start_pos; p < w; ++p) {
        const unsigned int shift = 2u * (unsigned int)p;
        for (unsigned int delta = 1u; delta <= 3u; ++delta) {
            visit_bins(c, key ^ (delta << shift), w, k_rem - 1, p + 1);
        }
    }
}

extern "C" __global__ void scan_count_prefilter_kernel(
    const unsigned long long* __restrict__ entries,
    const unsigned int* __restrict__ bin_off,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high,
    const unsigned int n_guides,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, const unsigned int bin_w,
    const unsigned int bin_start_bit, unsigned int* __restrict__ out_counts)
{
    const unsigned int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= n_guides) return;
    const unsigned long long g_low = guide_low[gi];
    const unsigned long long g_high = guide_high[gi];
    const unsigned int stride = max_mm + 1u;

    unsigned int local[8];
    for (unsigned int m = 0; m < stride; ++m) local[m] = 0u;

    ScanCtx c;
    c.entries = entries;
    c.bin_off = bin_off;
    c.e_low = g_low;
    c.e_high = g_high;
    c.mask_low = mask_low;
    c.mask_high = mask_high;
    c.max_mm = max_mm;
    c.local = local;

    const unsigned int key = bin_key_of(g_low, g_high, bin_start_bit, 2u * bin_w);
    visit_bins(&c, key, (int)bin_w, (int)max_mm, 0);

    for (unsigned int m = 0; m < stride; ++m) {
        out_counts[gi * stride + m] = local[m];
    }
}

// ---- Bin-prefilter HIT-EMIT path ----
// Same per-guide candidate-bin traversal as scan_count_prefilter_kernel, but
// each match is APPENDED to the output buffer (guide, entry, mm) via a global
// atomic counter — the fast genome-scale equivalent of the brute-force
// scan_kernel. The emitted hit set is identical to scan_kernel's (the prefilter
// only skips bins that cannot contain a match); it's just found ~40x faster by
// not touching the ~99% of non-candidate bins.
struct EmitCtx {
    const unsigned long long* entries;
    const unsigned int* bin_off;
    unsigned long long e_low, e_high, mask_low, mask_high;
    unsigned int max_mm;
    unsigned int gi;
    unsigned int* out_count;
    unsigned int out_cap;
    unsigned int cap;            // per-guide emit cap (on+off); 0 = uncapped
    unsigned int* out_guide;
    unsigned int* out_entry;
    unsigned int* out_mm;
};

__device__ void process_bin_emit(const EmitCtx* c, unsigned int key, unsigned int* emitted) {
    const unsigned int start = c->bin_off[2u * key];
    const unsigned int count = c->bin_off[2u * key + 1u];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int j = 0; j < count; ++j) {
        const unsigned long long idx = (unsigned long long)(start + j);
        const unsigned long long o_low = c->entries[3 * idx];
        const unsigned long long o_high = c->entries[3 * idx + 1];
        const unsigned long long cl = (o_low ^ c->e_low) & c->mask_low;
        const unsigned long long ch = (o_high ^ c->e_high) & c->mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= c->max_mm) {
            // Per-guide cap (on+off): once this guide has emitted `cap` hits,
            // stop — mirrors the CPU/CFD-kernel early-exit, bounds the device
            // buffer to n_guides * cap, and lets satellite guides quit early.
            if (c->cap != 0u && *emitted >= c->cap) return;
            *emitted += 1u;
            const unsigned int slot = atomicAdd(c->out_count, 1u);
            if (slot < c->out_cap) {
                c->out_guide[slot] = c->gi;
                c->out_entry[slot] = (unsigned int)idx;
                c->out_mm[slot] = mm;
            }
        }
    }
}

__device__ void visit_bins_emit(
    const EmitCtx* c, unsigned int key, int w, int k_rem, int start_pos, unsigned int* emitted)
{
    if (c->cap != 0u && *emitted >= c->cap) return;   // guide saturated: stop visiting bins
    process_bin_emit(c, key, emitted);
    if (k_rem == 0) return;
    for (int p = start_pos; p < w; ++p) {
        const unsigned int shift = 2u * (unsigned int)p;
        for (unsigned int delta = 1u; delta <= 3u; ++delta) {
            if (c->cap != 0u && *emitted >= c->cap) return;
            visit_bins_emit(c, key ^ (delta << shift), w, k_rem - 1, p + 1, emitted);
        }
    }
}

extern "C" __global__ void scan_prefilter_emit_kernel(
    const unsigned long long* __restrict__ entries,
    const unsigned int* __restrict__ bin_off,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high,
    const unsigned int n_guides,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, const unsigned int bin_w,
    const unsigned int bin_start_bit,
    unsigned int* __restrict__ out_count, const unsigned int out_cap,
    unsigned int* __restrict__ out_guide,
    unsigned int* __restrict__ out_entry,
    unsigned int* __restrict__ out_mm,
    const unsigned int cap_per_guide)
{
    const unsigned int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= n_guides) return;
    const unsigned long long g_low = guide_low[gi];
    const unsigned long long g_high = guide_high[gi];

    EmitCtx c;
    c.entries = entries;
    c.bin_off = bin_off;
    c.e_low = g_low;
    c.e_high = g_high;
    c.mask_low = mask_low;
    c.mask_high = mask_high;
    c.max_mm = max_mm;
    c.gi = gi;
    c.out_count = out_count;
    c.out_cap = out_cap;
    c.cap = cap_per_guide;
    c.out_guide = out_guide;
    c.out_entry = out_entry;
    c.out_mm = out_mm;

    unsigned int emitted = 0u;   // per-guide (per-thread) emitted count for the cap
    const unsigned int key = bin_key_of(g_low, g_high, bin_start_bit, 2u * bin_w);
    visit_bins_emit(&c, key, (int)bin_w, (int)max_mm, 0, &emitted);
}

// ---- Bin-prefilter CFD-ACCUMULATION path (Mode-1 specificity) ----
// Same per-guide candidate-bin traversal as scan_count_prefilter_kernel, but
// instead of bucketing counts it computes the SpCas9 CFD of each off-target on
// the fly and accumulates, per guide, the off-target CFD sum, the max CFD, and
// the on-target (mm==0) multiplicity — in registers, no atomics. Output is
// 3 scalars per guide regardless of hit volume, so NO hit list ever leaves the
// GPU: this is the only path that scales Mode-1 to a full genome of repeat-rich
// guides (hundreds of billions of hits would otherwise overflow the emit
// buffer). CFD itself is a verbatim port of crispr_score::Cfd::score_pair: a
// 23-bp Cas9 site lives entirely in the `low` word, so only `o_low`/`g_low`
// feed the score; the penalty[4][4][20] and pam[4][4] tables are uploaded once
// and passed as read-only pointers (flattened g*80+o*20+p and b1*4+b2).

__device__ __forceinline__ double cfd_score_pair(
    unsigned long long g_low, unsigned long long o_low,
    const double* __restrict__ penalty, const double* __restrict__ pam)
{
    double s = 1.0;
    for (unsigned int p = 0; p < 20u; ++p) {
        const unsigned int bit = (22u - p) * 2u;
        const unsigned int g = (unsigned int)((g_low >> bit) & 3ULL);
        const unsigned int o = (unsigned int)((o_low >> bit) & 3ULL);
        s *= penalty[g * 80u + o * 20u + p];
    }
    const unsigned int pam1 = (unsigned int)((o_low >> 2) & 3ULL);
    const unsigned int pam2 = (unsigned int)(o_low & 3ULL);
    return s * pam[pam1 * 4u + pam2];
}

struct CfdCtx {
    const unsigned long long* entries;
    const unsigned int* bin_off;
    const double* penalty;
    const double* pam;
    unsigned long long g_low, g_high, mask_low, mask_high;
    unsigned int max_mm;
    double off_sum;
    double max_cfd;
    unsigned int on_count;
    unsigned int off_count;
    unsigned int cap;        // 0 = unlimited; else saturate at on+off >= cap
    double off_sum_cap;      // 0 = disabled; else saturate once off_sum > this
                             //   (= 1/min_specificity - 1: spec provably below floor)
    unsigned int saturated;  // 1 once either cap was hit (partial sums; flagged)
};

// __noinline__ keeps this (heavy) leaf out of the recursive cfd_visit frame, so
// the per-thread stack that accumulates with recursion depth stays as small as
// the count kernel's — it is called and returns before each recursive descent.
__device__ __noinline__ void cfd_process_bin(CfdCtx* c, unsigned int key) {
    const unsigned int start = c->bin_off[2u * key];
    const unsigned int count = c->bin_off[2u * key + 1u];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int j = 0; j < count; ++j) {
        const unsigned long long idx = (unsigned long long)(start + j);
        const unsigned long long o_low = c->entries[3 * idx];
        const unsigned long long o_high = c->entries[3 * idx + 1];
        const unsigned long long cl = (o_low ^ c->g_low) & c->mask_low;
        const unsigned long long ch = (o_high ^ c->g_high) & c->mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= c->max_mm) {
            if (mm == 0u) {
                c->on_count += 1u;
            } else {
                const double v = cfd_score_pair(c->g_low, o_low, c->penalty, c->pam);
                c->off_sum += v;
                if (v > c->max_cfd) c->max_cfd = v;
                c->off_count += 1u;
            }
            // Bound per-guide work + flag provably non-specific guides. Stop on
            // EITHER: too many hits (count cap), or the off-target CFD sum so
            // high that specificity is already below the floor (off_sum cap —
            // the principled cut; exits the worst guides almost immediately).
            // Partial sums; host reports the guide as saturated (inexact).
            if ((c->cap != 0u && (c->on_count + c->off_count) >= c->cap) ||
                (c->off_sum_cap > 0.0 && c->off_sum > c->off_sum_cap)) {
                c->saturated = 1u;
                return;
            }
        }
    }
}

__device__ void cfd_visit(
    CfdCtx* c, unsigned int key, int w, int k_rem, int start_pos)
{
    if (c->saturated) return;
    cfd_process_bin(c, key);
    if (k_rem == 0 || c->saturated) return;
    for (int p = start_pos; p < w; ++p) {
        const unsigned int shift = 2u * (unsigned int)p;
        for (unsigned int delta = 1u; delta <= 3u; ++delta) {
            cfd_visit(c, key ^ (delta << shift), w, k_rem - 1, p + 1);
            if (c->saturated) return;
        }
    }
}

extern "C" __global__ void scan_cfd_prefilter_kernel(
    const unsigned long long* __restrict__ entries,
    const unsigned int* __restrict__ bin_off,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high,
    const unsigned int n_guides,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, const unsigned int bin_w,
    const unsigned int bin_start_bit,
    const double* __restrict__ penalty, const double* __restrict__ pam,
    const unsigned int cap, const double off_sum_cap,
    double* __restrict__ out_sum, double* __restrict__ out_max,
    unsigned int* __restrict__ out_on, unsigned int* __restrict__ out_off,
    unsigned int* __restrict__ out_sat)
{
    const unsigned int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= n_guides) return;

    CfdCtx c;
    c.entries = entries;
    c.bin_off = bin_off;
    c.penalty = penalty;
    c.pam = pam;
    c.g_low = guide_low[gi];
    c.g_high = guide_high[gi];
    c.mask_low = mask_low;
    c.mask_high = mask_high;
    c.max_mm = max_mm;
    c.off_sum = 0.0;
    c.max_cfd = 0.0;
    c.on_count = 0u;
    c.off_count = 0u;
    c.cap = cap;
    c.off_sum_cap = off_sum_cap;
    c.saturated = 0u;

    const unsigned int key = bin_key_of(c.g_low, c.g_high, bin_start_bit, 2u * bin_w);
    cfd_visit(&c, key, (int)bin_w, (int)max_mm, 0);

    out_sum[gi] = c.off_sum;
    out_max[gi] = c.max_cfd;
    out_on[gi] = c.on_count;
    out_off[gi] = c.off_count;
    out_sat[gi] = c.saturated;
}

// ---- Bin-prefilter SCREEN path (early-exit specificity pre-screen) ----
// One thread per guide; traverse candidate bins at mm<=t and BAIL the instant
// the guide is proven non-specific: a 0-mm DUPLICATE (>=2 zero-mismatch hits —
// the on-target plus at least one more copy) or ANY 1..=t-mm off-target. Most
// non-specific guides exit almost immediately (a near-match sits in their own
// d=0 bin), so this is far cheaper than counting every hit. Output: per guide,
// 1 = drop (non-specific), 0 = keep.
struct ScreenCtx {
    const unsigned long long* entries;
    const unsigned int* bin_off;
    unsigned long long g_low, g_high, mask_low, mask_high;
    unsigned int max_mm;     // = t
    unsigned int on_count;   // zero-mismatch hits seen so far
    unsigned int dropped;    // 1 once proven non-specific
};

__device__ __noinline__ void screen_process_bin(ScreenCtx* c, unsigned int key) {
    const unsigned int start = c->bin_off[2u * key];
    const unsigned int count = c->bin_off[2u * key + 1u];
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int j = 0; j < count; ++j) {
        const unsigned long long idx = (unsigned long long)(start + j);
        const unsigned long long o_low = c->entries[3 * idx];
        const unsigned long long o_high = c->entries[3 * idx + 1];
        const unsigned long long cl = (o_low ^ c->g_low) & c->mask_low;
        const unsigned long long ch = (o_high ^ c->g_high) & c->mask_high;
        const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
        const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
        const unsigned int mm = __popcll(il) + __popcll(ih);
        if (mm <= c->max_mm) {
            if (mm == 0u) {
                c->on_count += 1u;
                if (c->on_count >= 2u) { c->dropped = 1u; return; }
            } else {
                c->dropped = 1u; // any 1..=t-mm off-target ⇒ non-specific
                return;
            }
        }
    }
}

__device__ void screen_visit(
    ScreenCtx* c, unsigned int key, int w, int k_rem, int start_pos)
{
    if (c->dropped) return;
    screen_process_bin(c, key);
    if (k_rem == 0 || c->dropped) return;
    for (int p = start_pos; p < w; ++p) {
        const unsigned int shift = 2u * (unsigned int)p;
        for (unsigned int delta = 1u; delta <= 3u; ++delta) {
            screen_visit(c, key ^ (delta << shift), w, k_rem - 1, p + 1);
            if (c->dropped) return;
        }
    }
}

extern "C" __global__ void scan_screen_kernel(
    const unsigned long long* __restrict__ entries,
    const unsigned int* __restrict__ bin_off,
    const unsigned long long* __restrict__ guide_low,
    const unsigned long long* __restrict__ guide_high,
    const unsigned int n_guides,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, const unsigned int bin_w,
    const unsigned int bin_start_bit, unsigned int* __restrict__ out_dropped)
{
    const unsigned int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= n_guides) return;
    ScreenCtx c;
    c.entries = entries;
    c.bin_off = bin_off;
    c.g_low = guide_low[gi];
    c.g_high = guide_high[gi];
    c.mask_low = mask_low;
    c.mask_high = mask_high;
    c.max_mm = max_mm;
    c.on_count = 0u;
    c.dropped = 0u;
    const unsigned int key = bin_key_of(c.g_low, c.g_high, bin_start_bit, 2u * bin_w);
    screen_visit(&c, key, (int)bin_w, (int)max_mm, 0);
    out_dropped[gi] = c.dropped;
}

// ---- Bin-pair-tiled self-scan (all-vs-all) ----
// One block per guide bin B. Because every guide in B shares the bin-key B, the
// candidate-bin enumeration is done ONCE per bin (amortised over its ~66
// guides) instead of once per guide. The w*3 single-mutation subtrees are split
// across the block's threads; each visited candidate bin is compared against
// ALL of B's guides, accumulating per-guide counts via global atomics. Output
// is indexed by absolute entry index: out_counts[entry*stride + mm].

struct TiledCtx {
    const unsigned long long* entries;
    const unsigned int* bin_off;
    unsigned long long mask_low, mask_high;
    unsigned int max_mm, stride;
    unsigned int* out_counts;
    unsigned int g_start, g_count;
};

__device__ void tiled_process(const TiledCtx* c, unsigned int key) {
    const unsigned int o_start = c->bin_off[2u * key];
    const unsigned int o_count = c->bin_off[2u * key + 1u];
    if (o_count == 0u) return;
    const unsigned long long UPPER = 0xAAAAAAAAAAAAAAAAULL;
    for (unsigned int gi = 0; gi < c->g_count; ++gi) {
        const unsigned long long ge = (unsigned long long)(c->g_start + gi);
        const unsigned long long g_low = c->entries[3 * ge];
        const unsigned long long g_high = c->entries[3 * ge + 1];
        const unsigned int outbase = (c->g_start + gi) * c->stride;
        for (unsigned int oi = 0; oi < o_count; ++oi) {
            const unsigned long long oe = (unsigned long long)(o_start + oi);
            const unsigned long long o_low = c->entries[3 * oe];
            const unsigned long long o_high = c->entries[3 * oe + 1];
            const unsigned long long cl = (g_low ^ o_low) & c->mask_low;
            const unsigned long long ch = (g_high ^ o_high) & c->mask_high;
            const unsigned long long il = (cl & UPPER) | ((cl << 1) & UPPER);
            const unsigned long long ih = (ch & UPPER) | ((ch << 1) & UPPER);
            const unsigned int mm = __popcll(il) + __popcll(ih);
            if (mm <= c->max_mm) {
                atomicAdd(&c->out_counts[outbase + mm], 1u);
            }
        }
    }
}

__device__ void tiled_visit(
    const TiledCtx* c, unsigned int key, int w, int k_rem, int start_pos)
{
    tiled_process(c, key);
    if (k_rem == 0) return;
    for (int p = start_pos; p < w; ++p) {
        const unsigned int shift = 2u * (unsigned int)p;
        for (unsigned int delta = 1u; delta <= 3u; ++delta) {
            tiled_visit(c, key ^ (delta << shift), w, k_rem - 1, p + 1);
        }
    }
}

extern "C" __global__ void tiled_self_scan_kernel(
    const unsigned long long* __restrict__ entries,
    const unsigned int* __restrict__ bin_off,
    const unsigned int bin_lo, const unsigned int bin_hi,
    const unsigned long long mask_low, const unsigned long long mask_high,
    const unsigned int max_mm, const unsigned int bin_w,
    unsigned int* __restrict__ out_counts)
{
    const unsigned int B = bin_lo + blockIdx.x;
    if (B >= bin_hi) return;
    const unsigned int g_start = bin_off[2u * B];
    const unsigned int g_count = bin_off[2u * B + 1u];
    if (g_count == 0u) return;

    TiledCtx c;
    c.entries = entries;
    c.bin_off = bin_off;
    c.mask_low = mask_low;
    c.mask_high = mask_high;
    c.max_mm = max_mm;
    c.stride = max_mm + 1u;
    c.out_counts = out_counts;
    c.g_start = g_start;
    c.g_count = g_count;

    const int w = (int)bin_w;
    const unsigned int tid = threadIdx.x;
    // d=0 root (bin B against itself): one thread.
    if (tid == 0) tiled_process(&c, B);
    // Single-mutation subtrees, distributed across threads. Each subtree visits
    // the candidate bins whose lowest mutated position is p (delta), recursing
    // to higher positions — so every candidate bin is visited exactly once.
    if (max_mm >= 1u) {
        const unsigned int nbranch = (unsigned int)w * 3u;
        for (unsigned int branch = tid; branch < nbranch; branch += blockDim.x) {
            const unsigned int p = branch / 3u;
            const unsigned int delta = (branch % 3u) + 1u;
            const unsigned int key = B ^ (delta << (2u * p));
            tiled_visit(&c, key, w, (int)max_mm - 1, (int)p + 1);
        }
    }
}
"#;

/// Threads per block.
const BLOCK_DIM: u32 = 256;
/// Initial output-buffer capacity (hits). Doubled-and-retried on overflow.
/// 16 M hits × 3 × 4 B = 192 MB — trivial next to the resident index.
const INITIAL_OUTPUT_CAP: u32 = 16 << 20;

/// Raw kernel output: the total match count plus the (guide index, entry
/// index, mismatch) columns, truncated to whatever fit in the output buffer.
///
/// `guide` indices are relative to the `&[Guide]` slice passed to the scan;
/// `entry` indices address the resident database and decode to a `Site` +
/// `Position` via [`GpuScanner::entry`]. When `count > guide.len()` the buffer
/// overflowed and the caller should retry with a larger `out_cap`.
#[derive(Debug, Clone)]
pub struct ScanRaw {
    pub count: u32,
    pub guide: Vec<u32>,
    pub entry: Vec<u32>,
    pub mm: Vec<u32>,
}

/// Per-guide CFD aggregates from [`GpuScanner::scan_cfd_specificity`] — the
/// raw quantities needed to form a specificity score on the host, computed
/// entirely on the GPU without ever materializing a hit list.
///
/// `off_sum` = Σ CFD over off-targets (`mm > 0`); `max_cfd` = the largest
/// single off-target CFD; `on_count` = number of on-target (`mm == 0`)
/// positions; `off_count` = number of off-target (`mm > 0`) positions. The
/// host turns these into FlashFry `1/(1+off_sum)` or GuideScan
/// `1/(on_count+off_sum)` specificity.
///
/// `saturated` is set when the per-guide off-target cap was hit: the traversal
/// stopped early so the sums are **partial** (an upper bound on specificity —
/// such a guide is already specificity ≈ 0). Bounds the kernel's per-guide work
/// so a handful of homopolymer guides can't gate genome-scale wall time.
#[derive(Debug, Clone, Copy)]
pub struct CfdAgg {
    pub off_sum: f64,
    pub max_cfd: f64,
    pub on_count: u32,
    pub off_count: u32,
    pub saturated: bool,
}

/// Errors from GPU setup. `scan` itself is infallible (the `Scanner` trait
/// returns `Vec<Hit>`), so runtime CUDA failures there panic with context.
#[derive(Debug)]
pub enum GpuError {
    /// A CUDA driver call failed (no device, OOM, launch error, …).
    Driver(cudarc::driver::DriverError),
    /// The kernel failed to compile under NVRTC.
    Compile(cudarc::nvrtc::CompileError),
}

impl std::fmt::Display for GpuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Driver(e) => write!(f, "CUDA driver error: {e}"),
            Self::Compile(e) => write!(f, "NVRTC kernel compile error: {e}"),
        }
    }
}

impl std::error::Error for GpuError {}

impl From<cudarc::driver::DriverError> for GpuError {
    fn from(e: cudarc::driver::DriverError) -> Self {
        Self::Driver(e)
    }
}
impl From<cudarc::nvrtc::CompileError> for GpuError {
    fn from(e: cudarc::nvrtc::CompileError) -> Self {
        Self::Compile(e)
    }
}

/// GPU off-target scanner. Holds the database resident in VRAM; each `scan`
/// uploads the guides, runs the kernel, and copies the hits back.
pub struct GpuScanner {
    stream: Arc<CudaStream>,
    /// Hit-enumeration kernel (`scan_kernel`).
    func: CudaFunction,
    /// Per-guide brute-force count-aggregation kernel (`scan_count_kernel`).
    func_count: CudaFunction,
    /// Per-guide bin-prefilter count kernel (`scan_count_prefilter_kernel`).
    func_prefilter: CudaFunction,
    /// Per-guide bin-prefilter HIT-EMIT kernel (`scan_prefilter_emit_kernel`) —
    /// the fast genome-scale hit path used by streaming scored output.
    func_prefilter_emit: CudaFunction,
    /// Per-guide bin-prefilter CFD-accumulation kernel
    /// (`scan_cfd_prefilter_kernel`) — Mode-1 specificity with no hit download.
    func_cfd: CudaFunction,
    /// Per-guide early-exit pre-screen kernel (`scan_screen_kernel`) — flags
    /// guides with any off-target within a mismatch threshold, bailing at the
    /// first one.
    func_screen: CudaFunction,
    /// Bin-pair-tiled all-vs-all self-scan kernel (`tiled_self_scan_kernel`).
    func_tiled: CudaFunction,
    /// Database entries on the device (resident), as a flat `u64` array
    /// (3 per entry).
    d_entries: CudaSlice<u64>,
    /// Dense bin-offset table on the device: 2 `u32` per bin key
    /// (`[2*key]` = start index into entries, `[2*key+1]` = count). Empty bins
    /// have count 0. Used by the prefilter kernel to jump to candidate bins.
    d_bin_offsets: CudaSlice<u32>,
    /// The same entries on the host, to recover `Site` + `Position` for each
    /// returned hit by index.
    entries: Vec<PackedEntry>,
    n_entries: u32,
    compare_mask: SiteMask,
    /// Bin width in bases (for the prefilter's candidate-bin enumeration).
    bin_w: u32,
    /// Bit offset of the bin key within a `Site` (mirrors `MmapDb::bin_key`).
    bin_key_start_bit: u32,
    /// Resident CFD mismatch-penalty table, flattened `[g*80 + o*20 + p]`
    /// (320 f64), for `scan_cfd_prefilter_kernel`.
    d_cfd_penalty: CudaSlice<f64>,
    /// Resident CFD PAM-tail table, flattened `[b1*4 + b2]` (16 f64).
    d_cfd_pam: CudaSlice<f64>,
}

impl std::fmt::Debug for GpuScanner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuScanner")
            .field("n_entries", &self.n_entries)
            .field("compare_mask", &self.compare_mask)
            .finish_non_exhaustive()
    }
}

impl GpuScanner {
    /// Build a scanner from any [`BinSource`], uploading its entire entry set
    /// to CUDA device 0.
    ///
    /// # Errors
    /// Returns [`GpuError`] if the CUDA device/context can't be created, the
    /// kernel fails to compile, or the upload fails.
    ///
    /// # Panics
    /// Panics if the source holds more than `u32::MAX` entries (no realistic
    /// mammalian genome does).
    pub fn new(source: &dyn BinSource) -> Result<Self, GpuError> {
        Self::new_on_device(source, 0)
    }

    /// As [`new`](Self::new) but on a specific device ordinal.
    ///
    /// # Errors
    /// See [`new`](Self::new).
    ///
    /// # Panics
    /// See [`new`](Self::new).
    pub fn new_on_device(source: &dyn BinSource, ordinal: usize) -> Result<Self, GpuError> {
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let enzyme = source.enzyme();
        let compare_mask = enzyme.compare_mask;

        // Bin-key geometry (mirrors `MmapDb::bin_key`), needed by the prefilter
        // kernel to enumerate each guide's candidate bins.
        let bin_w = u32::from(source.bin_width());
        let bin_bits = bin_w * 2;
        let proto_len = u32::from(enzyme.protospacer_len);
        let pam_len = u32::from(enzyme.pam_len);
        let three_prime = matches!(enzyme.pam_side, crispr_encoding::PamSide::ThreePrime);
        let bin_top_from_lsb = if three_prime {
            (proto_len + pam_len) * 2
        } else {
            proto_len * 2
        };
        let bin_key_start_bit = bin_top_from_lsb - bin_bits;
        let bin_count = 1usize << bin_bits; // 4^bin_w

        // Flatten every bin's entries into one contiguous host array (the order
        // they'll occupy on the device) AND build a dense bin-offset table in
        // the same pass, so a guide's candidate bins map straight to entry
        // ranges. NOTE: at genome scale the dominant cost is paging the multi-GB
        // index in from storage on first touch (I/O-bound), not the copy — a
        // one-time cost amortised while the index stays VRAM-resident.
        let t_flat = std::time::Instant::now();
        let total: usize = source.iter_bins().map(|(_, e)| e.len()).sum();
        let mut entries: Vec<PackedEntry> = Vec::with_capacity(total);
        let mut bin_offsets: Vec<u32> = vec![0u32; bin_count * 2];
        for (bin_key, e) in source.iter_bins() {
            let k = bin_key as usize;
            bin_offsets[2 * k] = u32::try_from(entries.len()).expect("bin start fits u32");
            bin_offsets[2 * k + 1] = u32::try_from(e.len()).expect("bin count fits u32");
            entries.extend_from_slice(e);
        }
        let n_entries =
            u32::try_from(entries.len()).expect("entry count fits in u32 (< 4 G sites)");
        let flat_s = t_flat.elapsed().as_secs_f64();

        let ctx = CudaContext::new(ordinal)?;
        // The bin-prefilter kernels recurse over candidate bins (depth ≤ max_mm).
        // The CFD-accumulation kernel's frame is heavier than the count kernel's
        // (inlined `cfd_score_pair` with f64 accumulators), so the default
        // ~1 KB per-thread stack overruns → CUDA_ERROR_ILLEGAL_ADDRESS. Raise it
        // to a generous fixed size; a few-KB stack over the resident thread set
        // is negligible next to the multi-GB resident index.
        ctx.set_limit(
            cudarc::driver::sys::CUlimit::CU_LIMIT_STACK_SIZE,
            16 * 1024,
        )?;
        let stream = ctx.default_stream();

        // Compile for an Ampere baseline (compute_80); the driver JITs to the
        // actual SM (8.6 on A5500, 8.0 on A100) at module load.
        let t_compile = std::time::Instant::now();
        let opts = CompileOptions {
            arch: Some("compute_80"),
            ..Default::default()
        };
        let ptx = compile_ptx_with_opts(KERNEL_SRC, opts)?;
        let module = ctx.load_module(ptx)?;
        let func = module.load_function("scan_kernel")?;
        let func_count = module.load_function("scan_count_kernel")?;
        let func_prefilter = module.load_function("scan_count_prefilter_kernel")?;
        let func_prefilter_emit = module.load_function("scan_prefilter_emit_kernel")?;
        let func_cfd = module.load_function("scan_cfd_prefilter_kernel")?;
        let func_screen = module.load_function("scan_screen_kernel")?;
        let func_tiled = module.load_function("tiled_self_scan_kernel")?;
        let compile_s = t_compile.elapsed().as_secs_f64();

        // Flatten the CFD tables in the exact layout `scan_cfd_prefilter_kernel`
        // indexes (`penalty[g*80 + o*20 + p]`, `pam[b1*4 + b2]`) straight from
        // crispr-score, so the GPU scores byte-for-byte what the CPU does.
        let cfd = crispr_score::Cfd::new();
        let pen = cfd.penalty_table();
        let mut penalty_flat: Vec<f64> = Vec::with_capacity(320);
        for g in 0..4 {
            for o in 0..4 {
                for p in 0..20 {
                    penalty_flat.push(pen[g][o][p]);
                }
            }
        }
        let pam = cfd.pam_table();
        let mut pam_flat: Vec<f64> = Vec::with_capacity(16);
        for b1 in 0..4 {
            for b2 in 0..4 {
                pam_flat.push(pam[b1][b2]);
            }
        }

        // Upload the entries, reinterpreted as a flat u64 slice (PackedEntry
        // is `Pod`, 3×u64). This is the one-time cost amortised across every
        // query while the index stays VRAM-resident.
        let entries_u64: &[u64] = bytemuck::cast_slice(&entries);
        let t_up = std::time::Instant::now();
        let d_entries = stream.clone_htod(entries_u64)?;
        let d_bin_offsets = stream.clone_htod(&bin_offsets)?;
        let d_cfd_penalty = stream.clone_htod(&penalty_flat)?;
        let d_cfd_pam = stream.clone_htod(&pam_flat)?;
        stream.synchronize()?;
        let up_s = t_up.elapsed().as_secs_f64();
        if trace {
            eprintln!(
                "[trace] gpu setup: {n_entries} entries ({:.1} GB) + {bin_count} bins — \
                 flatten {flat_s:.2}s, nvrtc compile {compile_s:.2}s, H2D upload {up_s:.2}s",
                entries_u64.len() as f64 * 8.0 / 1e9,
            );
        }

        Ok(Self {
            stream,
            func,
            func_count,
            func_prefilter,
            func_prefilter_emit,
            func_cfd,
            func_screen,
            func_tiled,
            d_entries,
            d_bin_offsets,
            entries,
            n_entries,
            compare_mask,
            bin_w,
            bin_key_start_bit,
            d_cfd_penalty,
            d_cfd_pam,
        })
    }

    /// Number of database entries resident on the device.
    #[must_use]
    pub fn entry_count(&self) -> u32 {
        self.n_entries
    }

    /// Total number of bin keys (`4^bin_width`), i.e. the exclusive upper bound
    /// for `self_scan_counts`'s `bin_hi`.
    #[must_use]
    pub fn bin_count(&self) -> u32 {
        1u32 << (2 * self.bin_w)
    }

    /// Genome-scale path: instead of returning every hit, accumulate per-guide
    /// off-target **counts** by mismatch level on the GPU (`atomicAdd`),
    /// returning a compact `n_guides × (max_mismatches + 1)` row-major `u32`
    /// array indexed `out[gi * (max_mm + 1) + mm]`. This is the right output for
    /// genome-wide all-vs-all specificity over millions of guides, where a full
    /// hit list would be petabytes. The `mm == 0` bucket includes each guide's
    /// own on-target site. The index stays VRAM-resident across calls.
    ///
    /// # Panics
    /// Panics on any CUDA error (the operation is otherwise infallible).
    #[must_use]
    pub fn scan_counts(&self, guides: &[Guide], max_mismatches: u8) -> Vec<u32> {
        let stride = (u32::from(max_mismatches) + 1) as usize;
        if guides.is_empty() || self.n_entries == 0 {
            return vec![0; guides.len() * stride];
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let t = std::time::Instant::now();
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let stream = &self.stream;
        let d_glow = stream.clone_htod(&glow).expect("upload guide lows");
        let d_ghigh = stream.clone_htod(&ghigh).expect("upload guide highs");
        let mut d_counts = stream
            .alloc_zeros::<u32>(guides.len() * stride)
            .expect("alloc count buffer");

        let n_entries = self.n_entries;
        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mismatches);
        let cfg = LaunchConfig {
            grid_dim: (n_entries.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };
        {
            let mut lb = stream.launch_builder(&self.func_count);
            lb.arg(&self.d_entries);
            lb.arg(&n_entries);
            lb.arg(&d_glow);
            lb.arg(&d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&mut d_counts);
            // SAFETY: the 9 arguments match scan_count_kernel's parameter list.
            unsafe { lb.launch(cfg).expect("launch count kernel") };
        }
        stream.synchronize().expect("sync count kernel");
        let counts = stream.clone_dtoh(&d_counts).expect("download counts");
        stream.synchronize().expect("sync count download");
        if trace {
            eprintln!(
                "[trace] gpu scan_counts: {n_guides} guides, kernel+I/O {:.3}s (index resident)",
                t.elapsed().as_secs_f64()
            );
        }
        counts
    }

    /// Bin-prefilter variant of [`scan_counts`](Self::scan_counts): each guide
    /// visits only its candidate bins (within `max_mismatches` base-mismatches
    /// of its bin key), not all entries. Produces **identical** counts to
    /// `scan_counts` but with ~100×+ fewer comparisons at genome scale — the
    /// path for all-vs-all over millions of guides. Returns the same compact
    /// `n_guides × (max_mismatches + 1)` row-major `u32` array.
    ///
    /// # Panics
    /// Panics on any CUDA error (the operation is otherwise infallible).
    #[must_use]
    pub fn scan_counts_prefilter(&self, guides: &[Guide], max_mismatches: u8) -> Vec<u32> {
        let stride = (u32::from(max_mismatches) + 1) as usize;
        if guides.is_empty() || self.n_entries == 0 {
            return vec![0; guides.len() * stride];
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let t = std::time::Instant::now();
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let stream = &self.stream;
        let d_glow = stream.clone_htod(&glow).expect("upload guide lows");
        let d_ghigh = stream.clone_htod(&ghigh).expect("upload guide highs");
        let mut d_counts = stream
            .alloc_zeros::<u32>(guides.len() * stride)
            .expect("alloc count buffer");

        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mismatches);
        let bin_w = self.bin_w;
        let bin_start_bit = self.bin_key_start_bit;
        let cfg = LaunchConfig {
            grid_dim: (n_guides.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };
        {
            let mut lb = stream.launch_builder(&self.func_prefilter);
            lb.arg(&self.d_entries);
            lb.arg(&self.d_bin_offsets);
            lb.arg(&d_glow);
            lb.arg(&d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&bin_w);
            lb.arg(&bin_start_bit);
            lb.arg(&mut d_counts);
            // SAFETY: the 11 arguments match scan_count_prefilter_kernel.
            unsafe { lb.launch(cfg).expect("launch prefilter kernel") };
        }
        stream.synchronize().expect("sync prefilter kernel");
        let counts = stream.clone_dtoh(&d_counts).expect("download counts");
        stream.synchronize().expect("sync count download");
        if trace {
            eprintln!(
                "[trace] gpu scan_counts_prefilter: {n_guides} guides, kernel+I/O {:.3}s (resident)",
                t.elapsed().as_secs_f64()
            );
        }
        counts
    }

    /// Early-exit pre-screen: per guide, return `1` if it has any off-target
    /// within `threshold` mismatches (a 0-mm duplicate — `>= 2` zero-mismatch
    /// hits — or any `1..=threshold`-mm match), else `0`. Bails at the first
    /// such hit, so non-specific guides are far cheaper than counting all of
    /// them. No hit download; the index stays VRAM-resident across calls.
    ///
    /// # Panics
    /// Panics on any CUDA error.
    #[must_use]
    pub fn scan_screen(&self, guides: &[Guide], threshold: u8) -> Vec<u32> {
        if guides.is_empty() || self.n_entries == 0 {
            return vec![0; guides.len()];
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let t = std::time::Instant::now();
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let stream = &self.stream;
        let d_glow = stream.clone_htod(&glow).expect("upload guide lows");
        let d_ghigh = stream.clone_htod(&ghigh).expect("upload guide highs");
        let mut d_drop = stream
            .alloc_zeros::<u32>(guides.len())
            .expect("alloc screen buffer");

        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(threshold);
        let bin_w = self.bin_w;
        let bin_start_bit = self.bin_key_start_bit;
        let cfg = LaunchConfig {
            grid_dim: (n_guides.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };
        {
            let mut lb = stream.launch_builder(&self.func_screen);
            lb.arg(&self.d_entries);
            lb.arg(&self.d_bin_offsets);
            lb.arg(&d_glow);
            lb.arg(&d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&bin_w);
            lb.arg(&bin_start_bit);
            lb.arg(&mut d_drop);
            // SAFETY: the 11 arguments match scan_screen_kernel's parameter list.
            unsafe { lb.launch(cfg).expect("launch screen kernel") };
        }
        stream.synchronize().expect("sync screen kernel");
        let dropped = stream.clone_dtoh(&d_drop).expect("download screen flags");
        stream.synchronize().expect("sync screen download");
        if trace {
            eprintln!(
                "[trace] gpu scan_screen: {n_guides} guides, kernel+I/O {:.3}s (resident)",
                t.elapsed().as_secs_f64()
            );
        }
        dropped
    }

    /// Bin-pair-tiled all-vs-all self-scan: treats every index entry in bins
    /// `[bin_lo, bin_hi)` as a guide and counts its off-targets by mismatch
    /// level, amortising the candidate-bin enumeration over each bin's guides
    /// (one CUDA block per bin). Returns a compact `n_entries × (max_mm + 1)`
    /// row-major `u32` array indexed `out[entry * stride + mm]`; entries outside
    /// the bin range are left zero. Pass `bin_lo = 0`, `bin_hi = bin_count()`
    /// for the full genome self-scan, or a sub-range to benchmark a slice.
    ///
    /// # Panics
    /// Panics on any CUDA error.
    #[must_use]
    pub fn self_scan_counts(&self, bin_lo: u32, bin_hi: u32, max_mismatches: u8) -> Vec<u32> {
        let stride = (u32::from(max_mismatches) + 1) as usize;
        let n_out = self.n_entries as usize * stride;
        if n_out == 0 || bin_hi <= bin_lo {
            return vec![0; n_out];
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let stream = &self.stream;
        let mut d_counts = stream.alloc_zeros::<u32>(n_out).expect("alloc count buffer");

        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mismatches);
        let bin_w = self.bin_w;
        let cfg = LaunchConfig {
            grid_dim: (bin_hi - bin_lo, 1, 1),
            block_dim: (64, 1, 1),
            shared_mem_bytes: 0,
        };
        let t_k = std::time::Instant::now();
        {
            let mut lb = stream.launch_builder(&self.func_tiled);
            lb.arg(&self.d_entries);
            lb.arg(&self.d_bin_offsets);
            lb.arg(&bin_lo);
            lb.arg(&bin_hi);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&bin_w);
            lb.arg(&mut d_counts);
            // SAFETY: the 9 arguments match tiled_self_scan_kernel.
            unsafe { lb.launch(cfg).expect("launch tiled kernel") };
        }
        stream.synchronize().expect("sync tiled kernel");
        let kernel_s = t_k.elapsed().as_secs_f64();
        let counts = stream.clone_dtoh(&d_counts).expect("download counts");
        stream.synchronize().expect("sync count download");
        if trace {
            eprintln!(
                "[trace] gpu self_scan tiled: bins [{bin_lo},{bin_hi}), kernel {kernel_s:.3}s + download"
            );
        }
        counts
    }

    /// Run the kernel once with a given output capacity, returning a
    /// [`ScanRaw`]. If `count > out_cap` the columns are truncated to
    /// `out_cap` and the caller should retry with a larger cap.
    fn run_once(
        &self,
        d_glow: &CudaSlice<u64>,
        d_ghigh: &CudaSlice<u64>,
        n_guides: u32,
        max_mm: u8,
        out_cap: u32,
    ) -> Result<ScanRaw, GpuError> {
        let stream = &self.stream;
        let mut d_count = stream.alloc_zeros::<u32>(1)?;
        let mut d_guide = stream.alloc_zeros::<u32>(out_cap as usize)?;
        let mut d_entry = stream.alloc_zeros::<u32>(out_cap as usize)?;
        let mut d_mm = stream.alloc_zeros::<u32>(out_cap as usize)?;

        let n_entries = self.n_entries;
        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mm);
        let cfg = LaunchConfig {
            grid_dim: (n_entries.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };

        {
            let mut lb = stream.launch_builder(&self.func);
            lb.arg(&self.d_entries);
            lb.arg(&n_entries);
            lb.arg(d_glow);
            lb.arg(d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&mut d_count);
            lb.arg(&out_cap);
            lb.arg(&mut d_guide);
            lb.arg(&mut d_entry);
            lb.arg(&mut d_mm);
            // SAFETY: the 13 arguments match scan_kernel's parameter list in
            // count, type, and order.
            unsafe { lb.launch(cfg)? };
        }
        stream.synchronize()?;

        let count_v = stream.clone_dtoh(&d_count)?;
        let mut guide = stream.clone_dtoh(&d_guide)?;
        let mut entry = stream.clone_dtoh(&d_entry)?;
        let mut mm = stream.clone_dtoh(&d_mm)?;
        stream.synchronize()?;

        let count = count_v[0];
        let kept = count.min(out_cap) as usize;
        guide.truncate(kept);
        entry.truncate(kept);
        mm.truncate(kept);
        Ok(ScanRaw {
            count,
            guide,
            entry,
            mm,
        })
    }

    /// Upload `guides`, run a single bounded hit-emission pass against the
    /// resident index, and return the raw `(guide, entry, mismatch)` columns.
    ///
    /// `out_cap` bounds the device output buffer: if the true hit count
    /// exceeds it, [`ScanRaw::count`] reports the full count but the columns
    /// are truncated to `out_cap`. Size it via [`Self::scan_counts_prefilter`]
    /// (sum the per-guide counts) to emit every hit in one pass with no retry.
    ///
    /// This is the building block for streaming/batched scoring: hold the
    /// index resident, call `scan_raw` on each guide batch, decode + process
    /// its hits, then move on — peak host memory stays `O(batch hits)`.
    ///
    /// # Errors
    /// Returns [`GpuError`] on any CUDA failure (upload, launch, download).
    pub fn scan_raw(
        &self,
        guides: &[Guide],
        max_mismatches: u8,
        out_cap: u32,
    ) -> Result<ScanRaw, GpuError> {
        if guides.is_empty() || self.n_entries == 0 {
            return Ok(ScanRaw {
                count: 0,
                guide: Vec::new(),
                entry: Vec::new(),
                mm: Vec::new(),
            });
        }
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let d_glow = self.stream.clone_htod(&glow)?;
        let d_ghigh = self.stream.clone_htod(&ghigh)?;
        self.run_once(&d_glow, &d_ghigh, n_guides, max_mismatches, out_cap)
    }

    /// Decode a database entry index (as returned in [`ScanRaw::entry`]) into
    /// its packed `Site` + `Position` from the resident host copy. Cheap — a
    /// single `Vec` index, no device round-trip.
    ///
    /// # Panics
    /// Panics if `entry_index` is out of range (`>= entry_count()`).
    #[must_use]
    pub fn entry(&self, entry_index: u32) -> PackedEntry {
        self.entries[entry_index as usize]
    }

    /// Like [`Self::scan_raw`] but emits hits via the bin-prefilter kernel —
    /// each guide visits only its candidate bins, ~40× faster than the
    /// brute-force path at genome scale, with a byte-identical hit set. Size
    /// `out_cap` from [`Self::scan_counts_prefilter`] (sum the per-guide
    /// counts) to emit every hit in one pass with no retry.
    ///
    /// # Errors
    /// Returns [`GpuError`] on any CUDA failure (upload, launch, download).
    #[allow(clippy::cast_possible_truncation)]
    pub fn scan_prefilter_raw(
        &self,
        guides: &[Guide],
        max_mismatches: u8,
        out_cap: u32,
        cap_per_guide: u32,
    ) -> Result<ScanRaw, GpuError> {
        if guides.is_empty() || self.n_entries == 0 {
            return Ok(ScanRaw {
                count: 0,
                guide: Vec::new(),
                entry: Vec::new(),
                mm: Vec::new(),
            });
        }
        let stream = &self.stream;
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let d_glow = stream.clone_htod(&glow)?;
        let d_ghigh = stream.clone_htod(&ghigh)?;
        let mut d_count = stream.alloc_zeros::<u32>(1)?;
        let mut d_guide = stream.alloc_zeros::<u32>(out_cap as usize)?;
        let mut d_entry = stream.alloc_zeros::<u32>(out_cap as usize)?;
        let mut d_mm = stream.alloc_zeros::<u32>(out_cap as usize)?;

        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mismatches);
        let bin_w = self.bin_w;
        let bin_start_bit = self.bin_key_start_bit;
        let cfg = LaunchConfig {
            grid_dim: (n_guides.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };
        {
            let mut lb = stream.launch_builder(&self.func_prefilter_emit);
            lb.arg(&self.d_entries);
            lb.arg(&self.d_bin_offsets);
            lb.arg(&d_glow);
            lb.arg(&d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&bin_w);
            lb.arg(&bin_start_bit);
            lb.arg(&mut d_count);
            lb.arg(&out_cap);
            lb.arg(&mut d_guide);
            lb.arg(&mut d_entry);
            lb.arg(&mut d_mm);
            lb.arg(&cap_per_guide);
            // SAFETY: the 16 arguments match scan_prefilter_emit_kernel's
            // parameter list in count, type, and order.
            unsafe { lb.launch(cfg)? };
        }
        stream.synchronize()?;

        let count_v = stream.clone_dtoh(&d_count)?;
        let mut guide = stream.clone_dtoh(&d_guide)?;
        let mut entry = stream.clone_dtoh(&d_entry)?;
        let mut mm = stream.clone_dtoh(&d_mm)?;
        stream.synchronize()?;

        let count = count_v[0];
        let kept = count.min(out_cap) as usize;
        guide.truncate(kept);
        entry.truncate(kept);
        mm.truncate(kept);
        Ok(ScanRaw {
            count,
            guide,
            entry,
            mm,
        })
    }

    /// Mode-1 specificity path: for each guide, compute its off-target CFD
    /// **sum**, **max**, and on-target multiplicity entirely on the GPU
    /// (bin-prefilter traversal + on-the-fly `Cfd::score_pair`, register
    /// accumulation) and return one [`CfdAgg`] per guide. **No hit list is
    /// ever materialized or downloaded**, so memory and bandwidth are
    /// `O(#guides)` regardless of off-target volume — the only path that
    /// scales SpCas9 Mode-1 scoring to a full genome of repeat-rich guides.
    ///
    /// SpCas9-only: CFD is defined for 20-nt protospacer + NGG. For Cas12a or
    /// per-off-target output, use the [`scan_prefilter_raw`](Self::scan_prefilter_raw)
    /// emit path instead.
    ///
    /// # Errors
    /// Returns [`GpuError`] on any CUDA failure (upload, launch, download).
    pub fn scan_cfd_specificity(
        &self,
        guides: &[Guide],
        max_mismatches: u8,
        ot_cap: u32,
        off_sum_cap: f64,
    ) -> Result<Vec<CfdAgg>, GpuError> {
        if guides.is_empty() || self.n_entries == 0 {
            return Ok(vec![
                CfdAgg {
                    off_sum: 0.0,
                    max_cfd: 0.0,
                    on_count: 0,
                    off_count: 0,
                    saturated: false,
                };
                guides.len()
            ]);
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let t = std::time::Instant::now();
        let stream = &self.stream;
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let d_glow = stream.clone_htod(&glow)?;
        let d_ghigh = stream.clone_htod(&ghigh)?;
        let mut d_sum = stream.alloc_zeros::<f64>(guides.len())?;
        let mut d_max = stream.alloc_zeros::<f64>(guides.len())?;
        let mut d_on = stream.alloc_zeros::<u32>(guides.len())?;
        let mut d_off = stream.alloc_zeros::<u32>(guides.len())?;
        let mut d_sat = stream.alloc_zeros::<u32>(guides.len())?;

        let mask_low = self.compare_mask.low;
        let mask_high = self.compare_mask.high;
        let max_mm_u32 = u32::from(max_mismatches);
        let bin_w = self.bin_w;
        let bin_start_bit = self.bin_key_start_bit;
        let cfg = LaunchConfig {
            grid_dim: (n_guides.div_ceil(BLOCK_DIM), 1, 1),
            block_dim: (BLOCK_DIM, 1, 1),
            shared_mem_bytes: 0,
        };
        {
            let mut lb = stream.launch_builder(&self.func_cfd);
            lb.arg(&self.d_entries);
            lb.arg(&self.d_bin_offsets);
            lb.arg(&d_glow);
            lb.arg(&d_ghigh);
            lb.arg(&n_guides);
            lb.arg(&mask_low);
            lb.arg(&mask_high);
            lb.arg(&max_mm_u32);
            lb.arg(&bin_w);
            lb.arg(&bin_start_bit);
            lb.arg(&self.d_cfd_penalty);
            lb.arg(&self.d_cfd_pam);
            lb.arg(&ot_cap);
            lb.arg(&off_sum_cap);
            lb.arg(&mut d_sum);
            lb.arg(&mut d_max);
            lb.arg(&mut d_on);
            lb.arg(&mut d_off);
            lb.arg(&mut d_sat);
            // SAFETY: the 19 arguments match scan_cfd_prefilter_kernel's
            // parameter list in count, type, and order.
            unsafe { lb.launch(cfg)? };
        }
        stream.synchronize()?;
        let sum = stream.clone_dtoh(&d_sum)?;
        let max = stream.clone_dtoh(&d_max)?;
        let on = stream.clone_dtoh(&d_on)?;
        let off = stream.clone_dtoh(&d_off)?;
        let sat = stream.clone_dtoh(&d_sat)?;
        stream.synchronize()?;
        if trace {
            eprintln!(
                "[trace] gpu scan_cfd_specificity: {n_guides} guides, kernel+I/O {:.3}s (resident)",
                t.elapsed().as_secs_f64()
            );
        }
        Ok((0..guides.len())
            .map(|i| CfdAgg {
                off_sum: sum[i],
                max_cfd: max[i],
                on_count: on[i],
                off_count: off[i],
                saturated: sat[i] != 0,
            })
            .collect())
    }
}

impl Scanner for GpuScanner {
    fn scan(&self, guides: &[Guide], max_mismatches: u8) -> Vec<Hit> {
        if guides.is_empty() || self.n_entries == 0 {
            return Vec::new();
        }
        let trace = std::env::var_os("CRISPR_OTS_TRACE").is_some();
        let t_scan = std::time::Instant::now();
        let n_guides = u32::try_from(guides.len()).expect("guide count fits in u32");
        let glow: Vec<u64> = guides.iter().map(|g| g.site.low).collect();
        let ghigh: Vec<u64> = guides.iter().map(|g| g.site.high).collect();
        let d_glow = self.stream.clone_htod(&glow).expect("upload guide lows");
        let d_ghigh = self.stream.clone_htod(&ghigh).expect("upload guide highs");

        // Adaptive cap: re-run with a larger buffer if the true match count
        // exceeded it (a homopolymer guide alone can have ~10^6 off-targets).
        // The count is deterministic, so this converges in at most two passes.
        let mut cap = INITIAL_OUTPUT_CAP;
        loop {
            let raw = self
                .run_once(&d_glow, &d_ghigh, n_guides, max_mismatches, cap)
                .expect("GPU scan kernel");
            if raw.count <= cap {
                let mut hits = Vec::with_capacity(raw.guide.len());
                for ((&gi, &ei), &m) in raw.guide.iter().zip(&raw.entry).zip(&raw.mm) {
                    let e = self.entries[ei as usize];
                    hits.push(Hit {
                        guide_index: gi,
                        off_target: e.site(),
                        position: e.position(),
                        mismatches: u8::try_from(m).expect("mismatch fits in u8"),
                    });
                }
                if trace {
                    eprintln!(
                        "[trace] gpu scan: {n_guides} guides, {} hits, kernel+I/O {:.3}s \
                         (index resident, no re-upload)",
                        hits.len(),
                        t_scan.elapsed().as_secs_f64()
                    );
                }
                return hits;
            }
            cap = raw
                .count
                .checked_next_power_of_two()
                .expect("output cap fits in u32");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crispr_db::{read_fasta_from, BinTable, SiteFinder};
    use crispr_encoding::{Enzyme, Site};
    use crispr_scan::BinScanner;

    fn build_table(fasta: &[u8], enzyme: Enzyme) -> BinTable {
        let contigs = read_fasta_from(fasta).expect("parses FASTA");
        let finder = SiteFinder::new(enzyme.clone());
        let mut table = BinTable::for_enzyme(enzyme);
        for (i, contig) in contigs.iter().enumerate() {
            finder.scan(
                &contig.sequence,
                u32::try_from(i).expect("contig fits u32"),
                &mut table,
            );
        }
        table
    }

    /// The GPU kernel must find exactly the same hit set as the CPU
    /// `BinScanner` on the same database and guides, at every mismatch
    /// budget. Skips on hosts without a CUDA device (the login node and the
    /// CPU-only compute nodes), where the dispatch would never select GPU.
    #[test]
    #[allow(clippy::similar_names)]
    fn gpu_matches_cpu_when_device_present() {
        if CudaContext::new(0).is_err() {
            return;
        }
        // On-target plus off-targets at 1, 3 mismatches, and a very different
        // guide — exercises the popcount across the masked protospacer.
        let fasta = b">c\n\
            AAAAAAAAAAAAAAAAAAAAAGG\
            CCCCCCC\
            AAAAAAAAAAAAAAAAAAATAGG\
            CCCCCCC\
            AAAAAAAAAAAAAAAACAATAGG\
            CCCCCCC\
            TTTTTTTTTTTTTTTTTTTTAGG\n";
        let table = build_table(fasta, Enzyme::spcas9_ngg());
        let guides = vec![
            Guide {
                id: "g0".into(),
                site: Site::encode_ascii(b"AAAAAAAAAAAAAAAAAAAAAGG"),
            },
            Guide {
                id: "g1".into(),
                site: Site::encode_ascii(b"TTTTTTTTTTTTTTTTTTTTAGG"),
            },
        ];
        let gpu = GpuScanner::new(&table).expect("build GPU scanner");
        let cpu = BinScanner::new(&table);
        let key = |h: &Hit| {
            (
                h.guide_index,
                h.position.contig_id,
                h.position.offset,
                h.position.strand.as_char(),
                h.mismatches,
                h.off_target.low,
                h.off_target.high,
            )
        };
        for mm in [0u8, 1, 2, 3, 4] {
            let mut hits_gpu = gpu.scan(&guides, mm);
            let mut hits_cpu = cpu.scan(&guides, mm);
            hits_gpu.sort_by_key(key);
            hits_cpu.sort_by_key(key);
            let keys_gpu: Vec<_> = hits_gpu.iter().map(key).collect();
            let keys_cpu: Vec<_> = hits_cpu.iter().map(key).collect();
            assert_eq!(keys_gpu, keys_cpu, "GPU and CPU disagree at mm={mm}");
        }
    }
}
