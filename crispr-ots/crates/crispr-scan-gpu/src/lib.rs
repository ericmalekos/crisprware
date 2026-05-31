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
struct ScanRaw {
    count: u32,
    guide: Vec<u32>,
    entry: Vec<u32>,
    mm: Vec<u32>,
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
        let func_tiled = module.load_function("tiled_self_scan_kernel")?;
        let compile_s = t_compile.elapsed().as_secs_f64();

        // Upload the entries, reinterpreted as a flat u64 slice (PackedEntry
        // is `Pod`, 3×u64). This is the one-time cost amortised across every
        // query while the index stays VRAM-resident.
        let entries_u64: &[u64] = bytemuck::cast_slice(&entries);
        let t_up = std::time::Instant::now();
        let d_entries = stream.clone_htod(entries_u64)?;
        let d_bin_offsets = stream.clone_htod(&bin_offsets)?;
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
            func_tiled,
            d_entries,
            d_bin_offsets,
            entries,
            n_entries,
            compare_mask,
            bin_w,
            bin_key_start_bit,
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
