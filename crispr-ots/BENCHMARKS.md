# Benchmarks

Living comparison of `crispr-ots` against FlashFry and GuideScan2 on the
same workload. Updated whenever a substantive change lands in the engine.

## Workload

- **Reference:** human chromosome 22 (`chr22.fa.gz`, ~50 Mbp, ~5 M Cas9 NGG
  sites)
- **Query:** 1 000 random chr22 SpCas9 NGG protospacers + PAM, generated
  with seed 42 from the forward-strand NGG matches in chr22 (deterministic
  helper in `NOTES.md`)
- **Enzyme:** SpCas9 NGG, 20-nt protospacer + 3-nt PAM
- **Max mismatches:** 4
- **Scoring:** Doench 2016 CFD

## Methodology

All three tools run on the same machine. Wall and peak RSS are measured
via `/usr/bin/time -v`. Disk size is `du -b` of every file the tool needs
at query time. Each tool is allowed to use its default thread count (so
the comparison reflects out-of-the-box performance, not a forced
single-threaded budget). Discover-phase numbers are warm-cache; the
preceding index/build invocation has populated the page cache.

GuideScan2's PAM matching defaults to *exact* per-guide-row PAM. To get
NGG-equivalent semantics for the comparison, we pass
`--alt-pam AGG --alt-pam CGG --alt-pam GGG --alt-pam TGG`. See
`NOTES.md` for the convention difference and the spec-formula reconciliation.

## Hardware

- Intel Core i7-10870H @ 2.20 GHz, 16 logical cores
- 62 GiB RAM, NVMe SSD
- Linux 6.17.9

## Snapshot: commit `7c590d8` (Phase 4a + bin-prefix prefilter, single-thread scalar)

### Index / build phase (one-time, fresh state)

| Tool | Wall time | Peak RSS | On-disk size |
|---|---:|---:|---:|
| **FlashFry** `index` | 32.5 s | 1 623 MB | 42 MB |
| **GuideScan2** `index` | 14.6 s | 289 MB | 53 MB (query-needed)¹ |
| **crispr-ots** `build` | **3.9 s** | **268 MB** | 139 MB |

¹ GuideScan2 writes 154 MB during build (`chr22.fa.forward.dna` +
`chr22.fa.reverse.dna` + the `.gs2.*` FM-indexes), but only the 53 MB of
`.gs2.*` files are needed for `enumerate`. The `.dna` files can be deleted
after build.

### Discover / enumerate phase (1 000 chr22 NGG guides + CFD, warm cache)

| Tool | Wall time | Peak RSS | Threads used |
|---|---:|---:|---:|
| **FlashFry** `discover` + `score` | 6.23 s | 635 MB | 1–2² |
| **GuideScan2** `enumerate` (default) | 8.50 s | 319 MB | up to 16 (default `-n`) |
| **crispr-ots** `discover --score cfd` | **3.79 s** | **241 MB** | **1** |

² FlashFry's JVM uses background threads for JIT/GC but the scan itself
is single-threaded.

### Correctness against both reference tools

- vs. **FlashFry** (`cli_cfd_chr22_1000` integration test): 930 / 930
  scoreable guides match exactly. `otCount`: 0 mismatches. CFD
  specificity: max |Δ| = 1.11e-16 (one ULP of float64).
- vs. **GuideScan2** (`cli_guidescan2_equivalence` integration test):
  1 000 / 1 000 guides match. Off-target position count: 0 mismatches.
  CFD specificity reconciles via
  `spec_gs2 = 1 / (1/spec_ours + (n_on_target - 1))` to max |Δ| = 6e-7
  (GuideScan2's 6-significant-digit print precision floor).

### Headline takeaways

- **Index/build:** crispr-ots is 8.3× faster than FlashFry's index step
  and 3.7× faster than GuideScan2's, with 6× less peak RAM than FlashFry.
- **Live discover:** crispr-ots is fastest in wall time *and* lowest in
  peak RSS, despite being single-threaded against GuideScan2's
  multi-threaded default and FlashFry's optimized bin-scan.
- **Disk:** crispr-ots is the **worst** here at 139 MB (bincode is not
  space-optimal). FlashFry's BGZF-compressed format is 3.3× smaller,
  GuideScan2's FM-index is 2.6× smaller. Phase 4b (custom mmap-friendly
  format + optional zstd) is the cure.

## 10 000-guide workload (Phase 6 single-thread regime)

`random_10k.fasta`: 10 000 random chr22 NGG sites (seed 42; the 1 000-guide
fixture × 10). At single-thread the scan dominates (~85 % of wall),
which is the right regime to see SIMD impact clearly. Built with
`.cargo/config.toml` setting `target-cpu=x86-64-v3` so `wide`'s u64x4
lowers to AVX2 + BMI2 + POPCNT.

| Tool | Threads | Wall | Peak RSS |
|---|---:|---:|---:|
| FlashFry discover + score (10 K) | 1-2 | 40.6 s | ~1.6 GB |
| **crispr-ots** (w=9, AVX2, 1 thread) | 1 | 16.3 s | 1.28 GB |
| **crispr-ots** (w=9, AVX2, 16 threads) | 16 | **4.68 s** | 2.17 GB |

Phase 6 single-thread breakdown vs the Phase 5 scalar baseline:

| Phase | Scan time (10 K, w=9, 1 thread) | Total wall |
|---|---:|---:|
| Phase 5 (scalar `mismatches_masked`) | 16.96 s | 19.10 s |
| Phase 6 (`wide::u64x4` + SWAR popcount + AVX2) | 13.70 s | 16.33 s |

The ~24 % scan speedup is real but smaller than the ~4× a naive u64-lane
count would predict. The reason is that this benchmark machine
(i7-10870H, Comet Lake) has AVX2 but **not** AVX-512 VPOPCNTDQ, so each
SIMD lane has to do a 10-op SWAR popcount instead of a single hardware
instruction. Targets that do have VPOPCNTDQ (Ice Lake+, Tiger Lake+,
Zen 4+, Sapphire Rapids) should see closer to the full theoretical
speedup — adding a `target-feature=+avx512vpopcntdq` code path is a
clean future follow-on.

## Thread sweep (Phase 5: rayon over the bin loop)

Workload as above; same chr22 DB; warm cache. `RAYON_NUM_THREADS` set
explicitly. Total wall includes DB load + FASTA read + scan + TSV write.

| Threads | width 7 wall | width 9 wall |
|---:|---:|---:|
| 1 | 4.03 s | 2.20 s |
| 2 | 2.21 s | 1.39 s |
| 4 | 1.34 s | 1.00 s |
| 8 | 0.91 s | 0.83 s |
| 16 | 0.92 s | **0.76 s** |

Scan dominates at 1-2 threads; fixed overheads (DB-load, FASTA-read,
TSV-write) dominate at 8+. Width 9 + 16 threads is **8.2× faster than
FlashFry** (6.23 s) and **11.2× faster than GuideScan2** (8.50 s) on
this workload, with comparable or lower peak RSS.

Correctness preserved at every thread count: max |Δ spec| vs FlashFry
= 2.22e-16 across the sweep.

## Headline comparison (Phase 5, width 9, 16 threads)

| Tool | Discover wall | Peak RSS | Threads |
|---|---:|---:|---:|
| FlashFry discover + score | 6.23 s | 635 MB | 1-2 |
| GuideScan2 enumerate | 8.50 s | 319 MB | 16 (default) |
| **crispr-ots** discover (w=9, 16t) | **0.76 s** | 406 MB | 16 |

## Phase 4b: mmap-friendly v2 format

The on-disk format is now a hand-rolled header + bin-offsets table +
contiguous `PackedEntry` array (24 bytes / entry), all `repr(C) Pod`.
`crispr-ots discover` opens it via `memmap2` and casts the entries
region directly to `&[PackedEntry]` — zero deserialization. The same
trait surface (`BinSource`) is implemented by both the in-memory
`BinTable` (build path) and the mmap'd `MmapDb` (query path), so the
scanner is generic over either.

Chr22, width 9, single-thread, 10 K guides + CFD:

| Format | DB size | Build wall | Discover wall (warm) |
|---|---:|---:|---:|
| Phase 4a (bincode) | 139 MB | 3.91 s | 16.33 s |
| Phase 4b (mmap, v2) | **116 MB** | 5.02 s | 16.57 s |

**Why the disk shrank**: 24-byte `PackedEntry` records replace
bincode's ~28 byte (Site, Position) encoding (bincode pays for variant
tags on the `Strand` enum and for length prefixes on each Vec).

**Why the warm-cache discover wall is unchanged**: the OS keeps the
bincode payload in the page cache between runs, so the
"deserialization" cost approaches a memcpy in both cases. The mmap
advantage is on **first-run / cold-cache / huge-DB** scenarios — for
example, a future full-human-genome database (~6 GB) opens in
constant time via mmap, but would take tens of seconds to bincode-
deserialize into RAM. Empirical confirmation is queued for whenever
we have a human FASTA on hand.

Other wins from the refactor that the table above doesn't show:
- The 116 MB file is **completely self-describing**: header carries
  magic/version, enzyme metadata, contig table, and per-bin offsets.
  A corrupt or wrong-tool file fails at `MmapDb::open` with a clean
  error instead of a bincode decode panic.
- The mmap path uses ~zero heap. Process RSS during discover is
  dominated by hit accumulators + per-thread scratch, not by the
  database itself (the OS pages the mmap region in on demand).

## Full mouse genome (GRCm39 primary assembly)

The chr22 benchmark settles the question of correctness; it doesn't
settle the question of scale. Mouse GRCm39 is **~2.7 Gbp** and ~277 M
SpCas9 NGG sites (≈55× chr22), which is roughly the size of any real
target genome (human is ~3 Gbp). This is where the integration-time
indexes will be built.

### Setup
- **Reference:** `GRCm39.primary_assembly.genome.fa.gz` (739 MB compressed,
  ~2.6 GB decompressed)
- **Query:** 1 000 random forward-strand NGG sites sampled with seed 42
  across all contigs (`mm39-bench/make_guides.py`)
- **Enzyme:** SpCas9 NGG; max mismatches 4
- **crispr-ots bin width:** 11 (~4× chr22's width-9 sweet spot)
- **FlashFry parallelism:** 8 parallel JVM processes, query FASTA split
  evenly. FlashFry's discover is single-threaded internally; chunking
  across processes mirrors how a production pipeline would scale it.

### Index / build phase (one-time)

| Tool | Wall time | Peak RSS | On-disk size |
|---|---:|---:|---:|
| **FlashFry** `index` | 38 m 38 s | 10.3 GB | **2.1 GB** |
| **GuideScan2** `index` | 31 m 19 s | 25.4 GB | **2.6 GB** |
| **crispr-ots** `index` (`--bin-width 11`) | **6 m 26 s** | 12.4 GB | 6.3 GB |

crispr-ots's index step is **4.9× faster than GuideScan2** and
**6.0× faster than FlashFry** at mouse scale, and uses **half the peak
RSS of GuideScan2's FM-index construction**. On-disk size is the
remaining weak spot for the `.crot` format — the 24-byte `PackedEntry`
representation is uncompressed; a future zstd-compressed cold-storage
variant would close the gap.

### Enumerate phase (1 000 NGG guides + CFD, warm cache)

| Tool | Wall time | Peak RSS | Threads |
|---|---:|---:|---:|
| **GuideScan2** `enumerate` | 1 m 16 s | 7.2 GB | 16 |
| **FlashFry** chunked `discover + score` | 37.5 s | 1.1 GB / proc × 8 ≈ 8.6 GB | 8 procs |
| **crispr-ots** `enumerate` (cold, 1 thread) | 18.6 s | 8.0 GB | 1 |
| **crispr-ots** `enumerate` (warm, 16 threads) | **5.9 s** | 8.5 GB | 16 |

### Mouse-scale speedups

| Comparison | Speedup |
|---|---:|
| crispr-ots 16-thread vs GuideScan2 16-thread | **12.9× faster** |
| crispr-ots 16-thread vs FlashFry chunked (8 procs) | **6.4× faster** |
| crispr-ots 1-thread vs FlashFry chunked (8 procs) | 2.0× faster |
| crispr-ots index vs GuideScan2 index | **4.9× faster** |
| crispr-ots index vs FlashFry index | **6.0× faster** |

### Bin-width sweep at mouse scale

The chr22 sweep settled width 9 for ~5 M sites. Mouse has ~138 M NGG
sites (≈28×), so the prediction was width ≈ 11 (`log4(28) ≈ 2.4`
above chr22's 9). Empirical verification on the same 1 000-guide
workload, warm-cache, all measurements with `/usr/bin/time -v`:

| Width | Build wall | t=1 enum wall | t=16 enum wall | `.crot` size |
|---:|---:|---:|---:|---:|
| 9 | 5 m 08 s | 44.96 s | 8.94 s | 6.62 GB |
| 10 | 5 m 35 s | 23.89 s | 6.55 s | 6.63 GB |
| **11** | 6 m 26 s* | **18.79 s** | **6.00 s** | 6.69 GB |
| 12 | 6 m 11 s | 29.04 s | 7.59 s | 6.79 GB |

*Width 11's build wall is from the headline build; the sweep reused
that index so didn't re-time it.

**Width 11 is the empirical optimum on both axes** — single-thread
*and* 16-thread. Width 12 regresses because the per-bin fixed cost
(prefix-match check, prefetch, function-call boilerplate) starts to
dominate when each bin holds only ~8 entries (138 M / 4¹² = 8.2/bin)
versus ~33 entries at width 11 (138 M / 4¹¹). Width 9 is 7.5× slower
at one thread because the prefilter passes far too many candidates.

Disk and peak RSS are roughly flat across the sweep (`.crot` grows
1 % from width 9 to 12; RSS within 100 MB), so the choice is a
pure runtime trade-off. Driver: `mm39-bench/binwidth_sweep.sh`,
raw numbers in `mm39-bench/binwidth_sweep.tsv`.

### Correctness against both reference tools (mouse)

End-to-end correctness on the same 1 000-guide workload, scored against
both reference outputs via `mm39-bench/compare_specs.py`:

- **vs GuideScan2** (per-guide CFD specificity under
  `--spec-convention guidescan`): 1 000 / 1 000 common guides, max |Δ|
  = **1.0e-6** — GuideScan2's 6-significant-figure print precision floor.
  Bit-identical when both tools use the same formula.
- **vs FlashFry** (per-guide CFD specificity under
  `--spec-convention flashfry`): 830 / 830 common guides (FlashFry
  drops the 224 OVERFLOW guides at the score step), max |Δ| =
  **1.11e-16** — one ULP of `f64`. Bit-identical.

Both equivalences run automatically in `cli_guidescan2_equivalence`
(chr22 fixture) and `mm39-bench/compare_specs.py` (mouse), so the
guarantee survives further engine changes.

## Cas12a (Cpf1, TTTN) on mouse GRCm39

Same fixture protocol as SpCas9 but with `--enzyme cpf1`. The protospacer
is now **23-nt** (matching the 23-position published Cas12a activity
matrices; we diverge from FlashFry's 24-bp encoding cap). Total scan
length 27 bp — inside our 28-base `MAX_SITE_LEN` envelope. GuideScan2
doesn't ship Cas12a support, so the comparison is crispr-ots vs FlashFry
only and is limited to historical numbers from the 20-nt-protospacer run
(FlashFry can't go past 20 nt without bumping its bit-encoding scheme).

### Index / build phase

| Tool | Wall time | Peak RSS | On-disk size |
|---|---:|---:|---:|
| **FlashFry** `index --enzyme cpf1` (20 nt) | 23 m 09 s | 11.8 GB | **1.3 GB** |
| **crispr-ots** `index --enzyme cpf1 --bin-width 11` (20 nt) | 4 m 25 s | 8.9 GB | 4.0 GB |
| **crispr-ots** `index --enzyme cpf1 --bin-width 11` (23 nt) | **4 m 29 s** | 8.9 GB | 4.0 GB |

The 23-nt build is essentially free vs 20-nt — same site count, same
disk, +4 s scan time. ~5.2× faster index than FlashFry's 20-nt path.

### Discover / enumerate phase (1 000 TTTN guides + mm ≤ 4)

Mouse has ~88.5 M TTTN sites. Side-by-side with the 20-nt-protospacer
numbers (which were the historical comparison vs FlashFry) and the new
23-nt-protospacer default:

| Configuration | Wall time | Peak RSS | Notes |
|---|---:|---:|---|
| **FlashFry** `discover` (20 nt, 1 proc) | 2 m 11 s | 1.1 GB | reference baseline |
| **FlashFry** `discover` chunked (20 nt, 8 procs) | 34.6 s | ~9.4 GB total | fairer wall-clock |
| **crispr-ots** (20 nt protospacer, 16 t) | 9.75 s | 7.0 GB | apples-to-apples vs FlashFry |
| **crispr-ots** (23 nt protospacer, 16 t, cap=500) | **4.90 s** | 5.8 GB | new default |
| **crispr-ots** (23 nt protospacer, 16 t, cap=-1) | 5.27 s | 5.8 GB | no cap |

The 23-nt switch was a **2× speedup over the same workload**. The
explanation is a sharper bin-prefix prefilter: the 7-bp bin key
covers 30 % of the protospacer instead of 35 %, but the prefix
selectivity scales with the underlying genome complexity at that
length — 23-nt protospacers spread across more bins with fewer per-bin
collisions, so each guide's per-bin candidate set shrinks. Same `.crot`
file, same site count, no memory cost; pure scan-time win.

`--max-off-targets-per-bin` (default 500) caps distinct off-target
sequences per `(guide, mismatch)` bin. 500 × 5 mismatch bins = ≤ 2 500
sequences per guide. Modest 7 % speedup on the post-scan grouping +
scoring path (the bulk of cost is in the SIMD scan, which still runs
fully). The bigger win is *bounded output* — homopolymer guides used
to emit 800 k+ off-targets each; capped output is suitable for direct
ingestion by downstream tools. Pass `-1` to disable.

### Correctness against FlashFry (off-target set, not specificity)

CFD is **SpCas9-specific** (different mismatch-penalty table, different
PAM model), so per-pair CFD numbers don't transfer to Cas12a. The
comparison here is the **off-target enumeration itself**: per-guide
`otCount` and the full `{(sequence, count, mismatches)}` set returned
by each tool.

Numbers below are from the **20-nt-protospacer crispr-ots run** (the
only apples-to-apples vs FlashFry, since FlashFry can't do 23-nt at
all). 1 028 / 1 034 guides both tools emit:

| FlashFry `overflow` tag | Guides | Off-target set match |
|---|---:|---:|
| `OK` | 752 | **752 / 752 exact** (perfect on every off-target tuple) |
| `OVERFLOW` | 282 | 276 mismatched, all `crispr-ots ⊃ FlashFry` |

FlashFry's discover caps per-guide off-target enumeration around ~2 000
entries; OVERFLOW-tagged guides have their `offTargets` field truncated.
crispr-ots emits the true counts — e.g. the homopolymer
`TTTTTGTTTGTTTTTTGTTTTTTT…` guide actually has ~797 k off-targets at
mm ≤ 4 on mouse, all genuine. No guide had off-targets in FlashFry that
crispr-ots missed.

### Aggregated Cas12a specificity score

crispr-ots ships a dedicated Cas12a scorer behind `--score cfd-cas12a`
(alias `cfd-cas12a:2xnls`) and `--score cfd-cas12a:encas12a`, backed by
the two matrices bundled in `crispr-score/data/`. The matrices are
byte-for-byte copies of `crisprware/parasol_scripts/off_targ_2xNLS_Cas12a.csv`
and `off_targ_enCas12a.csv` — the same data crisprware's
`score_flashfry_cfd.py` pipeline uses today.

Per-pair score is the standard multiplicative CFD form: trim the 4-bp
PAM, walk the **23** protospacer positions, multiply the matrix penalty
for each mismatch. **No PAM weight tail** (unlike SpCas9's `pam_weight`).
The per-guide aggregate is the GuideScan2 convention,
`1 / Σ (cfd × count)`, in two flavors:

- `cas12a_spec_tttn` — sum over all off-targets.
- `cas12a_spec_tttv` — sum excluding TTTT-prefixed off-targets (the
  PAM variant Cas12a cleaves poorly, so excluding it tilts the score
  toward biologically credible off-targets).

The crispr-score crate's `cas12a_parasol_validate` integration test
checks per-pair scores against the parasol-script reference on five
hand-picked cases — including one that exercises PAM-distal positions
21-23 of the matrix; max |Δ| = 4e-11 (print-precision rounding). The
aggregation matches the parasol convention by construction.

## Cluster baseline: prism Genoa-X (2026-05-31, pre-acceleration)

First real numbers on the cluster hardware, replacing the i7-10870H laptop
figures as the reference point for the AVX-512 and GPU work. **This is the
current AVX2 / `x86-64-v3` engine — no AVX-512 path, no GPU yet** — i.e. the
apples-to-apples "before" for the two acceleration tasks.

### Hardware
- **Node:** phoenix-09, AMD **EPYC 9684X** ("Genoa-X", Zen 4, AVX-512 incl.
  `avx512_vpopcntdq`, 1152 MB V-Cache). 64 of 384 logical cores allocated.
- **Build:** `target-cpu=x86-64-v3` (AVX2 `wide::u64x4` + SWAR popcount).
- Reference `GRCm39.fa` and the `.crot` live on the ceph group FS (NFS).

### Index build + enumerate (1 000 NGG guides, mm ≤ 4, CFD, warm cache)

| Step | Wall | Peak RSS | Notes |
|---|---:|---:|---|
| `index --bin-width 11` | 5 m 07 s | 12.4 GB | 6.69 GB `.crot`; I/O-bound on the NFS FASTA read (~1.3× vs laptop — `index` isn't the hot loop) |
| `enumerate -t 1` | 19.58 s | 8.0 GB | scan ≈ 73 % of single-thread wall |
| `enumerate -t 16` | 5.59 s | 8.5 GB | 3.5× over t=1 |
| `enumerate -t 64` | 5.36 s | 9.2 GB | ≈ flat vs t=16 |

### Methodology note that drives the accelerator work
At 1 000 guides the engine **saturates by ~16 threads** (t=64 ≈ t=16). The
residual ~5 s is fixed overhead — mmap page-in of the 6.7 GB index plus
single-threaded hit grouping, scoring, and output — **not** the SIMD scan.
Consequences for the acceleration tasks:

- **AVX-512 VPOPCNTDQ** helps most where the scan dominates: single-/few-thread
  runs, scan-heavy enzymes (Cas12a), or larger guide sets. At t=64 on this
  1 K-guide SpCas9 workload its overall headroom is Amdahl-capped by the ~5 s
  fixed-overhead floor.
- **GPU** likewise needs a scan-bound workload to shine. Benchmark both
  accelerators on **10 K–100 K guides** (scan-dominated) and on Cas12a, with
  the index kept VRAM-resident across queries, to measure their real ceiling.
  A `random_10k.fasta` fixture should be generated for these comparisons; the
  1 K fixture stays for correctness pinning.

## AVX-512 VPOPCNTDQ scan path (2026-05-31)

Added a runtime-dispatched 8-wide scan kernel (`scan_bin_avx512` in
`crispr-scan/src/bin_scanner.rs`): it replaces the portable `wide::u64x4` +
10-op SWAR pair-popcount with `__m512i` and a single hardware
`_mm512_popcnt_epi64` per 64-bit lane, **eight guides per iteration**.
Selected at runtime via `is_x86_feature_detected!("avx512f" / "avx512vpopcntdq")`,
so one binary runs the AVX-512 path on Genoa / Genoa-X / Ice Lake-SP /
Sapphire Rapids and transparently falls back to the AVX2 path on Zen 2 (the
GPU hosts) and anything older. `CRISPR_OTS_NO_AVX512=1` forces the fallback
on the same binary (used for the A/B below). No new dependency (pure
`std::arch`); MSRV bumped to 1.89, the first stable rustc with the 512-bit
intrinsics.

### A/B on Genoa-X (phoenix-09, 1 000 NGG guides, mm ≤ 4, CFD, warm cache)

| Threads | AVX-512 | AVX2 fallback | Speedup | Regime |
|---:|---:|---:|---:|---|
| 1 | 12.35 s | 17.49 s | **1.42×** | scan-bound (scan ≈ 73 % of wall) |
| 16 | 9.69 s | 9.60 s | ~1.0× | memory-bandwidth-bound |

User (compute) time at t=1 is 9.49 s vs 12.85 s = **1.35× less CPU**; the
kernel-level scan speedup is ~1.7× (overall 1.42× after the ~5 s fixed
overhead). The gain vanishes at 16 threads because the scan there is
**memory-bandwidth-bound, not compute-bound** — the run sits at ~234 % CPU
(≈ 2.4 busy cores), and the 6.7 GB working set overflows even Genoa-X's
1152 MB V-Cache, so threads stall on memory and faster popcount can't help.

**Correctness:** AVX-512 output is **byte-identical** to the fallback on the
full mouse workload — **9 497 796 off-target rows match exactly** (`cmp`) —
on top of the `avx512_kernel_matches_x4_kernel` unit test that drives full
8-lane chunks + a remainder across both the SpCas9 and Cas12a masks.

**Takeaway:** a clean win for scan-bound use (single-/few-thread, scan-heavy
enzymes), free on capable hardware, harmless elsewhere. The high-thread
ceiling is set by memory bandwidth, not the kernel — which is exactly the
ceiling the GPU path must beat by moving the working set into the A5500's
~768 GB/s VRAM, resident across queries.

## GPU (CUDA / cudarc) scan path (2026-05-31)

A `crispr-scan-gpu` crate implements the same `Scanner` trait with a CUDA
kernel compiled at runtime via NVRTC. The off-target database — a flat array of
24-byte `PackedEntry` records — is uploaded once and kept VRAM-resident (the
mouse `.crot` is 6.7 GB, well inside the A5500's 24 GB). The kernel runs one
thread per database entry, inner-looping over every guide with the exact same
XOR / compare-mask / 2-bit-pair popcount as the CPU `BinScanner` (hardware
`__popcll`), and atomic-appends matches. No bin-prefix prefilter — brute force
finds the identical hit set and the GPU's throughput makes it cheap. Selected
with `--scanner gpu`; it lives behind an optional `gpu` cargo feature so the
default build is CUDA-free, and cudarc's `dynamic-loading` lets it build even
on the CUDA-less login node and load the real libraries on the GPU nodes.

### Correctness
GPU output is **byte-identical** to the CPU path on the full mouse workload:
1 000 NGG guides, mm ≤ 4, `--max-off-targets-per-bin -1` → **9 497 796 rows
match exactly** (`cmp`). Locked by the `gpu_matches_cpu_when_device_present`
unit test (runs wherever a CUDA device is present, skips elsewhere), which
checks GPU ≡ CPU hit sets at mm 0–4. (With the default cap the two can pick
*different* sequences to drop, because the cap keeps the first N in
hit-iteration order and the GPU's atomic-append order differs — so equivalence
is asserted with the cap disabled.)

### A5500 vs 16-thread CPU (phoenix-01, 1 000 NGG guides, mm ≤ 4, CFD, warm)

| Stage | GPU (A5500) | CPU (EPYC 7662, 16 t) |
|---|---:|---:|
| index load (one-time) | flatten 3.6 s + H2D 0.55 s | — (mmap, paged during scan) |
| **scan kernel, index resident** | **0.98 s** | 6.69 s |
| full `scan_on_gpu` (load + kernel) | 5.97 s | 6.69 s |
| total wall (incl. group/score/write) | 10.7 s | 11.1 s |
| peak RSS | 13.7 GB | 8.4 GB |

**Headline: the resident scan kernel is ~6.8× faster than the 16-thread CPU
scan** (0.98 s vs 6.69 s) over 9.5 M hits — the right metric for the canonical
index-once / query-many pattern (`GpuScanner` uploads the index once and serves
many `scan()` calls at ~1 s each).

**A single `enumerate` is roughly a wash end-to-end.** The GPU kernel is ~7×
faster, but (a) the 6.7 GB index load — single-threaded flatten + upload, ~4 s
warm / ~12 s cold off ceph — and (b) the ~4 s single-threaded
group/score/write tail (shared with the CPU path) dominate the wall. Higher GPU
RSS (13.7 vs 8.4 GB) is the host `Vec` copy of the entries used for upload and
hit reconstruction.

### Scaling to 10 K guides — the advantage *shrinks*, not grows

Replicating the 1 K set 10× (same index, warm) to probe scan-bound scaling.
Counter to the naive "more scan work favors the GPU" expectation, the lead
**erodes** at higher guide counts:

| guides | hits | GPU scan (resident) | CPU scan (16 t) | GPU total | CPU total |
|---:|---:|---:|---:|---:|---:|
| 1 000 | 9.5 M | 0.92 s | 2.6–6.7 s* | 10.8 s | 7.0 s |
| 10 000 | 95 M | 14.3 s | 16.6 s | 69.0 s | 58.7 s |

\* the CPU scan is host-page-cache-sensitive — ~2.6 s fully warm, ~6.7 s under
cache pressure. The GPU scan (~0.9 s at 1 K) is stable because the index is
VRAM-resident, immune to host cache eviction.

Why the GPU's lead shrinks as guides grow:
- **Hit volume is the GPU's weak point.** Output scales with hits (9.5 M →
  95 M); the GPU must atomic-append them, copy them host-ward, and reconstruct
  them on a single core — none of which the kernel accelerates. At 10 K that
  hit-handling dominates the "scan" number.
- **A prototype inefficiency compounds it:** the adaptive output buffer (16 M
  initial cap, grown on overflow) **re-runs the whole kernel** when 95 M > 16 M,
  so the 10 K scan effectively runs the kernel twice. A two-pass
  count-then-write — or, better, a kernel-side per-bin cap that emits far fewer
  hits — removes this; with it the 10 K GPU scan should be ~7 s, not 14 s.
- **The serial post-scan tail dominates both backends.** group/score/write is
  ~42 s of the ~59–69 s total at 10 K for *either* scanner, so end-to-end is
  tail-bound and the GPU's extra index-load + retry overhead pushes it slightly
  behind.

So the GPU win is real **in the regime where kernel compute dominates and hit
volume is modest** — the 1 K-guide, index-once / query-many case. For large
high-hit-count sweeps, the levers below (kernel-side cap, no-retry output,
parallel tail) are prerequisites, not nice-to-haves. (This corrected an earlier
guess that larger sets would simply favor the GPU — they don't, until the hit
pipeline and the tail are fixed.)

### All-vs-all (genome self-scan): on-GPU aggregation + the brute-force ceiling

The real target is scanning **every** genome protospacer against the genome —
millions of guides. Two changes make it even approachable:
- **Output must be per-guide aggregates, not hits.** `GpuScanner::scan_counts`
  accumulates per-guide off-target counts by mismatch level on the GPU
  (`atomicAdd` into an `n_guides × (max_mm+1)` array) — compact regardless of
  hit volume (a full mouse all-vs-all hit list would be ~10¹² rows).
- **Guides must be sampled uniformly** across the index (`allvsall_bench`
  example). The first-n entries are all low-complexity poly-A sites (~62 k
  off-targets each) and badly skew any measurement.

Measured on the A5500 (mm ≤ 4, uniform-sampled real genome sites, index
VRAM-resident):

| guides | scan_counts | hits | per-guide | implied all-vs-all (277 M sites) |
|---:|---:|---:|---:|---:|
| 10 000 | 5.98 s | 96 M | 598 µs | ~46 h |
| 100 000 | 61.0 s | 975 M | 610 µs | ~47 h |
| 1 000 000 | 621.8 s | 9.75 B | 622 µs | ~48 h |

Per-guide time is **flat (~0.6 ms) → cleanly linear**, and the rate
(~4.6×10¹¹ comparisons/s) is **compute-bound at ~peak integer throughput** for
this op mix (~15 ops: XOR/mask/2-bit-pair-popcount over both site words). So
brute force is **~46 h on one A5500**; 8× multi-GPU → ~6 h — still too slow for
routine use.

**The bin-prefilter is the essential next lever.** Brute force compares every
site against all 277 M; the prefilter restricts each site to its candidate
bins (those within k mismatches on the bin-width prefix), cutting comparisons
by **~100×+** → a projected **tens of minutes on one A5500**, then ÷8 with
multi-GPU. The catch is memory access: candidate-bin / bin-offset lookups
scatter, so the win hinges on a **bin-pair-tiled** kernel — load a bin and a
neighbor bin into shared memory, compare all pairs, reuse — and likely a
**larger bin width** than the w=11 used for single-guide queries (a sharper
prefix means fewer candidate bins). That kernel, plus multi-GPU, is the path
from "46 h" to "minutes", and is the next build.

### Bin-prefilter: ~1.3 h all-vs-all per A5500 (and a benchmarking trap)

`GpuScanner::scan_counts_prefilter` — each guide computes its bin key,
recursively enumerates its candidate bins (those within k base-mismatches of
the prefix), scans only those, and accumulates counts in registers (no atomics,
each guide owns its counters). **Validated byte-identical to the brute-force
kernel.** The number of *comparisons* drops ~130× vs brute force.

**The benchmark trap (and the real result).** The kernel's speed is dominated
by **L2 reuse of the candidate bins**, which depends entirely on how the guides
are *sampled*:

| sampling | per-guide | projected all-vs-all |
|---|---:|---:|
| **strided** across the whole index (10 k–1 M) | 278–1442 µs | ~21–111 h |
| **contiguous** middle run (100 k) | 28.2 µs | 2.2 h |
| **contiguous** middle run (1 M) | **16.3 µs** | 1.26 h |
| **contiguous** middle run (5 M) | 17.4 µs | **1.34 h** |

Strided sampling spreads each guide's candidate bins across the whole index →
L2 thrash → 278 µs/guide and a misleading ~21 h projection. **The real
all-vs-all processes entries in bin order — i.e. contiguous — so adjacent
guides share candidate bins and the rate is ~17 µs/guide (16× faster):
the full mouse self-scan is ~1.3 h on one A5500**, hence **~10 min across the
node's 8 GPUs** (shard the guide set; the 6.7 GB index fits resident on each),
and less on the A100. No new kernel required — just run the existing per-guide
prefilter on the index in order.

### Bin-pair tiling: tried as a separate kernel, shelved

`tiled_self_scan_kernel` / `self_scan_counts` processes guides **by bin** (one
block per bin) so the candidate-bin enumeration is shared by the bin's guides —
the idea being to amortise enumeration. **Validated byte-identical to the
per-guide prefilter**, but it is **~100× slower** (4 000 bins took 1 578 s vs
the per-guide kernel doing the same 596 k guides in 16 s). The design
**collapses parallelism**: one block per bin with only the `w*3 = 33`
top-level subtrees active per block (≈130 k active threads vs the per-guide
kernel's one-thread-per-guide), each doing large uneven serial work. With the
contiguous per-guide result above already at ~1.3 h, the tiling is moot and is
kept only as a documented dead-end; the win is the access pattern, not the
decomposition. A *finer-grained* tiling (distribute comparisons, not subtrees,
across a block) could still help but is not worth it over per-guide + sharding.

### Multi-GPU all-vs-all: full mouse genome in 24 min on 4 A5500s

`multi_gpu_allvsall` shards the entries into N contiguous chunks (one per GPU)
and runs the per-guide prefilter on each shard in parallel threads, each with
its own resident `GpuScanner`. Measured on phoenix-01:

| run | guides | wall | grand-total hits |
|---|---:|---:|---:|
| 1 GPU, 3% middle slice | 8.32 M | 177.8 s | 103,883,517,496 |
| 4 GPU, same slice | 8.32 M | 57.5 s | 103,883,517,496 |
| **4 GPU, full self-scan** | **277 M** | **1444.7 s (24.1 min)** | 2,715,565,331,281 |

The sharded grand total is **identical to the single-GPU reference** (disjoint
contiguous shards → sharding is correct). So the **full mouse all-vs-all — all
277 M sites scanned against the genome, ~2.7×10¹² off-target relationships at
mm ≤ 4 — completes in 24 minutes on 4 A5500s** (≈12 min on all 8, linear). Per-
GPU scan times span 1282–1430 s: contiguous sharding balances guide *count* but
not hit density, so the wall is bound by the densest shard; work-balanced
sharding (split by estimated comparisons, not guide count) would tighten it.
From a 46 h brute-force estimate to ~24 min — the whole genome-wide self-scan
now runs over coffee.

### What unlocks the ceiling (future work, in rough priority)
- **Expose index-once / query-many from the CLI** (persistent/server mode):
  realizes the 6.8× directly — the `GpuScanner` API already supports it.
- **Scan-bound workloads** (10 K–100 K guides, Cas12a): the 0.98 s kernel scales
  far better than the CPU's linear scan; the fixed tail amortizes.
- **Zero-copy / parallel index load**: upload straight from `MmapDb::entries()`
  (downcast) or read in parallel — kills the ~4 s flatten and the extra 6.7 GB
  RSS. The flatten is I/O-bound (paging the index in), not copy-bound.
- **Kernel-side per-bin cap + emit `(guide, site, position, mm)` directly**:
  stop emitting the 4 M dropped hits and drop the host entries array entirely.
- **Parallelize the post-scan tail** (group/score/write) — the ~4 s floor that
  bounds *both* backends end-to-end (see the AVX-512 section's same finding).
- **Multi-GPU** (8× A5500 / A100 per node) for linear throughput scaling.

Profiled on the A5500 (compute 8.6, 24 GB, ~768 GB/s); the A100 (40/80 GB,
~1.5–2 TB/s) should push the kernel faster still.

## Snapshot history

Track perf regressions and the impact of each optimization phase here.

| Commit | Date | Workload | Discover wall | Discover RSS | Notes |
|---|---|---|---:|---:|---|
| `ae6fb18` | 2026-05-30 | 1 K, 1 t | n/a | n/a | Initial scaffold (no CLI yet) |
| `aac081b` | 2026-05-30 | 1 K, 1 t | 7.11 s | n/a | Bin-prefix prefilter on top of scalar scan |
| `7c98dab` | 2026-05-30 | 1 K, 1 t | 3.70 s | 241 MB | Build/discover split, persistent DB |
| `c25e14d` | 2026-05-30 | 1 K, 1 t | 2.12 s | 261 MB | `--bin-width 9` (chr22 sweet spot) |
| `32b74f9` | 2026-05-30 | 1 K, 16 t | 0.76 s | 406 MB | Phase 5 rayon, width 9 |
| Phase 6 | 2026-05-30 | 10 K, 1 t | 16.33 s | 1.28 GB | wide::u64x4 SIMD, AVX2 baseline |
| Phase 6 | 2026-05-30 | 10 K, 16 t | **4.68 s** | 2.17 GB | Phase 6 with rayon |
| Phase 4b | 2026-05-30 | 10 K, 1 t | 16.57 s | 1.20 GB | mmap v2 format, DB shrunk 139→116 MB |
| Phase 4b | 2026-05-30 | 10 K, 16 t | 4.81 s | 1.93 GB | mmap v2 format, 16 threads |
| mm39-bench | 2026-05-30 | 1 K mouse, 1 t | 18.6 s | 8.0 GB | Full GRCm39, bin-width 11, cold-cache |
| mm39-bench | 2026-05-30 | 1 K mouse, 16 t | **5.9 s** | 8.5 GB | Full GRCm39, bin-width 11, warm-cache |
| mm39-bench | 2026-05-30 | 1 K mouse Cas12a, 16 t | **9.75 s** | 7.0 GB | Cpf1 TTTN, mouse, bin-width 11 |
| prism-base | 2026-05-31 | 1 K mouse, 1 t | 19.6 s | 8.0 GB | **Genoa-X EPYC 9684X**, AVX2/v3 build, warm; scan ≈73% of wall |
| prism-base | 2026-05-31 | 1 K mouse, 16 t | 5.59 s | 8.5 GB | Genoa-X; 3.5× over t=1 |
| prism-base | 2026-05-31 | 1 K mouse, 64 t | 5.36 s | 9.2 GB | Genoa-X; ≈flat vs t=16 → ~5 s fixed-overhead floor, scan no longer bottleneck |
| avx512 | 2026-05-31 | 1 K mouse, 1 t | **12.35 s** | 8.0 GB | Genoa-X VPOPCNTDQ; 1.42× vs AVX2 fallback (17.49 s); output bit-identical (9.5 M rows) |
| avx512 | 2026-05-31 | 1 K mouse, 16 t | 9.69 s | 8.5 GB | memory-bound; ≈ AVX2 fallback (9.60 s) |
| gpu-a5500 | 2026-05-31 | 1 K mouse, resident kernel | **0.98 s** | 13.7 GB | A5500 scan kernel only; 6.8× vs CPU-16t (6.69 s); output bit-identical (9.5 M rows) |
| gpu-a5500 | 2026-05-31 | 1 K mouse, 1 enumerate (warm) | 10.7 s | 13.7 GB | incl. ~4 s index load + ~4 s serial group/score/write (≈ CPU 11.1 s) |

## Reproducing

1. **Stand up the fixture** (one-time):
   ```bash
   cd /tmp/ffry-quickstart
   tar xzf /home/eric/Projects/FlashFry/FlashFry/test_data/quickstart_data.tar.gz
   gunzip -k chr22.fa.gz
   # Generate random_1000.fasta + random_1000.gs2.csv with the
   # Python helper in NOTES.md.
   ```

2. **Run each tool's build and time it**:
   ```bash
   # FlashFry index
   /usr/bin/time -v java -Xmx4g -jar FlashFry-assembly-1.15.jar index \
       --tmpLocation ./tmp --database chr22_cas9ngg_database \
       --reference chr22.fa.gz --enzyme spcas9ngg

   # GuideScan2 index
   /usr/bin/time -v guidescan index --index chr22.gs2 chr22.fa

   # crispr-ots build
   /usr/bin/time -v crispr-ots build \
       --reference chr22.fa.gz --enzyme spcas9ngg --output chr22.crot
   ```

3. **Run each tool's discover and time it**:
   ```bash
   # FlashFry discover+score
   /usr/bin/time -v bash -c '
       java -jar FlashFry-assembly-1.15.jar discover \
           --database chr22_cas9ngg_database --fasta random_1000.fasta \
           --output random_1000.output >/dev/null 2>&1
       java -jar FlashFry-assembly-1.15.jar score \
           --input random_1000.output --output random_1000.scored \
           --scoringMetrics doench2016cfd --database chr22_cas9ngg_database
   '

   # GuideScan2 enumerate (NGG alt-PAMs)
   /usr/bin/time -v guidescan enumerate \
       --kmers-file random_1000.gs2.csv --mismatches 4 \
       --alt-pam AGG --alt-pam CGG --alt-pam GGG --alt-pam TGG \
       --format csv --mode complete \
       --output random_1000.gs2.ngg.csv.out chr22.gs2

   # crispr-ots discover
   /usr/bin/time -v crispr-ots discover \
       --database chr22.crot --queries random_1000.fasta \
       --max-mismatches 4 --score cfd --output random_1000.ours.tsv
   ```

## Bin-width sweep on chr22

The `--bin-width` flag on `crispr-ots build` controls the size of the
bin-prefix key. Wider bins → smaller per-bin work but a larger bin-offset
table; past a point the prefilter outer loop dominates.

| Width | Build wall | Build RSS | DB size | Scan wall | Scan RSS |
|---:|---:|---:|---:|---:|---:|
| 6 | 3.72 s | 257 MB | 133 MB | 6.60 s | 240 MB |
| 7 (default) | 3.82 s | 267 MB | 133 MB | 3.82 s | 241 MB |
| 8 | 4.16 s | 283 MB | 134 MB | 2.34 s | 245 MB |
| **9** | 4.52 s | 300 MB | 136 MB | **2.12 s** | 261 MB |
| 10 | 4.89 s | 334 MB | 142 MB | 3.50 s | 286 MB |
| 11 | 5.43 s | 449 MB | 154 MB | 6.63 s | 336 MB |
| 12 | 6.00 s | 632 MB | 166 MB | 9.72 s | 419 MB |

Empirical sweet spot on chr22 (~5 M sites): **width 9**, 1.8× faster scan
than the default width 7. All widths produce bit-identical CFD scores;
correctness is preserved at every setting (max |Δ spec| vs FlashFry ≤
2.22e-16 across the sweep).

### Extrapolation to the full human genome

For ~250 M sites (50× chr22), the sweet spot shifts upward by roughly
`log₄(50) ≈ 2.8`. Expected useful range:

| Width | Bins | Sites/bin (~250 M sites) | Bin-table footprint |
|---:|---:|---:|---:|
| 9 | 256 K | ~975 | 4 MB |
| 10 | 1 M | ~244 | 16 MB |
| **11** | 4 M | ~62 | 64 MB |
| 12 | 16 M | ~16 | 256 MB |

A back-of-envelope model (prefilter cost = `bins × guides × ~1 ns`,
inner cost = `bins × pass_rate(w,4) × guides × sites/bin × ~3 ns`)
predicts **width 11** as the sweet spot for the full human genome with
~1 000 guides, beating width 10 by ~1.7× and width 12 by ~2×.

These predictions need empirical confirmation once we have a human DB
to test against — recording here so the actual measurement is a single
`build --bin-width 11 ...` away.

## Build-time tuning

The workspace ships a `.cargo/config.toml` that sets
`target-cpu=x86-64-v3` for x86_64 release builds. That guarantees AVX2,
BMI2, and POPCNT — enough for `wide::u64x4` to lower to true SIMD
instead of scalar emulation. The constraint is satisfied by any x86 CPU
made after 2013 (Haswell, Excavator, and later).

To squeeze more out of newer hardware (especially AVX-512 VPOPCNTDQ for
hardware lane-wise popcount):

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

Apple Silicon / aarch64 binaries are unaffected — NEON is part of the
AArch64 ABI.

## GPU feasibility analysis (2026-05-30)

Profile of the mouse Cas12a 16-thread enumerate (1 000 guides,
`--bin-width 11`, warm cache, 23-nt protospacer, cap=500):

```
[trace] scan 1038 guides: 3.85s (12,160,287 raw hits)
[trace] dropped 7,345,252 off-target sequences past --max-off-targets-per-bin cap
[trace] total: 4.67s
```

The bin-scan SIMD kernel is **77 %** of wall (3.85 s / 4.67 s). GPU
acceleration only matters if it speeds up that kernel by enough to
shrink total wall meaningfully — Amdahl's ceiling for this profile is
~4.3× even with an instant-fast scan.

### Hardware on this box
- **GPU**: NVIDIA GeForce GTX 1650 Ti Mobile (Turing, compute 7.5,
  1 024 CUDA cores, ~2.8 TFLOPS FP32, ~128 GB/s on-card memory bandwidth)
- **VRAM**: **4 GiB** (≈ 3.7 GiB free at idle)
- **CPU**: i7-10870H, 16 logical cores, AVX2 + BMI2 + POPCNT but no
  AVX-512 VPOPCNTDQ
- **CPU RAM**: 62 GiB; mouse `.crot` files mmap'd from NVMe

### VRAM is the deciding constraint
| Index | On-disk | Fits in 4 GiB VRAM? |
|---|---:|---:|
| chr22 SpCas9, width 9 | 116 MB | trivially |
| Mouse Cas12a (23 nt) | 4.0 GB | ❌ no headroom |
| Mouse SpCas9 | 6.3 GB | ❌ way over |

For mouse scale on this machine, we'd have to **stream bins from RAM
over PCIe** rather than keep the database resident. Laptop PCIe is
typically x4 Gen3 (≈ 4 GB/s); streaming the full 4-6 GB index per
query costs **1.0-1.5 s of PCIe just for the read**. That's a hard floor
before any kernel runs — bigger than the headroom we have over the
4.67 s CPU baseline.

### Theoretical kernel speedup (ignoring transfer)
The bin-scan inner loop is roughly:
```
mm = popcount_pairs((guide.high XOR ot.high) & mask.high)
   + popcount_pairs((guide.low  XOR ot.low ) & mask.low)
if mm <= max_mismatches: append_hit
```

~10 simple integer ops per (guide, bin_entry) pair. Mouse Cas12a at
mm ≤ 4: ~12 M raw hits over ~88 M sites × prefilter pass-rate ≈
several billion pair evaluations.

- GPU compute ceiling: 2.8 TIOPS / 10 ops = ~280 G pair/s → kernel
  finishes in ~20 ms.
- GPU memory-bandwidth ceiling: 128 GB/s. Each pair touches 16 B of
  bin-entry data; 4 GB of data → 30 ms.

So an *in-VRAM* kernel could run in **20-30 ms** vs the CPU's 3.85 s.
That's ~150× on the inner kernel, ~4× on overall wall (Amdahl).

### Practical bottleneck on this hardware: PCIe ingress
Adding PCIe-streaming bins:
- Total bytes to ship: 4-6 GB (one .crot pass per scan)
- Effective laptop x4 Gen3: 4 GB/s
- Streaming time: 1.0-1.5 s
- Total expected GPU wall: 1.5 s + 30 ms = **~1.5 s** vs CPU 4.67 s

≈ 3× overall speedup *if* the streaming and compute overlap perfectly,
which on this card requires careful tile pipelining. Not worth the
implementation cost given the existing 4.67 s baseline already beats
GuideScan2 by 13× and FlashFry chunked-8-procs by 3.5×.

### When GPU starts to win
A modern desktop GPU changes the calculus:
| GPU | VRAM | Bandwidth | Mouse `.crot` fits? |
|---|---:|---:|---:|
| RTX 3060 | 12 GB | 360 GB/s | yes |
| RTX 4070 | 12 GB | 504 GB/s | yes |
| RTX 4090 | 24 GB | 1008 GB/s | comfortably, plus headroom for two genomes |

On any of these, the .crot stays VRAM-resident across queries (the
canonical use-case: index once, run thousands of scans). PCIe is paid
only on the first scan after open; subsequent scans are memory-bound
on the GPU side at ~30-100 ms each. **Expected: 30-50× speedup on the
scan kernel, 4× overall.**

### Recommendation
Don't prototype on the 1650 Ti — the 4 GiB VRAM kills mouse-scale and
the chr22 numbers we already have (0.76 s at 16 threads) are too small
to evaluate GPU benefit cleanly. When the production hardware is
chosen (a workstation GPU with ≥ 12 GiB), a `wgpu` or `cudarc`
prototype of the bin-scan inner loop is a tractable 1-2-week effort and
should land 3-5× overall on mouse-scale workloads. Until then, the
**single highest-impact CPU follow-up** is the AVX-512 VPOPCNTDQ path
(listed below) — comparable gains on the hottest line of the existing
SIMD scan, no new dependency.

## Outstanding work that should move these numbers

- **Phase 4b** — mmap-friendly on-disk format. Target: shrink the
  crispr-ots DB from 139 MB to ≤ FlashFry's 42 MB (likely with a custom
  binary layout + optional zstd of the cold-storage form), and load it in
  ≤ 100 ms instead of bincode-deserializing the whole thing on startup.
- **Phase 5** — rayon parallelism on the bin-scan outer loop. Expected
  4–8× discover speedup on this 16-core box; brings discover under 1 s.
- **Phase 6** — SIMD-batched mismatch (AVX-512 / AVX2 / NEON via
  `multiversion`). Another 4–8× on the inner mismatch loop on top of
  rayon.
- **Phase 7** — adaptive bin width + GuideScan2-style threshold-then-
  enumerate. Smaller wins.

Each phase commit should append a row to the snapshot-history table above.
