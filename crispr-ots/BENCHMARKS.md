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
