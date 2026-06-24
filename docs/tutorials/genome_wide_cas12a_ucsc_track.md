# Genome-wide Cas12a UCSC track

This tutorial builds a **genome-wide Cas12a (Cpf1) track for the UCSC Genome Browser**
from a reference FASTA: every TTTV protospacer in the genome, scored for on-target
efficiency (four models) and off-target specificity against **two CFD matrices**
(EnCas12a and AsCas12a/2xNLS), packaged as a self-contained bigBed track.

It is written around the GRCh38 primary assembly (canonical chromosomes), but the exact
same commands work for any genome. The whole-genome scoring is a **single
`score_guides --ucscgb` command** that internally splits the work into chunks, fans the
GPU off-target passes across every GPU it is given, and merges the pieces into one track.

```{contents} On this page
:local:
:depth: 2
```

## What you get

Three files that together form one UCSC "bigDataUrl" track:

| File | What it is |
|------|------------|
| `cas12a.bb` | bigBed (bed9+) of every guide: position, strand, on-target color, and the per-guide score/specificity/efficiency columns shown on hover and on the details page. |
| `crisprDetails.tab.gz` (+ `.gzi`) | bgzip-compressed, byte-offset-indexed per-guide off-target lists (the long tables shown when you click a guide). |
| `cas12aTargets.as` | the autoSql schema describing the bigBed's columns. |

Each guide row carries: EnCas12a + AsCas12a TTTV/TTTN **specificity** (CFD), and four
**on-target** efficiencies (DeepCpf1, EnPAM-GB, EnCas12a-DeepCpf1, AsCas12a-DeepCpf1).

## Prerequisites

- **crisprware** installed (this repo), and its conda env active.
- The **GPU-built** crispr-ots engine (`crispr-ots` compiled `--features gpu`). The
  off-target pass JIT-compiles CUDA via nvrtc at runtime, so make `libcuda`/`libnvrtc`
  reachable (e.g. `export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH`). A
  CPU-only engine also works with `--ucscgb_scanner cpu`, just slower.
- UCSC `bedToBigBed` and `bgzip` on `PATH` (the track assembly shells out to them).
- One or more NVIDIA GPUs for genome-scale runs (a single workstation GPU is fine for a
  chromosome).

```{tip}
Confirm your engine is GPU-capable: `crispr-ots enumerate --scanner gpu ...` should not
print *"requires a binary built with `--features gpu`"*.
```

## Step 0 - choose and filter the reference

Decide which sequences belong in the track. For a clean human track you usually want the
**primary assembly's canonical chromosomes only** (no chrM, no unplaced/unlocalized
`GL*/KI*` scaffolds, no ALT/patch contigs):

```bash
GENOME=GRCh38.primary_assembly.genome.fa        # already excludes ALT/patch
samtools faidx "$GENOME" \
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
  chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
  > hg38_primary_chr.fna
samtools faidx hg38_primary_chr.fna
cut -f1,2 hg38_primary_chr.fna.fai > hg38.chrom.sizes      # chrom.sizes for bedToBigBed
```

The `chrom.sizes` file is required later by `bedToBigBed`.

```{note}
Off-targets are only searched within the sequences you index here. Dropping scaffolds
makes the track cleaner but means a guide's off-targets on those scaffolds are not
counted. Keep them if that matters for your use case.
```

## Step 1 - build the off-target index

```bash
crisprware index_genome \
  -f hg38_primary_chr.fna \
  --pam TTTN -l 23 --pam_5_prime --bin_width 14 \
  -o .
```

Why these flags:

- **`--pam TTTN`** (not TTTV). EnCas12a also cleaves the weak `TTTT` PAM, so the
  off-target index must contain all four `TTT[ACGT]` sites; otherwise `TTTT`
  off-targets are invisible. (The *guides* are still TTTV-only - Step 2.)
- **`-l 23`** - 23-nt protospacer (Cas12a).
- **`--pam_5_prime`** - Cas12a's PAM is 5' of the protospacer (unlike Cas9's 3' NGG).
- **`--bin_width 14`** - the database bin width; 14 is the measured speed/memory optimum
  and fits comfortably on a 24 GB GPU.

Output: `hg38_primary_chr_gscan2/hg38_primary_chr_gscan2.crot` (the index prefix you pass
to `-i` later).

## Step 2 - enumerate every TTTV protospacer

```bash
crisprware generate_guides \
  -f hg38_primary_chr.fna \
  --pam TTTV --pam_5_prime -5 19 -3 23 -l 23 -w 4 3 \
  --threads 16 \
  -o .
```

Why these flags:

- **`--pam TTTV`** - the canonical Cas12a guide PAM set (`V` = A/C/G; excludes the weak
  `TTTT`).
- **`-w 4 3`** emits a **34-nt context** per guide (4 upstream + 4-nt PAM + 23-nt
  protospacer + 3 downstream) - exactly the window the on-target models
  (DeepCpf1, EnPAM-GB, EnseqDeepCpf1) expect.
- **`--threads 16`** uses the window-parallel engine (large genome-wide speedup).

Output: `hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed` (one row per guide; this is the
`-b` input to scoring).

## Step 3 - score + build the track (the single command)

This one command runs the on-target models, the off-target pass against **both** CFD
matrices, and assembles the three track files. On a genome it chunks internally and fans
across GPUs:

```bash
crisprware score_guides \
  -b hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed \
  -i hg38_primary_chr_gscan2/hg38_primary_chr_gscan2 \
  --cas12a_scorer deepcpf1 enpam_gb enseq_deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a \
  --ucscgb track_hg38_cas12a \
  --chrom_sizes hg38.chrom.sizes \
  --crispr_ots_bin /path/to/target-gpu/release/crispr-ots \
  --ucscgb_scanner gpu \
  --ucscgb_cfd_threshold 0.023 \
  --ucscgb_list_cap 1000 \
  --ucscgb_blank_threshold 0 \
  --mismatches 4 \
  --threads 8 \
  --ucscgb_chunk_max 3000000 \
  -o .
```

What each part does:

- **`--ucscgb track_hg38_cas12a`** selects the UCSC-track output path (instead of the
  default TSV). It writes `cas12a.bb`, `crisprDetails.tab.gz` (+`.gzi`), and
  `cas12aTargets.as` into that directory.
- **`--cas12a_scorer deepcpf1 enpam_gb enseq_deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a`**
  computes the four on-target columns:
  - `deepcpf1` -> **DeepCpf1** (Kim 2018, wild-type AsCas12a/LbCas12a),
  - `enpam_gb` -> **EnPAM-GB** (en(As)Cas12a),
  - `enseq_deepcpf1` -> **EnCas12a-DeepCpf1** (Chen 2025),
  - `seq_deepcpf1variants --cas12a_variant AsCas12a` -> **AsCas12a-DeepCpf1** (the
    variant-specific column that drives the red/yellow/green item color).
- **off-target / two matrices.** The off-target pass always scores the **EnCas12a** CFD
  matrix; `--ucscgb_2xnls` (on by default) adds the **AsCas12a / 2xNLS** matrix in the
  same scan. You get `TTTV_enCas12a_specificity`, `TTTN_enCas12a_specificity`,
  `TTTV_AsCas12a_specificity`, `TTTN_AsCas12a_specificity`. Add `--no-ucscgb_2xnls` to
  score only EnCas12a.
- **`--ucscgb_cfd_threshold 0.023`** is the CFD floor for which off-targets are *listed*
  in `crisprDetails` (the per-bucket mismatch *counts* are unfloored).
- **`--ucscgb_list_cap 1000`** caps the stored off-target list per guide (top by CFD) to
  bound disk; **`--ucscgb_blank_threshold 0`** disables the "too many off-targets ->
  empty list" behavior (every guide keeps its capped list).
- **`--mismatches 4`** is the off-target search radius.
- **`--ucscgb_chunk_max 3000000`** turns on the chunked/multi-GPU path (next section).
  Omit it (or set `0`) to score the whole genome in a single pass on one GPU.

Output: `track_hg38_cas12a/` with the three files.

## The chunked / multi-GPU engine

A whole genome is ~100+ million guides - too much for one GPU pass in reasonable time,
and the single per-guide off-target file would be multiple TB. `--ucscgb_chunk_max`
makes `score_guides --ucscgb` do the fan-out for you:

1. **split** the guide BED into `<=N`-guide chunks,
2. **score** each chunk as an isolated single-pass `--ucscgb` run, one chunk per GPU,
   deleting that chunk's large off-target intermediates as soon as it assembles (so peak
   disk is only the handful of chunks running at once),
3. **merge** the per-chunk tracks into one, recomputing every `crisprDetails.tab` byte
   offset across the concatenated pieces and re-sorting the bigBed.

```bash
# Fan across specific GPUs:
crisprware score_guides ... --ucscgb_chunk_max 3000000 --ucscgb_gpus 0,1,2,3,4,5,6,7
```

- **`--ucscgb_gpus`** lists the devices to use, one persistent worker per GPU, so
  concurrent chunks never collide on a device. If you omit it, the GPUs are
  **auto-detected from `CUDA_VISIBLE_DEVICES`** - i.e. exactly the devices your SLURM job
  (or your shell) was given. On an 8-GPU allocation you simply do not pass the flag.
- **`--ucscgb_keep_chunks`** keeps the `<ucscgb>/chunks/` working directory after the
  merge (for inspection or resume); by default it is removed on success.

```{important}
**Percentiles are per-chunk.** The bigBed's numeric *score*, item *color*, and the
percentile prefixes (e.g. `63% (0.6275)`) are ranked **within each chunk**, not across
the whole genome (the merge does not re-rank). The raw values in parentheses - the actual
CFD and on-target model outputs - are exact and identical no matter how you chunk. If you
need globally-ranked percentiles, score in a single pass (`--ucscgb_chunk_max 0`).
```

### Equivalence and correctness

The chunked path is byte-for-byte equivalent to the manual
"`chunk_guides.py` + per-chunk `score_guides` + `merge_tracks.py`" workflow (it reuses the
same single-pass scorer and the same merge). On a chr22 control (50,000 guides, 2 chunks
vs. one pass) every population-invariant value matches exactly - coordinates, sequence,
unique flags, raw CFD specificity, raw on-target scores, and the entire off-target detail
content - while only the per-chunk percentile columns rebase, as expected. Every merged
`_offset` resolves to the correct off-target line.

## Running on a SLURM cluster

Genome-scale, on one 8-GPU node (GPUs auto-detected, so no `--ucscgb_gpus`):

```bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:8
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=48:00:00
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
crisprware score_guides ... --threads 8 --ucscgb_chunk_max 3000000   # 8 chunks x 8 threads = 64
```

Rules of thumb (A5500-class GPUs):

- `--threads` is *per chunk*; size `cpus-per-task` ~= `#GPUs x --threads`.
- Peak disk ~= `(#concurrent chunks) x (per-chunk off-target file)`. A 3M-guide chunk's
  intermediate is ~100+ GB, so 8 in flight is ~1 TB transient - put the output dir on a
  big filesystem.
- Memory ~= one scored chunk DataFrame + model + off-target file per worker; ~30 GB/GPU
  is a safe budget.

## Loading the track in UCSC

Host the three files at a public URL and add a `bigDataUrl` track (Cas12a uses
`bed9+` with an external detail table and the autoSql above). The viewer fetches the 27-nt
site (`[chromStart, chromStart+27)`), reverse-complements on the `-` strand, and reads the
4-nt PAM + 23-nt protospacer; clicking a guide pulls its off-target list out of
`crisprDetails.tab.gz` via the `.gzi` index using the row's `_offset`.

## Reference: the standalone scaling scripts

The same chunk-and-merge is available as standalone scripts under
`scripts/ucsc_gb_scripts/` (`chunk_guides.py`, `merge_tracks.py`) for cases the built-in
option does not cover - chiefly **multi-node** fan-out via a SLURM job array (the built-in
option uses the GPUs of a single node). `merge_tracks.py` and the built-in path share one
implementation (`crisprware.ucsc_track.merge_ucscgb_tracks`), so their output is identical.

## Worked example: GRCh38 primary assembly

- Reference filtered to 24 canonical chromosomes: **3,088,269,832 bp**.
- Guides (TTTV): **135,867,611** (window-parallel `generate_guides`, ~8.5 min).
- Track build: one `score_guides --ucscgb --ucscgb_chunk_max 3000000` command across
  **7x A5500** (GPUs auto-detected from the SLURM allocation), 46 chunks, **17.3 h** wall
  (peak job RSS ~157 GB; ~1.6 TB transient disk, bounded; the 46-piece merge took 72 min).
- **106,017,484 (78%)** of guides are uniquely placed and carry off-target detail lists;
  the rest have a perfect-match site elsewhere and store an empty list.
- Per-chunk timing: on-target ~2 min, off-target GPU enumerate ~30-42 min, track assembly
  (the wall-time bottleneck) ~90-150 min.
- Outputs in `track_hg38_cas12a/`: `cas12a.bb` (11.7 GB), `crisprDetails.tab.gz`
  (78.7 GB) + `.gzi`, `cas12aTargets.as`.
