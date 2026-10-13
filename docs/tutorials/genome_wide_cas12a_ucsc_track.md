# Genome-wide Cas12a UCSC track

This tutorial builds a **genome-wide Cas12a (Cpf1) track for the UCSC Genome Browser**
from a reference FASTA: every TTTV protospacer in the genome, scored for on-target
efficiency and off-target specificity, packaged as a self-contained bigBed track. By
default `score_guides --ucscgb` builds the **AsCas12a** track (`--ucscgb_enzyme ascas12a`):
one AsCas12a (WT/2xNLS) CFD specificity over a **TTTV-only** off-target index, plus two
AsCas12a on-target models. A separate **enCas12a** track is available via
`--ucscgb_enzyme encas12a` (covered at the end of this page).

It is written around the GRCh38 primary assembly (canonical chromosomes), but the exact
same commands work for any genome. The whole-genome scoring is a **single
`score_guides --ucscgb` command** that internally splits the work into chunks, fans the
GPU off-target passes across every GPU it is given, and merges the pieces into one track.

```{warning}
**Cas12a CFD correction (2026-07-18).** A strand-complement bug in the Cas12a CFD scorer made every
Cas12a specificity/CFD number produced before 2026-07-18 too low (up to ~23x). The worked-example
specificity values in this tutorial predate the fix and are illustrative only; rebuild with
`crispr-ots --version` 0.3.0+ for correct values. See `CAS12A_CFD_CORRECTION.md` at the repo root.
SpCas9 is unaffected.
```

```{contents} On this page
:local:
:depth: 2
```

## What you get

Three files that together form one UCSC "bigDataUrl" track (default `--ucscgb_enzyme ascas12a`):

| File | What it is |
|------|------------|
| `ascas12a.bb` | bigBed (`bed9+12`, 21 fields) of every guide: position, strand, three colour columns (`itemRgb`, `_crisprDeepCpf1Color`, `_specColor`, selectable in the browser via the trackDb `colorFields` setting), and the per-guide specificity and on-target columns shown on hover and on the details page. Each on-target model contributes a raw column plus a `...Perc` percentile sibling, matching the Cas9 track's `fusi`/`fusiPerc` layout. |
| `ascas12a_crisprDetails.tab.gz` (+ `.gzi`) | bgzip-compressed, byte-offset-indexed per-guide off-target lists (the long tables shown when you click a guide). |
| `ascas12aTargets.as` | the autoSql schema describing the bigBed's columns. |

Each guide row carries one **AsCas12a (WT/2xNLS) CFD specificity** (raw, 0-1, over the TTTV
off-target index) and two **on-target** efficiencies: DeepCpf1 (Kim 2018) and the AsCas12a
seq-DeepCpf1 variant (Chen 2025). (The `encas12a` and legacy `dual` tracks carry different
column sets; see the "Reading the track" pages.)

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
  --pam TTTV -l 23 --pam_5_prime --bin_width 14 \
  -o .
```

Why these flags:

- **`--pam TTTV`** (the strong canonical Cas12a PAM set; `V = A/C/G`). The AsCas12a track's
  off-target search is **TTTV-only**: WT/2xNLS AsCas12a does not cleave the weak `TTTT` PAM,
  so `TTTT` sites are not off-targets and must not count against a guide's specificity.
  Indexing `TTTV` (not `TTTN`) is what makes the index match the enzyme. (Build a broader
  `TTTN` or fa>0.3 index only for the legacy `dual` or the `encas12a` track; see the enCas12a
  section.)
- **`-l 23`** - 23-nt protospacer (Cas12a).
- **`--pam_5_prime`** - Cas12a's PAM is 5' of the protospacer (unlike Cas9's 3' NGG).
- **`--bin_width 14`** - the database bin width; 14 is the measured speed/memory optimum
  and fits comfortably on a 24 GB GPU.

Output: `hg38_primary_chr_crisprots/hg38_primary_chr_crisprots` (the index prefix you pass
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
  -i hg38_primary_chr_crisprots/hg38_primary_chr_crisprots \
  --cas12a_scorer deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a \
  --ucscgb track_hg38_cas12a \
  --ucscgb_enzyme ascas12a \
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
  default TSV). With the default `--ucscgb_enzyme ascas12a` it writes `ascas12a.bb`,
  `ascas12a_crisprDetails.tab.gz` (+`.gzi`), and `ascas12aTargets.as` into that directory.
- **`--ucscgb_enzyme ascas12a`** (the default; may be omitted) builds the single-enzyme
  **AsCas12a** track: one specificity over the TTTV off-target index, no enCas12a/TTTN
  columns. The other choices are `encas12a` (see the enCas12a section) and `dual` (legacy).
- **`--cas12a_scorer deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a`** computes the
  AsCas12a track's **two** on-target columns:
  - `deepcpf1` -> **DeepCpf1** (Kim 2018, wild-type AsCas12a/LbCas12a); its percentile drives
    the red/yellow/green item color,
  - `seq_deepcpf1variants --cas12a_variant AsCas12a` -> **AsCas12a seq-DeepCpf1** (Chen 2025).
- **off-target / one matrix.** For the AsCas12a track the off-target pass scores the single
  **AsCas12a (2xNLS)** CFD matrix against the TTTV-only index, yielding one raw
  `AsCas12a_specificity` (0-1) per guide. **No `--ucscgb_genome` is needed**: a TTTV index has
  no `TTTT` off-targets to exclude, so the AsCas12a-TTTV correction is intrinsic.
  (`--ucscgb_genome` is required only for `encas12a`, and is used by the legacy `dual` track to
  read each off-target's reference PAM and drop `TTTT` sites from its broader `TTTN` scan.)
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
**Percentiles are genome-wide, however you chunk.** Each chunk is scored in isolation, but the
merge accumulates a histogram over every chunk and re-ranks the `deepcpf1Perc` and
`ascas12aDeepcpf1Perc` columns, together with the two efficiency colour fields those percentiles
drive (`itemRgb` and `_crisprDeepCpf1Color`), against the whole genome before the final bigBed is
written. `--ucscgb_chunk_max` therefore has no effect on any published value, and no post-pass
re-ranking is required.

Ranking is on the value **as the track displays it** (rounded to 4 decimals), not the
full-precision float behind it, which is what makes the chunked and unchunked paths agree
exactly. `score`, the raw `AsCas12a_specificity` and `_specColor` are not ranks at all: they are
absolute functions of the raw CFD specificity, so they are identical under any chunking.
```

### Equivalence and correctness

The chunked path reuses the same single-pass scorer and the same merge as the manual
"`chunk_guides.py` + per-chunk `score_guides` + `merge_tracks.py`" workflow. Since the merge
re-ranks genome-wide, a chunked build now matches an unchunked one on the percentile and colour
columns as well, not just the population-invariant ones.

Verified on a full sacCer3 control (590,425 guides, 6 chunks vs. one pass): `itemRgb`,
`_crisprDeepCpf1Color`, `_specColor`, `score`, `AsCas12a_specificity`, `deepcpf1Perc`,
`ascas12aDeepcpf1Perc` and `_mouseOver` are identical in both builds. The only differences are 25
guides whose raw on-target floats move in the 4th decimal place, which is CNN batch-size
nondeterminism in the on-target models and not a chunking artefact. Every merged `_offset`
resolves to the correct off-target line.

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

## Building the enCas12a track instead (`--ucscgb_enzyme encas12a`)

The default track is AsCas12a. To build the **enCas12a** (enAsCas12a) track instead, four things
change: the guides use enAsCas12a's high-activity **fa>0.6 (15-PAM)** on-target set, the off-target
index becomes the broad **fa>0.3 (30-PAM)** set (enAsCas12a cleaves far more than `TTTV`),
`--ucscgb_genome` becomes **required**, and the on-target models switch to EnPAM-GB + the
enAsCas12a-HF1 variant.

**1. A fa>0.3 (30-PAM) off-target index.** enAsCas12a's accepted PAMs are not a single IUPAC
pattern, so build the index from an explicit PAM list (`--pam_list_file`) rather than `--pam`.
The 30 PAMs with frac-active > 0.3 are the `in_fa0.3` rows of the shipped weights table
(`crisprware/data/encas12a_offtarget_weights.tsv`):

```bash
# 30 enAsCas12a PAMs with frac-active > 0.3 (column 6 == 1 in the shipped weights table)
awk -F'\t' 'NR>1 && $6=="1" {print $1}' \
  crisprware/data/encas12a_offtarget_weights.tsv > encas12a_fa03.pams

crisprware index_genome \
  -f hg38_primary_chr.fna \
  --pam_list_file encas12a_fa03.pams -l 23 --pam_5_prime --bin_width 14 \
  -o idx_encas12a
```

Output prefix: `idx_encas12a/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots` (use a
separate `-o` so it does not collide with the TTTV index from Step 1).

**2. Guides over the fa>0.6 on-target set.** enCas12a *guides* use the high-activity **fa>0.6**
PAM set (the 15 PAMs enAsCas12a cleaves at >0.6 relative activity), not the broad fa>0.3 off-target
set, so weaker-PAM target sites are not enumerated as guides. These are the `in_fa0.6` rows of the
shipped weights table (a superset of `TTTV`):

```bash
# 15 enAsCas12a on-target PAMs with frac-active > 0.6 (column 5 == 1 in the shipped weights table)
awk -F'\t' 'NR>1 && $5=="1" {print $1}' \
  crisprware/data/encas12a_offtarget_weights.tsv > encas12a_fa06.pams

crisprware generate_guides \
  -f hg38_primary_chr.fna \
  --pam_list_file encas12a_fa06.pams --pam_5_prime -5 19 -3 23 -l 23 -w 4 3 \
  --threads 16 -o encas12a_guides
```

**3. Score with the enCas12a enzyme.** `--ucscgb_genome` is **required** here: each
off-target's reference PAM is reconstructed from the FASTA to weight the specificity by PAM
activity.

```bash
crisprware score_guides \
  -b encas12a_guides/hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed \
  -i idx_encas12a/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots \
  --cas12a_scorer enpam_gb seq_deepcpf1variants --cas12a_variant enAsCas12a-HF1 \
  --ucscgb track_hg38_encas12a \
  --ucscgb_enzyme encas12a \
  --ucscgb_genome hg38_primary_chr.fna \
  --chrom_sizes hg38.chrom.sizes \
  --crispr_ots_bin /path/to/target-gpu/release/crispr-ots \
  --ucscgb_scanner gpu \
  --ucscgb_cfd_threshold 0.023 \
  --ucscgb_list_cap 1000 --ucscgb_blank_threshold 0 \
  --mismatches 4 --threads 8 --ucscgb_chunk_max 3000000 \
  -o .
```

What is different from the AsCas12a build:

- **Files:** `encas12a.bb` (`bed9+11`, 20 fields), `encas12a_crisprDetails.tab.gz` (+`.gzi`),
  `encas12aTargets.as`.
- **Four specificities, all raw (0-1):** `enCas12a_spec_fa06`, `enCas12a_spec_fa06_wt`,
  `enCas12a_spec_fa03`, `enCas12a_spec_fa03_wt`. Each is the guidescan specificity
  `1 / (on_count + sum CFD[x weight])`, evaluated raw and PAM-activity-weighted (each
  off-target's CFD scaled by its PAM's frac-active weight), at the fa>0.6 (high-activity) and
  fa>0.3 (moderate-activity) PAM cutoffs. The scan runs **unfloored** so these sums are exact;
  `--ucscgb_cfd_threshold` then floors only the displayed popup list.
- **Score and color:** `enCas12a_spec_fa03_wt` (the fa>0.3 PAM-weighted specificity) drives
  the field-5 `score` and the low-specificity greying; the color efficiency ramp is the
  **EnPAM-GB** percentile. A guide with a perfect (0-mismatch) off-target elsewhere is blanked,
  marked `unique = False`, and greyed.
- **On-target models:** `enpam_gb` -> **EnPAM-GB** (Kim 2020), and
  `seq_deepcpf1variants --cas12a_variant enAsCas12a-HF1` -> the higher-fidelity
  **enAsCas12a-HF1** seq-DeepCpf1 variant (Chen 2025).

The chunked/multi-GPU engine and the UCSC hosting steps are identical to the AsCas12a build, with
one difference that matters: **the genome-wide percentile re-grade is AsCas12a-only.** A chunked
enCas12a build leaves its percentiles and efficiency-ramp colours ranked within each chunk and
prints `WARNING: no percentile re-grading rule for enzyme='encas12a'` when it merges. The legacy single-track **dual** build
(`--ucscgb_enzyme dual`, which carries both enzymes' TTTV/TTTN specificities in one 23-field
`cas12a.bb`) is retained for back-compat only.

## Loading the track in UCSC

Host the three files at a public URL and add a `bigDataUrl` track (Cas12a uses
`bed9+` with an external detail table and the autoSql above). The viewer fetches the 27-nt
site (`[chromStart, chromStart+27)`), reverse-complements on the `-` strand, and reads the
4-nt PAM + 23-nt protospacer; clicking a guide pulls its off-target list out of the track's
`*_crisprDetails.tab.gz` (e.g. `ascas12a_crisprDetails.tab.gz`) via the `.gzi` index using the
row's `_offset`.

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
- Outputs in `track_hg38_cas12a/`: `ascas12a.bb`, `ascas12a_crisprDetails.tab.gz` + `.gzi`,
  `ascas12aTargets.as`. (A comparable legacy dual-matrix run produced an 11.7 GB `.bb` and a
  78.7 GB details file; the 17-field AsCas12a track is somewhat smaller.)
