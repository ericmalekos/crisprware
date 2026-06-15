# ucsc_gb_scripts — genome-scale Cas12a UCSC track helpers

Helpers for building a genome-wide `score_guides --ucscgb` Cas12a track when the
guide set is too large to score in one pass: score per-chromosome (or per-chunk)
in parallel, then merge the per-piece tracks into a single track. Built for the
T2T-CHM13v2 genome-wide track (~141 M guides).

## Why two scripts

`score_guides --ucscgb` builds one self-contained track (`minimumCas12A.bb`,
`crisprDetails.tab.gz` + `.gzi`, `cas12aTargets.as`) per run, and the GPU
off-target `enumerate` is one-GPU-per-run. For a whole genome you therefore want to
**fan the scoring out** (per chromosome or finer) and **merge** the results — but the
bigBed `_offset` column points at uncompressed byte offsets in that piece's
`crisprDetails.tab`, so a merge must recompute every offset.

## Workflow

1. `generate_guides` over the whole genome → one guide BED. Split it by chromosome
   (e.g. `awk` on column 1, keeping the header) into `bed_by_chrom/<chrom>.bed`.
2. **`chunk_guides.py`** — split those per-chrom BEDs into `<=--chunk-max`-guide
   chunks (each keeps the header) plus a `manifest.txt`, so no single scoring job is
   unbounded (the per-chrom *assembly* is the wall-time bottleneck):
   ```bash
   python chunk_guides.py --bed-dir bed_by_chrom --out-dir . --chunk-max 3000000 \
       --exclude chr21 chr22          # e.g. skip chromosomes already scored
   ```
3. Score each chunk on one GPU, throttled with a SLURM array (`--array=0-(N-1)%C`;
   each task picks `sed -n "$((SLURM_ARRAY_TASK_ID+1))p" chunks/manifest.txt`), e.g.
   `score_guides -b chunks/<chunk>.bed -i <index> --ucscgb track_<chunk>
   --ucscgb_list_cap 1000 --ucscgb_scanner gpu ...`. Delete each chunk's
   `offtargets.csv.ot.tsv` after its assembly to bound peak disk.
4. **`merge_tracks.py`** — merge the per-chunk tracks into one, recomputing each
   guide's `crisprDetails.tab` byte offset across the concatenated pieces, then
   re-sort + `bedToBigBed` and `bgzip` + index:
   ```bash
   python merge_tracks.py --tracks track_chr1_p00 track_chr1_p01 ... \
       --out track_merged --chrom-sizes chrom.sizes \
       --as-file track_chr1_p00/cas12aTargets.as
   ```
   Merge order is arbitrary — offsets are assigned at merge time and the bigBed is
   sorted at the end. Requires `bedToBigBed` and `bgzip` on `PATH`.

## Track schema

`cas12aTargets.as` here is the canonical autoSql schema for the bigBed — 23 fields:
bed9 + `guideSeq`/`pam`, the unique flags, EnCas12a + AsCas12a TTTV/TTTN specificity,
the four on-target efficiency scores, and the required `_mouseOver`/`_offset` trailer.
It is emitted by `score_guides --ucscgb` (every per-chunk track's `cas12aTargets.as` is
an identical copy) and is what you pass to `merge_tracks.py --as-file` (or `bedToBigBed
-as`). Keep it in sync with `crisprware/ucsc_track.py:build_autosql`, which generates it.

## Requirements

Python 3 (stdlib only). `merge_tracks.py` shells out to `bgzip` and `bedToBigBed`.
