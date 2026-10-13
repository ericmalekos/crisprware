# Score and assemble the track

One command does everything: on-target scoring, off-target enumeration under the enzyme's CFD matrix, and
assembly of the three browser files. With `--ucscgb_chunk_max` it also **chunks the genome, fans the
off-target passes across every allocated GPU, and merges the pieces** into a single track. Pick the
nuclease with `--ucscgb_enzyme`.

::::{tab-set}

:::{tab-item} AsCas12a (default)
```bash
crisprware score_guides \
  -b $PROJ/hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed \
  -i $PROJ/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots \
  --ucscgb_enzyme ascas12a \
  --cas12a_scorer deepcpf1 seq_deepcpf1variants \
  --cas12a_variant AsCas12a \
  --ucscgb $PROJ/track_hg38_ascas12a \
  --chrom_sizes $PROJ/hg38.chrom.sizes \
  --crispr_ots_bin $BIN \
  --ucscgb_scanner gpu \
  --ucscgb_cfd_threshold 0.023 \
  --ucscgb_list_cap 1000 \
  --ucscgb_blank_threshold 0 \
  --mismatches 4 \
  --threads 8 \
  --ucscgb_chunk_max 3000000 \
  -o $PROJ
```
`ascas12a` is the default, so `--ucscgb_enzyme ascas12a` may be omitted. `-i` is the **`TTTV`** index.
**No `--ucscgb_genome`** is needed (it is ignored for AsCas12a): a TTTV index has no TTTT off-targets to
exclude. The two on-target columns are **DeepCpf1** and the **AsCas12a seq-DeepCpf1 variant**
(`--cas12a_scorer deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a`).
:::

:::{tab-item} enCas12a
```bash
crisprware score_guides \
  -b $PROJ/hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed \
  -i $PROJ/encas12a/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots \
  --ucscgb_enzyme encas12a \
  --cas12a_scorer enpam_gb seq_deepcpf1variants \
  --cas12a_variant enAsCas12a-HF1 \
  --ucscgb $PROJ/track_hg38_encas12a \
  --chrom_sizes $PROJ/hg38.chrom.sizes \
  --ucscgb_genome $PROJ/hg38_primary_chr.fna \
  --crispr_ots_bin $BIN \
  --ucscgb_scanner gpu \
  --ucscgb_cfd_threshold 0.023 \
  --ucscgb_list_cap 1000 \
  --ucscgb_blank_threshold 0 \
  --mismatches 4 \
  --threads 8 \
  --ucscgb_chunk_max 3000000 \
  -o $PROJ
```
`-i` is the **`fa>0.3`** index. enCas12a **requires** `--ucscgb_genome` (the FASTA the index was built
from): each off-target's PAM is reconstructed from it to weight the four guidescan specificities. The two
on-target columns are **EnPAM-GB** and the **seq-DeepCpf1 enAsCas12a-HF1** variant
(`--cas12a_scorer enpam_gb seq_deepcpf1variants --cas12a_variant enAsCas12a-HF1`).
:::

::::

## Flags, grouped

**Inputs**
: `-b` guide BED (from step 2). `-i` the off-target index (from step 1), given as the prefix without the
  `.crot` extension. `--chrom_sizes` for `bedToBigBed`. `--crispr_ots_bin` the GPU engine.

**Enzyme**
: `--ucscgb_enzyme {ascas12a, encas12a, dual}` picks the track. `ascas12a` (default) scores AsCas12a
  (2xNLS) against a `TTTV` index (one specificity, no TTTN/enCas12a columns). `encas12a` scores enCas12a
  against the broad `fa>0.3` index (four guidescan specificities) and requires `--ucscgb_genome`. `dual`
  is the legacy both-nucleases track.

**On-target scores** (two columns in the bigBed)
: **AsCas12a** -- `--cas12a_scorer deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a` produces
  **DeepCpf1** (Kim 2018) and the **AsCas12a-DeepCpf1** seq-variant (Chen 2025). **enCas12a** --
  `--cas12a_scorer enpam_gb seq_deepcpf1variants --cas12a_variant enAsCas12a-HF1` produces **EnPAM-GB**
  and the **enAsCas12a-HF1** seq-variant.

**Off-target scoring**
: `--ucscgb` is the switch that emits the UCSC track (its value is the output track dir). The pass runs
  the engine with **one CFD matrix** for the chosen enzyme (`cfd-cas12a:2xnls` for AsCas12a,
  `cfd-cas12a:encas12a` for enCas12a). `--mismatches 4` is the search radius. `--ucscgb_genome` (the
  FASTA the index was built from) is **required for enCas12a**: each off-target's reference PAM is read
  to weight the four guidescan specificities `1 / (on_count + sum of off-target CFD [x PAM-activity
  weight])`, at the `fa>0.6` and `fa>0.3` activity cutoffs, raw and PAM-weighted (the `fa>0.3`
  PAM-weighted column drives the score and color). For enCas12a the scan is run **unfloored**
  internally, so those specificity sums are exact; `--ucscgb_cfd_threshold` then floors only the
  displayed popup list. `--ucscgb_genome` is **ignored for AsCas12a** (a TTTV index needs no PAM lookup);
  a `.fai` is built automatically if missing.

**Off-target list shaping**
: `--ucscgb_cfd_threshold 0.023` floors which off-targets are *listed* (counts are unfloored).
  `--ucscgb_list_cap 1000` caps stored off-targets per guide. `--ucscgb_blank_threshold 0` disables the
  "too many off-targets, blank the list" behavior.

**Scale / hardware**
: `--ucscgb_scanner gpu` uses the GPU engine. `--threads 8` is the CPU thread count **per chunk**.
  `--ucscgb_chunk_max 3000000` turns on chunking (see below).

## What `--ucscgb_chunk_max` does

Without it, the whole genome is one pass: simplest, but peak memory and disk scale with the full guide
set, and only one GPU is used. With `--ucscgb_chunk_max N` the guide BED is split into `<= N`-guide
chunks, each scored as an isolated single-pass `--ucscgb` run, then the pieces are merged (recomputing
`crisprDetails.tab` byte offsets):

```{mermaid}
flowchart LR
    G["guide BED<br/>135.9 M"] --> S["split into<br/>46 x 3 M chunks"]
    S --> W0["GPU 0"]
    S --> W1["GPU 1"]
    S --> Wn["GPU n"]
    W0 --> M["merge<br/>(rebase offsets)"]
    W1 --> M
    Wn --> M
    M --> T["one track"]
    classDef s fill:#eaf2fb,stroke:#3a6ea5,color:#16314a;
    class G,S,W0,W1,Wn,M,T s;
```

- One persistent worker per GPU, each claiming a device from a shared queue, so chunks never share a GPU.
- GPUs are **auto-detected** from `CUDA_VISIBLE_DEVICES` (your SLURM allocation). Override with
  `--ucscgb_gpus 0,1,2`.
- Each chunk's large off-target file is deleted as soon as it assembles, so peak disk is bounded by the
  chunks in flight, not the whole genome.
- `--ucscgb_keep_chunks` keeps the per-chunk dirs for debugging.

```{note}
**Percentiles are genome-wide, including on the chunked path, on every enzyme.** Each chunk is scored in
isolation, but the merge accumulates a histogram of every chunk's cleavage-score values and re-ranks each
model's `...Perc` column, plus the efficiency colour field it drives, against the whole genome before
writing the final bigBed. No post-pass is needed, and `--ucscgb_chunk_max` does not change a published
percentile.

The rule is derived from the track's own autoSql rather than tabulated per enzyme
(`track_regrade_spec`), so a run that scored one cleavage model, or three, re-grades exactly the models
it has. Ranking is on the value **as the track displays it** (4 decimal places), not on the
full-precision float behind it, so the chunked and unchunked paths agree exactly. `score`, the raw
specificity columns and `_specColor` are never ranked at all: they are absolute functions of the raw CFD
specificity, and the off-target **counts** are exact and population-independent.

The legacy `dual` track is the exception; it packs percentile and raw into one string, has no re-grading
rule, and warns when a chunked build merges.
```

## Worked example (hg38, 7x A5500)

135.9 M guides, 46 chunks of 3 M, GPUs auto-detected from the allocation:

```text
46 chunks x ~2.1 h/chunk / 7 workers  ~=  17.3 h wall
  per chunk: on-target ~2 min | GPU off-target enumerate ~30-42 min | assembly (sort+bigBed) ~90 min
final merge (46 pieces):  72 min
peak job RSS ~157 GB | transient chunk disk ~1.6 TB (bounded)
```

The assembly sort, not the GPU scan, is the wall-time bottleneck (GPU duty cycle ~28%), so more GPUs help
the enumerate stage but not the assembly tail.
