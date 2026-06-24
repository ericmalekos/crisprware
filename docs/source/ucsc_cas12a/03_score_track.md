# Score and assemble the track

One command does everything: on-target scoring, off-target enumeration under both CFD matrices, and
assembly of the three browser files. With `--ucscgb_chunk_max` it also **chunks the genome, fans the
off-target passes across every allocated GPU, and merges the pieces** into a single track.

```bash
crisprware score_guides \
  -b $PROJ/hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed \
  -i $PROJ/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots \
  --cas12a_scorer deepcpf1 enpam_gb enseq_deepcpf1 seq_deepcpf1variants \
  --cas12a_variant AsCas12a \
  --ucscgb $PROJ/track_hg38_cas12a \
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

## Flags, grouped

**Inputs**
: `-b` guide BED (from step 2). `-i` the off-target index (from step 1). `--chrom_sizes` for `bedToBigBed`.
  `--crispr_ots_bin` the GPU engine.

**On-target scores** (four columns in the bigBed)
: `--cas12a_scorer deepcpf1 enpam_gb enseq_deepcpf1 seq_deepcpf1variants` with
  `--cas12a_variant AsCas12a` produces **DeepCpf1**, **EnPAM-GB**, **EnCas12a-DeepCpf1**, and
  **AsCas12a-DeepCpf1**.

**Off-target scoring**
: `--ucscgb` is the switch that emits the UCSC track (its value is the output track dir). The pass runs
  the engine with **two CFD matrices** (enCas12a + AsCas12a/2xNLS, the second is on by default), giving
  both nucleases' specificity. `--mismatches 4` is the search radius.

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

```{warning}
The bigBed **percentile** columns (score, color, the "NN%" display) are ranked **within each chunk**, the
merge does not recompute global percentiles. Raw specificity, raw on-target scores (the parenthesized
values), and the off-target **counts** are exact and population-independent. For genome-global
percentiles, re-rank the stored raw values in a post-pass; no rescoring is needed.
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
