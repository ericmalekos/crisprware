# Setup

Three things must be on hand: the `crisprware` Python env, the **GPU-built** off-target engine, and the
UCSC binaries (`bedToBigBed`, `bgzip`).

## Environment variables

Define these once; every later command reuses them.

```bash
export PROJ=/path/to/cas12a-track            # working directory for this build
export REPO=/path/to/crisprware              # cloned crisprware repo
export ENV=/path/to/conda_envs/crisprware    # env with crisprware + its deps
export BIN=$REPO/crispr-ots/target-gpu/release/crispr-ots   # GPU off-target engine

export PYTHONPATH=$REPO
export PATH=/path/to/conda_envs/ucsc/bin:$(dirname "$BIN"):$PATH   # bedToBigBed, bgzip, crispr-ots
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}  # nvrtc/libcuda for the GPU scanner
```

```{important}
The GPU path needs the engine built **with the `gpu` feature** (`target-gpu/release/crispr-ots`, a cudarc
build). A plain CPU build errors with *"requires a binary built with --features gpu"* the moment you pass
`--ucscgb_scanner gpu`. The `LD_LIBRARY_PATH` line is mandatory: the scanner JITs its CFD kernel via
`nvrtc` at runtime.
```

## The CLI

`crisprware` exposes the pipeline subcommands used below:

```bash
crisprware index_genome   ...   # build the off-target index
crisprware generate_guides ...  # enumerate protospacers
crisprware score_guides   ...   # score + (with --ucscgb) assemble the track
```

`crisprware <cmd>` is equivalent to `python -m crisprware.cli <cmd>`; the latter is handy when running
straight from a repo checkout without installing.

## Input genome

Start from a FASTA of just the sequences you want in the track. For the GRCh38 worked example we keep the
**24 canonical chromosomes** and drop `chrM` and all unplaced/unlocalized `GL*`/`KI*` scaffolds, then
build a `chrom.sizes` for `bedToBigBed`:

```bash
cd $PROJ
samtools faidx /path/to/GRCh38.primary_assembly.genome.fa \
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
  chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hg38_primary_chr.fna
samtools faidx hg38_primary_chr.fna
cut -f1,2 hg38_primary_chr.fna.fai > hg38.chrom.sizes
```

```text
chroms kept: 24    total bp: 3,088,269,832
```
