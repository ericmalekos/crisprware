# Quickstart

The whole pipeline on the bundled *C. elegans* chrIII test data (`tests/test_data/ce11/`): build gene
models, index the genome, find guides in coding sequence, score them, and pick the best per gene. Each
step is one `crisprware` subcommand; the per-module pages cover the options.

## 1. Gene models

```bash
crisprware preprocess_annotation -g tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
  -m metagene consensus longest shortest
```

Collapses multi-isoform genes into four model GTFs (pick one for the steps below):

```text
chrIII_ce11.ncbiRefSeq/
  chrIII_ce11.ncbiRefSeq_consensus.gtf
  chrIII_ce11.ncbiRefSeq_longestCDS.gtf
  chrIII_ce11.ncbiRefSeq_meta.gtf
  chrIII_ce11.ncbiRefSeq_shortestCDS.gtf
```

## 2. Off-target index

```bash
crisprware index_genome -f tests/test_data/ce11/chrIII_sequence.fasta -p NGG -l 20
```

Builds the `crispr-ots` off-target index into `chrIII_sequence_crisprots/` (the `.crot` mmap database).

## 3. Generate guides

```bash
crisprware generate_guides -f tests/test_data/ce11/chrIII_sequence.fasta \
  -k tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf --feature CDS
```

Every NGG protospacer whose cut site lands in a CDS -> `chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed`
(312,925 guides):

```text
#chr   start  stop   id,sequence,pam,chromosome,position,sense              context                         strand
chrIII 1482   1502   chrIII:1483:+,GGAATGTACTTCTTCCCAAA,NGG,chrIII,1483,+   TGTTGGAATGTACTTCTTCCCAAAAGGTTC  +
```

## 4. Score guides

```bash
crisprware score_guides -b chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed \
  -i chrIII_sequence_crisprots/chrIII_sequence_crisprots --tracr Chen2013 --threads 8
```

Adds an RS3 on-target column and a `crispr-ots` off-target specificity column (the engine is
auto-detected from the index) -> `chrIII_sequence_scoredgRNA/...bed`:

```text
#chr   start  stop   context                         strand sequence             RS3_score_Chen2013  specificity_chrIII_sequence_crisprots
chrIII 10569  10589  GCTGCCTACATGTACTTTTATTTGAGGGTC  +      CCTACATGTACTTTTATTTG -1.4272             1.0
chrIII 10590  10610  TTGAGGGTCCCCATGATCTTGAAGAGGAGA  +      GGGTCCCCATGATCTTGAAG -0.0695             1.0
```

RS3 is a z-score (centered at 0); specificity is 0-1 (1.0 = no off-targets).

## 5. Rank guides

```bash
crisprware rank_guides -k chrIII_sequence_scoredgRNA/...bed \
  -t tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf -f CDS \
  -c RS3_score_Chen2013 specificity_chrIII_sequence_crisprots -m 0 0.2 -p 5 65 \
  -r RS3_score_Chen2013 --output_all
```

Keeps guides with RS3 >= 0, specificity >= 0.2, in the 5th-65th CDS percentile, ranks by RS3, and
returns the picks per gene -> `..._rankedgRNA/...bed`:

```text
#chr   start   stop    strand sequence             RS3_score_Chen2013  specificity_...  target_id  RS3_..._normalized  combined_weighted
chrIII 574162  574182  +      TAAGCTTAGTACTTATGAGG 0.9041              1.0              B0353.1    0.8819              0.8819
```
