# score_guides

Add **on-target** activity and **off-target** specificity columns to a guide BED. On-target defaults to
RuleSet3 (`--tracr`); add deep models with `--cas9_scorer` / `--cas12a_scorer`. Off-target specificity
comes from each `-i` index (engine auto-detected).

```bash
crisprware score_guides -b chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed \
  -i chrIII_sequence_crisprots/chrIII_sequence_crisprots --tracr Chen2013 --threads 8
```

```text
#chr   start  stop   context                         strand sequence             RS3_score_Chen2013  specificity_chrIII_sequence_crisprots
chrIII 10569  10589  GCTGCCTACATGTACTTTTATTTGAGGGTC  +      CCTACATGTACTTTTATTTG -1.4272             1.0
chrIII 10590  10610  TTGAGGGTCCCCATGATCTTGAAGAGGAGA  +      GGGTCCCCATGATCTTGAAG -0.0695             1.0
```

One output column per scorer and per index. RuleSet3 is a z-score (centered at 0); specificity is 0-1
(1.0 = no off-targets).

## Options you reach for

- **On-target**: `--tracr Hsu2013|Chen2013|both` (RuleSet3); `--cas9_scorer deepspcas9 deephf_wt_u6 ...`
  or `--cas12a_scorer deepcpf1 enpam_gb ...` (`--cas12a_variant <name>` for the per-variant model).
  Multiple values allowed, one column each. See the on-target and off-target scoring pages for the menu.
- **Off-target**: one or more `-i` indices (crispr-ots and/or guidescan2); `--mismatches`, `--threshold`
  (discard guides with an off-target within this many mismatches), and -- guidescan2 only --
  `--rna_bulges` / `--dna_bulges`.
- **Speed / memory**: `--threads`, `--chunk_size` (guides held in memory at once), `--skip_rs3`,
  `--skip_gs2`, `--min_<scorer>` (drop low scorers before the off-target pass).
- **Genome-wide Cas12a track**: `--ucscgb <dir>` writes a UCSC Genome Browser Cas12a track instead of a
  scored BED (see the Cas12a UCSC tracks section).
