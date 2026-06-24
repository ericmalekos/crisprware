# rank_guides

Filter scored guides, then select and rank the best few per target (gene, transcript, or BED interval).

```bash
crisprware rank_guides -k chrIII_sequence_scoredgRNA/chrIII_sequence_scoredgRNA.bed \
  -t tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf -f CDS \
  -c RS3_score_Chen2013 specificity_chrIII_sequence_gscan2 -m 0 0.2 -p 5 65 \
  -r RS3_score_Chen2013 --output_all
```

Keeps guides with RS3 >= 0 **and** specificity >= 0.2, in the 5th-65th CDS percentile, ranks them by
RS3, and returns the picks per gene -> `..._rankedgRNA/...bed`:

```text
#chr   start   stop    strand sequence             RS3_score_Chen2013  specificity_...  target_id  RS3_..._normalized  combined_weighted
chrIII 574162  574182  +      TAAGCTTAGTACTTATGAGG 0.9041              1.0              B0353.1    0.8819              0.8819
```

## Options you reach for

- **Targets**: `-t <gtf/bed>`; for GTF/GFF set `--target_mode gene|tx` and `-f/--feature CDS` (or exon,
  5UTR, ...). `-n/--number_of_guides` per target, `--min_spacing` between picks.
- **Filter** (`-c` columns paired with `-m` minimums): e.g. `-c RS3_score_Chen2013 specificity_idx
  -m 0 0.2`. `-p/--percentile_range lo hi` restricts position within the feature.
- **Rank** (`-r` columns, optional `-w/--column_weights`): combined into a weighted score;
  `--normalize_columns` scales each to 0-1 first.
- `--output_all` / `--plot_histogram` -- write the intermediate TSV (and histograms) at each filter
  stage, useful for tuning thresholds.

```{note}
`-c`/`-m` and `-r`/`-w` are positional pairs -- the order of columns must match the order of values.
```
