# crisprScore integration (optional, R)

Since v0.2 every on-target scorer is built into the Python package, so the Bioconductor
[crisprScore](https://bioconductor.org/packages/crisprScore/) R bridge is **optional** -- reach for it
only for methods not ported to Python: **RuleSet1, Azimuth, Lindel, CRISPRscan, CRISPRater**, and the
DeepHF T7 variants.

It is a packaged R wrapper, `crisprscore_multi.R`, run in a separate R/conda env that has the
crisprScore package installed:

```bash
conda activate <crisprscore_r_env>

crisprscore_multi.R input.tsv <method_numbers> output_scored.tsv <enzyme> <5p_flank> <3p_flank>
# RuleSet1 + Azimuth + CRISPRater on Cas9 guides:
crisprscore_multi.R input.tsv 1,2,11 scored.tsv Cas9 13 29
```

Input is a TSV with a `context` column (everything else is preserved); the script trims `context` to
each method's required length. Methods:

| # | Method | Enzyme |
|---|--------|--------|
| 1 | RuleSet1 | Cas9 |
| 2 | Azimuth | Cas9 |
| 3-8 | DeepHF variants (U6 / T7) | Cas9 |
| 9 | Lindel | Cas9 |
| 10 | CRISPRscan | Cas9 |
| 11 | CRISPRater | Cas9 |
| 12 | DeepSpCas9 | Cas9 |
| 13-14 | RuleSet3 (Hsu / Chen) | Cas9 |
| 15-16 | DeepCpf1 | Cas12a |
| 17 | EnPAMGB | Cas12a |

```{note}
Many of these (12, 13-17, DeepHF U6) are already native in `score_guides` and don't need R. Use
crisprScore mainly for 1, 2, 9, 10, 11.
```
