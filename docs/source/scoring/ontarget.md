# On-target scoring

On-target scorers predict cleavage activity. RuleSet3 is the default (`--tracr`); the deep models are
added with `--cas9_scorer` (SpCas9) or `--cas12a_scorer` (Cas12a). Each adds one column;
`--min_<scorer>` drops low guides before the off-target pass.

## RuleSet3

`--tracr Hsu2013 | Chen2013 | both` -- **SpCas9** (default). Output is a z-score; the two values are
tracrRNA variants. [DeWeirdt et al. 2022](https://www.nature.com/articles/s41467-022-33024-2) (the
[`rs3`](https://github.com/gpp-rnd/rs3) package).

## DeepSpCas9

`--cas9_scorer deepspcas9` -- **SpCas9**, CNN on-target activity.
[Kim et al. 2019, Sci. Adv.](https://www.science.org/doi/10.1126/sciadv.aax9249)

## DeepHF

`--cas9_scorer deephf_wt_u6 | deephf_esp | deephf_hf` -- **SpCas9** (wildtype/U6, eSpCas9-1.1,
SpCas9-HF1); BiLSTM + biofeatures.
[Wang et al. 2019, Nat. Commun.](https://www.nature.com/articles/s41467-019-12281-8)

## DeepCpf1

`--cas12a_scorer deepcpf1` -- **Cas12a** (AsCas12a / LbCas12a), the original seq-DeepCpf1 CNN.
[Kim et al. 2018, Nat. Biotechnol.](https://www.nature.com/articles/nbt.4061)

## enPAM+GB

`--cas12a_scorer enpam_gb` -- **enAsCas12a**; gradient-boosted on RuleSet2 features.
[DeWeirdt et al. 2021, Nat. Biotechnol.](https://www.nature.com/articles/s41587-020-0600-6)

## enseq-DeepCpf1

`--cas12a_scorer enseq_deepcpf1` -- **Cas12a** (modern AsCas12a), sequence-only CNN.
[Chen et al. 2025, Nat. Commun.](https://www.nature.com/articles/s41467-025-57150-9)

## seq-DeepCpf1variants

`--cas12a_scorer seq_deepcpf1variants --cas12a_variant <name>` -- per-variant CNNs for **23 Cas12a
variants**: AsCas12a, AsCas12a_Ultra (2xNLS), enAsCas12a-HF1, HyperFi-AsCas12a, and the LbCas12a /
FnCas12a / iCas12a / Lb2Cas12a / Eb / Ce families. Names are case- and separator-insensitive.
[Chen et al. 2025, Nat. Commun.](https://www.nature.com/articles/s41467-025-57150-9)
