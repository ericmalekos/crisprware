# The bigBed columns

`score_guides --ucscgb` builds one of three Cas12a tracks, selected by `--ucscgb_enzyme`. Each is one
row per guide, defined by its own autoSql. The default is `ascas12a`:

| `--ucscgb_enzyme` | files | layout | fields |
|-------------------|-------|--------|--------|
| `ascas12a` (default) | `ascas12a.bb` / `ascas12aTargets.as` / `ascas12a_crisprDetails.tab.gz` | `bed9+12` | **21** |
| `encas12a` | `encas12a.bb` / `encas12aTargets.as` / `encas12a_crisprDetails.tab.gz` | `bed9+11` | **20** |
| `dual` (legacy) | `cas12a.bb` / `cas12aTargets.as` / `crisprDetails.tab.gz` | `bed9+14` | 23 |

The leading bed9 columns and the `_mouseOver`/`_offset` trailer are identical across all three. The
AsCas12a track additionally carries two extra color fields right after `itemRgb`, matching the Cas9
track (see [Color scheme mappings](02_colors.md)). Only the middle block (uniqueness + specificity + on-target efficiency) differs.

## AsCas12a track (default)

`bed9+12` = **21 fields**. Off-targets are searched against a **TTTV-only** index (AsCas12a does not
cleave a `TTTT` PAM), so there is a single specificity column and no enCas12a/TTTN column. A real row
(sacCer3, `chrII:593863`):

```text
chrom                   chrII
chromStart              593863                   0-based start of the 27-nt PAM+spacer site
chromEnd                593890                   chromStart + 27
name                    (blank)
score                   100                      round(AsCas12a-TTTV raw specificity x 100)
strand                  +
thickStart              593867                   spacer start (chromStart + 4); the thin 4 nt is the PAM
thickEnd                593890
itemRgb                 0,200,0                  default color: specificity greys, else DeepCpf1 ramp
_crisprDeepCpf1Color    0,200,0                  same scheme on the seq-DeepCpf1 AsCas12a model
_specColor              0,255,0                  specificity alone, no greys
guideSeq                AATATGGCTTCCCTGAGAGCTTT  23-nt spacer
pam                     TTTA                     4-nt PAM
unique                  True                     no identical protospacer elsewhere at a TTTV PAM
AsCas12a_specificity    1.0                      RAW AsCas12a (WT/2xNLS) CFD specificity, 0-1
deepcpf1                61.8948                  DeepCpf1 (Kim 2018), raw
deepcpf1Perc            75                       DeepCpf1, within-genome percentile
ascas12aDeepcpf1        0.7324                   AsCas12a seq-DeepCpf1 variant (Chen 2025), raw
ascas12aDeepcpf1Perc    69                       seq-DeepCpf1 variant, within-genome percentile
_mouseOver              AsCas12a Spec: 53, DeepCpf1: 75, AsCas12a-DeepCpf1: 69   (all percentiles)
_offset                 0                        byte offset into ascas12a_crisprDetails.tab (0 = no list)
```

Each on-target model contributes **two** fields, a raw value and a `...Perc` percentile, following the
Cas9 track's `fusi`/`fusiPerc` convention. Packing them into one `"75% (61.8948)"` string reads fine in
the popup but leaves neither value filterable in the Table Browser or the API.

```{note}
The `_mouseOver` specificity is a **percentile**, while `score` is the raw value x100. The row above
reads `score 100` but `Spec: 53`, because a raw specificity of 1.0 is a large tie whose mid-rank is the
53rd percentile. See [Color scheme mappings](02_colors.md).
```

## enCas12a track (`--ucscgb_enzyme encas12a`)

`bed9+11` = **20 fields**. Off-targets are searched against the broad **fa>0.3 (30-PAM)** index, and the
specificity is reported four ways (raw and PAM-activity-weighted, at the fa>0.6 and fa>0.3 PAM-activity
cutoffs). A real row:

```text
chrom                   chr1
chromStart              1000166
chromEnd                1000193
name                    (blank)
score                   98                       round(enCas12a fa>0.3 PAM-weighted specificity x 100)
strand                  +
thickStart              1000170
thickEnd                1000193
itemRgb                 0,200,0
guideSeq                CCACACTCGCCCCAGCCAATCGA
pam                     TTTC
unique                  True                     no identical protospacer elsewhere at a fa>0.3 PAM
enCas12a_spec_fa06      0.9958                   RAW CFD specificity; off-targets at high-activity fa>0.6 PAMs
enCas12a_spec_fa06_wt   0.9971                   RAW PAM-activity-weighted specificity; fa>0.6 PAMs
enCas12a_spec_fa03      0.9762                   RAW CFD specificity; off-targets at moderate-activity fa>0.3 PAMs
enCas12a_spec_fa03_wt   0.9810                   RAW PAM-weighted specificity; fa>0.3 PAMs (drives score/color)
enpam_gb_score          19% (0.7867)             EnPAM-GB (Kim 2020)
encas12a_hf1_score      61% (0.6454)             seq-DeepCpf1 enAsCas12a-HF1 variant (Chen 2025)
_mouseOver              enCas12a Spec (fa>0.3 wt): 0.981, EnPAM-GB: 19, enAsCas12a-HF1: 61
_offset                 27614949                 byte offset into encas12a_crisprDetails.tab (0 = no list)
```

The four specificities are the guidescan specificity `1 / (on_count + sum CFD[x weight])` evaluated two
ways (raw, and each off-target's CFD scaled by its PAM's frac-active weight) at two PAM-activity cutoffs
(fa>0.6 and fa>0.3). All four are RAW, on 0-1. `enCas12a_spec_fa03_wt` is the operative value: it drives
`score`, the browser score filter, and the low-specificity greying.

The legacy `dual` track (`--ucscgb_enzyme dual`) keeps both enzymes in one row: `unique_TTTV`/`unique_TTTN`,
four `{enCas12a, AsCas12a} x {TTTV, TTTN}` specificity columns, and four on-target columns
(`ascas12a_deepcpf1_score`, `enseq_deepcpf1_score`, `deepcpf1_score`, `enpam_gb_score`). It is retained only
for back-compat; new tracks should use `ascas12a` or `encas12a`.

## How the columns read

- **Specificity columns are RAW**, on 0-1 (`1.0` = perfectly specific, no off-targets). No percentile is
  shown: the vast majority of guides tie at `1.0`, so a percentile rank there is degenerate. AsCas12a has
  one such column (`AsCas12a_specificity`); enCas12a has four (`enCas12a_spec_fa06`, `_fa06_wt`, `_fa03`,
  `_fa03_wt`).
- **On-target efficiency columns read `percentile (raw)`**, e.g. `96% (0.8212)`. The **percentile** is
  ranked within the guide's chunk; the **raw** value in parentheses is exact and global (DeepCpf1's native
  output is ~0-100, the others are 0-1).

## score and color

**AsCas12a track**

- **score** (field 5, 0-1000): `round(AsCas12a-TTTV raw specificity x 100)`, e.g. `0.9958 -> 100`. The
  browser's score filter is therefore a specificity filter, on the **TTTV** aggregate (AsCas12a does not
  cleave a `TTTT` PAM).
- **itemRgb** (field 9): ramps on the **DeepCpf1 (Kim 2018) percentile**, using the Cas9 track's cuts of
  green `> 55`, yellow `> 30`, red `<= 30`, overridden (in this order) by grey for low **AsCas12a-TTTV**
  specificity (an **absolute** two-tier cutoff on the raw value: `< 0.20` darkest `50,50,50`,
  `0.20-0.50` lighter `90,90,90`), then `120,120,120` for a TTTV 1/2-mismatch off-target, then
  `150,150,150` non-unique. The row above is green because DeepCpf1 efficiency is high and the
  specificity (`1.0`) is well above `0.50`.
- **`_crisprDeepCpf1Color`** (field 10) and **`_specColor`** (field 11): two further pre-computed
  colorings of the same guide, selectable from the browser's **Color by:** dropdown. Full mapping in
  [Color scheme mappings](02_colors.md).

**enCas12a track**

- **score** (field 5, 0-1000): `round(enCas12a fa>0.3 PAM-weighted specificity x 100)` (the
  `enCas12a_spec_fa03_wt` column), e.g. `0.9810 -> 98`.
- **itemRgb** (field 9): ramps on the **EnPAM-GB percentile** with the same green/yellow/red thresholds,
  overridden by the same absolute two-tier low-specificity greys (on `enCas12a_spec_fa03_wt`), then
  close-off-target grey, then non-unique. enCas12a **does** cleave `TTTT`, so the close-off-target counts
  are unrestricted. A guide with a perfect (0-mismatch) off-target elsewhere is **non-unique**: its four
  specificity columns are blanked, `unique = False`, and it is greyed `150,150,150`.

```{warning}
The itemRgb efficiency **ramp** (green/yellow/red) and the `NN%` percentile prefixes on the on-target
columns are ranked **within each chunk** (see the [scoring step](../ucsc_cas12a/03_score_track.md)), not
genome-wide. The **raw** specificity columns, the parenthesized raw on-target values, and the
specificity-derived `score` are exact and comparable across the whole genome.
```
