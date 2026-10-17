# The bigBed columns

`score_guides --ucscgb` builds a single-enzyme Cas12a track, selected by `--ucscgb_enzyme`. Each is one
row per guide, defined by its own autoSql. The default is `ascas12a`:

| `--ucscgb_enzyme` | files | layout | fields |
|-------------------|-------|--------|--------|
| `ascas12a` (default) | `ascas12a.bb` / `ascas12aTargets.as` / `ascas12a_crisprDetails.tab.gz` | `bed9+12` | **21** |
| `encas12a` | `encas12a.bb` / `encas12aTargets.as` / `encas12a_crisprDetails.tab.gz` | `bed9+15` | **24** |

The leading bed9 columns, the two extra color fields right after `itemRgb`, and the
`_mouseOver`/`_offset` trailer are identical on both, matching the Cas9 track (see
[Color scheme mappings](02_colors.md)). Only the middle block (uniqueness + specificity + on-target
efficiency) differs, and it differs only in which specificity columns and which cleavage models the
enzyme carries.

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
deepcpf1                62                       DeepCpf1 (Kim 2018), raw (integer)
deepcpf1Perc            75                       DeepCpf1, genome-wide percentile
ascas12aDeepcpf1        0.73                     AsCas12a seq-DeepCpf1 variant (Chen 2025), raw (2 dp)
ascas12aDeepcpf1Perc    69                       seq-DeepCpf1 variant, genome-wide percentile
_mouseOver              AsCas12a Spec: 100, DeepCpf1: 75, AsCas12a-DeepCpf1: 69
_offset                 0                        byte offset into ascas12a_crisprDetails.tab (0 = no list)
```

Each on-target model contributes **two** fields, a raw value and a `...Perc` percentile, following the
Cas9 track's `fusi`/`fusiPerc` convention. Packing them into one `"75% (62)"` string reads fine in
the popup but leaves neither value filterable in the Table Browser or the API.

Each raw value is written at the precision its model supports, not at a uniform one: `deepcpf1` runs
roughly 1 to 115 and is an **integer**, while the 0-to-1 models (`ascas12aDeepcpf1`, `enpamGb`,
`encas12aHf1`) carry **2 decimals**. The percentile ranks that displayed value, so the two halves of a
pair always agree with each other.

```{note}
The `_mouseOver` specificity is the **raw value x100**, identical to `score`, as on the Cas9 track whose
hover reads `MIT Spec. Score: 80` with `score` 80. The two cleavage figures beside it *are* percentiles.
Specificity is deliberately not ranked: its distribution is dominated by the tie at a perfectly specific
1.0, so a rank would report the row above as the 53rd percentile while `score` says 100. See
[Color scheme mappings](02_colors.md).
```

## enCas12a track (`--ucscgb_enzyme encas12a`)

`bed9+15` = **24 fields**, the same shape as the AsCas12a track. Off-targets are searched against the
broad **fa>0.3 (30-PAM)** index, and the specificity is reported four ways (raw and PAM-activity-weighted,
at the fa>0.6 and fa>0.3 PAM-activity cutoffs). A real row:

```text
chrom                   chr1
chromStart              1000166
chromEnd                1000193
name                    (blank)
score                   98                       round(enCas12a fa>0.3 PAM-weighted specificity x 100)
strand                  +
thickStart              1000170
thickEnd                1000193
itemRgb                 0,200,0                  specificity grey, else the EnPAM-GB efficiency ramp
_crisprDeepCpf1Color    0,200,0                  same scheme on the second model (enAsCas12a-HF1)
_specColor              0,255,0                  specificity alone, no uniqueness greys
guideSeq                CCACACTCGCCCCAGCCAATCGA
pam                     TTTC
unique                  True                     no identical protospacer elsewhere at a fa>0.3 PAM
enCas12a_spec_fa06      0.9958                   RAW CFD specificity; off-targets at high-activity fa>0.6 PAMs
enCas12a_spec_fa06_wt   0.9971                   RAW PAM-activity-weighted specificity; fa>0.6 PAMs
enCas12a_spec_fa03      0.9762                   RAW CFD specificity; off-targets at moderate-activity fa>0.3 PAMs
enCas12a_spec_fa03_wt   0.9810                   RAW PAM-weighted specificity; fa>0.3 PAMs (drives score/color)
enpamGb                 0.79                     EnPAM-GB (Kim 2020), raw (2 dp)
enpamGbPerc             19                       EnPAM-GB, genome-wide percentile
encas12aHf1             0.65                     seq-DeepCpf1 enAsCas12a-HF1 variant (Chen 2025), raw (2 dp)
encas12aHf1Perc         61                       enAsCas12a-HF1, genome-wide percentile
_mouseOver              enCas12a Spec (fa>0.3 wt): 98, EnPAM-GB: 19, enAsCas12a-HF1: 61
_offset                 27614949                 byte offset into encas12a_crisprDetails.tab (0 = no list)
```

The four specificities are the guidescan specificity `1 / (on_count + sum CFD[x weight])` evaluated two
ways (raw, and each off-target's CFD scaled by its PAM's frac-active weight) at two PAM-activity cutoffs
(fa>0.6 and fa>0.3). All four are RAW, on 0-1. `enCas12a_spec_fa03_wt` is the operative value: it drives
`score`, the browser score filter, and the low-specificity greying.

## How the columns read

- **Specificity columns are RAW**, on 0-1 (`1.0` = perfectly specific, no off-targets). No percentile is
  shown: the vast majority of guides tie at `1.0`, so a percentile rank there is degenerate. AsCas12a has
  one such column (`AsCas12a_specificity`); enCas12a has four (`enCas12a_spec_fa06`, `_fa06_wt`, `_fa03`,
  `_fa03_wt`).
- **On-target efficiency columns** emit two columns per model, a raw value and a `...Perc` percentile
  sibling (`deepcpf1` / `deepcpf1Perc`, `enpamGb` / `enpamGbPerc`), following the Cas9 track's
  `fusi`/`fusiPerc` layout. This is the same on every enzyme. The percentile is ranked over **all guides
  in the genome**, and the raw value is exact (DeepCpf1's native output is ~0-100, the others are 0-1).

## score and color

**AsCas12a track**

- **score** (field 5, 0-1000): `round(AsCas12a-TTTV raw specificity x 100)`, e.g. `0.9958 -> 100`. The
  browser's score filter is therefore a specificity filter, on the **TTTV** aggregate (AsCas12a does not
  cleave a `TTTT` PAM).
- **itemRgb** (field 9): ramps on the **DeepCpf1 (Kim 2018) percentile**, using the Cas9 track's cuts of
  green `> 55`, yellow `> 30`, red `<= 30`, overridden (in this order) by grey for low **AsCas12a-TTTV**
  specificity (an **absolute** cutoff on the raw value: `< 0.20` -> `50,50,50`; a guide between
  0.20 and 0.50 is NOT greyed and keeps its ramp colour), then `120,120,120` for a TTTV
  1/2-mismatch off-target, then
  `150,150,150` non-unique. The row above is green because DeepCpf1 efficiency is high and the
  specificity (`1.0`) is well above `0.50`.
- **`_crisprDeepCpf1Color`** (field 10) and **`_specColor`** (field 11): two further pre-computed
  colorings of the same guide, selectable from the browser's **Color by:** dropdown. Full mapping in
  [Color scheme mappings](02_colors.md).

**enCas12a track**

- **score** (field 5, 0-1000): `round(enCas12a fa>0.3 PAM-weighted specificity x 100)` (the
  `enCas12a_spec_fa03_wt` column), e.g. `0.9810 -> 98`.
- **itemRgb** (field 9): ramps on the **EnPAM-GB percentile** with the same green/yellow/red thresholds,
  overridden by the same absolute `< 0.20` low-specificity grey (on `enCas12a_spec_fa03_wt`), then
  close-off-target grey, then non-unique. enCas12a **does** cleave `TTTT`, so the close-off-target counts
  are unrestricted. A guide with a perfect (0-mismatch) off-target elsewhere is **non-unique**: its four
  specificity columns are blanked, `unique = False`, and it is greyed `150,150,150`.

```{note}
**Percentiles are genome-wide on every track.** The `...Perc` columns and the efficiency ramps built from
them are re-ranked over all guides in the genome when the chunks merge, so they are comparable anywhere in
the genome regardless of how the build was chunked. See the
[scoring step](../ucsc_cas12a/03_score_track.md).

The **raw** specificity columns, the raw on-target values, and the specificity-derived `score` and
`_specColor` are absolute, never ranked, and comparable across the whole genome.

```
