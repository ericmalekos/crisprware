# The bigBed columns

`cas12a.bb` is `bed9+14` = **23 fields**, one row per guide, defined by `cas12aTargets.as`. A real row
(hg38, `chr1:1000166`):

```text
chrom        chr1
chromStart   1000166                    0-based start of the 27-nt PAM+spacer site
chromEnd     1000193                    chromStart + 27
name         (blank)
score        96                         drives the browser score filter (see below)
strand       +
thickStart   1000170                    spacer start (chromStart + 4); the thin 4 nt is the PAM
thickEnd     1000193
itemRgb      0,200,0                    display color (see below)
guideSeq     CCACACTCGCCCCAGCCAATCGA    23-nt spacer
pam          TTTC                       4-nt PAM
unique_TTTV  True                       no identical protospacer elsewhere at a TTTV PAM
unique_TTTN  True                       ... at a TTTN PAM
TTTV_enCas12a_specificity   90% (0.9377)
TTTN_enCas12a_specificity   85% (0.7001)
TTTV_AsCas12a_specificity   89% (0.9958)
TTTN_AsCas12a_specificity   87% (0.9628)
ascas12a_deepcpf1_score     96% (0.8212)
enseq_deepcpf1_score        61% (0.6454)
deepcpf1_score              57% (53.5321)
enpam_gb_score              19% (0.7867)
_mouseOver   EnCas12a Spec: 85, AsCas12a Spec: 87, EnCas12a-DeepCpf1: 61, ...
_offset      27614949                   byte offset into crisprDetails.tab (0 = no list)
```

## The eight score columns read `percentile (raw)`

- **percentile** ranked within the guide's chunk.
- **raw** is exact and global: specificity on 0-1 (`1.0` = no off-targets); DeepCpf1 is the
  exception, its native output is ~0-100.

Four **off-target specificity** columns = {enCas12a, AsCas12a} x {scanned vs the TTTV, TTTN index}. Four
**on-target efficiency** columns = AsCas12a-DeepCpf1, EnCas12a-DeepCpf1, DeepCpf1, EnPAM-GB.

## score and color

- **score** (field 5, 0-1000): `round(AsCas12a TTTN specificity_raw x 100)`, e.g. `0.9628 -> 96`. The
  browser's score filter is therefore a specificity filter.
- **itemRgb** (field 9): ramps on the **AsCas12a-DeepCpf1 percentile**, green `>75`, yellow `>50`, red
  `<=50`; overridden to grey for low specificity, then 1/2-mismatch, then non-unique. The row above is
  green because AsCas12a-DeepCpf1 = 96%.

```{warning}
Percentiles are per-chunk (see the [scoring step](../ucsc_cas12a/03_score_track.md)), so `score`,
`itemRgb`, and the `NN%` are ranked within a chunk, not genome-wide. The parenthesized **raw** values are
exact and comparable across the whole genome.
```
