# Overview

A UCSC CRISPR track is **three files**:

| File | What it is |
|------|------------|
| `cas12a.bb` | bigBed (`bed9+`) of on-target guides, one row per protospacer |
| `crisprDetails.tab.gz` | bgzip-compressed per-guide off-target lists |
| `crisprDetails.tab.gz.gzi` | the bgzip `.gzi` index (produced by `bgzip -i`) |

The column schema, `cas12aTargets.as`, is a fixed autoSql file kept **in the repo** and passed to
`bedToBigBed` at build time. It is not a per-genome output.

Each bigBed row stores a guide's coordinates, color, on-target scores, off-target specificity, and an
`_offset`: a byte offset into the uncompressed `crisprDetails.tab`. Clicking a guide makes the browser
seek that offset (via the `.gzi` index) and render the off-target list. `_offset = 0` means no list is
stored, the guide is either non-unique or already perfectly specific (see
[Off-target records](../reading_the_track/01_offtargets.md)).

## Cas12a specifics

- **PAM is 5', 4 nt** (`TTTV`), opposite end from Cas9's 3' `NGG`. `V = A/C/G`.
- A guide site spans **PAM(4) + spacer(23) = 27 nt**; the spacer is the bigBed "thick" region.
- We score with **two CFD off-target matrices** in one pass, enCas12a (WT) and AsCas12a/2xNLS.

## The off-target asymmetry that drives the index

The **guide** PAM is `TTTV` (the strong, canonical set). The **off-target index** uses `TTTN`, because
enCas12a also cleaves the weak `TTTT` PAM.

```{note}
Throughout, "specificity" is an aggregate over a guide's per-off-target CFD scores: 1.0 = perfectly
unique, lower = more / closer off-targets. On-target scores (DeepCpf1 and friends) are a separate axis
predicting cutting efficiency.
```

## References

- **Off-target CFD matrices** [DeWeirdt et al. 2021](https://www.nature.com/articles/s41587-020-0600-6)
  ([cas12a_manuscript code](https://github.com/PeterDeWeirdt/cas12a_manuscript))
- **On-target models** DeepCpf1 [Kim et al. 2018](https://www.nature.com/articles/nbt.4061);
  EnPAM-GB [DeWeirdt et al. 2021](https://www.nature.com/articles/s41587-020-0600-6);
  enseq-DeepCpf1 / seq-DeepCpf1variants [Chen et al. 2025](https://www.nature.com/articles/s41467-025-57150-9)
