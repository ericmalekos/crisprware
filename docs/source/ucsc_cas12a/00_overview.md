# Overview

`score_guides --ucscgb <outdir>` builds a Cas12a UCSC Genome Browser track in **one command**: it runs
the crispr-ots off-target pass, the on-target scoring, and the bigBed assembly internally (and, with
`--ucscgb_chunk_max` / `--ucscgb_gpus`, chunks the genome across multiple GPUs and merges the pieces).
Pick the nuclease with `--ucscgb_enzyme`:

- **`ascas12a`** (default) -- AsCas12a (WT / 2xNLS), scored against a **`TTTV`** off-target index.
- **`encas12a`** -- enCas12a (enAsCas12a), scored against the broad **`fa>0.3`** (30-PAM) off-target index.

Every track is single-enzyme: one CFD matrix, one off-target index, one specificity axis.

A UCSC CRISPR track is **three files**. For the default AsCas12a track:

| File | What it is |
|------|------------|
| `ascas12a.bb` | bigBed (`bed9+12`, **21 fields**) of on-target guides, one row per protospacer |
| `ascas12a_crisprDetails.tab.gz` | bgzip-compressed per-guide off-target lists |
| `ascas12a_crisprDetails.tab.gz.gzi` | the bgzip `.gzi` index (produced by `bgzip -i`) |

The enCas12a track writes the same three files under `encas12a.*` names, in the same shape
(`encas12a.bb` is `bed9+15`, **24 fields**: the extra columns are enCas12a's four specificities rather
than AsCas12a's one, and its own two cleavage models). The column schema -- `ascas12aTargets.as` or
`encas12aTargets.as` -- is an autoSql file written next to the bigBed and passed to `bedToBigBed` at
build time.

Track shape is universal across enzymes: bed9, then the same three colourings (`itemRgb`,
`_crisprDeepCpf1Color`, `_specColor`), then `guideSeq`/`pam`/`unique`, then that enzyme's specificity
columns, then a raw plus a `...Perc` column for each cleavage model it scored, then
`_mouseOver`/`_offset`. Which cleavage scores exist is the only thing that differs.

Each bigBed row stores a guide's coordinates, color, on-target scores, off-target specificity, and an
`_offset`: a byte offset into the uncompressed `*_crisprDetails.tab`. Clicking a guide makes the browser
seek that offset (via the `.gzi` index) and render the off-target list. `_offset = 0` means no list is
stored, the guide is either non-unique or already perfectly specific (see
[Off-target records](../reading_the_track/01_offtargets.md)).

## Cas12a specifics

- **PAM is 5', 4 nt** (`TTTV`), opposite end from Cas9's 3' `NGG`. `V = A/C/G`.
- A guide site spans **PAM(4) + spacer(23) = 27 nt**; the spacer is the bigBed "thick" region.
- Each track scores off-targets with **one CFD matrix**: `cfd-cas12a:2xnls` for AsCas12a,
  `cfd-cas12a:encas12a` for enCas12a.

## The two nucleases differ in what they cleave

The **guide** PAM is `TTTV` (the strong, canonical set) for both tracks. What differs is the **off-target
index**, because the two nucleases cleave different weak PAMs:

- **AsCas12a** does **not** cleave a `TTTT` PAM, so its off-target search uses a **`TTTV`** index. That is
  exactly the right search space, and no reference-PAM lookup is needed -- `--ucscgb_genome` is ignored
  for `ascas12a`. The track carries **one** specificity column, `AsCas12a_specificity`, the raw CFD
  aggregate on 0-1 (not a percentile).
- **enCas12a** cleaves a broad set of weak PAMs, so its off-target search uses the **`fa>0.3`** index: the
  30 four-mer PAMs enAsCas12a cleaves at >0.3 relative activity. enCas12a **requires** `--ucscgb_genome`
  -- each off-target's PAM is reconstructed from the reference to weight the specificity by that PAM's
  activity. The track carries **four** specificity columns: raw and PAM-activity-weighted, at the
  `fa>0.6` (high-activity) and `fa>0.3` (moderate) activity cutoffs. Each is a guidescan-style specificity,
  `1 / (on_count + sum of off-target CFD [x PAM-activity weight])`.

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
  seq-DeepCpf1 / seq-DeepCpf1variants [Chen et al. 2025](https://www.nature.com/articles/s41467-025-57150-9)
