# Outputs

## What you get

Each track is three files, named for the enzyme. For the default **AsCas12a** build (`bed9+8`, 17
fields); sizes are from a reference hg38 build and are approximate (they scale with the field count):

```text
track_hg38_ascas12a/
  ascas12a.bb                        bigBed, 135,867,611 rows x 17 fields
  ascas12a_crisprDetails.tab.gz      bgzip off-target lists
  ascas12a_crisprDetails.tab.gz.gzi  bgzip index
  ascas12aTargets.as                 autoSql column schema
```

The **enCas12a** build (`--ucscgb_enzyme encas12a`) writes the same layout under `encas12a.*` names, with
`encas12a.bb` at `bed9+11` (**20 fields**) and `encas12a_crisprDetails.tab.gz` (+ `.gzi`) /
`encas12aTargets.as`. The extra columns are the four guidescan specificities (raw + PAM-weighted, at the
`fa>0.6` and `fa>0.3` cutoffs) instead of AsCas12a's single `AsCas12a_specificity`.

For the hg38 reference build, **106,017,484 of 135,867,611 guides (78%)** carry a stored off-target list
(`_offset > 0`); the other 22% store none (`_offset = 0`): non-unique guides, or unique guides already
perfectly specific. Uniqueness itself is read from the single **`unique`** column, see
[Reading the track](../reading_the_track/01_offtargets.md).

## Reading a details record

`*_crisprDetails.tab` is two columns: per-bucket mismatch **counts**, then the `|`-separated off-target
list. Counts are raw absolute totals (`[0-mm, 1-mm, 2-mm, 3-mm, 4-mm]`), unaffected by chunking.

```text
_mismatchCounts   _crisprOfftargets
0,2,0,0,2         chr2;182255206-;466|chr15;101948968-;344|chr19;83578+;344|...
```

Each per-enzyme track scans a **single** CFD matrix, so each off-target is `chrom;pos±;cfd` -- one CFD
score, x1000 (AsCas12a's `cfd-cas12a:2xnls`, or enCas12a's `cfd-cas12a:encas12a`). A guide's bigBed row
points at one such line by uncompressed byte offset (`_offset`), which the browser resolves through the
`.gzi`.

```{note}
The legacy `dual` track differs here: its `cas12a.bb` is `bed9+14` (23 fields) with `unique_TTTV` /
`unique_TTTN` flags, and each off-target in `crisprDetails.tab` carries **two** CFDs,
`chrom;pos±;enCas12a_cfd;asCas12a_cfd`, because that track scans both matrices in one pass.
```
