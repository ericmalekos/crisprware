# Outputs

## What you get

```text
track_hg38_cas12a/
  cas12a.bb                  11.7 GB   bigBed, 135,867,611 rows x 23 fields
  crisprDetails.tab.gz       78.7 GB   bgzip off-target lists (243.5 GB uncompressed)
  crisprDetails.tab.gz.gzi   60 MB     bgzip index
```

For the hg38 example, **106,017,484 of 135,867,611 guides (78%)** carry a stored off-target list
(`_offset > 0`); the other 22% store none (`_offset = 0`): non-unique guides, or unique guides already
perfectly specific. Uniqueness itself is read from the `unique_TTTV` / `unique_TTTN` columns, see
[Reading the track](../reading_the_track/01_offtargets.md).

## Reading a details record

`crisprDetails.tab` is two columns: per-bucket mismatch **counts**, then the `|`-separated off-target
list. Counts are raw absolute totals (`[0-mm, 1-mm, 2-mm, 3-mm, 4-mm]`), unaffected by chunking.

```text
_mismatchCounts   _crisprOfftargets
0,2,0,0,2         chr2;182255206-;466;117|chr15;101948968-;344;94|chr19;83578+;344;94|...
```

Each off-target is `chrom;pos±;encas12a_cfd;ascas12a_cfd` (CFD scores x1000). A guide's bigBed row points
at one such line by uncompressed byte offset (`_offset`), which the browser resolves through the `.gzi`.
