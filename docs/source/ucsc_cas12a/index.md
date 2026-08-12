# The UCSC Cas12a track pipeline

```{warning}
**Cas12a CFD correction (2026-07-18).** A strand-complement bug in the Cas12a CFD scorer made every
Cas12a specificity number produced before 2026-07-18 too low (up to ~23x). The worked-example
specificity/CFD values and `_mouseOver` numbers on these pages predate the fix and are illustrative
only; rebuild with a binary reporting `crispr-ots --version` 0.3.0+ to get correct values. See
`CAS12A_CFD_CORRECTION.md` at the repo root. SpCas9 is unaffected. (The per-enzyme track structure
below -- separate AsCas12a and enCas12a tracks -- is unchanged by the fix.)
```

One command, `score_guides --ucscgb`, takes a genome FASTA to a browser-ready Cas12a track. Pick the
nuclease with `--ucscgb_enzyme`:

- **`ascas12a`** (default) -- AsCas12a (WT / 2xNLS), scored against a **`TTTV`** off-target index, one
  specificity column.
- **`encas12a`** -- enCas12a (enAsCas12a), scored against the broad **`fa>0.3`** (30-PAM) off-target
  index, four specificity columns; **requires** `--ucscgb_genome`.
- **`dual`** (legacy) -- one track carrying both nucleases' TTTV/TTTN specificities. Kept for
  back-compat only.

Four steps take a genome FASTA to that track:

1. **Setup** the environment, the GPU engine binary, and the UCSC tools.
2. **Index + guides** build the off-target index (a `TTTV` index for AsCas12a, or the broad `fa>0.3`
   index for enCas12a) and enumerate every `TTTV` protospacer.
3. **Score + assemble** run the one `score_guides --ucscgb` command: on-target + off-target scoring,
   internally chunked and fanned across GPUs, merged into one track.
4. **Outputs** inspect the three files the build produces.

```{toctree}
:maxdepth: 1

00_overview
01_setup
02_index_and_guides
03_score_track
04_outputs
```

## How guide colors are assigned

Every guide's `itemRgb` is decided by a fixed **priority cascade**: the first condition that matches
wins, so the greys (which flag a targeting problem) always override the efficiency colors below them.
The cascade below describes the **AsCas12a** (default) track; the enCas12a differences are noted after.

1. **Non-unique** (`150,150,150`, dark grey): an identical 23-nt protospacer exists elsewhere in the
   genome (a 0-mismatch duplicate), or the guide was dropped at the off-target pre-screen. Off-target
   analysis is skipped and the guide is not recommended.
2. **Close off-target** (`120,120,120`, grey): unique, but with one or more off-targets differing by
   only 1 or 2 mismatches across the 23-mer. For AsCas12a these are **TTTV** off-targets: TTTT sites are
   excluded, because AsCas12a cannot cleave them.
3. **Low AsCas12a-TTTV specificity**, an absolute cutoff on the raw 0-to-1 CFD specificity:
   - raw `< 0.20` -> `50,50,50` (very dark: low specificity),
   - `0.20 <= raw < 0.50` -> `90,90,90` (marginal specificity).
4. **Otherwise** (unique and specific) the color ramps on the on-target **DeepCpf1 (Kim 2018)**
   percentile (falling back to the AsCas12a seq-DeepCpf1 variant when DeepCpf1 is not scored):
   - `<= 50` -> red `255,100,100` (low predicted cleavage),
   - `> 50` and `<= 75` -> yellow `255,255,0` (medium),
   - `> 75` -> green `0,200,0` (high predicted cleavage).

Because the specificity greys are an absolute cutoff (rule 3), not a within-population rank, a genome of
uniformly specific guides is colored by efficiency (red/yellow/green) instead of greyed out. The browser
**score** (0 to 100) is that same AsCas12a-TTTV specificity times 100, so the score filter doubles as a
specificity filter. Only the efficiency-ramp percentile in rule 4 is ranked within each processing chunk
(relative); the raw specificity value and the parenthesized raw scores are exact and comparable
genome-wide.

**enCas12a** uses the same cascade with two substitutions: rules 2-3 key off the enCas12a **fa>0.3
PAM-activity-weighted** specificity (enCas12a *does* cleave TTTT, so its close-off-target grey and score
are not TTTV-restricted), and rule 4 ramps on the **EnPAM-GB** percentile. See [the column
reference](../reading_the_track/00_columns.md) for the full field list and the hover text.
