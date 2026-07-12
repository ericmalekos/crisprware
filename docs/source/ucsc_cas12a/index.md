# The UCSC Cas12a track pipeline

Four steps take a genome FASTA to a browser-ready track:

1. **Setup** the environment, the GPU engine binary, and the UCSC tools.
2. **Index + guides** build the off-target index (PAM `TTTN`) and enumerate every `TTTV` protospacer.
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

Every guide's `itemRgb` is decided by a fixed **priority cascade**: the first condition
that matches wins, so the greys (which flag a targeting problem) always override the
efficiency colors below them. Specificity is judged on the AsCas12a-**TTTV** aggregate
(it excludes TTTT off-targets, which AsCas12a cannot cleave) and its greys use an
**absolute** cutoff, not a within-population rank.

1. **Non-unique** (`150,150,150`, dark grey): an identical 23-nt protospacer exists
   elsewhere in the genome (a 0-mismatch duplicate), or the guide was dropped at the
   off-target pre-screen. Off-target analysis is skipped and the guide is not recommended.
2. **Close off-target** (`120,120,120`, grey): unique, but with one or more **TTTV**
   off-targets differing by only 1 or 2 mismatches across the 23-mer.
3. **Low AsCas12a-TTTV specificity**, an absolute cutoff on the raw 0-to-1 aggregate:
   - raw `< 0.20` -> `50,50,50` (very dark: low specificity),
   - `0.20 <= raw < 0.50` -> `90,90,90` (marginal specificity).
4. **Otherwise** (unique and specific) the color ramps on the on-target
   **AsCas12a-DeepCpf1** percentile:
   - `<= 50` -> red `255,100,100` (low predicted cleavage),
   - `> 50` and `<= 75` -> yellow `255,255,0` (medium),
   - `> 75` -> green `0,200,0` (high predicted cleavage).

Because the specificity greys are an absolute AsCas12a-TTTV cutoff (rule 3), a genome of
uniformly specific guides is colored by efficiency (red/yellow/green) instead of greyed
out. The browser **score** (0 to 1000) is that same AsCas12a-TTTV specificity times 100,
so the score filter doubles as a specificity filter. Only the efficiency-ramp percentile
in rule 4 is ranked within each processing chunk (relative); the specificity value and
the parenthesized raw scores are exact and comparable genome-wide. See
[the column reference](../reading_the_track/00_columns.md) for the full field list and
the hover text.
