# Color scheme mappings

The AsCas12a track ships **three** colors per guide, as the UCSC Cas9 track does. The browser shows
one at a time, chosen from the **Color by:** dropdown that the `colorFields` trackDb setting creates.
All three are pre-computed strings in the bigBed, so switching schemes redraws instantly and needs no
browser-side code.

| field | position | Cas9 counterpart | what it encodes |
|---|---|---|---|
| `itemRgb` | 9 | `reserved` | **Default.** Specificity greys, else the DeepCpf1 efficiency ramp |
| `_crisprDeepCpf1Color` | 10 | `_crisprScanColor` | Same scheme, ramped on the **second** on-target model |
| `_specColor` | 11 | `_specColor` | Specificity **alone**, no greys |

A field whose name starts with `_` is hidden from the details page by UCSC convention, which is why
the two alternative colors do not show up as rows in the popup.

## itemRgb: the blended default

Evaluated in this order, and **later rules win**. This is the Cas9 scheme: a guide that is hard to
target is greyed no matter how efficiently it is predicted to cut, so the efficiency ramp is only
ever visible on guides that are worth using.

| # | condition | color | swatch | meaning |
|---|---|---|---|---|
| 1 | DeepCpf1 percentile `<= 30` | `255,100,100` | <span style="color:#FF6464">red</span> | low predicted cleavage |
| 2 | DeepCpf1 percentile `30` to `55` | `255,255,0` | yellow | medium predicted cleavage |
| 3 | DeepCpf1 percentile `> 55` | `0,200,0` | green | high predicted cleavage |
| 4 | specificity `< 0.20` | `50,50,50` | dark grey | low specificity |
| 5 | `mm1 > 0` or `mm2 > 0` | `120,120,120` | mid grey | has a 1-mm or 2-mm off-target |
| 6 | not unique (`mm0 > 0`) or dropped | `150,150,150` | light grey | a perfect match exists elsewhere |

There is **no grey for a middle specificity band.** A guide between `0.20` and `0.50` keeps its
efficiency ramp colour: its specificity is poor but usable, and greying it hid the efficiency of a
guide somebody might still want. Read that guide's specificity from `score`, from `_specColor`, or
from the hover, all three of which agree with the default view that such a guide is usable:
`_specColor`'s red cut is `0.20`, not `0.50`, for exactly this reason. An earlier version of these
tracks greyed that band `90,90,90`; the colour is a function of the stored percentile, so existing
tracks were repainted in place without rescoring.

The ramp percentiles are ranked over **every guide in the genome**. A chunked build ranks within each
chunk first, then the merge re-ranks across all chunks before writing the bigBed, so the colors are
comparable between any two guides in the track. See [genome-wide percentiles](#genome-wide-percentiles).

Two conventions in that table are worth stating explicitly, because both differ from a naive reading:

- **The specificity cutoffs are absolute, not percentile.** They are applied to the raw 0 to 1
  AsCas12a specificity. A percentile would be useless here: in a genome of mostly unique guides, such
  as yeast, roughly 93% of guides would fall in the "bad" tail despite a median specificity of 1.0.
- **A missing value is treated as benign.** A `NaN` specificity counts as fully specific (not greyed)
  and a `NaN` efficiency percentile counts as 100 (green). In practice neither arises in these runs:
  both on-target models score every guide, and specificity is blank only for non-unique guides, which
  rule 7 has already greyed.

## _crisprDeepCpf1Color: the second on-target model

Identical logic to `itemRgb`, including every grey tier, with one substitution: the ramp runs on the
percentile of the **other** on-target model. The default ramps on DeepCpf1 (Kim 2018); this one ramps
on the AsCas12a seq-DeepCpf1 variant (Chen 2025).

The field exists because the two models disagree far more than the marginal color counts suggest. On
sacCer3, the two columns assign a **different color to 171,858 of 590,425 guides, 29.1%**, even
though each color's genome-wide total barely moves between them. Comparing the totals rather than
pairing per guide hides this almost completely.

Every disagreement sits inside the efficiency ramp; none involves a grey, since the greys are driven
by specificity and uniqueness and so are identical in both columns:

| default `itemRgb` | `_crisprDeepCpf1Color` | guides |
|---|---|---|
| green | yellow | 45,073 |
| yellow | green | 46,899 |
| red | yellow | 32,902 |
| yellow | red | 30,913 |
| green | red | 9,105 |
| red | green | 6,966 |

16,071 guides, 2.7%, invert outright between red and green. Cas9 fills this slot with Moreno-Mateos
and calls it `_crisprScanColor`; keeping that name here would have labelled a Cas12a model after a
Cas9 one, so the field is named for what it holds.

## _specColor: specificity alone

Purely a function of specificity. **No grey tier at all**, so a guide that `itemRgb` has already
greyed still reports its specificity here. That is the point of the column: the default scheme
collapses every kind of "hard to target" into grey, and this one recovers how bad the specificity
actually is.

| condition | color | swatch |
|---|---|---|
| specificity x100 `<= 20` | `255,0,0` | red |
| specificity x100 `20` to `70` | `128,128,0` | olive |
| specificity x100 `> 70` | `0,255,0` | green |

The red cut is `20`, the same number as the `itemRgb` grey, and **not** Cas9's `50`.
Carrying both meant the track drew two different "bad specificity" lines on one
quantity: a guide at `0.35` kept its efficiency ramp in the default view, which says it
is not hard to target, and read red here, which says it is. `0.20` is the line the
default scheme, the `score` filter and the description page already use, so this column
moved to it rather than the other way round.

Non-unique and dropped guides are forced to the worst tier (red) rather than left uncolored, since a
guide with a perfect match elsewhere has no meaningful specificity to report.

## Alignment with the Cas9 track

| aspect | Cas9 | AsCas12a | aligned? |
|---|---|---|---|
| efficiency ramp cuts | 30 / 55 percentile | 30 / 55 percentile | yes, but see the basis below |
| efficiency percentile basis | rank against a **fixed** mammalian reference distribution | rank **within this genome** | no, unavoidable |
| ramp colors | `255,100,100` / `255,255,0` / `0,200,0` | same | yes |
| non-unique grey | `150,150,150` (`#969696`) | `150,150,150` | yes |
| close-off-target grey | `120,120,120`, when alignment stopped (too many off-targets) | `120,120,120`, when a 1-mm or 2-mm off-target exists | same color, different trigger |
| low-specificity grey | one tier, `80,80,80` (`#505050`) | one tier, `50,50,50` | same idea, different shade |
| low-spec grey basis | MIT specificity `< 50` (grey to 49, never at 50) | raw specificity `< 0.20` | no, this track greys a narrower band |
| `_specColor` cuts | `<= 50` / `<= 70` / `> 70` | `<= 20` / `<= 70` / `> 70` | no, the red cut follows this track's grey |
| column 5 / score filter | MIT specificity, 0 to 100 | AsCas12a specificity x100 | yes |
| missing-efficiency color | `0,0,100` (`#000064`) | none | no, see below |

Five divergences are real:

1. **The percentile basis differs, and it is the one that changes how the colors read.** Cas9's
   efficiency percentile ranks a guide against a fixed reference distribution of mammalian guides, so
   its percentile is an absolute claim about predicted activity and the same raw score maps to the
   same percentile in every genome. Verified directly: across the hg38 and mm39 Cas9 tracks, all 87
   overlapping raw Doench values map to an identical percentile (raw 50 to the 50th, raw 60 to the
   81st, raw 70 to the 98th). This track ranks within the genome being displayed, because no
   comparable published reference distribution exists for Cas12a. So the shared 30 and 55 cut points
   do **not** mark the same predicted activity here as they do on the Cas9 track: they mark positions
   in this genome's own distribution. Closing this would mean choosing and shipping a Cas12a
   reference distribution, which is a scientific decision rather than a code change.
2. **A narrower low-specificity grey band.** Cas9 greys everything at or below MIT 50. This track
   greys only below 0.20, so a guide between 0.20 and 0.50 still shows its efficiency ramp here
   where Cas9 would hide it. That band is poor but usable, and `score` and `_specColor` already
   report the specificity.
3. **No "efficiency could not be calculated" color.** Cas9 reserves `0,0,100` for guides its
   efficiency model could not score. Both Cas12a models score every guide in these runs (0 blanks
   across sacCer3), so the case does not arise and no color is assigned to it.
4. **The second color field is renamed.** `_crisprScanColor` to `_crisprDeepCpf1Color`, as above.
5. **`_specColor`'s red cut is `20`, not `50`.** Cas9 can keep `50` because its default
   scheme greys everything at or below MIT 50, so its two schemes never disagree about
   which guides are bad. This track greys only below `0.20`, so a `50` red cut here
   contradicted the default view across the whole `0.20` to `0.50` band. The cut follows
   `SPEC_GREY_BAD`.

Cas9's own `_specColor` cuts are **undocumented but verified**: its track description page
documents only the `reserved` greys and ramp, and says nothing about `_specColor` or
`_crisprScanColor`. They were derived from real Cas9 rows and then checked against 165,992 scored
hg38 guides, where they hold with **zero** exceptions (no green at specificity <= 70, no olive outside
50 to 70, no red above 50), which is what establishes that Cas9's red cut really is `50` and so what
this track is diverging from. The same sample pins Cas9's other boundaries exactly: `itemRgb` is grey up
to specificity 49 and never at 50 or above, red runs to percentile 30, yellow from 32 to 53 and green
from 56 up. Note also
that the Cas9 page's own wording leaves percentile exactly 55 in neither the yellow nor the green
band (`> 30 and < 55`, then `> 55`); this track puts 55 in yellow.

## Genome-wide percentiles

Every percentile in the track, the `...Perc` columns, the hover numbers and the efficiency ramp that
`itemRgb` and `_crisprDeepCpf1Color` are built from, is ranked over **all guides in the genome**. The
numbers on this page are from the **AsCas12a** track, but the scheme is universal: every enzyme ships the
same three colourings on the same cuts, driven by whichever cleavage models it scored.

That takes an extra step on a large genome, because the build is chunked. `--ucscgb_chunk_max` splits
the guide BED into slices of roughly 3 million guides and scores each as an isolated run, so each
slice can only rank percentiles within itself. The guide BED is position-ordered, so a slice is a
contiguous stretch of one or two chromosomes: left alone, a percentile would describe a guide against
its regional neighbours, and the same raw score would land on different percentiles, and sometimes
different ramp colors, in different parts of the genome. On rn7 that drift reached 5 percentile
points for one model.

So `merge_ucscgb_tracks` re-ranks everything across all chunks before building the bigBed. The raw
score columns are exact per guide, which is what makes this recoverable with no rescan and no
re-scoring: the merge histograms the raw columns while it is already streaming each chunk's BED, then
rewrites the percentiles, the hover and the ramp colors in one further pass.

Two things are deliberately left alone, because neither is a rank:

- **`score`** (column 5) is the raw specificity times 100.
- **`_specColor`** uses absolute cutoffs on the raw specificity.

Grey cells are also preserved. The greys encode specificity and uniqueness, which `_item_rgb` applies
*after* the ramp and which are genome-wide to begin with, so a non-grey cell is always a pure ramp
color and is the only thing re-derived.

An unchunked build (`--ucscgb_chunk_max 0`, the default, and what the small genomes use) ranks over
every guide in one pass and skips the re-grade entirely. The two paths agree: `regrade_test` builds
sacCer3 both ways and asserts the bigBeds match on every percentile, color and hover.

## The hover specificity is the raw value, matching Cas9

The hover reports specificity as **raw x100**, identical to column 5, as the Cas9 track does with
`MIT Spec. Score: 48`. The two cleavage figures beside it *are* percentiles, and carry a `%` suffix
to say so:

```text
AsCas12a Spec. Score: 40, Kim 2018: 88%, Chen 2025: 69%
```

Specificity is deliberately never ranked. Its distribution is dominated by the tie at a perfectly
specific 1.0, so a rank reports the 53rd percentile for a guide whose score column says 100, and on
sacCer3 524,103 of 590,425 guides would read the same `53`. An earlier revision of this track did
rank it, which is what that number came from.

The cleavage models are **cited** rather than named because the obvious name is wrong: crisprware
loads `weights/Seq_deepCpf1_weights.h5`, the sequence-only model of Kim et al. 2018. DeepCpf1 proper
is the chromatin-accessibility variant and is not what runs. See
[Track release notes](../track_release_notes.md).

## Where this lives in the code

All of it is in `crisprware/ucsc_track.py`, and `build_track` is the single source of truth: there is
no post-hoc color converter, and a change here requires rebuilding the track.

| what | symbol |
|---|---|
| grey cutoff | `SPEC_GREY_BAD` (the 0.20 to 0.50 `SPEC_GREY_MARGINAL` tier is retired) |
| ramp cutoffs | `EFF_RAMP_LOW`, `EFF_RAMP_HIGH` |
| `_specColor` cutoffs | `SPEC_COLOR_BAD` (= `SPEC_GREY_BAD` x100), `SPEC_COLOR_MID` |
| ramp plus specificity grey | `_spec_eff_color` |
| full blended hierarchy | `_item_rgb` |
| specificity-only color | `_spec_color` |
| per-guide assignment | `build_track`, `enzyme="ascas12a"` branch |

The trackDb stanza that exposes the dropdown is `docs/trackhub/ascas12a.ra`, and the track
description page is `docs/trackhub/ascas12a.html`.
