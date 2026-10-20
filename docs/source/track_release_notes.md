# Track release notes

Every change that alters a **published** column value, colour or label, newest first.

The rule this file exists to enforce: if a push changes how a value is DERIVED, not
just how it is formatted, it goes here before the push, so it is read rather than
found by diffing two deployed tracks.

A change is "derivation" and belongs here if it moves a number for a guide whose
underlying scores did not change. Rounding a column that nothing ranks is formatting.
Rounding a column that a percentile ranks is derivation, because the percentile moves
with it.

## 2026-10-19

### Specificity columns display 3 decimals (formatting)

`AsCas12a_specificity` and the four `enCas12a_spec_*` columns round to 3 places
instead of 4, trailing zeros stripped but one kept, so an exact 1 reads `1.0`.

Display only, and safe to change, because **no percentile is derived from
specificity**: `score` (column 5) remains the x100 integer taken from the
full-precision value, the three colour fields still threshold on the full-precision
value, and `_offset` is untouched so existing details files stay valid. On sacCer3
this moved the displayed value for 29,946 of 590,425 guides (5.07%); nothing rounded
to zero. The displayed value now sits at most 0.0005 from the true one, i.e. 0.05 of
a `score` unit, well inside the tolerance `validate_track_bed` allows between the
two.

### Hover labels follow the Cas9 track, and cite the right model (labels)

`_mouseOver` changes shape, from

    AsCas12a Spec: 40, DeepCpf1: 88, AsCas12a-DeepCpf1: 69

to

    AsCas12a Spec. Score: 40, Kim 2018: 88%, Chen 2025: 69%

matching Cas9's `MIT Spec. Score: 48, Doench 2016: 18%, Moreno-Mateos: 72%`,
including the `%` suffix that marks the cleavage figures as percentiles.

**"DeepCpf1" was wrong, not merely verbose.** crisprware loads
`weights/Seq_deepCpf1_weights.h5`, the sequence-only model from Kim et al. 2018.
DeepCpf1 proper is the chromatin-accessibility variant and is not what runs. Tracks
shipped before this date carried the wrong model name in every row's hover and in the
autoSql. The fields are now cited rather than named, which avoids the problem.

The autoSql comments move to Cas9's `crispr.as` convention for the same reason:
`Efficiency: Kim et al. 2018 Score (raw)` and `(percentile)`, and
`AsCas12a_specificity` reads `AsCas12a aggregate CFD specificity`.

## 2026-10-18

### Cleavage score percentiles are re-ranked on the rounded value (DERIVATION)

Flagged late. This was applied in the 2026-10-18 push without a release note, and
was found by diffing. It is the reason this file exists.

`deepcpf1` became an integer and `ascas12aDeepcpf1` 2 decimals, to match each
model's real precision and the Cas9 track's integer score fields. The percentiles
were **recomputed from the rounded values** rather than carried over.

That was not asked for, and it is a derivation change:

- the displayed percentile moved for about 38.5% of guides, by up to 3 points;
- `itemRgb` moved for 0.79% of guides (0.46% by a within-chromosome estimate made
  before the run, 1.16% measured genome-wide on rn7 chr20 after it: the pre-run
  proxy understated it);
- the percentile **range** changed, from 0 to 100 down to 2 to 100 and 1 to 100,
  because with only 112 and 94 distinct values the lowest tie's mid-rank is no
  longer near zero.

It was kept because the alternative is worse: the percentile is defined as the rank
of the value the track displays, so rounding the display while keeping the old
percentiles would leave 68 of 112 displayed `deepcpf1` values carrying two or more
different percentiles. `validate_track_bed` rejects exactly that, and the property is
what makes a percentile reproducible from the published file and what keeps the
chunked and unchunked build paths in agreement.

The lesson recorded here: a rounding request on a ranked column is never
display-only, and the consequence belongs in this file at the time, not in a diff
afterwards.

## 2026-10-13

### The 0.20 to 0.50 specificity grey is retired (colour)

Guides with specificity in 0.20 to 0.50 were greyed `90,90,90`, which hid the
efficiency ramp of a guide whose specificity is poor but still usable. They now show
the ramp like any other specific guide, and their specificity is read from `score`,
`_specColor` and the hover. The three real greys are unchanged: `150,150,150`
non-unique, `120,120,120` close off-target, `50,50,50` specificity below 0.20.

Retrofitted onto the finished tracks from their stored percentiles, with no rescan.
