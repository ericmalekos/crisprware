# Track release notes

Every change that alters a **published** column value, colour or label, newest first.

The rule this file exists to enforce: if a push changes how a value is DERIVED, not
just how it is formatted, it goes here before the push, so it is read rather than
found by diffing two deployed tracks.

A change is "derivation" and belongs here if it moves a number for a guide whose
underlying scores did not change. Rounding a column that nothing ranks is formatting.
Rounding a column that a percentile ranks is derivation, because the percentile moves
with it.

## 2026-10-22

### AsCas12a off-targets are searched on TTTV, not TTTN (DERIVATION)

Every AsCas12a track is rebuilt against an off-target index that enumerates **TTTV**.
The eight indexes used until now enumerated TTTN, which includes TTTT, and AsCas12a does
not cleave a TTTT-prefixed protospacer. Every TTTT site in those indexes was therefore a
false off-target.

The scanner excluded TTTT from exactly **one** of its nine outputs:

| scanner output | consumed as | excluded TTTT? |
|---|---|---|
| `tttv_sum` to `specificity_tttv` | `score` (column 5), `AsCas12a_specificity`, hover, `_specColor`, low-specificity greys | **yes** |
| `off_sum` to `specificity_tttn` | the TTTN column, correct for enCas12a | no, by design |
| `mm_counts` | mismatch buckets, the close-off-target grey, `_mismatchCounts` | no |
| the off-target list floor | `_crisprOfftargets` in the details popup | no |
| `on_count` | the 0-mm screen behind `dropped` and `unique` | no |
| `off_count` | `off_target_count` | no |
| `max_cfd` | max CFD | no |
| the saturation cap | `saturated`, and the per-guide count and CFD-sum caps | no |
| the CSV `specificity` column | `specificity_<index>` in every ordinary `score_guides` run | no |

So the published specificity was already right, and everything built on the mismatch
counts was not. The visible consequence is the `120,120,120` close-off-target grey,
which fires on `mm1 > 0 or mm2 > 0` and is not a CFD cutoff at all, so a TTTT site with
a high CFD greys a guide whose specificity reads near 1.0:

| track | guides | close-off-target grey | of those, specificity >= 0.20 |
|---|---|---|---|
| sacCer3 | 590,425 | 7,266 (1.23%) | 6,907 (1.17% of all guides) |
| ce11 | 6,617,735 | 770,060 (11.64%) | 656,095 (9.91%) |
| dm6 | 6,611,240 | 264,913 (4.01%) | 214,220 (3.24%) |
| danRer11 chr1 | 3,241,084 | 842,303 (25.99%) | 521,813 (16.10%) |
| rn7 chr1 | 10,469,943 | 1,688,044 (16.12%) | 1,117,752 (10.68%) |
| mm39 chr1 | 8,322,179 | 1,316,878 (15.82%) | 849,099 (10.20%) |
| hg38 chr1 | 10,225,880 | 1,898,802 (18.57%) | 1,276,773 (12.49%) |
| hs1 chr1 | 10,673,054 | 1,908,085 (17.88%) | 1,274,982 (11.95%) |

0.20 is the threshold the track itself calls hard to target (`SPEC_GREY_BAD`), so every
guide in the last column was greyed as hard to target while its own published
specificity said it was not.

One row, sacCer3 `chrXII:1061684-`, shows the whole mechanism. `specificity_tttn`
0.237351 against `specificity_tttv` 0.971510; the published `score` is 97, correctly
from the TTTV figure. Its details popup lists 8 off-targets at CFD 844, 703, 703, 703,
77, 77, 77 and 29, and reading those coordinates back out of the FASTA shows **7 of the
8 carry a TTTT PAM**, including the single 1-mismatch hit at CFD 0.844 that triggers the
grey. Its `_mismatchCounts` of `0,1,4,3,0` counts them all. The row therefore ships
`itemRgb` `120,120,120` (hard to target) beside `_specColor` `0,255,0` (green) and a
specificity of 0.972.

TTTT was not a quarter of the index. Rebuilding sacCer3 on TTTV takes it from 945,116
sites to 590,426, so **37.5%** of every TTTN index was TTTT, because poly-T runs are
common and concentrate in the soft-masked repeats. The `unique` flag was affected too,
but barely: 118 of sacCer3's 32,837 non-unique guides are unique among TTTV sites,
0.02% of all guides.

Direction of the correction: specificity does not move, and guides are **released** from
grey rather than newly greyed, so anyone who used the greys as a filter was discarding
usable guides. Nothing previously shown as usable becomes unusable.

The fix is at the index, not at each consumer, because eight of the nine outputs above
would each have needed their own guard and any one of them could be forgotten again:
`index_genome --pam TTTV` leaves no TTTT site for them to count. What makes that
durable is `crisprware/utils/crot_index.py`, which reads the enumerated PAM set out of
the `.crot` header, and `score_guides --ucscgb_enzyme ascas12a`, which now **refuses** a
TTTT-bearing index before scanning anything. Three separate comments in `ucsc_track.py`
had asserted "TTTV-only index" in prose while all eight genomes were indexed on TTTN;
the assumption is now read off the file. `tests/test_crot_index.py` covers it, including
a test that drives the real track entry point rather than the helper, because a
unit-tested helper nobody calls is how this shipped.

enCas12a is unaffected and deliberately keeps a TTTN index: it does cleave TTTT.

### `_specColor`'s red cut moves from 0.50 to 0.20 (colour)

`SPEC_COLOR_BAD` goes from 50 to 20, so the "Off-target specificity (AsCas12a)" view
calls a guide bad at the same specificity the default view greys at.

The track was drawing two different "bad specificity" lines on one quantity. A guide at
0.35 kept its efficiency ramp in the default scheme, which says it is not hard to
target, and read red in the specificity scheme, which says it is. 0.20 is the line
`SPEC_GREY_BAD`, the `score` filter and the description page already use, so this column
moved to it. The olive band is now 0.20 to 0.70; green above 0.70 is unchanged.

On sacCer3 this moves 866 of 590,425 guides (0.15%) from red to olive. The 32,837
non-unique guides are still forced to the worst tier, as on the Cas9 track.

This is now a deliberate divergence from Cas9, which keeps a 50 red cut because its
default scheme greys everything at or below MIT 50 and so has no disagreement to
resolve. It is pure recolouring and could have been retrofitted onto a finished track
from its stored specificity with no rescan, but it ships with the TTTV rebuild rather
than as a second push.

### All eight genomes rescanned on the fixed Cas12a CFD binary (DERIVATION)

Every off-target was rescored and every track rebuilt. For five genomes this moves
specificity, `score` (column 5), the low-specificity greys and `_specColor` for guides
whose underlying sequence did not change, so it is derivation in the strongest sense
this file recognises.

The cause: danRer11, rn7, mm39, hg38 and hs1 were scanned with a GPU binary dated
2026-06-05, which predates `ad3ecc6` ("fix Cas12a off-target CFD matrix application:
complement the d-base"). About a third of Cas12a mismatch types were left at the default
penalty of 1.0, i.e. no penalty, which overestimates off-target activity and understates
specificity. It is not patchable after the fact, because specificity is
`1 / (1 + sum of off-target CFDs)`.

The symptom that identified it, off-targets scoring exactly 1000 despite carrying at
least one mismatch, which cannot legitimately happen:

| genome | binary | off-targets at exactly 1000 |
|---|---|---|
| sacCer3 | 2026-08-09 | 158 / 119,269 = 0.13% |
| ce11 | 2026-08-09 | 2,047 / 3,339,879 = 0.06% |
| dm6 | 2026-08-09 | 4,271 / 2,406,062 = 0.18% |
| danRer11 | 2026-06-05 | 892,025 / 15,084,841 = 5.91% |
| rn7 | 2026-06-05 | 218,585 / 6,046,467 = 3.62% |
| mm39 | 2026-06-05 | 331,614 / 8,724,856 = 3.80% |

A 20 to 90-fold split falling exactly on the binary boundary. sacCer3, dm6 and ce11 had
used the post-fix CPU binary and were already correct; they were rescanned anyway so all
eight tracks come from one verified binary.

Specificity therefore RISES for guides in the five affected genomes, and no guide's
specificity falls. Anyone who filtered or ranked on the deployed danRer11 or rn7 tracks
was reading conservative values.

How the stale binary survived: five successive snapshots each copied the previous one's
binaries, so the discipline that froze the Python carried a 2026-06-05 Rust binary
forward five times. `slurm/make_snapshot.sh` now refuses to carry a binary over unless
`git diff <base>..HEAD -- crispr-ots/` is empty, and `slurm/pipeline_shard.sh` checks the
scanner binary's md5 against an allow and deny list before scoring anything.

Verification of the rebuilt tracks, all eight by `slurm/verify_track.sh`: bigBed
`itemCount` equal to the source guide BED exactly, 21-field autoSql with three colour
columns, zero off-target invariant violations, and one autoSql md5
(`71138e14cfd62e5004eb1721fccd85c4`) across the campaign. The sharded GPU build path was
separately validated bit-exact against an unsharded reference on sacCer3 at whole-genome
scale before the other seven genomes were run.

**Never pushed to hgwdev, and superseded before it was.** `slurm/push_ascas12a_tracks.sh`
was repointed from the pre-rescan `*_rounded` directories to `*_rescan` later the same
day, and then the TTTV index defect above was found. The `*_rescan` tracks carry the
fixed CFD binary but were still scanned against a TTTN index, so they must NOT be
pushed; the deployable build is the TTTV one.

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
