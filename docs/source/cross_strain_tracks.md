# Cross-strain specificity tracks

The [Cas12a track](ucsc_cas12a/index.md) answers "how specific is this guide in the reference?".
This one answers **"does this guide stay specific in every other assembly of the species?"** Each
reference guide is one bigBed feature; an external details file gives that guide's off-target
specificity in *every* genome in the panel.

Written around the mouse panel, **mm39/GRCm39 plus 19 strain assemblies** (the Cell Genomics mouse
strain set, which includes two T2T co-assemblies and one *Mus spretus* outgroup), producing two
tracks: SpCas9 and Cas12a. The recipe works for any set of assemblies of one species.

```{contents} On this page
:local:
:depth: 1
```

## What you get

| track | guide PAM | length | specificity | features | `.bb` | details |
|-------|-----------|--------|-------------|----------|-------|---------|
| SpCas9 | `NGG` (3') | 20 nt | NGG | 268,496,809 | 22.4 GB | 14.9 GB |
| Cas12a (AsCas12a) | `TTTV` (5') | 23 nt | AsCas12a, TTTV off-target | 108,414,811 | 8.9 GB | 3.8 GB |

Three files per track: `<track>.bb` (bigBed `bed9+8`), `<track>_details.tab.gz` (+ `.gzi`), and the
autoSql `<track>.as`. Same two-file pattern as `score_guides --ucscgb`: the bigBed row carries an
`_offset`, an uncompressed byte offset into the details file, which the browser resolves via the
`.gzi`.

The Cas12a track is **AsCas12a** only, one specificity per genome. Guides are the canonical `TTTV`
design PAM, and the off-target CFD is summed over the `TTTV` PAM space (`specificity_tttv`): WT
AsCas12a is effective only against TTTV, so TTTT-PAM off-targets are EXCLUDED. The index is built
`--pam TTTN` (the broader set, so the search is complete) but only TTTV off-targets are counted.
(Caveat: the multimapper flag is on the TTTN `dropped` basis, conservative for ~0.4% of cells.)

```{note}
**Cas12a CFD correction (2026-07-18).** A strand-complement bug in the Cas12a CFD scorer deflated
every Cas12a specificity produced before 2026-07-18 (fixed; regenerate with `crispr-ots` 0.3.0+).
**This cross-strain panel is AsCas12a only.** enCas12a's accepted PAMs are far broader than TTTN
(TTYN, VTTV, TRTV, ...) and are not a single IUPAC pattern, so no enCas12a arm was carried across the
N per-strain indexes in this panel's design. (The main single-reference pipeline *does* score an
enCas12a track: `score_guides --ucscgb --ucscgb_enzyme encas12a` runs it against a broad fa>0.3
(30-PAM) off-target index built with `index_genome --pam_list_file`.)
```

## The design: one superset, N indexes

Scoring each genome's own guides against its own index would give N unrelated tables that cannot be
joined by row, and would silently mishandle absence. Instead, build **one guide superset** (the
deduplicated union of the guides found in all N genomes) and score *that same list* against *every*
genome's index.

Two things fall out of it:

- **Presence is a bitmask, not a score.** A guide absent from a strain has no perfect match there, and
  the scorer substitutes a synthetic on-target rather than reporting nothing, so an absent guide still
  gets a plausible-looking specificity. Absence must be read from the superset's **origin mask**, never
  from the score.
- **The N outputs are row-aligned.** Every genome is scored on the same guide list in the same order,
  so the N aggregate files combine by column-wise `paste`. No join, no sort.

## Step 1 - index and guides, per genome

Index every genome (off-target search space), and find guides in every genome (so the superset is the
union, not just the reference's guides).

```bash
for G in $GENOMES; do
  # SpCas9
  crispr-ots index --pam NGG  -l 20 -w 11 -o idx/$G  fa/$G.fa
  crisprware generate_guides -f fa/$G.fa -p NGG  -l 20 -t 8 -o guides/$G

  # Cas12a: guides TTTV, index TTTN (broad, complete); AsCas12a off-target PAM counted = TTTV only
  crispr-ots index --pam TTTN -l 23 --pam-5-prime -w 11 -o idx_cas12a/$G fa/$G.fa
  crisprware generate_guides -f fa/$G.fa -p TTTV --pam_5_prime -l 23 -t 8 -o guides_cas12a/$G
done
```

## Step 2 - superset and origin mask

Reduce each genome's guide BED to a sorted unique sequence list, then k-way streaming merge into one
superset. Each guide gets a stable id (`sup_N`) and an N-bit hex **origin mask** (bit *i* set = the
guide exists in genome *i*; bit 0 = the reference, so the representative locus prefers it).

```text
superset_guides.kmers.csv   id,sequence,pam,chromosome,position,sense
superset_origin.tsv         sup_0   0fffff   20            # id, origin mask, n_genomes
origin_legend.tsv           0  mm39 / 1  129S1_SvImJ / ...
```

Mouse panel: **356,572,107** SpCas9 and **160,465,870** Cas12a superset guides.

## Step 3 - score the superset against every genome

One `crispr-ots enumerate` per (genome, shard), as a SLURM array. Aggregated output mode gives one row
per guide, and `--keep-dropped` keeps multimappers as rows so the output stays row-aligned.

```bash
crispr-ots enumerate --output-mode aggregated --scanner gpu \
  --score cfd \                    # cfd-cas12a:2xnls (WT AsCas12a) for the Cas12a arm
  --threshold 0 -m 3 --max-off-targets 201 --min-specificity 0.001 --keep-dropped \
  --kmers-file shards/shard_0000.csv -o out/$G/shard_0000.agg  idx/$G
```

Identical parameters for every genome, or the cross-strain comparison is not apples to apples. Each
row yields the two facts the track needs: `specificity_tttv` (the Cas12a value; the single `cfd`
specificity for Cas9), and `dropped` (1 = the guide has 2 or more perfect on-targets, a multimapper).

## Step 4 - cross-strain matrix

Read the kmers file, the origin mask, and the N aggregate files **in lockstep** (they are row-aligned)
and emit one sequence-keyed row per reference guide. Each genome's cell resolves in priority order:

| cell | condition |
|------|-----------|
| `ABSENT` | origin bit unset: the guide does not exist in that genome |
| `MULTI` | `dropped = 1`: 2 or more perfect on-targets there, unusable |
| float | the specificity: `specificity_tttv` (Cas12a, AsCas12a) or `cfd` (SpCas9) |

Keep only rows whose origin bit 0 is set (guides that exist in the reference), and sort by sequence
for the join back to reference coordinates.

## Step 5 - assemble the track

Join the matrix back to the reference guide BED (by sequence, so a multimapping guide correctly yields
one feature per locus), then transform each row into one compact bigBed line plus one details line.

```bash
bedToBigBed -sort -tab -type=bed9+8 -as=<track>.as -extraIndex=name \
  <track>.compact.bed mm39.canon.chrom.sizes <track>.bb
bgzip -i <track>_details.tab
```

- **bigBed `bed9+8`**: bed9 (feature = protospacer + PAM, `thick` = protospacer, `name` = protospacer,
  `itemRgb` = color) plus `spec_ref`, `nStrainsPresent`, `nStrainsSpecific`, `nStrainsFragile`,
  `minStrainSpec`, `meanStrainSpec`, `_mouseOver`, `_offset`.
- **Details line**: `<protospacer>` TAB `<genome=value;genome=value;...>`, genomes in a fixed order on
  every row, `value` a single float or `ABSENT` or `MULTI` (both tracks; the Cas12a value is the
  AsCas12a `specificity_tttv`). Line 0 is the header `#_seq` TAB `_strainSpec`, so real offsets start
  above 0 and `_offset = 0` remains the "no detail" sentinel.
- **Color and summary fields** are computed over a **comparison subset** of the panel, here the 16
  *M. musculus* strains, excluding the two T2T co-assemblies (redundant with their primary assemblies)
  and *M. spretus* (a different species, so its divergence would swamp the signal). The details file
  still reports all 20 genomes.
- **Not unique in the reference** (a multimapper in GRCm39 primary **or** in C57BL/6J-T2T, which
  resolves duplications GRCm39 collapses) forces grey `150,150,150`, score 0, and no details line.
  Filter them out with any `scoreFilter > 0`.

```{tip}
`bedToBigBed` needs `-sort`: extending a feature to protospacer+PAM can pull a minus-strand
`chromStart` earlier than its neighbor's, which breaks an otherwise position-sorted BED.
```

## Reading a guide's per-genome scores

```bash
bgzip -b <OFFSET> -s 8192 mm39_cas9_crossStrain_details.tab.gz | head -1
# CCCCAGCCCCATGTCGTTGCATA   mm39=1.0000;C57BL_6J_T2T=ABSENT;C57BL_6NJ=1.0000;129S1_SvImJ=0.8631;...
```

## Loading in UCSC

```text
track type=bigBed name="mm39 Cas9 crossStrain" db=mm39 itemRgb=on labelFields=none hideItemName=on
  bigDataUrl=https://.../mm39_cas9_crossStrain.bb
```

`labelFields=none` + `hideItemName=on` hide the `name` field (it is the protospacer, used as the click
key to disambiguate two opposite-strand guides in the same window, not a label). The per-strain table
on the click page is rendered from `_strainSpec` by the `extFieldStrainSpec` hgc handler.

## Worked example: 20 mouse genomes

- Panel: mm39 + 19 strain assemblies. Superset 356.6 M (SpCas9) / 160.5 M (Cas12a) guides, each scored
  against all 20 indexes (`enumerate`, GPU for Cas9, CPU for the Cas12a arm).
- Reference guides in the tracks: 268,496,809 (SpCas9) and 108,414,811 (Cas12a), on chr1-19 and chrX.
- Assembly of one Cas12a track (matrix, join, transform, bigBed) is about 3 h on 64 CPU cores.

```{note}
**CFD correction.** Cas12a specificity here uses a binary that fixes a strand-complement bug in the
Cas12a CFD scorer (pre-2026-07-18 numbers were too low; rebuild with `crispr-ots` 0.3.0+). This panel
produces only the AsCas12a arm: enCas12a's PAM set (TTYN, VTTV, TRTV, ...) is far broader than TTTN and
is not a single IUPAC pattern, so it was not carried across the N per-strain indexes here. The main
single-reference pipeline scores enCas12a separately via `score_guides --ucscgb --ucscgb_enzyme encas12a`
against a broad fa>0.3 (30-PAM) index.
```
