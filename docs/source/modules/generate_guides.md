# generate_guides

Find every protospacer for a PAM in the genome (optionally restricted to regions), with the context
window the on-target scorers need.

```bash
crisprware generate_guides -f tests/test_data/ce11/chrIII_sequence.fasta \
  -k tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf --feature CDS
```

Every NGG protospacer whose cut site lands in a CDS -> `chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed`
(312,925 guides):

```text
#chr   start  stop   id,sequence,pam,chromosome,position,sense              context                         strand
chrIII 1482   1502   chrIII:1483:+,GGAATGTACTTCTTCCCAAA,NGG,chrIII,1483,+   TGTTGGAATGTACTTCTTCCCAAAAGGTTC  +
```

## Options you reach for

- `-p/--pam NGG` (any IUPAC code, auto-expanded), `-l/--sgRNA_length 20`.
- `-w/--context_window <up> <down>` -- flanking context for scoring (`4 6` for RuleSet3; Cas12a uses
  `-w 4 3` with `--pam_5_prime`).
- **Where**: `-k/--locations_to_keep <bed/gtf>` + `--feature CDS` (keep guides cutting here),
  `--locations_to_discard` (drop guides cutting here), `--join_operation merge|intersect` for multiple
  keep files.
- **Quality filters**: `--gc_range lo hi`, `--discard_poly_T` (PolIII), `--discard_poly_G`,
  `--restriction_patterns GCGGCCGC ...` (with `--flank_5`/`--flank_3` vector context).
- `--pam_5_prime`, `-5/-3` active-site offsets -- Cas12a geometry.

```{tip}
The `context` column (not the bare spacer) is what `score_guides` feeds to the on-target models, so
generate with the window the scorer expects.
```
