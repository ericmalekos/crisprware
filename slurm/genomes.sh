#!/bin/bash
# Per-genome paths for the Cas12a track campaign, in ONE place.
#
# Sourced by build_index.sh, submit_shard_chain.sh and pipeline_shard.sh. It lives on
# its own because the alternative was a 20-line case block copied into each script, and
# a copied assumption is exactly how eight genomes came to be scanned against a TTTN
# index while three separate comments asserted the index was TTTV.
#
# Usage:
#   . "$(dirname "${BASH_SOURCE[0]}")/genomes.sh"
#   resolve_genome yeast          # sets DB W STEM SIZES CHUNK FNA IDX_TTTV
#
# PROJ must already be set by the caller (the campaign data dir).

# The off-target PAM every AsCas12a index must use. AsCas12a does not cleave a
# TTTT-prefixed protospacer, so a TTTT site in the index is a false off-target.
# crisprware/utils/crot_index.py enforces this by reading the index header; this is
# the value that satisfies it.
ASCAS12A_PAM=TTTV
ASCAS12A_PROTO_LEN=23
ASCAS12A_BIN_WIDTH=14

resolve_genome() {
  local org="$1"
  : "${PROJ:?resolve_genome needs PROJ set}"
  # CHUNK defaults are not preferences, they are a record of what is ON DISK. The merge
  # recomputes the expected chunk count as ceil(guides / chunk_max) and compares it
  # against what it finds, so split and every later phase must agree on chunk_max or
  # the merge refuses, and a shard index past the recomputed count is out of range.
  # Change these only together with a re-split.
  W=""; STEM=""; SIZES=""
  case "$org" in
    yeast)     DB=sacCer3;  CHUNK=0 ;;
    worm)      DB=ce11;     CHUNK=0 ;;
    fly)       DB=dm6;      CHUNK=0 ;;
    zebrafish) DB=danRer11; CHUNK=2000000 ;;
    mouse)     DB=mm39;     CHUNK=1000000 ;;
    rat)       DB=rn7;      CHUNK=1000000 ;;
    # hg38 and hs1 were built before this project settled on <db>/<db>_canon_*, so they
    # keep their own directories and stems. Pointing at them in place is deliberate: the
    # index alone is 7 GB, and a copy or a symlink farm would be one more thing that can
    # drift from the tree the earlier tracks were actually built against.
    human)     DB=hg38; CHUNK=1000000
               W=/private/groups/carpenterlab/emalekos/tmp/hg38-cas12a
               STEM=hg38_primary_chr; SIZES=$W/hg38.chrom.sizes ;;
    # DB=hs1 is the UCSC name; the directory still says t2t.
    t2t)       DB=hs1;  CHUNK=1000000
               W=/private/groups/carpenterlab/emalekos/tmp/t2t-cas12a
               STEM=t2t_chm13v2; SIZES=$W/t2t.chrom.sizes ;;
    *) echo "unknown ORG=$org (use yeast|worm|fly|zebrafish|mouse|rat|human|t2t)" >&2; return 1 ;;
  esac
  W=${W:-$PROJ/$DB}
  STEM=${STEM:-${DB}_canon}
  SIZES=${SIZES:-$W/${DB}.chrom.sizes}
  FNA=$W/${STEM}.fna
  # The TTTV index sits in its own tttv/ subtree rather than replacing the TTTN one:
  # enCas12a legitimately needs TTTN, and nothing should overwrite an index a job may
  # still be reading.
  IDX_TTTV=$W/tttv/${STEM}_crisprots/${STEM}_crisprots
  GUIDES=$W/${STEM}_gRNA/${STEM}_gRNA.bed
  export DB W STEM SIZES CHUNK FNA IDX_TTTV GUIDES
}

ALL_ORGS="yeast worm fly zebrafish mouse rat human t2t"
