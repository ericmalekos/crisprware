#!/usr/bin/env bash
# Push the re-graded AsCas12a (bed9+12, 21-field) tracks to hgwdev.
# Run from mustard inside screen/tmux. ~135 GB; --partial makes it resumable.
#
# Destination is ascas12a_tracks, NOT the existing cas12a_tracks: this is a new
# schema with new filenames, and the old deployment stays untouched.
#
# Usage:
#   push_ascas12a_tracks.sh            push everything
#   push_ascas12a_tracks.sh --check    run the pre-flight only, touch nothing remote
set -euo pipefail

CHECK_ONLY=0
[ "${1:-}" = "--check" ] && CHECK_ONLY=1

SRC=/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a
DEST=emalekos@hgwdev.gi.ucsc.edu:/hive/users/emalekos/ascas12a_tracks
DESTPATH=/hive/users/emalekos/ascas12a_tracks

# ##########################################################################
# STOP: the DIRS map below is STALE as of 2026-10-22 and must not be pushed as it
# stands. It points at the *_rounded dirs, which are the PRE-RESCAN tracks. Five of
# them (danRer11, rn7, mm39, hg38, hs1) were scanned with a GPU binary dated
# 2026-06-05 that predates ad3ecc6, so about a third of Cas12a mismatch types sit at
# the default penalty of 1.0 and their specificity is overstated. That is the exact
# defect the October 2026 rescan existed to remove.
#
# The verified tracks are the *_rescan dirs. All eight passed slurm/verify_track.sh:
# bigBed itemCount equal to the source guide BED, 21-field autoSql, three colour
# columns, zero invariant violations, one autoSql md5 across the campaign.
#
# Repoint every line below from _rounded to _rescan before running this, and re-read
# the rsync notes above afterwards: the reround hardlink argument does NOT hold for
# the rescan, because every off-target was rescored, so the details .tab.gz is a new
# file for all eight genomes and nothing will be skipped.
# ##########################################################################
#
# db -> source track dir, ABSOLUTE. These are the *_rounded dirs: genome-wide
# percentiles, the retired 0.20-to-0.50 grey repainted, and the reviewed score
# rounding with its re-rank, each step validated by the pipeline's own validator.
#   recolor 41864284 / 42097870, validate 41865529 / 42097870, reround 42423247
#
# Only ascas12a.bb and ascas12aTargets.as actually changed in the reround: _offset
# never moved, so the details .tab.gz is the same inode as before (hardlinked) and
# rsync -a skips it for the five genomes already on hgwdev. It is still listed in
# FILES so the MD5SUMS manifest covers every shipped file, and so mm39/hg38/hs1,
# which have never been pushed, get their details on this run.
# hg38 and hs1 were built in their own scan dirs, so they sit outside $SRC.
# The older *_cas9schema / *_cas9schema2 / *_regraded dirs predate the colour
# fix; do not push those.
DIRS="
sacCer3  $SRC/sacCer3/track_sacCer3_cas12a_rounded
ce11     $SRC/ce11/track_ce11_cas12a_rounded
dm6      $SRC/dm6/track_dm6_cas12a_rounded
danRer11 $SRC/danRer11/track_danRer11_cas12a_rounded
rn7      $SRC/rn7/track_rn7_cas12a_rounded
mm39     $SRC/mm39/track_mm39_cas12a_rounded
hg38     /private/groups/carpenterlab/emalekos/tmp/hg38-cas12a/track_hg38_cas12a_rounded
hs1      /private/groups/carpenterlab/emalekos/tmp/t2t-cas12a/track_t2t_cas12a_rounded
"
FILES="ascas12a.bb ascas12a_crisprDetails.tab.gz ascas12a_crisprDetails.tab.gz.gzi ascas12aTargets.as"

# Pre-flight: every deliverable must be present before anything is sent, so a
# genome still being recolored aborts the run up front instead of leaving a
# half-pushed directory on hgwdev.
missing=0
while read -r db dir <&3; do
    [ -z "${db:-}" ] && continue
    for f in $FILES; do
      [ -s "$dir/$f" ] || { echo "MISSING: $db  $dir/$f" >&2; missing=1; }
    done
  done 3<<< "$DIRS"
if [ "$missing" != 0 ]; then
    echo "Nothing pushed. Wait for the reround job (42423247) to finish." >&2
    exit 1
fi
echo "pre-flight OK: $(grep -c '[^[:space:]]' <<< "$DIRS") genomes, $(du -shc $(while read -r db dir <&3; do [ -n "${db:-}" ] && echo "$dir"; done 3<<< "$DIRS") 2>/dev/null | tail -1 | cut -f1) total"
if [ "$CHECK_ONLY" = 1 ]; then
    echo "--check given: nothing pushed."
    exit 0
fi

MANIFEST=$(mktemp -d)
trap 'rm -rf "$MANIFEST"' EXIT


# Loop input comes from fd 3, not stdin: ssh and rsync both read stdin, and with the
# loop reading stdin they would swallow the remaining lines and silently push only the
# first genome. ssh also gets -n for the same reason.
while read -r db dir <&3; do
    [ -z "${db:-}" ] && continue
    echo "=== $db ==="
    ssh -n emalekos@hgwdev.gi.ucsc.edu "mkdir -p $DESTPATH/$db"
    for f in $FILES; do
      rsync -avP --partial "$dir/$f" "$DEST/$db/$f"
    done 
    # manifest goes to scratch, not into the track dir: these dirs are verified
    # deliverables and should stay exactly as the validator saw them
    ( cd "$dir" && md5sum $FILES ) > "$MANIFEST/$db.MD5SUMS"
    rsync -avP "$MANIFEST/$db.MD5SUMS" "$DEST/$db/MD5SUMS"
    ssh -n emalekos@hgwdev.gi.ucsc.edu "cd $DESTPATH/$db && md5sum -c MD5SUMS"
  done 3<<< "$DIRS"

# description page + trackDb stanza
HUB=/private/groups/carpenterlab/emalekos/crispr-ots/crisprware/docs/trackhub
ssh -n emalekos@hgwdev.gi.ucsc.edu "mkdir -p $DESTPATH/html"
rsync -avP "$HUB/ascas12a.html" "$HUB/ascas12a.ra" "$DEST/html/"

echo "ALL PUSHED AND CHECKSUM-VERIFIED"
