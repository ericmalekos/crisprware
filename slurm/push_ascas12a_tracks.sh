#!/usr/bin/env bash
# Push the TTTV-indexed AsCas12a (bed9+12, 21-field) tracks to hgwdev.
# Run from a login node inside screen/tmux. ~80 GB; --partial makes it resumable.
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

# db -> source track dir, ABSOLUTE. These are the *_tttv dirs: the eight-genome build
# on the fixed Cas12a CFD binary AND a TTTV off-target index. They are the only tracks
# that should ever be pushed.
#
# Two earlier maps must not come back, both for reasons in
# docs/source/track_release_notes.md, entry 2026-10-22:
#   *_rounded  danRer11, rn7, mm39, hg38 and hs1 were scanned with a GPU binary dated
#              2026-06-05, predating ad3ecc6, which left about a third of Cas12a
#              mismatch types at the default penalty of 1.0, so it overstated
#              off-target activity and understated specificity.
#   *_rescan    fixed that binary but was still scanned against a TTTN index. AsCas12a
#              cannot cleave a TTTT-prefixed protospacer, and only `specificity_tttv`
#              excluded TTTT: the mismatch buckets behind the close-off-target grey,
#              the off-target list, the off-target count, max CFD, the saturation cap
#              and the 0-mm screen behind `unique` all counted TTTT sites. On hg38
#              chr1 that greyed 12.49% of guides whose own published specificity said
#              they were fine. Never deployed.
#
# Every one of the eight must pass slurm/verify_track.sh before this runs: bigBed
# itemCount equal to the source guide BED exactly, 21-field autoSql, three colour
# columns, zero off-target invariant violations, and ONE autoSql md5 across all eight.
# The md5 changes with the TTTV rebuild (the autoSql comments are unchanged, but the
# _specColor cut moved from 0.50 to 0.20), so compare the eight against each other
# rather than against a value written here.
#
# Transfer size is ~80 GB and NOTHING will be skipped, unlike the reround push. That
# one changed only ascas12a.bb and ascas12aTargets.as, left _offset in place, and so
# hardlinked the details .tab.gz and let rsync skip it. The TTTV build rescored every
# off-target and dropped every TTTT one, so all four files are new for all eight.
#
# hg38 and hs1 were built in their own scan dirs, so they sit outside $SRC, and hs1's
# track dir is named track_hs1_* while its guide files use the t2t_chm13v2 stem.
# The older *_cas9schema / *_cas9schema2 / *_regraded / *_rounded / *_rescan dirs all
# predate the TTTV build; do not push those.
DIRS="
sacCer3  $SRC/sacCer3/track_sacCer3_cas12a_tttv
ce11     $SRC/ce11/track_ce11_cas12a_tttv
dm6      $SRC/dm6/track_dm6_cas12a_tttv
danRer11 $SRC/danRer11/track_danRer11_cas12a_tttv
rn7      $SRC/rn7/track_rn7_cas12a_tttv
mm39     $SRC/mm39/track_mm39_cas12a_tttv
hg38     /private/groups/carpenterlab/emalekos/tmp/hg38-cas12a/track_hg38_cas12a_tttv
hs1      /private/groups/carpenterlab/emalekos/tmp/t2t-cas12a/track_hs1_cas12a_tttv
"
FILES="ascas12a.bb ascas12a_crisprDetails.tab.gz ascas12a_crisprDetails.tab.gz.gzi ascas12aTargets.as"

# Pre-flight: every deliverable must be present before anything is sent, so a genome
# whose merge has not landed aborts the run up front instead of leaving a half-pushed
# directory on hgwdev.
#
# Presence is NOT the same as correctness, and this loop only checks presence. The
# correctness check is slurm/verify_track.sh, which compares the bigBed's itemCount
# against the source guide BED; run it per genome before pushing, because a track that
# dropped rows inside a chunk has all four files at plausible sizes.
missing=0
while read -r db dir <&3; do
    [ -z "${db:-}" ] && continue
    case "$dir" in
      *_tttv) ;;
      *) echo "REFUSING: $db -> $dir is not a *_tttv dir" >&2; missing=1 ;;
    esac
    for f in $FILES; do
      [ -s "$dir/$f" ] || { echo "MISSING: $db  $dir/$f" >&2; missing=1; }
    done
  done 3<<< "$DIRS"
if [ "$missing" != 0 ]; then
    echo "Nothing pushed." >&2
    exit 1
fi
# Size the FILES being shipped, not the track directories. The *_rescan dirs each
# carry a chunks/ tree of tens of GB that is deliberately kept for re-merging and is
# never pushed, so `du -sh` on the directory both overstates the transfer several-fold
# and takes minutes to walk. The *_rounded dirs had no chunks/, which is why the
# original line got away with it.
echo "pre-flight OK: $(grep -c '[^[:space:]]' <<< "$DIRS") genomes, $(
  while read -r db dir <&3; do
    [ -z "${db:-}" ] && continue
    for f in $FILES; do echo "$dir/$f"; done
  done 3<<< "$DIRS" | tr '\n' '\0' | du -cb --files0-from=- 2>/dev/null \
    | tail -1 | awk '{printf "%.1f GB", $1/1073741824}'
) to transfer"
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
