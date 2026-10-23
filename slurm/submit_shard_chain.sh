#!/bin/bash
# Submit one genome's rescan as split -> shard array -> merge.
#
# Why this shape: chunk work was always independent (each chunk is its own process with
# its own output dir and no shared state), but the in-process worker count is
# max(1, len(gpus)), which is 1 on the CPU scanner, so a 68-chunk genome ran strictly
# serially. The array runs the same chunks concurrently.
#
# Two properties worth knowing before using it:
#   * the split is the ONLY writer of chunk_*.bed, so it must finish before any shard
#     starts; N shards splitting concurrently would interleave their writes and produce
#     chunks that are individually well-formed and collectively wrong;
#   * the campaign is resumable at shard granularity. Chunks are kept, so a shard that
#     times out or fails costs that shard, not the genome: resubmit the array with just
#     that index and then rerun the merge, which refuses an incomplete set anyway.
#
# Usage:  submit_shard_chain.sh <org> [chunk_max] [shard_hours]
#   e.g.  submit_shard_chain.sh zebrafish 2000000 24
set -euo pipefail

ORG="${1:?usage: submit_shard_chain.sh <yeast|worm|fly|zebrafish|mouse|rat|human|t2t> [chunk_max] [shard_hours] [cpu|gpu]}"
CHUNK="${2:-2000000}"
SHARD_H="${3:-24}"
SCANNER="${4:-cpu}"
# SHARDS_ONLY=1 for a genome whose split has ALREADY run: skip the split job and submit
# the array with NO dependency. Both halves matter. Re-splitting would rewrite chunk
# BEDs that are already correct (60 GB across the four 1M genomes) and breaks the
# one-writer rule for no gain, since the split is a pure function of (guide BED,
# chunk_max). And `PriorityFlags` is EMPTY on this cluster, so the age factor runs from
# eligible time, not submit time: an array held on `Dependency` banks zero age. danRer11
# waited 18.3 h on Reason=Priority before placing, so a dependency-blocked array would
# have started that clock only after its split finished.
SHARDS_ONLY="${SHARDS_ONLY:-0}"
# THROTTLE caps concurrent array tasks (the %N suffix). danRer11 ran 37 shards
# unthrottled and self-limited to about 11 concurrent as GPUs freed, which was fine at
# that size. 501 shards is 13.5x that, and the gpu partition has 48 GPUs shared with
# other users, so the larger genomes get a cap. Empty means no cap.
THROTTLE="${THROTTLE:-}"
# PROJ is the campaign DATA dir (per-genome trees and logs), overridable for a later
# campaign. HERE is where these scripts live, which is a different thing: the worker is
# resolved next to this submitter rather than inside PROJ, so the pair stays together
# when the data dir moves. sbatch needs an absolute path, hence the pwd.
PROJ=${PROJ:-/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORKER=$HERE/pipeline_shard.sh
[ -x "$WORKER" ] || { echo "FATAL: no executable worker at $WORKER" >&2; exit 1; }
# SUFFIX names the track directory (track_<db>_cas12a_$SUFFIX), so it is how one
# campaign's output is kept off another's. `tttv` is the TTTV-index rebuild; `rescan`
# was the CFD-binary rescan, which was scanned against a TTTN index and must not be
# pushed. See docs/source/track_release_notes.md, 2026-10-22.
SUFFIX=${SUFFIX:-tttv}

# Per-genome paths come from genomes.sh, which pipeline_shard.sh also sources. They have
# to derive the same paths or the chain would split one tree and shard another, and hg38
# and hs1 predate the <db>/<db>_canon_* convention, so the mapping is not derivable.
. "$HERE/genomes.sh"
resolve_genome "$ORG" || exit 1
BED=$GUIDES
[ -s "$BED" ] || { echo "FATAL: no guide BED at $BED" >&2; exit 1; }
# Refuse before submitting anything: the array would otherwise fail per task, 500 times,
# on an index that was wrong before the first sbatch.
[ -s "$IDX_TTTV.crot" ] || {
  echo "FATAL: no TTTV off-target index at $IDX_TTTV.crot" >&2
  echo "  build it first:  sbatch -J idx_$DB -p gpu -c 4 --mem=48G -t 01:00:00 --export=ALL,ORG=$ORG $HERE/build_index.sh" >&2
  exit 1; }

# Must agree with _ucscgb_expected_chunks: ceil((lines - header) / chunk_max).
# One pass, not two: the hg38 and hs1 guide BEDs are 16 GB and 17 GB, so counting lines
# in Python and then again with `wc -l` for the echo would read 33 GB off the login node
# to print one number. `wc -l` alone is the fast path.
LINES=$(wc -l < "$BED")
N=$(python3 -c "
import math
print(max(1, math.ceil(($LINES - 1) / $CHUNK)))")
LAST=$((N - 1))
echo "$DB: $((LINES - 1)) guides -> $N chunks of <= $CHUNK"

# One writer per path. A chunks/ tree another job still holds open is the documented
# way to get a corrupt output of exactly the right size.
TRACK=$W/track_${DB}_cas12a_${SUFFIX}
if squeue -u "$USER" -h -o '%j' | grep -q "c12a_${DB}_"; then
  echo "FATAL: a c12a_${DB}_* job is already queued or running; it may own $TRACK/chunks" >&2
  squeue -u "$USER" -o '%.12i %.18j %.2t' | grep "c12a_${DB}_" >&2
  exit 1
fi

COMMON="--export=ALL,ORG=$ORG,SCANNER=$SCANNER,SUFFIX=$SUFFIX,CHUNK_OVERRIDE=$CHUNK"

# The shard array is the only phase whose hardware differs. `--gres=gpu:1` and never a
# model: naming one restricts the job to a subset of nodes and it then queues behind
# that subset while other GPUs sit idle. Cas12a does not care which GPU it gets, since
# the GPU only scans and the host does the scoring.
case "$SCANNER" in
  gpu) SHARD_SBATCH=(-p gpu --gres=gpu:1 -c 16 --mem=96G)
       echo "  scanner: gpu (1 GPU + 16 CPUs per shard; GPU-node CPUs are a separate pool from long)" ;;
  cpu) SHARD_SBATCH=(-p long -c 16 --mem=96G) ;;
  *)   echo "FATAL: SCANNER must be cpu or gpu, got $SCANNER" >&2; exit 1 ;;
esac

ARRAY_SPEC="0-$LAST${THROTTLE:+%$THROTTLE}"

# RESUME=1 submits only the chunks that are NOT already complete. The sharded design has
# always been resumable in principle, since chunks are kept and a failed shard costs that
# shard rather than the genome, but acting on it meant working out the index list by hand,
# which is where an off-by-one costs a whole merge.
#
# "Complete" uses the SAME criterion _ucscgb_merge enforces, matched by glob so it stays
# enzyme-agnostic: guides.bed, the bigBed, the autoSql, and the bgzipped details. The
# bigBed matters most because build_track writes it last, so its presence is what
# separates a finished chunk from one killed mid-assembly.
if [ "${RESUME:-0}" = "1" ]; then
  TODO=()
  for i in $(seq 0 "$LAST"); do
    d=$(printf '%s/chunks/track_%03d' "$TRACK" "$i")
    ok=1
    [ -s "$d/guides.bed" ] || ok=0
    compgen -G "$d/*.bb" >/dev/null || ok=0
    compgen -G "$d/*Targets.as" >/dev/null || ok=0
    compgen -G "$d/*_crisprDetails.tab.gz" >/dev/null || ok=0
    [ "$ok" -eq 1 ] || TODO+=("$i")
  done
  done_n=$(( N - ${#TODO[@]} ))
  echo "  resume: $done_n of $N chunks already complete, ${#TODO[@]} to run"
  if [ "${#TODO[@]}" -eq 0 ]; then
    echo "  all $N chunks complete; submitting the merge alone"
    ARRAY_SPEC=""
  else
    ARRAY_SPEC="$(IFS=,; echo "${TODO[*]}")${THROTTLE:+%$THROTTLE}"
  fi
fi

if [ "$SHARDS_ONLY" = "1" ]; then
  # Assert the split really did run, for THIS chunk_max. The count is recomputed from
  # the guide BED exactly as the merge recomputes it, so a chunk set written at a
  # different chunk_max is caught here rather than at merge time after N GPU hours.
  CHUNKS=$TRACK/chunks
  have=$(ls "$CHUNKS"/chunk_*.bed 2>/dev/null | wc -l)
  [ "$have" -eq "$N" ] || {
    echo "FATAL: SHARDS_ONLY=1 but $CHUNKS has $have chunk BEDs, expected $N" >&2
    echo "  either the split has not run at chunk_max $CHUNK, or it ran at a different one" >&2
    exit 1; }
  echo "  split   SKIPPED ($have chunk BEDs present, matching $N)"
  SPLIT_DEP=()
else
  SPLIT=$(sbatch --parsable -J "c12a_${DB}_split" -p medium -c 2 --mem=8G -t 00:30:00 \
    -o "$PROJ/logs/${SUFFIX}_${DB}_split_%j.log" \
    $COMMON,PHASE=split "$WORKER")
  echo "  split  $SPLIT"
  SPLIT_DEP=(--dependency="afterok:$SPLIT")
fi

if [ -n "$ARRAY_SPEC" ]; then
  ARR=$(sbatch --parsable -J "c12a_${DB}_shard" "${SHARD_SBATCH[@]}" -t "${SHARD_H}:00:00" \
    --array="$ARRAY_SPEC" "${SPLIT_DEP[@]+"${SPLIT_DEP[@]}"}" \
    -o "$PROJ/logs/${SUFFIX}_${DB}_shard_%A_%a.log" \
    $COMMON,PHASE=shard "$WORKER")
  echo "  shards $ARR [$ARRAY_SPEC]"
  MERGE_DEP=(--dependency="afterok:$ARR")
else
  # Nothing left to scan. The merge still runs, and it re-derives the expected chunk
  # count from the guide BED and refuses an incomplete set, so it remains the check on
  # whether "all complete" was true rather than trusting the resume scan above.
  ARR="(none)"
  MERGE_DEP=()
fi

# Merge partition and limit. The merge is pure I/O (reads each chunk's guides.bed and
# BGZF details, rebases offsets, re-ranks percentiles, writes the bigBed): no TensorFlow,
# no scan, so no GPU and no avx2/avx512 cleavage-score sensitivity. Observed 1394 s to
# 1942 s for 37 to 113 chunks, so 2 h not 8 h; an 8 h request has to find an 8 h window.
# MERGE_PART defaults to medium but gpu is often the right answer, because short/medium/
# long share one saturated CPU pool while the GPU nodes carry hundreds of idle CPUs.
MERGE=$(sbatch --parsable -J "c12a_${DB}_merge" -p "${MERGE_PART:-medium}" -c 16 --mem=64G -t "${MERGE_HOURS:-2}:00:00" \
  "${MERGE_DEP[@]+"${MERGE_DEP[@]}"}" \
  -o "$PROJ/logs/${SUFFIX}_${DB}_merge_%j.log" \
  $COMMON,PHASE=merge "$WORKER")
echo "  merge  $MERGE"
echo "chain: ${SPLIT:-(split skipped)} -> $ARR[$ARRAY_SPEC] -> $MERGE  ->  $TRACK"
