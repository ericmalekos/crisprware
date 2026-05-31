#!/usr/bin/env bash
# Run FlashFry's discover + score across `N_CHUNKS` parallel processes.
# Each chunk reads the shared mouse-genome DB (mmap'd page cache) and
# emits a per-chunk `.scored` file; we concatenate them at the end.
#
# Usage:
#   run_flashfry_chunked.sh <database-prefix> <query.fasta> <output.scored>

set -euo pipefail

DB=${1:?database prefix}
QUERY=${2:?query fasta}
OUTPUT=${3:?output scored path}
N_CHUNKS=${N_CHUNKS:-8}
JAR=${JAR:-/tmp/FlashFry-assembly-1.15.jar}
WORKDIR=$(mktemp -d -p "$(dirname "$OUTPUT")" ffchunk.XXXXXX)
trap 'rm -rf "$WORKDIR"' EXIT

# Split: each FASTA record is 2 lines (header + 23-bp). Round up.
N_RECORDS=$(grep -c '^>' "$QUERY")
LINES_PER_CHUNK=$(( (N_RECORDS + N_CHUNKS - 1) / N_CHUNKS ))
LINES_PER_CHUNK=$(( LINES_PER_CHUNK * 2 ))
echo "Splitting $N_RECORDS records into $N_CHUNKS chunks of <= $LINES_PER_CHUNK lines" >&2
split -d -l "$LINES_PER_CHUNK" "$QUERY" "$WORKDIR/chunk."

# Launch one discover+score per chunk.
PIDS=()
for chunk in "$WORKDIR"/chunk.*; do
    name=$(basename "$chunk")
    (
        java -Xmx4g -jar "$JAR" discover \
            --database "$DB" \
            --fasta "$chunk" \
            --output "$WORKDIR/$name.output" 2> "$WORKDIR/$name.discover.log"
        java -Xmx4g -jar "$JAR" score \
            --input "$WORKDIR/$name.output" \
            --output "$WORKDIR/$name.scored" \
            --scoringMetrics doench2016cfd \
            --database "$DB" 2> "$WORKDIR/$name.score.log"
    ) &
    PIDS+=($!)
done
echo "Launched ${#PIDS[@]} workers" >&2
for pid in "${PIDS[@]}"; do
    wait "$pid"
done
echo "All workers finished" >&2

# Concatenate, keeping the header from the first chunk only.
FIRST=1
> "$OUTPUT"
for f in "$WORKDIR"/chunk.*.scored; do
    if [ ! -f "$f" ]; then
        echo "Missing scored output: $f" >&2
        continue
    fi
    if [ "$FIRST" = "1" ]; then
        cat "$f" >> "$OUTPUT"
        FIRST=0
    else
        tail -n +2 "$f" >> "$OUTPUT"
    fi
done
echo "Wrote concatenated scored output to $OUTPUT" >&2
