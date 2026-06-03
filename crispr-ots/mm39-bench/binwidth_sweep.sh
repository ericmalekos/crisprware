#!/usr/bin/env bash
# Mouse-genome bin-width sweep.
#
# For each width: build the index, warm the page cache with a discard
# run, then measure 1-thread and 16-thread enumerate wall + RSS.
# Width 11 reuses the existing mm39.crispr-ots.crot if present so we
# don't rebuild it. Other widths' .crot files are deleted after their
# measurements.
#
# Run from a working directory that already contains the kmers CSV
# (see make_guides.py). Override REF / KMERS / WIDTHS via environment:
#
#   REF=/path/to/genome.fa.gz WIDTHS="9 10 11 12" ./binwidth_sweep.sh

set -euo pipefail

REF=${REF:?set REF to the reference genome FASTA (.gz ok)}
KMERS=${KMERS:-random_1000.kmers.csv}
RESULTS=${RESULTS:-binwidth_sweep.tsv}
read -r -a WIDTHS <<<"${WIDTHS:-9 10 11 12}"

echo -e "width\tbuild_wall_s\tbuild_rss_kb\tenum_t1_wall_s\tenum_t1_rss_kb\tenum_t16_wall_s\tenum_t16_rss_kb\tcrot_bytes" > "$RESULTS"

for W in "${WIDTHS[@]}"; do
    PREFIX="mm39.crispr-ots.w${W}"
    CROT="${PREFIX}.crot"
    BUILD_LOG="${PREFIX}.build.log"
    T1_LOG="${PREFIX}.t1.log"
    T16_LOG="${PREFIX}.t16.log"

    echo "=== width $W ===" >&2

    if [ "$W" = "11" ] && [ -f "mm39.crispr-ots.crot" ]; then
        echo "  reusing existing mm39.crispr-ots.crot (build skipped)" >&2
        PREFIX="mm39.crispr-ots"
        CROT="${PREFIX}.crot"
        BUILD_WALL="reused"
        BUILD_RSS="reused"
    else
        echo "  building..." >&2
        /usr/bin/time -v crispr-ots index --pam NGG -l 20 --bin-width "$W" \
            --index "$PREFIX" "$REF" 2> "$BUILD_LOG"
        BUILD_WALL=$(grep "Elapsed" "$BUILD_LOG" | awk '{print $NF}')
        BUILD_RSS=$(grep "Maximum resident" "$BUILD_LOG" | awk '{print $NF}')
    fi

    # Warmup: page-cache the .crot. Output discarded.
    echo "  warmup..." >&2
    crispr-ots enumerate --threads 1 --mismatches 4 --format csv \
        --spec-convention guidescan --kmers-file "$KMERS" \
        --output /tmp/discard.csv --alt-pam None --threshold -1 \
        "$PREFIX" >/dev/null 2>&1

    echo "  t1 measurement..." >&2
    /usr/bin/time -v crispr-ots enumerate --threads 1 --mismatches 4 --format csv \
        --spec-convention guidescan --kmers-file "$KMERS" \
        --output /tmp/discard.csv --alt-pam None --threshold -1 \
        "$PREFIX" 2> "$T1_LOG"
    T1_WALL=$(grep "Elapsed" "$T1_LOG" | awk '{print $NF}')
    T1_RSS=$(grep "Maximum resident" "$T1_LOG" | awk '{print $NF}')

    echo "  t16 measurement..." >&2
    /usr/bin/time -v crispr-ots enumerate --threads 16 --mismatches 4 --format csv \
        --spec-convention guidescan --kmers-file "$KMERS" \
        --output /tmp/discard.csv --alt-pam None --threshold -1 \
        "$PREFIX" 2> "$T16_LOG"
    T16_WALL=$(grep "Elapsed" "$T16_LOG" | awk '{print $NF}')
    T16_RSS=$(grep "Maximum resident" "$T16_LOG" | awk '{print $NF}')

    CROT_BYTES=$(stat -c '%s' "$CROT")

    echo -e "${W}\t${BUILD_WALL}\t${BUILD_RSS}\t${T1_WALL}\t${T1_RSS}\t${T16_WALL}\t${T16_RSS}\t${CROT_BYTES}" >> "$RESULTS"

    # Delete .crot for widths != 11 (the headline mouse index we want to keep).
    if [ "$W" != "11" ]; then
        echo "  cleaning up .crot for width $W..." >&2
        rm -f "$CROT" "${PREFIX}.forward" "${PREFIX}.reverse" "${PREFIX}.gs"
    fi

    rm -f /tmp/discard.csv
done

echo "" >&2
echo "Sweep complete. Results:" >&2
column -t "$RESULTS" >&2
