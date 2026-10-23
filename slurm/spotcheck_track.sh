#!/bin/bash
# Spot-check a finished AsCas12a track against the GENOME, not against itself.
#
# verify_track.sh proves the track is internally consistent: right item count, right
# schema, the off-target invariant holds. None of that can notice that every number was
# computed from the wrong off-target index, which is exactly what happened for eight
# genomes. So this script samples guides out of the finished bigBed, rescans just those
# guides UNFLOORED (--cfd-threshold 0, so the specificity sum is complete rather than
# trimmed to the popup's 0.023 floor), and hands the result to spotcheck_verify.py,
# which pulls every off-target's sequence out of the FASTA and recomputes the mismatch
# count, the CFD and the aggregate specificity from the published matrix.
#
# Usage:  spotcheck_track.sh <org> [suffix]      e.g.  spotcheck_track.sh human tttv
# Exit 0 only if every check passes.
set -uo pipefail

ORG="${1:?usage: spotcheck_track.sh <yeast|worm|fly|zebrafish|mouse|rat|human|t2t> [suffix]}"
SUFFIX="${2:-tttv}"
PROJ=${PROJ:-/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a}
REPO=${REPO:-/private/groups/carpenterlab/emalekos/crispr-ots/snapshots/crisprware-947524d}
UCSC=${UCSC:-/private/groups/carpenterlab/emalekos/conda_envs/ucsc/bin}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PATH=$UCSC:$PATH

# The build's own parameters. Change these together with pipeline_shard.sh or the popup
# comparisons below will fail for the wrong reason.
FLOOR=${FLOOR:-0.023}
LIST_CAP=${LIST_CAP:-1000}
MISMATCHES=${MISMATCHES:-4}
NWIN=${NWIN:-3}          # sampling windows
WINSZ=${WINSZ:-200000}   # bp per window

for t in bigBedToBed samtools bgzip; do
  command -v "$t" >/dev/null || { echo "FATAL: $t not on PATH (UCSC=$UCSC)" >&2; exit 1; }
done

GENOMES=${GENOMES:-$REPO/slurm/genomes.sh}
[ -r "$GENOMES" ] || { echo "FATAL: no genomes.sh at $GENOMES" >&2; exit 1; }
. "$GENOMES"
resolve_genome "$ORG" || exit 1
export ASCAS12A_PAM   # awk reads it through ENVIRON; genomes.sh only assigns it

BIN=$REPO/crispr-ots/target/x86_64-unknown-linux-gnu/release/crispr-ots
MATRIX=$REPO/crispr-ots/crates/crispr-score/data/cas12a_2xnls.csv
TRACK=$W/track_${DB}_cas12a_${SUFFIX}
BB=$TRACK/ascas12a.bb
DET=$TRACK/ascas12a_crisprDetails.tab.gz
D=$PROJ/${DB}_spotcheck_${SUFFIX}

for f in "$BB" "$DET" "$FNA" "$MATRIX" "$BIN" "$IDX_TTTV.crot" "$HERE/spotcheck_verify.py"; do
  [ -s "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
mkdir -p "$D" || exit 1

echo "=== $DB / $SUFFIX spot check"
echo "track : $BB"
echo "index : $IDX_TTTV"
echo "binary: $BIN  md5=$(md5sum "$BIN" | cut -d' ' -f1)"

# --- sample guides ----------------------------------------------------------------
# Windows on the largest contigs, each at a fixed fraction of the contig so the sample
# is reproducible rather than random. A big window is cheap (bigBed is indexed) and
# keeps small genomes from coming up empty.
POOL=$D/pool.bed
: > "$POOL"
sort -k2,2nr "$SIZES" | head -"$NWIN" | while read -r chrom len; do
  start=$(( len * 2 / 5 ))
  end=$(( start + WINSZ ))
  [ "$end" -gt "$len" ] && end=$len
  bigBedToBed "$BB" -chrom="$chrom" -start="$start" -end="$end" stdout >> "$POOL"
done
[ -s "$POOL" ] || { echo "FATAL: no guides sampled from $BB" >&2; exit 1; }
echo "sampled $(wc -l < "$POOL") guides from $NWIN window(s) of ${WINSZ}bp"

# Stratify: perfect, upper-mid, mid, low, and pre-screened (unique=False). A sample of
# only good guides would never exercise the saturation, floor or drop paths.
PICKED=$D/picked.bed
{
  awk -F'\t' '$14=="True" && $15+0==1'                      "$POOL" | head -2
  awk -F'\t' '$14=="True" && $15+0>=0.4 && $15+0<1'         "$POOL" | sort -t$'\t' -k15,15g | head -2
  awk -F'\t' '$14=="True" && $15+0>=0.1 && $15+0<0.4'       "$POOL" | sort -t$'\t' -k15,15g | head -3
  awk -F'\t' '$14=="True" && $15!="" && $15+0<0.1'          "$POOL" | sort -t$'\t' -k15,15g | head -2
  awk -F'\t' '$14=="False"'                                 "$POOL" | head -3
} | sort -u -t$'\t' -k1,1 -k2,2n > "$PICKED"
N=$(wc -l < "$PICKED")
[ "$N" -gt 0 ] || { echo "FATAL: stratified sample came out empty" >&2; exit 1; }
echo "picked $N guides across the specificity range"

{ echo "id,sequence,pam,chromosome,position,sense"
  awk -F'\t' '{printf "%s:%s:%s,%s,%s,%s,%s,%s\n",$1,$2,$6,$12,ENVIRON["ASCAS12A_PAM"],$1,$2,$6}' "$PICKED"
} > "$D/kmers.tsv"

# --- rescan: unfloored, plus a 0-mm pass for the pre-screened guides ---------------
# --cfd-threshold 0 keeps EVERY off-target, so the specificity sum re-derived from the
# genome is the same sum the scanner used. A 0-mm hit goes to on_count and returns
# before the off-target writer, so a dropped guide's multiplicity is only recoverable
# as 1/(on_count + 0) from a --mismatches 0 run; hence the second scan.
run() {
  srun -p short -c 8 --mem=48G -t 00:50:00 -J "spot_$DB" \
    "$BIN" enumerate --scanner cpu --threads 8 --score cfd-cas12a:2xnls \
    --output-mode both --ot-format tsv --cfd-threshold 0 "$@"
}
echo "--- unfloored rescan (mm<=$MISMATCHES)"
run --mismatches "$MISMATCHES" --threshold 0 --keep-dropped \
    --kmers-file "$D/kmers.tsv" --output "$D/out.csv" "$IDX_TTTV" || exit 1
echo "--- multiplicity pass (mm=0)"
run --mismatches 0 --kmers-file "$D/kmers.tsv" --output "$D/out0.csv" "$IDX_TTTV" || exit 1

# --- re-derive and compare --------------------------------------------------------
python3 -I "$HERE/spotcheck_verify.py" "$D" "$FNA" "$DET" "$MATRIX" "$FLOOR" "$LIST_CAP"
rc=$?
echo
[ $rc -eq 0 ] && echo "=== $DB SPOT CHECK PASSED ===" || echo "=== $DB SPOT CHECK FAILED ==="
exit $rc
