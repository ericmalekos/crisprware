#!/bin/bash
# Full RESCAN after the pre-fix Cas12a CFD binary was found in the published tracks.
#
# danRer11, rn7, mm39, hg38 and hs1 were scanned with a GPU binary dated 2026-06-05,
# which predates ad3ecc6 ("complement the d-base"). About a third of Cas12a mismatch
# types were left at the default penalty 1.0, overestimating off-target activity and
# understating specificity. That is not patchable posthoc: specificity is
# 1/(1 + sum of off-target CFDs), so every off-target has to be rescored.
# sacCer3, dm6 and ce11 used the post-fix CPU binary and were already correct; they
# are rescanned anyway so all eight come from one verified binary.
#
# Reuses the existing canonical FASTA, index and guide BED; runs ONLY step 3
# (score_guides --ucscgb). Writes to a NEW track dir so the published one is never
# clobbered.
#
# Submit:  sbatch -J c12a_<org> -p gpu --gres=gpu:A5500:<n> -c <cpus> --mem=<m> -t <t> \
#                 --export=ALL,ORG=<org> pipeline_trackonly.sh
# CPU instead:  --export=ALL,ORG=<org>,SCANNER=cpu   (fine for yeast/worm/fly)
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ORG="${ORG:?set ORG=yeast|worm|fly|zebrafish|mouse|rat}"
SCANNER="${SCANNER:-gpu}"
SUFFIX="${SUFFIX:-tttv}"   # track_<db>_cas12a_$SUFFIX; `rescan` was the TTTN-index build
# PHASE splits a chunked build across SLURM jobs. Chunk work was always independent
# (one process, own output dir, no shared state); the serial behaviour came from
# n_workers = max(1, len(gpus)), which is 1 on the CPU scanner. Fanning the same
# chunks across the CPU partition turns days of wall clock into hours.
#   split  -> write the per-chunk guide BEDs, score nothing
#   shard  -> score ONLY chunk $SHARD (set from SLURM_ARRAY_TASK_ID)
#   merge  -> merge every chunk, refusing an incomplete set
#   all    -> the original single-process behaviour
PHASE="${PHASE:-all}"
SHARD="${SHARD:-${SLURM_ARRAY_TASK_ID:-}}"

# PROJ is the campaign DATA dir (per-genome trees and logs); REPO is the frozen code
# snapshot. Both are overridable so a later campaign can point at its own data dir and
# its own snapshot without editing a script that a running array may be reading.
# REPO must stay a FROZEN snapshot rather than the live checkout: a live edit reaching a
# running job is how five genomes were published from a binary nobody had re-verified.
# See the snapshot's BINARIES.md for the binary it vouches for, and the md5 guard below.
PROJ=${PROJ:-/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a}
REPO=${REPO:-/private/groups/carpenterlab/emalekos/crispr-ots/snapshots/crisprware-947524d}
ENV=/private/groups/carpenterlab/emalekos/conda_envs/crisprware
MM=/private/home/emalekos/.local/bin/micromamba
if [ "$SCANNER" = "gpu" ]; then
  BIN=$REPO/crispr-ots/target-gpu/release/crispr-ots
else
  BIN=$REPO/crispr-ots/target/x86_64-unknown-linux-gnu/release/crispr-ots
fi
export PYTHONPATH=$REPO

# Refuse to run on a binary that cannot be trusted. This is the check whose absence
# let a 2026-06-05 pre-fix binary score five published genomes: nothing verified which
# binary the snapshot carried, and each snapshot inherited the last one's.
if [ ! -x "$BIN" ]; then
  echo "FATAL: no binary at $BIN" >&2
  echo "  SCANNER=$SCANNER needs a $([ "$SCANNER" = gpu ] && echo '--features gpu' || echo 'release') build in $REPO" >&2
  echo "  see $REPO/BINARIES.md" >&2
  exit 1
fi
BIN_MD5=$(md5sum "$BIN" | cut -d' ' -f1)
case "$BIN_MD5" in
  7e7ded6c0dc046b1a69804d583a71903|22e7030f06b1e994a6c81f0ab3d21807|\
  d9527c5e830b8e1d0b70897026afd4e7|7a938aa44311b79225a3f2f300acc9de|\
  4746ff52ffd4b85f44d54141b92c55f1)
    echo "FATAL: $BIN is a QUARANTINED pre-fix binary (md5 $BIN_MD5)." >&2
    echo "  Its Cas12a off-target CFDs are wrong; that is what this rescan exists to undo." >&2
    exit 1 ;;
esac
echo "binary: $BIN  md5=$BIN_MD5"

export PATH=/private/groups/carpenterlab/emalekos/conda_envs/ucsc/bin:$(dirname "$BIN"):/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}
export TMPDIR=/data/tmp/emalekos; mkdir -p "$TMPDIR"
THREADS=${THREADS:-16}

# Per-genome paths come from genomes.sh so this script and the chain submitter cannot
# disagree about where a genome lives or how large its chunks are.
. "$HERE/genomes.sh"
resolve_genome "$ORG" || exit 1

# CHUNK override. The small genomes default to CHUNK=0 (unchunked), so forcing a chunk
# size is how the sharded path gets validated: a chunked+sharded ce11 must reproduce
# the already-verified unchunked ce11 byte for byte, which is exactly the chunked ==
# unchunked equivalence the merge-time percentile re-rank exists to guarantee.
CHUNK="${CHUNK_OVERRIDE:-$CHUNK}"

# The TTTV off-target index, not the TTTN one the first eight tracks were built from.
# AsCas12a cannot cleave a TTTT-prefixed protospacer, so a TTTT site in the index is a
# false off-target. Only `specificity_tttv` ever excluded them; the mismatch buckets
# behind the close-off-target grey, the off-target list, the off-target count, max CFD,
# the saturation cap and the 0-mm screen behind `unique` all counted them. Indexing on
# TTTV removes them at the source instead of relying on a guard at each consumer.
IDX=$IDX_TTTV
TRACK=$W/track_${DB}_cas12a_${SUFFIX}   # NEW dir; published track_${DB}_cas12a untouched
mkdir -p "$TRACK" "$PROJ/logs"

for f in "$IDX.crot" "$GUIDES" "$SIZES"; do
  [ -s "$f" ] || { echo "MISSING PREREQ: $f"; exit 1; }
done
# Refuse a TTTT-bearing index up front. This is the check whose absence let eight
# genomes be scanned against TTTN while three comments asserted the index was TTTV.
"$MM" run -p "$ENV" python -c "
import sys
sys.path.insert(0, '$REPO')
from crisprware.utils.crot_index import assert_enzyme_pams
h = assert_enzyme_pams('$IDX', 'ascas12a')
print('index:', h.enzyme_name, '| PAMs', ' '.join(h.pams), '| sites', format(h.total_entries, ','))
" || { echo "FATAL: off-target index is not usable for AsCas12a" >&2; exit 1; }
echo "START $ORG ($DB) $(date) host=$(hostname) scanner=$SCANNER chunk=$CHUNK phase=$PHASE shard=${SHARD:-none} cuda=${CUDA_VISIBLE_DEVICES:-unset}"
echo "guides: $(($(wc -l < "$GUIDES") - 1))   ceph_used: $(getfattr --only-values -n ceph.dir.rbytes /private/groups/carpenterlab 2>/dev/null)"

SCORE_ARGS=( -b "$GUIDES" -i "$IDX"
  --cas12a_scorer deepcpf1 enpam_gb enseq_deepcpf1 seq_deepcpf1variants --cas12a_variant AsCas12a
  --ucscgb "$TRACK" --chrom_sizes "$SIZES" --crispr_ots_bin "$BIN"
  --ucscgb_scanner "$SCANNER" --ucscgb_cfd_threshold 0.023 --ucscgb_list_cap 1000 --ucscgb_blank_threshold 0
  --mismatches 4 --threads "$THREADS" -o "$W" )
[ "$CHUNK" -gt 0 ] && SCORE_ARGS+=( --ucscgb_chunk_max "$CHUNK" )
if [ "$PHASE" != "all" ]; then
  [ "$CHUNK" -gt 0 ] || { echo "FATAL: PHASE=$PHASE needs a chunked genome (CHUNK>0), got CHUNK=$CHUNK" >&2; exit 1; }
  SCORE_ARGS+=( --ucscgb_phase "$PHASE" --ucscgb_keep_chunks )
  if [ "$PHASE" = "shard" ]; then
    [ -n "$SHARD" ] || { echo "FATAL: PHASE=shard needs SHARD or SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    SCORE_ARGS+=( --ucscgb_shard "$SHARD" )
  fi
fi

# Pin TensorFlow's thread count. The cleavage scores are float32 TF output and depend
# on BOTH the host ISA group and the intra-op thread count, measured 2026-10-20 on
# 200,000 ce11 contexts:
#   * the avx2+fma nodes (every GPU node: phoenix-00 EPYC 7713, phoenix-01 to 05 EPYC
#     7662) agree BYTE FOR BYTE with each other, microarchitecture notwithstanding;
#   * the avx512 nodes (phoenix-09, 10, 22, 23, 24, EPYC 9684X) agree with each other
#     and differ from the avx2 group on 53% of values, up to 2.3e-05;
#   * on one host, 4 threads vs 8 threads differ on 3 of 200,000, up to 3.8e-06.
# Those low-bit differences flip the DISPLAYED value for about 3.3 ppm of guides.
# `nproc` does not respect the cgroup on every node (phoenix-23 reports 384 under a
# 4-CPU allocation), so TF's default thread count is not a function of the allocation
# and cannot be relied on. Pinning it to the requested CPU count makes the one
# remaining knob explicit and equal across every shard of a campaign.
export TF_NUM_INTRAOP_THREADS=${SLURM_CPUS_PER_TASK:-16}
export TF_NUM_INTEROP_THREADS=1
echo "tf threads: intra=$TF_NUM_INTRAOP_THREADS inter=$TF_NUM_INTEROP_THREADS"

TAG=$SUFFIX
[ "$PHASE" != "all" ] && TAG="${SUFFIX}_${PHASE}${SHARD:+_$SHARD}"
# `python -m` puts the CURRENT DIRECTORY at sys.path[0], AHEAD of PYTHONPATH. A job
# whose working directory happens to be the live checkout therefore imports the live
# crisprware and silently ignores the frozen snapshot. The ce11 shard test
# (42472185 to 42472187) ran live-repo Python for exactly this reason; it was harmless
# only because the working tree happened to sit on the snapshot's own commit. cd into
# the snapshot so cwd and PYTHONPATH name one tree, then prove which file was imported
# before any work starts. Same spirit as the binary md5 guard above: assert provenance,
# do not assume it.
cd "$REPO"
RESOLVED=$("$MM" run -p "$ENV" python -c 'import crisprware.score_guides as m; print(m.__file__)')
case "$RESOLVED" in
  "$REPO"/*) echo "crisprware: $RESOLVED" ;;
  *) echo "FATAL: crisprware resolved to $RESOLVED, outside the frozen snapshot $REPO" >&2
     exit 1 ;;
esac

/usr/bin/time -v -o "$PROJ/logs/time_${DB}_track_${TAG}.txt" \
  "$MM" run -p "$ENV" python -m crisprware.cli score_guides "${SCORE_ARGS[@]}"
rc=$?
grep -E 'Maximum resident|Elapsed \(wall' "$PROJ/logs/time_${DB}_track_${TAG}.txt" || true
[ $rc -eq 0 ] || { echo "SCORE_GUIDES FAILED rc=$rc"; exit $rc; }

# The checks below inspect a FINISHED track: its autoSql, its field count and its
# details file. Only the merge (or an unsharded "all" run) produces one. A split
# writes chunk BEDs and a shard writes one chunk's track, so running these after
# either would fail on a track that was never supposed to exist yet.
if [ "$PHASE" = "split" ] || [ "$PHASE" = "shard" ]; then
  echo "=== phase $PHASE complete${SHARD:+ (shard $SHARD)}; no finished track to verify yet ==="
  ls -la "$TRACK/chunks" | tail -5
  echo "DONE $ORG ($DB) phase=$PHASE $(date)"
  exit 0
fi

echo "=== track files ($DB) ==="; ls -la "$TRACK"

echo "=== schema check ($DB) ==="
AS=$(ls "$TRACK"/*Targets.as 2>/dev/null | head -1)
for f in itemRgb _crisprDeepCpf1Color _specColor deepcpf1 deepcpf1Perc ascas12aDeepcpf1 ascas12aDeepcpf1Perc; do
  grep -q "[ \t]$f;" "$AS" || { echo "SCHEMA FAIL: $f missing from $AS"; exit 1; }
done
ncol=$(grep -cE '^\s+(string|uint|char|bigint)' "$AS")
[ "$ncol" -eq 21 ] || { echo "SCHEMA FAIL: expected 21 fields, got $ncol"; exit 1; }
echo "schema OK: 21 fields, three colour columns present"

echo "=== invariant check ($DB) ==="
"$MM" run -p "$ENV" python - "$TRACK" <<'PY'
import gzip, glob, os, sys
hits = glob.glob(os.path.join(sys.argv[1], "*crisprDetails.tab.gz"))
if not hits:
    print("FAIL: no crisprDetails.tab.gz"); sys.exit(1)
ok = bad = over = 0; first = None
with gzip.open(hits[0], "rt") as f:
    next(f)
    for line in f:
        c, l = line.rstrip("\n").split("\t")
        s = sum(int(x) for x in c.split(",")); n = len(l.split("|")) if l else 0
        if l == "" and s > 0: over += 1
        elif s == n: ok += 1
        else:
            bad += 1
            if first is None: first = (c, n)
print(f"invariant holds     : {ok}")
print(f"over-cap empty list : {over}")
print(f"VIOLATIONS          : {bad}" + (f"   first={first}" if first else ""))
sys.exit(1 if bad else 0)
PY
inv=$?

# Drop the bulky per-off-target intermediates; ceph is near the 15 TB lab cap.
rm -f "$TRACK"/offtargets.csv.ot.tsv "$TRACK"/guides.kmers.csv
echo "invariant_exit=$inv  ceph_used_after: $(getfattr --only-values -n ceph.dir.rbytes /private/groups/carpenterlab 2>/dev/null)"
echo "DONE $ORG ($DB) $(date)"
exit $inv
