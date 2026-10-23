#!/bin/bash
#SBATCH --output=/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a/logs/%x_%j.log
# Build the AsCas12a (TTTV) off-target index for one genome.
#
# This is step 1 of a new genome, and it is a separate job from the scan because it is
# CPU and memory bound with no GPU use (peak RSS measured 2.3 GB on sacCer3, 11.9 GB on
# danRer11, 17.3 GB on mm39) and finishes in seconds to minutes.
#
# The PAM is TTTV, not TTTN. AsCas12a does not cleave a TTTT-prefixed protospacer, so a
# TTTT site in the index is a false off-target, and only ONE of the scanner's outputs
# (`specificity_tttv`) ever excluded them. The mismatch-bucket counts behind the
# close-off-target grey, the off-target list in the details popup, the off-target count,
# max CFD, the saturation cap and the 0-mm screen behind `unique` all counted them.
# Excluding TTTT at the index is the only fix that cannot be forgotten at a consumer.
# On sacCer3 TTTT was 37.5% of all TTTN sites (354,690 of 945,116), concentrated in
# soft-masked repeats, so this is not a 1-in-4 correction.
#
# enCas12a is different and DOES cleave TTTT: it keeps a TTTN index. This script is for
# AsCas12a only.
#
# Submit:  sbatch -J idx_<db> -p gpu -c 4 --mem=48G -t 01:00:00 --export=ALL,ORG=<org> build_index.sh
# (-p gpu with NO --gres: short, medium and long share one saturated CPU pool on prism,
#  while the GPU nodes carry a separate pool with idle CPUs. This job needs no GPU.)
set -uo pipefail
ORG="${ORG:?set ORG=yeast|worm|fly|zebrafish|mouse|rat|human|t2t}"

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJ=${PROJ:-/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a}
REPO=${REPO:-/private/groups/carpenterlab/emalekos/crispr-ots/snapshots/crisprware-d35c812}
ENV=/private/groups/carpenterlab/emalekos/conda_envs/crisprware
MM=/private/home/emalekos/.local/bin/micromamba
BIN=$REPO/crispr-ots/target/x86_64-unknown-linux-gnu/release/crispr-ots

[ -x "$BIN" ] || { echo "FATAL: no CPU binary at $BIN (see $REPO/BINARIES.md)" >&2; exit 1; }
BIN_MD5=$(md5sum "$BIN" | cut -d' ' -f1)
# Same quarantine list as pipeline_shard.sh. An index carries no CFD scores, so a
# pre-fix binary would not corrupt it, but refusing here keeps one answer to "which
# binary built this tree" instead of two.
case "$BIN_MD5" in
  7e7ded6c0dc046b1a69804d583a71903|22e7030f06b1e994a6c81f0ab3d21807|\
  d9527c5e830b8e1d0b70897026afd4e7|7a938aa44311b79225a3f2f300acc9de|\
  4746ff52ffd4b85f44d54141b92c55f1)
    echo "FATAL: $BIN is a QUARANTINED pre-fix binary (md5 $BIN_MD5)." >&2; exit 1 ;;
esac
echo "binary: $BIN  md5=$BIN_MD5"

export PATH=$(dirname "$BIN"):/usr/bin:$PATH
export PYTHONPATH=$REPO
export TMPDIR=/data/tmp/emalekos; mkdir -p "$TMPDIR"

. "$HERE/genomes.sh"
resolve_genome "$ORG" || exit 1
mkdir -p "$PROJ/logs"

[ -s "$FNA" ] || { echo "FATAL: no canonical FASTA at $FNA" >&2; exit 1; }

# One writer per path. An index another job may be reading is never overwritten; a
# rebuild is an explicit FORCE=1 plus a manual removal, so a half-written .crot cannot
# replace a good one.
if [ -s "$IDX_TTTV.crot" ]; then
  if [ "${FORCE:-0}" = 1 ]; then
    echo "FATAL: FORCE=1 will not overwrite in place. Remove $(dirname "$IDX_TTTV") yourself first." >&2
  else
    echo "already built: $IDX_TTTV.crot ($(stat -c %s "$IDX_TTTV.crot") bytes). Nothing to do." >&2
  fi
  exit 1
fi

OUT=$(dirname "$(dirname "$IDX_TTTV")")   # the tttv/ subtree, one level above <stem>_crisprots/
mkdir -p "$OUT"
echo "START idx $ORG ($DB) $(date) host=$(hostname) pam=$ASCAS12A_PAM"
echo "fasta: $FNA  out: $OUT"
echo "ceph_used: $(getfattr --only-values -n ceph.dir.rbytes /private/groups/carpenterlab 2>/dev/null)"

# cd to the snapshot: `python -m` puts the CURRENT DIRECTORY at sys.path[0], ahead of
# PYTHONPATH, so running from anywhere else can silently import a live checkout.
cd "$REPO" || exit 1
/usr/bin/time -v -o "$PROJ/logs/time_${DB}_index_tttv.txt" \
  "$MM" run -p "$ENV" python -m crisprware.cli index_genome \
  -f "$FNA" --pam "$ASCAS12A_PAM" -l "$ASCAS12A_PROTO_LEN" --pam_5_prime \
  --bin_width "$ASCAS12A_BIN_WIDTH" -o "$OUT"
rc=$?
grep -E 'Maximum resident|Elapsed \(wall' "$PROJ/logs/time_${DB}_index_tttv.txt" || true
[ "$rc" = 0 ] || { echo "FATAL: index_genome exited $rc" >&2; exit "$rc"; }

# Verify the thing that matters: the enumerated PAM set, read back out of the header.
# A build that succeeded with the wrong PAM is the failure this whole change exists to
# prevent, so it is checked here and not left to the scan.
"$MM" run -p "$ENV" python -c "
import sys
sys.path.insert(0, '$REPO')
from crisprware.utils.crot_index import assert_enzyme_pams
h = assert_enzyme_pams('$IDX_TTTV', 'ascas12a')
print('VERIFIED', h.enzyme_name, '| PAMs', ' '.join(h.pams), '| sites', format(h.total_entries, ','),
      '| proto', h.protospacer_len, '| pam_side', h.pam_side_label, '| bin_width', h.bin_width)
" || { echo "FATAL: index header check failed" >&2; exit 1; }
echo "DONE idx $ORG ($DB) $(date)"
