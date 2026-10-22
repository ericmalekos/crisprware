#!/bin/bash
# Verify a finished --ucscgb track, sharded or not.
#
# Why this exists as its own script: the build's own post-merge checks confirm the
# schema and the off-target invariant, and the merge refuses to run on an incomplete
# chunk set, but NONE of that can detect rows dropped INSIDE a chunk. The merge counts
# chunks, not guides. Comparing the bigBed's itemCount against the source guide BED
# does, so that is the load-bearing check here and it runs first.
#
# The other three checks are re-run rather than trusted from the build log, so a track
# can be verified long after the job that made it, and so a track built by hand or by a
# superseded driver can be held to the same standard.
#
# Usage:  verify_track.sh <org> [suffix]
#   e.g.  verify_track.sh human rescan
# Exit 0 only if every check passes.
set -uo pipefail

ORG="${1:?usage: verify_track.sh <yeast|worm|fly|zebrafish|mouse|rat|human|t2t> [suffix]}"
SUFFIX="${2:-rescan}"
PROJ=${PROJ:-/private/groups/carpenterlab/emalekos/tmp/multigenome-cas12a}
ENV=${ENV:-/private/groups/carpenterlab/emalekos/conda_envs/crisprware}
MM=${MM:-/private/home/emalekos/.local/bin/micromamba}
UCSC=${UCSC:-/private/groups/carpenterlab/emalekos/conda_envs/ucsc/bin}
export PATH=$UCSC:$PATH

# Assert the tool is reachable BEFORE any work. A verdict job that dies on a missing
# bigBedToBed after its schema check has already passed reads exactly like "the track
# does not match", and that misreading nearly cost a validated design once.
command -v bigBedInfo >/dev/null || {
  echo "FATAL: bigBedInfo not on PATH (looked in $UCSC)" >&2; exit 1; }

# Mirrors pipeline_shard.sh: hg38 and hs1 predate the <db>/<db>_canon_* convention and
# keep their own directories and stems, so every script has to derive the same paths.
case "$ORG" in
  yeast) DB=sacCer3 ;; worm) DB=ce11 ;; fly) DB=dm6 ;;
  zebrafish) DB=danRer11 ;; mouse) DB=mm39 ;; rat) DB=rn7 ;;
  human) DB=hg38; W=/private/groups/carpenterlab/emalekos/tmp/hg38-cas12a
         STEM=hg38_primary_chr ;;
  t2t)   DB=hs1;  W=/private/groups/carpenterlab/emalekos/tmp/t2t-cas12a
         STEM=t2t_chm13v2 ;;
  *) echo "FATAL: unknown ORG $ORG" >&2; exit 1 ;;
esac
W=${W:-$PROJ/$DB}
STEM=${STEM:-${DB}_canon}
BED=$W/${STEM}_gRNA/${STEM}_gRNA.bed
TRACK=$W/track_${DB}_cas12a_${SUFFIX}

[ -s "$BED" ]  || { echo "FATAL: no guide BED at $BED" >&2; exit 1; }
[ -d "$TRACK" ] || { echo "FATAL: no track dir at $TRACK" >&2; exit 1; }
BB=$(ls "$TRACK"/*.bb 2>/dev/null | head -1)
[ -n "$BB" ] || { echo "FATAL: no bigBed in $TRACK" >&2; exit 1; }

fail=0
echo "=== $DB  track=$TRACK"

# 1. Row count. The check nothing else makes.
echo "--- row count"
ITEMS=$(bigBedInfo "$BB" | awk '/itemCount/{gsub(",","",$2); print $2}')
SRC=$(( $(wc -l < "$BED") - 1 ))
printf '  bigBed itemCount : %s\n  source guides    : %s\n' "$ITEMS" "$SRC"
if [ "$ITEMS" = "$SRC" ]; then
  echo "  row count OK"
else
  echo "  ROW COUNT FAIL: $(( SRC - ITEMS )) guides missing from the track"; fail=1
fi

# 2. Schema. Same assertions the build makes, re-run against the file on disk.
echo "--- schema"
AS=$(ls "$TRACK"/*Targets.as 2>/dev/null | head -1)
if [ -z "$AS" ]; then
  echo "  SCHEMA FAIL: no autoSql in $TRACK"; fail=1
else
  for f in itemRgb _crisprDeepCpf1Color _specColor deepcpf1 deepcpf1Perc \
           ascas12aDeepcpf1 ascas12aDeepcpf1Perc; do
    grep -q "[ \t]$f;" "$AS" || { echo "  SCHEMA FAIL: $f missing from $AS"; fail=1; }
  done
  ncol=$(grep -cE '^\s+(string|uint|char|bigint)' "$AS")
  if [ "$ncol" -eq 21 ]; then
    echo "  schema OK: 21 fields, three colour columns present"
  else
    echo "  SCHEMA FAIL: expected 21 fields, got $ncol"; fail=1
  fi
  # Printed so tracks can be compared across genomes. One md5 across a campaign means
  # one schema, rather than N schemas that each independently pass a field count.
  echo "  autoSql md5      : $(md5sum "$AS" | cut -d' ' -f1)"
fi

# 3. Off-target invariant: sum(_mismatchCounts) == len(off-target list), except where
#    the list was blanked by the per-bin cap, which legitimately leaves counts with an
#    empty list.
echo "--- off-target invariant"
"$MM" run -p "$ENV" python - "$TRACK" <<'PY'
import gzip, glob, os, sys
hits = glob.glob(os.path.join(sys.argv[1], "*crisprDetails.tab.gz"))
if not hits:
    print("  FAIL: no crisprDetails.tab.gz"); sys.exit(1)
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
print(f"  invariant holds     : {ok}")
print(f"  over-cap empty list : {over}")
print(f"  VIOLATIONS          : {bad}" + (f"   first={first}" if first else ""))
sys.exit(1 if bad else 0)
PY
[ $? -eq 0 ] || fail=1

# 4. Side-car reachability. A bigBed with no usable details file still passes every
#    check above, and the track would then render with empty off-target popups.
echo "--- details side-car"
DET=$(ls "$TRACK"/*crisprDetails.tab.gz 2>/dev/null | head -1)
if [ -n "$DET" ] && [ -s "${DET}.gzi" ]; then
  echo "  bgzipped details plus .gzi index present"
else
  echo "  DETAILS FAIL: missing ${DET:-crisprDetails.tab.gz} or its .gzi index"; fail=1
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "=== $DB VERIFIED ==="
else
  echo "=== $DB FAILED VERIFICATION ==="
fi
exit $fail
