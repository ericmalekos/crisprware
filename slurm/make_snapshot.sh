#!/bin/bash
# Freeze the live crisprware checkout as a snapshot a SLURM campaign can run from.
#
# Exists because the hand-rolled version of this is what caused the damage: five
# successive snapshots each copied the previous one's binaries, carrying a 2026-06-05
# pre-fix GPU build forward while the Python was correctly re-frozen each time, and
# five published genomes got wrong Cas12a specificity. So a binary is only ever
# carried over when `git diff <base>..HEAD -- crispr-ots/` is EMPTY, which is checked
# here and recorded in the snapshot's BINARIES.md.
#
# Usage:  make_snapshot.sh <base-snapshot-sha>
#   e.g.  make_snapshot.sh 2d9a403
set -euo pipefail

BASE_SHA="${1:?usage: make_snapshot.sh <base-snapshot-sha whose binary to carry over>}"
LIVE=/private/groups/carpenterlab/emalekos/crispr-ots/crisprware
SNAPDIR=/private/groups/carpenterlab/emalekos/crispr-ots/snapshots
BASE=$SNAPDIR/crisprware-$BASE_SHA
BINREL=crispr-ots/target/x86_64-unknown-linux-gnu/release/crispr-ots

[ -d "$BASE" ] || { echo "FATAL: no base snapshot $BASE" >&2; exit 1; }
[ -x "$BASE/$BINREL" ] || { echo "FATAL: base snapshot has no binary at $BINREL" >&2; exit 1; }

cd "$LIVE"
# Freezing a dirty tree would produce a snapshot no commit describes.
if [ -n "$(git status --porcelain)" ]; then
  echo "FATAL: working tree is dirty; commit before snapshotting" >&2
  git status --short >&2
  exit 1
fi
SHA=$(git rev-parse --short HEAD)
NEW=$SNAPDIR/crisprware-$SHA
[ -e "$NEW" ] && { echo "snapshot $NEW already exists; nothing to do"; exit 0; }

# The binary may be carried over ONLY if no Rust source changed.
if [ -n "$(git diff --stat "$BASE_SHA..HEAD" -- crispr-ots/)" ]; then
  echo "FATAL: Rust source changed between $BASE_SHA and $SHA." >&2
  echo "  The binary cannot be carried over; build one in the new snapshot and vouch for it." >&2
  git diff --stat "$BASE_SHA..HEAD" -- crispr-ots/ >&2
  exit 1
fi

mkdir -p "$NEW"
git archive HEAD | tar -x -C "$NEW"
mkdir -p "$NEW/$(dirname "$BINREL")"
install -m 555 "$BASE/$BINREL" "$NEW/$BINREL"
cmp "$BASE/$BINREL" "$NEW/$BINREL"
M=$(md5sum "$NEW/$BINREL" | cut -d' ' -f1)

{
  echo "# Binaries in this snapshot"
  echo
  echo "| path | md5 | provenance |"
  echo "|---|---|---|"
  echo "| \`$BINREL\` | \`${M:0:12}\` | carried over from \`crisprware-$BASE_SHA\` by make_snapshot.sh on $(date +%F), after confirming \`git diff $BASE_SHA..$SHA -- crispr-ots/\` is empty, so this snapshot's Rust source is identical to the one that binary was built from |"
  echo
  echo "Full md5: \`$M\`"
  echo
  # carry over the base's scanner note and its verification claim, nothing else:
  # the closing rule is printed once, below.
  sed -n '/^\*\*/,/^A snapshot must/p' "$BASE/BINARIES.md" | sed '/^A snapshot must/d'
  echo "## What changed from \`crisprware-$BASE_SHA\`"
  echo
  echo "Python only. Commits:"
  echo
  git log --oneline "$BASE_SHA..HEAD" | sed 's/^/* /'
  echo
  echo "A snapshot must BUILD or explicitly vouch for its binaries, never silently"
  echo "inherit them. The carry-over above is vouched for by the empty Rust diff."
} > "$NEW/BINARIES.md"

chmod -R a-w "$NEW"
echo "snapshot: $NEW"
echo "binary md5: $M"
