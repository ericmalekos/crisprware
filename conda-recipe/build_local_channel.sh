#!/usr/bin/env bash
# Build crispr-ots into a local conda channel for offline / pre-bioconda use.
#
# Usage:
#   conda-recipe/build_local_channel.sh                   # writes to ./conda-channel
#   CHANNEL_DIR=/scratch/$USER/conda-channel \
#       conda-recipe/build_local_channel.sh               # custom location
#
# After build:
#   - The .tar.bz2 / .conda artifact lands in $CHANNEL_DIR/linux-64/
#   - `conda index $CHANNEL_DIR` regenerates the channel metadata
#   - Add `- file://$CHANNEL_DIR` near the top of `environment.yml`'s
#     `channels:` list to pick it up.
#
# The `git_url: ../..` source spec in meta.yaml means conda-build clones
# whatever's currently committed on the working branch — uncommitted
# changes are skipped. Iterate by: edit code, `git commit -a --amend`,
# rerun this script.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHANNEL_DIR="${CHANNEL_DIR:-$REPO_ROOT/conda-channel}"
RECIPE_DIR="$REPO_ROOT/conda-recipe/crispr-ots"

mkdir -p "$CHANNEL_DIR"

echo "Building crispr-ots conda package into $CHANNEL_DIR"
echo "  recipe:        $RECIPE_DIR"
echo "  source commit: $(git -C "$REPO_ROOT" rev-parse --short HEAD) ($(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD))"
echo

conda build "$RECIPE_DIR" --output-folder "$CHANNEL_DIR"

# Refresh channel index so conda can find the new package.
conda index "$CHANNEL_DIR"

echo
echo "Done. To use:"
echo
echo "  channels:"
echo "    - file://$CHANNEL_DIR"
echo "    - conda-forge"
echo "    - bioconda"
echo "  dependencies:"
echo "    - crispr-ots"
echo
