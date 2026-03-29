#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $(basename "$0") -l <tsv_list> -a <assembly>" >&2
  echo "  -l  path to list of TSV files" >&2
  echo "  -a  assembly name (e.g. hg38, mm10)" >&2
  exit 1
}

# Parse arguments
TSV_LIST=""
ASSEMBLY=""
while getopts "l:a:" opt; do
  case $opt in
    l) TSV_LIST="$OPTARG" ;;
    a) ASSEMBLY="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$TSV_LIST" || -z "$ASSEMBLY" ]] && usage
[[ ! -f "$TSV_LIST" ]] && { echo "Error: list file '$TSV_LIST' not found" >&2; exit 1; }

# Directories
SCRIPT_DIR="/hive/users/emalekos/${ASSEMBLY}/enpam_gb_scripts"
LOG_DIR="/hive/users/emalekos/${ASSEMBLY}/enpam_gb_logs"
SCORED_DIR="/hive/users/emalekos/${ASSEMBLY}/enpam_scored_gRNAs"
SPLIT_DIR="/hive/users/emalekos/${ASSEMBLY}/gRNA_split"

mkdir -p "$SCRIPT_DIR" "$LOG_DIR" "$SCORED_DIR"

while read -r p; do
  tsv="$(basename "$p")"
  out="${SCRIPT_DIR}/enpam_${tsv%.tsv}.sh"
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -euo pipefail' \
    '' \
    'export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/$USER/.cache}"' \
    'mkdir -p "$XDG_CACHE_HOME"' \
    "mkdir -p ${LOG_DIR}" \
    "mkdir -p ${SCORED_DIR}" \
    'export MAMBA_NO_LOCK=1' \
    'export OMP_NUM_THREADS=3' \
    'export MKL_NUM_THREADS=3' \
    'export OPENBLAS_NUM_THREADS=3' \
    'export NUMEXPR_NUM_THREADS=3' \
    '' \
    'cleanup() {' \
    '  echo "Timeout or interrupt — killing process group" >&2' \
    '  kill -- -$$' \
    '  wait' \
    '}' \
    'trap cleanup TERM INT' \
    '' \
    'timeout --kill-after=30s 1h \' \
    '  taskset -c 0-2 \' \
    '  /cluster/home/emalekos/bin/micromamba run \' \
    '    -p /hive/users/emalekos/enpamgb_env \' \
    "    /hive/users/emalekos/scripts/EnPAMGB.py -i ${SPLIT_DIR}/${tsv} -o ${SCORED_DIR}/enpam_${tsv} > ${LOG_DIR}/enpam_${tsv%.tsv}.log 2>&1 &" \
    '' \
    'wait $!' \
    'EXIT_CODE=$?' \
    '' \
    'if [ $EXIT_CODE -eq 124 ]; then' \
    '  echo "Script timed out" >&2' \
    'fi' \
    '' \
    'exit $EXIT_CODE' \
    > "$out"
  chmod +x "$out"
done < "$TSV_LIST"

echo "Done — scripts written to $SCRIPT_DIR"
