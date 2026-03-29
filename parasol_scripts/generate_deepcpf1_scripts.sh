#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $(basename "$0") -l <tsv_list> -a <assembly>" >&2
  echo "  -l  path to list of TSV files" >&2
  echo "  -a  assembly name (e.g. hg38, mm10, ce11)" >&2
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
SCRIPT_DIR="/hive/users/emalekos/${ASSEMBLY}/deep_cpf1_scripts"
LOG_DIR="/hive/users/emalekos/${ASSEMBLY}/deepcpf1_logs"
SCORED_DIR="/hive/users/emalekos/${ASSEMBLY}/enpam_scored_gRNAs"

mkdir -p "$SCRIPT_DIR" "$LOG_DIR" "$SCORED_DIR"

while read -r p; do
  tsv="$(basename "$p")"
  # Replace leading "enpam_" with "deepcpf1_" in the script name
  script_name="deepcpf1_${tsv#enpam_}"
  out="${SCRIPT_DIR}/${script_name%.tsv}.sh"
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'set -euo pipefail' \
    '' \
    'export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/$USER/.cache}"' \
    'mkdir -p "$XDG_CACHE_HOME"' \
    "mkdir -p ${LOG_DIR}" \
    "mkdir -p ${SCORED_DIR}" \
    'export MAMBA_NO_LOCK=1' \
    'export OMP_NUM_THREADS=1' \
    'export MKL_NUM_THREADS=1' \
    'export OPENBLAS_NUM_THREADS=1' \
    'export NUMEXPR_NUM_THREADS=1' \
    "export THEANO_FLAGS=\"base_compiledir=/tmp/\${USER}/.theano/${script_name%.tsv},openmp=True,openmp_elemwise_minsize=200000\"" \
    '' \
    'cleanup() {' \
    '  echo "Timeout or interrupt - killing process group" >&2' \
    '  kill -- -$$' \
    '  wait' \
    '}' \
    'trap cleanup TERM INT' \
    '' \
    'timeout --kill-after=30s 1h \' \
    '  taskset -c 0 \' \
    '  /cluster/home/emalekos/bin/micromamba run \' \
    '    -p /hive/users/emalekos/deepcpf1 \' \
    "    python /hive/users/emalekos/scripts/DeepCpf1_seq.py -i ${SCORED_DIR}/${tsv} -o ${SCORED_DIR}/${script_name} > ${LOG_DIR}/${script_name%.tsv}.log 2>&1 &" \
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

echo "Done - scripts written to $SCRIPT_DIR"
