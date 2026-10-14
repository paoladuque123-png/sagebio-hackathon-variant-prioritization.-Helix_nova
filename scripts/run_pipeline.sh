#!/bin/bash
# =============================================================
# Run the full pipeline in order, logging each step to LOG_DIR.
# Usage : bash scripts/run_pipeline.sh [config/config.sh] [first_step]
#   first_step: 0 fastqc | 1 fastp | 2 fastq_to_vcf | 3 vep | 4 prioritize | 5 clinvar_plp
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
FROM="${2:-0}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$LOG_DIR"
STAMP="$(date +%F_%H%M)"

STEPS=(00_fastqc.sh 01_fastp.sh 02_fastq_to_vcf.sh 03_vep_annotate.sh 04_prioritize.sh 05_clinvar_plp.sh)
for i in "${!STEPS[@]}"; do
    (( i < FROM )) && continue
    s="${STEPS[$i]}"
    log="${LOG_DIR}/${s%.sh}_${STAMP}.log"
    echo "[$(date +%T)] ${s} -> ${log}"
    bash "${SCRIPT_DIR}/${s}" "$CONFIG" > "$log" 2>&1
done

bash "${SCRIPT_DIR}/report_versions.sh" "$CONFIG" > "${LOG_DIR}/versions_${STAMP}.txt" 2>&1 || true
echo "[$(date +%T)] Pipeline finished."
