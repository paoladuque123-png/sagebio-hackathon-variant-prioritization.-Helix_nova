#!/bin/bash
# =============================================================
# Step 0 — Quality control of raw reads with FastQC (+ MultiQC summary)
# Input : <FASTQ_DIR>/<ID>*_1.fastq.gz, <ID>*_2.fastq.gz
# Output: <FASTQC_RAW_DIR>/<ID>*_fastqc.html / .zip
#         <FASTQC_RAW_DIR>/multiqc_report.html (if multiqc is available)
# Usage : bash scripts/00_fastqc.sh [config/config.sh]
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

if [[ -n "$ENV_QC" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    conda activate "$ENV_QC"
fi

[[ -f "$SAMPLES" ]] || { echo "ERROR: sample list not found: $SAMPLES"; exit 1; }
mkdir -p "$FASTQC_RAW_DIR"

while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue

    r1=$(find "$FASTQ_DIR" -type f -name "${sample}*_1.fastq.gz" | head -n 1)
    r2=$(find "$FASTQ_DIR" -type f -name "${sample}*_2.fastq.gz" | head -n 1)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "ERROR: paired FASTQ not found for $sample"
        exit 1
    fi

    echo "FastQC: $sample"
    fastqc --threads "$THREADS_FASTQC" --outdir "$FASTQC_RAW_DIR" "$r1" "$r2"

done < "$SAMPLES"

if command -v multiqc >/dev/null 2>&1; then
    multiqc --force --outdir "$FASTQC_RAW_DIR" "$FASTQC_RAW_DIR"
fi

echo "FastQC finished."
