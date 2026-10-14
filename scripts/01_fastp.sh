#!/bin/bash
# =============================================================
# Step 1 — Read QC and trimming with fastp
# Input : <FASTQ_DIR>/<ID>*_1.fastq.gz, <ID>*_2.fastq.gz
# Output: <TRIMMED_DIR>/<ID>_1.fq.gz, <ID>_2.fq.gz
#         <QC_DIR>/<ID>_report_fastp.html, <ID>_fastp.json
#         <FASTQC_TRIM_DIR>/<ID>_{1,2}_fastqc.html (FastQC on trimmed reads)
# Usage : bash scripts/01_fastp.sh [config/config.sh]
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

if [[ -n "$ENV_FASTP" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    conda activate "$ENV_FASTP"
fi

[[ -f "$SAMPLES" ]] || { echo "ERROR: sample list not found: $SAMPLES"; exit 1; }
mkdir -p "$TRIMMED_DIR" "$QC_DIR"

while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue

    r1=$(find "$FASTQ_DIR" -type f -name "${sample}*_1.fastq.gz" | head -n 1)
    r2=$(find "$FASTQ_DIR" -type f -name "${sample}*_2.fastq.gz" | head -n 1)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "ERROR: paired FASTQ not found for $sample"
        exit 1
    fi

    echo "fastp: $sample"

    fastp --thread "$THREADS_FASTP" \
        -i "$r1" -I "$r2" \
        -o "${TRIMMED_DIR}/${sample}_1.fq.gz" \
        -O "${TRIMMED_DIR}/${sample}_2.fq.gz" \
        --correction \
        --cut_front 25 --cut_tail 25 --cut_mean_quality 25 \
        --trim_front1 7 --trim_front2 7 \
        -l 70 \
        -h "${QC_DIR}/${sample}_report_fastp.html" \
        -j "${QC_DIR}/${sample}_fastp.json"

    # FastQC on trimmed reads (compare with the raw-read reports from step 0)
    if command -v fastqc >/dev/null 2>&1; then
        mkdir -p "$FASTQC_TRIM_DIR"
        fastqc --threads "$THREADS_FASTQC" --outdir "$FASTQC_TRIM_DIR" \
            "${TRIMMED_DIR}/${sample}_1.fq.gz" "${TRIMMED_DIR}/${sample}_2.fq.gz"
    fi

done < "$SAMPLES"

echo "fastp finished."
