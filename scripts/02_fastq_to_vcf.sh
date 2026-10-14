#!/bin/bash
# =============================================================
# Step 2 — Alignment and germline variant calling (GATK Best Practices)
# BWA-MEM -> samtools sort -> Picard MarkDuplicates
# -> Picard AddOrReplaceReadGroups -> GATK BaseRecalibrator
# -> GATK ApplyBQSR -> GATK HaplotypeCaller
# Reference: GRCh38 (GCA_000001405.15, no-alt analysis set)
# Usage : bash scripts/02_fastq_to_vcf.sh [config/config.sh]
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
# shellcheck source=/dev/null
source "$CONDA_SH"

mkdir -p "$BAM_DIR" "$BAM_RG_DIR" "$RECAL_DIR" "$VCF_DIR" "$LOG_DIR" "$TMP_DIR"

# -------------------------------------------------------------
# Input checks
# -------------------------------------------------------------
die() { echo "ERROR: $1"; [[ -n "$2" ]] && echo "Run: $2"; exit 1; }

DICT="${REF%.*}.dict"

[[ -f "$SAMPLES" ]]              || die "sample list not found: $SAMPLES"
[[ -d "$TRIMMED_DIR" ]]          || die "trimmed FASTQ dir not found: $TRIMMED_DIR"
[[ -f "$REF" ]]                  || die "reference not found: $REF"
[[ -f "${REF}.bwt" ]]            || die "BWA index missing" "bwa index $REF"
[[ -f "${REF}.fai" ]]            || die "FASTA index missing" "samtools faidx $REF"
[[ -f "$DICT" ]]                 || die "sequence dictionary missing" "picard CreateSequenceDictionary R=$REF O=$DICT"
[[ -f "$KNOWN_SITES_1" ]]        || die "known-sites not found: $KNOWN_SITES_1"
[[ -f "${KNOWN_SITES_1}.tbi" ]]  || die "tabix index missing" "tabix -p vcf $KNOWN_SITES_1"
[[ -f "$KNOWN_SITES_2" ]]        || die "known-sites not found: $KNOWN_SITES_2"
[[ -f "${KNOWN_SITES_2}.idx" ]]  || die "GATK index missing" "gatk IndexFeatureFile -I $KNOWN_SITES_2"

echo "Input checks passed."

# -------------------------------------------------------------
# Per-sample processing
# -------------------------------------------------------------
while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue

    echo "============================================================"
    echo "Sample: $sample"
    echo "============================================================"

    r1=$(find "$TRIMMED_DIR" -type f -name "${sample}*_1.fq.gz" | head -n 1)
    r2=$(find "$TRIMMED_DIR" -type f -name "${sample}*_2.fq.gz" | head -n 1)
    [[ -n "$r1" && -n "$r2" ]] || die "paired trimmed FASTQ not found for $sample"

    # 1. Alignment + coordinate sort
    echo "[1/6] BWA-MEM + samtools sort"
    conda deactivate || true
    conda activate "$ENV_ALIGN"
    bwa mem -t "$THREADS_BWA" "$REF" "$r1" "$r2" \
        | samtools sort -@ "$THREADS_BWA" -o "${BAM_DIR}/${sample}.bam" -
    samtools index "${BAM_DIR}/${sample}.bam"

    # 2. Duplicate marking
    echo "[2/6] Picard MarkDuplicates"
    conda deactivate || true
    conda activate "$ENV_GATK"
    picard -Xms10g -Xmx40g MarkDuplicates \
        I="${BAM_DIR}/${sample}.bam" \
        O="${BAM_DIR}/${sample}.dedup.bam" \
        M="${BAM_DIR}/${sample}.metrics.txt" \
        CREATE_INDEX=true \
        TMP_DIR="$TMP_DIR"

    # 3. Read groups
    echo "[3/6] Picard AddOrReplaceReadGroups"
    picard -Xms10g -Xmx30g AddOrReplaceReadGroups \
        I="${BAM_DIR}/${sample}.dedup.bam" \
        O="${BAM_RG_DIR}/${sample}.name.bam" \
        RGLB="$RG_LB" \
        RGPL="$RG_PL" \
        RGPU="$RG_PU" \
        RGSM="$sample" \
        CREATE_INDEX=true \
        TMP_DIR="$TMP_DIR"

    # 4. Base quality score recalibration — model
    echo "[4/6] GATK BaseRecalibrator"
    gatk --java-options "-Xmx60g -Xms40g -Djava.io.tmpdir=${TMP_DIR}" BaseRecalibrator \
        -I "${BAM_RG_DIR}/${sample}.name.bam" \
        -R "$REF" \
        --known-sites "$KNOWN_SITES_1" \
        --known-sites "$KNOWN_SITES_2" \
        -O "${RECAL_DIR}/${sample}_table"

    # 5. Base quality score recalibration — apply
    echo "[5/6] GATK ApplyBQSR"
    gatk --java-options "-Xmx30g -Djava.io.tmpdir=${TMP_DIR}" ApplyBQSR \
        -R "$REF" \
        -I "${BAM_RG_DIR}/${sample}.name.bam" \
        --bqsr-recal-file "${RECAL_DIR}/${sample}_table" \
        -O "${RECAL_DIR}/${sample}.recal.bam"

    conda deactivate || true
    conda activate "$ENV_ALIGN"
    samtools index "${RECAL_DIR}/${sample}.recal.bam"

    # 6. Germline variant calling
    echo "[6/6] GATK HaplotypeCaller"
    conda deactivate || true
    conda activate "$ENV_GATK"
    gatk --java-options "-Xmx30g -Xms15g -Djava.io.tmpdir=${TMP_DIR}" HaplotypeCaller \
        -R "$REF" \
        -I "${RECAL_DIR}/${sample}.recal.bam" \
        -O "${VCF_DIR}/${sample}.vcf.gz"

    echo "Done: $sample"
done < "$SAMPLES"

echo "Pipeline finished."
