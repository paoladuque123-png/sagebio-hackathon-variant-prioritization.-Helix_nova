#!/bin/bash
# Print versions of all tools used, for the methods description.
# Usage: bash scripts/report_versions.sh [config/config.sh] > versions.txt
CONFIG="${1:-config/config.sh}"
# shellcheck source=/dev/null
source "$CONFIG"
# shellcheck source=/dev/null
source "$CONDA_SH"

[[ -n "$ENV_QC" ]] && conda activate "$ENV_QC"
echo "fastqc:   $(fastqc --version 2>&1 | head -n 1)"
command -v multiqc >/dev/null 2>&1 && echo "multiqc:  $(multiqc --version 2>&1 | head -n 1)"

[[ -n "$ENV_FASTP" ]] && conda activate "$ENV_FASTP"
echo "fastp:    $(fastp --version 2>&1 | head -n 1)"

conda activate "$ENV_ALIGN"
echo "bwa:      $(bwa 2>&1 | grep -m1 Version)"
echo "samtools: $(samtools --version | head -n 1)"
echo "tabix:    $(tabix --version 2>&1 | head -n 1)"

conda activate "$ENV_GATK"
echo "gatk:     $(gatk --version 2>&1 | grep -m1 'The Genome Analysis Toolkit')"
echo "picard:   $(picard MarkDuplicates --version 2>&1 | tail -n 1)"
echo "java:     $(java -version 2>&1 | head -n 1)"

conda activate "$ENV_VEP"
echo "vep:      $(vep 2>&1 | grep -m1 'ensembl-vep')"
echo "bcftools: $(bcftools --version | head -n 1)"

echo "python:   $(python --version 2>&1)"
