#!/bin/bash
# =============================================================
# One-time setup of VEP annotation resources (requires internet on the server)
#  - Ensembl VEP cache (GRCh38; skipped if already present)
#  - VEP plugin modules (.pm) into VEP_PLUGINS_DIR
#  - AlphaMissense data file (public) + tabix index
#  - ClinVar GRCh38 VCF + tabix index
#  - HPO ontology (hp.obo) + gene-phenotype annotations (genes_to_phenotype.txt)
# Data for the remaining plugins (REVEL, ClinPred, CADD, SpliceAI, LOEUF,
# DosageSensitivity, UTRAnnotator) is downloaded manually; see README, "VEP plugin data".
# Usage : bash scripts/setup_vep_resources.sh [config/config.sh]
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
# shellcheck source=/dev/null
source "$CONFIG"
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$ENV_VEP"

mkdir -p "$VEP_CACHE_DIR" "$VEP_PLUGINS_DIR" "$VEP_PLUGIN_DATA" "$(dirname "$CLINVAR_VCF")"

# VEP cache (large download; --CONVERT builds tabix-indexed cache for faster lookups)
case "${VEP_TRANSCRIPT_SET:-ensembl}" in
    merged)  CACHE_SPECIES="homo_sapiens_merged" ;;
    refseq)  CACHE_SPECIES="homo_sapiens_refseq" ;;
    *)       CACHE_SPECIES="homo_sapiens" ;;
esac
if [[ -d "${VEP_CACHE_DIR}/${CACHE_SPECIES}/${VEP_CACHE_VERSION}_GRCh38" ]]; then
    echo "VEP cache already present: ${VEP_CACHE_DIR}/${CACHE_SPECIES}/${VEP_CACHE_VERSION}_GRCh38"
else
    vep_install -a c -s "$CACHE_SPECIES" -y GRCh38 -c "$VEP_CACHE_DIR" \
        --CACHE_VERSION "$VEP_CACHE_VERSION" --CONVERT --NO_UPDATE
fi

# VEP plugins (.pm files). Data files for each plugin are downloaded separately:
# the download and indexing instructions are in the header of each .pm file.
vep_install -a p -c "$VEP_CACHE_DIR" --PLUGINSDIR "$VEP_PLUGINS_DIR" \
    --PLUGINS NMD,AlphaMissense,REVEL,ClinPred,CADD,SpliceAI,LOEUF,DosageSensitivity,UTRAnnotator \
    --NO_UPDATE

# AlphaMissense (public, Google DeepMind)
cd "$VEP_PLUGIN_DATA"
wget -N https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -f -s 1 -b 2 -e 2 -S 1 AlphaMissense_hg38.tsv.gz

# ClinVar (record the download date in docs/prioritization.md)
cd "$(dirname "$CLINVAR_VCF")"
wget -N https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -N https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
echo "ClinVar file date: $(bcftools view -h clinvar.vcf.gz | grep -m1 '^##fileDate')"

# HPO ontology and gene-phenotype annotations (for prioritization; record the release)
mkdir -p "$(dirname "$HPO_OBO")" && cd "$(dirname "$HPO_OBO")"
wget -N https://purl.obolibrary.org/obo/hp.obo
wget -N https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt
echo "HPO release: $(grep -m1 '^data-version' hp.obo)"

echo "Annotation resources ready."
