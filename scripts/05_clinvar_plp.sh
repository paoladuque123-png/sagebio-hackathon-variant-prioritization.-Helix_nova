#!/bin/bash
# =============================================================
# Step 9 — ClinVar P/LP variants in disease-associated genes (independent of the ranking)
# Recorre el VCF anotado completo (sin filtro de frecuencia) con control de calidad.
# Input : <ANNOT_DIR>/<ID>.vep.vcf.gz
# Output: <ANNOT_DIR>/clinvar_plp/<ID>.clinvar_plp_disease_genes.tsv
# Usage : bash scripts/05_clinvar_plp.sh [config/config.sh]
# =============================================================
set -eo pipefail
CONFIG="${1:-config/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$ENV_VEP"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ID="$(grep -m1 -v '^$' "$SAMPLES")"
SEX="${PROBAND_SEX:-unknown}"
VEP_VCF="${ANNOT_DIR}/${ID}.vep.vcf.gz"
HPO_DIR="$(dirname "$HPO_G2P")"
HPO_TERMS="samples/${ID}.hpo.txt"
OUT="${ANNOT_DIR}/clinvar_plp"

for f in "$VEP_VCF" "$HPO_OBO" "$HPO_G2P" "$HPO_TERMS" "${SCRIPT_DIR}/05_clinvar_plp.py"; do
    [[ -f "$f" ]] || { echo "ERROR: falta $f"; exit 1; }
done
mkdir -p "$OUT"

FIELDS="SYMBOL,Consequence,HGVSc,HGVSp,MAX_AF,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN"
FMT="%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%GQ]"
for f in ${FIELDS//,/ }; do FMT="${FMT}\t%${f}"; done

echo "[1/2] Extracción de variantes con genotipo alternativo (DP>=${MIN_DP}, GQ>=${MIN_GQ})"
{
  echo -e "CHROM\tPOS\tREF\tALT\tGT\tDP\tGQ\t${FIELDS//,/\\t}"
  bcftools view -i "GT=\"alt\" && FMT/DP>=${MIN_DP} && FMT/GQ>=${MIN_GQ}" "$VEP_VCF" -Ou \
  | bcftools +split-vep -s worst -d -f "${FMT}\n" 2>/dev/null \
  | awk -F'\t' '$13 ~ /athogenic/'
} > "$OUT/${ID}.clinvar_any_plp.tmp.tsv"

echo "[2/2] Filtro P/LP en genes asociados a enfermedad"
python "${SCRIPT_DIR}/05_clinvar_plp.py" \
  --variants "$OUT/${ID}.clinvar_any_plp.tmp.tsv" \
  --g2p "$HPO_G2P" \
  --hpo-obo "$HPO_OBO" \
  --hpo-terms "$HPO_TERMS" \
  --sex "$SEX" \
  --out "$OUT/${ID}.clinvar_plp_disease_genes.tsv"
rm -f "$OUT/${ID}.clinvar_any_plp.tmp.tsv"

echo
echo "== Hallazgos (gen, HGVSc, HGVSp, cigosidad, MAX_AF, ClinVar, estrellas, MOI, HPO, interpretación) =="
awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  {print $c["gene"], $c["HGVSc"], $c["HGVSp"], $c["zygosity"], $c["MAX_AF"], $c["ClinVar_class"],
         $c["ClinVar_stars"]"*", $c["MOI"], $c["n_hpo_matched"], $c["interpretation"]}' \
  "$OUT/${ID}.clinvar_plp_disease_genes.tsv"
