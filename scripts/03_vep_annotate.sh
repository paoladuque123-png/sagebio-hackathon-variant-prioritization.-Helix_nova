#!/bin/bash
# =============================================================
# Step 7 — Normalization and local variant annotation with Ensembl VEP (offline cache)
# Input : <VCF_DIR>/<ID>.vcf.gz  (GATK HaplotypeCaller, GRCh38)
# Output: <ANNOT_DIR>/<ID>.norm.vcf.gz         (normalized, biallelic)
#         <ANNOT_DIR>/<ID>.vep.vcf.gz          (VEP-annotated VCF)
#         <ANNOT_DIR>/<ID>.vep_summary.html    (VEP run summary)
#         <ANNOT_DIR>/<ID>.vep.tsv             (one row per variant, most severe consequence)
#         <ANNOT_DIR>/<ID>.candidates.tsv      (quality + rarity + impact pre-filter)
#         <ANNOT_DIR>/<ID>.annotation_stats.txt (counts only, for the report)
# All processing is local; no variant data leaves the server.
# Usage : bash scripts/03_vep_annotate.sh [config/config.sh]
# =============================================================
set -eo pipefail

CONFIG="${1:-config/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$ENV_VEP"

die() { echo "ERROR: $1"; exit 1; }

[[ -f "$SAMPLES" ]]          || die "sample list not found: $SAMPLES"
[[ -f "$REF" ]]              || die "reference not found: $REF"
[[ -d "$VEP_CACHE_DIR" ]]    || die "VEP cache dir not found: $VEP_CACHE_DIR"
[[ -f "$CLINVAR_VCF" ]]      || die "ClinVar VCF not found: $CLINVAR_VCF"
[[ -f "${CLINVAR_VCF}.tbi" ]] || die "ClinVar index missing: tabix -p vcf $CLINVAR_VCF"
mkdir -p "$ANNOT_DIR"

# Optional plugins: each is enabled only if its data file(s) are configured and present
PLUGIN_ARGS=()
have() { for f in "$@"; do [[ -n "$f" && -f "$f" ]] || return 1; done; }

[[ "$USE_NMD" == "1" ]] && PLUGIN_ARGS+=(--plugin NMD)
have "$ALPHAMISSENSE_FILE"         && PLUGIN_ARGS+=(--plugin "AlphaMissense,file=${ALPHAMISSENSE_FILE}")
have "$REVEL_FILE"                 && PLUGIN_ARGS+=(--plugin "REVEL,file=${REVEL_FILE}")
have "$CLINPRED_FILE"              && PLUGIN_ARGS+=(--plugin "ClinPred,file=${CLINPRED_FILE}")
have "$CADD_SNV" "$CADD_INDEL"     && PLUGIN_ARGS+=(--plugin "CADD,snv=${CADD_SNV},indels=${CADD_INDEL}")
have "$SPLICEAI_SNV" "$SPLICEAI_INDEL" && PLUGIN_ARGS+=(--plugin "SpliceAI,snv=${SPLICEAI_SNV},indel=${SPLICEAI_INDEL}")
have "$LOEUF_FILE"                 && PLUGIN_ARGS+=(--plugin "LOEUF,file=${LOEUF_FILE},match_by=transcript")
have "$DOSAGE_FILE"                && PLUGIN_ARGS+=(--plugin "DosageSensitivity,file=${DOSAGE_FILE}")
have "$UTRANNOTATOR_FILE"          && PLUGIN_ARGS+=(--plugin "UTRAnnotator,file=${UTRANNOTATOR_FILE}")

if [[ ${#PLUGIN_ARGS[@]} -gt 0 ]]; then
    [[ -d "$VEP_PLUGINS_DIR" ]] || die "plugins configured but VEP_PLUGINS_DIR not found: $VEP_PLUGINS_DIR"
    PLUGIN_ARGS=(--dir_plugins "$VEP_PLUGINS_DIR" "${PLUGIN_ARGS[@]}")
    echo "Plugins enabled: $(printf '%s\n' "${PLUGIN_ARGS[@]}" | grep -v -- '--' | cut -d, -f1 | grep -v "$VEP_PLUGINS_DIR" | paste -sd ' ')"
fi

# Transcript set: ensembl (default), merged (Ensembl/GENCODE + RefSeq) or refseq
case "${VEP_TRANSCRIPT_SET:-ensembl}" in
    merged)  TRANSCRIPT_ARGS=(--merged) ;;
    refseq)  TRANSCRIPT_ARGS=(--refseq) ;;
    ensembl) TRANSCRIPT_ARGS=() ;;
    *) die "VEP_TRANSCRIPT_SET must be merged, ensembl or refseq" ;;
esac

while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue
    in_vcf="${VCF_DIR}/${sample}.vcf.gz"
    [[ -f "$in_vcf" ]] || die "VCF not found: $in_vcf"
    vcf_chr=$(bcftools view -H "$in_vcf" 2>/dev/null | head -n 1 | cut -f1 || true)
    grep -q -P "^${vcf_chr}\t" "${REF}.fai" \
        || die "contig '${vcf_chr}' of the VCF is not in ${REF}.fai (reference mismatch)"

    echo "============================================================"
    echo "VEP annotation: $sample"
    echo "============================================================"

    # 1. Normalize: split multiallelic sites and left-align indels
    echo "[1/4] bcftools norm"
    bcftools norm -m -any -f "$REF" -Oz -o "${ANNOT_DIR}/${sample}.norm.vcf.gz" "$in_vcf"
    tabix -f -p vcf "${ANNOT_DIR}/${sample}.norm.vcf.gz"

    # 2. VEP annotation (offline, Ensembl cache, ClinVar as custom annotation)
    echo "[2/4] VEP"
    vep \
        --input_file "${ANNOT_DIR}/${sample}.norm.vcf.gz" \
        --output_file "${ANNOT_DIR}/${sample}.vep.vcf.gz" \
        --stats_file "${ANNOT_DIR}/${sample}.vep_summary.html" \
        --vcf --compress_output bgzip --force_overwrite \
        --offline --cache --dir_cache "$VEP_CACHE_DIR" --cache_version "$VEP_CACHE_VERSION" \
        "${TRANSCRIPT_ARGS[@]}" \
        --species homo_sapiens --assembly GRCh38 \
        --fasta "$REF" \
        --fork "$THREADS_VEP" \
        --symbol --transcript_version --hgvs --hgvsg --canonical --mane --biotype --numbers --protein --uniprot \
        --variant_class --check_existing --pubmed --distance 5000 \
        --af --af_1kg --af_gnomade --af_gnomadg --max_af \
        --sift b --polyphen b \
        --flag_pick \
        --custom "file=${CLINVAR_VCF},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" \
        "${PLUGIN_ARGS[@]}"
    tabix -f -p vcf "${ANNOT_DIR}/${sample}.vep.vcf.gz"

    # 3. Flatten to a table (most severe consequence per variant), keeping genotype and quality
    echo "[3/4] bcftools +split-vep -> TSV"
    CSQ_FIELDS=$(bcftools +split-vep -l "${ANNOT_DIR}/${sample}.vep.vcf.gz" | cut -f2 | paste -sd ',')
    WANT="SYMBOL,Gene,Feature,MANE_SELECT,CANONICAL,BIOTYPE,Consequence,IMPACT,EXON,INTRON,HGVSc,HGVSp,HGVSg,Existing_variation,gnomADe_AF,gnomADg_AF,MAX_AF,SIFT,PolyPhen,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN"
    for opt in NMD am_class am_pathogenicity REVEL ClinPred CADD_PHRED \
               SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL SpliceAI_pred_SYMBOL \
               LOEUF DosageSensitivity_pHaplo DosageSensitivity_pTriplo \
               5UTR_consequence 5UTR_annotation PUBMED; do
        [[ ",${CSQ_FIELDS}," == *",${opt},"* ]] && WANT="${WANT},${opt}"
    done
    FMT="%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t[%GT\t%DP\t%GQ\t%AD]"
    for f in ${WANT//,/ }; do FMT="${FMT}\t%${f}"; done
    {
        echo -e "CHROM\tPOS\tREF\tALT\tFILTER\tQUAL\tGT\tDP\tGQ\tAD\t${WANT//,/\\t}"
        bcftools +split-vep -s worst -d -f "${FMT}\n" "${ANNOT_DIR}/${sample}.vep.vcf.gz"
    } > "${ANNOT_DIR}/${sample}.vep.tsv"

    # 4. Pre-filter candidates: non-reference genotype, basic quality, rare, and
    #    protein-altering/splice or ClinVar P/LP. Thresholds are set in config.
    echo "[4/4] Candidate pre-filter"
    bcftools view -i "GT!=\"ref\" && GT!=\"mis\" && FMT/DP>=${MIN_DP} && FMT/GQ>=${MIN_GQ}" \
        "${ANNOT_DIR}/${sample}.vep.vcf.gz" -Ou \
    | bcftools +split-vep -s worst -d -f "${FMT}\n" \
    | { echo -e "CHROM\tPOS\tREF\tALT\tFILTER\tQUAL\tGT\tDP\tGQ\tAD\t${WANT//,/\\t}"; cat; } \
    | awk -F'\t' -v OFS='\t' -v maxaf="$MAX_POP_AF" '
        NR==1 { for (i=1; i<=NF; i++) c[$i]=i; print; next }
        {
            af = $c["MAX_AF"]
            rare = (af == "." || af == "" || af + 0 < maxaf)
            impact = ($c["IMPACT"] == "HIGH" || $c["IMPACT"] == "MODERATE" || $c["Consequence"] ~ /splice/)
            clinvar = ($c["ClinVar_CLNSIG"] ~ /athogenic/ && $c["ClinVar_CLNSIG"] !~ /onflicting/)
            if (rare && (impact || clinvar)) print
        }' \
        > "${ANNOT_DIR}/${sample}.candidates.tsv"

    # 5. Run summary (counts only)
    {
        echo "Date                  : $(date +%F)"
        echo "VEP / cache           : $(vep 2>&1 | grep -m1 'ensembl-vep' | awk '{print $NF}') / ${VEP_CACHE_VERSION}"
        echo "ClinVar               : $(bcftools view -h "$CLINVAR_VCF" | grep -m1 '^##fileDate' | cut -d= -f2)"
        echo "Records in VCF        : $(bcftools view -H "$in_vcf" | wc -l)"
        echo "After normalization   : $(bcftools view -H "${ANNOT_DIR}/${sample}.norm.vcf.gz" | wc -l)"
        echo "Candidates            : $(($(wc -l < "${ANNOT_DIR}/${sample}.candidates.tsv") - 1))"
        echo "Candidates by IMPACT  :"
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}{n[$c["IMPACT"]]++} END{for(k in n) print "  "k"\t"n[k]}' \
            "${ANNOT_DIR}/${sample}.candidates.tsv" | sort
    } | tee "${ANNOT_DIR}/${sample}.annotation_stats.txt"
done < "$SAMPLES"

echo "VEP annotation finished."
