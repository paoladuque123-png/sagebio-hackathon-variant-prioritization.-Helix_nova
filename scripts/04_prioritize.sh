#!/bin/bash
# =============================================================
# Steps 8–9 — Prioritization and submission file
# Input : <ANNOT_DIR>/<ID>.candidates.tsv, samples/<ID>.hpo.txt, HPO ontology + annotations
# Output: <ANNOT_DIR>/<ID>.ranked.tsv
#         submission/<SUBMISSION_USER>_<APPROACH_NAME>.csv
# Usage : bash scripts/04_prioritize.sh [config/config.sh]
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
REPO_DIR="$(dirname "$SCRIPT_DIR")"
mkdir -p "${REPO_DIR}/submission"

[[ -f "$HPO_OBO" ]] || { echo "ERROR: hp.obo not found: $HPO_OBO"; exit 1; }
[[ -f "$HPO_G2P" ]] || { echo "ERROR: genes_to_phenotype.txt not found: $HPO_G2P"; exit 1; }
[[ "$SUBMISSION_USER" != "hf-username" && -n "$SUBMISSION_USER" ]] \
    || { echo "ERROR: set SUBMISSION_USER (Hugging Face username or team) in $CONFIG"; exit 1; }

while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue
    hpo_file="samples/${sample}.hpo.txt"
    [[ -f "$hpo_file" ]] || { echo "ERROR: HPO term file not found: $hpo_file"; exit 1; }

    echo "HPO release: $(grep -m1 '^data-version' "$HPO_OBO")"
    echo "== HPO terms checked against hp.obo =="
    grep '^HP:' "$hpo_file" | while read -r hp rest; do
        name=$(awk -v id="id: $hp" '$0==id{getline; sub(/^name: /,""); print; exit}' "$HPO_OBO")
        echo "$hp  given: $rest  |  hp.obo: ${name:-NOT FOUND}"
    done

    python "${SCRIPT_DIR}/04_prioritize.py" \
        --proband-id "$sample" \
        --candidates "${ANNOT_DIR}/${sample}.candidates.tsv" \
        --hpo-terms "$hpo_file" \
        --hpo-obo "$HPO_OBO" \
        --g2p "$HPO_G2P" \
        --sex "${PROBAND_SEX:-unknown}" \
        --max-rows "${MAX_SUBMISSION_ROWS:-10}" \
        --out-ranked "${ANNOT_DIR}/${sample}.ranked.tsv" \
        --out-submission "${REPO_DIR}/submission/${SUBMISSION_USER}_${APPROACH_NAME}.csv"
done < "$SAMPLES"
