#!/usr/bin/env python3
"""
Step 9 — ClinVar pathogenic / likely pathogenic variants in disease-associated genes.

Entrada
  --variants  tabla de variantes del probando (bcftools +split-vep output of 05_clinvar_plp.sh,
              no population-frequency filter)
  --g2p       genes_to_phenotype.txt (HPO): define los genes con enfermedad asociada
              (OMIM/ORPHA/DECIPHER) y su modo de herencia
  --hpo-obo   hp.obo
  --hpo-terms términos HPO del probando
Criterios
  1. ClinVar CLNSIG exactamente Pathogenic, Likely_pathogenic o Pathogenic/Likely_pathogenic
     (se excluyen conflictos y clasificaciones combinadas como risk_factor o low_penetrance)
  2. Gen con al menos una enfermedad asociada en HPO
Salida
  --out  tabla con gen, HGVS, cigosidad, frecuencia, ClinVar (clase, estrellas, condición),
         enfermedades del gen, modo de herencia, HPO coincidentes e interpretación
"""
import argparse
import csv
import re
from collections import defaultdict

MOI = {"HP:0000006": "AD", "HP:0000007": "AR", "HP:0001417": "XL",
       "HP:0001419": "XLR", "HP:0001423": "XLD", "HP:0001427": "MT"}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def load_parents(path):
    parents, term = defaultdict(set), None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                term = None
            elif line.startswith("id: HP:"):
                term = line[4:]
            elif line.startswith("is_a: HP:") and term:
                parents[term].add(line[6:16])
            elif line.startswith("alt_id: HP:") and term:
                parents[line[8:18]].add(term)
    return parents


def ancestors(t, parents, cache):
    if t in cache:
        return cache[t]
    seen, stack = {t}, [t]
    while stack:
        for p in parents.get(stack.pop(), ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    cache[t] = seen
    return seen


def stars(revstat):
    r = (revstat or "").lower()
    if "practice_guideline" in r:
        return 4
    if "expert_panel" in r:
        return 3
    if "multiple_submitters" in r and "no_conflicts" in r:
        return 2
    if "single_submitter" in r:
        return 1
    return 0


def zygosity(gt):
    al = re.split(r"[/|]", gt)
    alt = [a for a in al if a not in ("0", ".")]
    return "hom" if alt and len(alt) == len(al) else "het"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--g2p", required=True)
    ap.add_argument("--hpo-obo", required=True)
    ap.add_argument("--hpo-terms", required=True)
    ap.add_argument("--sex", default="unknown")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    parents, cache = load_parents(a.hpo_obo), {}
    patient = {l.split()[0] for l in open(a.hpo_terms) if l.strip().startswith("HP:")}

    gene_terms, gene_moi, gene_dis = defaultdict(set), defaultdict(set), defaultdict(set)
    with open(a.g2p, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            g, h, d = row["gene_symbol"], row["hpo_id"], row.get("disease_id", "")
            if d:
                gene_dis[g].add(d)
            if h in MOI:
                gene_moi[g].add(MOI[h])
            else:
                gene_terms[g] |= ancestors(h, parents, cache)

    n_plp, out = 0, []
    with open(a.variants, newline="") as fh:
        for v in csv.DictReader(fh, delimiter="\t"):
            if v.get("ClinVar_CLNSIG", "") not in PLP:
                continue
            n_plp += 1
            g = v.get("SYMBOL", "")
            if not gene_dis.get(g):
                continue
            zyg = zygosity(v["GT"])
            modes = gene_moi.get(g, set())
            matched = sorted(patient & gene_terms.get(g, set()))
            is_x = v["CHROM"].replace("chr", "") == "X"
            if zyg == "hom" or (is_x and a.sex == "male"):
                interp = "biallelic/hemizygous: consistent with disease"
            elif modes & {"AD", "XLD"}:
                interp = "heterozygous in dominant gene: consistent with disease"
            elif modes & {"AR", "XLR", "XL"}:
                interp = "heterozygous in recessive gene: carrier unless a second allele is present"
            else:
                interp = "heterozygous; inheritance not annotated"
            out.append({
                "gene": g, "chrom": v["CHROM"], "pos": v["POS"], "ref": v["REF"], "alt": v["ALT"],
                "HGVSc": v.get("HGVSc", ""), "HGVSp": v.get("HGVSp", ""), "consequence": v.get("Consequence", ""),
                "zygosity": zyg, "GT": v["GT"], "DP": v.get("DP", ""), "GQ": v.get("GQ", ""),
                "MAX_AF": v.get("MAX_AF", ""),
                "ClinVar_class": v["ClinVar_CLNSIG"], "ClinVar_stars": stars(v.get("ClinVar_CLNREVSTAT")),
                "ClinVar_condition": v.get("ClinVar_CLNDN", ""),
                "gene_diseases": ";".join(sorted(gene_dis[g])), "MOI": "/".join(sorted(modes)) or "unknown",
                "n_hpo_matched": len(matched), "matched_hpo": ";".join(matched),
                "interpretation": interp,
            })

    out.sort(key=lambda r: (-r["n_hpo_matched"], -r["ClinVar_stars"], r["gene"]))
    with open(a.out, "w", newline="") as fh:
        cols = list(out[0].keys()) if out else ["gene"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"ClinVar P/LP en el probando          : {n_plp}")
    print(f"... en genes asociados a enfermedad  : {len(out)}")
    print(f"... con coincidencia de HPO          : {sum(r['n_hpo_matched'] > 0 for r in out)}")
    for k in sorted({r['interpretation'] for r in out}):
        print(f"  {sum(r['interpretation'] == k for r in out):4d}  {k}")


if __name__ == "__main__":
    main()
