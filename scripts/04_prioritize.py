#!/usr/bin/env python3
"""
Step 8-9 — Phenotype- and evidence-based prioritization of candidate variants and
generation of the submission CSV.

Input
  --candidates  <ID>.candidates.tsv from 03_vep_annotate.sh
  --hpo-terms   text file with the proband's HPO terms (one HP:XXXXXXX per line)
  --hpo-obo     hp.obo (HPO ontology)
  --g2p         genes_to_phenotype.txt (HPO gene-phenotype annotations)
Output
  --out-ranked      full ranked table (all candidates and compound-het pairs), with scores
  --out-submission  submission CSV (top N rows, Hackathon format)

Scoring (all components in [0, 1])
  variant_score  evidence from ClinVar, predicted consequence, in silico predictors,
                 splicing predictions, gene constraint and population frequency
  pheno_score    fraction of the proband's HPO terms annotated to the gene,
                 counting annotations to the term or any of its descendants
  inheritance    consistency of the genotype with the gene's known mode(s) of inheritance
  EPCR = variant_score * (W_FLOOR + (1 - W_FLOOR) * pheno_score) * inheritance

Only the Python standard library is used.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

W_FLOOR = 0.3                 # weight kept by a variant whose gene has no phenotype match
LOF_TERMS = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
             "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation")
MOI = {                       # HPO mode-of-inheritance terms
    "HP:0000006": "AD", "HP:0000007": "AR", "HP:0001417": "XL",
    "HP:0001419": "XLR", "HP:0001423": "XLD", "HP:0001427": "MT",
}
SUBMISSION_FIELDS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                     "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]


# ---------------------------------------------------------------- helpers
def fnum(x, default=None):
    """First numeric value in a VEP field ('0.85', '0.85&0.10', 'deleterious(0.01)')."""
    if x is None:
        return default
    vals = [float(v) for v in re.findall(r"-?\d+\.?\d*(?:[eE]-?\d+)?", str(x))]
    return max(vals) if vals else default


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_hpo_parents(obo_path):
    parents, term = defaultdict(set), None
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                term = None
            elif line.startswith("id: HP:"):
                term = line[4:]
            elif line.startswith("is_a: HP:") and term:
                parents[term].add(line[6:16])
            elif line.startswith("alt_id: HP:") and term:
                parents[line[8:18]].add(term)          # map obsolete ids to current term
    return parents


def ancestors(term, parents, cache):
    if term in cache:
        return cache[term]
    seen, stack = {term}, [term]
    while stack:
        for p in parents.get(stack.pop(), ()):
            if p not in seen:
                seen.add(p)
                stack.append(p)
    cache[term] = seen
    return seen


def load_gene_annotations(g2p_path, parents):
    """gene symbol -> (set of HPO terms incl. ancestors, set of modes of inheritance)."""
    terms, moi, cache = defaultdict(set), defaultdict(set), {}
    with open(g2p_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sym_col = "gene_symbol" if "gene_symbol" in reader.fieldnames else reader.fieldnames[1]
        hpo_col = "hpo_id" if "hpo_id" in reader.fieldnames else reader.fieldnames[2]
        for row in reader:
            gene, hpo = row[sym_col], row[hpo_col]
            if hpo in MOI:
                moi[gene].add(MOI[hpo])
            else:
                terms[gene] |= ancestors(hpo, parents, cache)
    return terms, moi


# ---------------------------------------------------------------- scoring
def variant_evidence(v):
    """Return (score, list of evidence labels) for one variant."""
    csq = v.get("Consequence", "")
    impact = v.get("IMPACT", "")
    clnsig = v.get("ClinVar_CLNSIG", "")
    ev, score = [], 0.2

    if re.search(r"athogenic", clnsig) and not re.search(r"onflicting", clnsig):
        score = max(score, 1.0 if re.search(r"(^|[_&|,])Pathogenic", clnsig) else 0.9)
        ev.append(f"ClinVar {clnsig}")

    if any(t in csq for t in LOF_TERMS):
        loeuf = fnum(v.get("LOEUF"))
        nmd = v.get("NMD", "")
        s = 0.85
        if loeuf is not None and loeuf < 0.6:
            s = 0.95
            ev.append(f"LoF in constrained gene (LOEUF {loeuf:.2f})")
        else:
            ev.append("predicted LoF")
        if "NMD_escaping" in nmd:
            s -= 0.15
            ev.append("NMD-escaping")
        score = max(score, s)

    am = v.get("am_class", "")
    revel = fnum(v.get("REVEL"))
    clinpred = fnum(v.get("ClinPred"))
    if "missense_variant" in csq:
        hits = []
        if "likely_pathogenic" in am:
            hits.append("AlphaMissense LP")
        if revel is not None and revel >= 0.75:
            hits.append(f"REVEL {revel:.2f}")
        if clinpred is not None and clinpred >= 0.5:
            hits.append(f"ClinPred {clinpred:.2f}")
        s = {0: 0.35, 1: 0.6, 2: 0.7}.get(len(hits), 0.8)
        score = max(score, s)
        ev += hits

    ds = [fnum(v.get(k)) for k in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                                   "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")]
    ds = max([d for d in ds if d is not None], default=None)
    if ds is not None and ds >= 0.2:
        score = max(score, 0.8 if ds >= 0.5 else 0.5)
        ev.append(f"SpliceAI {ds:.2f}")

    if impact == "MODERATE":
        score = max(score, 0.35)

    max_af = fnum(v.get("MAX_AF"))
    if max_af is None:
        score = min(1.0, score + 0.05)
        ev.append("absent from gnomAD")
    elif max_af > 0.001:
        score *= 0.9

    return round(score, 4), ev


def zygosity(gt):
    alleles = re.split(r"[/|]", gt)
    alt = [a for a in alleles if a not in ("0", ".")]
    if len(alt) == 0:
        return "ref"
    if len(alt) == len(alleles):
        return "hom"
    return "het"


def inheritance_factor(modes, zyg, chrom, sex, n_qualifying_het):
    is_x = chrom.replace("chr", "") == "X"
    if not modes:
        return 0.8, "inheritance unknown"
    if zyg == "hom" or (is_x and sex == "male"):
        return 1.0, f"{'hemizygous' if is_x and sex == 'male' else 'homozygous'}; MOI {'/'.join(sorted(modes))}"
    # heterozygous
    if modes & {"AD", "XLD"}:
        return 1.0, f"heterozygous; MOI {'/'.join(sorted(modes))}"
    if modes & {"AR", "XLR", "XL"}:
        if n_qualifying_het >= 2:
            return 1.0, "candidate compound heterozygous (phase not determined)"
        return 0.35, "single heterozygous variant in recessive gene"
    return 0.8, f"MOI {'/'.join(sorted(modes))}"


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proband-id", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--hpo-terms", required=True)
    ap.add_argument("--hpo-obo", required=True)
    ap.add_argument("--g2p", required=True)
    ap.add_argument("--sex", choices=["male", "female", "unknown"], default="unknown")
    ap.add_argument("--max-rows", type=int, default=10)
    ap.add_argument("--min-het-score", type=float, default=0.5,
                    help="minimum variant_score for a heterozygous variant to count toward a compound-het pair")
    ap.add_argument("--out-ranked", required=True)
    ap.add_argument("--out-submission", required=True)
    a = ap.parse_args()

    parents = load_hpo_parents(a.hpo_obo)
    patient_terms = {l.split()[0] for l in open(a.hpo_terms) if l.strip().startswith("HP:")}
    if not patient_terms:
        sys.exit("ERROR: no HPO terms found in " + a.hpo_terms)
    gene_terms, gene_moi = load_gene_annotations(a.g2p, parents)

    variants = read_tsv(a.candidates)
    for v in variants:
        v["zyg"] = zygosity(v["GT"])
        v["variant_score"], v["evidence"] = variant_evidence(v)
        gene = v.get("SYMBOL", "")
        matched = patient_terms & gene_terms.get(gene, set())
        v["pheno_score"] = round(len(matched) / len(patient_terms), 4)
        v["matched_hpo"] = ";".join(sorted(matched))

    # qualifying heterozygous variants per gene (for compound-het pairing)
    het_by_gene = defaultdict(list)
    for v in variants:
        if v["zyg"] == "het" and v["variant_score"] >= a.min_het_score and v.get("SYMBOL"):
            het_by_gene[v["SYMBOL"]].append(v)

    # compound-het pairs: best two qualifying het variants per recessive-compatible gene.
    # The pair is scored as one finding (geometric mean of the two variant scores); its
    # members are kept in the ranked table but not submitted as separate rows.
    rows, pair_genes = [], set()
    for gene, hets in het_by_gene.items():
        modes = gene_moi.get(gene, set())
        if len(hets) < 2 or (modes and not modes & {"AR", "XLR", "XL"}) or (modes & {"AD", "XLD"}):
            continue
        h1, h2 = sorted(hets, key=lambda x: -x["variant_score"])[:2]
        vs = (h1["variant_score"] * h2["variant_score"]) ** 0.5
        factor = 1.0 if modes else 0.8
        note = "candidate compound heterozygous (phase not determined)"
        epcr = vs * (W_FLOOR + (1 - W_FLOOR) * h1["pheno_score"]) * factor
        rows.append({"type": "compound_het", "gene": gene, "v1": h1, "v2": h2, "epcr": epcr,
                     "vscore": round(vs, 4), "inheritance": note, "submit": True})
        pair_genes.add(gene)

    for v in variants:
        gene = v.get("SYMBOL", "")
        modes = gene_moi.get(gene, set())
        if gene in pair_genes and v["zyg"] == "het":
            inh, inh_note, submit = 0.35, "member of compound-het pair (see pair row)", False
        else:
            inh, inh_note = inheritance_factor(modes, v["zyg"], v["CHROM"], a.sex, 0)
            submit = True
        epcr = v["variant_score"] * (W_FLOOR + (1 - W_FLOOR) * v["pheno_score"]) * inh
        rows.append({"type": "single", "gene": gene, "v1": v, "v2": None, "epcr": epcr,
                     "vscore": v["variant_score"], "inheritance": inh_note, "submit": submit})

    for r in rows:
        r["epcr"] = round(min(1.0, max(r["epcr"], 0.001)), 4)
    rows.sort(key=lambda r: -r["epcr"])

    # full ranked table
    with open(a.out_ranked, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["rank", "type", "gene", "variant_1", "variant_2", "HGVSc_1", "HGVSp_1",
                    "zygosity_1", "variant_score", "pheno_score", "matched_hpo",
                    "inheritance", "evidence", "epcr"])
        for i, r in enumerate(rows, 1):
            v1, v2 = r["v1"], r["v2"]
            vid = lambda v: f'{v["CHROM"]}:{v["POS"]}:{v["REF"]}>{v["ALT"]}' if v else ""
            ev = v1["evidence"] + (v2["evidence"] if v2 else [])
            w.writerow([i, r["type"], r["gene"], vid(v1), vid(v2), v1.get("HGVSc", ""), v1.get("HGVSp", ""),
                        v1["zyg"], r["vscore"], v1["pheno_score"], v1["matched_hpo"], r["inheritance"],
                        "; ".join(ev), r["epcr"]])

    # submission CSV
    with open(a.out_submission, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUBMISSION_FIELDS)
        w.writeheader()
        submitted = [r for r in rows if r["submit"]][: a.max_rows]
        for r in submitted:
            v1, v2 = r["v1"], r["v2"]
            primary = v1["pheno_score"] > 0
            note = "; ".join([r["gene"], r["inheritance"]] + v1["evidence"] + (v2["evidence"] if v2 else []))
            if not primary:
                note = "phenotype not matched; " + note
            w.writerow({
                "proband_id": a.proband_id,
                "chrom_1": v1["CHROM"], "pos_1": v1["POS"], "ref_1": v1["REF"], "alt_1": v1["ALT"],
                "chrom_2": v2["CHROM"] if v2 else "", "pos_2": v2["POS"] if v2 else "",
                "ref_2": v2["REF"] if v2 else "", "alt_2": v2["ALT"] if v2 else "",
                "epcr": r["epcr"],
                "finding_type": "primary" if primary else "secondary",
                "notes": note[:300],
            })

    print(f"Candidates scored : {len(variants)}")
    print(f"Ranked rows       : {len(rows)} (incl. {sum(r['type'] == 'compound_het' for r in rows)} compound-het pairs)")
    print(f"Submission rows   : {len(submitted)} -> {a.out_submission}")


if __name__ == "__main__":
    main()
