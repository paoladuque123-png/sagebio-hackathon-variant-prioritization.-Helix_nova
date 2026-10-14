# Annotation and variant prioritization

Annotation is performed locally with Ensembl VEP (offline cache) by
`scripts/03\_vep\_annotate.sh`. No proband data are uploaded to external services.
Prioritization starts from `<ID>.candidates.tsv`.

## Analysis record

|Item|Value|
|-|-|
|Date of annotation|2026-10-13|
|VEP / cache version|115.2 / 115 (GRCh38, Ensembl)|
|ClinVar release (`##fileDate`)|2026-10-11|
|Plugins used|NMD, AlphaMissense (REVEL, ClinPred, CADD, SpliceAI, LOEUF, DosageSensitivity and UTRAnnotator supported but not installed)|
|Records in VCF / after normalization|5,007,018 / 5,100,885 (93,867 multiallelic sites split; 44,998 indels left-aligned; 0 removed)|
|Candidates after pre-filter|1,028 (HIGH 98; MODERATE 570; LOW 360, splice-region consequences; ClinVar P/LP 5)|
|Phenotype input (HPO terms)|HP:0004322 Short stature; HP:0000121 Nephrocalcinosis; HP:0003202 Skeletal muscle atrophy; HP:0002859 Rhabdomyosarcoma; HP:0001508 Failure to thrive; HP:0001622 Premature birth; HP:0200067 Recurrent spontaneous abortion; HP:0001518 Small for gestational age|
|HPO release (`data-version`)|hp/releases/2026-09-01|
|Sex|<!-- as provided in dataset -->|
|Inheritance models considered|AD/XLD (heterozygous), AR/XLR (homozygous or compound heterozygous), hemizygous X in males|
|Prioritization output|<!-- n candidates scored / ranked rows / compound-het pairs, from the step 8 log -->|

## Pre-filter (automated, `03\_vep\_annotate.sh`)

* Genotype: non-reference, called
* Quality: DP ≥ `MIN\_DP` (10), GQ ≥ `MIN\_GQ` (20)
* Population frequency: gnomAD `MAX\_AF` < `MAX\_POP\_AF` (0.01) or absent
* Effect: HIGH or MODERATE impact, splice-related consequence, or ClinVar
pathogenic/likely pathogenic (non-conflicting)

## Prioritization (automated, `04\_prioritize.py`)

Each candidate is scored with three components in \[0, 1]:

**variant\_score** — strongest evidence among:

|Evidence|Score|
|-|-|
|ClinVar Pathogenic (non-conflicting)|1.00|
|ClinVar Likely pathogenic (non-conflicting)|0.90|
|Predicted LoF in a constrained gene (LOEUF < 0.6)|0.95|
|Predicted LoF, other genes|0.85 (−0.15 if NMD-escaping)|
|Missense supported by ≥ 3 / 2 / 1 / 0 predictors (AlphaMissense LP, REVEL ≥ 0.75, ClinPred ≥ 0.5)|0.80 / 0.70 / 0.60 / 0.35|
|SpliceAI max Δ ≥ 0.5 / ≥ 0.2|0.80 / 0.50|
|Other MODERATE impact|0.35|
|Baseline|0.20|

Absence from gnomAD adds 0.05; `MAX\_AF` > 0.001 multiplies the score by 0.9.

**pheno\_score** — fraction of the proband's HPO terms annotated to the gene in
`genes\_to\_phenotype.txt`, counting annotations to the term itself or to any of its
descendants in the HPO hierarchy.

**inheritance** — consistency of the genotype with the gene's modes of inheritance in HPO:

|Situation|Factor|
|-|-|
|Homozygous / hemizygous (male, chrX), or heterozygous in an AD/XLD gene|1.00|
|Two qualifying heterozygous variants in an AR/XL gene (scored as one compound-het pair)|1.00|
|Mode of inheritance not annotated|0.80|
|Single heterozygous variant in an AR gene, or member of a pair (pair row is submitted instead)|0.35|

Compound-heterozygous pairs combine the two best heterozygous variants (variant\_score ≥ 0.5)
in a recessive gene; the pair's variant\_score is the geometric mean of both. Phase cannot be
determined from a single proband without parental samples or read-backed phasing, so pairs
are reported as candidates with phase not determined.

## Rationale and limitations of the scoring weights

The scoring weights are **heuristic**: they were set by the authors and were not fitted or
calibrated on data. Their relative order follows the strength of the ACMG/AMP evidence
categories; the exact numeric values are design choices. Because automated scoring uses only the
ranking of the submitted rows, the order of the evidence categories and the phenotype floor
determine the result more than the exact values.

|Element|Basis|
|-|-|
|Evidence order: ClinVar P/LP > predicted LoF > missense supported by predictors > other|ACMG/AMP evidence strength (Richards et al., 2015, *Genet Med*): PS/PP5 previously reported, PVS1 null variant, PP3 computational evidence|
|Bonus for absence from gnomAD|ACMG/AMP criterion PM2|
|Reduction for NMD-escaping LoF|PVS1 decision tree (Abou Tayoun et al., 2018, *Hum Mutat*)|
|AlphaMissense `likely\_pathogenic` class|Class defined by the model authors (Cheng et al., 2023, *Science*)|
|REVEL ≥ 0.75; SpliceAI Δ ≥ 0.2 / ≥ 0.5|Thresholds from the original publications (Ioannidis et al., 2016, *Am J Hum Genet*; Jaganathan et al., 2019, *Cell*); implemented but not active in this run (plugins not installed)|
|Numeric values (1.00, 0.90, 0.85, 0.60, 0.35, 0.20), +0.05 / ×0.9 frequency adjustments, phenotype floor 0.3 and weight 0.7, inheritance factors 0.80 and 0.35, geometric mean for compound-heterozygous pairs|Author choice, not calibrated|

Limitations of this run:

* Only the NMD and AlphaMissense plugins were active. Missense variants therefore reach at
most 0.60 (AlphaMissense only), predicted LoF variants cannot be up-weighted by gene constraint
(no LOEUF), and non-canonical splice variants remain at the baseline score (no SpliceAI).
* Phase of heterozygous variants in the same gene cannot be determined (single proband, no
parental samples, no read-backed phasing); compound-heterozygous pairs are candidates.
* HaplotypeCaller calls were not hard-filtered or recalibrated (VQSR); the DP/GQ pre-filter and
manual review of top candidates mitigate false positives.
* Gene–phenotype matching depends on HPO annotation completeness; genes with sparse annotation
are down-weighted but retained through the phenotype floor.

## From ranking to EPCR

`EPCR = variant\_score × (0.3 + 0.7 × pheno\_score) × inheritance`, bounded to (0, 1].
Findings are ranked by EPCR and the top 10 are written to the submission file. EPCR is a
relative ranking score, not a calibrated probability; only the ranking affects automated
scoring.

`finding\_type` is `primary` when the gene matches at least one HPO term of the proband and
`secondary` otherwise (for example, a ClinVar pathogenic variant unrelated to the phenotype).
The `notes` column summarizes gene, inheritance and supporting evidence.

Top candidates are reviewed manually before submission (see "Manual evaluation").

## ClinVar P/LP variants in disease-associated genes (step 9)

Independently of the phenotype-driven ranking, `05\_clinvar\_plp.py` lists every variant of the
proband (alternate genotype, DP ≥ 10, GQ ≥ 20; no population-frequency filter) that meets both:

1. ClinVar classification exactly `Pathogenic`, `Likely\_pathogenic` or
`Pathogenic/Likely\_pathogenic` (conflicting interpretations and combined classifications such
as risk factor or low penetrance are excluded);
2. the gene has at least one disease (OMIM, Orphanet or DECIPHER) annotated in HPO
`genes\_to\_phenotype.txt`.

Each variant is reported with zygosity, gnomAD `MAX\_AF`, ClinVar review stars and condition, the
gene's diseases and mode(s) of inheritance, the proband's HPO terms matched by the gene, and an
interpretation (biallelic/hemizygous; heterozygous in a dominant gene; heterozygous in a recessive
gene, i.e. carrier unless a second allele is present). Rows are sorted by number of matched HPO terms
and review stars. This list is used to review the ranking and to report secondary findings; it does
not modify EPCR.

## Manual evaluation

The candidate list generated by this prioritization (step 8 ranking and submission file) and the
list of ClinVar pathogenic / likely pathogenic variants in disease-associated genes (step 9) were
evaluated manually to verify the congruence of the process: consistency of each finding with the
proband's phenotype, the gene's mode of inheritance and the genotype, and agreement between the
automated evidence and the annotation sources. Any change made to the automated ranking as a
result of this review is recorded in this document.

## Data handling

Annotation outputs (`\*.norm.vcf.gz`, `\*.vep.vcf.gz`, `\*.vep.tsv`, `\*.candidates.tsv`)
the ranked table (`\*.ranked.tsv`) and the ClinVar list (`clinvar\_plp/`) contain per-variant genotypes of the proband. They are kept only on the analysis server,
are excluded from this repository (`.gitignore`), and are deleted together with the
FASTQ, BAM and VCF files at the end of the Hackathon, as required by the Data Transfer
Agreement. The ranked candidate list (submission CSV) and this documentation are retained.

