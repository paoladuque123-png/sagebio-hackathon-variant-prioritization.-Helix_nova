# VEP configuration

Annotation is run locally with Ensembl VEP 115.2 (offline cache, GRCh38) by
`scripts/03_vep_annotate.sh`. This document lists every annotation option, the
command-line flag that implements it, and the equivalent setting in the Ensembl VEP
web interface, so the configuration can be reviewed and reproduced in either form.

The proband VCF was not submitted to the VEP web interface; the web column is
provided only as a reference for readers who use the web tool with their own data.

## Input

| Setting | Web interface | Command line |
|---|---|---|
| Species / assembly | Homo_sapiens, GRCh38 | `--species homo_sapiens --assembly GRCh38` |
| Input | VCF file | `--input_file <ID>.norm.vcf.gz` (normalized with `bcftools norm -m -any`) |
| Transcript database | Ensembl/GENCODE transcripts | default Ensembl cache `homo_sapiens` (release 115); `--mane` adds the RefSeq NM_ ID of the MANE Select transcript |

## Identifiers

| Option | Web | Command line |
|---|---|---|
| Gene symbol | ✓ | `--symbol` |
| Transcript version | ✓ | `--transcript_version` |
| Protein | ✓ | `--protein` |
| UniProt | ✓ | `--uniprot` |
| HGVS | ✓ | `--hgvs --hgvsg` |
| MANE / canonical transcript flags | (shown by default) | `--mane --canonical` |

## Variants and frequency data

| Option | Web | Command line |
|---|---|---|
| Find co-located short variants | Yes | `--check_existing` |
| PubMed IDs for co-located variants | ✓ | `--pubmed` |
| 1000 Genomes global MAF | ✓ | `--af` |
| 1000 Genomes continental AF | ✓ | `--af_1kg` |
| gnomAD exomes AF | ✓ | `--af_gnomade` |
| gnomAD genomes AF | ✓ | `--af_gnomadg` |
| Maximum population AF | — | `--max_af` (used for filtering) |
| Variant synonyms, flagged variants, paralogue variants, Open Targets | ✗ | not used |
| ClinVar (significance, review status, condition) | via co-located variants | `--custom file=clinvar.vcf.gz,short_name=ClinVar,...,fields=CLNSIG%CLNREVSTAT%CLNDN` |

## Additional annotations

| Option | Web | Command line |
|---|---|---|
| DosageSensitivity | ✓ | `--plugin DosageSensitivity,file=...` |
| LOEUF | ✓ | `--plugin LOEUF,file=...,match_by=transcript` |
| Transcript biotype | ✓ | `--biotype` |
| Exon and intron numbers | ✓ | `--numbers` |
| Upstream/downstream distance | 5000 bp | `--distance 5000` |
| NMD | ✓ | `--plugin NMD` |
| UTRAnnotator | ✓ | `--plugin UTRAnnotator,file=...` |
| miRNA structure, RiboseqORFs, mutfunc, ProtVar, IntAct, MaveDB, Enformer | ✗ | not used |
| Regulatory region consequences | No | `--regulatory` not set |
| Phenotypes | ✓ | not run locally (requires a database connection to build the plugin file); gene–phenotype matching is done during prioritization |
| Gene Ontology | optional | not used |
| Geno2MP, Mastermind | ✗ | not used |

## Pathogenicity predictions

| Option | Web | Command line |
|---|---|---|
| SIFT | Prediction only | `--sift b` (prediction and score) |
| PolyPhen | Prediction only | `--polyphen b` (prediction and score) |
| AlphaMissense | ✓ | `--plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz` |
| CADD | Enabled (SNVs and InDels) | `--plugin CADD,snv=...,indels=...` |
| REVEL | ✓ | `--plugin REVEL,file=...` |
| ClinPred | ✓ | `--plugin ClinPred,file=...` |
| dbNSFP, AVI, EVE | ✗ | not used |

## Splicing and conservation

| Option | Web | Command line |
|---|---|---|
| SpliceAI | Enabled, masked scores (SNVs and indels) | `--plugin SpliceAI,snv=...,indel=...` (masked files) |
| MaxEntScan, BLOSUM62, ancestral allele | ✗ | not used |

## Filtering options

| Option | Web | Command line |
|---|---|---|
| Filter by frequency | No filtering | not set — frequency filtering is applied after annotation on allele-specific gnomAD `MAX_AF` |
| Coding regions only | ✗ | not set |
| Restrict results | Show all results | not set; `--flag_pick` marks one consequence per variant without removing the rest |

Rationale: filtering inside VEP uses the frequency of co-located variants, which at
multiallelic sites can remove a rare allele because a different allele at the same
position is common; restricting to the most severe consequence can drop the consequence
on the MANE Select transcript. Both steps are therefore deferred to the post-annotation
filter in `03_vep_annotate.sh`.

## Advanced options

| Option | Web | Command line |
|---|---|---|
| Buffer size | 5000 | default (`--buffer_size 5000`) |
| Right align variants | No | default (`--shift_3prime` not set; HGVS is 3'-shifted by default) |
| Parallel processes | — | `--fork $THREADS_VEP` |
| Output | — | `--vcf --compress_output bgzip` |
