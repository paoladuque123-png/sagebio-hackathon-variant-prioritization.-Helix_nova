# =============================================================
# Pipeline configuration (example)
# Copy to config/config.sh and edit paths for your system:
#   cp config/config.example.sh config/config.sh
# config/config.sh is git-ignored.
# =============================================================

# --- Conda ---------------------------------------------------
CONDA_SH="/path/to/miniconda3/etc/profile.d/conda.sh"
ENV_QC="qc"           # conda env with fastqc
ENV_FASTP="qc"        # conda env with fastp
ENV_ALIGN="align"    # conda env with bwa + samtools (original system: "ancestry")
ENV_GATK="gatk_vcf"   # conda env with picard + gatk4 + samtools
ENV_VEP="vep_env"     # conda env with ensembl-vep + bcftools + python

# --- Sample list (one sample ID per line) --------------------
SAMPLES="samples/samples.txt"

# --- Input / output directories ------------------------------
FASTQ_DIR="/path/to/0_fastq"        # raw paired FASTQ: <ID>*_1.fastq.gz / <ID>*_2.fastq.gz
FASTQC_RAW_DIR="/path/to/1_fastqc/raw"  # FastQC/MultiQC reports on raw reads
FASTQC_TRIM_DIR="/path/to/1_fastqc/trimmed"  # FastQC reports on trimmed reads
QC_DIR="/path/to/1_fastqc"          # fastp reports
TRIMMED_DIR="/path/to/2_trimmed"    # trimmed FASTQ
BAM_DIR="/path/to/3_bam"
BAM_RG_DIR="/path/to/3_bam_name"
RECAL_DIR="/path/to/4_recalibration"
VCF_DIR="/path/to/5_vcf"
ANNOT_DIR="/path/to/7_annotation"
LOG_DIR="/path/to/logs"
TMP_DIR="/path/to/tmp"

# --- Reference (GRCh38, no-alt analysis set) -----------------
REF="/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# --- Known sites for BQSR ------------------------------------
KNOWN_SITES_1="/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
KNOWN_SITES_2="/path/to/Homo_sapiens_assembly38.dbsnp138.vcf"

# --- Annotation (VEP, local/offline) ------------------------
VEP_CACHE_DIR="/path/to/.vep"
VEP_CACHE_VERSION=115                 # must match the installed ensembl-vep version
VEP_TRANSCRIPT_SET="ensembl"          # ensembl (Ensembl/GENCODE; MANE gives RefSeq NM_ IDs) | merged | refseq
CLINVAR_VCF="/path/to/clinvar/clinvar.vcf.gz"
VEP_PLUGINS_DIR="/path/to/vep_plugins"  # VEP plugin .pm files (vep_install -a p)
VEP_PLUGIN_DATA="/path/to/vep_plugin_data"
# Each plugin below is enabled only if its file(s) exist; leave "" to skip
USE_NMD=1                             # NMD plugin (no data file needed)
ALPHAMISSENSE_FILE=""                 # AlphaMissense_hg38.tsv.gz (tabix-indexed)
REVEL_FILE=""                         # new_tabbed_revel_grch38.tsv.gz (see REVEL.pm header)
CLINPRED_FILE=""                      # ClinPred tabix file (see ClinPred.pm header)
CADD_SNV=""                           # CADD whole_genome_SNVs.tsv.gz (GRCh38)
CADD_INDEL=""                         # CADD gnomad.genomes indel file (GRCh38)
SPLICEAI_SNV=""                       # SpliceAI masked SNV scores (GRCh38)
SPLICEAI_INDEL=""                     # SpliceAI masked indel scores (GRCh38)
LOEUF_FILE=""                         # gnomAD constraint metrics (see LOEUF.pm header)
DOSAGE_FILE=""                        # Collins et al. pHaplo/pTriplo (see DosageSensitivity.pm header)
UTRANNOTATOR_FILE=""                  # uORF_5UTR_GRCh38_PUBLIC.txt (UTRAnnotator)

# --- Prioritization (HPO) -----------------------------------
HPO_OBO="/path/to/hpo/hp.obo"
HPO_G2P="/path/to/hpo/genes_to_phenotype.txt"
PROBAND_SEX="unknown"                 # male | female | unknown (as provided in the dataset)

# --- Submission ----------------------------------------------
SUBMISSION_USER="hf-username"         # Hugging Face username or team name
APPROACH_NAME="gatk-vep-hpo"          # short name of the approach (used in the CSV filename)
MAX_SUBMISSION_ROWS=10

# --- Candidate pre-filter ------------------------------------
MIN_DP=10
MIN_GQ=20
MAX_POP_AF=0.01                       # max gnomAD AF across populations (VEP MAX_AF)

# --- Read group ----------------------------------------------
RG_LB="GP"
RG_PL="MGI"
RG_PU="unit1"

# --- Resources -----------------------------------------------
THREADS_FASTQC=8
THREADS_FASTP=40
THREADS_BWA=35
THREADS_VEP=16
