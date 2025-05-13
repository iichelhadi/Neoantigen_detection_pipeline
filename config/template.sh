# Configuration file for the neoantigen discovery pipeline
# This file defines paths, parameters, and resources used throughout the pipeline

# Project information
PROJECT_NAME="NeoantigenPipeline"
OUTPUT_DIR="/path/to/output"
RAW_DATA_DIR="/path/to/raw_data"
RNA_OUTPUT_DIR="${OUTPUT_DIR}/rna"

# Computational resources
THREADS=24

# Reference genomes and indices
REFERENCE_GENOME_FASTA="/path/to/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REFERENCE_GENOME_BWA_INDEX="/path/to/reference/bwa_index/GRCh38_bwaidx"
HISAT2_INDEX="/path/to/reference/hisat_index/hisat_index"
GTF_FILE="/path/to/reference/Homo_sapiens.GRCh38.gtf"

# Resources for variant calling
PANEL_OF_NORMALS="/path/to/resources/somatic-hg38_1000g_pon.hg38.vcf.gz"
GERMLINE_RESOURCE="/path/to/resources/somatic-hg38_af-only-gnomad.hg38.vcf.gz"

# VEP resources
VEP_FASTA="/path/to/vep/cache/homo_sapiens/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
VEP_CACHE_DIR="/path/to/vep/cache"
VEP_PLUGINS_DIR="/path/to/vep/cache/Plugins"
VEP_SYNONYMS="/path/to/vep/cache/homo_sapiens/GRCh38/chr_synonyms.txt"

# pVACtools settings
IEDB_DIR="/path/to/IEDB"
BINDING_THRESHOLD=500

# Paths to tools
TRIMGALORE_PATH="/path/to/TrimGalore"

# Conda environments
CONDA_ENV_PREPROCESS="biotools"
CONDA_ENV_ALIGN="biotools"
CONDA_ENV_RNA="biotools"
CONDA_ENV_OPTITYPE="optitype"
CONDA_ENV_VARIANT="biotools"
CONDA_ENV_VEP="VEP"
CONDA_ENV_PVACTOOLS="pvactools"
