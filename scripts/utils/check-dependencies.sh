#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# check_dependencies.sh - Check if all required dependencies are installed
# This script verifies that the required tools and resources are available

set -e  # Exit immediately if a command exits with a non-zero status

# Parse command line arguments
while getopts ":c:" opt; do
  case ${opt} in
    c )
      config_file=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done

if [[ -z "$config_file" ]]; then
    echo "Usage: $0 -c CONFIG_FILE"
    exit 1
fi

# Source the config file to get global variables
source "$config_file"

# Function to check if a command is available
check_command() {
    local cmd=$1
    local msg=$2
    
    if command -v "$cmd" &> /dev/null; then
        echo "✓ $cmd is available"
        return 0
    else
        echo "✗ $cmd is not available: $msg"
        return 1
    fi
}

# Function to check if a file exists
check_file() {
    local file=$1
    local msg=$2
    
    if [[ -f "$file" ]]; then
        echo "✓ File exists: $file"
        return 0
    else
        echo "✗ File not found: $file - $msg"
        return 1
    fi
}

# Function to check if a directory exists
check_dir() {
    local dir=$1
    local msg=$2
    
    if [[ -d "$dir" ]]; then
        echo "✓ Directory exists: $dir"
        return 0
    else
        echo "✗ Directory not found: $dir - $msg"
        return 1
    fi
}

# Function to check if a conda environment exists
check_conda_env() {
    local env=$1
    
    if conda info --envs | grep -q "$env"; then
        echo "✓ Conda environment exists: $env"
        return 0
    else
        echo "✗ Conda environment not found: $env"
        return 1
    fi
}

echo "Checking dependencies for the neoantigen discovery pipeline..."
echo "============================================================="

# Check required directories
check_dir "$OUTPUT_DIR" "Output directory not found. Create it or update the config file."
check_dir "$RAW_DATA_DIR" "Raw data directory not found. Create it or update the config file."

# Check required reference files
check_file "$REFERENCE_GENOME_FASTA" "Reference genome FASTA file not found. Download it or update the config file."
check_file "${REFERENCE_GENOME_BWA_INDEX}.amb" "BWA index files not found. Create them with 'bwa index'."
check_file "$GTF_FILE" "GTF annotation file not found. Download it or update the config file."

# Check resources for variant calling
check_file "$PANEL_OF_NORMALS" "Panel of normals not found. Download it or update the config file."
check_file "$GERMLINE_RESOURCE" "Germline resource not found. Download it or update the config file."

# Check VEP resources
check_file "$VEP_FASTA" "VEP FASTA file not found. Download it or update the config file."
check_dir "$VEP_CACHE_DIR" "VEP cache directory not found. Download VEP cache or update the config file."
check_dir "$VEP_PLUGINS_DIR" "VEP plugins directory not found. Install VEP plugins or update the config file."
check_file "$VEP_SYNONYMS" "VEP synonyms file not found. Download it or update the config file."

# Check required tools
check_command "fastqc" "Install with 'conda install -c bioconda fastqc'"
check_command "bwa" "Install with 'conda install -c bioconda bwa'"
check_command "samtools" "Install with 'conda install -c bioconda samtools'"
check_command "picard" "Install with 'conda install -c bioconda picard'"
check_command "gatk" "Install with 'conda install -c bioconda gatk'"
check_command "hisat2" "Install with 'conda install -c bioconda hisat2'"
check_command "stringtie" "Install with 'conda install -c bioconda stringtie'"
check_command "vep" "Install with 'conda install -c bioconda ensembl-vep'"
check_command "pvacseq" "Install with 'conda install -c bioconda pvactools'"

# Check IEDB directory
check_dir "$IEDB_DIR" "IEDB directory not found. Download IEDB tools or update the config file."

# Check conda environments
check_conda_env "$CONDA_ENV_PREPROCESS"
check_conda_env "$CONDA_ENV_ALIGN"
check_conda_env "$CONDA_ENV_RNA"
check_conda_env "$CONDA_ENV_OPTITYPE"
check_conda_env "$CONDA_ENV_VARIANT"
check_conda_env "$CONDA_ENV_VEP"
check_conda_env "$CONDA_ENV_PVACTOOLS"

echo "============================================================="
echo "Dependency check completed."
