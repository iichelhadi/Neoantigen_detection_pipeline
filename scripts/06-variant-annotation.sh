#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 06_variant_annotation.sh - Annotate variants with VEP
# This script annotates VCF files using Ensembl's Variant Effect Predictor

set -e  # Exit immediately if a command exits with a non-zero status

# Parse command line arguments
while getopts ":c:s:" opt; do
  case ${opt} in
    c )
      config_file=$OPTARG
      ;;
    s )
      sample_info=$OPTARG
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

if [[ -z "$config_file" || -z "$sample_info" ]]; then
    echo "Usage: $0 -c CONFIG_FILE -s SAMPLE_INFO"
    exit 1
fi

# Source the config file to get global variables
source "$config_file"

# Create output directory for annotated variants
mkdir -p "${OUTPUT_DIR}/vcf_files_annotated"

echo "Starting variant annotation at $(date)"

# Function to annotate variants for a single patient
annotate_variants() {
    local patient_id=$1
    
    # Check if filtered VCF exists
    if [[ ! -f "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.filtered.vcf.gz" ]]; then
        echo "ERROR: Filtered VCF file not found for patient ${patient_id}"
        return 1
    fi
    
    # Check if annotated VCF already exists
    if [[ -f "${OUTPUT_DIR}/vcf_files_annotated/${patient_id}_somatic.filtered.VEP.vcf" ]]; then
        echo "Annotated VCF for patient ${patient_id} already exists, skipping"
        return 0
    fi
    
    echo "Annotating variants for patient ${patient_id}"
    
    # Run VEP on the filtered VCF
    vep --input_file "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.filtered.vcf.gz" \
        --output_file "${OUTPUT_DIR}/vcf_files_annotated/${patient_id}_somatic.filtered.VEP.vcf" -v \
        --everything --check_existing --total_length --allele_number --xref_refseq --format vcf --vcf --terms SO \
        --fasta "${VEP_FASTA}" \
        --offline --cache --dir_cache "${VEP_CACHE_DIR}" --force_overwrite \
        --plugin Frameshift --plugin Wildtype --dir_plugins "${VEP_PLUGINS_DIR}" \
        --synonyms "${VEP_SYNONYMS}"
    
    # Check if the annotation was successful
    if [[ ! -f "${OUTPUT_DIR}/vcf_files_annotated/${patient_id}_somatic.filtered.VEP.vcf" ]]; then
        echo "ERROR: VEP annotation failed for patient ${patient_id}"
        return 1
    fi
    
    echo "Annotation successful for patient ${patient_id}"
    return 0
}

# Extract patient IDs from sample info file (skip header line)
patient_ids=$(tail -n +2 "$sample_info" | cut -f1 | sort | uniq)

# Process each patient
for patient_id in $patient_ids; do
    annotate_variants "$patient_id"
done

echo "Variant annotation completed at $(date)"
