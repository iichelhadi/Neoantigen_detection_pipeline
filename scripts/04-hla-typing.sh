#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 04_hla_typing.sh - HLA typing using OptiType
# This script performs HLA typing from preprocessed fastq files

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

# Create output directory for HLA typing
mkdir -p "${OUTPUT_DIR}/optitype"

echo "Starting HLA typing at $(date)"

# Function to perform HLA typing for a single sample
perform_hla_typing() {
    local sample=$1
    
    # Check if trimmed fastq files exist
    local r1_trim="${OUTPUT_DIR}/trim_fq/${sample}_1_val_1.fq.gz"
    local r2_trim="${OUTPUT_DIR}/trim_fq/${sample}_2_val_2.fq.gz"
    
    if [[ ! -f "$r1_trim" || ! -f "$r2_trim" ]]; then
        echo "ERROR: Could not find trimmed fastq files for sample ${sample}"
        echo "Expected: $r1_trim and $r2_trim"
        return 1
    fi

    # Check if OptiType result already exists
    if [[ -f "${OUTPUT_DIR}/optitype/${sample}__result.tsv" ]]; then
        echo "OptiType result for ${sample} already exists, skipping"
        return 0
    fi

    echo "Running OptiType for sample ${sample}"
    OptiTypePipeline.py \
        --input "$r1_trim" "$r2_trim" \
        --dna --enumerate 3 \
        --outdir "${OUTPUT_DIR}/optitype" \
        --prefix "${sample}_" -v

    return 0
}

# Extract normal samples from sample info file (typically use normal samples for HLA typing)
# We filter for normal samples which are usually in the second column
normal_samples=$(tail -n +2 "$sample_info" | cut -f2 | sort | uniq)

# Process each normal sample
for sample in $normal_samples; do
    perform_hla_typing "$sample"
done

# Summarize HLA typing results
echo "Summarizing HLA typing results"
output_file="${OUTPUT_DIR}/optitype/hla_typing_summary.tsv"
echo -e "Sample\tA1\tA2\tB1\tB2\tC1\tC2" > "$output_file"

for sample in $normal_samples; do
    result_file="${OUTPUT_DIR}/optitype/${sample}__result.tsv"
    if [[ -f "$result_file" ]]; then
        # Extract HLA types from the OptiType result file
        # Assuming the format is consistent with OptiType's output
        hla_types=$(tail -n 1 "$result_file" | cut -f2-7)
        echo -e "${sample}\t${hla_types}" >> "$output_file"
    else
        echo "WARNING: No HLA typing result found for sample ${sample}"
    fi
done

echo "HLA typing completed at $(date)"
