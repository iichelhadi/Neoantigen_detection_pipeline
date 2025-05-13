#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 02_dna_alignment.sh - Alignment of DNA/WES data to reference genome
# This script performs BWA-MEM alignment, marks duplicates, and generates BAM files

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

# Create output directories
mkdir -p "${OUTPUT_DIR}/SAM_files"
mkdir -p "${OUTPUT_DIR}/BAM_files"
mkdir -p "${OUTPUT_DIR}/metrics"

echo "Starting DNA alignment at $(date)"

# Function to align a single sample
align_sample() {
    local sample=$1
    
    # Check if trimmed fastq files exist
    local r1_trim="${OUTPUT_DIR}/trim_fq/${sample}_1_val_1.fq.gz"
    local r2_trim="${OUTPUT_DIR}/trim_fq/${sample}_2_val_2.fq.gz"
    
    if [[ ! -f "$r1_trim" || ! -f "$r2_trim" ]]; then
        echo "ERROR: Could not find trimmed fastq files for sample ${sample}"
        echo "Expected: $r1_trim and $r2_trim"
        return 1
    fi

    # Check if final BAM already exists
    if [[ -f "${OUTPUT_DIR}/BAM_files/${sample}_recal_reads.bam" ]]; then
        echo "Processed BAM file for ${sample} already exists, skipping alignment"
        return 0
    fi

    # Alignment with BWA MEM
    if [[ ! -f "${OUTPUT_DIR}/SAM_files/${sample}.sam" && ! -f "${OUTPUT_DIR}/SAM_files/${sample}.sam.gz" ]]; then
        echo "Aligning sample ${sample} with BWA-MEM"
        bwa mem -M -Y -t "${THREADS}" -R "@RG\\tID:${sample}\\tLB:${sample}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample}" \
            "${REFERENCE_GENOME_BWA_INDEX}" "$r1_trim" "$r2_trim" > "${OUTPUT_DIR}/SAM_files/${sample}.sam"
    else
        echo "SAM file for ${sample} already exists, skipping alignment"
    fi

    # Mark duplicates and sort
    if [[ ! -f "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" ]]; then
        echo "Marking duplicates for sample ${sample}"
        gatk MarkDuplicatesSpark \
            -I "${OUTPUT_DIR}/SAM_files/${sample}.sam" \
            -M "${OUTPUT_DIR}/metrics/${sample}_dedup_metrics.txt" \
            -O "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam"
    fi

    # Collect alignment metrics
    if [[ ! -f "${OUTPUT_DIR}/metrics/${sample}_alignment_metrics.txt" ]]; then
        echo "Collecting alignment metrics for sample ${sample}"
        gatk CollectAlignmentSummaryMetrics \
            -R "${REFERENCE_GENOME_FASTA}" \
            -I "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" \
            -O "${OUTPUT_DIR}/metrics/${sample}_alignment_metrics.txt"
    fi

    # Collect insert size metrics
    if [[ ! -f "${OUTPUT_DIR}/metrics/${sample}_insert_metrics.txt" ]]; then
        echo "Collecting insert size metrics for sample ${sample}"
        gatk CollectInsertSizeMetrics \
            --INPUT "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" \
            --OUTPUT "${OUTPUT_DIR}/metrics/${sample}_insert_metrics.txt" \
            --HISTOGRAM_FILE "${OUTPUT_DIR}/metrics/${sample}_insert_size_histogram.pdf"
    fi

    # Calculate depth
    if [[ ! -f "${OUTPUT_DIR}/metrics/${sample}_depth_out.txt" ]]; then
        echo "Calculating depth for sample ${sample}"
        samtools depth -a "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" > "${OUTPUT_DIR}/metrics/${sample}_depth_out.txt"
    fi

    # Build BAM index
    echo "Building BAM index for sample ${sample}"
    gatk BuildBamIndex --INPUT "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam"

    return 0
}

# Extract samples from sample info file (skip header line)
samples=$(tail -n +2 "$sample_info" | cut -f2,4 | tr '\t' '\n' | sort | uniq)

# Process each sample
for sample in $samples; do
    align_sample "$sample"
done

echo "DNA alignment completed at $(date)"
