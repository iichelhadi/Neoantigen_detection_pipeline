#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 05_variant_calling.sh - Somatic variant calling
# This script performs somatic variant calling using GATK Mutect2

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

# Create output directory for variant calls
mkdir -p "${OUTPUT_DIR}/vcf_files"

echo "Starting variant calling at $(date)"

# Function to call variants for a single sample
call_variants() {
    local tumor_sample=$2
    local normal_sample=$1
    local patient_id=$3
    
    # First process is individual variant calling for each sample
    for sample in "$normal_sample" "$tumor_sample"; do
        # Skip if recalibrated BAM already exists
        if [[ -f "${OUTPUT_DIR}/BAM_files/${sample}_recal_reads.bam" ]]; then
            continue
        fi
        
        # Check if sorted, deduplicated BAM exists
        if [[ ! -f "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" ]]; then
            echo "ERROR: Required input BAM file not found for sample ${sample}"
            return 1
        fi
        
        # If raw variants VCF doesn't exist, call variants
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_raw_variants.vcf" ]]; then
            echo "Calling variants for sample ${sample}"
            gatk HaplotypeCaller \
                -R "${REFERENCE_GENOME_FASTA}" \
                -I "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_raw_variants.vcf"
        fi
        
        # Extract SNPs
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_raw_snps.vcf" ]]; then
            echo "Extracting SNPs for sample ${sample}"
            gatk SelectVariants \
                -R "${REFERENCE_GENOME_FASTA}" \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_raw_variants.vcf" -select-type SNP \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_raw_snps.vcf"
        fi
        
        # Extract INDELs
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_raw_INDEL.vcf" ]]; then
            echo "Extracting INDELs for sample ${sample}"
            gatk SelectVariants \
                -R "${REFERENCE_GENOME_FASTA}" \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_raw_variants.vcf" -select-type INDEL \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_raw_INDEL.vcf"
        fi
        
        # Filter SNPs
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_filtered_snps.vcf" ]]; then
            echo "Filtering SNPs for sample ${sample}"
            gatk VariantFiltration \
                -R "${REFERENCE_GENOME_FASTA}" \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_raw_snps.vcf" \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_filtered_snps.vcf" \
                -filter-name "QD_filter" -filter "QD < 2.0" \
                -filter-name "FS_filter" -filter "FS > 60.0" \
                -filter-name "MQ_filter" -filter "MQ < 40.0" \
                -filter-name "SOR_filter" -filter "SOR > 3.0" \
                -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
                -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
        fi
        
        # Filter INDELs
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_filtered_indels.vcf" ]]; then
            echo "Filtering INDELs for sample ${sample}"
            gatk VariantFiltration \
                -R "${REFERENCE_GENOME_FASTA}" \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_raw_INDEL.vcf" \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_filtered_indels.vcf" \
                -filter-name "QD_filter" -filter "QD < 2.0" \
                -filter-name "FS_filter" -filter "FS > 200.0" \
                -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0"
        fi
        
        # Exclude filtered variants
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_snps.vcf" ]]; then
            echo "Excluding filtered SNPs for sample ${sample}"
            gatk SelectVariants \
                --exclude-filtered \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_filtered_snps.vcf" \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_snps.vcf"
        fi
        
        if [[ ! -f "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_indels.vcf" ]]; then
            echo "Excluding filtered INDELs for sample ${sample}"
            gatk SelectVariants \
                --exclude-filtered \
                -V "${OUTPUT_DIR}/vcf_files/${sample}_filtered_indels.vcf" \
                -O "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_indels.vcf"
        fi
        
        # Base Quality Score Recalibration (BQSR)
        mkdir -p "${OUTPUT_DIR}/recal"
        if [[ ! -f "${OUTPUT_DIR}/recal/${sample}_recal_data.table" ]]; then
            echo "Running BQSR for sample ${sample}"
            gatk BaseRecalibrator \
                -R "${REFERENCE_GENOME_FASTA}" \
                -I "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" \
                --known-sites "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_snps.vcf" \
                --known-sites "${OUTPUT_DIR}/vcf_files/${sample}_bqsr_indels.vcf" \
                -O "${OUTPUT_DIR}/recal/${sample}_recal_data.table"
        fi
        
        if [[ ! -f "${OUTPUT_DIR}/BAM_files/${sample}_recal_reads.bam" ]]; then
            echo "Applying BQSR for sample ${sample}"
            gatk ApplyBQSR \
                -R "${REFERENCE_GENOME_FASTA}" \
                -I "${OUTPUT_DIR}/BAM_files/${sample}_sorted_dedup_reads.bam" \
                -bqsr "${OUTPUT_DIR}/recal/${sample}_recal_data.table" \
                -O "${OUTPUT_DIR}/BAM_files/${sample}_recal_reads.bam"
        fi
    done
    
    # Prepare BAM files for Mutect2
    echo "Preparing BAM files for Mutect2 analysis"
    
    # Add/replace read groups for normal sample
    if [[ ! -f "${OUTPUT_DIR}/BAM_files/${normal_sample}_recal_reads_final.bam" ]]; then
        echo "Adding read groups to normal sample ${normal_sample}"
        gatk AddOrReplaceReadGroups \
            -I "${OUTPUT_DIR}/BAM_files/${normal_sample}_recal_reads.bam" \
            -O "${OUTPUT_DIR}/BAM_files/${normal_sample}_recal_reads_final.bam" \
            --RGID "${patient_id}_N" \
            --RGSM "${patient_id}_N" \
            --RGPL "illumina" \
            --RGLB "${normal_sample}" \
            --RGPU "${normal_sample}"
        
        samtools index "${OUTPUT_DIR}/BAM_files/${normal_sample}_recal_reads_final.bam"
    fi
    
    # Add/replace read groups for tumor sample
    if [[ ! -f "${OUTPUT_DIR}/BAM_files/${tumor_sample}_recal_reads_final.bam" ]]; then
        echo "Adding read groups to tumor sample ${tumor_sample}"
        gatk AddOrReplaceReadGroups \
            -I "${OUTPUT_DIR}/BAM_files/${tumor_sample}_recal_reads.bam" \
            -O "${OUTPUT_DIR}/BAM_files/${tumor_sample}_recal_reads_final.bam" \
            --RGID "${tumor_sample}" \
            --RGSM "${tumor_sample}" \
            --RGPL "illumina" \
            --RGLB "${tumor_sample}" \
            --RGPU "${tumor_sample}"
        
        samtools index "${OUTPUT_DIR}/BAM_files/${tumor_sample}_recal_reads_final.bam"
    fi
    
    # Call somatic variants with Mutect2
    if [[ ! -f "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.vcf.gz" ]]; then
        echo "Calling somatic variants for patient ${patient_id}"
        gatk Mutect2 \
            -R "${REFERENCE_GENOME_FASTA}" \
            -I "${OUTPUT_DIR}/BAM_files/${tumor_sample}_recal_reads_final.bam" \
            -I "${OUTPUT_DIR}/BAM_files/${normal_sample}_recal_reads_final.bam" \
            -normal "${patient_id}_N" \
            --panel-of-normals "${PANEL_OF_NORMALS}" \
            --germline-resource "${GERMLINE_RESOURCE}" \
            -O "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.vcf.gz"
    fi
    
    # Filter Mutect2 calls
    if [[ ! -f "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.filtered.vcf.gz" ]]; then
        echo "Filtering somatic variants for patient ${patient_id}"
        gatk FilterMutectCalls \
            -R "${REFERENCE_GENOME_FASTA}" \
            -V "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.vcf.gz" \
            -O "${OUTPUT_DIR}/vcf_files/${patient_id}_somatic.filtered.vcf.gz"
    fi
    
    return 0
}

# Extract sample information from sample info file (skip header line)
# Format: patient_id normal_sample normal_rg tumor_sample tumor_rg
sample_pairs=$(tail -n +2 "$sample_info" | cut -f1,2,3,4,5)

# Process each sample pair
while read -r patient_id normal_sample normal_rg tumor_sample tumor_rg; do
    echo "Processing patient: ${patient_id}"
    call_variants "$normal_sample" "$tumor_sample" "$patient_id"
done <<< "$sample_pairs"

echo "Variant calling completed at $(date)"
