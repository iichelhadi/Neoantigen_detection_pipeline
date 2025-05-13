#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 03_rna_alignment.sh - Alignment of RNA-seq data to reference genome
# This script performs RNA-seq alignment with HISAT2 and generates gene counts

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
mkdir -p "${RNA_OUTPUT_DIR}/SAM_files"
mkdir -p "${RNA_OUTPUT_DIR}/BAM_files"
mkdir -p "${RNA_OUTPUT_DIR}/feature_counts"
mkdir -p "${RNA_OUTPUT_DIR}/flagstat"
mkdir -p "${RNA_OUTPUT_DIR}/stringtie_assembly"
mkdir -p "${RNA_OUTPUT_DIR}/ballgown"

echo "Starting RNA-seq alignment at $(date)"

# Function to align a single RNA sample
align_rna_sample() {
    local sample=$1
    
    # Check if trimmed fastq files exist
    local r1_trim="${OUTPUT_DIR}/trim_fq/${sample}_1_val_1.fq.gz"
    local r2_trim="${OUTPUT_DIR}/trim_fq/${sample}_2_val_2.fq.gz"
    
    if [[ ! -f "$r1_trim" || ! -f "$r2_trim" ]]; then
        echo "ERROR: Could not find trimmed fastq files for sample ${sample}"
        echo "Expected: $r1_trim and $r2_trim"
        return 1
    fi

    # Align with HISAT2
    if [[ ! -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam" && ! -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam.gz" ]]; then
        echo "Aligning RNA sample ${sample} with HISAT2"
        hisat2 -p "${THREADS}" --dta -x "${HISAT2_INDEX}" \
            -1 "$r1_trim" -2 "$r2_trim" \
            -S "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam"
    else
        echo "SAM file for RNA sample ${sample} already exists, skipping alignment"
    fi

    # Run featureCounts
    if [[ ! -f "${RNA_OUTPUT_DIR}/feature_counts/${sample}_counts.txt" ]]; then
        echo "Running featureCounts for sample ${sample}"
        featureCounts -T "${THREADS}" -p -s 2 -t 'exon' -g 'gene_id' \
            -a "${GTF_FILE}" \
            -o "${RNA_OUTPUT_DIR}/feature_counts/${sample}_counts.txt" \
            "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam"
    fi

    # Convert SAM to BAM and sort
    if [[ ! -f "${RNA_OUTPUT_DIR}/BAM_files/${sample}_sorted.bam" ]]; then
        echo "Converting SAM to BAM and sorting for sample ${sample}"
        
        # Check if SAM is gzipped
        if [[ -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam.gz" && ! -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam" ]]; then
            gunzip "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam.gz"
        fi
        
        samtools view -Shb "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam" > "${RNA_OUTPUT_DIR}/BAM_files/${sample}.bam"
        samtools sort -@ "${THREADS}" -o "${RNA_OUTPUT_DIR}/BAM_files/${sample}_sorted.bam" "${RNA_OUTPUT_DIR}/BAM_files/${sample}.bam"
        rm "${RNA_OUTPUT_DIR}/BAM_files/${sample}.bam"  # Remove unsorted BAM
    fi

    # Generate flagstat report
    if [[ ! -f "${RNA_OUTPUT_DIR}/flagstat/flagstat_${sample}.txt" ]]; then
        echo "Generating flagstat report for sample ${sample}"
        samtools flagstat "${RNA_OUTPUT_DIR}/BAM_files/${sample}_sorted.bam" -@ "${THREADS}" > "${RNA_OUTPUT_DIR}/flagstat/flagstat_${sample}.txt"
    fi

    # Compress SAM file to save space
    if [[ -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam" && ! -f "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam.gz" ]]; then
        echo "Compressing SAM file for sample ${sample}"
        gzip "${RNA_OUTPUT_DIR}/SAM_files/${sample}.sam"
    fi

    # Generate transcript assembly with StringTie
    if [[ ! -f "${RNA_OUTPUT_DIR}/stringtie_assembly/${sample}_assembly.gtf" ]]; then
        echo "Generating transcript assembly with StringTie for sample ${sample}"
        stringtie "${RNA_OUTPUT_DIR}/BAM_files/${sample}_sorted.bam" \
            -l "${sample}" -p "${THREADS}" \
            -G "${GTF_FILE}" \
            -o "${RNA_OUTPUT_DIR}/stringtie_assembly/${sample}_assembly.gtf"
    fi

    return 0
}

# Extract RNA samples from sample info file (if different from DNA samples)
if [[ -z "$RNA_SAMPLES" ]]; then
    # If RNA_SAMPLES not explicitly defined, use the same samples as DNA
    RNA_SAMPLES=$(tail -n +2 "$sample_info" | cut -f2,4 | tr '\t' '\n' | sort | uniq)
fi

# Process each RNA sample
for sample in $RNA_SAMPLES; do
    align_rna_sample "$sample"
done

# Create merge list for StringTie
echo "Creating StringTie merge list"
merge_list="${RNA_OUTPUT_DIR}/stringtie_assembly/mergelist.txt"
> "$merge_list"  # Clear the file if it exists

for sample in $RNA_SAMPLES; do
    echo "${RNA_OUTPUT_DIR}/stringtie_assembly/${sample}_assembly.gtf" >> "$merge_list"
done

# Merge transcript assemblies
if [[ ! -f "${RNA_OUTPUT_DIR}/stringtie_assembly/stringtie_merged.gtf" ]]; then
    echo "Merging transcript assemblies with StringTie"
    stringtie --merge -p "${THREADS}" \
        -G "${GTF_FILE}" \
        -o "${RNA_OUTPUT_DIR}/stringtie_assembly/stringtie_merged.gtf" \
        "$merge_list"
fi

# Generate Ballgown input data
for sample in $RNA_SAMPLES; do
    ballgown_dir="${RNA_OUTPUT_DIR}/ballgown/${sample}"
    mkdir -p "$ballgown_dir"
    
    if [[ ! -f "${ballgown_dir}/${sample}.gtf" ]]; then
        echo "Generating Ballgown data for sample ${sample}"
        stringtie -e -B -p "${THREADS}" \
            -G "${RNA_OUTPUT_DIR}/stringtie_assembly/stringtie_merged.gtf" \
            -o "${ballgown_dir}/${sample}.gtf" \
            "${RNA_OUTPUT_DIR}/BAM_files/${sample}_sorted.bam"
    fi
done

echo "RNA-seq alignment completed at $(date)"
