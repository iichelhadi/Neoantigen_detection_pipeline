#!/bin/bash
# 01_preprocess.sh - Quality control and preprocessing of raw FASTQ files
# This script performs adapter trimming and quality filtering using TrimGalore
# and generates QC reports using FastQC

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
mkdir -p "${OUTPUT_DIR}/trim_fq"
mkdir -p "${OUTPUT_DIR}/trim_QC"

# Process samples
echo "Starting preprocessing at $(date)"

# Function to process a single sample
process_sample() {
    local sample=$1
    local r1_file="${RAW_DATA_DIR}/${sample}/${sample}_1.fastq.gz"
    local r2_file="${RAW_DATA_DIR}/${sample}/${sample}_2.fastq.gz"
    
    # Check if files exist
    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        echo "ERROR: Could not find raw fastq files for sample ${sample}"
        echo "Expected: $r1_file and $r2_file"
        return 1
    fi

    # Check if output already exists
    if [[ -f "${OUTPUT_DIR}/trim_fq/${sample}_1_val_1.fq.gz" && \
          -f "${OUTPUT_DIR}/trim_fq/${sample}_2_val_2.fq.gz" ]]; then
        echo "Trimmed fastq files for ${sample} already exist, skipping trimming"
    else
        echo "Running TrimGalore on sample ${sample}"
        ${TRIMGALORE_PATH}/trim_galore --paired \
            "$r1_file" "$r2_file" \
            -o "${OUTPUT_DIR}/trim_fq/"
    fi

    # Run FastQC on trimmed files
    if [[ ! -f "${OUTPUT_DIR}/trim_QC/${sample}_1_val_1_fastqc.zip" || \
          ! -f "${OUTPUT_DIR}/trim_QC/${sample}_2_val_2_fastqc.zip" ]]; then
        echo "Running FastQC on sample ${sample}"
        fastqc -o "${OUTPUT_DIR}/trim_QC/" "${OUTPUT_DIR}/trim_fq/${sample}_1_val_1.fq.gz" "${OUTPUT_DIR}/trim_fq/${sample}_2_val_2.fq.gz"
    else
        echo "FastQC reports for ${sample} already exist, skipping FastQC"
    fi

    return 0
}

# Extract samples from sample info file (skip header line)
samples=$(tail -n +2 "$sample_info" | cut -f2,4 | tr '\t' '\n' | sort | uniq)

# Process each sample
for sample in $samples; do
    process_sample "$sample"
done

# Run MultiQC to summarize FastQC reports
if command -v multiqc &> /dev/null; then
    echo "Running MultiQC to summarize QC reports"
    multiqc -f "${OUTPUT_DIR}/trim_QC/"*zip -o "${OUTPUT_DIR}/trim_QC/"
else
    echo "MultiQC not found, skipping summary report generation"
fi

echo "Preprocessing completed at $(date)"
