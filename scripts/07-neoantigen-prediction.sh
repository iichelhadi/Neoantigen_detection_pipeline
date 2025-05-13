#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# 07_neoantigen_prediction.sh - Predict neoantigens using pVACtools
# This script uses pVACseq to predict neoantigens from VEP-annotated VCFs

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

# Create output directory for neoantigen predictions
mkdir -p "${OUTPUT_DIR}/pvacseq"

echo "Starting neoantigen prediction at $(date)"

# Function to predict neoantigens for a single patient
predict_neoantigens() {
    local patient_id=$1
    local hla_string=$2  # Format: HLA-A*02:01,HLA-A*11:01,etc.
    
    # Check if annotated VCF exists
    if [[ ! -f "${OUTPUT_DIR}/vcf_files_annotated/${patient_id}_somatic.filtered.VEP.vcf" ]]; then
        echo "ERROR: Annotated VCF file not found for patient ${patient_id}"
        return 1
    fi
    
    # Check if output directory already exists and has results
    if [[ -d "${OUTPUT_DIR}/pvacseq/${patient_id}" && -f "${OUTPUT_DIR}/pvacseq/${patient_id}/MHC_Class_I/${patient_id}_T.filtered.tsv" ]]; then
        echo "Neoantigen prediction results for patient ${patient_id} already exist, skipping"
        return 0
    fi
    
    echo "Predicting neoantigens for patient ${patient_id} with HLA types: ${hla_string}"
    
    # Run pVACseq
    pvacseq run \
        --iedb-install-directory "${IEDB_DIR}" \
        --binding-threshold "${BINDING_THRESHOLD:-500}" \
        --n-threads "${THREADS}" \
        --net-chop-method "cterm" \
        --netmhc-stab \
        --exclude-NAs \
        --normal-sample-name "${patient_id}_N" \
        "${OUTPUT_DIR}/vcf_files_annotated/${patient_id}_somatic.filtered.VEP.vcf" \
        "${patient_id}_T" \
        "${hla_string}" \
        MHCflurry MHCnuggetsI NetMHC NetMHCpan \
        "${OUTPUT_DIR}/pvacseq/${patient_id}"
    
    # Check if the prediction was successful
    if [[ ! -f "${OUTPUT_DIR}/pvacseq/${patient_id}/MHC_Class_I/${patient_id}_T.filtered.tsv" ]]; then
        echo "WARNING: pVACseq prediction might have failed for patient ${patient_id}"
    else
        echo "Neoantigen prediction successful for patient ${patient_id}"
    fi
    
    return 0
}

# Extract patient HLA types from sample info file (skip header line)
# Format: patient_id hla_a hla_b hla_c
patient_hlas=$(tail -n +2 "$sample_info" | cut -f1,6,7,8)

# Process each patient
while read -r patient_id hla_a hla_b hla_c; do
    # Combine HLA types
    hla_string="${hla_a},${hla_b},${hla_c}"
    
    predict_neoantigens "$patient_id" "$hla_string"
done <<< "$patient_hlas"

echo "Neoantigen prediction completed at $(date)"

# Optional: Run additional pVACtools analyses (pVACbind, pVACvector)
# Uncomment and adapt as needed for your specific requirements

# echo "Running additional pVACbind analysis for specific epitopes"
# pvacbind run -e1 8,9,10,11 \
#     --iedb-install-directory "${IEDB_DIR}" \
#     -b "${BINDING_THRESHOLD:-500}" -t "${THREADS}" \
#     --net-chop-method "cterm" --netmhc-stab --exclude-NAs \
#     "${EPITOPE_FASTA}" "${EPITOPE_OUTPUT_PREFIX}" "${HLA_TYPES}" \
#     MHCflurry MHCnuggetsI NetMHC NetMHCpan "${OUTPUT_DIR}/pvacseq_additional"

# echo "Running pVACvector analysis for vector design"
# pvacvector run \
#     "${EPITOPE_FASTA}" \
#     "${SAMPLE_NAME}" "${HLA_TYPES}" \
#     MHCflurry MHCnuggetsI NetMHC NetMHCpan \
#     --iedb-install-directory "${IEDB_DIR}" \
#     -e1 8,9,10,11 \
#     -v "${OUTPUT_DIR}/vcf_files_annotated/${SAMPLE_NAME}_somatic.filtered.VEP.vcf" \
#     "${OUTPUT_DIR}/pvacvector/${SAMPLE_NAME}"
