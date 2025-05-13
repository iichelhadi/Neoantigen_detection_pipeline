#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# run_pipeline.sh - Master script to run the complete neoantigen discovery pipeline
# This script orchestrates the execution of all pipeline steps

set -e  # Exit immediately if a command exits with a non-zero status

# Parse command line arguments
while getopts ":c:s:p:" opt; do
  case ${opt} in
    c )
      config_file=$OPTARG
      ;;
    s )
      sample_info=$OPTARG
      ;;
    p )
      steps=$OPTARG
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
    echo "Usage: $0 -c CONFIG_FILE -s SAMPLE_INFO [-p STEPS]"
    echo ""
    echo "Options:"
    echo "  -c CONFIG_FILE    Path to the configuration file"
    echo "  -s SAMPLE_INFO    Path to the sample information file"
    echo "  -p STEPS          Comma-separated list of steps to run (default: all)"
    echo "                    Available steps: preprocess,dna_align,rna_align,hla_typing,variant_call,annotate,predict"
    exit 1
fi

# Set default steps if not specified
if [[ -z "$steps" ]]; then
    steps="preprocess,dna_align,rna_align,hla_typing,variant_call,annotate,predict"
fi

# Source the config file to get global variables
source "$config_file"

# Get the base directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create log directory
mkdir -p "${OUTPUT_DIR}/logs"
LOG_DIR="${OUTPUT_DIR}/logs"

# Function to run a pipeline step
run_step() {
    local step=$1
    local script_path=$2
    local log_file="${LOG_DIR}/${step}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "=========================================================="
    echo "Starting step: ${step} at $(date)"
    echo "Log file: ${log_file}"
    echo "=========================================================="
    
    # Run the script and redirect output to log file
    bash "${script_path}" -c "${config_file}" -s "${sample_info}" 2>&1 | tee "${log_file}"
    
    # Check exit status
    local status=${PIPESTATUS[0]}
    if [[ $status -eq 0 ]]; then
        echo "Step ${step} completed successfully"
    else
        echo "ERROR: Step ${step} failed with status ${status}"
        exit $status
    fi
}

# Run each requested step
IFS=',' read -ra STEP_ARRAY <<< "$steps"
for step in "${STEP_ARRAY[@]}"; do
    case $step in
        preprocess)
            run_step "preprocess" "${SCRIPT_DIR}/01_preprocess.sh"
            ;;
        dna_align)
            run_step "dna_align" "${SCRIPT_DIR}/02_dna_alignment.sh"
            ;;
        rna_align)
            run_step "rna_align" "${SCRIPT_DIR}/03_rna_alignment.sh"
            ;;
        hla_typing)
            run_step "hla_typing" "${SCRIPT_DIR}/04_hla_typing.sh"
            ;;
        variant_call)
            run_step "variant_call" "${SCRIPT_DIR}/05_variant_calling.sh"
            ;;
        annotate)
            run_step "annotate" "${SCRIPT_DIR}/06_variant_annotation.sh"
            ;;
        predict)
            run_step "predict" "${SCRIPT_DIR}/07_neoantigen_prediction.sh"
            ;;
        *)
            echo "WARNING: Unknown step '${step}', skipping"
            ;;
    esac
done

echo "=========================================================="
echo "Pipeline completed successfully at $(date)"
echo "=========================================================="
