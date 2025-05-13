#!/bin/bash
# Copyright (c) 2022 Elhadi Iich
# Licensed under the MIT License - see the LICENSE file in the root directory for details

# submit_pipeline.sh - Submit pipeline jobs to a cluster scheduler
# This script creates and submits job files for running the pipeline on an HPC cluster

# Parse command line arguments
while getopts ":c:s:p:q:m:t:" opt; do
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
    q )
      queue=$OPTARG
      ;;
    m )
      memory=$OPTARG
      ;;
    t )
      time=$OPTARG
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
    echo "Usage: $0 -c CONFIG_FILE -s SAMPLE_INFO [-p STEPS] [-q QUEUE] [-m MEMORY] [-t TIME]"
    echo ""
    echo "Options:"
    echo "  -c CONFIG_FILE    Path to the configuration file"
    echo "  -s SAMPLE_INFO    Path to the sample information file"
    echo "  -p STEPS          Comma-separated list of steps to run (default: all)"
    echo "                    Available steps: preprocess,dna_align,rna_align,hla_typing,variant_call,annotate,predict"
    echo "  -q QUEUE          Queue to submit jobs to (default: parallel24)"
    echo "  -m MEMORY         Memory to request in GB (default: 160)"
    echo "  -t TIME           Wall time to request in hours (default: 48)"
    exit 1
fi

# Set defaults
queue=${queue:-parallel24}
memory=${memory:-160}
time=${time:-48}
steps=${steps:-"preprocess,dna_align,rna_align,hla_typing,variant_call,annotate,predict"}

# Source the config file to get global variables
source "$config_file"

# Get the base directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Create job submission directory
mkdir -p "${OUTPUT_DIR}/jobs"
JOB_DIR="${OUTPUT_DIR}/jobs"

# Function to create and submit a job file
submit_job() {
    local step=$1
    local script_path=$2
    local job_name="neoantigen_${step}"
    local job_file="${JOB_DIR}/${job_name}.pbs"
    local cores=24  # Default cores
    local depends=""
    
    echo "Creating job file for step: ${step}"
    
    # Adjust cores and memory for different steps
    case $step in
        preprocess)
            cores=12
            ;;
        dna_align)
            cores=24
            depends="${JOB_IDS[preprocess]}"
            ;;
        rna_align)
            cores=24
            depends="${JOB_IDS[preprocess]}"
            ;;
        hla_typing)
            cores=12
            depends="${JOB_IDS[preprocess]}"
            ;;
        variant_call)
            cores=24
            depends="${JOB_IDS[dna_align]}"
            ;;
        annotate)
            cores=12
            depends="${JOB_IDS[variant_call]}"
            ;;
        predict)
            cores=24
            depends="${JOB_IDS[annotate]},${JOB_IDS[hla_typing]}"
            ;;
    esac
    
    # Add dependency if provided
    depend_opt=""
    if [[ -n "$depends" ]]; then
        depend_opt="#PBS -W depend=afterok:${depends}"
    fi
    
    # Create the PBS job file
    cat > "$job_file" <<EOF
#!/bin/bash

#PBS -P ${PROJECT_NAME}
#PBS -j oe
#PBS -N ${job_name}

#PBS -q ${queue}
#PBS -l select=1:ncpus=${cores}:mpiprocs=${cores}:mem=${memory}GB

#PBS -l walltime=${time}:00:00
${depend_opt}

cd \$PBS_O_WORKDIR;   ## this line is needed, do not delete and change.
np=\$( cat \${PBS_NODEFILE} |wc -l );  ### get number of CPUs, do not change

##--- Load modules and environment ---
source /etc/profile.d/rec_modules.sh

# Load required modules (customize as needed)
module load fastqc TrimGalore
module load bwa
module load picard
module load samtools
module load miniconda

# Set up for specific steps
case "${step}" in
    preprocess)
        module load fastqc TrimGalore
        source activate ${CONDA_ENV_PREPROCESS}
        ;;
    dna_align)
        module load bwa
        module load samtools
        module load picard
        module load gatk
        source activate ${CONDA_ENV_ALIGN}
        ;;
    rna_align)
        module load hisat2
        module load samtools
        module load stringtie
        source activate ${CONDA_ENV_RNA}
        ;;
    hla_typing)
        source activate ${CONDA_ENV_OPTITYPE}
        ;;
    variant_call)
        module load gatk
        module load samtools
        source activate ${CONDA_ENV_VARIANT}
        ;;
    annotate)
        source activate ${CONDA_ENV_VEP}
        ;;
    predict)
        source activate ${CONDA_ENV_PVACTOOLS}
        ;;
esac

# Run the script
bash ${script_path} -c ${config_file} -s ${sample_info}

mpirun -f \${PBS_NODEFILE} ./a.out
EOF
    
    # Submit the job
    echo "Submitting job for step: ${step}"
    job_id=$(qsub "$job_file" | cut -d '.' -f 1)
    echo "Job ID: $job_id"
    
    return $job_id
}

# Initialize job ID array
declare -A JOB_IDS

# Submit each requested step
IFS=',' read -ra STEP_ARRAY <<< "$steps"
for step in "${STEP_ARRAY[@]}"; do
    case $step in
        preprocess)
            job_id=$(submit_job "preprocess" "${SCRIPT_DIR}/01_preprocess.sh")
            JOB_IDS[preprocess]=$job_id
            ;;
        dna_align)
            job_id=$(submit_job "dna_align" "${SCRIPT_DIR}/02_dna_alignment.sh")
            JOB_IDS[dna_align]=$job_id
            ;;
        rna_align)
            job_id=$(submit_job "rna_align" "${SCRIPT_DIR}/03_rna_alignment.sh")
            JOB_IDS[rna_align]=$job_id
            ;;
        hla_typing)
            job_id=$(submit_job "hla_typing" "${SCRIPT_DIR}/04_hla_typing.sh")
            JOB_IDS[hla_typing]=$job_id
            ;;
        variant_call)
            job_id=$(submit_job "variant_call" "${SCRIPT_DIR}/05_variant_calling.sh")
            JOB_IDS[variant_call]=$job_id
            ;;
        annotate)
            job_id=$(submit_job "annotate" "${SCRIPT_DIR}/06_variant_annotation.sh")
            JOB_IDS[annotate]=$job_id
            ;;
        predict)
            job_id=$(submit_job "predict" "${SCRIPT_DIR}/07_neoantigen_prediction.sh")
            JOB_IDS[predict]=$job_id
            ;;
        *)
            echo "WARNING: Unknown step '${step}', skipping"
            ;;
    esac
done

echo "All jobs submitted successfully."
echo "Job dependencies:"
for step in "${!JOB_IDS[@]}"; do
    echo "  $step: ${JOB_IDS[$step]}"
done
