# Neoantigen Discovery Pipeline

A comprehensive bioinformatics pipeline for detecting neoantigens from paired tumor-normal samples using whole exome sequencing (WES) and RNA-seq data.

## Overview

This pipeline processes raw NGS data to identify potential neoantigens in cancer samples through the following steps:

1. Quality control and preprocessing of raw fastq files
2. Alignment to reference genome 
3. HLA typing
4. Somatic variant calling
5. Variant annotation
6. Neoantigen prediction using multiple prediction algorithms

## Prerequisites

The pipeline requires the following software to be installed:

- TrimGalore (v0.6.6 or later)
- FastQC
- BWA
- HISAT2
- GATK
- Samtools
- Picard
- OptiType
- VEP (Variant Effect Predictor)
- pVACtools (pVACseq, pVACbind, pVACvector)
- Conda

## Directory Structure

```
neoantigen-pipeline/
├── config/
│   └── config.yaml
├── data/
│   └── sample_info.tsv
├── scripts/
│   ├── 01_preprocess.sh
│   ├── 02_dna_alignment.sh
│   ├── 03_rna_alignment.sh
│   ├── 04_hla_typing.sh
│   ├── 05_variant_calling.sh
│   ├── 06_variant_annotation.sh
│   ├── 07_neoantigen_prediction.sh
│   └── utils/
│       ├── job_submission.sh
│       └── check_dependencies.sh
├── envs/
│   ├── alignment.yaml
│   ├── variant_calling.yaml
│   ├── vep.yaml
│   └── pvactools.yaml
└── README.md
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/neoantigen-pipeline.git
cd neoantigen-pipeline
```

2. Set up conda environments:
```bash
# Create and activate environments from provided YAML files
conda env create -f envs/alignment.yaml
conda env create -f envs/variant_calling.yaml
conda env create -f envs/vep.yaml
conda env create -f envs/pvactools.yaml
```

3. Configure the pipeline:
Edit the `config/config.yaml` file to specify reference genomes, paths to tools, and other parameters.

## Usage

1. Prepare your sample information file:
```
# Example sample_info.tsv format
patient_id	normal_sample	normal_rg	tumor_sample	tumor_rg	hla_a	hla_b	hla_c
patient1	patient1_PBMC	PBMC_RG	patient1_T	T_RG	HLA-A*02:01,HLA-A*11:01	HLA-B*35:01,HLA-B*40:01	HLA-C*03:04,HLA-C*04:01
```

2. Run the complete pipeline:
```bash
./scripts/run_pipeline.sh -c config/config.yaml -s data/sample_info.tsv
```

Or run individual steps:
```bash
# For preprocessing and QC
./scripts/01_preprocess.sh -c config/config.yaml -s data/sample_info.tsv

# For DNA alignment
./scripts/02_dna_alignment.sh -c config/config.yaml -s data/sample_info.tsv

# For HLA typing
./scripts/04_hla_typing.sh -c config/config.yaml -s data/sample_info.tsv

# For variant calling
./scripts/05_variant_calling.sh -c config/config.yaml -s data/sample_info.tsv

# For neoantigen prediction
./scripts/07_neoantigen_prediction.sh -c config/config.yaml -s data/sample_info.tsv
```

## Output

The pipeline generates the following output directories:

- `trim_fq/`: Quality-trimmed fastq files
- `trim_QC/`: FastQC reports on trimmed data
- `BAM_files/`: Aligned BAM files
- `metrics/`: Quality metrics for alignments
- `hla_typing/`: HLA typing results from OptiType
- `vcf_files/`: Variant calls in VCF format (raw, filtered, and annotated)
- `neoantigen_predictions/`: Final neoantigen prediction results

## References

If you use this pipeline in your research, please cite the following tools:

- TrimGalore: [citation]
- BWA-MEM: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM
- GATK: McKenna A, et al. (2010) The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data
- OptiType: Szolek A, et al. (2014) OptiType: precision HLA typing from next-generation sequencing data
- VEP: McLaren W, et al. (2016) The Ensembl Variant Effect Predictor
- pVACtools: Hundal J, et al. (2020) pVACtools: a computational toolkit to identify and visualize cancer neoantigens

## Citation

If you use this pipeline in your research, please cite:

Elhadi Iich. (2022). Neoantigen Detection Pipeline [Computer software]. https://github.com/iichelhadi/Neoantigen_detection_pipeline

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
