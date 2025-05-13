```markdown
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

## Documentation

- [Installation Guide](docs/installation.md) - How to set up the pipeline and its dependencies
- [Usage Guide](docs/usage.md) - How to run the pipeline with your data
- [Troubleshooting Guide](docs/troubleshooting.md) - Solutions for common issues

## Quick Start

```bash
# Clone the repository
git clone https://github.com/iichelhadi/Neoantigen_detection_pipeline.git
cd Neoantigen_detection_pipeline

# Set up conda environments
conda env create -f envs/alignment.yaml
conda env create -f envs/variant_calling.yaml
conda env create -f envs/vep.yaml
conda env create -f envs/pvactools.yaml

# Configure the pipeline
cp config/config.sh.template config/config.sh
# Edit config.sh with your specific paths and settings

# Run the pipeline
bash scripts/run_pipeline.sh -c config/config.sh -s data/sample_info.tsv
```

## Prerequisites

The pipeline requires the following software (installed via conda environments):

- TrimGalore and FastQC (quality control)
- BWA and HISAT2 (alignment)
- GATK and Samtools (variant calling)
- OptiType (HLA typing)
- VEP (Variant Effect Predictor)
- pVACtools (neoantigen prediction)

See the [Installation Guide](docs/installation.md) for detailed requirements and setup instructions.

## Directory Structure

```
neoantigen-pipeline/
├── config/
│   └── config.sh.template
├── data/
│   └── sample_info.tsv.template
├── scripts/
│   ├── 01_preprocess.sh
│   ├── 02_dna_alignment.sh
│   ├── 03_rna_alignment.sh
│   ├── 04_hla_typing.sh
│   ├── 05_variant_calling.sh
│   ├── 06_variant_annotation.sh
│   ├── 07_neoantigen_prediction.sh
│   ├── run_pipeline.sh
│   ├── submit_pipeline.sh
│   └── utils/
│       └── check_dependencies.sh
├── envs/
│   ├── alignment.yaml
│   ├── variant_calling.yaml
│   ├── vep.yaml
│   └── pvactools.yaml
├── docs/
│   ├── installation.md
│   ├── usage.md
│   └── troubleshooting.md
└── README.md
```

## Output

The pipeline generates the following output directories:

- `trim_fq/`: Quality-trimmed fastq files
- `BAM_files/`: Aligned BAM files
- `optitype/`: HLA typing results
- `vcf_files/`: Variant calls in VCF format
- `vcf_files_annotated/`: Annotated variants
- `pvacseq/`: Neoantigen predictions

<<<<<<< HEAD
For detailed information about running the pipeline and interpreting results, see the [Usage Guide](docs/usage.md).
=======
## References

If you use this pipeline in your research, please cite the following tools:

- [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
- BWA-MEM: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM
- GATK: McKenna A, et al. (2010) The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data
- OptiType: Szolek A, et al. (2014) OptiType: precision HLA typing from next-generation sequencing data
- VEP: McLaren W, et al. (2016) The Ensembl Variant Effect Predictor
- pVACtools: Hundal J, et al. (2020) pVACtools: a computational toolkit to identify and visualize cancer neoantigens
>>>>>>> 3e9b13110389cfed40cd07459da134d3b8cb1c67

## Citation

If you use this pipeline in your research, please cite:

Elhadi Iich. (2022). Neoantigen Detection Pipeline [Computer software]. https://github.com/iichelhadi/Neoantigen_detection_pipeline

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
```
