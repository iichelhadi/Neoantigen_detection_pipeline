```markdown
# Usage Guide

## Input Requirements
- Paired tumor-normal WES data (FASTQ format)
- Optional: RNA-seq data (FASTQ format)
- Sample information file (tab-delimited)

## Sample Information File Format
Create a file like `data/sample_info.tsv` with the following columns:
```
patient_id   normal_sample   normal_rg     tumor_sample   tumor_rg      hla_a                   hla_b                   hla_c
patient1     patient1_PBMC   PBMC_RG       patient1_T     T_RG          HLA-A*02:01,HLA-A*11:01 HLA-B*35:01,HLA-B*40:01 HLA-C*03:04,HLA-C*04:01
```

## Running the Pipeline

### Complete Pipeline
```bash
bash scripts/run_pipeline.sh -c config/config.sh -s data/sample_info.tsv
```

### Individual Steps
```bash
# Preprocessing
bash scripts/01_preprocess.sh -c config/config.sh -s data/sample_info.tsv

# DNA alignment
bash scripts/02_dna_alignment.sh -c config/config.sh -s data/sample_info.tsv

# RNA alignment
bash scripts/03_rna_alignment.sh -c config/config.sh -s data/sample_info.tsv

# HLA typing
bash scripts/04_hla_typing.sh -c config/config.sh -s data/sample_info.tsv

# Variant calling
bash scripts/05_variant_calling.sh -c config/config.sh -s data/sample_info.tsv

# Variant annotation
bash scripts/06_variant_annotation.sh -c config/config.sh -s data/sample_info.tsv

# Neoantigen prediction
bash scripts/07_neoantigen_prediction.sh -c config/config.sh -s data/sample_info.tsv
```

### Cluster Submission
```bash
bash scripts/submit_pipeline.sh -c config/config.sh -s data/sample_info.tsv
```

## Output Files
- `trim_fq/`: Trimmed FASTQ files
- `BAM_files/`: Aligned BAM files
- `optitype/`: HLA typing results
- `vcf_files/`: Variant calls
- `vcf_files_annotated/`: Annotated variants
- `pvacseq/`: Neoantigen predictions
```

