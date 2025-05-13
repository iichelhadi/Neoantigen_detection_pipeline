```markdown
# Installation Guide

## Prerequisites
- Linux system
- Git
- Conda

## Installation Steps

1. **Clone the repository**
   ```bash
   git clone https://github.com/iichelhadi/Neoantigen_detection_pipeline.git
   cd Neoantigen_detection_pipeline
   ```

2. **Create conda environments**
   ```bash
   conda env create -f envs/alignment.yaml
   conda env create -f envs/variant_calling.yaml
   conda env create -f envs/vep.yaml
   conda env create -f envs/pvactools.yaml
   ```

3. **Set up configuration**
   ```bash
   cp config/config.sh.template config/config.sh
   # Edit config.sh with your paths and settings
   ```

4. **Install reference data**
   - Download reference genome (GRCh38 recommended)
   - Create BWA index for the reference genome
   - Download HISAT2 index for RNA-seq alignment
   - Set up VEP cache
   - Download panel of normals for variant calling

5. **Verify installation**
   ```bash
   bash scripts/utils/check_dependencies.sh -c config/config.sh
   ```
```

