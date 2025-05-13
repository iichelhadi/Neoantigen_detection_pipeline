```markdown
# Troubleshooting Guide

This document provides solutions for common issues encountered when running the neoantigen pipeline.

## Installation Issues

### Conda Environment Creation Fails

**Problem**: Error when creating conda environments.

**Solution**:
- Ensure you have a stable internet connection
- Update conda: `conda update -n base conda`
- Check for conflicting channels in your `.condarc` file
- Try installing dependencies one by one to identify problematic packages
- Use `conda clean --all` to clear package caches before retrying

### Missing Dependencies After Installation

**Problem**: The dependency check script reports missing tools despite installation.

**Solution**:
- Ensure you've activated the right conda environment: `conda activate environment_name`
- Check if the tool is in your PATH: `which tool_name`
- Verify paths in config.sh match your actual installation locations
- Some tools may require manual installation steps after conda installation
- Reinstall problematic tools manually outside of conda if necessary

## Data Preprocessing Issues

### TrimGalore or FastQC Fails

**Problem**: Error during trimming or quality control.

**Solution**:
- Verify TrimGalore/FastQC installation: `which trim_galore`, `which fastqc`
- Check input FASTQ file integrity: `gzip -t your_file.fastq.gz`
- Ensure you have sufficient disk space for output files
- Check if input file paths in your commands are correct
- Try running with verbose flags for more detailed error messages

### Out of Memory During Preprocessing

**Problem**: Process killed due to insufficient memory.

**Solution**:
- Reduce number of threads used for processing
- Process samples one at a time rather than in parallel
- Use a machine with more RAM or configure your cluster job to request more memory
- Split large FASTQ files into smaller chunks and process separately

## Alignment Issues

### BWA Alignment Errors

**Problem**: BWA fails with error messages.

**Solution**:
- Verify reference genome path and index files exist: `ls -la /path/to/reference.*`
- Ensure the index files match your reference genome version
- Check if the reference genome was properly indexed: `bwa index reference.fa`
- Examine error messages carefully for clues (e.g., specific read errors)

### Low Alignment Rate

**Problem**: Very few reads aligning to the reference.

**Solution**:
- Check that you're using the correct reference genome version
- Verify sample species matches reference genome species
- Run FastQC to inspect read quality and contamination
- Look for adapter contamination or poor quality bases
- Consider using more sensitive alignment parameters

### Samtools/Picard Errors

**Problem**: Errors during SAM/BAM processing.

**Solution**:
- Verify input BAM/SAM files exist and aren't corrupted
- Check if you have the latest samtools/picard versions
- Ensure you have write permissions to output directories
- Try with increased memory allocation for Picard: `java -Xmx8g -jar picard.jar`
- Check for unsorted input when tools expect sorted files

## Variant Calling Issues

### Mutect2 Fails

**Problem**: GATK Mutect2 exits with errors.

**Solution**:
- Check if input BAM files are properly sorted and indexed
- Verify that read groups in BAM files match those in your command
- Increase Java heap size: `--java-options "-Xmx16g"`
- Ensure reference genome is properly formatted and indexed
- Check for problematic regions that might cause errors (e.g., excessive depth)

### No Variants Detected

**Problem**: Empty or nearly empty VCF files after variant calling.

**Solution**:
- Verify tumor purity (low purity can affect variant detection)
- Check coverage depth in your BAM files: `samtools depth`
- Ensure tumor and normal samples are correctly specified
- Try less stringent filtering parameters
- Manually examine BAM files in IGV to verify variants are present

### VCF Processing Errors

**Problem**: Errors during VCF filtering or manipulation.

**Solution**:
- Check if input VCF is properly formatted: `bcftools check input.vcf`
- Ensure VCF is sorted correctly for the reference genome
- Try bcftools to normalize variants: `bcftools norm`
- If compressed, ensure index files exist: `tabix -p vcf input.vcf.gz`

## Annotation Issues

### VEP Fails

**Problem**: Error during VEP annotation.

**Solution**:
- Verify VEP cache is properly installed and accessible
- Check VEP version compatibility with cache version
- Ensure you have all required plugins installed
- Try offline mode if online database access is problematic
- Check VCF format compatibility with VEP
- Increase memory allocation for VEP

### Missing or Incomplete Annotations

**Problem**: Annotations are incomplete or missing key information.

**Solution**:
- Verify you're using the `--everything` flag with VEP for comprehensive annotation
- Check if the correct plugins are specified and installed
- Ensure reference genome version matches your annotation database version
- Try updating to the latest VEP cache
- Check if variants fall in regions covered by the annotation database

## Neoantigen Prediction Issues

### pVACseq Fails

**Problem**: pVACseq errors out during prediction.

**Solution**:
- Verify VEP annotation format is correct and includes all required fields
- Check HLA format (should be HLA-A*02:01, not A*02:01 or HLA-A02:01)
- Ensure IEDB tools are properly installed and accessible
- Check Python version compatibility with pVACtools
- Verify all necessary prediction algorithms are installed

### No Neoantigens Predicted

**Problem**: Pipeline completes but no neoantigens are identified.

**Solution**:
- Check if somatic variants were detected (prerequisite for neoantigens)
- Verify binding threshold isn't too stringent (try increasing from 500nM to 1000nM)
- Ensure mutations lead to protein changes (missense, frameshift, etc.)
- Check HLA types are correct and supported by the prediction algorithms
- Look for errors in intermediate files in the pVACseq output directory

## HLA Typing Issues

### OptiType Fails

**Problem**: OptiType fails to complete HLA typing.

**Solution**:
- Check if input FASTQ files exist and have sufficient coverage
- Verify razers3 is properly installed and working
- Ensure you have the OptiType reference data available
- Try with different coverage or enumeration parameters
- Check Python version compatibility (OptiType typically requires Python 2.7)

### Inconsistent HLA Types

**Problem**: Different HLA typing tools give different results.

**Solution**:
- This is normal - different tools have different sensitivities
- Use high-coverage data for more accurate typing
- Consider running multiple HLA typing tools and taking a consensus
- Focus on tools with highest accuracy for your data type (WGS, WES, RNA-seq)

## Configuration and System Issues

### Path and Environment Issues

**Problem**: Tools not found despite being installed.

**Solution**:
- Check PATH environment variable: `echo $PATH`
- Ensure conda environment is activated
- Verify config.sh has correct paths to all tools and resources
- Use absolute paths instead of relative paths in commands
- Check file permissions for executables: `chmod +x script.sh`

### Resource Limitations

**Problem**: Jobs fail due to memory or time limits.

**Solution**:
- For cluster jobs, increase requested resources: `-l mem=32gb,walltime=48:00:00`
- Monitor resource usage to identify bottlenecks: `top`, `htop`, or cluster monitoring
- Split large jobs into smaller subtasks
- Consider using SSD storage for temporary files to reduce I/O bottlenecks
- Reduce thread count for memory-intensive steps

## Getting Help

If you encounter issues not covered in this guide:

1. **Check Logs**: Examine log files in the `logs/` directory for specific error messages
2. **Documentation**: Refer to the original tool documentation for specific error codes
3. **Search Online**: Many bioinformatics errors have been encountered by others
4. **GitHub Issues**: Open an issue on our GitHub repository with:
   - Description of the problem
   - Relevant log snippets
   - Your config.sh (with sensitive information removed)
   - Software versions: `conda list` output
   - System information: OS version, CPU, RAM
5. **Contact Developers**: For urgent issues, contact the pipeline maintainers directly
```

