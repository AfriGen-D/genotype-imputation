# H3ABionet ChipImputation Pipeline - Complete Usage Guide

## Table of Contents
1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Pipeline Overview](#pipeline-overview)
4. [Input Requirements](#input-requirements)
5. [Running the Pipeline](#running-the-pipeline)
6. [Output Structure](#output-structure)
7. [Advanced Configuration](#advanced-configuration)
8. [Troubleshooting](#troubleshooting)

## Quick Start

```bash
# Basic run with Docker
nextflow run main_simple.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile docker

# Run with Singularity on HPC
nextflow run main_simple.nf \
    --input samplesheet.csv \
    --outdir /scratch/results \
    -profile singularity,slurm
```

## Installation

### Prerequisites
- Nextflow (>=23.04.0)
- Docker/Singularity/Conda
- At least 8GB RAM for small datasets
- 100GB+ disk space for reference panels

### Setup
```bash
# Clone the repository
git clone https://github.com/h3abionet/chipimputation.git
cd chipimputation

# Test installation
nextflow run main_simple.nf --help
```

## Pipeline Overview

The pipeline performs the following steps:

### 1. Quality Control (Preprocessing)
- **CHECK_FILES**: Validates input VCF files
- **CHECK_CHROMOSOME**: Verifies chromosome naming consistency
- **QC_DUPL**: Removes duplicate variants
- **SPLIT_MULTI_ALLELIC**: Splits multi-allelic variants
- **FILTER_MIN_AC**: Filters by minimum allele count
- **TARGET_QC**: Comprehensive QC using BCFtools
- **QC_SITE_MISSINGNESS**: Filters by site missingness
- **SITES_ONLY**: Extracts sites for phasing

### 2. Phasing
- **EAGLE_PHASING**: Reference-based phasing using Eagle2

### 3. Imputation
- **IMPUTE_MINIMAC4**: Genotype imputation using Minimac4
- **EXTRACT_IMPUTE_INFO**: Extracts imputation quality metrics
- **COMBINE_IMPUTE**: Combines imputed chunks
- **COMBINE_INFO**: Merges info files

### 4. Reporting
- **FILTER_INFO_BY_TARGET**: Filters by R² threshold
- **REPORT_WELL_IMPUTED**: Reports well-imputed variants
- **PLOT_PERFORMANCE**: Generates performance plots
- **REPORT_ACCURACY**: Calculates accuracy metrics
- **PLOT_ACCURACY**: Creates accuracy visualizations
- **PLOT_R2_MAF**: R² vs MAF plots
- **PLOT_FREQ_COMPARISON**: Frequency distribution plots

## Input Requirements

### Sample Sheet Format
Create a CSV file with the following columns:

```csv
sample,vcf,population,study
sample1,/path/to/sample1.vcf.gz,AFR,study1
sample2,/path/to/sample2.vcf.gz,EUR,study1
sample3,/path/to/sample3.vcf.gz,EAS,study2
```

### VCF Requirements
- BGZipped and indexed VCF files (`.vcf.gz` and `.vcf.gz.tbi`)
- Chromosome naming consistent with reference panel
- Quality scores in QUAL field (optional)

### Reference Panels
Specify reference panels in the configuration:

```groovy
params.ref_panels = [
    ['1000G_phase3', '/path/to/1000g.m3vcf.gz', '/path/to/1000g.vcf.gz'],
    ['CAAPA', '/path/to/caapa.m3vcf.gz', '/path/to/caapa.vcf.gz']
]
```

## Running the Pipeline

### Basic Usage

```bash
nextflow run main_simple.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

### With Specific Chromosomes

```bash
nextflow run main_simple.nf \
    --input samplesheet.csv \
    --chromosomes 21,22 \
    --outdir results \
    -profile singularity
```

### Full Production Run

```bash
nextflow run main_simple.nf \
    --input samplesheet.csv \
    --ref_panels "[[\"1000G\",\"/ref/1000g.m3vcf.gz\",\"/ref/1000g.vcf.gz\"]]" \
    --eagle_genetic_map /ref/genetic_map_hg19.txt.gz \
    --reference_genome /ref/human_g1k_v37.fasta \
    --chromosomes ALL \
    --chunk_size 5000000 \
    --min_ac 2 \
    --site_missingness 0.05 \
    --r2_threshold 0.3 \
    --outdir /scratch/results \
    -profile singularity,slurm \
    -resume
```

### Test Run

```bash
# Quick test with stub mode (creates empty files)
nextflow run main_simple.nf \
    -profile test_complete \
    -stub-run

# Full test with sample data
./run_test.sh
```

## Output Structure

```
results/
├── qc/
│   ├── {sample}/
│   │   ├── *.nodup.vcf.gz          # Deduplicated VCF
│   │   ├── *.split.vcf.gz          # Multi-allelic split
│   │   ├── *.minac.vcf.gz          # AC filtered
│   │   ├── *.qc.vcf.gz             # QC passed variants
│   │   ├── *.stats                 # BCFtools statistics
│   │   ├── *.frq                   # Allele frequencies
│   │   └── *.missing               # Missingness report
├── phased/
│   ├── {sample}/
│   │   ├── *.phased.vcf.gz         # Phased genotypes
│   │   └── *.phasing.log           # Eagle log
├── imputed/
│   ├── {sample}/
│   │   ├── *.dose.vcf.gz           # Imputed dosages
│   │   ├── *.info                  # Imputation info
│   │   └── combined/
│   │       ├── *.combined.vcf.gz   # Combined chunks
│   │       └── *.combined.info.gz  # Combined info
├── reports/
│   ├── {sample}/
│   │   ├── *.well_imputed.txt      # Well imputed summary
│   │   ├── *.accuracy.txt          # Accuracy report
│   │   └── *.accuracy.tsv          # Accuracy data
├── plots/
│   ├── {sample}/
│   │   ├── *.performance.png       # Performance plot
│   │   ├── *.accuracy.png          # Accuracy plot
│   │   ├── *.r2_maf.png           # R² vs MAF
│   │   └── *.freq_comparison.png   # Frequency plots
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── pipeline_dag.html
```

## Advanced Configuration

### Custom Profiles

Create a custom profile in `conf/custom.config`:

```groovy
profiles {
    my_cluster {
        process {
            executor = 'slurm'
            queue = 'genomics'
            module = 'singularity/3.8.0'
        }
        singularity {
            enabled = true
            autoMounts = true
            cacheDir = '/shared/singularity_cache'
        }
    }
}
```

### Resource Allocation

Adjust resources in `conf/base.config`:

```groovy
process {
    withName: 'IMPUTE_MINIMAC4' {
        cpus = 8
        memory = 64.GB
        time = 48.h
    }
}
```

### Container Configuration

Use specific container versions:

```groovy
process {
    withName: 'EAGLE_PHASING' {
        container = 'mamana/eagle-vcf-processing:2.4.1'
    }
    withName: 'PLOT_.*' {
        container = 'mamana/python-plotting:1.0.0'
    }
}
```

## Parameters Reference

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to samplesheet CSV | Required |
| `--outdir` | Output directory | `./results` |
| `--chromosomes` | Chromosomes to process | `ALL` |
| `--chunk_size` | Chunk size for processing | `5000000` |
| `--min_ac` | Minimum allele count | `2` |
| `--site_missingness` | Max site missing rate | `0.05` |
| `--r2_threshold` | R² threshold for filtering | `0.3` |
| `--minRatio` | Minimac4 min ratio | `0.01` |
| `--eagle_pbwt_iters` | Eagle PBWT iterations | `2` |
| `--qc_plots` | Generate QC plots | `true` |
| `--max_cpus` | Maximum CPUs | `16` |
| `--max_memory` | Maximum memory | `128.GB` |
| `--max_time` | Maximum time | `240.h` |

## Troubleshooting

### Common Issues

#### 1. Memory Errors
```bash
# Increase memory allocation
nextflow run main_simple.nf \
    --max_memory 256.GB \
    --max_cpus 32
```

#### 2. Resume Failed Run
```bash
# Resume from last successful step
nextflow run main_simple.nf \
    -resume \
    --input samplesheet.csv
```

#### 3. Container Issues
```bash
# Pull containers manually
singularity pull docker://mamana/python-plotting:1.0.0
```

#### 4. Work Directory Space
```bash
# Use scratch directory
export NXF_WORK=/scratch/users/$USER/work
nextflow run main_simple.nf ...
```

### Debug Mode

```bash
# Enable debug output
nextflow run main_simple.nf \
    -with-trace \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-dag dag.png
```

### Getting Help

```bash
# Show help
nextflow run main_simple.nf --help

# Check version
nextflow -version

# View pipeline configuration
nextflow config -show-profiles
```

## Support

- GitHub Issues: https://github.com/h3abionet/chipimputation/issues
- Documentation: https://h3abionet.github.io/chipimputation
- Slack: https://h3africa.slack.com/channels/chipimputation

## Citation

If you use this pipeline, please cite:

> H3ABioNet Consortium. (2023). H3ABioNet Genotype Imputation Pipeline. 
> https://github.com/h3abionet/chipimputation

## License

This pipeline is released under the MIT License.