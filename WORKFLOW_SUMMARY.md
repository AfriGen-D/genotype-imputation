# Workflow Summary

## Recent Updates (August 2025)

### Major Fixes and Improvements

#### 1. Output Format Update
- **Issue**: Minimac4 was producing SAV (Zstandard compressed) format files when using MSAV reference panels
- **Solution**: Changed output format to BCF (Binary VCF) for better compatibility
- **Modified**: `modules/impute.nf` - Updated `impute_minimac4` process to use `--output-format bcf`

#### 2. Reference Panel Path Resolution
- **Issue**: Reference panel paths contained `%s` placeholders that weren't being replaced with actual chromosome names
- **Solution**: Modified workflow to properly substitute chromosome names in reference panel paths
- **Modified**: `main.nf` - Updated `generate_frequency` workflow section to handle chromosome-specific paths

#### 3. R Analysis Container
- **Added**: New container `r-analysis:1.0` with R 4.3.2 and required packages
- **Packages**: dplyr, data.table, ggplot2, tidyr, readr, stringr, purrr, tibble, and visualization libraries
- **Purpose**: Support for pipeline's analysis and visualization processes that require R

#### 4. Configuration Updates
- **Added**: Test configuration for v6 phased data (`v6_chr21_phased.config`)
- **Updated**: Container specifications in `nextflow.config` for R analysis processes
- **Modified**: minRatio parameter adjusted to 0.01 for better variant ratio handling

### Container Registry

All containers are automatically built and published to:
- GitHub Container Registry: `ghcr.io/mamana/*`
- Docker Hub: `mamana/*`

### Available Containers

| Container | Version | Purpose |
|-----------|---------|---------|
| minimac4 | minimac4-4.1.6 | Imputation with Minimac4 |
| eagle-phasing | eagle-2.4.1 | Phasing with Eagle2 |
| vcf-processing | bcftools-1.20 | VCF manipulation with BCFtools |
| r-analysis | 1.0 | R-based analysis and visualization |
| imputation | minimac4-4.1.6 | Combined imputation tools |
| phasing | eagle-2.4.1 | Combined phasing tools |

## Pipeline Output Structure

```
output/
├── imputed/
│   └── <ref_panel>/
│       └── <target_dataset>/
│           ├── *.imputed        # BCF format imputed genotypes
│           └── *.imputed.info   # Imputation quality metrics
├── phased/
│   └── <target_dataset>/
│       └── *.phased.vcf.gz     # Phased VCF files
├── frqs/
│   └── <ref_panel>/
│       ├── *.frq               # Allele frequencies
│       └── *.vcf.gz            # Frequency annotated VCFs
├── plots/
│   └── <ref_panel>/
│       ├── *_r2_SNPcount.png   # R² vs SNP count plots
│       ├── *_MAF_r2.png        # MAF vs R² plots
│       └── *_performance.png   # Performance metrics
└── reports/
    └── <ref_panel>/
        ├── *_well_imputed.tsv  # Well-imputed variants report
        └── *_accuracy.tsv      # Imputation accuracy metrics
```

## Running the Pipeline

### Basic Usage

```bash
# With test data
nextflow run main.nf -profile test,singularity

# With custom config
nextflow run main.nf -c your_config.config -profile slurm,singularity
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `chromosomes` | 'ALL' | Chromosomes to impute (ALL, specific list, or range) |
| `minRatio` | 0.01 | Minimum ratio of target to reference variants |
| `chunk_size` | 10000000 | Chunk size for imputation (in base pairs) |
| `buffer_size` | 500000 | Buffer size for chunk overlap |
| `impute_info_cutoff` | 0.8 | R² cutoff for well-imputed variants |
| `outDir` | 'output' | Output directory path |

### Profiles

- `singularity`: Use Singularity containers
- `docker`: Use Docker containers
- `slurm`: Submit jobs to SLURM cluster
- `test`: Run with test dataset

## Troubleshooting

### Common Issues

1. **Missing R packages**: Ensure the r-analysis container is available
2. **Chromosome naming**: Check that reference panel uses consistent naming (chr1 vs 1)
3. **Memory issues**: Increase memory allocation in config for large datasets
4. **Output format**: BCF files can be converted to VCF using `bcftools view`

### Monitoring Progress

Check the pipeline status with:
```bash
# View running jobs
squeue -u $USER

# Check logs
tail -f .nextflow.log

# View work directory for specific process
ls -la work/<hash>/
```

## Support

For issues or questions:
- GitHub Issues: https://github.com/h3abionet/chipimputation/issues
- Documentation: https://github.com/h3abionet/chipimputation/wiki