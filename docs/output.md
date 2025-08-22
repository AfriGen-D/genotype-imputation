# h3achipimputation: Output

## Table of Contents

* [Introduction](#introduction)
* [Pipeline Output Structure](#pipeline-output-structure)
* [Hierarchical Results Organization](#hierarchical-results-organization)
* [Quality Control Plots](#quality-control-plots)
* [Data Files and Reports](#data-files-and-reports)
* [Interpreting Results](#interpreting-results)
* [Troubleshooting Output Issues](#troubleshooting-output-issues)

## Introduction

The h3achipimputation pipeline produces a comprehensive set of outputs organized in a hierarchical structure that facilitates analysis at multiple scales: chunk-level, chromosome-level, and genome-wide. The pipeline generates imputed VCF files, extensive quality control visualizations, and detailed performance metrics.

## Pipeline Output Structure

### Default Output Directory
By default, the pipeline writes results to `./output/`, but this can be changed using the `--outdir` parameter.

### Hierarchical Organization (Updated 2025-08-21)
Results are organized in a three-tier hierarchical structure:

```
output/
└── {dataset}/
    └── {reference_panel}/
        ├── chunks/
        │   └── chr{N}/
        │       ├── plots/              # Chunk-level QC plots
        │       ├── frequency/          # Allele frequency analysis
        │       └── warnings/           # QC warnings and logs
        ├── chromosome/
        │   ├── plots/                  # Chromosome-level plots (22 files)
        │   ├── *.chr_summary.json      # Chromosome aggregated metrics
        │   └── *.chr_stats.txt         # Chromosome statistics
        └── genome/
            ├── plots/                  # Genome-wide plots (8 files)
            ├── *.genome_summary.json   # Genome-wide aggregated metrics
            └── *.genome_stats.txt      # Genome-wide statistics
```

**Example Structure:**
```
output/
└── awigen_500_b38/
    └── H3AR6x/
        ├── chunks/
        │   ├── chr1/plots/
        │   ├── chr2/plots/
        │   └── ...
        ├── chromosome/plots/    # 22 performance plots (~41KB each)
        └── genome/plots/        # 8 comprehensive plots (~2MB total)
```

## Quality Control Plots

### Chunk-Level Plots
Generated for each chromosome and genomic chunk:

* **Performance plots**: R² metrics, imputation accuracy, variant quality
* **Frequency analysis**: Allele frequency comparisons between target and reference
* **MAF distributions**: Minor allele frequency patterns
* **Warning reports**: QC flags and potential issues

### Chromosome-Level Plots (22 files)
Each chromosome generates a comprehensive performance dashboard:

**`{dataset}_chr{N}_{ref_panel}.chr_performance.pdf`** (~41KB each)
- **Chunk performance comparison**: R², concordance, well-imputed rates across chunks
- **Variant distribution**: Genotyped vs. well-imputed vs. other variants
- **MAF distribution**: Frequency spectrum analysis
- **R² by MAF categories**: Quality metrics stratified by frequency
- **Summary statistics**: Key metrics dashboard

### Genome-Wide Plots (8 files, ~2MB total)

#### 1. Comprehensive Summary
**`{dataset}_{ref_panel}.genome_plots.pdf`** (57KB)
- Multi-panel overview of genome-wide imputation performance
- Cross-chromosome comparisons and statistical summaries

#### 2. Performance Analysis
**`{dataset}_{ref_panel}.genome_performance.pdf`** (38KB)
- Mean R² by chromosome with genome-wide average
- Total variants per chromosome (millions)
- Well-imputed rates by chromosome
- Performance statistics dashboard

#### 3. MAF Analysis
**`{dataset}_{ref_panel}.genome_maf_analysis.pdf`** (34KB)
- Genome-wide MAF distribution (linear and log scale)
- Mean R² by MAF category with sample sizes
- Rare vs. common variant breakdown
- MAF-stratified quality metrics

#### 4. R² Distribution Analysis
**`{dataset}_{ref_panel}.genome_r2_distribution.pdf`** (42KB)
- Genome-wide R² histogram with quality thresholds
- R² quartiles by chromosome (median ± IQR)
- Cumulative R² distribution with percentile markers
- Quality category breakdown (high/medium/low)

#### 5. Manhattan Plot
**`{dataset}_{ref_panel}.genome_r2_position.pdf`** (1.3MB)
- **R² values by genomic position** across all chromosomes
- Quality threshold lines (R² = 0.3, 0.8, genome mean)
- Chromosome boundaries and labels
- R² distribution histogram

#### 6. MAF vs R² Relationship
**`{dataset}_{ref_panel}.genome_maf_r2.pdf`** (543KB)
- Scatter plot showing MAF vs R² colored by chromosome
- Trend analysis with correlation coefficients
- MAF distribution and R² by frequency categories
- Population genetics insights for African genomics

#### 7. SNP Count Analysis
**`{dataset}_{ref_panel}.genome_r2_snpcount.pdf`** (48KB)
- R² vs total SNP count by chromosome
- R² distribution and cumulative analysis
- Variants per chunk analysis
- Quality categories summary

#### 8. Chromosome Comparison
**`{dataset}_{ref_panel}.genome_chr_comparison.pdf`** (39KB)
- R² vs chromosome size scatter analysis
- Normalized metrics comparison
- Chromosome heterogeneity analysis
- Quality distribution across chromosomes

## Data Files and Reports

### JSON Summary Files
**Chromosome Level:** `{dataset}_chr{N}_{ref_panel}.chr_summary.json`
```json
{
  "chromosome": "chr1",
  "chunks_processed": 10,
  "total_variants": 4650164,
  "well_imputed_variants": 2112597,
  "mean_r2": 0.42982433477614984,
  "maf_bins": {...},
  "chunk_details": [...]
}
```

**Genome Level:** `{dataset}_{ref_panel}.genome_summary.json`
```json
{
  "dataset": "awigen_500_b38",
  "ref_name": "H3AR6x",
  "chromosomes_processed": 22,
  "total_variants": 58655394,
  "well_imputed_variants": 26956088,
  "mean_r2": 0.43508862449888236,
  "chromosome_details": [...]
}
```

### Imputed VCF Files
**Location:** Distributed across chunk and chromosome directories
**Format:** Standard VCF 4.2 with imputation quality scores
**Content:** 
- Imputed genotypes with dosage information
- R² quality scores in INFO field
- Allele frequency annotations

### Statistical Reports
**Text summaries:** `*.stats.txt` files containing:
- Variant counts by category
- Quality metric distributions
- Imputation performance summaries
- Reference panel coverage statistics

## Interpreting Results

### Quality Thresholds
- **R² ≥ 0.8**: High quality imputation, suitable for most analyses
- **R² ≥ 0.3**: Well-imputed threshold, commonly used cutoff
- **R² < 0.3**: Poor imputation quality, consider filtering

### MAF Considerations (African Genomics)
- **Very rare (MAF < 0.01)**: Limited imputation accuracy due to reference panel size
- **Rare (0.01-0.05)**: Moderate quality, population-specific patterns
- **Low frequency (0.05-0.2)**: Good imputation quality
- **Common (≥ 0.2)**: Excellent imputation accuracy

### Population-Specific Patterns
The pipeline is optimized for African genomics with H3Africa reference panels:
- MAF distributions reflect African population structure
- R² thresholds calibrated for African genetic diversity
- Reference panel coverage patterns specific to African populations

## Plot Quality Indicators

### File Size Diagnostics
- **Large files (40KB-1.3MB)**: Real data visualizations with detailed graphics
- **Small files (~13KB)**: Potential placeholder plots or minimal data
- **Multi-panel layouts**: Comprehensive scientific visualizations

### Visual Quality Checks
- **Manhattan plots**: Should show clear chromosomal patterns
- **Scatter plots**: Should display meaningful data clouds with trends
- **Histograms**: Should show realistic distribution shapes
- **Statistical panels**: Should contain numerical summaries

## Troubleshooting Output Issues

### Empty or Missing Plots
If plots appear empty or show placeholder text:
1. Check file sizes - small files (~13KB) may be placeholders
2. Use the plot regeneration tool: `python3 bin/regenerate_all_plots.py`
3. Verify container execution: ensure `mamana/python-plotting:1.1.0` is available
4. Check input data: ensure JSON summary files contain valid data

### Plot Regeneration
```bash
# Regenerate all plots with proper data-driven scripts
cd /path/to/results
python3 /path/to/pipeline/bin/regenerate_all_plots.py \
  --dataset your_dataset \
  --ref-name your_reference \
  --results-dir /path/to/results
```

### Common Issues
- **Container access**: Ensure Singularity/Docker can access plotting containers
- **Data format**: Verify VCF files follow Minimac4 output format
- **File permissions**: Check read/write access to output directories
- **Memory limits**: Large genome-wide plots may require sufficient memory

### Getting Help
For technical issues or questions about output interpretation:
1. Check the troubleshooting section in the main documentation
2. Review the CLAUDE.md file for implementation details
3. Submit issues to the project repository with example output files