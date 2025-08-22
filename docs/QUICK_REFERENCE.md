# h3achipimputation: Quick Reference Guide

## Output Structure at a Glance

```
output/{dataset}/{reference_panel}/
├── chunks/chr{N}/plots/     # Detailed chunk-level QC
├── chromosome/plots/        # 22 chromosome summaries (~41KB each)  
└── genome/plots/           # 8 genome-wide analyses (~2MB total)
```

## Key Quality Metrics

### R² Thresholds
- **≥ 0.8**: High quality - suitable for most analyses
- **≥ 0.3**: Well-imputed threshold - commonly used cutoff  
- **< 0.3**: Poor quality - consider filtering

### MAF Categories (African Genomics)
- **< 0.01**: Very rare - limited accuracy
- **0.01-0.05**: Rare - moderate quality
- **0.05-0.2**: Low frequency - good quality
- **≥ 0.2**: Common - excellent accuracy

## Plot Quick Reference

### Genome-Wide Plots (8 files)
1. **`genome_plots.pdf`** (57KB) - Overview dashboard
2. **`genome_performance.pdf`** (38KB) - Cross-chromosome metrics  
3. **`genome_maf_analysis.pdf`** (34KB) - Frequency analysis
4. **`genome_r2_distribution.pdf`** (42KB) - Quality distributions
5. **`genome_r2_position.pdf`** (1.3MB) - Manhattan plot 🌟
6. **`genome_maf_r2.pdf`** (543KB) - Frequency vs quality
7. **`genome_r2_snpcount.pdf`** (48KB) - Variant density analysis  
8. **`genome_chr_comparison.pdf`** (39KB) - Chromosome comparisons

### Chromosome Plots (22 files)
- **`chr{N}_performance.pdf`** (~41KB each) - Multi-panel QC dashboard

## Troubleshooting

### Empty Plots?
```bash
# Check file sizes - small files (~13KB) may be placeholders
ls -la output/*/plots/

# Regenerate with proper scripts
python3 bin/regenerate_all_plots.py --dataset your_data --ref-name your_ref
```

### Quality Diagnostics
- ✅ **Large files (>40KB)**: Real data visualization
- ❌ **Small files (~13KB)**: Likely placeholders
- ✅ **Multi-panel layouts**: Comprehensive analysis
- ❌ **Static text**: Placeholder content

## JSON Data Access

### Chromosome Summary
```bash
# Quick stats for chromosome 1
cat output/dataset/ref/chromosome/dataset_chr1_ref.chr_summary.json | jq '.mean_r2'
```

### Genome Summary  
```bash
# Overall statistics
cat output/dataset/ref/genome/dataset_ref.genome_summary.json | jq '.total_variants'
```

## African Genomics Notes

- MAF patterns reflect population structure
- R² thresholds calibrated for H3Africa panels
- Reference coverage varies by chromosome
- Population-specific variant frequencies expected

## Container Usage

All plots generated using:
```bash
singularity exec mamana/python-plotting:1.1.0 python3 script.py
```

## Common Commands

```bash
# View plot directory structure
tree output/*/plots/ -L 2

# Check plot file sizes  
find output -name "*.pdf" -exec ls -lh {} \; | sort -k5

# Count total plots
find output -name "*.pdf" | wc -l

# Regenerate specific genome plots
cd output/dataset/ref/genome
singularity exec container python3 /path/to/script.py --genome-summary file.json
```

For detailed documentation: [docs/output.md](output.md)