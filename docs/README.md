# h3achipimputation: Documentation

The h3achipimputation documentation is split into the following files:

## Core Documentation
1. [Installation](installation.md)
2. [Pipeline configuration](config_files.md)  
    2.1. [Configuration files](configs.md)  
    2.2. [Software requirements](soft_requirements.md)  
    2.3. [Other clusters](other_clusters.md)  
3. [Running the pipeline](usage.md)
4. [Output and how to interpret the results](output.md) ⭐ **Updated 2025-08-21**
5. [Troubleshooting](troubleshooting.md)

## Quick References
- [Quick Reference Guide](QUICK_REFERENCE.md) 🆕 **New**
- [Documentation Changelog](CHANGELOG.md) 🆕 **New**

## Key Updates (2025-08-21)

### Comprehensive Output Documentation
The pipeline now generates **30+ plot types** across three analysis levels:
- **Chunk-level**: Detailed QC for genomic regions
- **Chromosome-level**: 22 performance dashboards (~41KB each)  
- **Genome-wide**: 8 comprehensive analyses (~2MB total)

### New Plot Types
- **Manhattan plots**: R² by genomic position across all chromosomes
- **MAF vs R² analysis**: Frequency-quality relationships for African genomics
- **SNP count analysis**: Variant density and distribution patterns
- **Cross-chromosome comparisons**: Quality metrics across the genome

### Enhanced Features
- **Hierarchical organization**: Results structured by dataset → reference_panel → analysis_level
- **Quality diagnostics**: File size indicators for plot validation
- **African genomics focus**: Population-specific interpretation guidelines
- **Troubleshooting tools**: Plot regeneration and quality checking utilities

For quick access to key information, see the [Quick Reference Guide](QUICK_REFERENCE.md).
