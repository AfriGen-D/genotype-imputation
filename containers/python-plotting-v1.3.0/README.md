# Python Plotting Container v1.3.0

## Overview
This container provides a Python environment optimized for generating QC plots and reports for the genotype imputation pipeline. It includes all necessary libraries for scientific visualization, VCF processing, and statistical analysis.

## Version Information
- **Container Version**: 1.3.0
- **Base Image**: python:3.9-slim
- **Namespace**: mamana/python-plotting

## Included Packages

### Core Plotting Libraries
- **matplotlib** (3.7.1): Comprehensive plotting library
- **seaborn** (0.12.2): Statistical data visualization
- **plotly** (5.15.0): Interactive visualizations

### Data Processing
- **pandas** (2.0.3): Data manipulation and analysis
- **numpy** (1.24.3): Numerical computing
- **scipy** (1.11.1): Scientific computing

### Genomics Tools
- **pysam** (0.21.0): SAM/BAM/VCF file processing
- **cyvcf2** (0.30.22): Fast VCF parsing

### Machine Learning
- **scikit-learn** (1.3.0): Machine learning algorithms for QC metrics

## Usage

### Building the Container
```bash
# Build locally
./build.sh

# Build and push to Docker Hub
./build.sh --push
```

### Running the Container

#### With Docker
```bash
docker run -it -v $(pwd):/work mamana/python-plotting:1.3.0 python script.py
```

#### With Singularity
```bash
singularity exec docker://mamana/python-plotting:1.3.0 python script.py
```

#### In Nextflow Pipeline
```nextflow
process PLOT_QC_METRICS {
    container 'mamana/python-plotting:1.3.0'
    
    input:
    path vcf_file
    
    output:
    path "*.pdf"
    
    script:
    """
    python plot_qc_metrics.py ${vcf_file}
    """
}
```

## Key Features

### VCF Processing Capabilities
- Efficient parsing of large VCF/BCF files
- Support for compressed formats (gzip, bgzip)
- Streaming processing for memory efficiency

### Plot Types Supported
- Manhattan plots
- MAF vs R² scatter plots
- Chromosome-level QC dashboards
- SNP density distributions
- Heterozygosity analysis
- HWE deviation plots
- Dosage distributions
- Cross-validation metrics

### Performance Optimizations
- Memory-efficient VCF streaming with cyvcf2
- Vectorized operations with numpy
- Multiprocessing support via Python's built-in libraries

## Integration with Genotype Imputation Pipeline

This container is designed to work seamlessly with the H3ABioNet genotype imputation pipeline:

1. **Chunk-level Analysis**: Process individual imputation chunks
2. **Chromosome Aggregation**: Combine chunk metrics per chromosome
3. **Genome-wide Summary**: Generate whole-genome visualizations
4. **Report Generation**: Create comprehensive PDF reports

## Troubleshooting

### Common Issues

1. **Memory Errors with Large VCFs**
   - Use streaming mode with cyvcf2
   - Process in chunks using pandas chunking
   - Increase container memory allocation

2. **Missing System Libraries**
   - Container includes procps for process monitoring
   - All compilation tools included for pip installations

3. **Plot Display Issues**
   - Container uses Agg backend by default (no display required)
   - Outputs saved directly to PDF/PNG files

## Changelog

### v1.3.0 (Current)
- Updated all packages to latest stable versions
- Added cyvcf2 for faster VCF processing
- Improved memory efficiency for large datasets
- Added plotly for interactive visualizations

### v1.2.0
- Added scikit-learn for advanced QC metrics
- Updated matplotlib and seaborn versions
- Fixed compatibility issues with Python 3.9

### v1.1.0
- Initial production release
- Core plotting libraries included
- Basic VCF processing with pysam

## Support

For issues or questions:
- GitHub Issues: https://github.com/h3abionet/chipimputation/issues
- Container Registry: https://hub.docker.com/r/mamana/python-plotting

## License

This container is part of the H3ABioNet Chip Imputation pipeline and is distributed under the MIT License.