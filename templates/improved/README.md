# Improved Imputation Pipeline Scripts

This directory contains modernized, improved versions of the imputation pipeline scripts with better error handling, type hints, logging, and Python 3 compatibility.

## Features

### General Improvements
- ✅ **Python 3 Compatible**: All scripts use Python 3 with type hints
- ✅ **Better Error Handling**: Comprehensive exception handling with informative error messages
- ✅ **Logging**: Structured logging with configurable verbosity levels
- ✅ **Input Validation**: Robust validation of input files and parameters
- ✅ **Modular Design**: Object-oriented design with reusable components
- ✅ **Documentation**: Comprehensive docstrings and help messages
- ✅ **Performance**: Optimized for large-scale genomic data processing
- ✅ **Testing**: Includes validation and statistics reporting

## Scripts Overview

### 1. `filter_info_minimac.py`
Filters imputation info files based on quality thresholds.

**Key Features:**
- Process multiple info files in batch
- Configurable R² threshold for well-imputed variants
- MAF filtering
- Concordance data support
- Detailed summary statistics
- Progress logging

**Usage:**
```bash
python filter_info_minimac.py \
    --infoFiles file1.info.gz,file2.info.gz \
    --datasets dataset1,dataset2 \
    --out_prefix output_prefix \
    --infoCutoff 0.8
```

### 2. `generate_chunks.py`
Generates genomic chunks for parallel imputation processing.

**Key Features:**
- Flexible chunk size configuration
- Support for specific chromosomes or regions
- Natural chromosome sorting (1, 2, ..., 22, X, Y, MT)
- BED format output option
- Statistics reporting
- Region validation against actual data

**Usage:**
```bash
# Generate chunks for all chromosomes
python generate_chunks.py input.map chunks.txt --chunk-size 5000000

# Generate chunks for specific chromosomes
python generate_chunks.py input.map chunks.txt --chromosomes 1,2,3

# Generate chunks for specific region
python generate_chunks.py input.map chunks.txt --region chr1:1000000-5000000

# Get statistics about the data
python generate_chunks.py input.map chunks.txt --stats
```

### 3. `report_well_imputed.py`
Generates comprehensive reports on well-imputed variants.

**Key Features:**
- Multi-dataset support
- MAF and R² distribution analysis
- Customizable thresholds
- Multiple output formats (text, JSON)
- Variant list generation
- Per-dataset and combined statistics

**Usage:**
```bash
python report_well_imputed.py \
    --info-files file1.info,file2.info \
    --datasets data1,data2 \
    --output-prefix report \
    --rsq-threshold 0.8 \
    --maf-threshold 0.01 \
    --json \
    --variant-list
```

**Outputs:**
- `.summary.txt`: Human-readable summary report
- `.variants.txt`: List of well-imputed variant IDs
- `.stats.json`: Machine-readable statistics

### 4. `report_accuracy_by_maf.py`
Analyzes imputation accuracy stratified by minor allele frequency.

**Key Features:**
- Standard MAF bins for imputation studies
- Custom MAF bin support
- Mean and standard deviation of R²
- Empirical R² support
- Publication-ready tables
- JSON export for downstream analysis

**Usage:**
```bash
python report_accuracy_by_maf.py \
    --info-files file1.info,file2.info \
    --datasets data1,data2 \
    --output-prefix accuracy \
    --rsq-threshold 0.8 \
    --summary-table \
    --json

# With custom MAF bins
python report_accuracy_by_maf.py \
    --info-files file.info \
    --datasets data1 \
    --output-prefix accuracy \
    --custom-maf-bins 0,0.01,0.05,0.1,0.5
```

**Outputs:**
- `.report.txt`: Detailed accuracy report
- `.summary.tsv`: Tab-separated summary table
- `.metrics.json`: JSON format metrics

### 5. `extract_region_from_refs.py`
Extracts specific genomic regions from reference panel files.

**Key Features:**
- BCFtools integration
- Automatic indexing
- Sample filtering (include/exclude)
- Allele frequency filtering
- Multiple region extraction
- Region merging
- Support for VCF/BCF formats

**Usage:**
```bash
# Extract single region
python extract_region_from_refs.py \
    reference.vcf.gz \
    chr1:1000000-2000000 \
    output.vcf.gz

# Extract with sample filtering
python extract_region_from_refs.py \
    reference.vcf.gz \
    chr1:1000000-2000000 \
    output.vcf.gz \
    --samples sample1,sample2,sample3

# Extract with MAF filtering
python extract_region_from_refs.py \
    reference.vcf.gz \
    chr1:1000000-2000000 \
    output.vcf.gz \
    --min-af 0.01 \
    --max-af 0.5

# Extract multiple regions from file
python extract_region_from_refs.py \
    reference.vcf.gz \
    regions.txt \
    output_dir/ \
    --merge
```

## Installation Requirements

```bash
# Python packages (optional but recommended)
pip install numpy  # For statistical calculations

# System tools (required for extract_region_from_refs.py)
# Install bcftools and tabix
conda install -c bioconda bcftools tabix
# or
apt-get install bcftools tabix
```

## Integration with Nextflow Pipeline

To use these improved scripts in the Nextflow pipeline, update the template references:

```nextflow
process filter_info_by_target {
    script:
    template 'improved/filter_info_minimac.py'
}

process generate_chunks_vcf {
    script:
    template 'improved/generate_chunks.py'
}
```

## Error Handling

All scripts include comprehensive error handling:

- **File not found**: Clear error messages with file paths
- **Invalid input**: Validation with helpful error messages
- **Processing errors**: Graceful handling with recovery options
- **Progress tracking**: Regular status updates for long-running operations

## Performance Considerations

### Large Files
- Scripts use streaming processing where possible
- Memory-efficient data structures
- Progress logging for long operations

### Parallel Processing
- Scripts are designed to be run in parallel
- Thread-safe file operations
- Atomic write operations for output files

## Testing

Test the scripts with sample data:

```bash
# Generate test data
echo -e "chr1\t1000000\nchr1\t2000000\nchr1\t3000000" > test.map

# Test chunk generation
python generate_chunks.py test.map test_chunks.txt --chunk-size 1000000 --stats

# Test with actual info files
python filter_info_minimac.py \
    --infoFiles test.info \
    --datasets test \
    --out_prefix test_filtered \
    --infoCutoff 0.3
```

## Logging Levels

All scripts support different logging levels:

```bash
# Default (INFO level)
python script.py [arguments]

# Verbose output (DEBUG level)
python script.py [arguments] --verbose

# Suppress non-error output
python script.py [arguments] 2>/dev/null
```

## Output Formats

### JSON Output
Machine-readable format for integration with other tools:
```json
{
  "parameters": {
    "rsq_threshold": 0.8,
    "maf_threshold": 0.01
  },
  "datasets": {
    "dataset1": {
      "total_variants": 1000000,
      "well_imputed": 850000,
      "well_imputed_rate": 0.85
    }
  }
}
```

### TSV Output
Tab-separated format for spreadsheet import and R/Python analysis:
```
Dataset	MAF_Range	N_Variants	Mean_Rsq	SD_Rsq	Well_Imputed_Pct
data1	0-0.1%	1000	0.3421	0.1234	15.32
data1	0.1-0.5%	5000	0.5678	0.1456	42.15
```

## Advantages Over Original Scripts

1. **Reliability**: Better error handling prevents pipeline crashes
2. **Maintainability**: Clean, documented code with type hints
3. **Performance**: Optimized algorithms and data structures
4. **Usability**: Better command-line interfaces with help text
5. **Compatibility**: Python 3 support ensures long-term viability
6. **Extensibility**: Modular design allows easy feature additions
7. **Debugging**: Comprehensive logging aids troubleshooting

## Contributing

When modifying these scripts:

1. Maintain Python 3 compatibility
2. Add type hints for new functions
3. Include docstrings for documentation
4. Handle errors gracefully
5. Add logging for important operations
6. Update this README with new features

## License

These scripts are part of the h3achipimputation pipeline and follow the same license terms.