# Samplesheet Format

The pipeline requires a CSV samplesheet as input that specifies the VCF files to be processed.

## Format

The samplesheet must be a comma-separated file with the following columns:

| Column      | Description                                          | Required |
|-------------|------------------------------------------------------|----------|
| `sample`    | Unique sample identifier (no spaces allowed)        | Yes      |
| `vcf`       | Full path to VCF file (.vcf or .vcf.gz)            | Yes      |
| `population`| Population code (AFR/AMR/EAS/EUR/SAS/Unknown)      | No       |
| `sex`       | Sample sex (Male/Female/Unknown)                    | No       |

## Example

```csv
sample,vcf,population,sex
sample1,/path/to/sample1.vcf.gz,EUR,Female
sample2,/path/to/sample2.vcf.gz,AFR,Male
sample3,/path/to/sample3.vcf.gz,AMR,Unknown
```

## Templates

- **Template**: `assets/samplesheet_template.csv`
- **Schema**: `assets/schema_input.json`
- **Example for chr21**: `samplesheet_chr21.csv`

## Creating a Samplesheet

### From existing VCF files

If you have VCF files following the old `target_datasets` format, create a samplesheet like this:

```bash
# Create header
echo "sample,vcf,population,sex" > samplesheet.csv

# Add your samples
echo "my_sample,/path/to/my_sample.vcf.gz,AFR,Unknown" >> samplesheet.csv
```

### From the old config format

If migrating from the old pipeline format where you had:
```groovy
target_datasets = [
  ["sample_name", "/path/to/vcf.gz"]
]
```

Convert to:
```csv
sample,vcf,population,sex
sample_name,/path/to/vcf.gz,Unknown,Unknown
```

## Validation

The pipeline will validate your samplesheet and check:
- Sample names are unique and contain no spaces
- VCF files exist at the specified paths
- Population codes are valid (if provided)
- Sex values are valid (if provided)

## Multiple Samples

You can process multiple samples in a single run:

```csv
sample,vcf,population,sex
batch1_sample1,/data/batch1/sample1.vcf.gz,EUR,Female
batch1_sample2,/data/batch1/sample2.vcf.gz,EUR,Male
batch2_sample1,/data/batch2/sample1.vcf.gz,AFR,Female
batch2_sample2,/data/batch2/sample2.vcf.gz,AFR,Male
```

## Using with the Pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --reference_panels panels.json \
  -profile singularity \
  -c v6_chr21_phased_v2.config
```