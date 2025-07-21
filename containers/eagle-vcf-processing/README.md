# Eagle + VCF Processing Container

A comprehensive bioinformatics container that combines Eagle haplotype phasing with VCF processing tools.

## Tools Included

- **Eagle v2.4.1** - Haplotype phasing software from the Broad Institute
- **BCFtools v1.20** - VCF/BCF manipulation and analysis
- **VCFtools v0.1.16** - VCF file processing and analysis
- **HTSlib v1.20** - High-throughput sequencing library (includes tabix)
- **tabix** - Indexing and querying of tab-delimited genome position files

## Usage

### Pull the container
```bash
docker pull mamana/eagle-vcf-processing:eagle-2.4.1
```

### Run with Singularity
```bash
singularity pull docker://mamana/eagle-vcf-processing:eagle-2.4.1
```

### Test the container
```bash
# Check all tools
docker run --rm mamana/eagle-vcf-processing:eagle-2.4.1 eagle_vcf_tools_info

# Test individual tools
docker run --rm mamana/eagle-vcf-processing:eagle-2.4.1 eagle --version
docker run --rm mamana/eagle-vcf-processing:eagle-2.4.1 tabix --version
docker run --rm mamana/eagle-vcf-processing:eagle-2.4.1 bcftools --version
```

## Nextflow Integration

Add to your `nextflow.config`:

```groovy
process {
  withLabel: 'phasing_eagle' {
    container = 'mamana/eagle-vcf-processing:eagle-2.4.1'
  }
}
```

## Use Cases

This container is specifically designed for:
- **minimac4_phasing_eagle** process in genotype imputation workflows
- **impute5_phasing_eagle** process in genotype imputation workflows
- Any workflow requiring both Eagle phasing and VCF processing tools

## Base Image

Alpine Linux 3.19 for minimal size and security.

## License

The container is provided under the same license as the individual tools:
- Eagle: MIT License
- BCFtools: MIT License  
- VCFtools: LGPL v3
- HTSlib: MIT License 