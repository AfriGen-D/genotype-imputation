- not building docker locally, commit and push to GH, GH Actions will build the container.

## Recent Work (Last Updated: 2025-08-13)

### Completed
- Migrated all plotting from R to Python for better compatibility
- Fixed container configurations for Python processes in Nextflow workflow
- Added procps package to python-plotting container for process monitoring
- Implemented chromosome-level reporting for imputation QC
- Simplified GitHub Actions by removing unnecessary security scanning
- Properly assigned containers to processes (imputation, reporting, VCF subsetting)

### Current State
- Python plotting using matplotlib/seaborn in dedicated python-plotting container
- Workflow uses different containers for different process types:
  - `quay.io/h3abionet_org/py3plink` for imputation processes
  - `mamana/python-plotting:latest` for plotting/QC visualization
  - `quay.io/biocontainers/bcftools:1.11--h7c999a4_0` for bcftools operations
  - `sickleinafrica/r-analysis:latest` for any remaining R processes

### Known Issues/TODOs
- Monitor container builds on GitHub Actions after commits
- Verify all Python scripts work correctly with their assigned containers
- Always do deep analysis and research
- run test with nextflow run main.nf --input samplesheet_chr21.csv -c v6_chr21_phased_v2.config
      -profile slurm,singularity -resume