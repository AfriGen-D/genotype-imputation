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
- remember to add it at the end when wverything is working fine

## Instructions for Claude

### Bioinformatics Best Practices
- **Scientific Rigor**: This is a bioinformatics project requiring evidence-based decisions
  - Research relevant literature and established methods before implementing solutions
  - Reference bioinformatics best practices and gold-standard tools
  - Consider biological implications of computational choices
  - Verify parameter choices against published benchmarks when available
  
- **Domain-Specific Considerations**:
  - Understand the biological context (genotype imputation, QC metrics, population genetics)
  - Check for field-specific standards (e.g., GATK best practices, 1000 Genomes guidelines)
  - Consider computational resources and scalability for genomic data
  - Ensure compatibility with standard bioinformatics file formats (VCF, BCF, PLINK, etc.)
  - Account for population-specific considerations in imputation accuracy

### General Development Practices
- Always think deeply and thoroughly analyze problems before providing solutions
- Do comprehensive research and exploration before making changes
- Consider edge cases and potential impacts of any modifications
- Verify assumptions by checking actual code and configurations
- Use extended thinking for complex problems requiring careful analysis

### Workflow Modification Protocol
- **IMPORTANT**: Before editing any Nextflow process, module, or workflow file:
  1. First show the proposed changes with clear explanation of what will be modified
  2. Explain the rationale and potential impacts
  3. Wait for explicit approval before proceeding with the actual edits
  4. This applies to all .nf files, config files, and any workflow-related modifications