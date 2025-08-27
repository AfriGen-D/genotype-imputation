- not building docker locally, commit and push to GH, GH Actions will build the container.

## Recent Work (Last Updated: 2025-08-21)

### Completed
- Migrated all plotting from R to Python for better compatibility
- Fixed container configurations for Python processes in Nextflow workflow
- Added procps package to python-plotting container for process monitoring
- Implemented chromosome-level reporting for imputation QC
- Simplified GitHub Actions by removing unnecessary security scanning
- Properly assigned containers to processes (imputation, reporting, VCF subsetting)
- **Fixed critical chromosome mismatch bug** where Eagle was using chr1 reference for chr2 targets (removed .first() operators)
- **Resolved TARGET_QC configuration warnings** by removing references from all config files
- **Fixed Nextflow channel deadlock** in preprocess.nf when all chunks pass overlap check
- **Optimized SLURM parallelization** - pipeline now submits ~500 jobs efficiently
- **Implemented comprehensive plotting system** - Complete genome-wide and chromosome-level QC visualization suite
- **Fixed empty plot issue** - Replaced placeholder scripts with data-driven scientific visualizations
- **Created hierarchical results organization** - dataset → reference_panel → (chunks, chromosome, genome, reports)

### Current State
- Python plotting using matplotlib/seaborn in dedicated python-plotting container
- Workflow uses different containers for different process types:
  - `quay.io/h3abionet_org/py3plink` for imputation processes
  - `mamana/python-plotting:1.1.0` for plotting/QC visualization
  - `quay.io/biocontainers/bcftools:1.11--h7c999a4_0` for bcftools operations
  - `sickleinafrica/r-analysis:latest` for any remaining R processes
- **Pipeline successfully processes large datasets** (awigen_500_b38 with ~554 chunks)
- **SLURM optimization active**: queueSize=500, submitRateLimit=200/min, process maxForks configured
- **Complete QC visualization system** with 30+ plot types across chunk, chromosome, and genome levels

### Plotting System Architecture (2025-08-21)

#### **Comprehensive Plot Types**
- **Chunk-level plots** (22 chromosomes): Performance metrics, MAF analysis, variant distributions per chromosome
- **Genome-wide plots** (8 types): Manhattan plots, MAF vs R² analysis, SNP count distributions, cross-chromosome comparisons
- **Hierarchical organization**: Results organized by dataset → reference_panel → analysis_level structure

#### **Key Plot Implementations**
1. **Manhattan Plot** (`genome_r2_position.pdf`): R² values by genomic position across all chromosomes
2. **MAF vs R² Scatter** (`genome_maf_r2.pdf`): Frequency-quality relationship analysis with population genetics insights
3. **SNP Count Analysis** (`genome_r2_snpcount.pdf`): R² distribution and quality metrics by variant density
4. **Chromosome Performance** (`chr_performance.pdf`): Multi-panel QC dashboards per chromosome
5. **Performance Summaries** (`genome_performance.pdf`): Cross-chromosome quality comparisons

#### **Results Directory Structure**
```
/scratch3/users/mamana/results/
└── awigen_500_b38/
    └── H3AR6x/
        ├── chunks/
        │   └── chr{N}/
        │       ├── plots/          # Chunk-level QC plots
        │       ├── frequency/      # Allele frequency analysis
        │       └── warnings/       # QC warnings and logs
        ├── chromosome/
        │   ├── plots/              # 22 chromosome performance plots (~41KB each)
        │   └── *.chr_summary.json  # Chromosome aggregated metrics
        └── genome/
            ├── plots/              # 8 genome-wide plots (total ~2MB)
            └── *.genome_summary.json # Genome-wide aggregated metrics
```

#### **Plot Generation Tools**
- **`regenerate_all_plots.py`**: Systematic plot regeneration using proper data-driven scripts
- **Data-driven scripts**: 330+ line implementations with comprehensive statistical analysis
- **Container-based execution**: All plots generated in `mamana/python-plotting:1.1.0` container
- **Scientific visualization**: Matplotlib/seaborn with publication-quality outputs

### Performance Optimizations (v6_chr21_phased_nfcore.config)
```nextflow
// SLURM executor settings
executor.queueSize = 500           // Allow 500 concurrent jobs
executor.submitRateLimit = '200/1min'  // Submit up to 200 jobs/minute

// Process-specific parallelization
withName: 'EAGLE_PHASING' { maxForks = 100 }
withName: 'IMPUTE_MINIMAC4' { maxForks = 100 }
withName: 'QC_DUPL|SPLIT_MULTI_ALLELIC|FILTER_MIN_AC' { maxForks = 50 }
```

### Known Issues/TODOs
- **COMPARE_PRE_POST_IMPUTATION module** - disabled due to syntax errors (needs fixing)
- **PLOT_R2_GENOMIC_WINDOWS module** - disabled due to Unicode R² character issues
- Monitor container builds on GitHub Actions after commits
- All plotting scripts verified working with assigned containers (✅ Completed)

### Plotting Infrastructure Details

#### **Script Categories**
- **Real plotting scripts** (300+ lines): Comprehensive data-driven visualizations with multi-panel layouts
  - Examples: `plot_chr_performance.py`, `plot_genome_summary.py`, `plot_genome_maf_r2.py`
  - Generate publication-quality scientific plots with statistical analysis
- **Placeholder scripts** (56 lines): Simple text-based outputs (legacy, mostly replaced)
  - Generate basic "Placeholder Plot" PDFs with minimal content

#### **Plot Quality Indicators**
- **Large file sizes** (40KB-1.3MB): Indicate real data visualization with detailed graphics
- **Small file sizes** (~13KB): Usually placeholder plots with static text
- **Multi-panel layouts**: Real scripts generate 2x2 or 2x3 subplot arrangements with statistical summaries

#### **Data Flow**
1. **Chunk reports** → JSON summaries → **Chromosome aggregation** → **Genome aggregation**
2. **VCF parsing** → Minimac4 format handling → R²/MAF extraction → Plot generation
3. **Container execution** → Python matplotlib/seaborn → PDF output → Organized directory structure

#### **Regeneration Process**
- Use `regenerate_all_plots.py` for systematic plot creation
- Automatically identifies and runs proper data-driven scripts
- Handles both chromosome-level (22 plots) and genome-level (8 plots) generation
- Container-based execution ensures reproducible environments

#### **African Genomics Considerations**
- MAF distributions reflect African population structure
- Reference panel coverage patterns specific to H3Africa data
- R² thresholds calibrated for population-specific imputation accuracy
- Frequency spectrum analysis accounts for African variant diversity

## Recent Updates (Last Updated: 2025-08-27)

### Pipeline Testing Results
- **Core imputation workflow functional**: Successfully tested with chr21 test data
- **Phasing and imputation modules working**: EAGLE and MINIMAC4 completed successfully
- **Reporting modules have issues**: Several modules expect different input formats
  - PLOT_CALIBRATION, PLOT_CONCORDANCE_MAF, PLOT_CROSS_VALIDATION disabled (require validation data)
  - PLOT_AGGREGATED_R2_DASHBOARD expects .info files but pipeline produces .sites.vcf.gz
  - Added gzip support to plotting scripts for compressed VCF handling
- **Main outputs successfully created**: 
  - Phased VCFs (.phased.vcf.gz)
  - Imputed dosage VCFs (.dose.vcf.gz) 
  - Sites info files (.sites.vcf.gz)

## Instructions for Claude

### Documentation Maintenance
- **CRITICAL**: After major code updates, commits, or significant changes, ALWAYS update this CLAUDE.md file
- Document new patterns, conventions, or project-specific requirements discovered during work
- Update instructions when project structure or dependencies change significantly  
- Add warnings about known issues or areas requiring special attention
- Record important architectural decisions and their rationale
- Keep track of coding standards and style preferences observed in the codebase
- Before making commits, review if CLAUDE.md needs updating with any discoveries or fixes

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