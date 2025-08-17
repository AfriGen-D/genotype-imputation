# Migration to nf-core Standards - Complete Guide

## Overview
This document summarizes the migration of the h3achipimputation pipeline to follow nf-core community standards and best practices.

## ✅ Completed Migration Steps

### 1. Directory Structure (nf-core compliant)
```
.
├── .github/
│   └── workflows/
│       └── ci.yml                 # nf-core CI/CD workflow
├── assets/
│   ├── schema_input.json          # Input validation schema
│   ├── samplesheet.csv            # Example samplesheet
│   └── samplesheet_test.csv       # Test samplesheet
├── bin/                           # Custom scripts (existing)
├── conf/
│   ├── base.config                # Base configuration
│   ├── modules.config             # Module-specific settings
│   ├── test_nfcore.config         # nf-core test profile
│   └── ...                        # Other profiles
├── docs/                          # Documentation
├── lib/
│   ├── NfcoreTemplate.groovy     # nf-core utilities
│   ├── WorkflowMain.groovy       # Workflow utilities
│   └── Utils.groovy               # Helper functions
├── modules/
│   └── local/                     # Custom modules
│       ├── impute/
│       │   └── impute_minimac4.nf
│       ├── phasing/
│       │   └── eagle_phasing.nf
│       ├── qc/
│       │   └── check_files.nf
│       └── subset_vcf/
│           └── generate_chunks_vcf.nf
├── subworkflows/
│   ├── local/                     # Custom subworkflows
│   │   ├── preprocess.nf
│   │   ├── phase.nf
│   │   ├── impute.nf
│   │   └── report.nf
│   └── nf-core/                   # nf-core subworkflows (future)
├── workflows/
│   └── chipimputation.nf          # Main workflow
├── main_nfcore.nf                 # nf-core entry point
├── nextflow_nfcore.config         # nf-core configuration
├── CITATIONS.md                   # Tool citations
└── README_nfcore.md               # nf-core README
```

### 2. Key Files Created

#### Main Workflow Files
- **main_nfcore.nf**: Entry point with nf-core validation and parameter handling
- **workflows/chipimputation.nf**: Main workflow logic following nf-core patterns
- **nextflow_nfcore.config**: Configuration with nf-core defaults and profiles

#### Subworkflows (DSL2)
- **subworkflows/local/preprocess.nf**: QC and preprocessing steps
- **subworkflows/local/phase.nf**: Phasing with Eagle2
- **subworkflows/local/impute.nf**: Imputation with Minimac4
- **subworkflows/local/report.nf**: Reporting and visualization

#### Modules (DSL2)
- **modules/local/impute/impute_minimac4.nf**: Minimac4 imputation process
- **modules/local/phasing/eagle_phasing.nf**: Eagle phasing process
- **modules/local/qc/check_files.nf**: File validation process
- **modules/local/subset_vcf/generate_chunks_vcf.nf**: Chunking process

#### Configuration
- **conf/modules.config**: Process-specific configurations
- **conf/test_nfcore.config**: Test profile for CI/CD
- **assets/schema_input.json**: Input validation schema
- **assets/samplesheet.csv**: Example input format

#### CI/CD
- **.github/workflows/ci.yml**: GitHub Actions workflow for testing

#### Documentation
- **CITATIONS.md**: Comprehensive tool citations
- **README_nfcore.md**: nf-core-style documentation

### 3. nf-core Features Implemented

✅ **DSL2 Syntax**
- Modular workflow design
- Reusable processes and subworkflows
- Clear separation of concerns

✅ **Container Support**
- Docker/Singularity/Conda profiles
- BioContainers integration
- Custom containers preserved

✅ **Input Validation**
- JSON schema for samplesheet validation
- nf-validation plugin support
- Parameter validation

✅ **Testing Framework**
- Test profile with minimal dataset
- CI/CD with GitHub Actions
- Multiple parameter testing

✅ **Documentation**
- Comprehensive README
- Citations file
- Parameter documentation

✅ **Best Practices**
- Semantic versioning
- Process labels for resource allocation
- Consistent naming conventions
- Error handling and logging

## 🔄 Migration Path

### For Users
1. Use the new entry point: `nextflow run main_nfcore.nf`
2. Provide input via samplesheet: `--input samplesheet.csv`
3. Use nf-core profiles: `-profile docker,test`

### For Developers
1. Add new processes in `modules/local/`
2. Create subworkflows in `subworkflows/local/`
3. Update `conf/modules.config` for process settings
4. Follow nf-core module naming conventions

## 📝 Usage Examples

### Basic Usage
```bash
# Run with test data
nextflow run main_nfcore.nf -profile test,docker

# Run with real data
nextflow run main_nfcore.nf \
    --input samplesheet.csv \
    --genome GRCh37 \
    --outdir results \
    -profile docker
```

### Custom Parameters
```bash
nextflow run main_nfcore.nf \
    --input samplesheet.csv \
    --chromosomes 21 \
    --chunk_size 3000000 \
    --minRatio 0.01 \
    --qc_plots true \
    -profile singularity
```

## 🔧 Configuration

### Profile Selection
- `docker`: Run with Docker
- `singularity`: Run with Singularity
- `conda`: Run with Conda/Mamba
- `test`: Use test dataset
- `<institute>`: Use institutional configs

### Resource Allocation
Configure in `conf/base.config`:
```groovy
process {
    cpus   = { check_max( 2    * task.attempt, 'cpus'   ) }
    memory = { check_max( 6.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h  * task.attempt, 'time'   ) }
}
```

## 🚀 Benefits of nf-core Migration

1. **Standardization**: Follows community best practices
2. **Maintainability**: Modular design easier to maintain
3. **Reproducibility**: Better container integration
4. **Testing**: Automated CI/CD testing
5. **Documentation**: Comprehensive and standardized
6. **Community**: Access to nf-core tools and support
7. **Portability**: Runs on any nf-core compatible infrastructure

## 📊 Compatibility

### Preserved Features
- All original functionality maintained
- Existing containers still used
- Parameter names unchanged (where possible)
- Output structure preserved

### Enhanced Features
- Input validation
- MultiQC integration ready
- Better error handling
- Improved logging
- Resource optimization

## 🔮 Future Enhancements

1. **Add nf-core modules**:
   - Replace custom modules with nf-core versions where available
   - Contribute custom modules to nf-core

2. **MultiQC Integration**:
   - Add MultiQC custom content
   - Generate comprehensive reports

3. **Wave/Fusion Support**:
   - Enable Wave containers
   - Add Fusion file system support

4. **Enhanced Testing**:
   - Add unit tests for modules
   - Implement full test dataset
   - Add regression testing

5. **Documentation**:
   - Add parameter descriptions
   - Create usage tutorials
   - Add troubleshooting guide

## 📚 Resources

- [nf-core website](https://nf-co.re/)
- [nf-core tools](https://github.com/nf-core/tools)
- [nf-core modules](https://github.com/nf-core/modules)
- [Nextflow documentation](https://www.nextflow.io/docs/latest/)

## ✨ Summary

The pipeline has been successfully migrated to follow nf-core standards while preserving all original functionality. The modular structure makes it easier to maintain, test, and extend. The pipeline is now ready for:

- Community contribution
- Automated testing
- Deployment on any nf-core compatible infrastructure
- Integration with nf-core tools and workflows

For questions or issues, please open an issue on the GitHub repository.