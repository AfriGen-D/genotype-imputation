# Nextflow Pipeline Comparison: ImputationServer2 vs H3ABioNet Genotype Imputation

## Executive Summary

This document provides a detailed comparison between the ImputationServer2 Nextflow pipeline (v2.0.8) and the H3ABioNet genotype imputation pipeline, focusing on their Nextflow implementations, architectural patterns, and workflow orchestration approaches. Both pipelines use Nextflow DSL2 but with significantly different design philosophies and capabilities.

---

## Table of Contents

1. [Architecture Philosophy](#architecture-philosophy)
2. [Workflow Structure](#workflow-structure)
3. [Module Organization](#module-organization)
4. [Quality Control Implementation](#quality-control-implementation)
5. [Phasing Approaches](#phasing-approaches)
6. [Imputation Implementation](#imputation-implementation)
7. [Configuration Management](#configuration-management)
8. [Resource Management](#resource-management)
9. [Container Strategy](#container-strategy)
10. [Error Handling & Recovery](#error-handling--recovery)
11. [Reporting & Outputs](#reporting--outputs)
12. [Extended Features](#extended-features)
13. [Development Practices](#development-practices)
14. [Performance Optimization](#performance-optimization)
15. [Use Case Comparison](#use-case-comparison)

---

## 1. Architecture Philosophy

### ImputationServer2 Nextflow

**Design Philosophy**: **Simplified Pipeline-as-a-Service**
- Mirrors web service functionality in Nextflow
- Fixed workflow with optional steps
- Minimal configuration exposure
- Standardized processing paths
- Service-oriented architecture

```nextflow
// ImputationServer2 approach - simplified, linear
workflow {
    INPUT_VALIDATION()
    QUALITY_CONTROL()
    if (params.mode == 'imputation') {
        PHASING()
        IMPUTATION()
        if (params.merge_results) {
            ENCRYPTION()
        }
    }
}
```

### H3ABioNet Pipeline

**Design Philosophy**: **Flexible Research Framework**
- Modular, composable workflows
- Extensive parameterization
- Research-oriented flexibility
- Population-specific adaptations
- Framework architecture

```nextflow
// H3ABioNet approach - modular, adaptive
workflow CHIPIMPUTATION {
    PREPROCESS(ch_input)  // Multi-stage QC
    PHASE(PREPROCESS.out.vcf)  // Multiple methods
    IMPUTE(PHASE.out.phased)  // Multi-panel support
    REPORT(IMPUTE.out.imputed)  // 30+ modules
    REPORT_AGG(REPORT.out.chunk_reports)  // Hierarchical
}
```

### Philosophical Comparison

| Aspect | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Primary Goal** | Standardization | Flexibility |
| **Target Users** | Service users | Researchers |
| **Complexity** | Hidden | Exposed |
| **Customization** | Limited | Extensive |
| **Learning Curve** | Minimal | Significant |

---

## 2. Workflow Structure

### ImputationServer2

**Linear Workflow Pattern**:
```
┌──────────────┐
│    INPUT     │
│  VALIDATION  │
└──────┬───────┘
       ↓
┌──────────────┐
│   QUALITY    │
│   CONTROL    │
└──────┬───────┘
       ↓
┌──────────────┐
│   PHASING    │ (Optional)
└──────┬───────┘
       ↓
┌──────────────┐
│  IMPUTATION  │
└──────┬───────┘
       ↓
┌──────────────┐
│  ENCRYPTION  │ (Optional)
└──────┬───────┘
       ↓
┌──────────────┐
│   ANCESTRY   │ (Optional)
└──────┬───────┘
       ↓
┌──────────────┐
│     PGS      │ (Optional)
└──────────────┘
```

**Workflow Components**:
- 7 main modules
- Sequential execution
- Limited branching
- Optional steps via flags

### H3ABioNet Pipeline

**Modular Workflow Pattern**:
```
┌─────────────────────────────────────┐
│         CHIPIMPUTATION              │
├─────────────────────────────────────┤
│  ┌──────────────────────────────┐  │
│  │      PREPROCESS              │  │
│  │  ┌─────┐ ┌─────┐ ┌─────┐   │  │
│  │  │CHECK│→│CHUNK│→│  QC  │   │  │
│  │  └─────┘ └─────┘ └─────┘   │  │
│  └──────────────────────────────┘  │
│  ┌──────────────────────────────┐  │
│  │         PHASE                │  │
│  │  ┌─────────┐ ┌──────────┐   │  │
│  │  │ EAGLE   │ │ SHAPEIT4 │   │  │
│  │  └─────────┘ └──────────┘   │  │
│  └──────────────────────────────┘  │
│  ┌──────────────────────────────┐  │
│  │        IMPUTE                │  │
│  │  ┌─────────┐ ┌──────────┐   │  │
│  │  │MINIMAC4 │ │ IMPUTE5  │   │  │
│  │  └─────────┘ └──────────┘   │  │
│  └──────────────────────────────┘  │
│  ┌──────────────────────────────┐  │
│  │        REPORT                │  │
│  │  30+ visualization modules   │  │
│  └──────────────────────────────┘  │
└─────────────────────────────────────┘
```

**Workflow Components**:
- 4 main subworkflows
- 50+ individual processes
- Complex branching
- Dynamic execution paths

### Workflow Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Subworkflows** | 7 modules | 4 major + many minor |
| **Total Processes** | ~15-20 | 50+ |
| **Execution Pattern** | Linear | DAG with branches |
| **Conditional Logic** | Simple flags | Complex conditions |
| **Data Flow** | Sequential | Parallel channels |

---

## 3. Module Organization

### ImputationServer2

**Module Structure**:
```
modules/
└── local/
    ├── ancestry_estimation/
    │   └── ancestry_estimation.nf
    ├── compression/
    │   └── compress_vcf.nf
    ├── imputation/
    │   ├── impute_minimac5.nf
    │   └── merge_imputed_vcf.nf
    ├── input_validation/
    │   └── validate_input.nf
    ├── pgs_calculation/
    │   └── calculate_pgs.nf
    ├── phasing/
    │   ├── phase_eagle.nf
    │   └── phase_beagle.nf
    └── quality_control/
        ├── quality_control_report.nf
        └── quality_control_vcf.nf
```

**Module Characteristics**:
- Simple, flat structure
- Single-purpose modules
- Minimal interdependencies
- Standard naming conventions

### H3ABioNet Pipeline

**Module Structure**:
```
modules/
├── local/
│   ├── qc/
│   │   ├── check_files.nf
│   │   ├── check_chromosome.nf
│   │   ├── check_genome_build.nf
│   │   ├── check_mismatch.nf
│   │   ├── check_global_mismatch.nf
│   │   ├── check_overlap.nf
│   │   ├── filter_min_ac.nf
│   │   ├── generate_chunk_map.nf
│   │   ├── merge_adjacent_chunks.nf
│   │   ├── qc_dupl.nf
│   │   ├── qc_site_missingness.nf
│   │   └── split_multi_allelic.nf
│   ├── phasing/
│   │   ├── eagle_phasing.nf
│   │   ├── shapeit4_phasing.nf
│   │   └── beagle5_phasing.nf
│   ├── impute/
│   │   ├── impute_minimac4.nf
│   │   ├── impute_minimac4_safe.nf
│   │   ├── impute_impute5.nf
│   │   └── validate_imputation_chunk.nf
│   ├── report/
│   │   ├── plot_performance.nf
│   │   ├── plot_r2_maf.nf
│   │   ├── plot_accuracy.nf
│   │   ├── generate_chunk_json.nf
│   │   ├── plot_genome_summary.nf
│   │   └── [25+ more plotting modules]
│   └── subset_vcf/
│       ├── generate_chunks_vcf.nf
│       └── split_target_to_chunk.nf
└── nf-core/  # Optional nf-core modules
```

**Module Characteristics**:
- Deep hierarchical structure
- Multi-purpose modules
- Complex interdependencies
- Domain-specific organization

### Module Comparison

| Aspect | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Module Count** | ~10 | 50+ |
| **Structure Depth** | 2 levels | 3-4 levels |
| **Modularity** | Basic | Advanced |
| **Reusability** | Limited | High |
| **Testing** | Simple | Comprehensive |

---

## 4. Quality Control Implementation

### ImputationServer2

**QC Approach**:
```nextflow
// Simple, fixed QC
process QUALITY_CONTROL_VCF {
    input:
    tuple val(meta), path(vcf)
    
    script:
    """
    # Fixed thresholds
    bcftools +fill-tags ${vcf} -- -t AF,AC,AN |
    bcftools view \\
        --min-ac 1 \\
        --max-alleles 2 \\
        --min-alleles 2 \\
        -i 'F_MISSING<0.05' \\
        -O z -o ${meta.id}_qc.vcf.gz
    """
}
```

**QC Parameters**:
- Fixed thresholds
- Standard filters only
- Binary pass/fail
- No population adjustment

### H3ABioNet Pipeline

**QC Approach**:
```nextflow
// Adaptive, multi-stage QC
process CHECK_MISMATCH {
    input:
    tuple val(meta), path(vcf), path(ref)
    
    script:
    def max_mismatch = params.max_mismatch_rate ?: 0.2
    def min_variants = params.min_matching_variants ?: 20
    """
    # Adaptive thresholds
    if [[ "${meta.population}" == "African" ]]; then
        max_mismatch=0.25
        min_variants=15
    fi
    
    # Multi-tier assessment
    check_mismatch.py \\
        --vcf ${vcf} \\
        --ref ${ref} \\
        --max-mismatch \${max_mismatch} \\
        --min-variants \${min_variants} \\
        --output-status ${meta.id}.status
    
    # Handle warnings
    status=\$(cat ${meta.id}.status)
    if [[ "\${status}" == "WARN" ]]; then
        echo "Warning: Marginal quality for ${meta.id}"
    fi
    """
}

// Global QC check
process CHECK_GLOBAL_MISMATCH {
    input:
    val all_statuses
    
    script:
    """
    # Aggregate QC across all chunks
    fail_rate=\$(calculate_failure_rate.py ${all_statuses})
    if (( \$(echo "\${fail_rate} > 0.5" | bc -l) )); then
        echo "ERROR: Too many chunks failing QC"
        exit 1
    fi
    """
}
```

**QC Features**:
- Configurable thresholds
- Population-specific adjustments
- Three-tier system (PASS/WARN/FAIL)
- Global QC assessment
- Recovery mechanisms

### QC Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **QC Stages** | 2 | 12+ |
| **Threshold Type** | Fixed | Adaptive |
| **Population-Aware** | No | Yes |
| **Warning System** | No | Yes |
| **Recovery** | None | Chunk merging |
| **Global Checks** | No | Yes |

---

## 5. Phasing Approaches

### ImputationServer2

**Phasing Implementation**:
```nextflow
// Fixed phasing with limited options
process PHASE_EAGLE {
    label 'process_high'
    
    input:
    tuple val(meta), path(vcf)
    path reference
    path genetic_map
    
    script:
    """
    eagle \\
        --vcfTarget=${vcf} \\
        --vcfRef=${reference} \\
        --geneticMapFile=${genetic_map} \\
        --outPrefix=${meta.id} \\
        --numThreads=${task.cpus}
    """
}

// Alternative: Beagle only
process PHASE_BEAGLE {
    input:
    tuple val(meta), path(vcf)
    
    script:
    """
    beagle \\
        gt=${vcf} \\
        out=${meta.id}.phased \\
        nthreads=${task.cpus}
    """
}
```

**Phasing Options**:
- Eagle or Beagle
- No parameter tuning
- Fixed genetic maps
- Standard execution

### H3ABioNet Pipeline

**Phasing Implementation**:
```nextflow
// Flexible phasing with multiple methods
process EAGLE_PHASING {
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path genetic_map
    path ref_panel
    path ref_panel_index
    
    script:
    def args = task.ext.args ?: ''
    def chrm = meta.contig
    
    // Adaptive parameters
    def kpbwt = meta.sample_size > 5000 ? 50000 : 20000
    def pbwt_iters = params.phasing_iterations ?: 2
    def expect_ibd = params.expect_ibd_cm ?: 3.0
    
    // Population-specific settings
    def pop_args = ""
    if (meta.population == "African") {
        pop_args = "--expectIBDcM 5.0 --histFactor 1.0"
    }
    
    """
    eagle \\
        --vcfTarget=${vcf} \\
        --vcfRef=${ref_panel} \\
        --geneticMapFile=${genetic_map} \\
        --outPrefix=${meta.id}.phased \\
        --chrom=${chrm}:${meta.start}-${meta.end} \\
        --Kpbwt=${kpbwt} \\
        --pbwtIters=${pbwt_iters} \\
        --numThreads=${task.cpus} \\
        ${pop_args} \\
        ${args}
    """
}

// Alternative methods
process SHAPEIT4_PHASING {
    input:
    tuple val(meta), path(vcf)
    path reference
    path genetic_map
    
    script:
    def window = params.shapeit4_window ?: 2
    def pbwt_depth = params.pbwt_depth ?: 8
    """
    shapeit4 \\
        --input ${vcf} \\
        --reference ${reference} \\
        --map ${genetic_map} \\
        --region ${meta.contig}:${meta.start}-${meta.end} \\
        --window ${window} \\
        --pbwt-depth ${pbwt_depth} \\
        --thread ${task.cpus} \\
        --output ${meta.id}.phased.vcf.gz
    """
}
```

### Phasing Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Methods** | Eagle, Beagle | Eagle, SHAPEIT4, Beagle5 |
| **Parameter Control** | None | Full |
| **Population-Specific** | No | Yes |
| **Adaptive Settings** | No | Yes |
| **Custom Maps** | Limited | Full support |

---

## 6. Imputation Implementation

### ImputationServer2

**Imputation Process**:
```nextflow
process IMPUTE_MINIMAC5 {
    label 'process_medium'
    
    input:
    tuple val(meta), path(phased_vcf)
    path reference
    
    script:
    def window = params.imputation_window ?: 500000
    """
    minimac5 \\
        --haps ${phased_vcf} \\
        --refHaps ${reference} \\
        --prefix ${meta.id} \\
        --window ${window} \\
        --format GT,DS \\
        --threads ${task.cpus}
    """
}
```

**Features**:
- Minimac5 only
- Fixed parameters
- Standard output formats
- No validation

### H3ABioNet Pipeline

**Imputation Process**:
```nextflow
process IMPUTE_MINIMAC4_SAFE {
    input:
    tuple val(meta), path(phased_vcf)
    tuple val(ref_meta), val(ref_name), path(ref_m3vcf)
    
    script:
    def args = task.ext.args ?: ''
    def min_ratio = params.minRatio ?: 0.01
    def format = params.output_format ?: 'GT,DS'
    def window = params.impute_window ?: 500000
    
    """
    # Pre-validation
    validate_imputation_chunk.py \\
        --vcf ${phased_vcf} \\
        --ref ${ref_m3vcf} \\
        --min-ratio ${min_ratio} \\
        --region ${meta.contig}:${meta.start}-${meta.end} || {
        
        echo "Chunk ${meta.id} has insufficient variants"
        create_empty_output.sh ${meta.id}_${ref_name}
        exit 0
    }
    
    # Safe imputation with error handling
    minimac4 \\
        --refHaps ${ref_m3vcf} \\
        --haps ${phased_vcf} \\
        --prefix ${meta.id}_${ref_name} \\
        --format ${format} \\
        --minRatio ${min_ratio} \\
        --window ${window} \\
        --cpus ${task.cpus} \\
        ${args} || {
        
        echo "Imputation failed for ${meta.id}"
        # Attempt recovery
        if [[ -f "${phased_vcf}.backup" ]]; then
            retry_with_relaxed_params.sh
        else
            create_empty_output.sh ${meta.id}_${ref_name}
        fi
    }
    
    # Post-validation
    validate_output.py \\
        --vcf ${meta.id}_${ref_name}.dose.vcf.gz \\
        --min-variants 100
    """
}

// Alternative imputation method
process IMPUTE_IMPUTE5 {
    input:
    tuple val(meta), path(phased_vcf)
    path reference
    path genetic_map
    
    script:
    """
    impute5 \\
        --h ${reference} \\
        --m ${genetic_map} \\
        --g ${phased_vcf} \\
        --r ${meta.contig}:${meta.start}-${meta.end} \\
        --o ${meta.id}.imputed.vcf.gz \\
        --threads ${task.cpus}
    """
}
```

### Imputation Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Tools** | Minimac5 | Minimac4, IMPUTE5 |
| **Validation** | None | Pre & Post |
| **Error Recovery** | None | Multiple strategies |
| **Safe Mode** | No | Yes |
| **Multi-Panel** | No | Yes |

---

## 7. Configuration Management

### ImputationServer2

**Configuration Structure**:
```groovy
// nextflow.config - Simple, fixed
params {
    // Fixed parameters
    project = 'my-project'
    build = 'hg38'
    mode = 'imputation'
    
    // Limited options
    chunk_size = 20000000
    min_samples = 20
    max_samples = 50000
    
    // Phasing
    phasing_engine = 'eagle'  // or 'beagle'
    
    // Output
    merge_results = true
    encrypt_results = false
}

// Profiles
profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }
    
    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
    }
    
    slurm {
        process.executor = 'slurm'
    }
}
```

### H3ABioNet Pipeline

**Configuration Structure**:
```groovy
// nextflow.config - Extensive, flexible
params {
    // Input/Output
    input = null
    outdir = './results'
    
    // Processing options
    chromosomes = 'ALL'
    chunk_size = 25000000  // Adaptive
    buffer_size = 500000
    
    // QC Parameters (all configurable)
    site_miss = 0.05
    hwe = 1e-5
    mac = 1
    maf_thresh = 0.1
    r2_threshold = 0.3
    
    // Mismatch handling
    max_mismatch_rate = 0.2
    min_matching_variants = 20
    warn_mismatch_rate = 0.3
    global_failure_threshold = 0.5
    
    // Phasing options
    phasing_method = 'eagle'  // eagle, shapeit4, beagle5
    phasing_iterations = 2
    expect_ibd_cm = 3.0
    shapeit4_window = 2
    pbwt_depth = 8
    
    // Imputation options
    impute_method = 'minimac4'  // minimac4, impute5
    minRatio = 0.01
    output_format = 'GT,DS'
    impute_window = 500000
    use_safe_imputation = true
    validate_chunks = false
    
    // Reference panels (multiple)
    ref_panels = [
        ["1000G", "1kg_p3", "/refs/1000g/chr%s.bcf"],
        ["H3Africa", "h3a_v2", "/refs/h3a/chr%s.bcf"]
    ]
    
    // Population-specific
    population = null  // Triggers adaptive settings
    
    // Reporting
    generate_all_plots = true
    hierarchical_reports = true
    plot_formats = ['pdf', 'png']
    
    // Performance
    max_cpus = 16
    max_memory = '128.GB'
    max_time = '240.h'
    
    // Advanced options
    enable_gpu = false
    cache_dir = './work/cache'
    cleanup = false
}

// Complex profiles
profiles {
    standard {
        includeConfig 'conf/base.config'
    }
    
    african {
        params.max_mismatch_rate = 0.25
        params.chunk_size = 30000000
        params.expect_ibd_cm = 5.0
    }
    
    clinical {
        params.max_mismatch_rate = 0.10
        params.validate_chunks = true
        params.use_safe_imputation = true
    }
    
    test {
        includeConfig 'conf/test.config'
        params.max_cpus = 2
        params.max_memory = '8.GB'
    }
    
    slurm {
        includeConfig 'conf/slurm.config'
        executor.queueSize = 500
        executor.submitRateLimit = '200/1min'
    }
}
```

### Configuration Comparison

| Aspect | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Parameters** | ~15 | 50+ |
| **Profiles** | 3 basic | 10+ specialized |
| **Population-Specific** | No | Yes |
| **Resource Control** | Basic | Advanced |
| **Customization** | Minimal | Extensive |

---

## 8. Resource Management

### ImputationServer2

**Resource Allocation**:
```nextflow
// Simple, fixed resources
process {
    withLabel: process_low {
        cpus = 2
        memory = '4.GB'
        time = '2.h'
    }
    
    withLabel: process_medium {
        cpus = 4
        memory = '8.GB'
        time = '4.h'
    }
    
    withLabel: process_high {
        cpus = 8
        memory = '16.GB'
        time = '8.h'
    }
}
```

### H3ABioNet Pipeline

**Resource Allocation**:
```nextflow
// Adaptive, fine-grained resources
process {
    // Dynamic allocation based on input
    withName: EAGLE_PHASING {
        cpus = { meta.sample_size > 5000 ? 16 : 8 }
        memory = { meta.sample_size > 5000 ? '64.GB' : '32.GB' }
        time = { meta.chunk_size > 10000000 ? '8.h' : '4.h' }
        maxRetries = 3
        errorStrategy = { task.exitStatus in [104, 134, 139] ? 'retry' : 'finish' }
    }
    
    withName: IMPUTE_MINIMAC4 {
        cpus = 4
        memory = { 16.GB * task.attempt }
        time = { 4.h * task.attempt }
        maxForks = 200  // Parallel execution limit
    }
    
    withName: 'CHECK_MISMATCH|CHECK_OVERLAP' {
        cpus = 1
        memory = '2.GB'
        maxForks = 50
    }
    
    // Population-specific resources
    withName: 'PLOT_.*' {
        cpus = 2
        memory = '8.GB'
        container = 'mamana/python-plotting:1.1.0'
    }
}

// SLURM-specific optimizations
executor {
    $slurm {
        queueSize = 500
        submitRateLimit = '200/1min'
        exitReadTimeout = '120 min'
        killBatchSize = 50
    }
}
```

### Resource Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Allocation Type** | Static | Dynamic |
| **Retry Logic** | Basic | Advanced |
| **Parallelization** | Limited | Extensive |
| **Fork Control** | No | Yes |
| **Queue Management** | Simple | Optimized |

---

## 9. Container Strategy

### ImputationServer2

**Single Container Approach**:
```dockerfile
# Single monolithic container
FROM quay.io/genepi/imputationserver2:v2.0.8
# Contains all tools: eagle, beagle, minimac5, bcftools
```

**Container Usage**:
```nextflow
process.container = 'quay.io/genepi/imputationserver2:v2.0.8'
```

### H3ABioNet Pipeline

**Multi-Container Strategy**:
```nextflow
// Process-specific containers
process {
    withName: 'EAGLE_PHASING|SHAPEIT4_PHASING' {
        container = 'quay.io/h3abionet_org/py3plink'
    }
    
    withName: 'IMPUTE_MINIMAC4.*' {
        container = 'mamana/minimac4:latest'
    }
    
    withName: 'CHECK_.*|QC_.*' {
        container = 'quay.io/biocontainers/bcftools:1.11'
    }
    
    withName: 'PLOT_.*|REPORT_.*' {
        container = 'mamana/python-plotting:1.1.0'
    }
    
    withName: 'CONVERT_.*' {
        container = 'mamana/vcf-processing:bcftools-1.20'
    }
}
```

**Container Features**:
- Specialized containers per task
- Smaller image sizes
- Version pinning
- Local registry support

### Container Comparison

| Aspect | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Strategy** | Monolithic | Microservices |
| **Container Count** | 1 | 5+ |
| **Image Size** | Large (~2GB) | Small (~200-500MB each) |
| **Flexibility** | Low | High |
| **Maintenance** | Simple | Complex |

---

## 10. Error Handling & Recovery

### ImputationServer2

**Basic Error Handling**:
```nextflow
process IMPUTE_MINIMAC5 {
    errorStrategy 'terminate'
    
    script:
    """
    minimac5 ... || exit 1
    """
}
```

### H3ABioNet Pipeline

**Advanced Error Recovery**:
```nextflow
process IMPUTE_MINIMAC4_SAFE {
    maxRetries 3
    errorStrategy { task.exitStatus in [104, 134, 139, 140, 143, 137] ? 'retry' : 
                   task.exitStatus == 255 ? 'ignore' : 'finish' }
    
    script:
    """
    # Multiple recovery strategies
    set -euo pipefail
    
    # Pre-check
    if ! validate_chunk.py ${phased_vcf}; then
        echo "Invalid chunk - attempting recovery"
        fix_chunk.py ${phased_vcf} > fixed.vcf
        phased_vcf=fixed.vcf
    fi
    
    # Main process with fallback
    minimac4 ... || {
        exit_code=\$?
        case \$exit_code in
            137|143)  # Memory issues
                echo "Memory error - will retry with more memory"
                exit 137
                ;;
            255)  # Empty chunk
                echo "Empty chunk - creating placeholder"
                create_empty_output.sh
                exit 0
                ;;
            *)
                echo "Unknown error: \$exit_code"
                exit \$exit_code
                ;;
        esac
    }
    
    # Post-check
    if [[ ! -f output.vcf.gz ]] || [[ \$(bcftools view -H output.vcf.gz | wc -l) -eq 0 ]]; then
        echo "Output validation failed"
        exit 1
    fi
    """
}
```

### Error Handling Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Retry Logic** | Basic | Advanced |
| **Error Codes** | Not handled | Specific handling |
| **Recovery** | None | Multiple strategies |
| **Validation** | Post-only | Pre & Post |
| **Fallbacks** | None | Empty outputs, fixes |

---

## 11. Reporting & Outputs

### ImputationServer2

**Basic Reporting**:
```nextflow
process QUALITY_CONTROL_REPORT {
    publishDir "${params.output}/reports"
    
    input:
    path qc_stats
    
    output:
    path "qc_report.html"
    
    script:
    """
    generate_qc_report.py ${qc_stats} > qc_report.html
    """
}
```

**Outputs**:
- Basic QC report
- Imputed VCFs
- Info scores
- Simple statistics

### H3ABioNet Pipeline

**Comprehensive Reporting**:
```nextflow
// Hierarchical reporting system
workflow REPORT {
    take:
    ch_imputed
    
    main:
    // Chunk-level reports (22 types)
    PLOT_PERFORMANCE(ch_imputed)
    PLOT_R2_MAF(ch_imputed)
    PLOT_ACCURACY(ch_imputed)
    GENERATE_CHUNK_JSON(ch_imputed)
    
    // Chromosome-level aggregation
    COMBINE_FREQ_BY_CHR(EXTRACT_FREQ.out)
    PLOT_CHR_PERFORMANCE(COMBINE_FREQ_BY_CHR.out)
    
    // Genome-level summaries
    COMBINE_FREQ_GENOME(COMBINE_FREQ_BY_CHR.out)
    PLOT_GENOME_SUMMARY(COMBINE_FREQ_GENOME.out)
    PLOT_IMPUTATION_ACCURACY_MAF_BINS(COMBINE_FREQ_GENOME.out)
    PLOT_AGGREGATED_R2_DASHBOARD(COMBINE_FREQ_GENOME.out)
    
    // Generate final report
    GENERATE_SUMMARY_REPORT(
        PLOT_PERFORMANCE.out.plots.collect(),
        GENERATE_CHUNK_JSON.out.json.collect()
    )
}
```

**Output Structure**:
```
results/
├── imputed/
│   ├── chunks/
│   ├── chromosomes/
│   └── merged/
├── qc/
│   ├── preprocessing/
│   ├── mismatch/
│   └── global/
├── plots/
│   ├── chunk_level/
│   ├── chromosome_level/
│   └── genome_level/
├── reports/
│   ├── html/
│   ├── pdf/
│   └── json/
└── logs/
```

### Reporting Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Report Types** | 2-3 | 30+ |
| **Visualization** | Basic | Comprehensive |
| **Hierarchical** | No | Yes (3 levels) |
| **Formats** | HTML | HTML, PDF, JSON |
| **Interactive** | No | Dashboard options |

---

## 12. Extended Features

### ImputationServer2

**Additional Features**:
- Ancestry estimation
- PGS calculation
- Result encryption
- Email notifications

```nextflow
// Optional ancestry estimation
workflow ANCESTRY_ESTIMATION {
    take:
    ch_vcf
    
    main:
    ESTIMATE_ANCESTRY(ch_vcf)
    
    emit:
    ancestry = ESTIMATE_ANCESTRY.out
}

// Optional PGS
workflow PGS_CALCULATION {
    take:
    ch_imputed
    
    main:
    CALCULATE_PGS(ch_imputed, params.pgs_catalog)
    
    emit:
    scores = CALCULATE_PGS.out
}
```

### H3ABioNet Pipeline

**Extended Capabilities**:
- Multi-panel comparison
- Population stratification
- Method benchmarking
- Safe imputation mode
- Chunk recovery
- Custom reporting

```nextflow
// Multi-panel comparison
workflow COMPARE_PANELS {
    take:
    ch_phased
    
    main:
    // Run imputation with multiple panels
    ch_panels = Channel.fromList(params.ref_panels)
    
    ch_imputed = ch_panels.combine(ch_phased)
        .map { panel_info, phased ->
            def (name, version, path) = panel_info
            [phased, name, file(path)]
        }
        .set { ch_multi_panel }
    
    IMPUTE_MULTI(ch_multi_panel)
    
    // Compare results
    COMPARE_IMPUTATION_QUALITY(IMPUTE_MULTI.out.groupTuple())
    PLOT_PANEL_COMPARISON(COMPARE_IMPUTATION_QUALITY.out)
    
    emit:
    comparison = PLOT_PANEL_COMPARISON.out
}
```

### Extended Features Comparison

| Feature | ImputationServer2 | H3ABioNet |
|---------|------------------|-----------|
| **Ancestry** | Yes | Via external tools |
| **PGS** | Yes | Via external tools |
| **Multi-Panel** | No | Yes |
| **Benchmarking** | No | Yes |
| **Safe Mode** | No | Yes |
| **Custom Analysis** | Limited | Extensive |

---

## 13. Development Practices

### ImputationServer2

**Development Approach**:
- Version tagged releases
- Basic testing
- Docker-centric
- Service-oriented updates

```bash
# Testing approach
nextflow run main.nf -profile test,docker
```

### H3ABioNet Pipeline

**Development Approach**:
- Git flow branching
- Comprehensive testing
- CI/CD integration
- Module unit tests
- Integration tests

```bash
# Comprehensive testing
# Unit tests
pytest tests/modules/

# Integration tests
nextflow run main.nf -profile test,docker --test_data tiny
nextflow run main.nf -profile test,singularity --test_data small
nextflow run main.nf -profile test,slurm --test_data full

# Continuous integration
.github/workflows/ci.yml
```

### Development Comparison

| Aspect | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Testing** | Basic | Comprehensive |
| **CI/CD** | Limited | Full pipeline |
| **Documentation** | User-focused | Developer + User |
| **Modularity** | Low | High |
| **Contribution Model** | Closed | Open |

---

## 14. Performance Optimization

### ImputationServer2

**Performance Strategy**:
- Fixed chunking (20Mb)
- Standard parallelization
- Simple resource allocation

```nextflow
// Fixed performance settings
params.chunk_size = 20000000
process.cpus = 8
process.memory = '16.GB'
```

### H3ABioNet Pipeline

**Performance Strategy**:
- Adaptive chunking
- Dynamic parallelization
- Resource optimization
- Cache utilization

```nextflow
// Adaptive performance
params {
    // Dynamic chunking based on density
    chunk_size = { vcf ->
        def density = calculate_variant_density(vcf)
        density > 1000 ? 15000000 : 25000000
    }
}

// Parallelization optimization
process {
    withName: IMPUTE_MINIMAC4 {
        maxForks = { 
            def available_slots = executor.queueSize - executor.queued
            Math.min(200, available_slots)
        }
    }
}

// Cache strategy
process {
    cache = 'deep'  // Content-based caching
    storeDir = { "${params.cache_dir}/${task.process}/${task.hash}" }
}
```

### Performance Comparison

| Metric | ImputationServer2 | H3ABioNet |
|--------|------------------|-----------|
| **Chunking** | Fixed 20Mb | Adaptive 15-30Mb |
| **Parallelization** | Static | Dynamic |
| **Caching** | Basic | Advanced |
| **Resource Usage** | Fixed | Optimized |
| **Scalability** | Limited | Excellent |

---

## 15. Use Case Comparison

### Best Use Cases for ImputationServer2 Nextflow

1. **Standardized Production Pipelines**
   - Need consistent, reproducible results
   - Standard populations and panels
   - Minimal customization required

2. **Quick Deployment Scenarios**
   - Rapid setup needed
   - Limited bioinformatics expertise
   - Standard QC requirements

3. **Service Deployment**
   - Building imputation services
   - Multi-user environments
   - Simplified maintenance

**Example Configuration**:
```nextflow
// ImputationServer2 for production service
params {
    mode = 'imputation'
    build = 'hg38'
    phasing_engine = 'eagle'
    merge_results = true
    encrypt_results = true
}
```

### Best Use Cases for H3ABioNet Pipeline

1. **Research & Development**
   - Method comparison needed
   - Novel populations
   - Custom reference panels
   - Extensive QC requirements

2. **Population Genomics**
   - African populations
   - Admixed populations
   - Rare variant analysis
   - Population-specific optimization

3. **Clinical & Precision Medicine**
   - Strict QC requirements
   - Audit trails needed
   - Custom reporting
   - Safe mode operation

**Example Configuration**:
```nextflow
// H3ABioNet for African genomics
params {
    population = 'African'
    max_mismatch_rate = 0.25
    phasing_method = 'shapeit4'
    use_safe_imputation = true
    ref_panels = [
        ["H3Africa", "h3a_v2", "/refs/h3africa/chr%s.bcf"],
        ["1000G_AFR", "1kg_afr", "/refs/1000g_afr/chr%s.bcf"]
    ]
    generate_all_plots = true
}
```

---

## Summary Matrix

| Category | ImputationServer2 | H3ABioNet | Winner |
|----------|------------------|-----------|---------|
| **Simplicity** | ⭐⭐⭐⭐⭐ | ⭐⭐ | ImputationServer2 |
| **Flexibility** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **QC Depth** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Performance** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Population Support** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Error Recovery** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Reporting** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Maintenance** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ImputationServer2 |
| **Scalability** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |
| **Research Use** | ⭐⭐ | ⭐⭐⭐⭐⭐ | H3ABioNet |

---

## Conclusions

### ImputationServer2 Nextflow Pipeline
**Strengths**:
- Simple deployment and maintenance
- Standardized workflow
- Good for production services
- Minimal configuration needed
- Stable and tested

**Weaknesses**:
- Limited flexibility
- Basic QC and reporting
- No population-specific features
- Limited error recovery
- Fixed resource allocation

### H3ABioNet Pipeline
**Strengths**:
- Extreme flexibility and customization
- Comprehensive QC and reporting
- Population-specific optimizations
- Advanced error recovery
- Superior scalability

**Weaknesses**:
- Complex configuration
- Steeper learning curve
- Higher maintenance burden
- Requires more expertise

### Final Recommendation

**Choose ImputationServer2 Nextflow when**:
- You need a simple, standardized pipeline
- Quick deployment is priority
- Limited customization is acceptable
- Building a service platform

**Choose H3ABioNet Pipeline when**:
- Research flexibility is crucial
- Working with diverse populations
- Need comprehensive QC and reporting
- Require method comparison capabilities
- Building production genomics pipelines

Both pipelines represent valid approaches to genotype imputation, with ImputationServer2 prioritizing simplicity and standardization, while H3ABioNet emphasizes flexibility and comprehensive analysis capabilities. The choice depends on specific use case requirements and available expertise.

---

*Document Version: 1.0*  
*Last Updated: 2025*  
*Nextflow Pipeline Implementation Comparison*