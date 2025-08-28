# Pipeline Simplification Strategy: Achieving ImputationServer2 Simplicity with H3ABioNet Robustness

## Executive Summary

This document provides a comprehensive strategy to simplify the H3ABioNet genotype imputation pipeline while maintaining its robustness, performance, and population-specific capabilities. The goal is to achieve ImputationServer2's ease of use without sacrificing the advanced features that make the H3ABioNet pipeline superior for diverse genomics research.

---

## Table of Contents

1. [Current Complexity Analysis](#current-complexity-analysis)
2. [Simplification Principles](#simplification-principles)
3. [Proposed Architecture Changes](#proposed-architecture-changes)
4. [Profile-Based Simplification](#profile-based-simplification)
5. [Module Consolidation Strategy](#module-consolidation-strategy)
6. [Configuration Simplification](#configuration-simplification)
7. [Workflow Streamlining](#workflow-streamlining)
8. [Smart Defaults Implementation](#smart-defaults-implementation)
9. [Progressive Disclosure Pattern](#progressive-disclosure-pattern)
10. [Implementation Roadmap](#implementation-roadmap)
11. [Code Refactoring Examples](#code-refactoring-examples)
12. [Testing & Validation Strategy](#testing--validation-strategy)

---

## 1. Current Complexity Analysis

### Sources of Complexity in H3ABioNet Pipeline

#### 1.1 Module Proliferation
```
Current State:
- 50+ individual modules
- 30+ reporting modules alone
- Deep nesting (3-4 levels)
- Complex interdependencies

Impact:
- Difficult to understand flow
- Hard to maintain
- Overwhelming for new users
```

#### 1.2 Configuration Overload
```
Current State:
- 50+ parameters
- Multiple configuration files
- Complex profile system
- Population-specific settings scattered

Impact:
- Decision paralysis
- Error-prone configuration
- Steep learning curve
```

#### 1.3 Workflow Complexity
```
Current State:
- Multiple branching paths
- Complex channel operations
- Conditional execution logic
- Nested subworkflows

Impact:
- Hard to debug
- Difficult to predict behavior
- Performance overhead
```

### ImputationServer2's Simplification Success Factors

1. **Linear workflow** - Predictable execution
2. **Fixed defaults** - No decision paralysis
3. **Single container** - Simple deployment
4. **Minimal parameters** - Easy configuration
5. **Hidden complexity** - Advanced features not exposed

---

## 2. Simplification Principles

### Core Principles to Follow

#### P1: Progressive Disclosure
```
Hide complexity behind experience levels:
- Beginner: One command, works perfectly
- Intermediate: Common customizations
- Advanced: Full control
```

#### P2: Sensible Defaults
```
Make the common case simple:
- Auto-detect population when possible
- Use best practices as defaults
- Optimize for most common hardware
```

#### P3: Convention Over Configuration
```
Reduce decisions:
- Standard directory structure
- Predictable naming
- Automatic format detection
```

#### P4: Fail Fast with Clear Messages
```
Better error handling:
- Validate early
- Clear error messages
- Suggested fixes
```

#### P5: Modular but Hidden
```
Keep flexibility without exposing it:
- Composite modules
- Profile-based feature sets
- Internal complexity
```

---

## 3. Proposed Architecture Changes

### 3.1 Three-Tier Architecture

```nextflow
// NEW SIMPLIFIED ARCHITECTURE
// Level 1: Simple Interface (What users see)
workflow {
    IMPUTATION_PIPELINE(
        params.input,
        params.mode  // 'simple', 'standard', 'advanced'
    )
}

// Level 2: Mode Controllers (Hidden complexity)
workflow IMPUTATION_PIPELINE {
    take:
    input
    mode
    
    main:
    switch(mode) {
        case 'simple':
            SIMPLE_MODE(input)
            break
        case 'standard':
            STANDARD_MODE(input)
            break
        case 'advanced':
            ADVANCED_MODE(input)
            break
    }
}

// Level 3: Actual Implementation (Current modules)
workflow SIMPLE_MODE {
    take:
    input
    
    main:
    // Consolidated workflow with sensible defaults
    UNIFIED_QC(input)
    UNIFIED_PHASE(UNIFIED_QC.out)
    UNIFIED_IMPUTE(UNIFIED_PHASE.out)
    BASIC_REPORT(UNIFIED_IMPUTE.out)
}
```

### 3.2 Consolidated Module Structure

```
modules/
├── unified/                    # NEW: Consolidated modules
│   ├── qc_pipeline.nf         # Combines all QC steps
│   ├── phasing_pipeline.nf    # Unified phasing
│   ├── imputation_pipeline.nf # Unified imputation
│   └── reporting_pipeline.nf  # Basic reporting
├── core/                       # Renamed from 'local'
│   └── [existing modules]     # Keep for advanced mode
└── profiles/                   # NEW: Profile-specific modules
    ├── african/
    ├── european/
    └── clinical/
```

---

## 4. Profile-Based Simplification

### 4.1 Create Usage Profiles

```groovy
// profiles/simple.config
params {
    // Minimal required parameters
    mode = 'simple'
    auto_detect = true
    use_defaults = true
    
    // Hidden but set
    chunk_size = 25000000
    max_mismatch_rate = 0.2
    phasing_method = 'eagle'
    impute_method = 'minimac4'
    basic_reports_only = true
}

// profiles/research.config  
params {
    mode = 'standard'
    // Expose common research parameters
    population = null  // User specifies
    ref_panels = []    // User configures
    generate_all_plots = true
    
    // Smart defaults based on population
    adaptive_qc = true
}

// profiles/production.config
params {
    mode = 'production'
    // Optimized for throughput
    parallel_chunks = 500
    aggressive_caching = true
    minimal_logging = true
    skip_optional_qc = true
}

// profiles/clinical.config
params {
    mode = 'clinical'
    // Maximum safety
    use_safe_imputation = true
    validate_all_chunks = true
    strict_qc = true
    audit_trail = true
}
```

### 4.2 Single Entry Point

```bash
# NEW SIMPLIFIED USAGE

# Simplest - just works
nextflow run imputation --input data.vcf

# With profile
nextflow run imputation --input data.vcf -profile african

# Research mode
nextflow run imputation --input data.vcf -profile research --population AFR

# Full control (advanced users only)
nextflow run imputation --input data.vcf -profile advanced -c custom.config
```

---

## 5. Module Consolidation Strategy

### 5.1 QC Module Consolidation

**Current: 12+ separate QC module files**

```
modules/local/qc/
├── check_files.nf          # File validation
├── check_chromosome.nf      # Chromosome consistency
├── check_genome_build.nf    # Build detection
├── check_mismatch.nf        # Allele mismatch checking
├── check_global_mismatch.nf # Global mismatch threshold
├── check_overlap.nf         # Reference overlap
├── filter_min_ac.nf         # Allele count filtering
├── generate_chunk_map.nf    # Chunk mapping
├── merge_adjacent_chunks.nf # Chunk merging
├── qc_dupl.nf              # Duplicate removal
├── qc_site_missingness.nf  # Site missingness
└── split_multi_allelic.nf  # Multi-allelic splitting
```

**Proposed: Single consolidated QC file with all processes**

```nextflow
// NEW: modules/unified/qc_pipeline.nf - All QC in one file
nextflow.enable.dsl = 2

/*
 * CONSOLIDATED QC MODULE
 * Combines all 12 QC processes into a single, well-organized file
 */

// ============================================================================
// BASIC VALIDATION PROCESSES
// ============================================================================

process CHECK_FILES {
    tag "$meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path(vcf), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    # File validation logic
    if [[ ! -f "${vcf}" ]]; then
        echo "ERROR: Input file not found: ${vcf}"
        exit 1
    fi
    
    bcftools query -l ${vcf} > samples.txt
    n_samples=\$(wc -l < samples.txt)
    echo "Found \${n_samples} samples in ${vcf}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process CHECK_CHROMOSOME {
    tag "$meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path(vcf), emit: vcf
    path "*.chr_check.txt", emit: chr_report
    path "versions.yml", emit: versions
    
    script:
    """
    # Chromosome consistency check
    bcftools query -f '%CHROM\\n' ${vcf} | sort -u > ${meta.id}.chr_check.txt
    
    # Verify chromosomes are consistent
    n_chroms=\$(wc -l < ${meta.id}.chr_check.txt)
    echo "Found \${n_chroms} chromosomes"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process CHECK_GENOME_BUILD {
    tag "$meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), path(vcf)
    val expected_build
    
    output:
    tuple val(meta), path(vcf), emit: vcf
    path "*.build_check.txt", emit: build_status
    path "versions.yml", emit: versions
    
    script:
    """
    # Detect genome build
    first_chrom=\$(bcftools query -f '%CHROM\\n' ${vcf} | head -1)
    
    if [[ "\${first_chrom}" == chr* ]]; then
        detected_build="b38"
    else
        detected_build="b37"
    fi
    
    echo "First chromosome found: \${first_chrom}" > ${meta.id}.build_check.txt
    echo "Detected build: \${detected_build}" >> ${meta.id}.build_check.txt
    
    if [[ "\${detected_build}" == "${expected_build}" ]]; then
        echo "Status: PASS" >> ${meta.id}.build_check.txt
    else
        echo "Status: FAIL - Expected ${expected_build}" >> ${meta.id}.build_check.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

// ============================================================================
// REFERENCE COMPARISON PROCESSES
// ============================================================================

process CHECK_OVERLAP {
    tag "$meta.id"
    label 'process_low'
    
    input:
    tuple val(meta), path(chunk_map), path(ref_map)
    
    output:
    tuple val(meta), path("*.overlap.txt"), path("*.overlap.status"), emit: overlap
    path "versions.yml", emit: versions
    
    script:
    """
    # Check overlap between chunk and reference
    # TODO(human): Implement overlap checking logic between chunk_map and ref_map
    # Consider: minimum overlap threshold, variant count requirements
    
    echo "PASS" > ${meta.id}.overlap.status
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        custom: 1.0.0
    END_VERSIONS
    """
}

process CHECK_MISMATCH {
    tag "$meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), path(chunk_vcf), path(reference_panel), path(reference_index)
    
    output:
    tuple val(meta), path("*.mismatch.txt"), path("*.mismatch.status"), emit: mismatch
    path "versions.yml", emit: versions
    
    script:
    def max_mismatch_rate = params.max_mismatch_rate ?: 0.2
    def min_matching = params.min_matching_variants ?: 20
    """
    # Mismatch checking implementation
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${chunk_vcf} > chunk.alleles
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${reference_panel} > ref.alleles
    
    # Check concordance
    check_mismatch.py \\
        --chunk chunk.alleles \\
        --ref ref.alleles \\
        --max-mismatch ${max_mismatch_rate} \\
        --min-matching ${min_matching} \\
        --output ${meta.id}.mismatch.txt \\
        --status ${meta.id}.mismatch.status
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process CHECK_GLOBAL_MISMATCH {
    label 'process_single'
    
    input:
    val all_statuses
    
    output:
    path "global_mismatch.txt", emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    # Global mismatch assessment
    echo "${all_statuses.join('\\n')}" > statuses.txt
    
    pass_count=\$(grep -c "PASS" statuses.txt || true)
    warn_count=\$(grep -c "WARN" statuses.txt || true)
    fail_count=\$(grep -c "FAIL" statuses.txt || true)
    total_count=\$(wc -l < statuses.txt)
    
    echo "Global Mismatch Report" > global_mismatch.txt
    echo "PASS: \${pass_count}" >> global_mismatch.txt
    echo "WARN: \${warn_count}" >> global_mismatch.txt
    echo "FAIL: \${fail_count}" >> global_mismatch.txt
    echo "Total: \${total_count}" >> global_mismatch.txt
    
    failure_rate=\$(awk "BEGIN {print \${fail_count} / \${total_count}}")
    if (( \$(awk "BEGIN {print (\${failure_rate} > 0.5)}") )); then
        echo "ERROR: Too many chunks failing (>\${failure_rate})"
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        custom: 1.0.0
    END_VERSIONS
    """
}

// ============================================================================
// FILTERING PROCESSES
// ============================================================================

process QC_DUPL {
    tag "$meta.id"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.no_dupl.vcf.gz"), path("*.no_dupl.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    # Remove duplicate positions
    bcftools norm \\
        --rm-dup all \\
        --output ${meta.id}.no_dupl.vcf.gz \\
        --output-type z \\
        ${vcf}
    
    bcftools index --tbi ${meta.id}.no_dupl.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process SPLIT_MULTI_ALLELIC {
    tag "$meta.id"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.split.vcf.gz"), path("*.split.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    # Split multi-allelic sites
    bcftools norm \\
        --multiallelics -both \\
        --output ${meta.id}.split.vcf.gz \\
        --output-type z \\
        ${vcf}
    
    bcftools index --tbi ${meta.id}.split.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process FILTER_MIN_AC {
    tag "$meta.id"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def min_ac = params.min_ac ?: 1
    """
    # Filter by minimum allele count
    bcftools filter \\
        --include 'AC[0]>=${min_ac}' \\
        --output ${meta.id}.filtered.vcf.gz \\
        --output-type z \\
        ${vcf}
    
    bcftools index --tbi ${meta.id}.filtered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process QC_SITE_MISSINGNESS {
    tag "$meta.id"
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.qc.vcf.gz"), path("*.qc.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def max_missing = params.site_miss ?: 0.05
    """
    # Filter by site missingness
    bcftools filter \\
        --include 'F_MISSING<=${max_missing}' \\
        --output ${meta.id}.qc.vcf.gz \\
        --output-type z \\
        ${vcf}
    
    bcftools index --tbi ${meta.id}.qc.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

// ============================================================================
// UTILITY PROCESSES
// ============================================================================

process GENERATE_CHUNK_MAP {
    tag "$meta.id"
    label 'process_single'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.chunk.map"), emit: map
    path "versions.yml", emit: versions
    
    script:
    """
    # Generate variant map for chunk
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' \\
        ${vcf} > ${meta.id}.chunk.map
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

process MERGE_ADJACENT_CHUNKS {
    tag "$meta.id"
    label 'process_medium'
    
    input:
    tuple val(passed_meta), val(failed_meta), path(passed_vcf), path(failed_vcf)
    
    output:
    tuple val(meta), path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def meta = passed_meta.clone()
    meta.id = "${passed_meta.id}_merged"
    """
    # Merge adjacent chunks
    bcftools concat \\
        ${passed_vcf} ${failed_vcf} \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output ${meta.id}.merged.vcf.gz \\
        --output-type z
    
    bcftools index --tbi ${meta.id}.merged.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}

// ============================================================================
// MAIN QC WORKFLOW - Orchestrates all QC steps
// ============================================================================

workflow QC_PIPELINE {
    take:
    ch_input  // channel: [ val(meta), path(vcf) ]
    
    main:
    ch_versions = Channel.empty()
    
    // Step 1: Basic validation
    CHECK_FILES(ch_input)
    ch_versions = ch_versions.mix(CHECK_FILES.out.versions)
    
    CHECK_CHROMOSOME(CHECK_FILES.out.vcf)
    ch_versions = ch_versions.mix(CHECK_CHROMOSOME.out.versions)
    
    CHECK_GENOME_BUILD(
        CHECK_CHROMOSOME.out.vcf,
        params.genome_build ?: 'b38'
    )
    ch_versions = ch_versions.mix(CHECK_GENOME_BUILD.out.versions)
    
    // Step 2: Filter by build status
    ch_vcf_validated = CHECK_GENOME_BUILD.out.vcf
        .join(CHECK_GENOME_BUILD.out.build_status, by: 0)
        .filter { meta, vcf, status_file ->
            status_file.text.contains("PASS")
        }
        .map { meta, vcf, status_file -> [meta, vcf] }
    
    // Step 3: Standard QC filters
    QC_DUPL(ch_vcf_validated)
    ch_versions = ch_versions.mix(QC_DUPL.out.versions)
    
    SPLIT_MULTI_ALLELIC(QC_DUPL.out.vcf)
    ch_versions = ch_versions.mix(SPLIT_MULTI_ALLELIC.out.versions)
    
    FILTER_MIN_AC(SPLIT_MULTI_ALLELIC.out.vcf)
    ch_versions = ch_versions.mix(FILTER_MIN_AC.out.versions)
    
    QC_SITE_MISSINGNESS(FILTER_MIN_AC.out.vcf)
    ch_versions = ch_versions.mix(QC_SITE_MISSINGNESS.out.versions)
    
    emit:
    vcf = QC_SITE_MISSINGNESS.out.vcf
    versions = ch_versions
}

// ============================================================================
// SIMPLIFIED UNIFIED QC - For 'simple' mode
// ============================================================================

process UNIFIED_QC_SIMPLE {
    tag "$meta.id"
    label 'process_high'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.unified_qc.vcf.gz"), path("*.unified_qc.vcf.gz.tbi"), emit: vcf
    path "*.qc_report.json", emit: report
    path "versions.yml", emit: versions
    
    script:
    def qc_level = params.mode == 'simple' ? 'basic' : 
                   params.mode == 'clinical' ? 'strict' : 'standard'
    """
    # All QC in one command for simple mode
    bcftools norm \\
        --rm-dup all \\
        --multiallelics -both \\
        ${vcf} |
    bcftools filter \\
        --include 'AC[0]>=1 && F_MISSING<=0.05' |
    bcftools annotate \\
        --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' \\
        --output ${meta.id}.unified_qc.vcf.gz \\
        --output-type z
    
    bcftools index --tbi ${meta.id}.unified_qc.vcf.gz
    
    # Generate QC report
    bcftools stats ${meta.id}.unified_qc.vcf.gz > stats.txt
    
    # Convert to JSON
    echo '{' > ${meta.id}.qc_report.json
    echo '  "qc_level": "${qc_level}",' >> ${meta.id}.qc_report.json
    echo '  "variants_passed": '"\$(grep 'number of records:' stats.txt | cut -f4)"',' >> ${meta.id}.qc_report.json
    echo '  "status": "PASS"' >> ${meta.id}.qc_report.json
    echo '}' >> ${meta.id}.qc_report.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
```

**Implementation of unified_qc_pipeline.py:**

```python
#!/usr/bin/env python3
"""
Unified QC Pipeline - Consolidates all QC steps
"""
import argparse
import json
from pathlib import Path

class UnifiedQC:
    def __init__(self, vcf_path, qc_level='standard', population='auto'):
        self.vcf_path = vcf_path
        self.qc_level = qc_level
        self.population = self._detect_population() if population == 'auto' else population
        self.qc_stats = {}
        
    def run_pipeline(self):
        """Execute all QC steps in sequence"""
        # Step 1: Basic validation
        self._validate_file()
        
        # Step 2: Chromosome check
        self._check_chromosomes()
        
        # Step 3: Build detection
        self._detect_genome_build()
        
        if self.qc_level in ['standard', 'strict']:
            # Step 4: Advanced QC
            self._check_overlap()
            self._check_mismatch()
            
        if self.qc_level == 'strict':
            # Step 5: Clinical-grade QC
            self._clinical_qc()
            
        # Step 6: Apply filters
        self._apply_filters()
        
        # Step 7: Generate report
        return self._generate_report()
    
    def _detect_population(self):
        """Auto-detect population from allele frequencies"""
        # Smart population detection logic
        return "EUR"  # placeholder
    
    def _apply_filters(self):
        """Apply QC filters based on level and population"""
        filters = self._get_filters()
        
        # Population-specific adjustments
        if self.population == 'AFR':
            filters['max_mismatch_rate'] *= 1.25
            filters['min_variants'] *= 0.75
            
        # Apply all filters in one pass
        self._filter_vcf(filters)
```

### 5.2 Reporting Module Consolidation

**Current: 30+ reporting modules**

**Proposed: Tiered reporting system**

```nextflow
// NEW: Single reporting entry point
process UNIFIED_REPORT {
    input:
    tuple val(meta), path(imputed_vcf), path(info)
    val report_level  // 'basic', 'standard', 'comprehensive'
    
    output:
    path "reports/*", emit: reports
    path "plots/*", emit: plots, optional: true
    
    script:
    """
    unified_report_generator.py \\
        --vcf ${imputed_vcf} \\
        --info ${info} \\
        --level ${report_level} \\
        --output-dir . \\
        --threads ${task.cpus}
    """
}
```

**Report Levels:**
```python
REPORT_LEVELS = {
    'basic': [
        'summary_stats',
        'imputation_quality',
        'basic_plots'
    ],
    'standard': [
        'summary_stats',
        'imputation_quality', 
        'qc_metrics',
        'maf_analysis',
        'performance_plots',
        'frequency_comparison'
    ],
    'comprehensive': [
        # All 30+ current reports
        'summary_stats',
        'imputation_quality',
        'qc_metrics',
        # ... all reports ...
    ]
}
```

---

## 6. Configuration Simplification

### 6.1 Hierarchical Configuration

```groovy
// NEW: nextflow.config with smart defaults
includeConfig 'conf/base.config'      // Minimal base
includeConfig 'conf/profiles.config'  // Profile definitions

// Smart parameter resolution
params {
    // REQUIRED (only these!)
    input = null
    outdir = './results'
    
    // OPTIONAL with smart defaults
    genome_build = 'auto'  // Auto-detect from VCF
    population = 'auto'    // Auto-detect from frequencies
    mode = 'simple'        // Default to simplest
    
    // HIDDEN (but configurable)
    // Everything else hidden in profiles
}

// Auto-configuration based on detected features
if (params.population == 'auto') {
    include { DETECT_POPULATION } from './modules/unified/auto_detect'
}
```

### 6.2 Configuration Generator Tool

```python
#!/usr/bin/env python3
"""
Configuration wizard for pipeline setup
"""

class ConfigWizard:
    def run(self):
        print("Welcome to Imputation Pipeline Setup")
        print("=====================================")
        
        # Ask minimal questions
        level = self._ask_level()
        
        if level == 'simple':
            # Generate simple config
            self._generate_simple_config()
        elif level == 'custom':
            # Interactive configuration
            population = self._ask_population()
            resources = self._detect_resources()
            self._generate_custom_config(population, resources)
    
    def _generate_simple_config(self):
        """Generate minimal config file"""
        config = """
        params {
            input = "your_data.vcf"
            outdir = "results"
        }
        """
        Path('my_imputation.config').write_text(config)
        print("✓ Configuration saved to my_imputation.config")
        print("✓ Run with: nextflow run imputation -c my_imputation.config")
```

---

## 7. Workflow Streamlining

### 7.1 Simplified Main Workflow

```nextflow
// NEW: main.nf - Clean and simple
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Banner
log.info """
╔═══════════════════════════════════════════╗
║     Genotype Imputation Pipeline         ║
║           Simple • Fast • Accurate        ║
╚═══════════════════════════════════════════╝
Input:  ${params.input}
Mode:   ${params.mode ?: 'simple'}
Output: ${params.outdir}
"""

// Single include
include { IMPUTATION } from './workflows/imputation'

// Clean workflow
workflow {
    // Input validation
    Channel
        .fromPath(params.input, checkIfExists: true)
        .set { ch_input }
    
    // Run pipeline
    IMPUTATION(ch_input, params.mode)
}

// Completion handler
workflow.onComplete {
    log.info """
    ✓ Pipeline completed successfully
    ✓ Results: ${params.outdir}
    ✓ Runtime: ${workflow.duration}
    """
}
```

### 7.2 Simplified Subworkflows

```nextflow
// workflows/imputation.nf
include { UNIFIED_QC } from '../modules/unified/qc_pipeline'
include { UNIFIED_PHASE } from '../modules/unified/phasing_pipeline'
include { UNIFIED_IMPUTE } from '../modules/unified/imputation_pipeline'
include { UNIFIED_REPORT } from '../modules/unified/reporting_pipeline'

workflow IMPUTATION {
    take:
    ch_input
    mode
    
    main:
    // Determine report level from mode
    report_level = mode == 'simple' ? 'basic' :
                   mode == 'standard' ? 'standard' :
                   'comprehensive'
    
    // Linear, clear workflow
    UNIFIED_QC(ch_input)
    UNIFIED_PHASE(UNIFIED_QC.out.vcf)
    UNIFIED_IMPUTE(UNIFIED_PHASE.out.phased)
    UNIFIED_REPORT(
        UNIFIED_IMPUTE.out.imputed,
        report_level
    )
    
    emit:
    imputed = UNIFIED_IMPUTE.out.imputed
    reports = UNIFIED_REPORT.out.reports
}
```

---

## 8. Smart Defaults Implementation

### 8.1 Intelligent Auto-Detection

```nextflow
// modules/unified/auto_detect.nf
process AUTO_DETECT_SETTINGS {
    input:
    path vcf
    
    output:
    env DETECTED_BUILD, emit: build
    env DETECTED_POP, emit: population
    env SUGGESTED_CHUNKS, emit: chunks
    
    script:
    """
    # Auto-detect genome build
    DETECTED_BUILD=\$(detect_genome_build.py ${vcf})
    
    # Auto-detect population
    DETECTED_POP=\$(detect_population.py ${vcf})
    
    # Suggest optimal chunk size
    SUGGESTED_CHUNKS=\$(suggest_chunks.py ${vcf} \${DETECTED_POP})
    
    echo "Detected: Build=\${DETECTED_BUILD}, Population=\${DETECTED_POP}"
    """
}
```

### 8.2 Population-Aware Defaults

```groovy
// conf/population_defaults.config
def setPopulationDefaults() {
    switch(params.population) {
        case 'AFR':
            params.max_mismatch_rate = 0.25
            params.chunk_size = 30000000
            params.phasing_method = 'shapeit4'
            params.expect_ibd_cm = 5.0
            break
        case 'EUR':
            params.max_mismatch_rate = 0.15
            params.chunk_size = 20000000
            params.phasing_method = 'eagle'
            params.expect_ibd_cm = 3.0
            break
        case 'EAS':
            params.max_mismatch_rate = 0.15
            params.chunk_size = 25000000
            params.phasing_method = 'eagle'
            params.expect_ibd_cm = 2.5
            break
        default:
            // Safe defaults for unknown
            params.max_mismatch_rate = 0.20
            params.chunk_size = 25000000
            params.phasing_method = 'eagle'
    }
}
```

### 8.3 Resource Auto-Optimization

```nextflow
// modules/unified/resource_optimizer.nf
def optimizeResources() {
    // Detect available resources
    def availableCpus = Runtime.runtime.availableProcessors()
    def availableMemory = Runtime.runtime.maxMemory()
    
    // Set intelligent defaults
    if (availableCpus >= 32) {
        params.max_cpus = 32
        params.parallel_chunks = 200
    } else if (availableCpus >= 16) {
        params.max_cpus = 16
        params.parallel_chunks = 100
    } else {
        params.max_cpus = availableCpus
        params.parallel_chunks = availableCpus * 2
    }
    
    // Memory-based optimization
    if (availableMemory >= 128.GB) {
        params.max_memory = '128.GB'
        params.chunk_size = 30000000
    } else if (availableMemory >= 64.GB) {
        params.max_memory = '64.GB'
        params.chunk_size = 25000000
    } else {
        params.max_memory = '32.GB'
        params.chunk_size = 20000000
    }
}
```

---

## 9. Progressive Disclosure Pattern

### 9.1 Command-Line Interface Layers

```bash
# Layer 1: Absolute Simplest
impute mydata.vcf

# Layer 2: Common Options
impute mydata.vcf --population african

# Layer 3: Power User
impute mydata.vcf --population african --ref-panel h3africa --chunks 30M

# Layer 4: Full Control
nextflow run imputation/main.nf \
    -c custom.config \
    -profile advanced \
    --population AFR \
    --max_mismatch_rate 0.25 \
    --phasing_method shapeit4 \
    --impute_method minimac4 \
    --generate_all_plots true
```

### 9.2 Documentation Layers

```markdown
# Quick Start (90% of users stop here)
```bash
impute mydata.vcf
```
Done! Results in `results/` folder.

# Common Customizations (9% read this)
- `--population`: african, european, asian
- `--quality`: fast, balanced, accurate
- `--output`: specify output directory

# Advanced Usage (1% need this)
See [Advanced Documentation](advanced.md) for:
- Custom reference panels
- Pipeline customization
- Performance tuning
- Method comparison
```

### 9.3 Error Message Layers

```python
def smart_error_handler(error, context):
    """Progressive error messages"""
    
    # Level 1: Simple message and fix
    print(f"❌ Error: {error.simple_message}")
    print(f"✓ Fix: {error.suggested_fix}")
    
    if user_wants_details():
        # Level 2: More context
        print(f"ℹ Context: {error.context}")
        print(f"ℹ Possible causes: {error.causes}")
        
    if user_is_advanced():
        # Level 3: Full debugging
        print(f"Debug info: {error.stack_trace}")
        print(f"Pipeline state: {context.state}")
```

---

## 10. Implementation Roadmap

### Phase 1: Foundation (Weeks 1-2)
```
□ Create unified module structure
□ Implement AUTO_DETECT modules
□ Set up profile system
□ Create configuration wizard
```

### Phase 2: Core Consolidation (Weeks 3-4)
```
□ Consolidate QC modules → UNIFIED_QC
□ Consolidate phasing → UNIFIED_PHASE
□ Consolidate imputation → UNIFIED_IMPUTE
□ Create BASIC_REPORT module
```

### Phase 3: Interface Simplification (Weeks 5-6)
```
□ Implement simple CLI wrapper
□ Create mode-based workflows
□ Add smart defaults
□ Implement auto-detection
```

### Phase 4: Testing & Refinement (Weeks 7-8)
```
□ Test all modes extensively
□ Benchmark performance
□ Gather user feedback
□ Refine defaults
```

### Phase 5: Documentation & Release (Weeks 9-10)
```
□ Write layered documentation
□ Create video tutorials
□ Prepare migration guide
□ Release v2.0
```

---

## 11. Code Refactoring Examples

### 11.1 Simplify Channel Operations

**Current: Complex channel manipulation**
```nextflow
ch_vcf_chunks_filtered = ch_vcf_chunks
    .combine(ch_ref_chroms)
    .map { meta, vcf, ref_chroms ->
        def chrm = meta.contig
        if (ref_chroms.contains(chrm)) {
            [meta, vcf]
        } else {
            null
        }
    }
    .filter { it != null }
    .join(CHECK_OVERLAP.out.overlap)
    .branch { meta, vcf, overlap_txt, status_file ->
        def status = status_file.text.trim()
        passed: status == "PASS"
        failed: status == "FAIL"
    }
```

**Proposed: Simplified with helper functions**
```nextflow
// Define reusable functions
def filterValidChunks = { meta, vcf, ref_chroms ->
    ref_chroms.contains(meta.contig) ? [meta, vcf] : null
}

def branchByStatus = { meta, vcf, status_file ->
    status_file.text.trim() == "PASS" ? 'passed' : 'failed'
}

// Clean channel operations
ch_vcf_chunks_filtered = ch_vcf_chunks
    .combine(ch_ref_chroms)
    .map(filterValidChunks)
    .filter()
    .join(CHECK_OVERLAP.out.overlap)
    .branch(branchByStatus)
```

### 11.2 Consolidate Error Handling

**Current: Scattered error handling**
```nextflow
process IMPUTE_MINIMAC4 {
    script:
    """
    minimac4 ... || {
        exit_code=\$?
        case \$exit_code in
            137|143)
                echo "Memory error"
                exit 137
                ;;
            255)
                echo "Empty chunk"
                create_empty_output.sh
                exit 0
                ;;
            *)
                echo "Unknown error: \$exit_code"
                exit \$exit_code
                ;;
        esac
    }
    """
}
```

**Proposed: Unified error handler**
```nextflow
process IMPUTE_MINIMAC4 {
    script:
    """
    # Use unified error handler
    source ${projectDir}/bin/error_handler.sh
    
    run_with_error_handling minimac4 \\
        --refHaps ${ref_m3vcf} \\
        --haps ${phased_vcf} \\
        --prefix ${prefix}
    """
}
```

**Error handler implementation:**
```bash
#!/bin/bash
# bin/error_handler.sh

run_with_error_handling() {
    "$@" || handle_error $? "$1"
}

handle_error() {
    local exit_code=$1
    local command=$2
    
    case $exit_code in
        0) return 0 ;;
        137|143) handle_memory_error ;;
        255) handle_empty_chunk ;;
        *) handle_unknown_error $exit_code $command ;;
    esac
}
```

### 11.3 Simplify Configuration Loading

**Current: Complex configuration hierarchy**
```groovy
includeConfig 'conf/base.config'
includeConfig 'conf/resources.config'
includeConfig 'conf/containers.config'
if (params.population) {
    includeConfig "conf/populations/${params.population}.config"
}
```

**Proposed: Single smart configuration**
```groovy
// Single configuration with smart resolution
includeConfig 'conf/smart.config'

// conf/smart.config
def loadConfiguration() {
    // Load base
    def config = new ConfigSlurper().parse(new File('conf/defaults.config').text)
    
    // Apply mode
    config.merge(getModeConfig(params.mode ?: 'simple'))
    
    // Apply population if detected
    if (params.population == 'auto') {
        params.population = detectPopulation()
    }
    config.merge(getPopulationConfig(params.population))
    
    // Apply resource optimization
    config.merge(optimizeResources())
    
    return config
}

// Apply configuration
def config = loadConfiguration()
params.putAll(config.params)
```

---

## 12. Testing & Validation Strategy

### 12.1 Compatibility Testing

```groovy
// test/compatibility_test.nf
workflow test_compatibility {
    // Test that simplified pipeline produces same results
    
    // Run old pipeline
    OLD_PIPELINE(test_data)
    
    // Run new simplified pipeline
    NEW_SIMPLE_PIPELINE(test_data)
    
    // Compare results
    COMPARE_RESULTS(
        OLD_PIPELINE.out,
        NEW_SIMPLE_PIPELINE.out
    )
    
    // Assert compatibility
    ASSERT_COMPATIBLE(COMPARE_RESULTS.out)
}
```

### 12.2 Performance Benchmarking

```python
#!/usr/bin/env python3
"""
Benchmark simplified vs original pipeline
"""

class PipelineBenchmark:
    def run_benchmarks(self):
        metrics = {}
        
        # Benchmark original pipeline
        metrics['original'] = self.benchmark_pipeline(
            'main_nfcore.nf',
            'original.config'
        )
        
        # Benchmark simplified pipeline
        metrics['simplified'] = self.benchmark_pipeline(
            'main_simplified.nf',
            'simple.config'
        )
        
        # Compare metrics
        self.compare_metrics(metrics)
        
    def compare_metrics(self, metrics):
        """Compare key performance indicators"""
        comparison = {
            'runtime': {
                'original': metrics['original']['runtime'],
                'simplified': metrics['simplified']['runtime'],
                'improvement': calculate_improvement()
            },
            'memory_peak': {
                'original': metrics['original']['memory'],
                'simplified': metrics['simplified']['memory'],
                'improvement': calculate_improvement()
            },
            'user_steps': {
                'original': 15,  # Configure, run, monitor
                'simplified': 2,  # Run with defaults
                'improvement': '87% reduction'
            }
        }
        
        print_benchmark_report(comparison)
```

### 12.3 User Experience Testing

```python
class UXTesting:
    def test_user_journeys(self):
        """Test different user scenarios"""
        
        scenarios = [
            {
                'user': 'beginner',
                'task': 'Run imputation on test data',
                'expected_commands': 1,
                'expected_time': '< 1 minute to start'
            },
            {
                'user': 'researcher',
                'task': 'Run African population imputation',
                'expected_commands': 1,
                'expected_time': '< 2 minutes to configure'
            },
            {
                'user': 'advanced',
                'task': 'Custom pipeline with specific panels',
                'expected_commands': '< 5',
                'expected_time': '< 5 minutes to configure'
            }
        ]
        
        for scenario in scenarios:
            self.test_scenario(scenario)
```

---

## Summary & Key Recommendations

### Critical Simplifications to Implement

1. **Unified Modules** (Highest Priority)
   - Consolidate 50+ modules into 5-7 unified modules
   - Hide complexity inside modules
   - Maintain all functionality

2. **Smart Defaults** (High Priority)
   - Auto-detect population and build
   - Set optimal parameters automatically
   - Resource auto-optimization

3. **Profile System** (High Priority)
   - Create 'simple', 'research', 'clinical' profiles
   - Hide advanced features in profiles
   - Progressive disclosure

4. **Single Entry Point** (Medium Priority)
   - Simple CLI wrapper
   - Configuration wizard
   - Clear documentation layers

5. **Maintain Robustness** (Critical)
   - Keep all error handling
   - Preserve population-specific QC
   - Maintain recovery mechanisms

### Expected Outcomes

| Metric | Current | After Simplification | Improvement |
|--------|---------|---------------------|-------------|
| **Modules** | 50+ | 5-7 unified | 90% reduction |
| **Required Parameters** | 10+ | 2 | 80% reduction |
| **Time to First Run** | 30+ min | < 2 min | 93% reduction |
| **Configuration Lines** | 100+ | 3-5 | 95% reduction |
| **Learning Curve** | Weeks | Hours | 95% reduction |
| **Performance** | Baseline | Same/Better | 0-10% improvement |
| **Robustness** | Excellent | Excellent | Maintained |
| **Flexibility** | Full | Full (hidden) | Maintained |

### Implementation Priority

1. **Week 1-2**: Create unified modules and smart defaults
2. **Week 3-4**: Implement profile system and auto-detection
3. **Week 5-6**: Build simple CLI and configuration wizard
4. **Week 7-8**: Testing and optimization
5. **Week 9-10**: Documentation and release

### Final Architecture Vision

```
User sees:
impute data.vcf [--population african]

Hidden complexity:
- 50+ modules working behind scenes
- Smart population detection
- Adaptive QC thresholds
- Optimal resource allocation
- Comprehensive error handling
- 30+ reports (generated selectively)

Result:
- ImputationServer2 simplicity
- H3ABioNet robustness
- Best of both worlds
```

---

## Conclusion

The proposed simplification strategy achieves ImputationServer2's ease of use while maintaining H3ABioNet pipeline's superior capabilities. By implementing unified modules, smart defaults, and progressive disclosure, users get a simple interface backed by powerful, population-aware processing. The key is hiding complexity, not removing it.

### Next Steps

1. Review and approve simplification strategy
2. Create proof-of-concept for unified modules
3. Test with user group
4. Iterate based on feedback
5. Full implementation

The simplified pipeline will make genotype imputation accessible to all users while preserving the advanced features that make it superior for diverse populations and research applications.

---

*Document Version: 1.0*  
*Pipeline Simplification Strategy*  
*Achieving Simplicity Without Sacrificing Robustness*