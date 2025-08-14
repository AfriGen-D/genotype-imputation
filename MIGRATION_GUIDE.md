# Complete Migration Guide: v1.x to v2.0

## Table of Contents
1. [Overview](#overview)
2. [Pre-Migration Checklist](#pre-migration-checklist)
3. [Step-by-Step Migration](#step-by-step-migration)
4. [Parameter Mapping](#parameter-mapping)
5. [Testing Migration](#testing-migration)
6. [Rollback Plan](#rollback-plan)
7. [FAQ](#faq)

## Overview

This guide provides comprehensive instructions for migrating from the v1.x pipeline to the refactored v2.0 architecture.

### Key Changes
- **Modular architecture**: 5 focused subworkflows instead of 9 mixed workflows
- **Unified reporting**: Single reporting system for all output types
- **Modern DSL2**: Following Nextflow and nf-core best practices
- **Improved configuration**: Cleaner parameter organization

### Migration Timeline
- **Phase 1** (Weeks 1-2): Test v2.0 with sample data
- **Phase 2** (Weeks 3-4): Run parallel validation
- **Phase 3** (Week 5): Production switchover
- **Phase 4** (Week 6+): Decommission v1.x

## Pre-Migration Checklist

Before starting migration, ensure you have:

- [ ] Backed up current configuration files
- [ ] Documented custom modifications to v1.x
- [ ] Listed all active workflows and their parameters
- [ ] Verified container/singularity image availability
- [ ] Allocated time for testing (2-4 hours)
- [ ] Access to test datasets

## Step-by-Step Migration

### Step 1: Backup Current Setup

```bash
# Create backup directory
mkdir -p migration_backup_$(date +%Y%m%d)
cd migration_backup_$(date +%Y%m%d)

# Backup configurations
cp ../nextflow.config ./
cp ../v6_chr21_phased.config ./
cp -r ../conf ./

# Document current version
git rev-parse HEAD > current_version.txt
nextflow log > pipeline_runs.txt
```

### Step 2: Install v2.0 Pipeline

```bash
# Clone or switch to v2.0 branch
cd /users/mamana/genotype-imputation
git pull origin master

# Verify installation
cd workflows
bash bin/validate_pipeline.sh
```

### Step 3: Convert Configuration

```bash
# Automatic conversion
cd workflows
python bin/migrate_config.py \
  ../v6_chr21_phased.config \
  -o v6_chr21_phased_v2.config \
  --create-samplesheet \
  --verbose

# Review generated files
cat v6_chr21_phased_v2.config
cat v6_chr21_phased_v2_samplesheet.csv
```

### Step 4: Update File Paths

Edit the migrated configuration to update paths:

```groovy
// v6_chr21_phased_v2.config
params {
    // Update these paths to your actual locations
    input = '/actual/path/to/samplesheet.csv'
    reference_panels = [
        ['H3AR6x', 
         '/cbio/dbs/refpanels/h3a/v6/H3AR6x.m3vcf.gz',
         '/cbio/dbs/refpanels/h3a/v6/H3AR6x.vcf.gz']
    ]
    eagle_genetic_map = '/actual/path/to/genetic_map.txt.gz'
    reference_genome = '/actual/path/to/reference.fa'
}
```

### Step 5: Create Samplesheet

Convert your target datasets to CSV format:

```csv
sample,vcf,population,sex
sample1,/path/to/sample1.vcf.gz,AFR,M
sample2,/path/to/sample2.vcf.gz,EUR,F
sample3,/path/to/sample3.vcf.gz,AMR,U
```

### Step 6: Test Migration

```bash
# Dry run with stub
cd workflows
nextflow run main.nf \
  -c v6_chr21_phased_v2.config \
  -profile test,singularity \
  -stub

# Small test run
nextflow run main.nf \
  -c v6_chr21_phased_v2.config \
  -profile singularity \
  --chromosomes "22" \
  -resume
```

## Parameter Mapping

### Input/Output Parameters

| v1.x Parameter | v2.0 Parameter | Type | Notes |
|----------------|----------------|------|-------|
| `target_datasets` | `input` | String | Now points to CSV samplesheet |
| `ref_panels` | `reference_panels` | Array | JSON array format |
| `outDir` | `outdir` | String | Lowercase convention |
| `output_dir` | `outdir` | String | Standardized |

### QC Parameters

| v1.x Parameter | v2.0 Parameter | Default | Notes |
|----------------|----------------|---------|-------|
| `site_miss` | `qc_max_missing` | 0.05 | Inverted logic (max vs min) |
| `hwe` | `qc_hwe_pvalue` | 1e-6 | Same threshold |
| `mac` | `qc_min_ac` | 2 | Renamed for clarity |
| `min_ac` | `qc_min_ac` | 2 | Consolidated |
| `min_alleles` | - | - | Removed (always 2) |

### Imputation Parameters

| v1.x Parameter | v2.0 Parameter | Default | Notes |
|----------------|----------------|---------|-------|
| `NE` | `impute_ne` | 20000 | Effective population size |
| `impute_iter` | - | - | Tool-specific, removed |
| `impute_burnin` | - | - | Tool-specific, removed |
| `impute_info_cutoff` | `impute_info_cutoff` | 0.3 | Unchanged |
| `chunk_size` | `chunk_size` | 5000000 | For parallelization |
| `buffer_size` | `impute_buffer` | 250000 | Renamed |

### New Parameters in v2.0

| Parameter | Default | Description |
|-----------|---------|-------------|
| `phasing_tool` | 'eagle' | Choose: eagle, shapeit, beagle |
| `imputation_tool` | 'minimac4' | Choose: minimac4, impute5, beagle5 |
| `report_level` | 'detailed' | Choose: summary, detailed, full |
| `genome_build` | 'b38' | Choose: b37, b38 |
| `concordance_analysis` | false | Enable concordance checks |

## Testing Migration

### 1. Validate Configuration

```bash
cd workflows
cat > test_migration.nf << 'EOF'
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Test parameter loading
workflow {
    log.info "Testing parameter migration..."
    log.info "Input: ${params.input}"
    log.info "Reference panels: ${params.reference_panels}"
    log.info "QC min AC: ${params.qc_min_ac}"
    log.info "Imputation tool: ${params.imputation_tool}"
}
EOF

nextflow run test_migration.nf -c v6_chr21_phased_v2.config
```

### 2. Compare Outputs

```bash
# Run both pipelines on test data
# v1.x pipeline
nextflow run main.nf \
  -c v6_chr21_phased.config \
  -profile singularity \
  --chromosomes 22 \
  --outDir output_v1

# v2.0 pipeline
cd workflows
nextflow run main.nf \
  -c v6_chr21_phased_v2.config \
  -profile singularity \
  --chromosomes 22 \
  --outdir output_v2

# Compare outputs
diff -r ../output_v1/imputed output_v2/imputed
```

### 3. Validate Results

```bash
# Check imputation quality metrics
bcftools stats output_v2/imputed/*.vcf.gz > v2_stats.txt
bcftools stats ../output_v1/imputed/*.vcf.gz > v1_stats.txt
diff v1_stats.txt v2_stats.txt

# Compare INFO scores
zcat output_v2/metrics/*.info.gz | head
zcat ../output_v1/info/*.info.gz | head
```

## Rollback Plan

If issues arise, you can quickly rollback:

### Immediate Rollback

```bash
# Stop v2.0 pipeline
nextflow stop [run_name]

# Return to v1.x
cd /users/mamana/genotype-imputation
nextflow run main.nf -c v6_chr21_phased.config -resume
```

### Configuration Rollback

```bash
# Restore backed up configs
cp migration_backup_*/nextflow.config ./
cp migration_backup_*/v6_chr21_phased.config ./
```

### Data Recovery

```bash
# v2.0 outputs are in separate directory
# Original data remains untouched
ls -la output/       # v1.x outputs
ls -la workflows/results/  # v2.0 outputs
```

## FAQ

### Q: Can I run both versions in parallel?
**A:** Yes! The v2.0 pipeline is in the `workflows/` directory and uses different output paths. You can run both simultaneously for validation.

### Q: Will my existing outputs be affected?
**A:** No. The v2.0 pipeline writes to new output directories. Your existing results remain untouched.

### Q: How do I handle custom modifications?
**A:** Document your modifications and contact the development team. Most customizations can be ported to v2.0's modular structure more cleanly.

### Q: What if the migration script fails?
**A:** You can manually create the configuration using the template in `workflows/conf/test/test.config` as a starting point.

### Q: How long does migration take?
**A:** Configuration migration: 5-10 minutes. Testing: 1-2 hours. Full validation: 2-4 hours depending on data size.

### Q: Can I migrate incrementally?
**A:** Yes! Start with test datasets, then pilot projects, then full production. The pipelines can coexist.

## Support

For migration assistance:
1. Check the [REFACTORING_SUMMARY.md](REFACTORING_SUMMARY.md)
2. Review the [workflows/README.md](workflows/README.md)
3. Run the validation script: `bash workflows/bin/validate_pipeline.sh`
4. Open an issue on GitHub with the `migration` tag

## Migration Checklist Summary

- [ ] Backup current setup
- [ ] Install v2.0 pipeline
- [ ] Convert configuration files
- [ ] Create CSV samplesheet
- [ ] Update file paths
- [ ] Run validation tests
- [ ] Compare outputs with v1.x
- [ ] Document any issues
- [ ] Plan production switchover
- [ ] Monitor initial production runs
- [ ] Decommission v1.x (after 30 days)

## Success Criteria

Your migration is successful when:
- ✅ All tests pass validation
- ✅ Output VCFs have comparable metrics
- ✅ INFO scores match within 0.01
- ✅ Runtime is similar or better
- ✅ No errors in Nextflow logs
- ✅ Reports generate correctly

---

*Last updated: August 2024*
*Version: 2.0.0*