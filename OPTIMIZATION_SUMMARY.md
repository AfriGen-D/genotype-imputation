# Pipeline Resource Optimization Summary

## Date: 2025-08-16

### Overview
Comprehensive resource optimization for genotype imputation pipeline to improve performance, reliability, and throughput.

## Key Optimizations Implemented

### 1. Process-Specific Resource Allocations

#### Phasing Processes
- **EAGLE_PHASE**: 16 CPUs, 64GB RAM (retry on memory errors)
- **SHAPEIT_PHASE**: 16 CPUs, 64GB RAM with 3 retries
- **BEAGLE_PHASE**: 16 CPUs, 64GB RAM with 3 retries

#### Imputation Processes
- **MINIMAC4_IMPUTE**: 16 CPUs, 64GB RAM, 48h time limit, 3 retries
- **IMPUTE5_IMPUTE**: 8 CPUs, 32GB RAM with threading support
- **BEAGLE5_IMPUTE**: 12 CPUs, 48GB RAM with extended time limits

#### QC Processes
- **REMOVE_DUPLICATES**: 8 CPUs, 32GB RAM with error retry
- **SPLIT_MULTIALLELIC**: 8 CPUs, 32GB RAM
- **FILTER_VARIANTS**: 8 CPUs, 32GB RAM
- **CHECK_VCF_FORMAT**: 8 CPUs, 32GB RAM

#### Chunking & Merging
- **CREATE_CHUNK_LIST**: 2 CPUs, 8GB RAM (lightweight)
- **EXTRACT_VCF_CHUNK**: 4 CPUs, 16GB RAM with 3 retries
- **CHECK_CHUNK_OVERLAP**: 4 CPUs, 16GB RAM
- **MERGE_PHASED_CHUNKS**: 12 CPUs, 48GB RAM
- **MERGE_IMPUTED_CHUNKS**: 12 CPUs, 48GB RAM
- **MERGE_LOW_OVERLAP_CHUNKS**: 8 CPUs, 32GB RAM

#### Reporting & Visualization
- **CALCULATE_METRICS**: 4 CPUs, 16GB RAM
- **PLOT_INFO_DISTRIBUTION**: 2 CPUs, 8GB RAM
- **CREATE_HTML_REPORT**: 4 CPUs, 16GB RAM
- **AGGREGATE_OVERLAP_REPORTS**: 4 CPUs, 16GB RAM

### 2. Memory Scaling Strategy

Implemented exponential memory scaling on retry:
```groovy
memory = { check_max(BASE_MEMORY * Math.pow(2, task.attempt - 1), 'memory') }
```

This ensures:
- 1st attempt: Base memory
- 2nd attempt: 2x base memory
- 3rd attempt: 4x base memory

### 3. Error Handling

Added intelligent error handling based on exit codes:
- **137, 139, 143, 247**: Memory-related errors → retry with increased resources
- Other errors: Process-specific handling (ignore/finish/terminate)

### 4. Global Resource Limits

Increased from baseline:
- **max_cpus**: 20 → 32
- **max_memory**: 128GB → 256GB
- **max_time**: 48h → 72h

### 5. SLURM Optimization

Enhanced SLURM executor configuration:
- **queueSize**: 20 → 50 (more parallel jobs)
- **submitRateLimit**: 10/min → 20/min (faster job submission)
- **exitReadTimeout**: Added 60 min timeout
- **killBatchSize**: 50 jobs (efficient cleanup)
- **perJobMemLimit**: Enforced per-job memory limits

### 6. Container-Specific Optimizations

Each process now has:
- Dedicated container assignment
- Process-appropriate resource allocation
- Optimized thread usage via ext.args

### 7. Publishing Strategy

Organized output directories by process type:
- `/chunks` - Chunking outputs
- `/qc` - QC results
- `/phased` - Phasing outputs
- `/imputed` - Imputation results
- `/plots` - Visualizations
- `/metrics` - Performance metrics
- `/stats` - Statistical summaries

## Performance Benefits

1. **Reduced Failures**: Progressive retry with memory scaling reduces OOM errors
2. **Better Parallelization**: Increased queue size allows more concurrent jobs
3. **Faster Recovery**: Smart error handling prevents unnecessary pipeline failures
4. **Improved Throughput**: Optimized CPU/memory allocation per process type
5. **Resource Efficiency**: Process-specific allocations prevent over-provisioning

## Monitoring Recommendations

1. Track memory usage patterns per process
2. Monitor retry rates for memory-intensive processes
3. Analyze job queue times vs execution times
4. Review SLURM efficiency reports

## Future Optimization Opportunities

1. Implement dynamic chunk sizing based on variant density
2. Add GPU acceleration for compatible tools
3. Implement adaptive resource allocation based on input size
4. Add caching for reference panel operations
5. Optimize I/O patterns for large VCF files

## Configuration Files Modified

- `/conf/modules.config` - Process-specific configurations
- `/v6_chr21_phased_v2.config` - Project-specific settings
- SLURM profile settings - Executor optimizations

## Testing Command

```bash
nextflow run main.nf \
  --input samplesheet_chr21.csv \
  -c v6_chr21_phased_v2.config \
  -profile slurm,singularity \
  -resume
```

## Notes

- All optimizations maintain backward compatibility
- Settings can be overridden in project-specific configs
- Resource limits respect cluster maximums via check_max()
- Error strategies are non-destructive (retry before fail)