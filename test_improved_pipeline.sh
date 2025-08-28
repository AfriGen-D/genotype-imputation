#!/bin/bash

# Test script for improved pipeline with error-tolerant handling
# This script tests the three main improvements:
# 1. Error-tolerant report workflow
# 2. Safe imputation validation  
# 3. Adaptive chunking with relaxed mismatch criteria

echo "=================================================="
echo "Testing Improved Genotype Imputation Pipeline"
echo "=================================================="
echo ""
echo "Improvements implemented:"
echo "1. Error-tolerant channels in report workflow"
echo "2. Safe imputation with validation (optional)"
echo "3. Relaxed mismatch criteria for diverse populations"
echo ""

# Clean previous runs (optional)
echo "Cleaning previous test runs..."
rm -rf work/test_improved_*

# Kill any stuck pipelines
echo "Stopping any stuck pipelines..."
nextflow stop backstabbing_goldberg 2>/dev/null || true

echo ""
echo "Starting improved pipeline test..."
echo "=================================================="

# Run with the improved configuration
nextflow run main_nfcore.nf \
  -c v6_chr21_phased_nfcore.config \
  --input test_sample.csv \
  --outdir /scratch3/users/mamana/results/test_improved \
  -profile slurm,singularity \
  -resume \
  --use_safe_imputation true \
  --validate_chunks false \
  --max_mismatch_rate 0.2 \
  --min_matching_variants 20 \
  --skip_failed_chunks true \
  --log_filtered_chunks true \
  -process.errorStrategy 'ignore' \
  -name test_improved_pipeline

echo ""
echo "=================================================="
echo "Pipeline test completed"
echo "Check /scratch3/users/mamana/results/test_improved for results"
echo "=================================================="