#!/bin/bash
# Test script for pre/post-imputation comparison and R2 genomic window analyses
# This script demonstrates how to run the newly enabled reporting modules

set -e

echo "======================================"
echo "Testing Imputation Reporting Modules"
echo "======================================"

# Configuration
WORK_DIR="/scratch3/users/mamana/nextflow-work"
RESULTS_DIR="/scratch3/users/mamana/results"
CONFIG_FILE="v6_chr21_phased_nfcore.config"

# Create results directory if it doesn't exist
mkdir -p ${RESULTS_DIR}/reports/{pre_post_comparison,imputation_quality}

echo ""
echo "1. Running pipeline with reporting modules enabled..."
echo "   - COMPARE_PRE_POST_IMPUTATION: Compares pre and post-imputation VCFs"
echo "   - PLOT_R2_GENOMIC_WINDOWS: Analyzes R2 quality in genomic windows"

# Test with a small dataset first (chr21 or chr22)
echo ""
echo "2. Test run with chromosome 21 data..."

nextflow run main_nfcore.nf \
    -c ${CONFIG_FILE} \
    -profile slurm,singularity \
    -resume \
    --input test_sample.csv \
    --output ${RESULTS_DIR} \
    --chromosome chr21 \
    --max_memory 16.GB \
    --max_cpus 8 \
    --max_time 24.h \
    -with-report ${RESULTS_DIR}/reports/execution_report.html \
    -with-trace ${RESULTS_DIR}/reports/execution_trace.txt \
    -with-timeline ${RESULTS_DIR}/reports/execution_timeline.html \
    -with-dag ${RESULTS_DIR}/reports/pipeline_dag.svg

echo ""
echo "3. Expected outputs:"
echo "   Pre/Post Comparison Reports:"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/comparison_stats.txt"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/variant_counts.png"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/maf_distribution.png"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/coverage_improvement.png"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/variant_gain.png"
echo "   - ${RESULTS_DIR}/reports/pre_post_comparison/*/comparison_report.html"
echo ""
echo "   R2 Genomic Window Analysis:"
echo "   - ${RESULTS_DIR}/reports/imputation_quality/*/r2_windows.png"
echo "   - ${RESULTS_DIR}/reports/imputation_quality/*/poor_regions.txt"
echo "   - ${RESULTS_DIR}/reports/imputation_quality/*/r2_heatmap.png"
echo "   - ${RESULTS_DIR}/reports/imputation_quality/*/r2_distribution.png"
echo "   - ${RESULTS_DIR}/reports/imputation_quality/*/window_stats.txt"

echo ""
echo "4. To run with both reference panels (v6 and v7) for comparison:"
echo "   Modify the config to use both panels and the pipeline will generate"
echo "   comparative reports showing imputation quality differences."

echo ""
echo "======================================"
echo "Module Features:"
echo "======================================"
echo ""
echo "COMPARE_PRE_POST_IMPUTATION:"
echo "  - Quantifies variant gain from imputation"
echo "  - Shows MAF distribution changes"
echo "  - Identifies which frequency categories benefit most"
echo "  - Generates comprehensive HTML report"
echo ""
echo "PLOT_R2_GENOMIC_WINDOWS:"
echo "  - Identifies poorly imputed genomic regions"
echo "  - Uses sliding windows (default 1Mb)"
echo "  - Highlights regions below R2 threshold (default 0.3)"
echo "  - Creates heatmaps and distribution plots"
echo ""
echo "======================================"

# Optional: Run specific analyses on existing data
if [ -d "${WORK_DIR}" ]; then
    echo ""
    echo "5. Running standalone analysis on existing imputation results..."
    
    # Find recent info files for R2 analysis
    INFO_FILES=$(find ${WORK_DIR} -name "*.info" -mtime -1 2>/dev/null | head -5)
    
    if [ ! -z "$INFO_FILES" ]; then
        echo "   Found info files for R2 analysis:"
        echo "$INFO_FILES" | head -3
        
        # You could run the Python scripts directly here if needed
        # for info_file in $INFO_FILES; do
        #     python3 bin/analyze_r2_windows.py --info $info_file --output ${RESULTS_DIR}/r2_analysis
        # done
    fi
fi

echo ""
echo "Script complete. Check the results directory for outputs."
echo "======================================"