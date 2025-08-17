#!/bin/bash
###############################################
# Production run script for ChipImputation   #
###############################################

set -e

# Set environment
export NXF_WORK="/scratch3/users/mamana/nextflow-work"
export NXF_TEMP="/scratch3/users/mamana/nextflow-temp"
export NXF_SINGULARITY_CACHEDIR="/users/mamana/genotype-imputation/singularity_cache"

# Create directories if needed
mkdir -p $NXF_WORK
mkdir -p $NXF_TEMP
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Parse command line arguments
PROFILE="singularity"
RESUME=""
STUB=""
INPUT=""
OUTDIR="/scratch3/users/mamana/results"

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        --stub)
            STUB="-stub-run"
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --input FILE     Input samplesheet CSV file (required)"
            echo "  --outdir DIR     Output directory (default: /scratch3/users/mamana/results)"
            echo "  --profile PROF   Execution profile (default: singularity)"
            echo "  --resume         Resume from previous run"
            echo "  --stub           Run in stub mode (testing)"
            echo "  --help           Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Test run with stub mode"
            echo "  $0 --input test_sample.csv --stub"
            echo ""
            echo "  # Production run"
            echo "  $0 --input samples.csv --profile singularity,slurm"
            echo ""
            echo "  # Resume failed run"
            echo "  $0 --input samples.csv --resume"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 --help' for usage information"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$INPUT" ]; then
    echo "ERROR: Input file is required (--input)"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file does not exist: $INPUT"
    exit 1
fi

# Display run configuration
echo "=================================================="
echo "H3ABionet ChipImputation Pipeline"
echo "=================================================="
echo "Input:      $INPUT"
echo "Output:     $OUTDIR"
echo "Profile:    $PROFILE"
echo "Work dir:   $NXF_WORK"
echo "Resume:     $([ -n "$RESUME" ] && echo "Yes" || echo "No")"
echo "Stub mode:  $([ -n "$STUB" ] && echo "Yes" || echo "No")"
echo "=================================================="
echo ""

# Run the pipeline with clean configuration
# Using main_simple.nf to avoid configuration conflicts
nextflow run main_simple.nf \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --eagle_genetic_map null \
    --reference_genome null \
    --ref_panels "[]" \
    -profile "$PROFILE" \
    $RESUME \
    $STUB \
    -with-report "${OUTDIR}/pipeline_info/execution_report.html" \
    -with-timeline "${OUTDIR}/pipeline_info/execution_timeline.html" \
    -with-trace "${OUTDIR}/pipeline_info/execution_trace.txt"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Pipeline completed successfully!"
    echo "Results are in: $OUTDIR"
    echo "=================================================="
else
    echo ""
    echo "=================================================="
    echo "Pipeline failed. Check the logs for details."
    echo "=================================================="
    exit 1
fi