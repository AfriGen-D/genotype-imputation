#!/bin/bash
#############################################
# Test script for h3abionet/chipimputation #
#############################################

set -e

echo "=================================================="
echo "H3ABionet ChipImputation Pipeline - Test Run"
echo "=================================================="
echo ""

# Set up environment
export NXF_WORK="/scratch3/users/mamana/nextflow-work"
export NXF_TEMP="/scratch3/users/mamana/nextflow-temp"

# Create directories if they don't exist
mkdir -p $NXF_WORK
mkdir -p $NXF_TEMP

# Function to run tests
run_test() {
    local test_name=$1
    local test_cmd=$2
    
    echo "-------------------------------------------"
    echo "Running: $test_name"
    echo "-------------------------------------------"
    echo "Command: $test_cmd"
    echo ""
    
    eval $test_cmd
    
    if [ $? -eq 0 ]; then
        echo "✓ $test_name completed successfully"
    else
        echo "✗ $test_name failed"
        exit 1
    fi
    echo ""
}

# Test 1: Help message
run_test "Help Message Test" \
    "nextflow run main_simple.nf --help"

# Test 2: Stub run with test config
run_test "Stub Run Test" \
    "nextflow run main_simple.nf \
        -profile test_complete \
        -stub-run \
        --outdir /scratch3/users/mamana/test_stub"

# Test 3: Real run with minimal data (if test data exists)
if [ -f "test_data/chunk_chr21_1_5000000.vcf.gz" ]; then
    run_test "Real Data Test" \
        "nextflow run main_simple.nf \
            -profile test_complete,singularity \
            --input test_sample.csv \
            --outdir /scratch3/users/mamana/test_real \
            -resume"
else
    echo "Skipping real data test - test data not found"
fi

echo "=================================================="
echo "All tests completed!"
echo "=================================================="