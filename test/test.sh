#!/bin/bash
###############################################
# Test script for ChipImputation Pipeline    #
###############################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}==================================================${NC}"
echo -e "${GREEN}ChipImputation Pipeline - Test Run${NC}"
echo -e "${GREEN}==================================================${NC}"

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if test files exist in test directory
if [ ! -f "$SCRIPT_DIR/test_sample.csv" ]; then
    echo -e "${RED}Error: Test sample file not found at $SCRIPT_DIR/test_sample.csv${NC}"
    echo "Please ensure test files are in the test/ directory"
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/v6_chr21_phased_nfcore.config" ]; then
    echo -e "${RED}Error: Test config not found at $SCRIPT_DIR/v6_chr21_phased_nfcore.config${NC}"
    exit 1
fi

# Set test parameters (using test directory files)
INPUT="$SCRIPT_DIR/test_sample.csv"
CONFIG="$SCRIPT_DIR/v6_chr21_phased_nfcore.config"
OUTDIR="/scratch3/users/mamana/results/test_$(date +%Y%m%d_%H%M%S)"
PROFILE="singularity,slurm"

# Parse command line options
while [[ $# -gt 0 ]]; do
    case $1 in
        --stub)
            STUB="-stub-run"
            echo -e "${YELLOW}Running in stub mode (no actual processing)${NC}"
            shift
            ;;
        --resume)
            RESUME="-resume"
            echo -e "${YELLOW}Resuming from previous run${NC}"
            shift
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --stub           Run in stub mode (for testing)"
            echo "  --resume         Resume from previous run"
            echo "  --profile PROF   Execution profile (default: singularity,slurm)"
            echo "  --outdir DIR     Output directory (default: /scratch3/users/mamana/results/test_TIMESTAMP)"
            echo "  --help           Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Quick test with stub mode"
            echo "  test/test.sh --stub"
            echo ""
            echo "  # Full test run"
            echo "  test/test.sh"
            echo ""
            echo "  # Resume a failed test"
            echo "  test/test.sh --resume"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Run './test.sh --help' for usage"
            exit 1
            ;;
    esac
done

# Display test configuration
echo ""
echo "Test Configuration:"
echo "-------------------"
echo "Input:    $INPUT"
echo "Config:   $CONFIG"
echo "Output:   $OUTDIR"
echo "Profile:  $PROFILE"
echo ""

# Confirm before running
read -p "Proceed with test? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Test cancelled${NC}"
    exit 0
fi

# Run the test (from parent directory)
echo ""
echo -e "${GREEN}Starting test run...${NC}"
echo ""

# Change to parent directory to run pipeline
cd "$SCRIPT_DIR/.."

nextflow run main_nfcore.nf \
    -c "$CONFIG" \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    -profile "$PROFILE" \
    $RESUME \
    $STUB \
    -with-report "${OUTDIR}/pipeline_info/execution_report.html" \
    -with-timeline "${OUTDIR}/pipeline_info/execution_timeline.html" \
    -with-trace "${OUTDIR}/pipeline_info/execution_trace.txt"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}==================================================${NC}"
    echo -e "${GREEN}Test completed successfully!${NC}"
    echo -e "${GREEN}Results: $OUTDIR${NC}"
    echo -e "${GREEN}==================================================${NC}"
else
    echo ""
    echo -e "${RED}==================================================${NC}"
    echo -e "${RED}Test failed. Check the logs for details.${NC}"
    echo -e "${RED}==================================================${NC}"
    exit 1
fi