#!/bin/bash
#
# Validation script for refactored imputation pipeline
# Checks that all required components are in place and functional

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Counters
TESTS_PASSED=0
TESTS_FAILED=0

echo "========================================="
echo "Imputation Pipeline v2.0 Validation"
echo "========================================="
echo ""

# Function to check if file exists
check_file() {
    local file=$1
    local description=$2
    
    if [ -f "$WORKFLOW_DIR/$file" ]; then
        echo -e "${GREEN}✓${NC} $description exists"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗${NC} $description missing: $file"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Function to check if directory exists
check_dir() {
    local dir=$1
    local description=$2
    
    if [ -d "$WORKFLOW_DIR/$dir" ]; then
        echo -e "${GREEN}✓${NC} $description exists"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗${NC} $description missing: $dir"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Function to validate Nextflow syntax
check_nf_syntax() {
    local file=$1
    local description=$2
    
    if nextflow config -flat "$WORKFLOW_DIR/$file" >/dev/null 2>&1; then
        echo -e "${GREEN}✓${NC} $description syntax valid"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗${NC} $description has syntax errors"
        ((TESTS_FAILED++))
        return 1
    fi
}

echo "1. Checking directory structure..."
echo "-----------------------------------"
check_dir "subworkflows/local" "Subworkflows directory"
check_dir "modules/local" "Modules directory"
check_dir "lib" "Library directory"
check_dir "conf" "Configuration directory"
check_dir "bin" "Binary directory"
echo ""

echo "2. Checking main workflow files..."
echo "-----------------------------------"
check_file "main.nf" "Main workflow"
check_file "README.md" "Documentation"
echo ""

echo "3. Checking subworkflows..."
echo "-----------------------------------"
check_file "subworkflows/local/input_validation.nf" "Input validation subworkflow"
check_file "subworkflows/local/quality_control.nf" "Quality control subworkflow"
check_file "subworkflows/local/phasing.nf" "Phasing subworkflow"
check_file "subworkflows/local/imputation.nf" "Imputation subworkflow"
check_file "subworkflows/local/reporting.nf" "Reporting subworkflow"
echo ""

echo "4. Checking configuration files..."
echo "-----------------------------------"
check_file "conf/base.config" "Base configuration"
check_file "conf/modules.config" "Modules configuration"
check_file "conf/test/test.config" "Test configuration"
echo ""

echo "5. Checking key modules..."
echo "-----------------------------------"
check_file "modules/local/qc/check_vcf_format.nf" "VCF format checker"
check_file "modules/local/qc/remove_duplicates.nf" "Duplicate remover"
check_file "modules/local/qc/split_multiallelic.nf" "Multiallelic splitter"
check_file "modules/local/qc/filter_variants.nf" "Variant filter"
check_file "modules/local/phasing/eagle.nf" "Eagle phasing"
check_file "modules/local/imputation/minimac4.nf" "Minimac4 imputation"
check_file "modules/local/visualization/calculate_metrics.nf" "Metrics calculator"
echo ""

echo "6. Checking helper scripts..."
echo "-----------------------------------"
check_file "lib/WorkflowMain.groovy" "Main workflow library"
check_file "bin/migrate_config.py" "Migration script"
echo ""

echo "7. Validating Nextflow syntax..."
echo "-----------------------------------"
if command -v nextflow &> /dev/null; then
    check_nf_syntax "main.nf" "Main workflow"
    check_nf_syntax "conf/base.config" "Base configuration"
else
    echo -e "${YELLOW}⚠${NC} Nextflow not found - skipping syntax validation"
fi
echo ""

echo "8. Checking for required containers..."
echo "-----------------------------------"
CONTAINERS=(
    "mamana/minimac4:minimac4-4.1.6"
    "mamana/eagle-vcf-processing:eagle-2.4.1"
    "mamana/python-plotting:1.0.0"
    "mamana/vcf-processing:bcftools-1.20"
)

for container in "${CONTAINERS[@]}"; do
    echo -e "${YELLOW}ℹ${NC} Container required: $container"
done
echo ""

echo "========================================="
echo "Validation Summary"
echo "========================================="
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All validation checks passed!${NC}"
    echo ""
    echo "The pipeline is ready to run. Test with:"
    echo "  cd $WORKFLOW_DIR"
    echo "  nextflow run main.nf -profile test,docker"
    exit 0
else
    echo -e "${RED}✗ Validation failed with $TESTS_FAILED errors${NC}"
    echo ""
    echo "Please fix the issues above before running the pipeline."
    exit 1
fi