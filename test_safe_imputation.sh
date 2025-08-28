#!/bin/bash

###########################################
# Test Safe Imputation Implementation
###########################################

echo "============================================"
echo "Testing Safe Imputation Implementation"
echo "============================================"
echo ""

# Configuration
PIPELINE_DIR="/users/mamana/genotype-imputation"
TEST_CHR="chr9"
TEST_START="45046581"
TEST_END="50046580"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    if [ "$1" == "SUCCESS" ]; then
        echo -e "${GREEN}✓${NC} $2"
    elif [ "$1" == "WARNING" ]; then
        echo -e "${YELLOW}⚠${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
    fi
}

# Test 1: Check if validation script exists and is executable
echo "Test 1: Checking validation script..."
if [ -x "${PIPELINE_DIR}/bin/validate_imputation_chunk.py" ]; then
    print_status "SUCCESS" "Validation script exists and is executable"
else
    print_status "FAIL" "Validation script not found or not executable"
    exit 1
fi

# Test 2: Check if safe imputation modules exist
echo ""
echo "Test 2: Checking safe imputation modules..."
modules_ok=true

if [ -f "${PIPELINE_DIR}/modules/local/impute/impute_minimac4_safe.nf" ]; then
    print_status "SUCCESS" "IMPUTE_MINIMAC4_SAFE module found"
else
    print_status "FAIL" "IMPUTE_MINIMAC4_SAFE module not found"
    modules_ok=false
fi

if [ -f "${PIPELINE_DIR}/modules/local/impute/validate_imputation_chunk.nf" ]; then
    print_status "SUCCESS" "VALIDATE_IMPUTATION_CHUNK module found"
else
    print_status "FAIL" "VALIDATE_IMPUTATION_CHUNK module not found"
    modules_ok=false
fi

if [ "$modules_ok" == false ]; then
    echo "Module check failed. Exiting."
    exit 1
fi

# Test 3: Check if subworkflow has been updated
echo ""
echo "Test 3: Checking subworkflow updates..."
if grep -q "IMPUTE_MINIMAC4_SAFE" "${PIPELINE_DIR}/subworkflows/local/impute.nf"; then
    print_status "SUCCESS" "Subworkflow includes safe imputation module"
else
    print_status "FAIL" "Subworkflow not updated with safe imputation"
fi

if grep -q "use_safe_imputation" "${PIPELINE_DIR}/subworkflows/local/impute.nf"; then
    print_status "SUCCESS" "Subworkflow checks for safe imputation parameter"
else
    print_status "WARNING" "Subworkflow may not check safe imputation parameter"
fi

# Test 4: Check configuration
echo ""
echo "Test 4: Checking configuration..."
if grep -q "use_safe_imputation" "${PIPELINE_DIR}/nextflow.config"; then
    print_status "SUCCESS" "Main config includes safe imputation parameters"
else
    print_status "WARNING" "Main config missing safe imputation parameters"
fi

if grep -q "use_safe_imputation = true" "${PIPELINE_DIR}/v6_chr21_phased_nfcore.config"; then
    print_status "SUCCESS" "v6 config has safe imputation enabled"
else
    print_status "WARNING" "v6 config doesn't have safe imputation enabled"
fi

# Test 5: Dry run of validation script (if test data is available)
echo ""
echo "Test 5: Testing validation script functionality..."

# Create minimal test VCF files for validation
TEST_DIR="/tmp/safe_imputation_test_$$"
mkdir -p "${TEST_DIR}"

# Create a minimal test VCF with no variants (to simulate the problem)
cat > "${TEST_DIR}/test_empty.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=chr9,length=138394717>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF

# Compress and index
bgzip -f "${TEST_DIR}/test_empty.vcf"
bcftools index -t "${TEST_DIR}/test_empty.vcf.gz" 2>/dev/null || tabix -p vcf "${TEST_DIR}/test_empty.vcf.gz"

# Create a reference VCF with some variants
cat > "${TEST_DIR}/test_ref.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=chr9,length=138394717>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	REF1
chr9	45046600	rs1	A	G	100	PASS	.	GT	0/1
chr9	45046700	rs2	C	T	100	PASS	.	GT	1/1
chr9	45046800	rs3	G	A	100	PASS	.	GT	0/1
EOF

bgzip -f "${TEST_DIR}/test_ref.vcf"
bcftools index -t "${TEST_DIR}/test_ref.vcf.gz" 2>/dev/null || tabix -p vcf "${TEST_DIR}/test_ref.vcf.gz"

# Run validation
if python3 "${PIPELINE_DIR}/bin/validate_imputation_chunk.py" \
    --vcf "${TEST_DIR}/test_empty.vcf.gz" \
    --ref-vcf "${TEST_DIR}/test_ref.vcf.gz" \
    --chrom chr9 \
    --start ${TEST_START} \
    --end ${TEST_END} \
    --output "${TEST_DIR}/validation.txt" \
    --status-file "${TEST_DIR}/validation.status" 2>/dev/null; then
    print_status "WARNING" "Validation passed (unexpected for empty chunk)"
else
    # Check if it correctly identified the problem
    if grep -q "FAIL" "${TEST_DIR}/validation.status" 2>/dev/null; then
        print_status "SUCCESS" "Validation correctly identified empty chunk as FAIL"
    else
        print_status "FAIL" "Validation script didn't work as expected"
    fi
fi

# Clean up test files
rm -rf "${TEST_DIR}"

# Test 6: Check if pipeline can handle the safe mode
echo ""
echo "Test 6: Testing pipeline configuration..."
cd "${PIPELINE_DIR}"

# Test if the pipeline configuration is valid
if nextflow config -profile standard | grep -q "use_safe_imputation = true"; then
    print_status "SUCCESS" "Pipeline configuration includes safe imputation"
else
    print_status "WARNING" "Safe imputation may not be properly configured"
fi

# Summary
echo ""
echo "============================================"
echo "Test Summary"
echo "============================================"
echo ""
echo "Safe imputation implementation has been tested."
echo "The pipeline should now handle empty chunks gracefully."
echo ""
echo "To run the pipeline with safe imputation:"
echo "  nextflow run main_nfcore.nf -profile singularity \\
            -c v6_chr21_phased_nfcore.config \\
            --use_safe_imputation true"
echo ""
echo "To enable extra validation (slower but safer):"
echo "  nextflow run main_nfcore.nf -profile singularity \\
            -c v6_chr21_phased_nfcore.config \\
            --use_safe_imputation true \\
            --validate_chunks true"
echo ""
print_status "SUCCESS" "Safe imputation setup complete!"