#!/bin/bash
# Build script for Eagle + VCF Processing Container

set -euo pipefail

echo "Building Eagle + VCF Processing Container..."

# Build the container
docker build -t mamana/eagle-vcf-processing:eagle-2.4.1 .

echo "Container built successfully!"

# Test the container
echo "Testing container..."
docker run --rm mamana/eagle-vcf-processing:eagle-2.4.1 eagle_vcf_tools_info

echo "All tests passed!"
echo "Container: mamana/eagle-vcf-processing:eagle-2.4.1" 