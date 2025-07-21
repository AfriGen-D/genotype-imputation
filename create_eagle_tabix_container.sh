#!/bin/bash
# Script to create a container with both Eagle and tabix by copying eagle from eagle container to minimac4 container

set -euo pipefail

echo "Creating combined Eagle + tabix container..."

# Copy the minimac4 container (which has tabix) to create our new container
cp singularity_cache/mamana-imputation-minimac4-4.1.6.img singularity_cache/mamana-eagle-tabix-combined.img

# Extract eagle binary from eagle container and copy to new container
echo "Extracting eagle binary from eagle container..."
singularity exec singularity_cache/mamana-eagle-phasing-eagle-2.4.1.img cp /usr/local/bin/eagle /tmp/eagle

# Copy eagle binary to the new container
echo "Copying eagle binary to combined container..."
singularity exec --writable singularity_cache/mamana-eagle-tabix-combined.img cp /tmp/eagle /usr/local/bin/eagle
singularity exec --writable singularity_cache/mamana-eagle-tabix-combined.img chmod +x /usr/local/bin/eagle

# Clean up
rm -f /tmp/eagle

# Test the combined container
echo "Testing combined container..."
echo "Eagle version:"
singularity exec singularity_cache/mamana-eagle-tabix-combined.img eagle --version
echo "Tabix version:"
singularity exec singularity_cache/mamana-eagle-tabix-combined.img tabix --version

echo "Combined container created successfully: singularity_cache/mamana-eagle-tabix-combined.img" 