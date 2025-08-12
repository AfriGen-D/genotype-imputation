#!/bin/bash

# Build script for Python plotting container
# Usage: ./build.sh [tag_name]

set -e  # Exit on error

# Set default tag if not provided
TAG=${1:-"mamana/python-plotting:1.0.0"}
LATEST_TAG="mamana/python-plotting:latest"

echo "Building Python plotting container..."
echo "Tag: $TAG"

# Build the Docker image
docker build -t "$TAG" -t "$LATEST_TAG" .

if [ $? -eq 0 ]; then
    echo "✓ Successfully built Docker image: $TAG"
    echo "✓ Also tagged as: $LATEST_TAG"
    
    # Display image info
    echo ""
    echo "Image information:"
    docker images | grep "mamana/python-plotting"
    
    # Test the container
    echo ""
    echo "Testing container..."
    docker run --rm "$TAG" python -c "
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
print('✓ NumPy version:', np.__version__)
print('✓ Pandas version:', pd.__version__)
print('✓ Matplotlib version:', matplotlib.__version__)
print('✓ Seaborn version:', sns.__version__)
print('✓ All packages imported successfully!')
"
    
    echo ""
    echo "Container is ready to use!"
    echo ""
    echo "To push to registry:"
    echo "  docker push $TAG"
    echo "  docker push $LATEST_TAG"
    
else
    echo "✗ Failed to build Docker image"
    exit 1
fi