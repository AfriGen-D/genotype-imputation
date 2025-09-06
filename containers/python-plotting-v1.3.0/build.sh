#!/bin/bash

# Build script for python-plotting container v1.3.0
# Usage: ./build.sh [--push]

set -e

# Configuration
IMAGE_NAME="mamana/python-plotting"
VERSION="1.3.0"
FULL_TAG="${IMAGE_NAME}:${VERSION}"
LATEST_TAG="${IMAGE_NAME}:latest"

echo "Building Python Plotting Container v${VERSION}"
echo "=========================================="

# Build the Docker image
echo "Building Docker image: ${FULL_TAG}"
docker build -t ${FULL_TAG} .

# Tag as latest
echo "Tagging as latest..."
docker tag ${FULL_TAG} ${LATEST_TAG}

# Push to Docker Hub if requested
if [ "$1" == "--push" ]; then
    echo "Pushing to Docker Hub..."
    docker push ${FULL_TAG}
    docker push ${LATEST_TAG}
    echo "Successfully pushed ${FULL_TAG} and ${LATEST_TAG}"
else
    echo "Build complete. To push to Docker Hub, run: ./build.sh --push"
fi

echo "Done!"