#!/bin/bash

# Multi-platform build script for minimac4-all
# GA4GH Hackathon 2025 - African Genomics Team

set -e

echo "🔧 Building minimac4-all with multi-platform support..."

# Step 1: Set up Docker buildx for multi-platform builds
echo "📦 Setting up Docker buildx..."
docker buildx create --name multiplatform-builder --use 2>/dev/null || true
docker buildx inspect --bootstrap

# Step 2: Build multi-platform image
echo "🏗️ Building for linux/amd64 and linux/arm64..."
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    --tag ghcr.io/afrigen-d/minimac4-all:multiplatform \
    --tag ghcr.io/afrigen-d/minimac4-all:latest-local \
    --load \
    .

echo "✅ Multi-platform build completed!"

# Step 3: Test the image on current platform
echo "🧪 Testing the image on current platform ($(uname -m))..."
docker run --rm ghcr.io/afrigen-d/minimac4-all:latest-local test-azure-tools

echo "🎉 Multi-platform minimac4-all container ready!"
echo "📋 Supports:"
echo "   ✅ linux/amd64 (Intel/AMD x86_64)"
echo "   ✅ linux/arm64 (Apple Silicon M1/M2/M3)"
echo ""
echo "🚀 Usage:"
echo "   docker run --rm ghcr.io/afrigen-d/minimac4-all:latest-local test-azure-tools"
echo "   # No more --platform flags needed!"

echo ""
echo "💡 To push to registry:"
echo "   docker buildx build --platform linux/amd64,linux/arm64 --push --tag your-registry/minimac4-all:latest ." 