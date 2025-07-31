# Minimac4-All Container Build Instructions

## Overview

This directory contains the `minimac4-all` container - a comprehensive genotype imputation environment optimized for Azure deployments. We've created two versions to handle different scenarios:

### Available Dockerfiles

1. **`Dockerfile`** (Debian-based)
   - **Base**: `debian:bookworm-slim`
   - **Size**: ~2.5GB 
   - **Build Time**: 30-45 minutes
   - **Use Case**: Standard production environments, full compatibility

2. **`Dockerfile.alpine`** (Alpine-based) ⭐ **RECOMMENDED**
   - **Base**: `alpine:3.19`
   - **Size**: ~1.8GB
   - **Build Time**: 20-30 minutes
   - **Use Case**: Azure deployments, lightweight & fast

## Quick Start

### Option 1: Using Build Script (Recommended)

```bash
# Navigate to containers directory
cd /Users/mamana/projects/genotype-imputation/containers

# Build using the custom build script (includes minimac4-all)
./build-custom.sh

# Or build with specific tag
./build-custom.sh v2.0
```

### Option 2: Direct Docker Build

```bash
# Alpine version (recommended for Azure)
docker build --platform linux/amd64 \
  -t ghcr.io/afrigen-d/minimac4-all:alpine \
  -f ./minimac4-all/Dockerfile.alpine \
  ./minimac4-all/

# Debian version (if you need full compatibility)
docker build --platform linux/amd64 \
  -t ghcr.io/afrigen-d/minimac4-all:latest \
  ./minimac4-all/
```

## Included Tools

The `minimac4-all` container includes:

### Core Imputation Tools
- **Minimac4 v4.1.6** - Primary genotype imputation (built from source)
- **Eagle v2.4.1** - Haplotype phasing
- **PLINK2** - Genetic data analysis (Alpine version only)

### Genomics Utilities
- **SAMtools v1.19** - Sequence alignment/map manipulation
- **BCFtools v1.19** - Variant calling and manipulation  
- **HTSlib v1.19** - High-throughput sequencing data processing
- **VCFtools v0.1.16** - VCF file manipulation and analysis
- **BEDtools v2.31.1** - Genomic interval operations

### Data Analysis Environment
- **Python 3** with scientific packages (NumPy, Pandas, SciPy, Matplotlib, Seaborn)
- **Bash** shell environment
- **Non-root user** (`impuser`) for security

## Testing the Container

```bash
# Quick test
docker run --rm ghcr.io/afrigen-d/minimac4-all:alpine test-minimac4-all

# Interactive session
docker run -it --rm \
  -v $(pwd)/data:/data \
  ghcr.io/afrigen-d/minimac4-all:alpine

# Test specific tool
docker run --rm ghcr.io/afrigen-d/minimac4-all:alpine minimac4 --help
```

## Azure Deployment Examples

### Azure Container Instances

```bash
az container create \
  --resource-group myResourceGroup \
  --name minimac4-imputation \
  --image ghcr.io/afrigen-d/minimac4-all:alpine \
  --cpu 4 \
  --memory 16 \
  --environment-variables TZ=UTC \
  --azure-file-volume-account-name mystorageaccount \
  --azure-file-volume-account-key $STORAGE_KEY \
  --azure-file-volume-share-name data \
  --azure-file-volume-mount-path /data
```

### Azure Batch

```yaml
pool:
  id: minimac4-pool
  vmSize: Standard_D4s_v3
  containerConfiguration:
    containerImageNames:
      - ghcr.io/afrigen-d/minimac4-all:alpine
```

## Build Status Updates

### Recent Improvements (Latest Version)
- ✅ Fixed Ubuntu repository hash sum mismatch issues
- ✅ Added Alpine Linux version for better reliability
- ✅ Optimized build process and container size
- ✅ Added comprehensive testing scripts
- ✅ Improved Azure compatibility
- ✅ Built all tools from source for maximum compatibility

### Known Issues & Solutions

#### Issue: Ubuntu/Debian Repository Errors
**Solution**: Use the Alpine version (`Dockerfile.alpine`)

```bash
# Instead of this (problematic)
docker build -f Dockerfile ./minimac4-all/

# Use this (reliable)
docker build -f Dockerfile.alpine ./minimac4-all/
```

#### Issue: Build Takes Too Long
**Solution**: Use the Alpine version and enable BuildKit

```bash
export DOCKER_BUILDKIT=1
docker build -f Dockerfile.alpine ./minimac4-all/
```

#### Issue: Platform Architecture Warnings
**Solution**: Explicitly specify platform for Azure compatibility

```bash
docker build --platform linux/amd64 -f Dockerfile.alpine ./minimac4-all/
```

## Container Registry

Once built, containers are tagged as:
- `ghcr.io/afrigen-d/minimac4-all:latest` (Debian version)
- `ghcr.io/afrigen-d/minimac4-all:alpine` (Alpine version)

## Performance Recommendations

### For Azure Deployments
- Use the **Alpine version** for faster startup and smaller storage requirements
- Allocate **16GB+ RAM** for large imputation jobs
- Use **Premium SSD storage** for `/data/temp` directory
- Consider **proximity placement groups** for multi-node jobs

### For Development/Testing
- Mount local directories to `/data/input`, `/data/output`, `/data/reference`
- Use `test-minimac4-all` command to verify all tools are working
- Check container logs if builds fail

## Next Steps

1. **Build the container** using the Alpine version
2. **Test locally** with your data
3. **Deploy to Azure** using the provided examples
4. **Monitor performance** and adjust resources as needed

## Support

For container-specific issues:
1. Run the test suite: `docker run --rm <image> test-minimac4-all`
2. Check tool versions and availability
3. Verify Azure storage mount points
4. Review build logs for any compilation errors

The Alpine version is recommended for most Azure use cases due to its reliability, smaller size, and faster build times. 