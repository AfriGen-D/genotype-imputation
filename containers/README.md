# Container Architecture for Genotype Imputation Pipeline

## Overview

This pipeline uses a modular container architecture with specialized containers for different tools and workflows. The containers are available on Docker Hub for easy access and optimal performance.

## Container Registry

All containers are available on Docker Hub under the `mamana` organization:
- **Registry**: `docker.io/mamana/`
- **Alternative**: `ghcr.io/afrigen-d/` (GitHub Container Registry)

## Container Catalog

### 🔧 Core Tool Containers

#### 1. Minimac4 Imputation
```bash
docker pull mamana/minimac4:minimac4-4.1.6
```
- **Purpose**: Genotype imputation using Minimac4
- **Tools**: Minimac4 4.1.6
- **Base**: Alpine Linux 3.18
- **Size**: ~300MB

#### 2. Eagle Phasing
```bash
docker pull mamana/eagle-phasing:eagle-2.4.1
```
- **Purpose**: Haplotype phasing using Eagle
- **Tools**: Eagle 2.4.1, HTSlib 1.20
- **Base**: Alpine Linux 3.18
- **Size**: ~280MB

#### 3. VCF Processing
```bash
docker pull mamana/vcf-processing:bcftools-1.20
```
- **Purpose**: VCF file processing and quality control
- **Tools**: BCFtools 1.20, VCFtools 0.1.16, tabix, HTSlib 1.20
- **Base**: Alpine Linux 3.18
- **Size**: ~250MB

### 🚀 Workflow Containers

#### 4. Comprehensive Imputation
```bash
docker pull mamana/imputation:minimac4-4.1.6
```
- **Purpose**: Complete imputation workflow
- **Tools**: Minimac4 4.1.6, BCFtools 1.20, VCFtools 0.1.16, tabix
- **Base**: Alpine Linux 3.18
- **Size**: ~450MB

#### 5. Comprehensive Phasing
```bash
docker pull mamana/phasing:eagle-2.4.1
```
- **Purpose**: Complete phasing workflow
- **Tools**: Eagle 2.4.1, HTSlib 1.20, tabix
- **Base**: Alpine Linux 3.18
- **Size**: ~350MB

## Usage in Nextflow

### Process-Specific Container Assignment

The pipeline automatically assigns the appropriate container based on the process:

```nextflow
// Minimac4 imputation processes
withName: 'impute_minimac4' {
  container = 'mamana/minimac4:minimac4-4.1.6'
}

// Eagle phasing processes
withName: 'minimac4_phasing_eagle' {
  container = 'mamana/eagle-phasing:eagle-2.4.1'
}

// VCF processing and QC
withName: 'check_chromosome' {
  container = 'mamana/vcf-processing:bcftools-1.20'
}
```

### Docker Compose Profiles

Use Docker Compose for local development and testing:

```bash
# Production containers (optimized, lightweight)
docker-compose --profile production up

# All containers (including legacy)
docker-compose --profile all up

# Specific service
docker-compose up imputation
```

## Performance Optimizations

### 🏔️ Alpine Linux Base
- **Smaller size**: 40-60% smaller than Ubuntu-based containers
- **Faster startup**: Reduced image pull and container start time
- **Security**: Minimal attack surface with musl libc

### 🔧 Tool-Specific Containers
- **Targeted**: Only includes necessary tools for specific processes
- **Efficient**: No unnecessary dependencies
- **Reliable**: Fewer conflicts between tools

### 🚀 Multi-Stage Builds
- **Optimized**: Only runtime dependencies in final image
- **Smaller**: Build tools removed from final container
- **Cached**: Intermediate stages cached for faster rebuilds

## Container Mapping

| Process Category | Container | Tools Included |
|------------------|-----------|----------------|
| **Imputation** | `mamana/minimac4:minimac4-4.1.6` | Minimac4 4.1.6 |
| **Phasing** | `mamana/eagle-phasing:eagle-2.4.1` | Eagle 2.4.1, HTSlib 1.20 |
| **VCF Processing** | `mamana/vcf-processing:bcftools-1.20` | BCFtools 1.20, VCFtools 0.1.16, tabix |
| **Quality Control** | `mamana/vcf-processing:bcftools-1.20` | BCFtools 1.20, VCFtools 0.1.16 |
| **File Operations** | `mamana/vcf-processing:bcftools-1.20` | BCFtools 1.20, tabix |
| **Comprehensive** | `mamana/imputation:minimac4-4.1.6` | All imputation tools |

## Running the Pipeline

### Option 1: Docker Profile
```bash
nextflow run h3abionet/chipimputation -profile docker
```

### Option 2: Singularity Profile
```bash
nextflow run h3abionet/chipimputation -profile singularity
```

### Option 3: Custom Container
```bash
nextflow run h3abionet/chipimputation \
  --container mamana/imputation:minimac4-4.1.6 \
  -profile docker
```

## Migration from Legacy Containers

### Old Container References
```bash
# OLD (deprecated)
quay.io/h3abionet_org/imputation_tools
quay.io/mypandos/impute5:latest
```

### New Container References
```bash
# NEW (recommended)
mamana/imputation:minimac4-4.1.6          # Comprehensive imputation
mamana/minimac4:minimac4-4.1.6           # Minimac4 only
mamana/eagle-phasing:eagle-2.4.1         # Eagle phasing
mamana/vcf-processing:bcftools-1.20      # VCF processing
mamana/phasing:eagle-2.4.1              # Comprehensive phasing
```

## Advanced Usage

### Manual Container Pull
```bash
# Pull all containers
docker pull mamana/minimac4:minimac4-4.1.6
docker pull mamana/eagle-phasing:eagle-2.4.1
docker pull mamana/imputation:minimac4-4.1.6
docker pull mamana/phasing:eagle-2.4.1
docker pull mamana/vcf-processing:bcftools-1.20
```

### Container Testing
```bash
# Test Minimac4 container
docker run --rm mamana/minimac4:minimac4-4.1.6 minimac4 --help

# Test Eagle container
docker run --rm mamana/eagle-phasing:eagle-2.4.1 eagle --help

# Test VCF processing
docker run --rm mamana/vcf-processing:bcftools-1.20 bcftools --version
```

### Development with Docker Compose
```bash
# Start development environment
docker-compose --profile production up -d

# Run analysis interactively
docker-compose exec imputation bash

# View logs
docker-compose logs -f imputation
```

## Troubleshooting

### Common Issues

1. **Container Not Found**
   ```bash
   # Ensure you're using the correct registry
   docker pull mamana/minimac4:minimac4-4.1.6
   ```

2. **Permission Issues**
   ```bash
   # Run with proper user mapping
   docker run --rm -u $(id -u):$(id -g) mamana/minimac4:minimac4-4.1.6
   ```

3. **Memory Issues**
   ```bash
   # Increase Docker memory limit
   docker run --rm -m 8g mamana/imputation:minimac4-4.1.6
   ```

### Performance Tips

1. **Use SSD storage** for Docker volumes
2. **Increase Docker memory** allocation (8GB+ recommended)
3. **Use local registry** for frequently used containers
4. **Enable BuildKit** for faster builds:
   ```bash
   export DOCKER_BUILDKIT=1
   ```

## Container Maintenance

### Updating Containers
```bash
# Pull latest versions
docker pull mamana/minimac4:minimac4-4.1.6
docker pull mamana/eagle-phasing:eagle-2.4.1
docker pull mamana/vcf-processing:bcftools-1.20

# Clean up old images
docker image prune -f
```

### Security Scanning
```bash
# Scan for vulnerabilities
docker scan mamana/minimac4:minimac4-4.1.6
```

## Support

For container-related issues:
- Check the [GitHub Issues](https://github.com/h3abionet/chipimputation/issues)
- Review the [Docker Hub repositories](https://hub.docker.com/u/mamana)
- Contact the maintainers via [AfriGen-D](mailto:info@afrigen.org)

---

**Note**: This container architecture is optimized for HPC environments and standard x86_64 systems. ARM64 support is not currently available. 