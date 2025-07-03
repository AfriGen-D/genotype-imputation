# Genotype Imputation Containers

This directory contains Docker containers optimized for genotype imputation workflows, leveraging existing containers where possible and only building custom containers for missing tools.

## 🐳 Container Architecture

### Existing Containers (Ready to Use)
- **BCFtools/SAMtools/HTSlib**: `ghcr.io/wtsi-npg/samtools:latest`
- **VCFtools**: `quay.io/biocontainers/vcftools:latest`
- **Python Environment**: `python:3.11-slim`
- **R Environment**: `bioconductor/bioconductor_docker:latest`

### Custom Containers (We Build)
- **Eagle Phasing**: `ghcr.io/afrigen-d/eagle-phasing:latest`
- **Minimac4**: `ghcr.io/afrigen-d/minimac4:latest`
- **PLINK 2.0**: `ghcr.io/afrigen-d/plink2:latest`
- **All-in-One**: `ghcr.io/afrigen-d/genotype-imputation:latest`

## 🔧 Usage Examples

### Using Existing Containers
```bash
# VCF processing with BCFtools
docker run --rm -v $(pwd):/data ghcr.io/wtsi-npg/samtools:latest \
  bcftools view -H /data/input.vcf.gz | wc -l

# VCF filtering with VCFtools
docker run --rm -v $(pwd):/data quay.io/biocontainers/vcftools:latest \
  vcftools --gzvcf /data/input.vcf.gz --maf 0.05 --recode --stdout
```

### Using Custom Containers
```bash
# Eagle phasing
docker run --rm -v $(pwd):/data ghcr.io/afrigen-d/eagle-phasing:latest \
  eagle --vcf=/data/input.vcf.gz --geneticMap=/data/genetic_map.txt

# Minimac4 imputation
docker run --rm -v $(pwd):/data ghcr.io/afrigen-d/minimac4:latest \
  minimac4 --refHaps /data/reference.m3vcf.gz --haps /data/target.vcf.gz
```

## 📊 Container Comparison

| Container | Size | Build Time | Use Case |
|-----------|------|------------|----------|
| WTSI Samtools | ~200MB | 0min (pull) | VCF processing |
| BioContainers VCFtools | ~150MB | 0min (pull) | VCF analysis |
| Custom Eagle | ~300MB | 15min | Phasing |
| Custom Minimac4 | ~250MB | 20min | Imputation |
| Custom All-in-One | ~800MB | 35min | Complete workflow |

## 🚀 Quick Start

### 1. Use Docker Compose
```bash
# Start all services
docker-compose up -d

# Run specific workflow
docker-compose run --rm eagle-phasing \
  eagle --vcf=/data/input.vcf.gz --out=/data/phased
```

### 2. Build Custom Containers Only
```bash
# Build only what we need
./build-custom.sh

# Or build specific container
docker build -t genotype-imputation/eagle-phasing eagle-phasing/
```

### 3. Building on Ubuntu VM (Recommended for x86_64)
```bash
# Clone your repository
git clone https://github.com/YourUsername/genotype-imputation.git
cd genotype-imputation/containers

# Build all custom containers
./build-custom.sh

# Or build individual containers
docker build -t eagle-phasing:latest ./eagle-phasing
docker build -t minimac4:latest ./minimac4
docker build -t plink2:latest ./plink2
docker build -t all-in-one:latest ./all-in-one

# Test with docker-compose
docker-compose --profile custom up -d
```

### 3. Use in Nextflow
```nextflow
process EAGLE_PHASING {
    container 'ghcr.io/afrigen-d/eagle-phasing:latest'
    
    input:
    path vcf
    
    output:
    path "*.phased.vcf.gz"
    
    script:
    """
    eagle --vcf=${vcf} --geneticMap=genetic_map.txt --outPrefix=phased
    """
}
```

## 🔄 Migration from Monolithic

If migrating from the original monolithic container:
1. **VCF Processing**: Replace with `ghcr.io/wtsi-npg/samtools:latest`
2. **VCF Analysis**: Replace with `quay.io/biocontainers/vcftools:latest`
3. **Phasing**: Use new `ghcr.io/afrigen-d/eagle-phasing:latest`
4. **Imputation**: Use new `ghcr.io/afrigen-d/minimac4:latest`

## 🏗️ Development

### Adding New Tools
1. Check if container already exists in:
   - [BioContainers](https://biocontainers.pro)
   - [WTSI Containers](https://github.com/wtsi-npg)
   - [Docker Hub](https://hub.docker.com)
2. If not available, create new container in appropriate directory
3. Update CI/CD pipeline to build new container

### Testing
```bash
# Test existing containers
./test-existing.sh

# Test custom containers
./test-custom.sh

# Run full integration test
./test-integration.sh
```

## 📈 Performance Benefits

Using existing containers provides:
- **Faster deployment**: No compilation time for existing tools
- **Smaller total size**: Avoid duplicate dependencies
- **Better maintenance**: Upstream updates handled by maintainers
- **Proven stability**: Widely tested containers

## 🔧 Recent Improvements

### Fixed Ubuntu Repository Issues
The containers have been updated to resolve widespread Ubuntu repository hash sum mismatches:
- **Root Cause**: Ubuntu 20.04 and 22.04 repositories experiencing hash verification failures
- **Solution**: Migrated to stable Alpine Linux 3.19 base images
- **Impact**: Containers now build reliably without package installation failures

### Architecture Compatibility
- **Target Platform**: linux/amd64 for HPC compatibility
- **Build Environment**: Optimized for Ubuntu VM deployment
- **Cross-platform**: Works on both x86_64 and ARM64 hosts (via emulation)

### Container Optimizations
- **Size Reduction**: Alpine-based containers are 60-80% smaller than Ubuntu equivalents
- **Build Time**: Faster builds due to smaller base images
- **Security**: Minimal attack surface with Alpine's security-focused design
- **Tool Versions**: Updated to latest stable versions of all bioinformatics tools

### What's Changed
```bash
# Previous (Ubuntu-based, broken)
FROM ubuntu:20.04  # ❌ Hash sum mismatch errors
RUN apt-get update && apt-get install -y bcftools  # ❌ Fails to install

# Current (Alpine-based, working)
FROM alpine:3.19   # ✅ Stable, secure, lightweight
RUN apk add --no-cache bcftools  # ✅ Reliable installation
```

### Build Script Improvements
- **Parallel Building**: Build multiple containers simultaneously
- **Error Handling**: Comprehensive error reporting and recovery
- **Platform Support**: Automatic platform detection and optimization
- **External Container Check**: Verify availability of upstream containers

## 🛠️ Troubleshooting

### Architecture Issues

**Problem**: `rosetta error: failed to open elf at /lib64/ld-linux-x86-64.so.2`
**Solution**: You're running x86_64 containers on Apple Silicon (ARM64)

```bash
# Option 1: Build for your native architecture
docker build --platform linux/arm64 -t eagle-phasing:arm64 ./eagle-phasing

# Option 2: Use Ubuntu VM for x86_64 builds (recommended)
# The containers are optimized for linux/amd64 (HPC environments)
```

### Permission Issues

**Problem**: `Permission denied` when pushing to GitHub
**Solution**: 
1. Fork the repository: https://github.com/AfriGen-D/genotype-imputation
2. Update your remote: `git remote set-url origin https://github.com/YourUsername/genotype-imputation.git`
3. Push: `git push origin master`

### Container Build Failures

**Problem**: Package installation fails with hash sum mismatch
**Solution**: The containers now use Alpine Linux to avoid Ubuntu repository issues

**Problem**: Tool not found after container build
**Solution**: 
```bash
# Check if binary exists
docker run --rm container-name which tool-name

# Check container contents
docker run --rm -it container-name sh
```

### Docker Compose Issues

**Problem**: `image not found` errors
**Solution**: 
```bash
# Build custom containers first
./build-custom.sh

# Or pull existing containers
docker-compose pull

# Use specific profile
docker-compose --profile custom up -d
```

## 🤝 Contributing

When contributing new containers:
1. First check if tool already exists in container registries
2. If building custom, optimize for size and build time
3. Use multi-stage builds where possible
4. Include comprehensive tests
5. Document usage examples

## 📚 References

- [WTSI Samtools Container](https://github.com/wtsi-npg/samtools_container)
- [BioContainers Registry](https://biocontainers.pro)
- [Eagle Phasing Documentation](https://alkesgroup.broadinstitute.org/Eagle/)
- [Minimac4 Documentation](https://genome.sph.umich.edu/wiki/Minimac4) 