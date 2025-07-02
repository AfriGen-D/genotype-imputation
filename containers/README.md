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