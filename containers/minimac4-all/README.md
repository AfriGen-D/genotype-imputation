# Minimac4-All Container

A comprehensive genotype imputation container optimized for Azure environments with Minimac4 and all related genomics tools.

## Overview

This container provides a complete environment for genotype imputation workflows, including:

- **Minimac4 v4.1.6** - Primary genotype imputation software (built from source)
- **SAMtools v1.19** - Sequence alignment/map manipulation
- **BCFtools v1.19** - Variant calling and manipulation
- **HTSlib v1.19** - High-throughput sequencing data processing
- **VCFtools v0.1.16** - VCF file manipulation and analysis
- **BEDtools v2.31.1** - Genomic interval operations
- **PLINK2** - Genetic data analysis and manipulation
- **Python 3** with data science packages (NumPy, Pandas, SciPy, Matplotlib, Seaborn)

## Key Features

- **Azure Optimized**: Built specifically for Azure cloud environments
- **Source Compiled**: All tools built from source for maximum compatibility
- **Comprehensive**: All necessary tools in one container
- **Security**: Runs as non-root user (`impuser`)
- **Well-Organized**: Structured data directories (`/data/input`, `/data/output`, `/data/reference`)

## Quick Start

### Build the Container

```bash
# Build using the custom build script
./build-custom.sh

# Or build manually
docker build --platform linux/amd64 -t ghcr.io/afrigen-d/minimac4-all:latest ./minimac4-all/
```

### Run the Container

```bash
# Interactive session
docker run -it --rm \
  -v $(pwd)/data:/data \
  ghcr.io/afrigen-d/minimac4-all:latest

# Run with specific command
docker run --rm \
  -v $(pwd)/input:/data/input \
  -v $(pwd)/output:/data/output \
  -v $(pwd)/reference:/data/reference \
  ghcr.io/afrigen-d/minimac4-all:latest \
  minimac4 --help
```

### Test the Installation

```bash
# Run comprehensive test suite
docker run --rm ghcr.io/afrigen-d/minimac4-all:latest test-minimac4-all
```

## Directory Structure

- `/data/input` - Input genomic data files
- `/data/output` - Output imputation results
- `/data/reference` - Reference panels and genetic maps
- `/data/temp` - Temporary working files
- `/scripts` - Custom analysis scripts
- `/logs` - Log files and reports

## Usage Examples

### Basic Minimac4 Imputation

```bash
docker run --rm \
  -v $(pwd)/input:/data/input \
  -v $(pwd)/output:/data/output \
  -v $(pwd)/reference:/data/reference \
  ghcr.io/afrigen-d/minimac4-all:latest \
  minimac4 \
    --refHaps /data/reference/panel.m3vcf.gz \
    --haps /data/input/target.vcf.gz \
    --prefix /data/output/imputed
```

### VCF Processing Pipeline

```bash
# Inside the container
bcftools view -q 0.01 input.vcf.gz | \
bcftools norm -d both | \
bcftools filter -i 'INFO/R2>0.8' > filtered.vcf

# Convert to PLINK format
plink2 --vcf filtered.vcf --make-bed --out dataset
```

### Quality Control

```bash
# Calculate statistics
vcftools --vcf input.vcf --freq --out frequency_stats
vcftools --vcf input.vcf --missing-indv --out missing_stats

# Create BED intervals
bedtools intersect -a variants.vcf -b regions.bed > filtered_variants.vcf
```

## Azure Deployment

### Azure Container Instances (ACI)

```bash
az container create \
  --resource-group myResourceGroup \
  --name minimac4-imputation \
  --image ghcr.io/afrigen-d/minimac4-all:latest \
  --cpu 4 \
  --memory 16 \
  --azure-file-volume-account-name mystorageaccount \
  --azure-file-volume-account-key $STORAGE_KEY \
  --azure-file-volume-share-name data \
  --azure-file-volume-mount-path /data
```

### Azure Batch

```yaml
# batch-pool.yaml
pool:
  id: minimac4-pool
  vmSize: Standard_D4s_v3
  containerConfiguration:
    containerImageNames:
      - ghcr.io/afrigen-d/minimac4-all:latest
```

### Azure Kubernetes Service (AKS)

```yaml
# minimac4-job.yaml
apiVersion: batch/v1
kind: Job
metadata:
  name: minimac4-imputation
spec:
  template:
    spec:
      containers:
      - name: minimac4
        image: ghcr.io/afrigen-d/minimac4-all:latest
        command: ["minimac4"]
        args: ["--refHaps", "/data/ref/panel.m3vcf.gz", "--haps", "/data/input/target.vcf.gz", "--prefix", "/data/output/imputed"]
        volumeMounts:
        - name: data-volume
          mountPath: /data
      restartPolicy: Never
      volumes:
      - name: data-volume
        persistentVolumeClaim:
          claimName: genotype-data-pvc
```

## Environment Variables

- `PATH` - Includes `/usr/local/bin` for all compiled tools
- `LD_LIBRARY_PATH` - Includes `/usr/local/lib` for shared libraries
- `DEBIAN_FRONTEND=noninteractive` - Prevents interactive prompts
- `TZ=UTC` - Sets timezone for consistent logging

## Troubleshooting

### Permission Issues
```bash
# Fix ownership if running as root
docker run --rm -v $(pwd)/data:/data ghcr.io/afrigen-d/minimac4-all:latest \
  sudo chown -R impuser:impuser /data
```

### Memory Issues
```bash
# Increase container memory
docker run --memory=16g --rm ghcr.io/afrigen-d/minimac4-all:latest
```

### Tool Not Found
```bash
# Verify installation
docker run --rm ghcr.io/afrigen-d/minimac4-all:latest which minimac4
docker run --rm ghcr.io/afrigen-d/minimac4-all:latest test-minimac4-all
```

## Performance Optimization

### For Large Datasets
- Use SSD storage for `/data/temp`
- Allocate sufficient memory (16GB+ recommended)
- Use multiple CPU cores with `--cpus` parameter

### Azure Specific
- Use Premium Storage for reference panels
- Enable accelerated networking on VMs
- Consider proximity placement groups for multi-node jobs

## Container Size & Build Time

- **Size**: ~2.5GB (compressed)
- **Build Time**: ~30-45 minutes (depending on CPU cores)
- **Base**: Ubuntu 20.04 LTS (stable, Azure-compatible)

## Support

For issues specific to this container:
1. Run the test suite: `test-minimac4-all`
2. Check tool versions and installations
3. Verify input/output directory permissions
4. Check Azure storage mount points

For tool-specific issues, refer to the respective documentation:
- [Minimac4 Documentation](https://genome.sph.umich.edu/wiki/Minimac4)
- [SAMtools Documentation](http://www.htslib.org/)
- [BCFtools Documentation](http://samtools.github.io/bcftools/) 