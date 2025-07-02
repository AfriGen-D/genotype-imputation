# Genotype Imputation Containers

This directory contains a modular container architecture for the chipimputation workflow, optimized based on actual tool usage patterns in Nextflow processes.

## 🏗️ Container Architecture

### **Process-Based Tool Groupings**

Our containers are organized by **how tools are actually used together** in workflow processes:

```
containers/
├── vcf-processing/     # BCFtools + tabix + VCFtools
├── phasing/           # Eagle + tabix  
├── imputation/        # Minimac4 + BCFtools + VCFtools
├── analysis/          # BCFtools + Python + R
├── all-in-one/        # Complete environment (original Dockerfile)
├── legacy/            # Original containers (crossmap, impute5, etc.)
└── docker-compose.yml # Container orchestration
```

---

## 📋 Container Details

### **1. VCF Processing Container** (`vcf-processing/`)
**Tools:** BCFtools 1.22, HTSlib 1.22, VCFtools 0.1.16, tabix  
**Used in processes:**
- `qc_dupl` - Remove duplicate variants
- `split_multi_allelic` - Split multi-allelic sites
- `fill_tags_vcf` - Add allele frequency tags
- `filter_min_ac` - Filter by allele count
- `sites_only` - Extract sites-only VCF
- `combine_vcfs` - Concatenate VCF files

**Usage:**
```bash
# Build container
docker build -t afrigen-d/vcf-processing ./vcf-processing

# Run interactively
docker run -it -v $(pwd)/data:/data afrigen-d/vcf-processing

# Example process
bcftools norm --rm-dup both input.vcf.gz -Ob -o output.bcf
tabix output.bcf
```

### **2. Phasing Container** (`phasing/`)
**Tools:** Eagle 2.4.1, tabix  
**Used in processes:**
- `minimac4_phasing_eagle` - Phase with reference panel
- `impute5_phasing_eagle` - Phase for IMPUTE5
- `phasing_vcf_no_ref_chunk` - Phase without reference

**Usage:**
```bash
# Build container
docker build -t afrigen-d/phasing ./phasing

# Example process
eagle --vcfTarget=target.vcf.gz --vcfRef=ref.vcf.gz --geneticMapFile=map.txt --outPrefix=phased
```

### **3. Imputation Container** (`imputation/`)
**Tools:** Minimac4 1.0.2, BCFtools 1.22, VCFtools 0.1.16  
**Used in processes:**
- `impute_minimac4` - Main imputation process
- `impute_minimac4_1` - Alternative imputation process
- `combineImpute` - Combine imputed chunks

**Usage:**
```bash
# Build container  
docker build -t afrigen-d/imputation ./imputation

# Example process
minimac4 --refHaps ref.m3vcf.gz --haps phased.vcf.gz --format GT,DS --prefix imputed
```

### **4. Analysis Container** (`analysis/`)
**Tools:** BCFtools 1.22, Python 3, R, statistical packages  
**Used in processes:**
- `generate_frequency` - Calculate allele frequencies
- `plot_freq_comparison` - Compare frequencies
- `filter_info_by_target` - Filter imputation info
- All reporting and plotting processes

**Usage:**
```bash
# Build container
docker build -t afrigen-d/analysis ./analysis

# Example process
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' imputed.vcf.gz > frequencies.txt
Rscript plot_results.R
```

### **5. All-in-One Container** (`all-in-one/`)
**Tools:** All tools in one comprehensive environment  
**Used for:** Complete workflows, testing, development

**Usage:**
```bash
# Build container
docker build -t afrigen-d/genotype-imputation ./all-in-one

# Run Nextflow pipeline
docker run -it -v $(pwd):/workspace afrigen-d/genotype-imputation
```

---

## 🚀 Quick Start

### **Option 1: Use Pre-built Containers from GitHub Actions (Recommended)**
```bash
# Pull pre-built containers from GitHub Container Registry
docker pull ghcr.io/afrigen-d/vcf-processing:latest
docker pull ghcr.io/afrigen-d/phasing:latest
docker pull ghcr.io/afrigen-d/imputation:latest
docker pull ghcr.io/afrigen-d/analysis:latest
docker pull ghcr.io/afrigen-d/genotype-imputation:latest

# Run directly
docker run -it -v $(pwd)/data:/data ghcr.io/afrigen-d/vcf-processing:latest
```

### **Option 2: Use Docker Compose**
```bash
# Build all containers locally
cd containers/
docker-compose build

# Start specific container
docker-compose up vcf-processing

# Run command in container
docker-compose exec vcf-processing bcftools --version

# Stop all containers
docker-compose down
```

### **Option 3: Use Build Script (Local Development)**
```bash
# Build all containers locally with parallel builds
cd containers/
./build-all.sh

# Build with custom tag
./build-all.sh v1.0.0

# Build sequentially (if parallel fails)
PARALLEL=false ./build-all.sh
```

### **Option 2: Use Individual Containers**
```bash
# Build specific container
docker build -t afrigen-d/vcf-processing ./vcf-processing

# Run with data mounting
docker run -it \
  -v $(pwd)/data:/data \
  -v $(pwd)/input:/input \
  -v $(pwd)/output:/output \
  afrigen-d/vcf-processing
```

### **Option 3: Use with Nextflow**
Update your `nextflow.config`:
```groovy
process {
    withLabel: 'vcf_processing' {
        container = 'afrigen-d/vcf-processing:latest'
    }
    withLabel: 'phasing' {
        container = 'afrigen-d/phasing:latest'  
    }
    withLabel: 'imputation' {
        container = 'afrigen-d/imputation:latest'
    }
    withLabel: 'analysis' {
        container = 'afrigen-d/analysis:latest'
    }
}
```

---

## 🔄 Workflow Integration

### **Process Labels**
Add these labels to your Nextflow processes:

```groovy
process qc_dupl {
    label 'vcf_processing'
    // ... 
}

process minimac4_phasing_eagle {
    label 'phasing'
    // ...
}

process impute_minimac4 {
    label 'imputation' 
    // ...
}

process generate_frequency {
    label 'analysis'
    // ...
}
```

---

## 📊 Benefits of Modular Architecture

1. **🚀 Faster Builds**: Smaller, focused containers build quicker
2. **💾 Efficient Storage**: Avoid downloading unused tools
3. **🔧 Easy Maintenance**: Update tools independently  
4. **📈 Better Caching**: Docker layers cache more effectively
5. **🧪 Easier Testing**: Test individual workflow stages
6. **⚡ Parallel Builds**: Build containers simultaneously

---

## 🤖 GitHub Actions CI/CD

### **Automated Builds**
All containers are automatically built and pushed to GitHub Container Registry on:
- **Push to main/master**: Builds `latest` and `YYYY-MM-DD` tags
- **Pull Requests**: Builds for testing (not pushed)
- **Releases**: Builds with semantic version tags (`v1.0.0`, `v1.0`, `v1`)
- **Weekly Schedule**: Keeps containers updated (Sundays 2 AM UTC)

### **Available Images**
```bash
# Latest versions (updated on every push)
ghcr.io/afrigen-d/vcf-processing:latest
ghcr.io/afrigen-d/phasing:latest
ghcr.io/afrigen-d/imputation:latest
ghcr.io/afrigen-d/analysis:latest
ghcr.io/afrigen-d/genotype-imputation:latest

# Date-tagged versions (for reproducibility)
ghcr.io/afrigen-d/vcf-processing:2025-01-02
ghcr.io/afrigen-d/phasing:2025-01-02
# ... etc
```

### **Monitoring Builds**
- Check build status: `https://github.com/AfriGen-D/genotype-imputation/actions`
- Build matrix creates 5 containers in parallel
- Each container supports `linux/amd64` and `linux/arm64` platforms
- Artifacts include security attestations and SBOM

---

## 🛠️ Development

### **Adding New Tools**
1. Identify which processes use the tool together
2. Add to appropriate existing container, or create new one
3. Update Docker Compose and documentation

### **Testing Containers**
```bash
# Test individual container
docker build -t test-container ./vcf-processing
docker run --rm test-container bcftools --version

# Test full orchestration
docker-compose build --parallel
docker-compose up -d
```

### **Building for Production**
```bash
# Build optimized images
docker-compose build --parallel --no-cache

# Tag for registry
docker tag afrigen-d/vcf-processing:latest ghcr.io/afrigen-d/vcf-processing:latest

# Push to registry
docker push ghcr.io/afrigen-d/vcf-processing:latest
```

---

## 🏷️ Container Size Comparison

| Container | Tools | Estimated Size | Use Case |
|-----------|-------|----------------|----------|
| **vcf-processing** | BCFtools + VCFtools | ~200MB | QC workflows |
| **phasing** | Eagle + tabix | ~150MB | Phasing only |
| **imputation** | Minimac4 + BCFtools | ~250MB | Imputation workflows |
| **analysis** | R + Python + BCFtools | ~800MB | Analysis/plotting |
| **all-in-one** | Everything | ~1.2GB | Complete workflows |

---

## 📝 Legacy Containers

The `legacy/` directory contains the original containers:
- `bcftools/` - Original BCFtools container
- `crossmap/` - CrossMap container  
- `impute5/` - IMPUTE5 container
- `imputation_tools/` - Original comprehensive container

These are preserved for backward compatibility but should be migrated to the new modular architecture.

---

## 🤝 Contributing

1. Analyze tool usage in Nextflow processes
2. Group tools by actual co-usage patterns
3. Create focused, minimal containers
4. Update Docker Compose configuration  
5. Test with real workflows
6. Update documentation

---

For questions or issues, please contact the AfriGen-D development team. 