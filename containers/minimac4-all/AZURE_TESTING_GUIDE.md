# Azure Testing Guide for minimac4-all Container

This guide provides comprehensive testing approaches to ensure the `minimac4-all` container works perfectly on Azure.

## 🧪 Testing Approaches

### 1. Local Azure Simulation Testing

**Script**: `test-azure-compatibility.sh`

Tests Azure-specific conditions locally without requiring an Azure account:

```bash
# Run comprehensive local Azure simulation
./test-azure-compatibility.sh
```

**What it tests:**
- ✅ Azure Container Instances simulation
- ✅ Resource constraints (CPU/memory limits)
- ✅ Non-root security requirements
- ✅ Volume mount permissions
- ✅ Environment variables
- ✅ Network isolation
- ✅ Container registry compatibility
- ✅ Long-running job simulation
- ✅ Secrets management
- ✅ Monitoring and logging

### 2. Real Azure Deployment Testing

**Script**: `deploy-azure-test.sh`

Tests on actual Azure services (requires Azure CLI and account):

```bash
# Prerequisites
az login
az account set --subscription "your-subscription-id"

# Run real Azure tests
./deploy-azure-test.sh
```

**What it tests:**
- 🌐 Azure Container Instances deployment
- 📁 Azure Files integration  
- 📦 Azure Container Registry
- 🔐 Authentication and permissions
- 📊 Performance metrics
- 🧹 Resource cleanup

## 🎯 Azure Service-Specific Testing

### Azure Container Instances (ACI)

**Template**: `azure-templates/aci-genomics.yaml`

Deploy with:
```bash
az container create --resource-group myRG --file aci-genomics.yaml
```

**Features:**
- 4 CPU cores, 8GB memory
- Azure Files volume mount
- Environment variables
- Genomics-optimized configuration

### Azure Batch

**Template**: `azure-templates/batch-job.json`

Deploy with:
```bash
az batch job create --json-file batch-job.json
```

**Features:**
- Parallel chromosome processing
- Job manager for coordination
- Resource file management
- Output file upload
- Task retry logic

## 🔍 Key Azure Compatibility Checks

### Security Requirements
- ✅ **Non-root user**: Container runs as `impuser` (UID 1000)
- ✅ **Read-only filesystem**: Application data in `/data`
- ✅ **No privileged access**: Standard container permissions

### Resource Management
- ✅ **CPU limits**: Tested with 0.5-4 cores
- ✅ **Memory limits**: Tested with 1-8GB RAM
- ✅ **Disk space**: Efficient Alpine base (~254MB)

### Storage Integration
- ✅ **Azure Files**: SMB volume mounts work
- ✅ **Azure Blob**: Compatible with blob storage
- ✅ **Persistent volumes**: Data survives container restarts

### Networking
- ✅ **Private networks**: Works in isolated subnets
- ✅ **Offline execution**: No internet dependencies
- ✅ **Port exposure**: No incoming ports needed

## 🚀 Quick Start Commands

### Test Container Locally (with Azure simulation)
```bash
# Pull the latest container
docker pull --platform linux/amd64 ghcr.io/afrigen-d/minimac4-all:latest

# Run Azure compatibility test
docker run --platform linux/amd64 --rm ghcr.io/afrigen-d/minimac4-all:latest test-azure-tools

# Test with resource constraints (simulating Azure)
docker run --platform linux/amd64 --rm --memory=2g --cpus=1 \
  ghcr.io/afrigen-d/minimac4-all:latest \
  bash -c "echo 'Testing with Azure-like constraints' && minimac4 --help"
```

### Deploy to Azure Container Instances
```bash
# Quick deployment
az container create \
  --resource-group myRG \
  --name minimac4-job \
  --image ghcr.io/afrigen-d/minimac4-all:latest \
  --cpu 2 --memory 4 \
  --command-line "test-azure-tools"

# Check logs
az container logs --resource-group myRG --name minimac4-job --follow
```

### Deploy to Azure Batch
```bash
# Create Batch account and pool first, then:
az batch job create --json-file azure-templates/batch-job.json
az batch task list --job-id minimac4-imputation-job
```

## 📊 Performance Expectations

### Startup Time
- **ACI**: ~30-60 seconds (depending on region)
- **Batch**: ~1-2 minutes (including pool allocation)
- **AKS**: ~10-30 seconds (with pre-pulled images)

### Resource Usage
- **Base memory**: ~50MB (Alpine + tools)
- **Minimac4 memory**: Varies by dataset (typically 1-8GB)
- **CPU utilization**: Scales with thread count
- **Disk I/O**: Optimized for genomics file formats

### Cost Optimization
- **Small jobs**: Use ACI (pay per second)
- **Large jobs**: Use Batch (dedicated pools)
- **High availability**: Use AKS (managed clusters)

## 🔧 Troubleshooting

### Common Issues

**1. Platform Architecture Error**
```bash
# Issue: linux/arm64 vs linux/amd64
# Solution: Always specify platform
docker run --platform linux/amd64 ghcr.io/afrigen-d/minimac4-all:latest
```

**2. Permission Denied on Azure Files**
```bash
# Issue: Volume mount permissions
# Solution: Ensure storage account key is correct
az storage account keys list --account-name mystorageaccount
```

**3. Container Startup Timeout**
```bash
# Issue: Long container startup
# Solution: Increase timeout in Azure template
# Add: "timeoutInSeconds": 300
```

### Validation Commands
```bash
# Check container health
docker run --platform linux/amd64 --rm ghcr.io/afrigen-d/minimac4-all:latest \
  bash -c "minimac4 --help && python3 --version && whoami"

# Check Azure CLI setup
az account show
az container list --resource-group myRG
```

## 📈 Monitoring and Logging

### Azure Monitor Integration
```bash
# Enable container insights
az container create \
  --log-analytics-workspace myWorkspace \
  --log-analytics-workspace-key myKey \
  ...
```

### Log Analysis
```bash
# View container logs
az container logs --resource-group myRG --name minimac4-job

# Stream live logs  
az container logs --resource-group myRG --name minimac4-job --follow

# Export logs for analysis
az monitor log-analytics query \
  --workspace myWorkspace \
  --analytics-query "ContainerInstanceLog_CL | where ContainerGroup_s == 'minimac4-job'"
```

## 🎉 Success Criteria

Your container is Azure-ready when all these tests pass:

- ✅ **Local simulation**: `test-azure-compatibility.sh` completes successfully
- ✅ **ACI deployment**: Container runs and completes on Azure Container Instances  
- ✅ **Storage access**: Can read/write Azure Files and Blob storage
- ✅ **Resource limits**: Handles CPU/memory constraints gracefully
- ✅ **Security compliance**: Runs as non-root user
- ✅ **Monitoring**: Logs appear in Azure Monitor
- ✅ **Cost efficiency**: Completes within expected time/cost budgets

## 📚 Additional Resources

- [Azure Container Instances Documentation](https://docs.microsoft.com/en-us/azure/container-instances/)
- [Azure Batch for HPC Workloads](https://docs.microsoft.com/en-us/azure/batch/)
- [Azure Files Integration](https://docs.microsoft.com/en-us/azure/storage/files/)
- [Container Security Best Practices](https://docs.microsoft.com/en-us/azure/security/fundamentals/container-security)

---

🧬 **Ready to run genomics workloads on Azure!** 🚀 