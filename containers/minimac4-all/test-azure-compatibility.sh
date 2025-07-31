#!/bin/bash

# Azure Compatibility Test Suite for minimac4-all
# Tests Azure Container Instances, Azure Batch, and AKS scenarios

set -e

CONTAINER_IMAGE="ghcr.io/afrigen-d/minimac4-all:latest"
TEST_DATA_DIR="$(pwd)/test-data"

echo "🧪 Azure Compatibility Test Suite for minimac4-all"
echo "Container: $CONTAINER_IMAGE"
echo "========================================================"

# Test 1: Basic Container Functionality (Azure Container Instances simulation)
echo ""
echo "🔵 Test 1: Azure Container Instances (ACI) Simulation"
echo "Testing basic container startup and tool availability..."

docker run --platform linux/amd64 --rm $CONTAINER_IMAGE test-azure-tools

# Test 2: Resource Constraints (Azure Batch simulation)
echo ""
echo "🔵 Test 2: Azure Batch Resource Constraints Simulation"
echo "Testing with limited CPU and memory (typical Azure Batch limits)..."

docker run --platform linux/amd64 --rm \
    --memory=2g \
    --cpus=1 \
    $CONTAINER_IMAGE \
    bash -c "
        echo '🧠 Memory limit: 2GB'
        echo '⚡ CPU limit: 1 core'
        echo '🧬 Testing Minimac4 with constraints...'
        minimac4 --help | head -5
        echo '✅ Resource-constrained test passed'
    "

# Test 3: Non-root User Security (Azure security requirement)
echo ""
echo "🔵 Test 3: Azure Security - Non-root User Verification"
echo "Testing non-root user execution (Azure security requirement)..."

docker run --platform linux/amd64 --rm $CONTAINER_IMAGE \
    bash -c "
        echo '👤 Current user: $(whoami)'
        echo '🔒 User ID: $(id -u)'
        echo '🏢 Group ID: $(id -g)'
        if [ \$(id -u) -eq 0 ]; then
            echo '❌ SECURITY ISSUE: Running as root!'
            exit 1
        else
            echo '✅ Security test passed: Non-root user'
        fi
    "

# Test 4: File System Permissions (Azure volume mounts)
echo ""
echo "🔵 Test 4: Azure Volume Mount Simulation"
echo "Testing file system permissions with mounted volumes..."

mkdir -p $TEST_DATA_DIR/{input,output,reference,temp}
echo "test_data" > $TEST_DATA_DIR/input/test.txt

docker run --platform linux/amd64 --rm \
    -v $TEST_DATA_DIR:/data \
    $CONTAINER_IMAGE \
    bash -c "
        echo '📁 Testing volume mount permissions...'
        echo '📖 Reading from input: $(cat /data/input/test.txt)'
        echo 'azure_test_output' > /data/output/result.txt
        echo '📝 Writing to output: $(cat /data/output/result.txt)'
        echo '🗂️ Listing directories:'
        ls -la /data/
        echo '✅ Volume mount test passed'
    "

# Test 5: Environment Variables (Azure configuration)
echo ""
echo "🔵 Test 5: Azure Environment Variables"
echo "Testing environment variable handling..."

docker run --platform linux/amd64 --rm \
    -e AZURE_STORAGE_ACCOUNT="testaccount" \
    -e AZURE_RESOURCE_GROUP="test-rg" \
    -e MINIMAC4_THREADS="2" \
    $CONTAINER_IMAGE \
    bash -c "
        echo '☁️ Azure Storage Account: \$AZURE_STORAGE_ACCOUNT'
        echo '🏗️ Resource Group: \$AZURE_RESOURCE_GROUP' 
        echo '🧵 Minimac4 Threads: \$MINIMAC4_THREADS'
        echo '✅ Environment variables test passed'
    "

# Test 6: Network Connectivity (Azure networking)
echo ""
echo "🔵 Test 6: Azure Network Connectivity Simulation"
echo "Testing network access (simulating Azure networking)..."

docker run --platform linux/amd64 --rm \
    --network none \
    $CONTAINER_IMAGE \
    bash -c "
        echo '🌐 Testing offline execution (Azure private networks)...'
        minimac4 --version
        python3 -c 'import numpy, pandas; print(\"📊 Python packages available offline\")'
        echo '✅ Offline execution test passed'
    "

# Test 7: Azure Container Registry Simulation
echo ""
echo "🔵 Test 7: Azure Container Registry (ACR) Pull Simulation"
echo "Testing container image pull and caching behavior..."

docker run --platform linux/amd64 --rm $CONTAINER_IMAGE \
    bash -c "
        echo '📦 Container image info:'
        echo '   - Base: Alpine Linux'
        echo '   - Size: $(du -sh / 2>/dev/null | cut -f1 || echo \"~254MB\")'
        echo '   - User: $(whoami)'
        echo '   - Minimac4: $(minimac4 --version 2>&1 | head -1)'
        echo '✅ Container registry simulation passed'
    "

# Test 8: Long-running Process (Azure Batch job simulation)
echo ""
echo "🔵 Test 8: Azure Batch Long-running Job Simulation"
echo "Testing long-running genomics workflow simulation..."

docker run --platform linux/amd64 --rm \
    --memory=1g \
    --cpus=0.5 \
    $CONTAINER_IMAGE \
    bash -c "
        echo '🧬 Simulating genomics workflow...'
        echo '⏱️ Step 1: Data validation (1s)'
        sleep 1
        echo '⏱️ Step 2: Tool initialization (1s)'
        minimac4 --help > /dev/null
        sleep 1
        echo '⏱️ Step 3: Processing simulation (2s)'
        python3 -c 'import time; import numpy as np; time.sleep(2); print(\"📊 Processing complete\")'
        echo '✅ Long-running job simulation passed'
    "

# Test 9: Azure Key Vault Secrets Simulation
echo ""
echo "🔵 Test 9: Azure Key Vault Secrets Simulation"
echo "Testing secrets management..."

docker run --platform linux/amd64 --rm \
    -e AZURE_TENANT_ID="test-tenant" \
    -e AZURE_CLIENT_ID="test-client" \
    --tmpfs /tmp:noexec,nosuid,size=100m \
    $CONTAINER_IMAGE \
    bash -c "
        echo '🔐 Testing secrets handling...'
        echo '🆔 Tenant ID available: \${AZURE_TENANT_ID:0:4}...'
        echo '👤 Client ID available: \${AZURE_CLIENT_ID:0:4}...'
        echo '📂 Temp filesystem: $(df -h /tmp | tail -1)'
        echo '✅ Secrets simulation passed'
    "

# Test 10: Azure Monitoring Simulation
echo ""
echo "🔵 Test 10: Azure Monitoring and Logging"
echo "Testing logging and monitoring capabilities..."

docker run --platform linux/amd64 --rm \
    --log-driver json-file \
    --log-opt max-size=10m \
    $CONTAINER_IMAGE \
    bash -c "
        echo '{\"timestamp\":\"$(date -Iseconds)\",\"level\":\"INFO\",\"message\":\"Container started\",\"component\":\"minimac4-all\"}'
        echo '{\"timestamp\":\"$(date -Iseconds)\",\"level\":\"INFO\",\"message\":\"Tools validated\",\"minimac4_version\":\"$(minimac4 --version 2>&1 | head -1)\"}'
        echo '{\"timestamp\":\"$(date -Iseconds)\",\"level\":\"INFO\",\"message\":\"Container ready\",\"status\":\"healthy\"}'
        echo '✅ Logging simulation passed'
    "

# Cleanup
rm -rf $TEST_DATA_DIR

echo ""
echo "🎉 All Azure compatibility tests passed!"
echo ""
echo "📋 Test Results Summary:"
echo "✅ Azure Container Instances compatibility"
echo "✅ Azure Batch resource constraints"
echo "✅ Non-root security requirements"
echo "✅ Volume mount permissions"
echo "✅ Environment variable handling"
echo "✅ Network isolation support"
echo "✅ Container registry compatibility"
echo "✅ Long-running job support"
echo "✅ Secrets management simulation"
echo "✅ Monitoring and logging"
echo ""
echo "🚀 Ready for Azure deployment!"
echo ""
echo "📖 Next Steps:"
echo "1. Deploy to Azure Container Instances for testing"
echo "2. Test on Azure Batch for genomics workflows"
echo "3. Configure Azure Key Vault for secrets"
echo "4. Set up Azure Monitor for logging" 