#!/bin/bash

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
PARALLEL=${PARALLEL:-true}
TAG=${1:-latest}
REGISTRY=${REGISTRY:-ghcr.io/afrigen-d}

# Custom containers to build (only ones not available externally)
CUSTOM_CONTAINERS=(
    "eagle-phasing"
    "minimac4" 
    "minimac4-all"
    "plink2"
    "all-in-one"
)

# Function to print colored output
print_step() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to build a container
build_container() {
    local container=$1
    local tag=$2
    
    print_step "Building ${container}:${tag}..."
    
    if [[ ! -d "$container" ]]; then
        print_error "Directory $container does not exist"
        return 1
    fi
    
    # Build the container (x86_64/amd64 only for HPC compatibility)
    if docker build \
        --platform linux/amd64 \
        --tag "${REGISTRY}/${container}:${tag}" \
        --tag "${REGISTRY}/${container}:latest" \
        "./${container}"; then
        print_success "Built ${container}:${tag}"
        return 0
    else
        print_error "Failed to build ${container}:${tag}"
        return 1
    fi
}

# Function to check if external containers are available
check_external_containers() {
    print_step "Checking availability of external containers..."
    
    local external_containers=(
        "ghcr.io/wtsi-npg/samtools:latest"
        "quay.io/biocontainers/vcftools:latest"
        "python:3.11-slim"
        "bioconductor/bioconductor_docker:latest"
    )
    
    for container in "${external_containers[@]}"; do
        if docker manifest inspect "$container" >/dev/null 2>&1; then
            print_success "✓ $container is available"
        else
            print_warning "⚠ $container may not be available (check your internet connection)"
        fi
    done
}

# Function to pull external containers
pull_external_containers() {
    print_step "Pulling external containers..."
    
    local external_containers=(
        "ghcr.io/wtsi-npg/samtools:latest"
        "quay.io/biocontainers/vcftools:latest"
        "python:3.11-slim"
        "bioconductor/bioconductor_docker:latest"
    )
    
    for container in "${external_containers[@]}"; do
        print_step "Pulling $container..."
        if docker pull "$container"; then
            print_success "✓ Pulled $container"
        else
            print_warning "⚠ Failed to pull $container (may not be available)"
        fi
    done
}

# Function to build containers in parallel
build_parallel() {
    local tag=$1
    local pids=()
    local failed_containers=()
    
    print_step "Building custom containers in parallel (linux/amd64)..."
    
    # Start builds in background
    for container in "${CUSTOM_CONTAINERS[@]}"; do
        build_container "$container" "$tag" &
        pids+=($!)
    done
    
    # Wait for all builds to complete
    for i in "${!pids[@]}"; do
        local pid=${pids[$i]}
        local container=${CUSTOM_CONTAINERS[$i]}
        
        if wait $pid; then
            print_success "✓ ${container} completed successfully"
        else
            print_error "✗ ${container} failed"
            failed_containers+=("$container")
        fi
    done
    
    # Report results
    if [[ ${#failed_containers[@]} -eq 0 ]]; then
        print_success "All custom containers built successfully!"
        return 0
    else
        print_error "Failed to build: ${failed_containers[*]}"
        return 1
    fi
}

# Function to build containers sequentially
build_sequential() {
    local tag=$1
    local failed_containers=()
    
    print_step "Building custom containers sequentially (linux/amd64)..."
    
    for container in "${CUSTOM_CONTAINERS[@]}"; do
        if build_container "$container" "$tag"; then
            print_success "✓ ${container} completed successfully"
        else
            print_error "✗ ${container} failed"
            failed_containers+=("$container")
        fi
    done
    
    # Report results
    if [[ ${#failed_containers[@]} -eq 0 ]]; then
        print_success "All custom containers built successfully!"
        return 0
    else
        print_error "Failed to build: ${failed_containers[*]}"
        return 1
    fi
}

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 [TAG] [OPTIONS]

Build custom containers for genotype imputation workflow.
Uses existing containers where available, builds only missing tools.
Optimized for HPC environments (linux/amd64 only).

Arguments:
  TAG           Docker tag to use (default: latest)

Environment Variables:
  PARALLEL      Build in parallel (true/false, default: true)
  REGISTRY      Docker registry prefix (default: ghcr.io/afrigen-d)

Examples:
  $0                    # Build with latest tag
  $0 v1.0.0            # Build with specific tag
  PARALLEL=false $0    # Build sequentially
  
Custom containers built (linux/amd64):
  - eagle-phasing      (Eagle 2.4.1 for haplotype phasing)
  - minimac4          (Minimac4 1.0.2 for genotype imputation)
  - minimac4-all      (Comprehensive Minimac4 + genomics tools for Azure)
  - plink2            (PLINK 2.0 for genetic analysis)
  - all-in-one        (Complete environment)

External containers used:
  - ghcr.io/wtsi-npg/samtools:latest (BCFtools, SAMtools, HTSlib)
  - quay.io/biocontainers/vcftools:latest (VCFtools)
  - python:3.11-slim (Python environment)
  - bioconductor/bioconductor_docker:latest (R environment)

EOF
}

# Main function
main() {
    # Parse arguments
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        show_usage
        exit 0
    fi
    
    print_step "Starting genotype imputation custom container builds..."
    print_step "Platform: linux/amd64 (HPC optimized)"
    print_step "Registry: $REGISTRY"
    print_step "Tag: $TAG"
    print_step "Parallel: $PARALLEL"
    
    # Check Docker is available
    if ! docker --version >/dev/null 2>&1; then
        print_error "Docker is not available. Please install Docker first."
        exit 1
    fi
    
    # Check external containers
    check_external_containers
    
    # Optionally pull external containers
    read -p "Do you want to pull external containers now? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        pull_external_containers
    fi
    
    # Build custom containers
    if [[ "$PARALLEL" == "true" ]]; then
        if build_parallel "$TAG"; then
            print_success "Custom container build completed successfully!"
        else
            print_error "Some custom containers failed to build"
            exit 1
        fi
    else
        if build_sequential "$TAG"; then
            print_success "Custom container build completed successfully!"
        else
            print_error "Some custom containers failed to build"
            exit 1
        fi
    fi
    
    # Show final status
    print_step "Container build summary:"
    echo "Custom containers (built locally, linux/amd64):"
    for container in "${CUSTOM_CONTAINERS[@]}"; do
        echo "  ✓ ${REGISTRY}/${container}:${TAG}"
    done
    
    echo
    echo "External containers (use directly):"
    echo "  ✓ ghcr.io/wtsi-npg/samtools:latest"
    echo "  ✓ quay.io/biocontainers/vcftools:latest"
    echo "  ✓ python:3.11-slim"
    echo "  ✓ bioconductor/bioconductor_docker:latest"
    
    print_success "Build completed! Use 'docker-compose up' to start services."
}

# Run main function
main "$@" 