#!/bin/bash

#######################################################################
# Build All Genotype Imputation Containers
# This script builds all modular containers locally for development
#######################################################################

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Container definitions (name:context:description)
CONTAINERS=(
    "vcf-processing:./vcf-processing:VCF processing tools (BCFtools, VCFtools, tabix)"
    "phasing:./phasing:Haplotype phasing tools (Eagle, tabix)"  
    "imputation:./imputation:Genotype imputation tools (Minimac4, BCFtools, VCFtools)"
    "analysis:./analysis:Analysis and reporting tools (BCFtools, Python, R)"
    "genotype-imputation:./all-in-one:Complete genotype imputation environment"
)

# Configuration
REGISTRY="ghcr.io/afrigen-d"
TAG="${1:-latest}"
PARALLEL="${PARALLEL:-true}"
BUILD_ARGS="--no-cache"

echo -e "${BLUE}🏗️  Building Genotype Imputation Containers${NC}"
echo -e "${BLUE}===============================================${NC}"
echo "Registry: $REGISTRY"
echo "Tag: $TAG"
echo "Parallel: $PARALLEL"
echo ""

# Function to build a single container
build_container() {
    local container_def="$1"
    IFS=':' read -r name context description <<< "$container_def"
    
    echo -e "${YELLOW}📦 Building $name...${NC}"
    echo "   Context: $context"
    echo "   Description: $description"
    
    if docker build $BUILD_ARGS -t "$REGISTRY/$name:$TAG" "$context"; then
        echo -e "${GREEN}✅ $name built successfully${NC}"
        
        # Test the container
        echo -e "${BLUE}🧪 Testing $name...${NC}"
        if docker run --rm "$REGISTRY/$name:$TAG" /bin/bash -c "echo 'Container test successful'"; then
            echo -e "${GREEN}✅ $name tested successfully${NC}"
        else
            echo -e "${RED}❌ $name test failed${NC}"
            return 1
        fi
    else
        echo -e "${RED}❌ $name build failed${NC}"
        return 1
    fi
    echo ""
}

# Function to build all containers in parallel
build_parallel() {
    echo -e "${BLUE}🚀 Building containers in parallel...${NC}"
    
    # Start all builds in background
    for container_def in "${CONTAINERS[@]}"; do
        build_container "$container_def" &
    done
    
    # Wait for all builds to complete
    wait
    
    echo -e "${GREEN}🎉 All parallel builds completed!${NC}"
}

# Function to build all containers sequentially
build_sequential() {
    echo -e "${BLUE}🔄 Building containers sequentially...${NC}"
    
    for container_def in "${CONTAINERS[@]}"; do
        build_container "$container_def"
    done
    
    echo -e "${GREEN}🎉 All sequential builds completed!${NC}"
}

# Main execution
main() {
    # Check if we're in the containers directory
    if [[ ! -f "docker-compose.yml" ]]; then
        echo -e "${RED}❌ Error: Must run from containers/ directory${NC}"
        echo "Usage: cd containers && ./build-all.sh [tag]"
        exit 1
    fi
    
    # Check if Docker is available
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}❌ Error: Docker not found${NC}"
        exit 1
    fi
    
    echo "Starting build process..."
    start_time=$(date +%s)
    
    if [[ "$PARALLEL" == "true" ]]; then
        build_parallel
    else
        build_sequential
    fi
    
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    echo ""
    echo -e "${GREEN}🎯 Build Summary${NC}"
    echo -e "${GREEN}===============${NC}"
    echo "Total time: ${duration} seconds"
    echo ""
    echo "Built containers:"
    for container_def in "${CONTAINERS[@]}"; do
        IFS=':' read -r name context description <<< "$container_def"
        echo "  - $REGISTRY/$name:$TAG"
    done
    echo ""
    echo -e "${BLUE}💡 Usage Examples:${NC}"
    echo "  # Run VCF processing container"
    echo "  docker run -it -v \$(pwd)/data:/data $REGISTRY/vcf-processing:$TAG"
    echo ""
    echo "  # Use with Docker Compose"
    echo "  docker-compose up vcf-processing"
    echo ""
    echo "  # Push to registry"
    echo "  docker push $REGISTRY/vcf-processing:$TAG"
}

# Help function
show_help() {
    echo "Build All Genotype Imputation Containers"
    echo ""
    echo "Usage: $0 [tag] [options]"
    echo ""
    echo "Arguments:"
    echo "  tag                Tag for the containers (default: latest)"
    echo ""
    echo "Environment Variables:"
    echo "  PARALLEL=false     Build sequentially instead of parallel"
    echo ""
    echo "Examples:"
    echo "  $0                 Build all containers with 'latest' tag"
    echo "  $0 v1.0.0          Build all containers with 'v1.0.0' tag"
    echo "  PARALLEL=false $0  Build sequentially"
}

# Parse arguments
case "${1:-}" in
    -h|--help)
        show_help
        exit 0
        ;;
    *)
        main "$@"
        ;;
esac 