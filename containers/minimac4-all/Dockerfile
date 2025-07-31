FROM alpine:3.19

# Minimac4-all Azure Simple: Core tools for Azure
LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="4.1.6-azure-simple"
LABEL description="Azure-compatible genotype imputation container with core tools (simplified build)"

# Set timezone and working directory
ENV TZ=UTC
WORKDIR /opt

# Install system dependencies and build tools
RUN apk add --no-cache \
    # Core build tools
    build-base \
    cmake \
    autoconf \
    automake \
    libtool \
    pkgconfig \
    git \
    # System libraries
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    ncurses-dev \
    # Download tools
    wget \
    curl \
    unzip \
    bzip2 \
    tar \
    # Runtime dependencies
    ca-certificates \
    bash \
    # Python for scripts
    python3 \
    python3-dev \
    py3-pip

# Install Python packages for data processing
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    pandas \
    scipy \
    matplotlib \
    seaborn

# Build and install HTSlib (required for many genomics tools)
ENV HTSLIB_VERSION=1.19
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf htslib-${HTSLIB_VERSION}*

# Build and install SAMtools
ENV SAMTOOLS_VERSION=1.19
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VERSION}*

# Build and install BCFtools
ENV BCFTOOLS_VERSION=1.19
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*

# Install VCFtools
ENV VCFTOOLS_VERSION=0.1.16
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    tar -xzf vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    cd vcftools-${VCFTOOLS_VERSION} && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf vcftools-${VCFTOOLS_VERSION}*

# Install BEDtools (for genomic interval operations)
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xzf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make -j$(nproc) && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf bedtools-2.31.1.tar.gz bedtools2

# Install Eagle (phasing tool) - binary download
ENV EAGLE_VERSION=2.4.1
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v${EAGLE_VERSION}.tar.gz && \
    tar -xzf Eagle_v${EAGLE_VERSION}.tar.gz && \
    mv Eagle_v${EAGLE_VERSION}/eagle /usr/local/bin/ && \
    chmod +x /usr/local/bin/eagle && \
    rm -rf Eagle_v${EAGLE_VERSION}*

# Create Minimac4 environment setup (for manual installation)
RUN echo '#!/bin/bash' > /usr/local/bin/setup-minimac4 && \
    echo 'echo "🔧 Setting up Minimac4 environment for Azure..."' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "You can install Minimac4 by:"' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "1. Mount your own minimac4 binary to /usr/local/bin/minimac4"' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "2. Or download pre-built binary:"' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "   wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh"' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "   chmod +x minimac4-4.1.6-Linux-x86_64.sh"' >> /usr/local/bin/setup-minimac4 && \
    echo 'echo "   ./minimac4-4.1.6-Linux-x86_64.sh --skip-license --prefix=/usr/local"' >> /usr/local/bin/setup-minimac4 && \
    chmod +x /usr/local/bin/setup-minimac4

# Create working directories for Azure
RUN mkdir -p /data/{input,output,reference,temp} && \
    mkdir -p /scripts && \
    mkdir -p /logs

# Create non-root user for Azure security compliance
RUN adduser -D -s /bin/bash impuser && \
    chown -R impuser:impuser /data /scripts /logs

# Set environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Create Azure-compatible testing script
RUN echo '#!/bin/bash' > /usr/local/bin/test-azure-tools && \
    echo 'echo "=== Azure-Compatible Genomics Tools Test ==="' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🏠 Container: $(hostname)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "👤 User: $(whoami)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "📁 Working directory: $(pwd)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo ""' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🧬 Available Genomics Tools:"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  SAMtools: $(samtools --version | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  BCFtools: $(bcftools --version | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  VCFtools: $(vcftools --version 2>&1 | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  HTSlib (tabix): $(tabix 2>&1 | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  BEDtools: $(bedtools --version)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  Eagle: $(eagle 2>&1 | head -3 | tail -1 || echo \"eagle available\")"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo ""' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🐍 Python Environment:"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  Python: $(python3 --version)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  NumPy: $(python3 -c \"import numpy; print(f\\\"numpy {numpy.__version__}\\\")\"")"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "  Pandas: $(python3 -c \"import pandas; print(f\\\"pandas {pandas.__version__}\\\")\"")"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo ""' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "💡 For Minimac4, run: setup-minimac4"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo ""' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "✅ Azure-compatible genomics environment ready!"' >> /usr/local/bin/test-azure-tools && \
    chmod +x /usr/local/bin/test-azure-tools

# Create welcome message for Azure
RUN echo 'echo "🧬 Azure Genomics Environment"' > /etc/profile.d/welcome.sh && \
    echo 'echo "Core tools ready for Azure deployments"' >> /etc/profile.d/welcome.sh && \
    echo 'echo "Run test-azure-tools to verify installations"' >> /etc/profile.d/welcome.sh && \
    echo 'echo "Run setup-minimac4 for Minimac4 setup instructions"' >> /etc/profile.d/welcome.sh && \
    echo 'echo "Working directory: /data"' >> /etc/profile.d/welcome.sh

# Set default working directory and user for Azure
WORKDIR /data
USER impuser

# Set default command
CMD ["/bin/bash", "-l"] 