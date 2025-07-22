############################################################
# Eagle + VCF Processing Container
# Tools: Eagle 2.4.1 + BCFtools + tabix + VCFtools
# Used in: minimac4_phasing_eagle, impute5_phasing_eagle
# Based on Ubuntu for better compatibility with Eagle
############################################################

FROM ubuntu:22.04 AS eagle-builder

# Install build dependencies for Eagle
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    tar \
    gzip \
    ca-certificates \
    libgomp1 \
    libgcc-s1 \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Download and extract Eagle
WORKDIR /tmp
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xzf Eagle_v2.4.1.tar.gz

# Production stage
FROM ubuntu:22.04

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.4.1"
LABEL org.opencontainers.image.description="Eagle v2.4.1 haplotype phasing + VCF processing toolkit. Alpine Linux based container with Eagle, BCFtools v1.20, VCFtools v0.1.16, and HTSlib v1.20."

# Install system dependencies and build tools
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    bash \
    wget \
    curl \
    ca-certificates \
    unzip \
    bzip2 \
    tar \
    gzip \
    libgomp1 \
    libgcc-s1 \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tools

# HTSlib 1.20 (foundation for BCFtools and tabix)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# BCFtools 1.20
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# VCFtools 0.1.16
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xzf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf vcftools-0.1.16*

# Minimac4 4.1.6 (latest stable - for imputation)
RUN wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux.tar.gz && \
    tar -xzf minimac4-4.1.6-Linux.tar.gz && \
    mv minimac4-4.1.6-Linux/minimac4 /usr/local/bin/minimac4 && \
    chmod +x /usr/local/bin/minimac4 && \
    rm -rf minimac4-4.1.6-Linux*

# Copy Eagle binary from builder
COPY --from=eagle-builder /tmp/Eagle_v2.4.1/eagle /usr/local/bin/eagle

# Make Eagle executable and verify
RUN chmod +x /usr/local/bin/eagle && \
    ls -la /usr/local/bin/eagle && \
    /usr/local/bin/eagle --help >/dev/null 2>&1 || echo "Eagle binary verification failed"

# Create working directory
WORKDIR /data

# Test the installation
RUN test -x /usr/local/bin/eagle && echo "Eagle binary is ready" && \
    bcftools --version >/dev/null && echo "BCFtools is ready" && \
    vcftools --version >/dev/null 2>&1 && echo "VCFtools is ready" && \
    tabix --version >/dev/null 2>&1 && echo "tabix is ready" && \
    minimac4 --version >/dev/null 2>&1 && echo "Minimac4 is ready"

# Set default command
CMD ["eagle", "--help"] 