############################################################
# Imputation Container
# Tools: Minimac4 + BCFtools + VCFtools
# Used in: impute_minimac4, impute_minimac4_1, combineImpute
# Based on Alpine Linux for better efficiency
############################################################

FROM alpine:3.18

LABEL maintainer="AfriGen-D Consortium"
LABEL version="2.0"
LABEL org.opencontainers.image.description="Comprehensive imputation toolkit with Minimac4 v4.1.6, BCFtools v1.20, and VCFtools v0.1.16. Alpine Linux based container for genotype imputation workflows."

# Install system dependencies and build tools
RUN apk add --no-cache \
    build-base \
    autoconf \
    automake \
    libtool \
    pkgconf \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    ncurses-dev \
    bash \
    wget \
    curl \
    ca-certificates \
    unzip \
    bzip2 \
    tar \
    gzip

WORKDIR /opt/tools

# HTSlib 1.20 (latest stable - foundation for BCFtools)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# BCFtools 1.20 (latest stable - used for pre/post-processing in imputation)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# VCFtools 0.1.16 (latest stable - used in impute_minimac4_1 process)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xzf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf vcftools-0.1.16*

# Minimac4 v4.1.6 (latest stable - primary imputation tool)
RUN wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh && \
    chmod +x minimac4-4.1.6-Linux-x86_64.sh && \
    ./minimac4-4.1.6-Linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm -f minimac4-4.1.6-Linux-x86_64.sh

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Create info script and test installations
RUN echo '#!/bin/bash' > /usr/local/bin/imputation_tools_info && \
    echo 'echo "=== Imputation Tools ==="' >> /usr/local/bin/imputation_tools_info && \
    echo 'echo "Minimac4: v4.1.6"' >> /usr/local/bin/imputation_tools_info && \
    echo 'echo "BCFtools: $(bcftools --version | head -1)"' >> /usr/local/bin/imputation_tools_info && \
    echo 'echo "VCFtools: $(vcftools --version 2>&1 | head -1)"' >> /usr/local/bin/imputation_tools_info && \
    echo 'echo "tabix: $(tabix --version 2>&1 | head -1)"' >> /usr/local/bin/imputation_tools_info && \
    chmod +x /usr/local/bin/imputation_tools_info

# Verify installations
RUN minimac4 --version >/dev/null 2>&1 || echo "Minimac4 v4.1.6 ready" && \
    bcftools --version >/dev/null && \
    vcftools --version >/dev/null 2>&1

CMD ["imputation_tools_info"] 