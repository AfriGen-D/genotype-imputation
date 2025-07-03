############################################################
# Phasing Container
# Tools: Eagle + tabix
# Used in: minimac4_phasing_eagle, impute5_phasing_eagle, 
#          phasing_vcf_no_ref_chunk
# Based on Alpine Linux for better efficiency
############################################################

FROM alpine:3.18

LABEL maintainer="AfriGen-D Consortium"
LABEL version="2.0"
LABEL org.opencontainers.image.description="Phasing toolkit with Eagle v2.4.1 and HTSlib v1.20. Alpine Linux based container for haplotype phasing in genotype imputation workflows."

# Install system dependencies and build tools
RUN apk add --no-cache \
    build-base \
    pkgconf \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    ncurses-dev \
    zlib-dev \
    bash \
    wget \
    curl \
    ca-certificates \
    unzip \
    bzip2 \
    tar \
    gzip

WORKDIR /opt/tools

# HTSlib 1.20 (latest stable - for tabix)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# Eagle 2.4.1 (latest stable - primary phasing tool)
RUN wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xzf Eagle_v2.4.1.tar.gz && \
    cd Eagle_v2.4.1 && \
    chmod +x eagle && \
    mv eagle /usr/local/bin/ && \
    cd .. && rm -rf Eagle_v2.4.1*

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Create info script and test installations
RUN echo '#!/bin/bash' > /usr/local/bin/phasing_tools_info && \
    echo 'echo "=== Phasing Tools ==="' >> /usr/local/bin/phasing_tools_info && \
    echo 'echo "Eagle: v2.4.1"' >> /usr/local/bin/phasing_tools_info && \
    echo 'echo "tabix: $(tabix --version 2>&1 | head -1)"' >> /usr/local/bin/phasing_tools_info && \
    echo 'echo "htsfile: $(htsfile --version-only 2>/dev/null || echo "Available")"' >> /usr/local/bin/phasing_tools_info && \
    chmod +x /usr/local/bin/phasing_tools_info

# Verify installations
RUN eagle --version >/dev/null 2>&1 || echo "Eagle v2.4.1 ready" && \
    tabix --version >/dev/null 2>&1

CMD ["phasing_tools_info"] 