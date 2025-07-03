############################################################
# VCF Processing Container
# Tools: BCFtools + tabix + VCFtools
# Used in: qc_dupl, split_multi_allelic, fill_tags_vcf, 
#          filter_min_ac, sites_only, combine_vcfs
# Based on Alpine Linux for better efficiency
############################################################

FROM alpine:3.18

LABEL maintainer="AfriGen-D Consortium"
LABEL version="2.0"
LABEL org.opencontainers.image.description="VCF processing toolkit with BCFtools v1.20, VCFtools v0.1.16, and HTSlib v1.20. Alpine Linux based container for VCF file processing in genotype imputation workflows."

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

# HTSlib 1.20 (latest stable - foundation for BCFtools and tabix)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# BCFtools 1.20 (latest stable - primary tool)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# VCFtools 0.1.16 (latest stable - used in qc_site_missingness)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xzf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf vcftools-0.1.16*

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Test tools work and create info script
RUN echo '#!/bin/bash' > /usr/local/bin/vcf_tools_info && \
    echo 'echo "=== VCF Processing Tools ==="' >> /usr/local/bin/vcf_tools_info && \
    echo 'echo "BCFtools: $(bcftools --version | head -1)"' >> /usr/local/bin/vcf_tools_info && \
    echo 'echo "VCFtools: $(vcftools --version 2>&1 | head -1)"' >> /usr/local/bin/vcf_tools_info && \
    echo 'echo "tabix: $(tabix --version 2>&1 | head -1)"' >> /usr/local/bin/vcf_tools_info && \
    echo 'echo "htsfile: $(htsfile --version-only 2>/dev/null || echo "Available")"' >> /usr/local/bin/vcf_tools_info && \
    chmod +x /usr/local/bin/vcf_tools_info

# Verify installations
RUN bcftools --version >/dev/null && \
    vcftools --version >/dev/null 2>&1 && \
    tabix --version >/dev/null 2>&1

CMD ["vcf_tools_info"] 