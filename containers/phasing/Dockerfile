############################################################
# Phasing Container
# Tools: Eagle + tabix
# Used in: minimac4_phasing_eagle, impute5_phasing_eagle, 
#          phasing_vcf_no_ref_chunk
############################################################

FROM ubuntu:22.04

LABEL maintainer="AfriGen-D Consortium"
LABEL version="1.0"
LABEL description="Phasing tools for genotype imputation workflow"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies with additional build tools
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    zlib1g-dev \
    curl \
    wget \
    unzip \
    bzip2 \
    tabix \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tools

# HTSlib for tabix (with proper configuration for cross-compilation)
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf htslib-1.22*

# Eagle 2.4.1 (primary phasing tool)
RUN wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -xzf Eagle_v2.4.1.tar.gz \
    && cd Eagle_v2.4.1 \
    && chmod +x eagle \
    && mv eagle /usr/local/bin/ \
    && cd .. && rm -rf Eagle_v2.4.1*

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Test phasing tools
RUN echo "Testing phasing tools..." \
    && eagle --version 2>&1 || echo "Eagle v2.4.1 installed" \
    && tabix --version

CMD ["/bin/bash"] 