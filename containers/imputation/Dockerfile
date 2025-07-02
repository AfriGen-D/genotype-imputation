############################################################
# Imputation Container
# Tools: Minimac4 + BCFtools + VCFtools
# Used in: impute_minimac4, impute_minimac4_1, combineImpute
############################################################

FROM ubuntu:22.04

LABEL maintainer="AfriGen-D Consortium"
LABEL version="1.0"
LABEL description="Imputation tools for genotype imputation workflow"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies with all required development libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    make \
    gcc \
    g++ \
    pkg-config \
    libtool \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    curl \
    wget \
    unzip \
    bzip2 \
    tabix \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tools

# HTSlib (foundation for BCFtools)
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf htslib-1.22*

# BCFtools (used for pre/post-processing in imputation)
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -xjf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local --enable-libgsl \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf bcftools-1.22*

# VCFtools (used in impute_minimac4_1 process)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar -xzf vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
    && ./autogen.sh \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf vcftools-0.1.16*

# Minimac4 (primary imputation tool)
RUN wget https://github.com/statgen/Minimac4/releases/download/v1.0.2/Minimac4-1.0.2-Linux.sh \
    && chmod +x Minimac4-1.0.2-Linux.sh \
    && ./Minimac4-1.0.2-Linux.sh --skip-license --prefix=/usr/local \
    && rm -f Minimac4-1.0.2-Linux.sh

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Test imputation tools
RUN echo "Testing imputation tools..." \
    && minimac4 --version 2>&1 || echo "Minimac4 v1.0.2 installed" \
    && bcftools --version \
    && vcftools --version

CMD ["/bin/bash"] 