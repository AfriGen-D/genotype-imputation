############################################################
# VCF Processing Container
# Tools: BCFtools + tabix + VCFtools
# Used in: qc_dupl, split_multi_allelic, fill_tags_vcf, 
#          filter_min_ac, sites_only, combine_vcfs
############################################################

FROM ubuntu:22.04

LABEL maintainer="AfriGen-D Consortium"
LABEL version="1.0"
LABEL description="VCF processing tools for genotype imputation workflow"

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

# HTSlib 1.22 (foundation for BCFtools and tabix)
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf htslib-1.22*

# BCFtools 1.22 (primary tool)
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -xjf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local --enable-libgsl \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf bcftools-1.22*

# VCFtools 0.1.16 (used in qc_site_missingness)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
    && tar -xzf vcftools-0.1.16.tar.gz \
    && cd vcftools-0.1.16 \
    && ./autogen.sh \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf vcftools-0.1.16*

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Test tools work together
RUN echo "Testing VCF processing tools..." \
    && bcftools --version \
    && vcftools --version \
    && tabix --version

CMD ["/bin/bash"] 