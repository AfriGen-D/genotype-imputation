############################################################
# Dockerfile to build bcftools, samtools, htslib, vcftools, bedtools
# Based on Alpine Linux for better stability
############################################################

# Set the base image to Alpine Linux
FROM alpine:3.18

# File Author / Maintainer
LABEL maintainer="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"
LABEL description="BCFtools and genomics tools container - Latest versions"

################## BEGIN INSTALLATION ######################

# Install dependencies
RUN apk add --no-cache \
    build-base \
    autoconf \
    automake \
    git \
    wget \
    curl \
    tar \
    bzip2 \
    xz \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    ncurses-dev \
    ca-certificates \
    bash \
    perl \
    python3 \
    python3-dev

# Install htslib (latest version)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xvf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# Install samtools (latest version)
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.20*

# Install bcftools (latest version)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xvf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# Install VCFTools (latest version)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xvf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf vcftools-0.1.16*

# Install bedtools (latest version)
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/bedtools /usr/local/bin/ && \
    cd .. && rm -rf bedtools2 && rm -f bedtools-2.31.1.tar.gz

# Create working directory
WORKDIR /data

# Add /usr/local/bin to PATH
ENV PATH="/usr/local/bin:$PATH"

# Test installations
RUN bcftools --version && \
    samtools --version && \
    bedtools --version && \
    vcftools --version

CMD ["/bin/bash"]