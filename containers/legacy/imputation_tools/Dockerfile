############################################################
# Dockerfile to build comprehensive genotype imputation tools
# Based on Alpine Linux for better stability and efficiency
############################################################

# Set the base image to Alpine Linux
FROM alpine:3.18

# File Author / Maintainer
LABEL maintainer="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"
LABEL description="Comprehensive genotype imputation tools suite - Latest versions"
LABEL version="2024.1"

################## BEGIN INSTALLATION ######################

# Install build dependencies and system tools
RUN apk add --no-cache \
    build-base \
    autoconf \
    automake \
    libtool \
    git \
    wget \
    curl \
    bash \
    ca-certificates \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    ncurses-dev \
    libcurl \
    python3 \
    python3-dev \
    py3-pip \
    unzip \
    tar \
    gzip

# Install latest htslib (1.20)
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xvf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# Install latest samtools (1.20)
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xvf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.20*

# Install latest bcftools (1.20)  
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xvf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# Install latest VCFtools (still 0.1.16)
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xvf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf vcftools-0.1.16*

# Install latest bedtools (2.31.1)
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools2*

# Install IMPUTE2 (prebuilt binary)
RUN wget https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
    mv impute_v2.3.2_x86_64_static/impute2 /usr/local/bin/impute2 && \
    chmod +x /usr/local/bin/impute2 && \
    mkdir -p /opt/impute2/example && \
    mv impute_v2.3.2_x86_64_static/Example/* /opt/impute2/example/ && \
    rm -rf impute_v2.3.2_x86_64_static*

# Install PLINK2 (latest)
RUN wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O plink.zip && \
    unzip plink.zip -d /usr/local/bin/ && \
    rm -f plink.zip

# Install Eagle (latest v2.4.1)
RUN wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xzf Eagle_v2.4.1.tar.gz && \
    mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
    chmod +x /usr/local/bin/eagle && \
    rm -rf Eagle_v2.4.1*

# Install Minimac4 (prebuilt binary)
RUN wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh && \
    chmod +x minimac4-4.1.6-Linux-x86_64.sh && \
    ./minimac4-4.1.6-Linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm minimac4-4.1.6-Linux-x86_64.sh

# Note: Minimac3 is deprecated in favor of Minimac4
# Note: SnpEff/SnpSift require Java which would significantly increase image size
# Users can install these separately if needed

# Create working directory
WORKDIR /data

# Create wrapper scripts for common tasks
RUN echo '#!/bin/bash' > /usr/local/bin/imputation_help && \
    echo 'echo "=== Imputation Tools Suite ==="' >> /usr/local/bin/imputation_help && \
    echo 'echo "Available tools:"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - htslib: $(htsfile --version-only 2>/dev/null || echo unknown)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - samtools: $(samtools --version | head -1)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - bcftools: $(bcftools --version | head -1)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - vcftools: $(vcftools --version 2>&1 | head -1)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - bedtools: $(bedtools --version)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - impute2: Available"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - plink2: $(plink2 --version | head -1)"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - eagle: Available"' >> /usr/local/bin/imputation_help && \
    echo 'echo "  - minimac4: Available"' >> /usr/local/bin/imputation_help && \
    chmod +x /usr/local/bin/imputation_help

# Set default command
CMD ["imputation_help"]