############################################################
# Analysis Container
# Tools: BCFtools + Python + R
# Used in: generate_frequency, plot_freq_comparison, 
#          filter_info_by_target, all reporting processes
# Based on Alpine Linux for better efficiency
############################################################

FROM alpine:3.18

LABEL maintainer="AfriGen-D Consortium"
LABEL version="2.0"
LABEL org.opencontainers.image.description="Analysis and reporting toolkit with BCFtools v1.20, Python 3.11, R 4.3.0, and data analysis packages. Alpine Linux based container for genotype imputation analysis workflows."

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
    gzip \
    python3 \
    python3-dev \
    py3-pip \
    py3-setuptools \
    py3-wheel \
    gcc \
    gfortran \
    musl-dev \
    linux-headers \
    lapack-dev \
    openblas-dev

# Install R and development packages
RUN apk add --no-cache \
    R \
    R-dev

WORKDIR /opt/tools

# HTSlib 1.20 (latest stable)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && \
    ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.20*

# BCFtools 1.20 (latest stable)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar -xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bcftools-1.20*

# Install Python packages for analysis (use latest compatible versions)
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel && \
    pip3 install --no-cache-dir \
        numpy \
        pandas \
        scipy \
        matplotlib \
        seaborn \
        biopython

# Install R packages for plotting and analysis
RUN R -e "install.packages(c('data.table', 'ggplot2', 'dplyr', 'tidyr', 'readr', 'optparse', 'ggsci'), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Create info script and test installations
RUN echo '#!/bin/bash' > /usr/local/bin/analysis_tools_info && \
    echo 'echo "=== Analysis Tools ==="' >> /usr/local/bin/analysis_tools_info && \
    echo 'echo "BCFtools: $(bcftools --version | head -1)"' >> /usr/local/bin/analysis_tools_info && \
    echo 'echo "Python: $(python3 --version)"' >> /usr/local/bin/analysis_tools_info && \
    echo 'echo "R: $(R --version | head -1)"' >> /usr/local/bin/analysis_tools_info && \
    echo 'echo "Testing Python packages..."' >> /usr/local/bin/analysis_tools_info && \
    echo 'python3 -c "import pandas, numpy, matplotlib, seaborn, scipy, Bio; print(\"Python packages: OK\")"' >> /usr/local/bin/analysis_tools_info && \
    echo 'echo "Testing R packages..."' >> /usr/local/bin/analysis_tools_info && \
    echo 'R -e "library(ggplot2); library(data.table); library(dplyr); cat(\"R packages: OK\\n\")" --slave' >> /usr/local/bin/analysis_tools_info && \
    chmod +x /usr/local/bin/analysis_tools_info

# Verify core installations
RUN bcftools --version >/dev/null && \
    python3 --version >/dev/null && \
    R --version >/dev/null

CMD ["analysis_tools_info"] 