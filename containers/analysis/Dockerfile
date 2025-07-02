############################################################
# Analysis Container
# Tools: BCFtools + Python + R
# Used in: generate_frequency, plot_freq_comparison, 
#          filter_info_by_target, all reporting processes
############################################################

FROM ubuntu:22.04

LABEL maintainer="AfriGen-D Consortium"
LABEL version="1.0"
LABEL description="Analysis and reporting tools for genotype imputation workflow"

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
    python3 \
    python3-pip \
    python3-dev \
    r-base \
    r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tools

# HTSlib for BCFtools
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 \
    && tar -xjf htslib-1.22.tar.bz2 \
    && cd htslib-1.22 \
    && ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf htslib-1.22*

# BCFtools (used for frequency analysis and VCF queries)
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -xjf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local --enable-libgsl \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. && rm -rf bcftools-1.22*

# Python packages for analysis
RUN pip3 install --no-cache-dir \
    numpy==1.24.3 \
    pandas==2.0.3 \
    scipy==1.11.1 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    biopython==1.81

# R packages for plotting and analysis
RUN R -e "install.packages(c('data.table', 'ggplot2', 'dplyr', 'tidyr', 'readr', 'optparse', 'sm', 'ggsci'), repos='https://cran.rstudio.com/')"

# Create data directories
RUN mkdir -p /data /input /output

WORKDIR /data

# Test analysis tools
RUN echo "Testing analysis tools..." \
    && bcftools --version \
    && python3 --version \
    && R --version \
    && python3 -c "import pandas, numpy, matplotlib; print('Python packages OK')" \
    && R -e "library(ggplot2); library(data.table); cat('R packages OK\n')"

CMD ["/bin/bash"] 