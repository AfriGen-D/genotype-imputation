# Python Plotting Container v1.3.0
# For genotype imputation pipeline QC plots and reports

FROM python:3.9-slim

LABEL maintainer="H3ABioNet Genotype Imputation Team"
LABEL version="1.3.0"
LABEL description="Python plotting container for genotype imputation QC"

# Install system dependencies needed for compiling packages
RUN apt-get update && apt-get install -y \
    procps \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhts-dev \
    autoconf \
    automake \
    && rm -rf /var/lib/apt/lists/*

# Install numpy first (required by many other packages)
RUN pip install --no-cache-dir numpy==1.24.3

# Install core scientific packages
RUN pip install --no-cache-dir \
    pandas==2.0.3 \
    scipy==1.11.1 \
    scikit-learn==1.3.0

# Install plotting libraries
RUN pip install --no-cache-dir \
    matplotlib==3.7.1 \
    seaborn==0.12.2 \
    plotly==5.15.0

# Install genomics packages (these need compilation)
RUN pip install --no-cache-dir \
    pysam==0.21.0

# Install cyvcf2 separately as it's most likely to fail
RUN pip install --no-cache-dir \
    cyvcf2==0.30.22

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]