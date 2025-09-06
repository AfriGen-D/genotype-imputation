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
    plotly==5.15.0 \
    Pillow==10.0.0

# Install utility packages
RUN pip install --no-cache-dir \
    psutil==5.9.5 \
    PyYAML==6.0.1

# Note: Removed pysam and cyvcf2 as they're not used in the pipeline
# If genomics packages are needed later, consider a separate container

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]