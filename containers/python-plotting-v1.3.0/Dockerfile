# Python Plotting Container v1.3.0
# For genotype imputation pipeline QC plots and reports

FROM python:3.9-slim

LABEL maintainer="H3ABioNet Genotype Imputation Team"
LABEL version="1.3.0"
LABEL description="Python plotting container for genotype imputation QC"

# Install system dependencies
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
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for plotting and data processing
RUN pip install --no-cache-dir \
    matplotlib==3.7.1 \
    seaborn==0.12.2 \
    pandas==2.0.3 \
    numpy==1.24.3 \
    scipy==1.11.1 \
    scikit-learn==1.3.0 \
    plotly==5.15.0 \
    pysam==0.21.0 \
    cyvcf2==0.30.22

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]