############################################################
# Dockerfile for Genotype Imputation Environment
# Updated: July 2025
# Base: Debian 12 (Bookworm) - More stable than Ubuntu
############################################################

FROM debian:bookworm

LABEL maintainer="Genotype Imputation Team"
LABEL version="2.0"
LABEL org.opencontainers.image.description="Comprehensive genotype imputation analysis environment with BCFtools v1.22, SAMtools v1.22, Minimac4 v1.0.2, Eagle v2.4.1, PLINK2, Python 3, and R. Debian-based container with complete bioinformatics toolkit."

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Set working directory
WORKDIR /opt

############################################################
# System dependencies and basic tools
############################################################

RUN apt-get update && apt-get install -y \
    # Build essentials
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    ca-certificates \
    # System libraries for cross-compilation
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    # GSL library for BCFtools (required for polysomy command)
    libgsl-dev \
    # Python and R
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    r-base \
    r-base-dev \
    # Utilities
    wget \
    curl \
    unzip \
    bzip2 \
    git \
    vim \
    less \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

############################################################
# Install modern bioinformatics tools
############################################################

# Set tool versions (latest as of Jan 2025)
ENV HTSLIB_VERSION=1.22
ENV SAMTOOLS_VERSION=1.22
ENV BCFTOOLS_VERSION=1.22
ENV BEDTOOLS_VERSION=2.31.1
ENV VCFTOOLS_VERSION=0.1.16
ENV PLINK2_VERSION=20250627
ENV MINIMAC4_VERSION=1.0.2
ENV EAGLE_VERSION=2.4.1

# Install HTSlib, SAMtools, and BCFtools (same source, compatible versions)
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local --enable-libcurl && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf htslib-${HTSLIB_VERSION}*

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VERSION}*

RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --enable-libgsl && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf bcftools-${BCFTOOLS_VERSION}*

# Install BEDtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary && \
    chmod +x bedtools.static.binary && \
    mv bedtools.static.binary /usr/local/bin/bedtools

# Install VCFtools
RUN wget https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    tar -xzf vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    cd vcftools-${VCFTOOLS_VERSION} && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf vcftools-${VCFTOOLS_VERSION}*

# Install PLINK 2.0 (latest alpha)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_${PLINK2_VERSION}.zip && \
    unzip plink2_linux_x86_64_${PLINK2_VERSION}.zip && \
    mv plink2 /usr/local/bin/ && \
    rm plink2_linux_x86_64_${PLINK2_VERSION}.zip

# Install PLINK 1.9 (still useful for some legacy analyses)
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/plink1.9 && \
    rm plink_linux_x86_64_20231211.zip

# Install Minimac4
RUN wget https://github.com/statgen/Minimac4/releases/download/v${MINIMAC4_VERSION}/Minimac4-${MINIMAC4_VERSION}-Linux.sh && \
    chmod +x Minimac4-${MINIMAC4_VERSION}-Linux.sh && \
    ./Minimac4-${MINIMAC4_VERSION}-Linux.sh --skip-license --prefix=/usr/local && \
    rm Minimac4-${MINIMAC4_VERSION}-Linux.sh

# Install Eagle
RUN wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v${EAGLE_VERSION}.tar.gz && \
    tar -xzf Eagle_v${EAGLE_VERSION}.tar.gz && \
    mv Eagle_v${EAGLE_VERSION}/eagle /usr/local/bin/ && \
    rm -rf Eagle_v${EAGLE_VERSION}*

############################################################
# Install modern variant analysis tools
############################################################

# Install Beagle 5.4 (latest)
RUN wget https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar && \
    mv beagle.28Jun21.220.jar /usr/local/bin/beagle.jar

# Install SHAPEIT4 (modern phasing tool)
RUN wget https://github.com/odelaneau/shapeit4/archive/master.zip && \
    unzip master.zip && \
    cd shapeit4-master && \
    make -j$(nproc) && \
    cp bin/shapeit4 /usr/local/bin/ && \
    cd .. && \
    rm -rf shapeit4-master master.zip

# Install IMPUTE5 (latest imputation software)
RUN wget https://www.dropbox.com/s/mwnceyhir8yze13/impute5_v1.2.0.zip && \
    unzip impute5_v1.2.0.zip && \
    mv impute5_v1.2.0/impute5_1.2.0_static /usr/local/bin/impute5 && \
    chmod +x /usr/local/bin/impute5 && \
    rm -rf impute5_v1.2.0*

############################################################
# Install Python packages for analysis
############################################################

RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    jupyter \
    scikit-learn \
    pysam \
    cyvcf2 \
    pandas-plink \
    plotly \
    dash

############################################################
# Install R packages for genetic analysis
############################################################

RUN R -e "install.packages(c('data.table', 'ggplot2', 'dplyr', 'tidyr', 'optparse', 'R.utils', 'BH'), repos='https://cloud.r-project.org/')"

############################################################
# Create analysis directories and set up environment
############################################################

RUN mkdir -p /data/{input,output,reference,temp} && \
    mkdir -p /scripts && \
    mkdir -p /logs

# Create non-root user for security
RUN useradd -m -s /bin/bash analyst && \
    chown -R analyst:analyst /data /scripts /logs

# Set up environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Add helpful aliases
RUN echo 'alias ll="ls -la"' >> /home/analyst/.bashrc && \
    echo 'alias plink2="/usr/local/bin/plink2"' >> /home/analyst/.bashrc && \
    echo 'alias plink1="/usr/local/bin/plink1.9"' >> /home/analyst/.bashrc

############################################################
# Final setup and validation
############################################################

# Verify installations
RUN samtools --version && \
    bcftools --version && \
    bedtools --version && \
    vcftools --version && \
    plink2 --version && \
    python3 --version && \
    R --version

# Set working directory and user
WORKDIR /data
USER analyst

# Default command
CMD ["/bin/bash"]

############################################################
# Usage examples and documentation
############################################################

# Build: docker build -t genotype-imputation:latest .
# Run: docker run -it -v /path/to/data:/data genotype-imputation:latest
# With Jupyter: docker run -p 8888:8888 genotype-imputation:latest jupyter notebook --ip=0.0.0.0 --allow-root
