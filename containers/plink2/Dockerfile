FROM ubuntu:22.04

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.00a6LM"
LABEL description="PLINK 2.0 genetic data analysis software container"

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Download and install PLINK 2.0 (x86_64 only)
WORKDIR /tmp
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20240804.zip \
    && unzip plink2_linux_x86_64_20240804.zip \
    && mv plink2 /usr/local/bin/plink2 \
    && chmod +x /usr/local/bin/plink2 \
    && rm -rf /tmp/*

# Create working directory
WORKDIR /data

# Test the installation
RUN plink2 --help > /dev/null

# Set default command
CMD ["plink2", "--help"] 