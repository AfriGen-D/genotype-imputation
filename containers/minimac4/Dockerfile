FROM ubuntu:22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Build Minimac4
WORKDIR /tmp
RUN git clone https://github.com/statgen/Minimac4.git \
    && cd Minimac4 \
    && git checkout tags/v4.1.6 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc)

# Production stage
FROM ubuntu:22.04

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="4.1.6"
LABEL description="Minimac4 genotype imputation software container"

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4-gnutls-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy Minimac4 binary from builder
COPY --from=builder /tmp/Minimac4/build/minimac4 /usr/local/bin/minimac4

# Make executable
RUN chmod +x /usr/local/bin/minimac4

# Create working directory
WORKDIR /data

# Test the installation
RUN minimac4 --help > /dev/null

# Set default command
CMD ["minimac4", "--help"] 