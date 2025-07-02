FROM ubuntu:22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Download and install Eagle
WORKDIR /tmp
RUN wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -xzf Eagle_v2.4.1.tar.gz

# Production stage
FROM ubuntu:22.04

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.4.1"
LABEL description="Eagle haplotype phasing software container"

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy Eagle binary from builder
COPY --from=builder /tmp/Eagle_v2.4.1/eagle /usr/local/bin/eagle

# Make executable
RUN chmod +x /usr/local/bin/eagle

# Create working directory
WORKDIR /data

# Test the installation
RUN eagle --help > /dev/null

# Set default command
CMD ["eagle", "--help"] 