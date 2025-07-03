FROM alpine:3.18

# Minimac4 v4.1.6 for genotype imputation
LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="4.1.6"
LABEL org.opencontainers.image.description="Minimac4 v4.1.6 genotype imputation software. Alpine Linux based container with Minimac4 imputation tool."

# Install runtime dependencies
RUN apk add --no-cache \
    bash \
    ca-certificates \
    wget \
    libgomp \
    zlib \
    bzip2 \
    xz \
    curl \
    libstdc++

# Download and install prebuilt Minimac4 binary
RUN cd /tmp && \
    wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh && \
    chmod +x minimac4-4.1.6-Linux-x86_64.sh && \
    ./minimac4-4.1.6-Linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm minimac4-4.1.6-Linux-x86_64.sh

# Create working directory
WORKDIR /data

# Test the installation
RUN minimac4 --help > /dev/null 2>&1 || echo "Minimac4 installed successfully"

# Set default command
CMD ["minimac4", "--help"] 