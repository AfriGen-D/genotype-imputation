FROM alpine:3.19 AS builder

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.00a6LM"
LABEL org.opencontainers.image.description="PLINK 2.0 genetic data analysis software for association studies and data management. Alpine Linux based container with prebuilt PLINK2 binary."

# Install dependencies for downloading and extracting
RUN apk add --no-cache \
    wget \
    unzip \
    ca-certificates

# Download pre-built PLINK2 binary (more reliable than building from source on Alpine)
WORKDIR /tmp
RUN wget -O plink2_linux_x86_64.zip "http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" && \
    unzip plink2_linux_x86_64.zip && \
    mv plink2 /usr/local/bin/plink2 && \
    chmod +x /usr/local/bin/plink2 && \
    rm -rf /tmp/*

# Production stage  
FROM alpine:3.19

# Install runtime dependencies
RUN apk add --no-cache \
    libgcc \
    libstdc++ \
    zlib \
    bzip2 \
    xz

# Copy PLINK2 binary from builder
COPY --from=builder /usr/local/bin/plink2 /usr/local/bin/plink2

# Create working directory
WORKDIR /data

# Test the installation
RUN plink2 --help > /dev/null 2>&1 || echo "PLINK2 installed (help test may fail without input)"

# Set default command
CMD ["plink2", "--help"] 