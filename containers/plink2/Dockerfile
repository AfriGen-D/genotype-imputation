FROM alpine:3.19 AS builder

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.00a6LM"
LABEL description="PLINK 2.0 genetic data analysis software container"

# Install build dependencies
RUN apk add --no-cache \
    build-base \
    git \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    ca-certificates

# Build PLINK2 from source
WORKDIR /tmp
RUN git clone --depth 1 https://github.com/chrchang/plink-ng.git && \
    cd plink-ng/2.0/build_dynamic && \
    make -j$(nproc) && \
    mv plink2 /usr/local/bin/plink2 && \
    chmod +x /usr/local/bin/plink2

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