FROM alpine:3.19 AS builder

# Install build dependencies
RUN apk add --no-cache \
    build-base \
    wget \
    tar \
    gzip \
    ca-certificates

# Download and extract Eagle
WORKDIR /tmp
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xzf Eagle_v2.4.1.tar.gz

# Production stage
FROM alpine:3.19

LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="2.4.1"
LABEL org.opencontainers.image.description="Eagle v2.4.1 haplotype phasing software for genotype imputation workflows. Alpine Linux based container with Eagle phasing tool."

# Install runtime dependencies
RUN apk add --no-cache \
    libgomp \
    libgcc \
    libstdc++

# Copy Eagle binary from builder
COPY --from=builder /tmp/Eagle_v2.4.1/eagle /usr/local/bin/eagle

# Make executable
RUN chmod +x /usr/local/bin/eagle

# Create working directory
WORKDIR /data

# Test the installation (Eagle may need input files, so just check if binary exists and is executable)
RUN test -x /usr/local/bin/eagle && echo "Eagle binary is ready"

# Set default command
CMD ["eagle", "--help"] 