FROM alpine:3.18

# Minimac4-all: Multi-platform Azure-compatible container
LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="4.1.6-multiplatform"
LABEL description="Multi-platform Minimac4 v4.1.6 with Azure compatibility"
LABEL platforms="linux/amd64,linux/arm64"

# Set timezone for Azure
ENV TZ=UTC

# Install runtime dependencies for multi-platform support
RUN apk add --no-cache \
    bash \
    ca-certificates \
    wget \
    curl \
    # Multi-platform compatible packages
    build-base \
    cmake \
    git \
    autoconf \
    automake \
    libtool \
    pkgconfig \
    # Runtime libraries
    libgomp \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    libstdc++ \
    # Azure tools
    python3 \
    py3-pip

# Install minimal Python packages for Azure
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    pandas

# Build Minimac4 from source for multi-platform compatibility
# This ensures it works on both AMD64 and ARM64
ENV MINIMAC4_VERSION=4.1.6
RUN echo "Building Minimac4 v${MINIMAC4_VERSION} for $(uname -m) architecture..." && \
    git clone --depth 1 --branch v${MINIMAC4_VERSION} https://github.com/statgen/Minimac4.git && \
    cd Minimac4 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j$(nproc) && \
    make install && \
    echo "✅ Minimac4 built successfully for $(uname -m)" && \
    cd ../.. && \
    rm -rf Minimac4

# Create Azure-compatible directories
RUN mkdir -p /data/{input,output,reference,temp} && \
    mkdir -p /scripts && \
    mkdir -p /logs

# Create non-root user for Azure security compliance
RUN adduser -D -s /bin/bash impuser && \
    chown -R impuser:impuser /data /scripts /logs

# Set environment variables
ENV PATH="/usr/local/bin:${PATH}"

# Test the installation
RUN minimac4 --help > /dev/null 2>&1 && echo "✅ Minimac4 installation verified"

# Create multi-platform Azure test script
RUN echo '#!/bin/bash' > /usr/local/bin/test-azure-tools && \
    echo 'echo "=== Minimac4-All Multi-Platform Azure Test ==="' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🏠 Container: $(hostname)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "👤 User: $(whoami)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🏗️ Architecture: $(uname -m)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🐧 OS: $(cat /etc/os-release | grep PRETTY_NAME | cut -d= -f2 | tr -d \")"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "📁 Working directory: $(pwd)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🧬 Minimac4: $(minimac4 --version 2>&1 | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🐍 Python: $(python3 --version)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "✅ Azure-compatible Minimac4 ready on $(uname -m)!"' >> /usr/local/bin/test-azure-tools && \
    chmod +x /usr/local/bin/test-azure-tools

# Set Azure-compatible defaults
WORKDIR /data
USER impuser

# Set default command
CMD ["/bin/bash", "-l"] 