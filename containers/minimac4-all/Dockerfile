FROM alpine:3.18

# Minimac4-all: Azure-compatible container (same as minimac4 + Azure features)
LABEL maintainer="AfriGen-D <info@afrigen.org>"
LABEL version="4.1.6-azure"
LABEL description="Minimac4 v4.1.6 with Azure compatibility features"

# Set timezone for Azure
ENV TZ=UTC

# Install runtime dependencies (same as basic minimac4 + Azure tools)
RUN apk add --no-cache \
    bash \
    ca-certificates \
    wget \
    libgomp \
    zlib \
    bzip2 \
    xz \
    curl \
    libstdc++ \
    # Azure tools
    python3 \
    py3-pip

# Install minimal Python packages for Azure
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    pandas

# Download and install prebuilt Minimac4 binary (reliable approach)
ENV MINIMAC4_VERSION=4.1.6
RUN cd /tmp && \
    wget https://github.com/statgen/Minimac4/releases/download/v${MINIMAC4_VERSION}/minimac4-${MINIMAC4_VERSION}-Linux-x86_64.sh && \
    chmod +x minimac4-${MINIMAC4_VERSION}-Linux-x86_64.sh && \
    ./minimac4-${MINIMAC4_VERSION}-Linux-x86_64.sh --skip-license --prefix=/usr/local && \
    rm minimac4-${MINIMAC4_VERSION}-Linux-x86_64.sh

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

# Create Azure test script  
RUN echo '#!/bin/bash' > /usr/local/bin/test-azure-tools && \
    echo 'echo "=== Minimac4-All Azure Test ==="' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🏠 Container: $(hostname)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "👤 User: $(whoami)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "📁 Working directory: $(pwd)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🧬 Minimac4: $(minimac4 --version 2>&1 | head -1)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "🐍 Python: $(python3 --version)"' >> /usr/local/bin/test-azure-tools && \
    echo 'echo "✅ Azure-compatible Minimac4 ready!"' >> /usr/local/bin/test-azure-tools && \
    chmod +x /usr/local/bin/test-azure-tools

# Set Azure-compatible defaults
WORKDIR /data
USER impuser

# Set default command
CMD ["/bin/bash", "-l"] 