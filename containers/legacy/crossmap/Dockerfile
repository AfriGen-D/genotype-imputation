############################################################
# Dockerfile to build CrossMap coordinate conversion tool
# Based on Alpine Linux for better stability
############################################################

# Set the base image to Alpine Linux
FROM alpine:3.18

# File Author / Maintainer
LABEL maintainer="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"
LABEL description="CrossMap coordinate conversion tool - Latest version"
LABEL version="0.7.3"

################## BEGIN INSTALLATION ######################

# Install system dependencies
RUN apk add --no-cache \
    python3 \
    py3-pip \
    py3-setuptools \
    py3-wheel \
    bash \
    wget \
    ca-certificates \
    build-base \
    python3-dev \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    gcc \
    musl-dev

# Install Python dependencies step by step
RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install numpy
RUN pip3 install pysam
RUN pip3 install bx-python
RUN pip3 install CrossMap==0.7.3 --verbose

# List what was installed
RUN pip3 list | grep -i cross

# Create working directory
WORKDIR /data

# Create a wrapper script that points to the actual CrossMap executable
RUN echo '#!/bin/bash' > /usr/local/bin/crossmap && \
    echo '/usr/bin/CrossMap "$@"' >> /usr/local/bin/crossmap && \
    chmod +x /usr/local/bin/crossmap

# Set default command
CMD ["/usr/bin/CrossMap"]