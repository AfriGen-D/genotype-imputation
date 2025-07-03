############################################################
# Dockerfile to build IMPUTE5 genotype imputation tool
# Based on Alpine Linux for better stability
############################################################

# Set the base image to Alpine Linux
FROM alpine:3.18

# File Author / Maintainer
LABEL maintainer="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"
LABEL description="IMPUTE5 genotype imputation tool container"
LABEL version="1.1.5"

################## BEGIN INSTALLATION ######################

# Install runtime dependencies
RUN apk add --no-cache \
    bash \
    wget \
    ca-certificates \
    libgomp \
    libstdc++ \
    curl

# Create working directory
WORKDIR /data

# Note: IMPUTE5 is commercial software requiring licensing
# Users should download the appropriate binary for their system
# from the official Oxford website and mount it into the container
# 
# Example usage:
# docker run -v /path/to/impute5:/usr/local/bin/impute5 impute5:latest
#
# Alternative: Add your licensed binary to the image during build

# Create a placeholder script
RUN echo '#!/bin/bash' > /usr/local/bin/impute5 && \
    echo 'echo "IMPUTE5 is commercial software requiring licensing."' >> /usr/local/bin/impute5 && \
    echo 'echo "Please download the licensed binary from Oxford and mount it:"' >> /usr/local/bin/impute5 && \
    echo 'echo "docker run -v /path/to/your/impute5:/usr/local/bin/impute5 impute5:latest"' >> /usr/local/bin/impute5 && \
    chmod +x /usr/local/bin/impute5

# Set default command
CMD ["impute5", "--help"]
