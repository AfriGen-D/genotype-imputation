# Use official Python base image with slim variant for smaller size
FROM python:3.11-slim-bullseye

# Set metadata
LABEL maintainer="mamana@cbio"
LABEL description="Python container for genotype imputation QC plotting"
LABEL version="1.1.0"

# Set working directory
WORKDIR /app

# Install system dependencies required for scientific Python packages
# Including wkhtmltopdf for PDF generation
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libfreetype6-dev \
    libpng-dev \
    pkg-config \
    ca-certificates \
    curl \
    procps \
    wkhtmltopdf \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file
COPY requirements.txt /app/

# Upgrade pip and install Python packages
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Create a non-root user for running the container
RUN useradd -m -s /bin/bash plotuser && \
    chown -R plotuser:plotuser /app

# Set environment variables for matplotlib
ENV MPLBACKEND=Agg
ENV PYTHONUNBUFFERED=1

# Add /app to PYTHONPATH
ENV PYTHONPATH="/app:${PYTHONPATH}"

# Create directories for input/output
RUN mkdir -p /data /output && \
    chown -R plotuser:plotuser /data /output

# Switch to non-root user
USER plotuser

# Set default command
CMD ["python", "--version"]