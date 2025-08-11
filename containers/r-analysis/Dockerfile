# Dockerfile for R analysis container with required packages
FROM rocker/r-ver:4.3.2

# Install system dependencies including Python 3
RUN apt-get update && apt-get install -y \
    python3 \
    python3-minimal \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install R packages required by the pipeline
RUN R -e "install.packages(c('dplyr', 'data.table', 'ggplot2', 'tidyr', 'readr', 'stringr', 'purrr', 'tibble', 'R.utils'), repos='https://cloud.r-project.org/')"

# Install additional packages that might be needed for plotting
RUN R -e "install.packages(c('scales', 'gridExtra', 'cowplot', 'viridis', 'RColorBrewer'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /data

# Default command
CMD ["R"]