FROM ubuntu:22.04

LABEL maintainer="ChiPImputation Pipeline"
LABEL description="LaTeX and Python container for comprehensive report generation"
LABEL version="1.0.0"

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies and TeX Live
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    # Basic utilities
    wget \
    curl \
    ca-certificates \
    gnupg \
    lsb-release \
    git \
    make \
    # Python and pip
    python3.10 \
    python3-pip \
    python3-dev \
    # Build essentials for Python packages
    build-essential \
    gcc \
    g++ \
    # Image processing libraries
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    # PDF manipulation
    ghostscript \
    poppler-utils \
    # TeX Live full installation for comprehensive LaTeX support
    texlive-full \
    # Additional LaTeX packages
    texlive-latex-extra \
    texlive-fonts-extra \
    texlive-science \
    texlive-publishers \
    texlive-xetex \
    latexmk \
    # Clean up
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and setuptools
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python scientific packages
RUN python3 -m pip install --no-cache-dir \
    numpy==1.24.3 \
    pandas==2.0.3 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    scipy==1.11.1 \
    scikit-learn==1.3.0 \
    plotly==5.15.0 \
    bokeh==3.2.1 \
    statsmodels==0.14.0 \
    jupyter==1.0.0 \
    notebook==7.0.0

# Install bioinformatics packages
RUN python3 -m pip install --no-cache-dir \
    biopython==1.81 \
    pysam==0.21.0 \
    cyvcf2==0.30.22 \
    pyvcf==0.6.8

# Install reporting and documentation packages
RUN python3 -m pip install --no-cache-dir \
    reportlab==4.0.4 \
    pypdf==3.12.2 \
    Pillow==10.0.0 \
    jinja2==3.1.2 \
    markdown==3.4.3 \
    python-docx==0.8.11 \
    xlsxwriter==3.1.2 \
    openpyxl==3.1.2

# Install additional visualization tools
RUN python3 -m pip install --no-cache-dir \
    plotnine==0.12.2 \
    altair==5.0.1 \
    pygal==3.0.0 \
    wordcloud==1.9.2

# Create working directory
WORKDIR /work

# Add useful LaTeX templates directory
RUN mkdir -p /usr/local/share/latex-templates

# Create a LaTeX template for reports
RUN cat > /usr/local/share/latex-templates/chipimputation.cls << 'EOF'
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{chipimputation}[2025/01/01 ChiPImputation Report Class]

% Base class
\LoadClass[11pt,a4paper]{article}

% Required packages
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{geometry}
\RequirePackage{graphicx}
\RequirePackage{xcolor}
\RequirePackage{hyperref}
\RequirePackage{fancyhdr}
\RequirePackage{booktabs}
\RequirePackage{longtable}
\RequirePackage{float}
\RequirePackage{caption}
\RequirePackage{subcaption}
\RequirePackage{amsmath,amssymb}
\RequirePackage{tikz}
\RequirePackage{pgfplots}

% Page geometry
\geometry{
    left=25mm,
    right=25mm,
    top=30mm,
    bottom=30mm
}

% Color definitions
\definecolor{chipblue}{RGB}{44, 62, 80}
\definecolor{chipgreen}{RGB}{39, 174, 96}
\definecolor{chiporange}{RGB}{230, 126, 34}
\definecolor{chipred}{RGB}{192, 57, 43}

% Header and footer
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small ChiPImputation Report}
\fancyhead[R]{\small \today}
\fancyfoot[C]{\small Page \thepage}

% Hyperref setup
\hypersetup{
    colorlinks=true,
    linkcolor=chipblue,
    filecolor=chipblue,
    urlcolor=chipblue,
    citecolor=chipblue
}
EOF

# Create helper scripts
RUN cat > /usr/local/bin/compile-report << 'EOF'
#!/bin/bash
# Helper script to compile LaTeX reports with proper error handling

if [ $# -eq 0 ]; then
    echo "Usage: compile-report <tex-file>"
    exit 1
fi

TEX_FILE="$1"
BASE_NAME="${TEX_FILE%.tex}"

echo "Compiling LaTeX report: $TEX_FILE"

# Run pdflatex multiple times for references
for i in 1 2 3; do
    echo "Pass $i..."
    pdflatex -interaction=nonstopmode -file-line-error "$TEX_FILE" || {
        echo "LaTeX compilation failed on pass $i"
        cat "${BASE_NAME}.log" | grep -A 5 "^!"
        exit 1
    }
done

# Run bibtex if .bib file exists
if [ -f "${BASE_NAME}.bib" ]; then
    echo "Running BibTeX..."
    bibtex "$BASE_NAME"
    pdflatex -interaction=nonstopmode "$TEX_FILE"
    pdflatex -interaction=nonstopmode "$TEX_FILE"
fi

# Clean up auxiliary files
echo "Cleaning up..."
rm -f "${BASE_NAME}.aux" "${BASE_NAME}.log" "${BASE_NAME}.out" \
      "${BASE_NAME}.toc" "${BASE_NAME}.lof" "${BASE_NAME}.lot" \
      "${BASE_NAME}.bbl" "${BASE_NAME}.blg" "${BASE_NAME}.nav" \
      "${BASE_NAME}.snm" "${BASE_NAME}.vrb"

echo "Report compiled successfully: ${BASE_NAME}.pdf"
EOF

RUN chmod +x /usr/local/bin/compile-report

# Create a Python script for testing
RUN cat > /usr/local/bin/test-environment.py << 'EOF'
#!/usr/bin/env python3
import sys
print("Python environment test:")
print(f"Python version: {sys.version}")

try:
    import numpy as np
    print(f"NumPy: {np.__version__}")
except ImportError:
    print("NumPy: NOT INSTALLED")

try:
    import pandas as pd
    print(f"Pandas: {pd.__version__}")
except ImportError:
    print("Pandas: NOT INSTALLED")

try:
    import matplotlib
    print(f"Matplotlib: {matplotlib.__version__}")
except ImportError:
    print("Matplotlib: NOT INSTALLED")

try:
    import seaborn as sns
    print(f"Seaborn: {sns.__version__}")
except ImportError:
    print("Seaborn: NOT INSTALLED")

try:
    import reportlab
    print(f"ReportLab: {reportlab.__version__}")
except ImportError:
    print("ReportLab: NOT INSTALLED")

import subprocess
try:
    result = subprocess.run(['pdflatex', '--version'], capture_output=True, text=True)
    print("LaTeX: " + result.stdout.split('\n')[0])
except:
    print("LaTeX: NOT AVAILABLE")

print("\nEnvironment test complete!")
EOF

RUN chmod +x /usr/local/bin/test-environment.py

# Set Python3 as default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# Test the environment
RUN python /usr/local/bin/test-environment.py

# Default command
CMD ["/bin/bash"]