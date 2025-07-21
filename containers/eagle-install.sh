#!/bin/bash
# Script to install Eagle in a container that has VCF processing tools

set -euo pipefail

echo "Installing Eagle 2.4.1..."

# Download Eagle
cd /tmp
wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
    https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz

# Extract Eagle
tar -xzf Eagle_v2.4.1.tar.gz

# Copy Eagle binary to /usr/local/bin
cp Eagle_v2.4.1/eagle /usr/local/bin/eagle
chmod +x /usr/local/bin/eagle

# Clean up
rm -rf Eagle_v2.4.1*

echo "Eagle 2.4.1 installed successfully!"
eagle --version 