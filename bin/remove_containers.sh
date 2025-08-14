#!/bin/bash
# Script to remove hardcoded container definitions from process files

echo "Removing hardcoded container definitions from process files..."

# Find all .nf files in modules/local and remove container lines
find /users/mamana/genotype-imputation/modules/local -name "*.nf" -type f | while read file; do
    # Check if file has a container definition
    if grep -q "^    container '" "$file"; then
        echo "Processing: $file"
        # Remove the container line
        sed -i "/^    container '/d" "$file"
    fi
done

echo "Done! All container definitions have been removed from process files."
echo "Containers are now defined in conf/modules.config"