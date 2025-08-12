#!/usr/bin/env python3
"""
Calculate average R-squared value from imputation info file
"""

import pandas as pd
import numpy as np
import gzip

# Parse arguments from Nextflow template
file_infos = "${file_infos}"
meanr2_out = "${meanr2_out}"

# Helper function to read potentially gzipped files
def read_file(filename):
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

# Read the info file
info_file = read_file(file_infos)

# Filter out missing values and calculate mean
filtered_data = info_file[info_file['Rsq'] != '-'].copy()
filtered_data['Rsq'] = pd.to_numeric(filtered_data['Rsq'], errors='coerce')
filtered_data = filtered_data.dropna(subset=['Rsq'])

# Calculate average R-squared
average = filtered_data['Rsq'].mean()

# Write the result to output file
with open(meanr2_out, 'w') as f:
    f.write(str(average) + "\\n")