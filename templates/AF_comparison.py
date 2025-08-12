#!/usr/bin/env python3
"""
Allele Frequency Comparison Plot
Correlation between reference allele freq and target allele freq
Creates scatter plot colored based on the SNP r-squared values
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip

# Parse arguments from Nextflow template
info_file = "${info}"
target_file = "${target}"
frq_file = "${frq}"
output_color = "${outputcolor}"
rsq_threshold = 0  # Default value, can be modified if needed
subset_size = 20000  # Number of SNPs to display

# Helper function to read potentially gzipped files
def read_file(filename):
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

# Read input files
info = read_file(info_file)
frq = read_file(target_file)
ref_frq = read_file(frq_file)

# Modify SNP ID format if needed
if '_' in frq['SNP'].iloc[0]:
    frq[['CHR_tmp', 'Position', 'REF', 'ALT']] = frq['SNP'].str.split('_', expand=True)
    frq['SNP'] = frq['CHR'].astype(str) + ':' + frq['POS'].astype(str) + ':' + frq['REF'] + ':' + frq['ALT']

# Select relevant columns
frq = frq[['CHR', 'POS', 'SNP']]

# Merge tables
full = pd.merge(frq, info[['SNP', 'ALT_Frq', 'Rsq', 'Genotyped']], on='SNP', how='inner')
full = pd.merge(full, ref_frq, on=['CHR', 'POS'], how='inner')

# Filter for imputed SNPs and calculate frequency difference
imputed = full[full['Genotyped'] == 'Imputed'].copy()
imputed['diff'] = np.abs(imputed['ALT_Frq'].astype(float) - imputed['AF'].astype(float))

# Remove SNPs with missing R-squared values
imputed = imputed[(imputed['Rsq'] != '-') & (imputed['Rsq'].notna())]
imputed['Rsq'] = pd.to_numeric(imputed['Rsq'], errors='coerce')
imputed['ALT_Frq'] = pd.to_numeric(imputed['ALT_Frq'], errors='coerce')
imputed['AF'] = pd.to_numeric(imputed['AF'], errors='coerce')

# Apply R-squared threshold if specified
if rsq_threshold > 0:
    imputed = imputed[imputed['Rsq'] > rsq_threshold]

# Subset data if too many SNPs
if len(imputed) > subset_size:
    n = len(imputed) // subset_size
    imputed = imputed.iloc[::n, :]

# Create the plot colored by R-squared values
fig, ax = plt.subplots(figsize=(7, 7))

# Create scatter plot with color based on R-squared
scatter = ax.scatter(imputed['ALT_Frq'], imputed['AF'], 
                    c=imputed['Rsq'], cmap='coolwarm',
                    s=3, alpha=0.6, edgecolors='none')

# Add diagonal reference line
ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=0.5)

# Add color bar
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('R-squared', rotation=270, labelpad=15)

# Set axis labels and limits
ax.set_xlabel('Ref Allele Frequency (Uploaded Samples)', fontsize=11)
ax.set_ylabel('Ref Allele Frequency (Reference Panel)', fontsize=11)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xticks(np.arange(0, 1.1, 0.2))
ax.set_yticks(np.arange(0, 1.1, 0.2))

# Add text annotations
ax.text(0.2, 1.02, f'R-squared threshold: {rsq_threshold}', 
        transform=ax.transAxes, fontsize=10)
ax.text(0.2, 0.95, f'{len(imputed):,} SNPs', 
        transform=ax.transAxes, fontsize=10)

# Add grid for better readability
ax.grid(True, alpha=0.3)

# Style improvements
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('Allele Frequency Comparison', fontsize=12, pad=20)
plt.tight_layout()

# Save the plot
plt.savefig(output_color, dpi=150, bbox_inches='tight')
plt.close()