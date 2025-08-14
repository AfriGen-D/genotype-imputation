#!/usr/bin/env python3
"""
R-squared - Position - Plot
Plots imputation accuracy (R-squared) vs genomic position
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gzip

# Set style for professional-looking plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Parse arguments (these come from Nextflow template variables)
info_file = "${info}"
target_file = "${target}"
output_file = "${output}"
maf_thresh = float("${maf_thresh}") if "${maf_thresh}" != "" else 0

# Helper function to read potentially gzipped files
def read_file(filename):
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

# Read in the files
info = read_file(info_file)
info = info[['SNP', 'Rsq', 'Genotyped', 'MAF', 'ALT_Frq']]

frq = read_file(target_file)
frq = frq[['SNP', 'CHR', 'POS']]

# Modify SNP ID of frq file to match info file format
# Split SNP column if it's in the format CHR_POS_REF_ALT
if '_' in frq['SNP'].iloc[0]:
    frq[['CHR_tmp', 'Position', 'REF', 'ALT']] = frq['SNP'].str.split('_', expand=True)
    frq['SNP'] = frq['CHR'].astype(str) + ':' + frq['POS'].astype(str) + ':' + frq['REF'] + ':' + frq['ALT']

# Merge tables
full = pd.merge(frq, info, on='SNP', how='inner')

# Check if merge produced any results
if len(full) == 0:
    print(f"Warning: No matching SNPs found between frequency file and info file.")
    print(f"  Frequency file has {len(frq)} SNPs")
    print(f"  Info file has {len(info)} SNPs")
    if len(frq) > 0 and len(info) > 0:
        print(f"  Sample SNP IDs from frequency: {frq['SNP'].head(3).tolist()}")
        print(f"  Sample SNP IDs from info: {info['SNP'].head(3).tolist()}")

# Change chromosome names for display
full['CHR'] = 'Chromosome ' + full['CHR'].astype(str)

# Extract imputed SNPs and filter out SNPs with missing R-squared values
imputed = full[full['Genotyped'] == 'Imputed'].copy()
imputed = imputed[(imputed['Rsq'] != '-') & (imputed['Rsq'].notna())]
imputed['Rsq'] = pd.to_numeric(imputed['Rsq'], errors='coerce')

genotyped = full[full['Genotyped'] == 'Genotyped'].copy()

# Display only ~50,000 Imputed SNPs in the r2 - position plot
if len(imputed) > 50000:
    N = len(imputed) // 50000
    imputed = imputed.iloc[::N, :]

# Display only ~1,000 Genotyped SNPs as ticks
if len(genotyped) > 1000:
    N2 = len(genotyped) // 1000
    genotyped = genotyped.iloc[::N2, :]

# Filter out MAF below a given value
AF_thresh = maf_thresh / 100 if maf_thresh > 0 else 0
imputed['MAF'] = pd.to_numeric(imputed['MAF'], errors='coerce')
imputed = imputed[(imputed['MAF'] > AF_thresh) & (imputed['MAF'] != 0)]

# Categorize MAF levels
def categorize_maf(maf):
    if 0 < maf <= 0.001:
        return "extreme rare (0,0.001]"
    elif 0.001 < maf <= 0.01:
        return "moderate rare (0.001,0.01]"
    elif 0.01 < maf <= 0.02:
        return "rare (0.01,0.02]"
    elif 0.02 < maf <= 0.05:
        return "moderate (0.02,0.05]"
    elif 0.05 < maf <= 0.2:
        return "common (0.05,0.2]"
    elif 0.2 < maf <= 0.5:
        return "extreme common (0.2,0.5]"
    else:
        return "unknown"

imputed['MAF_category'] = imputed['MAF'].apply(categorize_maf)

# Get unique chromosomes
chromosomes = sorted(imputed['CHR'].unique())

# Check if there's data to plot
if len(imputed) == 0 or len(chromosomes) == 0:
    print("Warning: No imputed data found after filtering. Creating empty plot.")
    fig, ax = plt.subplots(1, 1, figsize=(10, 5.5))
    ax.text(0.5, 0.5, 'No imputed data available\\nafter filtering', 
            ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.set_xlabel('Position [bp]', fontsize=10)
    ax.set_ylabel('Imputation accuracy (r-squared)', fontsize=10)
    ax.set_title('No Data', fontsize=10)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    import sys
    sys.exit(0)

# Create the plot
n_chr = len(chromosomes)
fig, axes = plt.subplots(1, n_chr, figsize=(10, 5.5), sharey=True)

if n_chr == 1:
    axes = [axes]

# Define colors for MAF categories
maf_colors = {
    "extreme rare (0,0.001]": "#1f77b4",
    "moderate rare (0.001,0.01]": "#ff7f0e",
    "rare (0.01,0.02]": "#2ca02c",
    "moderate (0.02,0.05]": "#d62728",
    "common (0.05,0.2]": "#9467bd",
    "extreme common (0.2,0.5]": "#8c564b"
}

# Plot for each chromosome
for idx, (ax, chr_name) in enumerate(zip(axes, chromosomes)):
    chr_data = imputed[imputed['CHR'] == chr_name]
    chr_genotyped = genotyped[genotyped['CHR'] == chr_name]
    
    # Plot points colored by MAF category
    for category, color in maf_colors.items():
        cat_data = chr_data[chr_data['MAF_category'] == category]
        if not cat_data.empty:
            ax.scatter(cat_data['POS'], cat_data['Rsq'], 
                      c=color, s=1, alpha=0.6, label=category)
    
    # Add horizontal line at R-squared = 0.3
    ax.axhline(y=0.3, color='red', linestyle='--', alpha=0.5, linewidth=0.5)
    
    # Add rug plot for genotyped SNPs
    if not chr_genotyped.empty:
        ax.scatter(chr_genotyped['POS'], np.zeros(len(chr_genotyped)), 
                  marker='|', s=10, c='black', alpha=0.3)
    
    # Formatting
    ax.set_xlabel('Position [bp]' if idx == n_chr//2 else '', fontsize=10)
    ax.set_title(chr_name, fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.grid(True, alpha=0.3)
    
    if idx == 0:
        ax.set_ylabel('Imputation accuracy (r-squared)', fontsize=10)

# Add legend
handles, labels = axes[0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), 
              title='MAF Category', fontsize=8)

plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()