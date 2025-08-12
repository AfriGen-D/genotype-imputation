#!/usr/bin/env python3
"""
Concordance Rate vs MAF Plot
For masked/known variants, shows concordance between true and imputed genotypes across MAF bins
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip

def read_file(filename):
    """Read potentially gzipped files"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

# Parse arguments from Nextflow template
info_file = "${info_file}"
masked_file = "${masked_file}"  # Optional file with masked variant concordance data
output_file = "${output_file}"

# Read data
data = read_file(info_file)

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Check if we have empirical R-squared (concordance proxy) or actual concordance data
has_concordance = False

if masked_file and masked_file != "":
    # Read masked variant data if available
    try:
        masked_data = read_file(masked_file)
        if 'concordance' in masked_data.columns and 'MAF' in masked_data.columns:
            has_concordance = True
            concordance_data = masked_data
    except:
        pass

if not has_concordance and 'EmpRsq' in data.columns and 'MAF' in data.columns:
    # Use empirical R-squared as concordance proxy
    concordance_data = data[data['EmpRsq'] != '-'].copy()
    concordance_data['concordance'] = pd.to_numeric(concordance_data['EmpRsq'], errors='coerce')
    concordance_data['MAF'] = pd.to_numeric(concordance_data['MAF'], errors='coerce')
    concordance_data = concordance_data.dropna(subset=['concordance', 'MAF'])
    has_concordance = True

if has_concordance:
    # Define MAF bins
    maf_bins = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
    maf_labels = ['Rare\n(0-0.1%)', 'Very Low\n(0.1-1%)', 'Low\n(1-5%)', 
                  'Medium\n(5-10%)', 'Common\n(10-20%)', 'Very Common\n(20-50%)']
    
    concordance_data['MAF_bin'] = pd.cut(concordance_data['MAF'], bins=maf_bins, labels=maf_labels)
    
    # Plot 1: Bar plot with error bars
    grouped = concordance_data.groupby('MAF_bin')['concordance'].agg(['mean', 'std', 'count', 'min', 'max'])
    grouped = grouped.dropna()
    
    x = range(len(grouped))
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(grouped)))
    
    bars = ax1.bar(x, grouped['mean'], yerr=grouped['std'], capsize=5, 
                   color=colors, edgecolor='black', alpha=0.7)
    
    ax1.set_xticks(x)
    ax1.set_xticklabels(grouped.index, rotation=45, ha='right')
    ax1.set_xlabel('MAF Bin', fontsize=11)
    ax1.set_ylabel('Concordance Rate', fontsize=11)
    ax1.set_title('Concordance by MAF Category', fontsize=13)
    ax1.set_ylim(0, 1.05)
    
    # Add horizontal lines
    ax1.axhline(y=0.95, color='g', linestyle='--', alpha=0.5, label='95% concordance')
    ax1.axhline(y=0.90, color='orange', linestyle='--', alpha=0.5, label='90% concordance')
    ax1.axhline(y=0.80, color='r', linestyle='--', alpha=0.5, label='80% concordance')
    ax1.legend(loc='lower right', fontsize=9)
    
    # Add sample sizes and values on bars
    for i, (bar, (idx, row)) in enumerate(zip(bars, grouped.iterrows())):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{row["mean"]:.3f}\n(n={int(row["count"])})',
                ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Scatter plot with trend
    ax2.scatter(concordance_data['MAF'], concordance_data['concordance'], 
               alpha=0.3, s=10, c='blue', edgecolors='none')
    
    # Add binned means
    maf_ranges = np.linspace(0, 0.5, 21)
    bin_centers = []
    bin_means = []
    for i in range(len(maf_ranges)-1):
        mask = (concordance_data['MAF'] >= maf_ranges[i]) & (concordance_data['MAF'] < maf_ranges[i+1])
        if mask.sum() > 0:
            bin_centers.append((maf_ranges[i] + maf_ranges[i+1])/2)
            bin_means.append(concordance_data.loc[mask, 'concordance'].mean())
    
    if bin_centers:
        ax2.plot(bin_centers, bin_means, 'r-', linewidth=2, label='Binned mean')
    
    ax2.set_xlabel('Minor Allele Frequency', fontsize=11)
    ax2.set_ylabel('Concordance Rate', fontsize=11)
    ax2.set_title('Concordance vs MAF (Continuous)', fontsize=13)
    ax2.set_xlim(-0.01, 0.51)
    ax2.set_ylim(0, 1.05)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Add statistics
    overall_concordance = concordance_data['concordance'].mean()
    ax2.text(0.98, 0.02, f'Overall concordance: {overall_concordance:.3f}\nN = {len(concordance_data):,}',
            transform=ax2.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

else:
    # No concordance data available - show informative message
    for ax in [ax1, ax2]:
        ax.text(0.5, 0.5, 'Concordance data not available\n\n' +
                'To generate concordance metrics:\n' +
                '1. Mask known variants before imputation\n' +
                '2. Compare imputed vs true genotypes\n' +
                '3. Calculate concordance rates',
                ha='center', va='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_title('Concordance Analysis', fontsize=13)

plt.suptitle('Imputation Concordance by MAF', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Concordance vs MAF plot saved to {output_file}")