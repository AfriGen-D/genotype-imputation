#!/usr/bin/env python3
"""
INFO Score Distribution Plot
Histogram showing the distribution of INFO scores across all variants
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
output_file = "${output_file}"

# Read data
data = read_file(info_file)

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Get INFO score column (could be 'INFO', 'Rsq', or 'R2')
info_col = None
for col in ['INFO', 'Rsq', 'R2']:
    if col in data.columns:
        info_col = col
        break

if info_col:
    # Clean and convert data
    info_scores = pd.to_numeric(data[info_col].replace('-', np.nan), errors='coerce').dropna()
    
    # Plot 1: Histogram
    counts, bins, patches = ax1.hist(info_scores, bins=50, edgecolor='black', alpha=0.7)
    
    # Color code by quality
    for i, patch in enumerate(patches):
        if bins[i] < 0.3:
            patch.set_facecolor('#d62728')  # Red for poor
        elif bins[i] < 0.8:
            patch.set_facecolor('#ff7f0e')  # Orange for moderate
        else:
            patch.set_facecolor('#2ca02c')  # Green for good
    
    ax1.axvline(x=0.3, color='r', linestyle='--', alpha=0.7, linewidth=2, label='Poor quality threshold')
    ax1.axvline(x=0.8, color='g', linestyle='--', alpha=0.7, linewidth=2, label='High quality threshold')
    ax1.set_xlabel(f'{info_col} Score', fontsize=11)
    ax1.set_ylabel('Number of Variants', fontsize=11)
    ax1.set_title(f'{info_col} Score Distribution', fontsize=13)
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Calculate statistics
    mean_score = info_scores.mean()
    median_score = info_scores.median()
    q25 = info_scores.quantile(0.25)
    q75 = info_scores.quantile(0.75)
    
    poor_quality = (info_scores < 0.3).sum()
    moderate_quality = ((info_scores >= 0.3) & (info_scores < 0.8)).sum()
    high_quality = (info_scores >= 0.8).sum()
    total = len(info_scores)
    
    # Add statistics box
    stats_text = f'Statistics:\n'
    stats_text += f'Mean: {mean_score:.3f}\n'
    stats_text += f'Median: {median_score:.3f}\n'
    stats_text += f'Q25-Q75: {q25:.3f}-{q75:.3f}\n'
    stats_text += f'\nQuality Distribution:\n'
    stats_text += f'Poor (<0.3): {poor_quality:,} ({100*poor_quality/total:.1f}%)\n'
    stats_text += f'Moderate: {moderate_quality:,} ({100*moderate_quality/total:.1f}%)\n'
    stats_text += f'High (≥0.8): {high_quality:,} ({100*high_quality/total:.1f}%)\n'
    stats_text += f'\nTotal: {total:,} variants'
    
    ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # Plot 2: Cumulative distribution
    sorted_scores = np.sort(info_scores)
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    
    ax2.plot(sorted_scores, cumulative, linewidth=2, color='navy')
    ax2.fill_between(sorted_scores, 0, cumulative, alpha=0.3, color='skyblue')
    
    # Add reference lines
    ax2.axvline(x=0.3, color='r', linestyle='--', alpha=0.7, linewidth=1.5)
    ax2.axvline(x=0.8, color='g', linestyle='--', alpha=0.7, linewidth=1.5)
    ax2.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, linewidth=1)
    
    # Mark percentiles
    percentiles = [0.1, 0.25, 0.5, 0.75, 0.9]
    for p in percentiles:
        val = info_scores.quantile(p)
        idx = np.searchsorted(sorted_scores, val)
        if idx < len(sorted_scores):
            ax2.plot(val, p, 'ro', markersize=5)
            ax2.text(val, p, f' {int(p*100)}%', fontsize=8, ha='left')
    
    ax2.set_xlabel(f'{info_col} Score', fontsize=11)
    ax2.set_ylabel('Cumulative Proportion', fontsize=11)
    ax2.set_title('Cumulative Distribution', fontsize=13)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.grid(True, alpha=0.3)
    
    # Add proportion annotations
    prop_03 = (info_scores <= 0.3).mean()
    prop_08 = (info_scores <= 0.8).mean()
    ax2.text(0.3, 0.05, f'{prop_03:.1%} ≤ 0.3', fontsize=9, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.text(0.8, 0.05, f'{prop_08:.1%} ≤ 0.8', fontsize=9, ha='center',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

else:
    ax1.text(0.5, 0.5, 'INFO score data not available', 
            ha='center', va='center', fontsize=12)
    ax2.text(0.5, 0.5, 'INFO score data not available', 
            ha='center', va='center', fontsize=12)

plt.suptitle('INFO Score Quality Assessment', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"INFO score distribution plot saved to {output_file}")