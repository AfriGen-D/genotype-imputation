#!/usr/bin/env python3
"""
Dosage Distribution Plot
Shows certainty of imputed genotypes (dosages closer to 0, 1, or 2 indicate higher certainty)
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
fig, ax = plt.subplots(figsize=(8, 6))

# Check for dosage-related columns
if 'AvgCall' in data.columns:
    # AvgCall represents average call score (dosage certainty)
    avg_call = pd.to_numeric(data['AvgCall'], errors='coerce').dropna()
    
    # Plot histogram
    ax.hist(avg_call, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax.set_xlabel('Average Call Score', fontsize=11)
    ax.set_ylabel('Number of Variants', fontsize=11)
    ax.set_title('Dosage Certainty Distribution', fontsize=13)
    
    # Add vertical lines for thresholds
    ax.axvline(x=0.9, color='g', linestyle='--', alpha=0.7, label='High certainty (>0.9)')
    ax.axvline(x=0.7, color='orange', linestyle='--', alpha=0.7, label='Moderate certainty (>0.7)')
    ax.axvline(x=0.5, color='r', linestyle='--', alpha=0.7, label='Low certainty (<0.5)')
    
    # Add statistics
    mean_score = avg_call.mean()
    median_score = avg_call.median()
    high_cert = (avg_call > 0.9).sum()
    total = len(avg_call)
    
    stats_text = f'Mean: {mean_score:.3f}\nMedian: {median_score:.3f}\n'
    stats_text += f'High certainty: {high_cert:,} ({100*high_cert/total:.1f}%)\n'
    stats_text += f'Total variants: {total:,}'
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
elif 'Rsq' in data.columns:
    # Use R-squared as proxy for dosage certainty
    rsq = pd.to_numeric(data['Rsq'].replace('-', np.nan), errors='coerce').dropna()
    
    # Create bins for certainty levels
    certainty_bins = [0, 0.3, 0.5, 0.7, 0.9, 1.0]
    certainty_labels = ['Very Low\n(0-0.3)', 'Low\n(0.3-0.5)', 'Moderate\n(0.5-0.7)', 
                       'High\n(0.7-0.9)', 'Very High\n(0.9-1.0)']
    
    rsq_binned = pd.cut(rsq, bins=certainty_bins, labels=certainty_labels)
    counts = rsq_binned.value_counts().sort_index()
    
    # Plot bar chart
    colors = ['#d62728', '#ff7f0e', '#ffbb78', '#98df8a', '#2ca02c']
    bars = ax.bar(range(len(counts)), counts.values, color=colors, edgecolor='black', alpha=0.7)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=0)
    ax.set_xlabel('Imputation Certainty (R²)', fontsize=11)
    ax.set_ylabel('Number of Variants', fontsize=11)
    ax.set_title('Imputation Certainty Distribution', fontsize=13)
    
    # Add value labels on bars
    for bar, val in zip(bars, counts.values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:,}\n({100*val/len(rsq):.1f}%)',
                ha='center', va='bottom', fontsize=9)
    
    # Add total count
    ax.text(0.98, 0.98, f'Total variants: {len(rsq):,}', transform=ax.transAxes, 
            fontsize=10, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

else:
    ax.text(0.5, 0.5, 'Dosage certainty data not available', 
            ha='center', va='center', fontsize=12)

ax.legend(loc='best')
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Dosage distribution plot saved to {output_file}")