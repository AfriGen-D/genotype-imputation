#!/usr/bin/env python3
"""
SNP counts vs Mean Imputation Quality Score - Histogram
Creates histograms of imputation quality scores for each reference panel
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip
import seaborn as sns

# Set style for professional-looking plots
plt.style.use('seaborn-v0_8-whitegrid')

# Parse arguments from Nextflow template
infos = "${infos}"
ref_panels = "${ref_panels}"
plot_out = "${plot_out}"
impute_info_cutoff = float("${impute_info_cutoff}")

# Helper function to read potentially gzipped files
def read_file(filename):
    try:
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                return pd.read_csv(f, sep='\t')
        else:
            return pd.read_csv(filename, sep='\t')
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return pd.DataFrame()

# Parse input file lists
info_files = infos.split(',')
panel_names = ref_panels.split(',')

# Read and combine info files from all reference panels
all_data = []
for panel_name, info_file in zip(panel_names, info_files):
    if info_file and info_file != '':
        panel_data = read_file(info_file)
        if not panel_data.empty:
            # Select relevant columns
            panel_data = panel_data[['SNP', 'MAF', 'Rsq', 'Genotyped']]
            panel_data['R_Panel'] = panel_name
            all_data.append(panel_data)

# Combine all panels
if all_data:
    full = pd.concat(all_data, ignore_index=True)
else:
    print("No data to plot")
    # Create empty plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
    plt.savefig(plot_out, dpi=150, bbox_inches='tight')
    exit()

# Filter for imputed SNPs and remove missing values
imputed = full[full['Genotyped'] == 'Imputed'].copy()
imputed = imputed[(imputed['Rsq'] != '-') & (imputed['Rsq'].notna())]

# Convert to numeric
imputed['Rsq'] = pd.to_numeric(imputed['Rsq'], errors='coerce')

# Remove any remaining NaN values
imputed = imputed.dropna(subset=['Rsq'])

# Determine number of panels
n_panels = len(panel_names)

# Create subplots
fig, axes = plt.subplots(1, n_panels, figsize=(4*n_panels, 5), sharey=True)

# If only one panel, make axes a list
if n_panels == 1:
    axes = [axes]

# Use a color palette
colors = sns.color_palette("husl", n_panels)

# Create histogram for each reference panel
for idx, (panel_name, ax) in enumerate(zip(panel_names, axes)):
    panel_data = imputed[imputed['R_Panel'] == panel_name]
    
    if not panel_data.empty:
        # Create histogram
        ax.hist(panel_data['Rsq'], bins=21, range=(0, 1), 
                color=colors[idx], alpha=0.7, edgecolor='black', linewidth=0.5)
        
        # Add vertical line at cutoff
        ax.axvline(x=impute_info_cutoff, color='red', linestyle='--', 
                  alpha=0.7, linewidth=1.5)
        
        # Set title and labels
        ax.set_title(panel_name, fontsize=11)
        ax.set_xlabel('Imputation Quality Score' if idx == n_panels//2 else '', fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1.1, 0.2))
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        
        # Style improvements
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if idx == 0:
            ax.set_ylabel('SNP Count', fontsize=11)
            # Format y-axis with comma separator
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

# Add overall title
fig.suptitle('Distribution of Imputation Quality Scores', fontsize=13, y=1.02)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig(plot_out, dpi=150, bbox_inches='tight')
plt.close()