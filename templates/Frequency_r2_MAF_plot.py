#!/usr/bin/env python3
"""
Imputation quality score (r2-values) and SNP count plotted against MAF bin
Generates a dual-axis scatterplot: x = MAF bin,
y-axis left: SNPs count, y-axis right: mean r-squared
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
imputed = imputed[(imputed['MAF'] != '-') & (imputed['MAF'].notna())]

# Convert to numeric
imputed['Rsq'] = pd.to_numeric(imputed['Rsq'], errors='coerce')
imputed['MAF'] = pd.to_numeric(imputed['MAF'], errors='coerce')

# Remove any remaining NaN values
imputed = imputed.dropna(subset=['Rsq', 'MAF'])

# Round MAF to 2 decimal places for binning
imputed['MAF_bin'] = np.round(imputed['MAF'], 2)

# Calculate mean Rsq and count for each MAF bin and panel
grouped = imputed.groupby(['R_Panel', 'MAF_bin']).agg(
    Rsq_mean=('Rsq', 'mean'),
    N=('MAF_bin', 'count')
).reset_index()

# Create the plot with dual y-axes
fig, ax1 = plt.subplots(figsize=(10, 6))

# Use a color palette
colors = sns.color_palette("husl", len(panel_names))

# Plot for each reference panel
for idx, panel_name in enumerate(panel_names):
    panel_data = grouped[grouped['R_Panel'] == panel_name]
    if not panel_data.empty:
        # Sort by MAF_bin for proper line connection
        panel_data = panel_data.sort_values('MAF_bin')
        
        # Plot SNP count on primary y-axis
        line1 = ax1.plot(panel_data['MAF_bin'], panel_data['N'], 
                        marker='s', markersize=6, label=f'{panel_name} (Count)',
                        color=colors[idx], linewidth=1.5, alpha=0.8)

# Create secondary y-axis
ax2 = ax1.twinx()

# Plot mean R-squared on secondary y-axis
for idx, panel_name in enumerate(panel_names):
    panel_data = grouped[grouped['R_Panel'] == panel_name]
    if not panel_data.empty:
        # Sort by MAF_bin for proper line connection
        panel_data = panel_data.sort_values('MAF_bin')
        
        # Plot mean R-squared on secondary y-axis
        line2 = ax2.plot(panel_data['MAF_bin'], panel_data['Rsq_mean'], 
                        marker='D', markersize=4, label=f'{panel_name} (R²)',
                        color=colors[idx], linewidth=1.5, alpha=0.8, linestyle='--')

# Set labels and formatting for primary y-axis
ax1.set_xlabel('Minor Allele Frequency (MAF) Bin', fontsize=11)
ax1.set_ylabel('SNP Count', fontsize=11, color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

# Set labels and formatting for secondary y-axis
ax2.set_ylabel('Mean Imputation Quality Score (R²)', fontsize=11, color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.set_ylim(0, 1)
ax2.set_yticks(np.arange(0, 1.1, 0.2))

# Set x-axis limits
ax1.set_xlim(0, 0.5)

# Add grid
ax1.grid(True, alpha=0.3)

# Get all lines and labels for legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Create combined legend
ax1.legend(lines1 + lines2, labels1 + labels2, 
          loc='upper right', frameon=True, fancybox=True, fontsize=9)

# Add title
plt.title('MAF vs SNP Count and Imputation Quality', fontsize=12, pad=15)

# Style improvements
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.tight_layout()

# Save the plot
plt.savefig(plot_out, dpi=150, bbox_inches='tight')
plt.close()