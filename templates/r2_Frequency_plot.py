#!/usr/bin/env python3
"""
SNP counts vs Mean Imputation Quality Score
Plots the frequency of SNPs at different imputation quality levels
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

# Round Rsq to 2 decimal places for binning
imputed['Rsq_bin'] = np.round(imputed['Rsq'], 2)

# Calculate mean Rsq and count for each bin and panel
grouped = imputed.groupby(['R_Panel', 'Rsq_bin']).agg(
    Rsq_mean=('Rsq_bin', 'mean'),
    N=('Rsq_bin', 'count')
).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(8, 5))

# Use a color palette
colors = sns.color_palette("husl", len(panel_names))

# Plot for each reference panel
for idx, panel_name in enumerate(panel_names):
    panel_data = grouped[grouped['R_Panel'] == panel_name]
    if not panel_data.empty:
        # Sort by Rsq_mean for proper line connection
        panel_data = panel_data.sort_values('Rsq_mean')
        ax.plot(panel_data['Rsq_mean'], panel_data['N'], 
               marker='o', markersize=4, label=panel_name,
               color=colors[idx], linewidth=1.5, alpha=0.8)

# Add vertical line at cutoff
ax.axvline(x=impute_info_cutoff, color='red', linestyle='--', alpha=0.5, linewidth=1)

# Formatting
ax.set_xlabel('Mean Imputation Quality Score', fontsize=11)
ax.set_ylabel('SNP Count', fontsize=11)
ax.set_xlim(0, 1)
ax.set_xticks(np.arange(0, 1.1, 0.2))

# Format y-axis with comma separator for thousands
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))

# Add grid
ax.grid(True, alpha=0.3)

# Style improvements
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add legend
ax.legend(title='Reference Panel', loc='best', frameon=True, fancybox=True)

plt.title('SNP Count vs Imputation Quality', fontsize=12, pad=15)
plt.tight_layout()

# Save the plot
plt.savefig(plot_out, dpi=150, bbox_inches='tight')
plt.close()