#!/usr/bin/env python3
"""
Plot results by MAF bins
Plots accuracy/concordance by MAF from TSV summary table
Version: 2.0 - Fixed TSV format handling
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Set style for professional-looking plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Define color palette similar to the R script
superpop_plot_colours = [
    "#8A2BE2", "#CC0066", "#6495ED", "#66CDAA", "#A52A2A", 
    "#CDAA7D", "#66CD00", "#7AC5CD", "#CD5B45", "#CDC8B1", 
    "#00CDCD", "#CD950C", "#8B7500", "#800000", "#808000", 
    "#008000", "#800080", "#008080", "#000080"
]

# Read the data - the file should be a TSV with columns from the summary report
try:
    data = pd.read_csv("${report}", sep='\t')
except Exception as e:
    print(f"Error reading file: {e}", file=sys.stderr)
    # Create empty plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, 'Error reading data file', ha='center', va='center', fontsize=14)
    plt.savefig("${plot_by_maf}", dpi=150, bbox_inches='tight')
    sys.exit(0)

# Check if data is empty
if data.empty:
    print("Warning: No data to plot", file=sys.stderr)
    # Create empty plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=14)
    ax.set_xlabel("${xlab}", fontsize=12)
    ax.set_ylabel("${ylab}", fontsize=12)
    plt.savefig("${plot_by_maf}", dpi=150, bbox_inches='tight')
    sys.exit(0)

# Determine the type of data based on columns
if 'Dataset' in data.columns and 'MAF_Range' in data.columns:
    # This is a summary table format from report_accuracy_by_maf
    # Group column is Dataset, x-axis is MAF_Range
    group_col = 'Dataset'
    maf_col = 'MAF_Range'
    
    # Determine value column to plot
    if 'Mean_Rsq' in data.columns:
        value_col = 'Mean_Rsq'
        # Convert to percentage
        data['plot_value'] = pd.to_numeric(data[value_col], errors='coerce') * 100
    elif 'Well_Imputed_Pct' in data.columns:
        value_col = 'Well_Imputed_Pct'
        data['plot_value'] = pd.to_numeric(data[value_col], errors='coerce')
    else:
        print("Error: No suitable value column found", file=sys.stderr)
        sys.exit(1)
    
    # Get unique groups and MAF ranges
    groups = data[group_col].unique()
    maf_ranges = data[maf_col].unique()
    
elif 'Dataset' in data.columns:
    # This is the format from report_well_imputed generate_plot_data
    # Data has Dataset column and MAF bin columns
    # Need to melt the data
    group_col = 'Dataset'
    
    # Get MAF bin columns (all columns except Dataset)
    maf_bin_cols = [col for col in data.columns if col != 'Dataset']
    
    # Melt the data
    data = pd.melt(data, id_vars=['Dataset'], var_name='MAF_Range', value_name='plot_value')
    # Convert values to numeric
    data['plot_value'] = pd.to_numeric(data['plot_value'], errors='coerce')
    
    groups = data['Dataset'].unique()
    maf_ranges = maf_bin_cols  # Use original column order
    maf_col = 'MAF_Range'
    
else:
    # Try old format with hardcoded columns
    print("Warning: Unexpected data format, trying legacy format", file=sys.stderr)
    # Assume first column is group identifier
    group_col = data.columns[0]
    # Rest are MAF bins
    data = pd.melt(data, id_vars=[group_col], var_name='MAF_Range', value_name='plot_value')
    data['plot_value'] = pd.to_numeric(data['plot_value'], errors='coerce')
    groups = data[group_col].unique()
    maf_ranges = data['MAF_Range'].unique()
    maf_col = 'MAF_Range'

# Create the plot
fig, ax = plt.subplots(figsize=(12, 6))

# Calculate bar positions
n_groups = len(groups)
n_maf_bins = len(maf_ranges)

if n_groups == 0 or n_maf_bins == 0:
    print("Error: No groups or MAF bins to plot", file=sys.stderr)
    ax.text(0.5, 0.5, 'No data to plot', ha='center', va='center', fontsize=14)
else:
    bar_width = 0.8 / n_groups if n_groups > 0 else 0.8
    x = range(n_maf_bins)
    
    # Plot bars for each group
    for i, group in enumerate(groups):
        group_data = data[data[group_col] == group]
        
        # Create values array for this group, maintaining MAF bin order
        values = []
        for maf_range in maf_ranges:
            maf_data = group_data[group_data[maf_col] == maf_range]
            if not maf_data.empty:
                val = maf_data['plot_value'].values[0]
                values.append(val if pd.notna(val) else 0)
            else:
                values.append(0)
        
        # Calculate x positions for this group's bars
        x_pos = [j + i * bar_width - (n_groups - 1) * bar_width / 2 for j in x]
        
        # Plot bars
        ax.bar(x_pos, values, bar_width, 
               label=group, 
               color=superpop_plot_colours[i % len(superpop_plot_colours)],
               edgecolor='black', linewidth=0.5)
    
    # Customize the plot
    ax.set_xlabel("${xlab}", fontsize=12)
    ax.set_ylabel("${ylab}", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(maf_ranges, rotation=45, ha='right')
    
    # Set y-axis limits based on data type
    if 'Mean_Rsq' in data.columns or 'accuracy' in "${ylab}".lower():
        ax.set_ylim(0, 105)  # For percentage data
        ax.axhline(y=100, color='gray', linestyle='--', alpha=0.5)
    
    # Add legend
    if n_groups > 1:
        ax.legend(title="${group}", bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.legend(title=group_col, bbox_to_anchor=(1.05, 1), loc='upper left')

# Add grid for better readability
ax.grid(True, alpha=0.3, axis='y')

# Add title
plt.title(f"Imputation Quality by MAF Bin", fontsize=14, fontweight='bold')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig("${plot_by_maf}", dpi=150, bbox_inches='tight')
plt.close()

print(f"Plot saved to ${plot_by_maf}")