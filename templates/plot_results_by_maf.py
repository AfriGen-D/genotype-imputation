#!/usr/bin/env python3
"""
Plot results by MAF bins
Converts R script to Python for plotting performance/accuracy by MAF
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

# Read the data
data = pd.read_csv("${report}", sep='\t')

# Rename columns to match MAF bins
data.columns = ["${group}", "(0,0.001]", "(0.001,0.01]", "(0.01,0.02]", 
                "(0.02,0.05]", "(0.05,0.2]", "(0.2,0.5]"]

# Melt the data for plotting
mdata = pd.melt(data, id_vars=["${group}"], var_name='MAF_bin', value_name='value')

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Get unique groups for color mapping
groups = mdata["${group}"].unique()
n_groups = len(groups)

# Create bar plot
x_labels = mdata['MAF_bin'].unique()
x = range(len(x_labels))
width = 0.8 / n_groups

for i, group in enumerate(groups):
    group_data = mdata[mdata["${group}"] == group]
    values = group_data.groupby('MAF_bin')['value'].first().reindex(x_labels, fill_value=0)
    x_pos = [j + width * i for j in x]
    ax.bar(x_pos, values, width, label=group, 
           color=superpop_plot_colours[i % len(superpop_plot_colours)])

# Customize the plot
ax.set_xlabel("${xlab}", fontsize=12)
ax.set_ylabel("${ylab}", fontsize=12)
ax.set_xticks([j + width * (n_groups - 1) / 2 for j in x])
ax.set_xticklabels(x_labels, rotation=45, ha='right')
ax.legend(title="${group}", bbox_to_anchor=(1.05, 1), loc='upper left')

# Add grid for better readability
ax.grid(True, alpha=0.3, axis='y')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig("${plot_by_maf}", dpi=150, bbox_inches='tight')
plt.close()