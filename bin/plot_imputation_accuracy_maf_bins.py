#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot imputation accuracy (R²) across MAF bins - similar to Terra All of Us figures
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import json

def parse_info_file(info_file):
    """Parse Minimac4 info file to extract R² and MAF values"""
    data = []
    with open(info_file, 'r') as f:
        header = f.readline().strip().split('\t')
        maf_idx = header.index('MAF') if 'MAF' in header else -1
        r2_idx = header.index('Rsq') if 'Rsq' in header else header.index('R2') if 'R2' in header else -1
        
        if maf_idx == -1 or r2_idx == -1:
            print(f"Warning: Could not find MAF or R² columns in {info_file}")
            return pd.DataFrame()
        
        for line in f:
            parts = line.strip().split('\t')
            try:
                maf = float(parts[maf_idx])
                r2 = float(parts[r2_idx])
                data.append({'MAF': maf, 'R2': r2})
            except (ValueError, IndexError):
                continue
    
    return pd.DataFrame(data)

def create_maf_bins(df, bins=None):
    """Create MAF bins for analysis"""
    if bins is None:
        bins = [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    
    df['MAF_bin'] = pd.cut(df['MAF'], bins=bins, include_lowest=True)
    return df

def plot_accuracy_by_maf(df, output_prefix, title="Imputation Accuracy by MAF"):
    """Create Terra-style imputation accuracy plot"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Define MAF bins
    maf_bins = [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    df = create_maf_bins(df, maf_bins)
    
    # Color palette similar to Terra
    colors = sns.color_palette("viridis", n_colors=len(maf_bins)-1)
    
    # 1. Box plot of R² by MAF bins
    ax = axes[0, 0]
    maf_order = df['MAF_bin'].cat.categories
    box_data = [df[df['MAF_bin'] == bin]['R2'].dropna() for bin in maf_order]
    bp = ax.boxplot(box_data, patch_artist=True, notch=True)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Minor Allele Frequency Bins', fontsize=12)
    ax.set_ylabel('Imputation R²', fontsize=12)
    ax.set_title('Imputation Accuracy by MAF Bins', fontsize=14, fontweight='bold')
    ax.set_xticklabels([f'{b:.3f}' if b < 0.1 else f'{b:.1f}' 
                        for b in maf_bins[1:]], rotation=45)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1])
    
    # Add median line
    medians = [np.median(d) if len(d) > 0 else 0 for d in box_data]
    ax.plot(range(1, len(medians)+1), medians, 'r-', linewidth=2, label='Median R²')
    ax.legend()
    
    # 2. Scatter plot with trend line
    ax = axes[0, 1]
    sample_size = min(10000, len(df))
    df_sample = df.sample(n=sample_size) if len(df) > sample_size else df
    
    scatter = ax.scatter(df_sample['MAF'], df_sample['R2'], 
                        c=df_sample['MAF'], cmap='viridis', 
                        alpha=0.3, s=1)
    
    # Add binned means
    bin_means = df.groupby('MAF_bin').agg({'MAF': 'mean', 'R2': 'mean'}).reset_index()
    ax.plot(bin_means['MAF'], bin_means['R2'], 'r-', linewidth=2, label='Mean R²')
    ax.scatter(bin_means['MAF'], bin_means['R2'], color='red', s=50, zorder=5)
    
    ax.set_xlabel('Minor Allele Frequency', fontsize=12)
    ax.set_ylabel('Imputation R²', fontsize=12)
    ax.set_title('R² vs MAF Scatter Plot', fontsize=14, fontweight='bold')
    ax.set_xlim([0, 0.5])
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.colorbar(scatter, ax=ax, label='MAF')
    
    # 3. Histogram of R² distribution by MAF categories
    ax = axes[1, 0]
    
    # Define MAF categories
    maf_categories = {
        'Rare (MAF<0.01)': df[df['MAF'] < 0.01],
        'Low (0.01≤MAF<0.05)': df[(df['MAF'] >= 0.01) & (df['MAF'] < 0.05)],
        'Common (MAF≥0.05)': df[df['MAF'] >= 0.05]
    }
    
    for i, (label, data) in enumerate(maf_categories.items()):
        if len(data) > 0:
            ax.hist(data['R2'], bins=50, alpha=0.5, label=label, 
                   color=['red', 'orange', 'green'][i])
    
    ax.set_xlabel('Imputation R²', fontsize=12)
    ax.set_ylabel('Number of Variants', fontsize=12)
    ax.set_title('R² Distribution by Variant Frequency', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Cumulative accuracy plot
    ax = axes[1, 1]
    
    for label, data in maf_categories.items():
        if len(data) > 0:
            sorted_r2 = np.sort(data['R2'])
            cumulative = np.arange(1, len(sorted_r2) + 1) / len(sorted_r2)
            ax.plot(sorted_r2, cumulative, linewidth=2, label=label)
    
    # Add threshold lines
    for threshold in [0.3, 0.5, 0.8]:
        ax.axvline(threshold, color='gray', linestyle='--', alpha=0.5)
        ax.text(threshold, 0.02, f'R²={threshold}', rotation=90, 
               verticalalignment='bottom', fontsize=8)
    
    ax.set_xlabel('Imputation R²', fontsize=12)
    ax.set_ylabel('Cumulative Proportion of Variants', fontsize=12)
    ax.set_title('Cumulative R² Distribution', fontsize=14, fontweight='bold')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'{title}\\nn = {len(df):,} variants', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    # Save figure
    plt.savefig(f'{output_prefix}_maf_accuracy.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_maf_accuracy.png', dpi=150, bbox_inches='tight')
    
    # Generate summary statistics
    generate_summary_stats(df, output_prefix)

def generate_summary_stats(df, output_prefix):
    """Generate summary statistics for the imputation"""
    stats = {
        'total_variants': len(df),
        'mean_r2': df['R2'].mean(),
        'median_r2': df['R2'].median(),
        'well_imputed_0.3': (df['R2'] >= 0.3).sum(),
        'well_imputed_0.5': (df['R2'] >= 0.5).sum(),
        'well_imputed_0.8': (df['R2'] >= 0.8).sum(),
        'prop_well_imputed_0.3': (df['R2'] >= 0.3).mean(),
        'prop_well_imputed_0.5': (df['R2'] >= 0.5).mean(),
        'prop_well_imputed_0.8': (df['R2'] >= 0.8).mean()
    }
    
    # MAF-specific statistics
    maf_categories = [
        ('rare', df[df['MAF'] < 0.01]),
        ('low_freq', df[(df['MAF'] >= 0.01) & (df['MAF'] < 0.05)]),
        ('common', df[df['MAF'] >= 0.05])
    ]
    
    for cat_name, cat_df in maf_categories:
        if len(cat_df) > 0:
            stats[f'{cat_name}_count'] = len(cat_df)
            stats[f'{cat_name}_mean_r2'] = cat_df['R2'].mean()
            stats[f'{cat_name}_median_r2'] = cat_df['R2'].median()
            stats[f'{cat_name}_prop_well_imputed_0.3'] = (cat_df['R2'] >= 0.3).mean()
            stats[f'{cat_name}_prop_well_imputed_0.8'] = (cat_df['R2'] >= 0.8).mean()
    
    # Save statistics
    with open(f'{output_prefix}_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Create text summary
    with open(f'{output_prefix}_summary.txt', 'w') as f:
        f.write("=" * 80 + "\\n")
        f.write("IMPUTATION ACCURACY SUMMARY\\n")
        f.write("=" * 80 + "\\n\\n")
        
        f.write(f"Total variants analyzed: {stats['total_variants']:,}\\n")
        f.write(f"Mean R²: {stats['mean_r2']:.4f}\\n")
        f.write(f"Median R²: {stats['median_r2']:.4f}\\n\\n")
        
        f.write("Well-imputed variants (R² ≥ threshold):\\n")
        f.write(f"  R² ≥ 0.3: {stats['well_imputed_0.3']:,} ({stats['prop_well_imputed_0.3']:.1%})\\n")
        f.write(f"  R² ≥ 0.5: {stats['well_imputed_0.5']:,} ({stats['prop_well_imputed_0.5']:.1%})\\n")
        f.write(f"  R² ≥ 0.8: {stats['well_imputed_0.8']:,} ({stats['prop_well_imputed_0.8']:.1%})\\n\\n")
        
        f.write("By allele frequency:\\n")
        for cat_name, cat_label in [('rare', 'Rare (MAF<1%)'), 
                                    ('low_freq', 'Low frequency (1%≤MAF<5%)'),
                                    ('common', 'Common (MAF≥5%)')]:
            if f'{cat_name}_count' in stats:
                f.write(f"\\n  {cat_label}:\\n")
                f.write(f"    Count: {stats[f'{cat_name}_count']:,}\\n")
                f.write(f"    Mean R²: {stats[f'{cat_name}_mean_r2']:.4f}\\n")
                f.write(f"    R² ≥ 0.3: {stats[f'{cat_name}_prop_well_imputed_0.3']:.1%}\\n")
                f.write(f"    R² ≥ 0.8: {stats[f'{cat_name}_prop_well_imputed_0.8']:.1%}\\n")

def main():
    parser = argparse.ArgumentParser(description='Plot imputation accuracy by MAF bins')
    parser.add_argument('info_file', help='Minimac4 info file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--title', default='Imputation Accuracy Analysis', 
                       help='Plot title')
    parser.add_argument('--sample-name', help='Sample name for title')
    
    args = parser.parse_args()
    
    # Parse info file
    df = parse_info_file(args.info_file)
    
    if df.empty:
        print(f"Error: No data found in {args.info_file}")
        sys.exit(1)
    
    # Create title
    title = args.title
    if args.sample_name:
        title = f"{args.sample_name} - {title}"
    
    # Generate plots
    plot_accuracy_by_maf(df, args.output_prefix, title)
    
    print(f"Plots saved to {args.output_prefix}_maf_accuracy.pdf")
    print(f"Statistics saved to {args.output_prefix}_stats.json")

if __name__ == '__main__':
    main()