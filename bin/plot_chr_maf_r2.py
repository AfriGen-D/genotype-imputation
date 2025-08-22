#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import argparse
import sys
import os

# Set matplotlib backend for headless environments
import matplotlib
matplotlib.use('Agg')

def load_chromosome_summary(json_file):
    """Load and validate chromosome summary JSON."""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        required_fields = ['chromosome', 'maf_bins', 'mean_r2_by_maf']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        return data
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        sys.exit(1)

def maf_category_to_numeric(maf_cat):
    """Convert MAF category string to numeric midpoint."""
    try:
        if '-' in maf_cat:
            start, end = maf_cat.split('-')
            return (float(start) + float(end)) / 2
        else:
            return float(maf_cat)
    except:
        return 0.25  # Default fallback

def create_chr_maf_r2_plot(chr_summary, output_file):
    """Create comprehensive chromosome-level MAF vs R-squared visualization."""
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with 6 subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Extract data
    chromosome = chr_summary['chromosome']
    maf_bins = chr_summary['maf_bins']
    r2_by_maf = chr_summary['mean_r2_by_maf']
    total_variants = chr_summary.get('total_variants', 0)
    
    # Prepare data for plotting
    maf_categories = list(maf_bins.keys())
    variant_counts = list(maf_bins.values())
    
    # Convert MAF categories to numeric values for scatter plots
    maf_numeric = [maf_category_to_numeric(cat) for cat in maf_categories]
    
    # Get R� statistics
    r2_means = []
    r2_medians = []
    r2_stds = []
    
    for cat in maf_categories:
        if cat in r2_by_maf:
            r2_means.append(r2_by_maf[cat]['mean'])
            r2_medians.append(r2_by_maf[cat]['median'])
            r2_stds.append(r2_by_maf[cat]['std'])
        else:
            r2_means.append(0)
            r2_medians.append(0)
            r2_stds.append(0)
    
    # Subplot 1: MAF vs R� scatter plot with error bars
    ax1 = plt.subplot(3, 2, 1)
    
    # Create scatter plot with point sizes proportional to variant count
    max_count = max(variant_counts) if variant_counts else 1
    point_sizes = [100 + (count/max_count) * 400 for count in variant_counts]
    
    scatter = ax1.scatter(maf_numeric, r2_means, s=point_sizes, c=r2_means, 
                         cmap='RdYlGn', alpha=0.7, edgecolors='black', linewidth=1)
    
    # Add error bars
    ax1.errorbar(maf_numeric, r2_means, yerr=r2_stds, fmt='none', 
                ecolor='black', alpha=0.5, capsize=3)
    
    # Add trend line
    if len(maf_numeric) > 1:
        z = np.polyfit(maf_numeric, r2_means, 1)
        p = np.poly1d(z)
        ax1.plot(maf_numeric, p(maf_numeric), "r--", alpha=0.8, linewidth=2, 
                label=f'Trend: R� = {z[0]:.2f}�MAF + {z[1]:.3f}')
    
    ax1.set_xlabel('Minor Allele Frequency')
    ax1.set_ylabel('Mean R�')
    ax1.set_title(f'{chromosome}: MAF vs R� Relationship')
    ax1.set_xlim(0, 0.5)
    ax1.set_ylim(0, 1)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('R� Value')
    
    # Subplot 2: MAF vs R� with confidence intervals
    ax2 = plt.subplot(3, 2, 2)
    
    # Plot means with confidence intervals
    ax2.plot(maf_numeric, r2_means, 'o-', linewidth=2, markersize=8, 
            color='darkblue', label='Mean R�')
    ax2.fill_between(maf_numeric, 
                    [max(0, m-s) for m, s in zip(r2_means, r2_stds)],
                    [min(1, m+s) for m, s in zip(r2_means, r2_stds)],
                    alpha=0.3, color='blue', label='�1 SD')
    
    # Add median line
    ax2.plot(maf_numeric, r2_medians, 's--', linewidth=2, markersize=6, 
            color='red', label='Median R�')
    
    # Add quality thresholds
    ax2.axhline(y=0.3, color='orange', linestyle=':', alpha=0.7, label='Well-imputed (0.3)')
    ax2.axhline(y=0.8, color='green', linestyle=':', alpha=0.7, label='High-quality (0.8)')
    
    ax2.set_xlabel('Minor Allele Frequency')
    ax2.set_ylabel('R� Value')
    ax2.set_title(f'{chromosome}: R� Confidence by MAF')
    ax2.set_xlim(0, 0.5)
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Variant density vs MAF (log scale)
    ax3 = plt.subplot(3, 2, 3)
    
    bars = ax3.bar(maf_numeric, variant_counts, width=0.02, 
                  color=sns.color_palette("viridis", len(maf_categories)), 
                  alpha=0.7, edgecolor='black', linewidth=1)
    
    ax3.set_xlabel('Minor Allele Frequency')
    ax3.set_ylabel('Number of Variants (log scale)')
    ax3.set_title(f'{chromosome}: Variant Density by MAF')
    ax3.set_xlim(0, 0.5)
    ax3.set_yscale('log')
    ax3.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for i, (bar, count, maf) in enumerate(zip(bars, variant_counts, maf_numeric)):
        height = bar.get_height()
        ax3.text(maf, height * 1.5, f'{count:,}', 
                ha='center', va='bottom', fontsize=9, rotation=45)
    
    # Subplot 4: R� efficiency by MAF (R� per variant)
    ax4 = plt.subplot(3, 2, 4)
    
    # Calculate R� efficiency (mean R� weighted by variant count)
    r2_efficiency = []
    for r2_mean, count in zip(r2_means, variant_counts):
        if count > 0:
            r2_efficiency.append(r2_mean * np.log10(count + 1))  # Log-weighted
        else:
            r2_efficiency.append(0)
    
    bars4 = ax4.bar(maf_numeric, r2_efficiency, width=0.02,
                   color=sns.color_palette("plasma", len(maf_categories)), 
                   alpha=0.7, edgecolor='black', linewidth=1)
    
    ax4.set_xlabel('Minor Allele Frequency')
    ax4.set_ylabel('R� Efficiency (R� � log��(variants))')
    ax4.set_title(f'{chromosome}: Imputation Efficiency by MAF')
    ax4.set_xlim(0, 0.5)
    ax4.grid(True, alpha=0.3)
    
    # Subplot 5: MAF category performance comparison
    ax5 = plt.subplot(3, 2, 5)
    
    # Create a comparison chart showing different metrics
    x_pos = range(len(maf_categories))
    width = 0.25
    
    # Normalize values for comparison
    norm_r2 = [r / max(r2_means) if max(r2_means) > 0 else 0 for r in r2_means]
    norm_counts = [c / max(variant_counts) if max(variant_counts) > 0 else 0 for c in variant_counts]
    norm_stds = [s / max(r2_stds) if max(r2_stds) > 0 else 0 for s in r2_stds]
    
    ax5.bar([x - width for x in x_pos], norm_r2, width, 
           label='Normalized R�', alpha=0.8, color='lightblue')
    ax5.bar(x_pos, norm_counts, width, 
           label='Normalized Variant Count', alpha=0.8, color='lightgreen')
    ax5.bar([x + width for x in x_pos], [1-s for s in norm_stds], width, 
           label='Stability (1-norm_std)', alpha=0.8, color='lightcoral')
    
    ax5.set_xlabel('MAF Categories')
    ax5.set_ylabel('Normalized Values')
    ax5.set_title(f'{chromosome}: MAF Category Performance')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels(maf_categories, rotation=45)
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Subplot 6: Summary statistics table
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    # Create summary table
    summary_data = []
    summary_data.append(['Total Variants', f"{total_variants:,}"])
    summary_data.append(['MAF Categories', f"{len(maf_categories)}"])
    
    # Overall R� statistics
    weighted_r2 = sum(r * c for r, c in zip(r2_means, variant_counts)) / sum(variant_counts) if sum(variant_counts) > 0 else 0
    summary_data.append(['Weighted Mean R�', f"{weighted_r2:.4f}"])
    
    # Best and worst MAF categories
    if r2_means:
        best_idx = np.argmax(r2_means)
        worst_idx = np.argmin(r2_means)
        summary_data.append(['Best MAF Category', f"{maf_categories[best_idx]} (R�={r2_means[best_idx]:.3f})"])
        summary_data.append(['Worst MAF Category', f"{maf_categories[worst_idx]} (R�={r2_means[worst_idx]:.3f})"])
    
    # Most and least common MAF categories
    if variant_counts:
        most_common_idx = np.argmax(variant_counts)
        least_common_idx = np.argmin(variant_counts)
        summary_data.append(['Most Common MAF', f"{maf_categories[most_common_idx]} ({variant_counts[most_common_idx]:,} variants)"])
        summary_data.append(['Least Common MAF', f"{maf_categories[least_common_idx]} ({variant_counts[least_common_idx]:,} variants)"])
    
    # Correlation between MAF and R�
    if len(maf_numeric) > 1:
        correlation = np.corrcoef(maf_numeric, r2_means)[0, 1]
        summary_data.append(['MAF-R� Correlation', f"{correlation:.3f}"])
    
    # Quality metrics
    high_quality_vars = sum(c for r, c in zip(r2_means, variant_counts) if r >= 0.8)
    well_imputed_vars = sum(c for r, c in zip(r2_means, variant_counts) if r >= 0.3)
    summary_data.append(['High Quality (e0.8)', f"{high_quality_vars:,} ({high_quality_vars/total_variants*100:.1f}%)"])
    summary_data.append(['Well Imputed (e0.3)', f"{well_imputed_vars:,} ({well_imputed_vars/total_variants*100:.1f}%)"])
    
    table = ax6.table(cellText=summary_data,
                     colLabels=['Metric', 'Value'],
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.8)
    
    # Style the table
    for i in range(len(summary_data) + 1):
        for j in range(2):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#40466e')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')
    
    ax6.set_title(f'{chromosome}: MAF-R� Analysis Summary', pad=20, fontweight='bold')
    
    # Overall title and layout
    fig.suptitle(f'Chromosome {chromosome.replace("chr", "")}: Minor Allele Frequency vs R� Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.3, wspace=0.3)
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created chromosome MAF vs R� plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level MAF vs R� plots')
    parser.add_argument('--chr-summary', required=True, 
                       help='Path to chromosome summary JSON file')
    parser.add_argument('--output', required=True,
                       help='Output PDF file path')
    
    args = parser.parse_args()
    
    # Load data
    chr_summary = load_chromosome_summary(args.chr_summary)
    
    # Create plot
    create_chr_maf_r2_plot(chr_summary, args.output)

if __name__ == '__main__':
    main()