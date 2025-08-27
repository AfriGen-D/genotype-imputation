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
        
        required_fields = ['chromosome', 'maf_bins', 'mean_r2_by_maf', 'total_variants']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        return data
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        sys.exit(1)

def create_chr_maf_analysis_plot(chr_summary, output_file):
    """Create comprehensive chromosome-level MAF analysis visualization."""
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with 6 subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Extract data
    chromosome = chr_summary['chromosome']
    maf_bins = chr_summary['maf_bins']
    r2_by_maf = chr_summary['mean_r2_by_maf']
    total_variants = chr_summary['total_variants']
    
    # Prepare data for plotting
    maf_categories = list(maf_bins.keys())
    variant_counts = list(maf_bins.values())
    
    # Get R² statistics
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
    
    # Subplot 1: MAF distribution (bar plot)
    ax1 = plt.subplot(3, 2, 1)
    bars = ax1.bar(range(len(maf_categories)), variant_counts, 
                   color=sns.color_palette("viridis", len(maf_categories)))
    ax1.set_xlabel('MAF Categories')
    ax1.set_ylabel('Number of Variants')
    ax1.set_title(f'{chromosome}: Variant Distribution by MAF')
    ax1.set_xticks(range(len(maf_categories)))
    ax1.set_xticklabels(maf_categories, rotation=45)
    
    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, variant_counts)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=9)
    
    # Subplot 2: R² by MAF (line plot with error bars)
    ax2 = plt.subplot(3, 2, 2)
    x_pos = range(len(maf_categories))
    ax2.errorbar(x_pos, r2_means, yerr=r2_stds, 
                marker='o', linewidth=2, markersize=8, capsize=5)
    ax2.plot(x_pos, r2_medians, 'r--', marker='s', 
             linewidth=2, markersize=6, label='Median R²')
    ax2.set_xlabel('MAF Categories')
    ax2.set_ylabel('R² Value')
    ax2.set_title(f'{chromosome}: Imputation Quality by MAF')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(maf_categories, rotation=45)
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: MAF distribution (log scale)
    ax3 = plt.subplot(3, 2, 3)
    ax3.bar(range(len(maf_categories)), variant_counts, 
            color=sns.color_palette("plasma", len(maf_categories)))
    ax3.set_xlabel('MAF Categories')
    ax3.set_ylabel('Number of Variants (log scale)')
    ax3.set_title(f'{chromosome}: Variant Distribution (Log Scale)')
    ax3.set_xticks(range(len(maf_categories)))
    ax3.set_xticklabels(maf_categories, rotation=45)
    ax3.set_yscale('log')
    
    # Subplot 4: R² distribution violin plot
    ax4 = plt.subplot(3, 2, 4)
    
    # Create violin plot data (simulated from means and stds)
    violin_data = []
    labels = []
    for i, (cat, mean, std) in enumerate(zip(maf_categories, r2_means, r2_stds)):
        if std > 0:  # Only if we have valid data
            # Simulate data points for violin plot
            data_points = np.random.normal(mean, std, 100)
            data_points = np.clip(data_points, 0, 1)  # Clip to valid R² range
            violin_data.append(data_points)
            labels.append(cat)
    
    if violin_data:
        parts = ax4.violinplot(violin_data, positions=range(len(violin_data)), 
                              showmeans=True, showmedians=True)
        ax4.set_xlabel('MAF Categories')
        ax4.set_ylabel('R² Distribution')
        ax4.set_title(f'{chromosome}: R² Distribution by MAF')
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45)
        ax4.set_ylim(0, 1)
    else:
        ax4.text(0.5, 0.5, 'No R² distribution data available', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title(f'{chromosome}: R² Distribution by MAF')
    
    # Subplot 5: Chunk-level variability (if available)
    ax5 = plt.subplot(3, 2, 5)
    
    if 'chunk_details' in chr_summary and chr_summary['chunk_details']:
        chunk_r2_values = [chunk['mean_r2'] for chunk in chr_summary['chunk_details']]
        chunk_ids = [chunk['chunk_id'].split('_')[-2:] for chunk in chr_summary['chunk_details']]  # Get position info
        
        ax5.scatter(range(len(chunk_r2_values)), chunk_r2_values, 
                   alpha=0.7, s=60, color='darkblue')
        ax5.axhline(y=np.mean(chunk_r2_values), color='red', linestyle='--', 
                   label=f'Mean R² = {np.mean(chunk_r2_values):.3f}')
        ax5.set_xlabel('Chunk Index')
        ax5.set_ylabel('Mean R²')
        ax5.set_title(f'{chromosome}: R² Variability Across Chunks')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
    else:
        ax5.text(0.5, 0.5, 'No chunk-level data available', 
                ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title(f'{chromosome}: R² Variability Across Chunks')
    
    # Subplot 6: Summary statistics table
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    # Create summary table
    summary_data = []
    summary_data.append(['Total Variants', f"{total_variants:,}"])
    mean_r2 = chr_summary.get('mean_r2', 0)
    if mean_r2 is None:
        mean_r2 = 0
    summary_data.append(['Mean R²', f"{mean_r2:.4f}"])
    
    if 'info_score_stats' in chr_summary:
        stats = chr_summary['info_score_stats']
        summary_data.append(['Well-imputed (R²≥0.3)', f"{stats.get('above_0.3', 0):,}"])
        summary_data.append(['High-quality (R²≥0.8)', f"{stats.get('above_0.8', 0):,}"])
        summary_data.append(['Median INFO', f"{stats.get('median', 0):.4f}"])
    
    # Add MAF category with most variants
    max_maf_cat = max(maf_bins.items(), key=lambda x: x[1])
    summary_data.append(['Largest MAF category', f"{max_maf_cat[0]} ({max_maf_cat[1]:,} variants)"])
    
    # Best performing MAF category
    if r2_means:
        best_maf_idx = np.argmax(r2_means)
        best_maf_cat = maf_categories[best_maf_idx]
        summary_data.append(['Best R² category', f"{best_maf_cat} (R²={r2_means[best_maf_idx]:.3f})"])
    
    table = ax6.table(cellText=summary_data,
                     colLabels=['Metric', 'Value'],
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2)
    
    # Style the table
    for i in range(len(summary_data) + 1):
        for j in range(2):
            cell = table[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#40466e')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')
    
    ax6.set_title(f'{chromosome}: Summary Statistics', pad=20, fontweight='bold')
    
    # Overall title and layout
    fig.suptitle(f'Chromosome {chromosome.replace("chr", "")}: Minor Allele Frequency Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.3, wspace=0.3)
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created chromosome MAF analysis plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level MAF analysis plots')
    parser.add_argument('--chr-summary', required=True, 
                       help='Path to chromosome summary JSON file')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--chromosome', required=True, help='Chromosome name')
    
    args = parser.parse_args()
    
    # Load data
    chr_summary = load_chromosome_summary(args.chr_summary)
    
    # Create plot
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_maf_analysis.pdf"
    create_chr_maf_analysis_plot(chr_summary, output_file)

if __name__ == '__main__':
    main()