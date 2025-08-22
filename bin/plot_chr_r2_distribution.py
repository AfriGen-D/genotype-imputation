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
        
        required_fields = ['chromosome', 'info_score_stats', 'mean_r2_by_maf']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        return data
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        sys.exit(1)

def create_chr_r2_distribution_plot(chr_summary, output_file):
    """Create comprehensive chromosome-level R² distribution visualization."""
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with 6 subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Extract data
    chromosome = chr_summary['chromosome']
    info_stats = chr_summary['info_score_stats']
    r2_by_maf = chr_summary['mean_r2_by_maf']
    
    # Subplot 1: R² distribution histogram
    ax1 = plt.subplot(3, 2, 1)
    
    # Create simulated R² distribution from statistics
    mean_r2 = info_stats.get('mean', 0.5)
    # Use interquartile range to estimate distribution shape
    q25 = info_stats.get('q25', 0.1)
    q75 = info_stats.get('q75', 0.8)
    
    # Create bins for histogram
    r2_bins = np.linspace(0, 1, 21)  # 20 bins from 0 to 1
    
    # Simulate data points based on quartiles
    np.random.seed(42)  # For reproducibility
    n_points = 10000
    
    # Create a mixture of distributions to match quartiles
    low_r2 = np.random.beta(0.5, 2, n_points//3)  # Skewed towards 0
    mid_r2 = np.random.normal(0.5, 0.2, n_points//3)  # Normal around 0.5
    high_r2 = np.random.beta(2, 0.5, n_points//3)  # Skewed towards 1
    
    # Combine and clip
    r2_data = np.concatenate([low_r2, mid_r2, high_r2])
    r2_data = np.clip(r2_data, 0, 1)
    
    # Scale to match actual statistics
    r2_data = np.sort(r2_data)
    target_q25_idx = int(len(r2_data) * 0.25)
    target_q75_idx = int(len(r2_data) * 0.75)
    
    # Adjust to match quartiles
    scale_factor = (q75 - q25) / (r2_data[target_q75_idx] - r2_data[target_q25_idx])
    r2_data = q25 + (r2_data - r2_data[target_q25_idx]) * scale_factor
    r2_data = np.clip(r2_data, 0, 1)
    
    hist, bins, _ = ax1.hist(r2_data, bins=r2_bins, alpha=0.7, 
                            color='steelblue', edgecolor='black', linewidth=0.5)
    ax1.axvline(mean_r2, color='red', linestyle='--', linewidth=2, label=f'Mean = {mean_r2:.3f}')
    ax1.axvline(q25, color='orange', linestyle=':', linewidth=2, label=f'Q25 = {q25:.3f}')
    ax1.axvline(q75, color='orange', linestyle=':', linewidth=2, label=f'Q75 = {q75:.3f}')
    ax1.set_xlabel('R² Value')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'{chromosome}: R² Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Subplot 2: Cumulative distribution
    ax2 = plt.subplot(3, 2, 2)
    sorted_r2 = np.sort(r2_data)
    cumulative = np.arange(1, len(sorted_r2) + 1) / len(sorted_r2)
    ax2.plot(sorted_r2, cumulative, linewidth=2, color='darkgreen')
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.7, label='Median')
    ax2.axvline(info_stats.get('median', 0.5), color='red', linestyle='--', alpha=0.7)
    ax2.set_xlabel('R² Value')
    ax2.set_ylabel('Cumulative Probability')
    ax2.set_title(f'{chromosome}: Cumulative R² Distribution')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Subplot 3: Quality thresholds bar chart
    ax3 = plt.subplot(3, 2, 3)
    
    total_variants = chr_summary.get('total_variants', 1)
    above_03 = info_stats.get('above_0.3', 0)
    above_08 = info_stats.get('above_0.8', 0)
    
    below_03 = total_variants - above_03
    between_03_08 = above_03 - above_08
    
    categories = ['Poor\n(R² < 0.3)', 'Well-imputed\n(0.3 ≤ R² < 0.8)', 'High-quality\n(R² ≥ 0.8)']
    counts = [below_03, between_03_08, above_08]
    colors = ['lightcoral', 'gold', 'lightgreen']
    
    bars = ax3.bar(categories, counts, color=colors, edgecolor='black', linewidth=1)
    ax3.set_ylabel('Number of Variants')
    ax3.set_title(f'{chromosome}: Quality Distribution')
    
    # Add percentage labels
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        percentage = (count / total_variants) * 100
        ax3.text(bar.get_x() + bar.get_width()/2., height + total_variants*0.01,
                f'{count:,}\n({percentage:.1f}%)', 
                ha='center', va='bottom', fontsize=10)
    
    # Subplot 4: R² by MAF box plot
    ax4 = plt.subplot(3, 2, 4)
    
    maf_categories = list(r2_by_maf.keys())
    r2_means = [r2_by_maf[cat]['mean'] for cat in maf_categories]
    r2_medians = [r2_by_maf[cat]['median'] for cat in maf_categories]
    r2_stds = [r2_by_maf[cat]['std'] for cat in maf_categories]
    
    # Create box plot data (simulated from statistics)
    box_data = []
    for mean, std in zip(r2_means, r2_stds):
        if std > 0:
            data = np.random.normal(mean, std, 100)
            data = np.clip(data, 0, 1)
            box_data.append(data)
        else:
            box_data.append([mean] * 10)  # Constant values if no std
    
    bp = ax4.boxplot(box_data, labels=maf_categories, patch_artist=True)
    
    # Color the boxes
    colors = sns.color_palette("viridis", len(maf_categories))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax4.set_xlabel('MAF Categories')
    ax4.set_ylabel('R² Value')
    ax4.set_title(f'{chromosome}: R² Distribution by MAF')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(True, alpha=0.3)
    
    # Subplot 5: Chunk-level R² variability
    ax5 = plt.subplot(3, 2, 5)
    
    if 'chunk_details' in chr_summary and chr_summary['chunk_details']:
        chunk_r2_values = [chunk['mean_r2'] for chunk in chr_summary['chunk_details']]
        
        # Create histogram of chunk R² values
        ax5.hist(chunk_r2_values, bins=15, alpha=0.7, color='purple', 
                edgecolor='black', linewidth=0.5)
        ax5.axvline(np.mean(chunk_r2_values), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean = {np.mean(chunk_r2_values):.3f}')
        ax5.axvline(np.median(chunk_r2_values), color='orange', linestyle=':', 
                   linewidth=2, label=f'Median = {np.median(chunk_r2_values):.3f}')
        ax5.set_xlabel('Chunk Mean R²')
        ax5.set_ylabel('Number of Chunks')
        ax5.set_title(f'{chromosome}: Chunk R² Variability')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
    else:
        ax5.text(0.5, 0.5, 'No chunk-level data available', 
                ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title(f'{chromosome}: Chunk R² Variability')
    
    # Subplot 6: Summary statistics table
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    # Create summary table
    summary_data = []
    summary_data.append(['Total Variants', f"{total_variants:,}"])
    summary_data.append(['Mean R²', f"{info_stats.get('mean', 0):.4f}"])
    summary_data.append(['Median R²', f"{info_stats.get('median', 0):.4f}"])
    summary_data.append(['Min R²', f"{info_stats.get('min', 0):.4f}"])
    summary_data.append(['Max R²', f"{info_stats.get('max', 0):.4f}"])
    summary_data.append(['Well-imputed (≥0.3)', f"{above_03:,} ({(above_03/total_variants)*100:.1f}%)"])
    summary_data.append(['High-quality (≥0.8)', f"{above_08:,} ({(above_08/total_variants)*100:.1f}%)"])
    
    # Best MAF category
    if r2_means:
        best_maf_idx = np.argmax(r2_means)
        best_maf_cat = maf_categories[best_maf_idx]
        summary_data.append(['Best MAF category', f"{best_maf_cat} (R²={r2_means[best_maf_idx]:.3f})"])
    
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
    
    ax6.set_title(f'{chromosome}: R² Statistics Summary', pad=20, fontweight='bold')
    
    # Overall title and layout
    fig.suptitle(f'Chromosome {chromosome.replace("chr", "")}: R² Distribution Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.3, wspace=0.3)
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created chromosome R² distribution plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level R² distribution plots')
    parser.add_argument('--chr-summary', required=True, 
                       help='Path to chromosome summary JSON file')
    parser.add_argument('--output', required=True,
                       help='Output PDF file path')
    
    args = parser.parse_args()
    
    # Load data
    chr_summary = load_chromosome_summary(args.chr_summary)
    
    # Create plot
    create_chr_r2_distribution_plot(chr_summary, args.output)

if __name__ == '__main__':
    main()