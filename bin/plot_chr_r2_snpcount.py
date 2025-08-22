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
        
        required_fields = ['chromosome', 'chunk_details', 'total_variants']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        return data
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        sys.exit(1)

def create_chr_r2_snpcount_plot(chr_summary, output_file):
    """Create comprehensive chromosome-level R-squared vs SNP count visualization."""
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with 6 subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Extract data
    chromosome = chr_summary['chromosome']
    chunk_details = chr_summary['chunk_details']
    total_variants = chr_summary['total_variants']
    
    if not chunk_details:
        # Create empty plot with message
        ax = plt.subplot(1, 1, 1)
        ax.text(0.5, 0.5, f'No chunk data available for {chromosome}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(f'{chromosome}: R� vs SNP Count Analysis')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Extract chunk-level data
    variant_counts = [chunk['variants'] for chunk in chunk_details]
    r2_values = [chunk['mean_r2'] for chunk in chunk_details]
    well_imputed_counts = [chunk['well_imputed'] for chunk in chunk_details]
    
    # Convert to arrays
    variant_counts = np.array(variant_counts)
    r2_values = np.array(r2_values)
    well_imputed_counts = np.array(well_imputed_counts)
    
    # Calculate additional metrics
    well_imputed_fractions = well_imputed_counts / variant_counts
    chunk_ids = [chunk['chunk_id'] for chunk in chunk_details]
    
    # Subplot 1: R� vs SNP count scatter plot
    ax1 = plt.subplot(3, 2, 1)
    
    # Color points by well-imputed fraction
    scatter = ax1.scatter(variant_counts, r2_values, c=well_imputed_fractions, 
                         cmap='viridis', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Add trend line
    if len(variant_counts) > 1:
        z = np.polyfit(variant_counts, r2_values, 1)
        p = np.poly1d(z)
        ax1.plot(variant_counts, p(variant_counts), "r--", alpha=0.8, linewidth=2, 
                label=f'Trend: R� = {z[0]:.2e}�SNP + {z[1]:.3f}')
        
        # Calculate correlation
        correlation = np.corrcoef(variant_counts, r2_values)[0, 1]
        ax1.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=ax1.transAxes, bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.8))
    
    ax1.set_xlabel('Number of Variants per Chunk')
    ax1.set_ylabel('Mean R�')
    ax1.set_title(f'{chromosome}: R� vs Variant Count')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Well-imputed Fraction')
    
    # Subplot 2: SNP count distribution
    ax2 = plt.subplot(3, 2, 2)
    
    # Histogram of variant counts
    hist, bins, _ = ax2.hist(variant_counts, bins=15, alpha=0.7, color='steelblue', 
                            edgecolor='black', linewidth=0.5)
    
    # Add statistics lines
    mean_count = np.mean(variant_counts)
    median_count = np.median(variant_counts)
    ax2.axvline(mean_count, color='red', linestyle='--', linewidth=2, 
               label=f'Mean = {mean_count:,.0f}')
    ax2.axvline(median_count, color='orange', linestyle=':', linewidth=2, 
               label=f'Median = {median_count:,.0f}')
    
    ax2.set_xlabel('Variants per Chunk')
    ax2.set_ylabel('Number of Chunks')
    ax2.set_title(f'{chromosome}: Variant Count Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: R� distribution by variant count quartiles
    ax3 = plt.subplot(3, 2, 3)
    
    # Create quartile bins
    quartiles = np.percentile(variant_counts, [25, 50, 75])
    
    # Assign chunks to quartiles
    q1_r2 = r2_values[variant_counts <= quartiles[0]]
    q2_r2 = r2_values[(variant_counts > quartiles[0]) & (variant_counts <= quartiles[1])]
    q3_r2 = r2_values[(variant_counts > quartiles[1]) & (variant_counts <= quartiles[2])]
    q4_r2 = r2_values[variant_counts > quartiles[2]]
    
    # Box plot
    box_data = [q1_r2, q2_r2, q3_r2, q4_r2]
    box_labels = [f'Q1\n(d{quartiles[0]:,.0f})', f'Q2\n({quartiles[0]:,.0f}-{quartiles[1]:,.0f})', 
                 f'Q3\n({quartiles[1]:,.0f}-{quartiles[2]:,.0f})', f'Q4\n(>{quartiles[2]:,.0f})']
    
    bp = ax3.boxplot(box_data, labels=box_labels, patch_artist=True)
    
    # Color the boxes
    colors = ['lightcoral', 'lightblue', 'lightgreen', 'gold']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax3.set_xlabel('Variant Count Quartiles')
    ax3.set_ylabel('R� Value')
    ax3.set_title(f'{chromosome}: R� by Variant Density')
    ax3.grid(True, alpha=0.3)
    
    # Subplot 4: Well-imputed efficiency
    ax4 = plt.subplot(3, 2, 4)
    
    # Plot well-imputed count vs total count
    scatter4 = ax4.scatter(variant_counts, well_imputed_counts, c=r2_values, 
                          cmap='RdYlGn', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Add diagonal line for perfect imputation
    max_variants = max(variant_counts)
    ax4.plot([0, max_variants], [0, max_variants], 'k--', alpha=0.5, 
            label='Perfect imputation')
    
    # Add trend line
    if len(variant_counts) > 1:
        z4 = np.polyfit(variant_counts, well_imputed_counts, 1)
        p4 = np.poly1d(z4)
        ax4.plot(variant_counts, p4(variant_counts), "r-", alpha=0.8, linewidth=2, 
                label=f'Actual trend (slope={z4[0]:.3f})')
    
    ax4.set_xlabel('Total Variants per Chunk')
    ax4.set_ylabel('Well-imputed Variants per Chunk')
    ax4.set_title(f'{chromosome}: Imputation Efficiency')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar4 = plt.colorbar(scatter4, ax=ax4)
    cbar4.set_label('Mean R�')
    
    # Subplot 5: Chunk performance ranking
    ax5 = plt.subplot(3, 2, 5)
    
    # Sort chunks by R� for ranking
    sorted_indices = np.argsort(r2_values)
    sorted_r2 = r2_values[sorted_indices]
    sorted_counts = variant_counts[sorted_indices]
    
    # Create ranking plot
    chunk_ranks = range(len(sorted_r2))
    bars = ax5.bar(chunk_ranks, sorted_r2, color=plt.cm.RdYlGn(sorted_r2), 
                  alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add quality thresholds
    ax5.axhline(y=0.3, color='orange', linestyle='--', alpha=0.7, label='Well-imputed (e0.3)')
    ax5.axhline(y=0.8, color='green', linestyle='--', alpha=0.7, label='High-quality (e0.8)')
    
    ax5.set_xlabel('Chunk Rank (sorted by R�)')
    ax5.set_ylabel('R� Value')
    ax5.set_title(f'{chromosome}: Chunk Performance Ranking')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Subplot 6: Summary statistics table
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('tight')
    ax6.axis('off')
    
    # Calculate summary statistics
    mean_r2 = np.mean(r2_values)
    median_r2 = np.median(r2_values)
    std_r2 = np.std(r2_values)
    min_r2 = np.min(r2_values)
    max_r2 = np.max(r2_values)
    
    mean_variants = np.mean(variant_counts)
    median_variants = np.median(variant_counts)
    total_well_imputed = np.sum(well_imputed_counts)
    overall_fraction = total_well_imputed / total_variants if total_variants > 0 else 0
    
    # Correlation metrics
    snp_r2_corr = np.corrcoef(variant_counts, r2_values)[0, 1] if len(variant_counts) > 1 else 0
    
    # Performance categories
    high_quality_chunks = np.sum(r2_values >= 0.8)
    well_imputed_chunks = np.sum(r2_values >= 0.3)
    poor_chunks = np.sum(r2_values < 0.3)
    
    # Create summary table
    summary_data = []
    summary_data.append(['Total Chunks', f"{len(chunk_details)}"])
    summary_data.append(['Total Variants', f"{total_variants:,}"])
    summary_data.append(['Mean Variants/Chunk', f"{mean_variants:,.0f}"])
    summary_data.append(['Median Variants/Chunk', f"{median_variants:,.0f}"])
    summary_data.append(['Mean R�', f"{mean_r2:.4f}"])
    summary_data.append(['R� Range', f"{min_r2:.3f} - {max_r2:.3f}"])
    summary_data.append(['R� Std Dev', f"{std_r2:.4f}"])
    summary_data.append(['SNP-R� Correlation', f"{snp_r2_corr:.3f}"])
    summary_data.append(['Overall Well-imputed', f"{overall_fraction:.1%}"])
    summary_data.append(['High Quality Chunks', f"{high_quality_chunks}/{len(chunk_details)} ({high_quality_chunks/len(chunk_details)*100:.1f}%)"])
    summary_data.append(['Poor Quality Chunks', f"{poor_chunks}/{len(chunk_details)} ({poor_chunks/len(chunk_details)*100:.1f}%)"])
    
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
    
    ax6.set_title(f'{chromosome}: SNP Count Analysis Summary', pad=20, fontweight='bold')
    
    # Overall title and layout
    fig.suptitle(f'Chromosome {chromosome.replace("chr", "")}: R� vs SNP Count Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.3, wspace=0.3)
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created chromosome R� vs SNP count plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level R� vs SNP count plots')
    parser.add_argument('--chr-summary', required=True, 
                       help='Path to chromosome summary JSON file')
    parser.add_argument('--output', required=True,
                       help='Output PDF file path')
    
    args = parser.parse_args()
    
    # Load data
    chr_summary = load_chromosome_summary(args.chr_summary)
    
    # Create plot
    create_chr_r2_snpcount_plot(chr_summary, args.output)

if __name__ == '__main__':
    main()