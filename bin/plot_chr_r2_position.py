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
        
        required_fields = ['chromosome', 'chunk_details']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"Missing required field: {field}")
        
        return data
    except Exception as e:
        print(f"Error loading JSON file {json_file}: {e}", file=sys.stderr)
        sys.exit(1)

def extract_positions_from_chunk_id(chunk_id):
    """Extract start and end positions from chunk ID."""
    try:
        parts = chunk_id.split('_')
        # Format: dataset_chr_start_end
        start_pos = int(parts[-2])
        end_pos = int(parts[-1])
        return start_pos, end_pos
    except (ValueError, IndexError):
        return None, None

def create_chr_r2_position_plot(chr_summary, output_file):
    """Create chromosome-level R² by genomic position visualization."""
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with 4 subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Extract data
    chromosome = chr_summary['chromosome']
    chunk_details = chr_summary['chunk_details']
    
    if not chunk_details:
        # Create empty plot with message
        ax = plt.subplot(1, 1, 1)
        ax.text(0.5, 0.5, f'No chunk data available for {chromosome}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(f'{chromosome}: R² by Genomic Position')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Process chunk data
    positions = []
    r2_values = []
    variant_counts = []
    well_imputed_counts = []
    
    for chunk in chunk_details:
        start_pos, end_pos = extract_positions_from_chunk_id(chunk['chunk_id'])
        if start_pos is not None and end_pos is not None:
            mid_pos = (start_pos + end_pos) / 2
            positions.append(mid_pos)
            r2_values.append(chunk['mean_r2'])
            variant_counts.append(chunk['variants'])
            well_imputed_counts.append(chunk['well_imputed'])
    
    if not positions:
        # Create empty plot with message
        ax = plt.subplot(1, 1, 1)
        ax.text(0.5, 0.5, f'Could not extract position data for {chromosome}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(f'{chromosome}: R² by Genomic Position')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Convert to arrays for easier handling
    positions = np.array(positions)
    r2_values = np.array(r2_values)
    variant_counts = np.array(variant_counts)
    well_imputed_counts = np.array(well_imputed_counts)
    
    # Subplot 1: R² by position (main plot)
    ax1 = plt.subplot(2, 2, 1)
    
    # Color points by R² value
    scatter = ax1.scatter(positions/1e6, r2_values, c=r2_values, 
                         cmap='RdYlGn', s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Add horizontal lines for quality thresholds
    ax1.axhline(y=0.3, color='orange', linestyle='--', alpha=0.7, label='Well-imputed (R²≥0.3)')
    ax1.axhline(y=0.8, color='red', linestyle='--', alpha=0.7, label='High-quality (R²≥0.8)')
    
    ax1.set_xlabel('Genomic Position (Mb)')
    ax1.set_ylabel('Mean R²')
    ax1.set_title(f'{chromosome}: R² by Genomic Position')
    ax1.set_ylim(0, 1)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('R² Value')
    
    # Subplot 2: Variant density by position
    ax2 = plt.subplot(2, 2, 2)
    
    # Bar plot of variant counts
    bar_width = (positions[1] - positions[0]) / 1e6 * 0.8 if len(positions) > 1 else 1
    bars = ax2.bar(positions/1e6, variant_counts, width=bar_width, 
                  alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    
    ax2.set_xlabel('Genomic Position (Mb)')
    ax2.set_ylabel('Number of Variants')
    ax2.set_title(f'{chromosome}: Variant Density by Position')
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Well-imputed fraction by position
    ax3 = plt.subplot(2, 2, 3)
    
    # Calculate well-imputed fraction
    well_imputed_fraction = well_imputed_counts / variant_counts
    
    # Color points by fraction
    scatter3 = ax3.scatter(positions/1e6, well_imputed_fraction, 
                          c=well_imputed_fraction, cmap='viridis', 
                          s=variant_counts/50, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    ax3.axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label='50% well-imputed')
    ax3.set_xlabel('Genomic Position (Mb)')
    ax3.set_ylabel('Well-imputed Fraction (R²≥0.3)')
    ax3.set_title(f'{chromosome}: Imputation Success by Position')
    ax3.set_ylim(0, 1)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar3 = plt.colorbar(scatter3, ax=ax3)
    cbar3.set_label('Well-imputed Fraction')
    
    # Subplot 4: Position-based statistics table and trend analysis
    ax4 = plt.subplot(2, 2, 4)
    
    # Calculate position-based statistics
    pos_mb = positions / 1e6
    pos_range = f"{pos_mb.min():.1f} - {pos_mb.max():.1f} Mb"
    mean_r2 = np.mean(r2_values)
    median_r2 = np.median(r2_values)
    r2_range = f"{r2_values.min():.3f} - {r2_values.max():.3f}"
    total_variants = np.sum(variant_counts)
    total_well_imputed = np.sum(well_imputed_counts)
    overall_fraction = total_well_imputed / total_variants if total_variants > 0 else 0
    
    # Check for trends
    correlation_r2_pos = np.corrcoef(positions, r2_values)[0, 1]
    correlation_density_pos = np.corrcoef(positions, variant_counts)[0, 1]
    
    # Create summary table
    summary_data = []
    summary_data.append(['Position Range', pos_range])
    summary_data.append(['Number of Chunks', f"{len(positions)}"])
    summary_data.append(['Total Variants', f"{total_variants:,}"])
    summary_data.append(['Mean R²', f"{mean_r2:.4f}"])
    summary_data.append(['Median R²', f"{median_r2:.4f}"])
    summary_data.append(['R² Range', r2_range])
    summary_data.append(['Well-imputed Overall', f"{overall_fraction:.1%}"])
    summary_data.append(['R² vs Position Corr.', f"{correlation_r2_pos:.3f}"])
    summary_data.append(['Density vs Position Corr.', f"{correlation_density_pos:.3f}"])
    
    # Best and worst performing chunks
    best_chunk_idx = np.argmax(r2_values)
    worst_chunk_idx = np.argmin(r2_values)
    
    summary_data.append(['Best Chunk R²', f"{r2_values[best_chunk_idx]:.3f} @ {pos_mb[best_chunk_idx]:.1f}Mb"])
    summary_data.append(['Worst Chunk R²', f"{r2_values[worst_chunk_idx]:.3f} @ {pos_mb[worst_chunk_idx]:.1f}Mb"])
    
    ax4.axis('tight')
    ax4.axis('off')
    
    table = ax4.table(cellText=summary_data,
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
    
    ax4.set_title(f'{chromosome}: Position Analysis Summary', pad=20, fontweight='bold')
    
    # Overall title and layout
    fig.suptitle(f'Chromosome {chromosome.replace("chr", "")}: R² by Genomic Position Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, hspace=0.3, wspace=0.3)
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created chromosome R² by position plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level R² by position plots')
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
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_r2_position.pdf"
    create_chr_r2_position_plot(chr_summary, output_file)

if __name__ == '__main__':
    main()