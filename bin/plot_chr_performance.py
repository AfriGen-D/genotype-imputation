#!/usr/bin/env python3
"""
Generate chromosome-level performance plots aggregating chunk data.
"""

import json
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


def load_chr_summary(summary_path):
    """Load chromosome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)


def plot_chunk_performance_comparison(ax, chr_summary):
    """Plot performance metrics across chunks."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Extract metrics per chunk
    chunk_names = []
    r2_values = []
    concordance_values = []
    well_imputed_rates = []
    
    for chunk in chunk_details:
        chunk_id = chunk.get('chunk_id', 'Unknown')
        # Extract just the position part for cleaner labels
        if '_' in chunk_id:
            parts = chunk_id.split('_')
            if len(parts) >= 4:
                chunk_names.append(f"{parts[-2]}-{parts[-1]}")
            else:
                chunk_names.append(chunk_id)
        else:
            chunk_names.append(chunk_id)
        
        r2_values.append(chunk.get('mean_r2', 0) or 0)
        concordance_values.append(chunk.get('mean_concordance', 0) or 0)
        
        total_vars = chunk.get('total_variants', 0) or 1
        well_imp = chunk.get('well_imputed', 0) or 0
        well_imputed_rates.append((well_imp / total_vars) if total_vars > 0 else 0)
    
    x = np.arange(len(chunk_names))
    width = 0.25
    
    # Create grouped bar chart
    bars1 = ax.bar(x - width, r2_values, width, label='Mean R²', alpha=0.8, color='steelblue')
    bars2 = ax.bar(x, concordance_values, width, label='Concordance', alpha=0.8, color='coral')
    bars3 = ax.bar(x + width, well_imputed_rates, width, label='Well-Imputed Rate', alpha=0.8, color='green')
    
    ax.set_xlabel('Chunk Position')
    ax.set_ylabel('Score')
    ax.set_title(f'Performance Metrics by Chunk - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticks(x)
    ax.set_xticklabels(chunk_names, rotation=45, ha='right')
    ax.legend()
    ax.set_ylim(0, 1.1)
    
    # Add chromosome mean lines
    if chr_summary.get('mean_r2'):
        ax.axhline(y=chr_summary['mean_r2'], color='blue', linestyle='--', alpha=0.5, linewidth=1)
    if chr_summary.get('mean_concordance'):
        ax.axhline(y=chr_summary['mean_concordance'], color='red', linestyle='--', alpha=0.5, linewidth=1)


def plot_variant_distribution(ax, chr_summary):
    """Plot variant distribution across chunks."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    chunk_names = []
    total_variants = []
    well_imputed = []
    genotyped = []
    
    for chunk in chunk_details:
        chunk_id = chunk.get('chunk_id', 'Unknown')
        if '_' in chunk_id:
            parts = chunk_id.split('_')
            if len(parts) >= 4:
                chunk_names.append(f"{parts[-2]}-{parts[-1]}")
            else:
                chunk_names.append(chunk_id)
        else:
            chunk_names.append(chunk_id)
        
        total_variants.append(chunk.get('total_variants', 0) or 0)
        well_imputed.append(chunk.get('well_imputed', 0) or 0)
        genotyped.append(chunk.get('genotyped', 0) or 0)
    
    x = np.arange(len(chunk_names))
    
    # Create stacked bar chart
    bars1 = ax.bar(x, genotyped, label='Genotyped', alpha=0.8, color='blue')
    bars2 = ax.bar(x, well_imputed, bottom=genotyped, label='Well-Imputed', alpha=0.8, color='green')
    remaining = [t - w - g for t, w, g in zip(total_variants, well_imputed, genotyped)]
    bars3 = ax.bar(x, remaining, bottom=[w + g for w, g in zip(well_imputed, genotyped)], 
                   label='Other', alpha=0.8, color='lightgray')
    
    ax.set_xlabel('Chunk Position')
    ax.set_ylabel('Number of Variants')
    ax.set_title(f'Variant Distribution - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticks(x)
    ax.set_xticklabels(chunk_names, rotation=45, ha='right')
    ax.legend()


def plot_maf_distribution(ax, chr_summary):
    """Plot MAF distribution for chromosome."""
    maf_bins = chr_summary.get('maf_bins', {})
    
    if not maf_bins:
        ax.text(0.5, 0.5, 'No MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    bin_names = sorted(maf_bins.keys())
    counts = [maf_bins[bin_name] for bin_name in bin_names]
    total = sum(counts)
    percentages = [(c/total)*100 if total > 0 else 0 for c in counts]
    
    # Create bar chart
    bars = ax.bar(bin_names, percentages, alpha=0.8, color='teal')
    
    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Percentage of Variants (%)')
    ax.set_title(f'MAF Distribution - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticklabels(bin_names, rotation=45, ha='right')
    
    # Add count labels
    for bar, count, pct in zip(bars, counts, percentages):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
               f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=8)


def plot_r2_by_maf(ax, chr_summary):
    """Plot R² by MAF for chromosome."""
    r2_by_maf = chr_summary.get('mean_r2_by_maf', {})
    
    if not r2_by_maf:
        ax.text(0.5, 0.5, 'No R² by MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    bin_names = sorted(r2_by_maf.keys())
    mean_values = []
    std_values = []
    n_variants = []
    
    for bin_name in bin_names:
        data = r2_by_maf[bin_name]
        mean_values.append(data.get('mean', 0) or 0)
        std_values.append(data.get('std', 0) or 0)
        n_variants.append(data.get('n_variants', 0) or 0)
    
    x = np.arange(len(bin_names))
    
    # Create line plot with error bars
    ax.errorbar(x, mean_values, yerr=std_values, fmt='o-', linewidth=2, 
                markersize=8, capsize=5, capthick=2, color='darkblue', label='Mean R² ± SD')
    
    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Mean R²')
    ax.set_title(f'R² by MAF Category - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha='right')
    ax.set_ylim(0, max(mean_values) * 1.2 if mean_values else 1)
    ax.grid(True, alpha=0.3)
    
    # Add sample size annotations
    for i, (val, n) in enumerate(zip(mean_values, n_variants)):
        ax.text(i, val + std_values[i] + 0.02, f'n={n:,}', 
                ha='center', va='bottom', fontsize=8)
    
    # Add overall mean line
    if chr_summary.get('mean_r2'):
        ax.axhline(y=chr_summary['mean_r2'], color='red', linestyle='--', 
                  alpha=0.7, label=f'Overall Mean: {chr_summary["mean_r2"]:.3f}')
    
    ax.legend()


def plot_summary_stats(ax, chr_summary):
    """Plot summary statistics dashboard."""
    ax.axis('off')
    
    # Create dashboard display
    stats = [
        ('Chromosome', chr_summary.get('chromosome', 'Unknown')),
        ('Chunks Processed', chr_summary.get('chunks_processed', 0)),
        ('Total Variants', f"{chr_summary.get('total_variants', 0):,}"),
        ('Well-Imputed', f"{chr_summary.get('well_imputed_variants', 0):,}"),
        ('Mean R²', f"{(chr_summary.get('mean_r2', 0) or 0):.4f}"),
        ('Mean Concordance', f"{(chr_summary.get('mean_concordance', 0) or 0):.4f}"),
        ('Dataset', chr_summary.get('dataset', 'Unknown')),
        ('Reference', chr_summary.get('ref_name', 'Unknown')),
    ]
    
    # Create grid layout
    n_cols = 2
    n_rows = (len(stats) + n_cols - 1) // n_cols
    
    for i, (label, value) in enumerate(stats):
        row = i // n_cols
        col = i % n_cols
        
        x = 0.1 + col * 0.45
        y = 0.85 - row * 0.2
        
        # Draw box
        rect = plt.Rectangle((x-0.05, y-0.08), 0.35, 0.15, 
                            facecolor='lightsteelblue', alpha=0.5,
                            transform=ax.transAxes)
        ax.add_patch(rect)
        
        # Add text
        ax.text(x, y, label, fontsize=10, weight='bold',
               transform=ax.transAxes)
        ax.text(x, y-0.04, str(value), fontsize=14,
               transform=ax.transAxes)
    
    ax.set_title('Chromosome Summary Statistics', fontsize=14, weight='bold', pad=20)


def create_chromosome_performance_plots(chr_summary, output_prefix):
    """Create comprehensive chromosome-level performance plots."""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 12))
    
    # Add main title
    fig.suptitle(f'Chromosome {chr_summary.get("chromosome", "?")} Performance Summary - Dataset: {chr_summary.get("dataset", "Unknown")}', 
                 fontsize=16, fontweight='bold')
    
    # Create 2x3 grid
    ax1 = plt.subplot(2, 3, 1)
    plot_chunk_performance_comparison(ax1, chr_summary)
    
    ax2 = plt.subplot(2, 3, 2)
    plot_variant_distribution(ax2, chr_summary)
    
    ax3 = plt.subplot(2, 3, 3)
    plot_maf_distribution(ax3, chr_summary)
    
    ax4 = plt.subplot(2, 3, 4)
    plot_r2_by_maf(ax4, chr_summary)
    
    ax5 = plt.subplot(2, 3, 5)
    plot_summary_stats(ax5, chr_summary)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    
    return fig


def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level performance plots')
    parser.add_argument('--chr-summary', required=True,
                       help='Chromosome summary JSON file')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    parser.add_argument('--chromosome', required=True,
                       help='Chromosome name')
    
    args = parser.parse_args()
    
    # Load chromosome summary
    logger.info(f"Loading chromosome summary from {args.chr_summary}")
    chr_summary = load_chr_summary(args.chr_summary)
    
    # Create plots
    logger.info(f"Creating performance plots for chromosome {args.chromosome}")
    fig = create_chromosome_performance_plots(chr_summary, args.output_prefix)
    
    # Save to PDF
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_performance.pdf"
    logger.info(f"Saving plots to {output_file}")
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
        
        # Add metadata
        d = pdf.infodict()
        d['Title'] = f'Chromosome {args.chromosome} Performance Summary'
        d['Author'] = 'Genotype Imputation Pipeline'
        d['Subject'] = f'Chromosome Performance Metrics for {args.dataset}'
        d['Keywords'] = f'Imputation Performance Chromosome {args.chromosome}'
    
    plt.close(fig)
    
    logger.info(f"Performance plots saved successfully to {output_file}")
    
    # Print summary
    print(f"\nChromosome Performance Summary:")
    print(f"  Chromosome: {args.chromosome}")
    print(f"  Dataset: {args.dataset}")
    print(f"  Chunks: {chr_summary.get('chunks_processed', 0)}")
    print(f"  Output: {output_file}")


if __name__ == '__main__':
    main()