#!/usr/bin/env python3
"""
Generate chromosome-level plots for imputation quality metrics.
Creates comprehensive visualizations of imputation performance across chunks.
"""

import json
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


def load_chromosome_summary(summary_path):
    """Load chromosome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)


def plot_chunk_performance(ax, chr_summary):
    """Plot chunk-level performance metrics."""
    if not chr_summary.get('chunk_details'):
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    chunks = chr_summary['chunk_details']
    chunk_ids = [c['chunk_id'] for c in chunks]
    mean_r2 = [c.get('mean_r2', 0) for c in chunks]
    well_imputed = [c.get('well_imputed', 0) for c in chunks]
    
    # Create bar plot
    x = np.arange(len(chunks))
    width = 0.35
    
    # Normalize well_imputed for dual axis
    max_well_imputed = max(well_imputed) if well_imputed and max(well_imputed) > 0 else 1
    well_imputed_norm = [w/max_well_imputed for w in well_imputed]
    
    bars1 = ax.bar(x - width/2, mean_r2, width, label='Mean R²', alpha=0.8)
    bars2 = ax.bar(x + width/2, well_imputed_norm, width, label='Well-imputed (normalized)', alpha=0.8)
    
    ax.set_xlabel('Chunk')
    ax.set_ylabel('Score')
    ax.set_title(f'Chunk Performance - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticks(x)
    
    # Rotate labels if many chunks
    if len(chunks) > 10:
        ax.set_xticklabels([f"C{i+1}" for i in range(len(chunks))], rotation=45, ha='right')
    else:
        ax.set_xticklabels([c.split('_')[-1][:8] for c in chunk_ids], rotation=45, ha='right')
    
    ax.legend()
    ax.set_ylim(0, 1.1)
    
    # Add value labels on bars if few chunks
    if len(chunks) <= 20:
        for bar, val in zip(bars1, mean_r2):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=8)


def plot_maf_distribution(ax, chr_summary):
    """Plot MAF distribution."""
    maf_bins = chr_summary.get('maf_bins', {})
    
    if not maf_bins:
        ax.text(0.5, 0.5, 'No MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort bins for consistent ordering
    bin_names = sorted(maf_bins.keys())
    counts = [maf_bins[bin_name] for bin_name in bin_names]
    
    # Check if all counts are zero
    if all(c == 0 for c in counts):
        ax.text(0.5, 0.5, 'No variants in MAF bins', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'MAF Distribution - Chromosome {chr_summary.get("chromosome", "?")}')
        return
    
    # Create pie chart
    colors = plt.cm.Set3(np.linspace(0, 1, len(bin_names)))
    wedges, texts, autotexts = ax.pie(counts, labels=bin_names, colors=colors,
                                       autopct='%1.1f%%', startangle=90)
    
    ax.set_title(f'MAF Distribution - Chromosome {chr_summary.get("chromosome", "?")}')
    
    # Make percentage text more readable
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_weight('bold')
        autotext.set_fontsize(9)


def plot_r2_by_maf(ax, chr_summary):
    """Plot R² by MAF bins."""
    r2_by_maf = chr_summary.get('mean_r2_by_maf', {})
    
    if not r2_by_maf:
        ax.text(0.5, 0.5, 'No R² by MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort bins for consistent ordering
    bin_names = sorted(r2_by_maf.keys())
    mean_values = []
    std_values = []
    n_values = []
    
    for bin_name in bin_names:
        stats = r2_by_maf[bin_name]
        mean_values.append(stats.get('mean', 0))
        std_values.append(stats.get('std', 0))
        n_values.append(stats.get('n', 0))
    
    # Create bar plot with error bars
    x = np.arange(len(bin_names))
    bars = ax.bar(x, mean_values, yerr=std_values, capsize=5, alpha=0.8,
                  color='steelblue', edgecolor='navy', linewidth=1.5)
    
    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Mean R²')
    ax.set_title(f'R² by MAF - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha='right')
    ax.set_ylim(0, max(mean_values) * 1.2 if mean_values else 1)
    
    # Add sample size annotations
    for i, (bar, n) in enumerate(zip(bars, n_values)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + std_values[i] + 0.02,
               f'n={n:,}', ha='center', va='bottom', fontsize=8)
    
    # Add horizontal line for overall mean
    if chr_summary.get('mean_r2'):
        ax.axhline(y=chr_summary['mean_r2'], color='red', linestyle='--', 
                  alpha=0.7, label=f'Overall Mean: {chr_summary["mean_r2"]:.3f}')
        ax.legend()


def plot_info_score_distribution(ax, chr_summary):
    """Plot info score distribution."""
    info_stats = chr_summary.get('info_score_stats', {})
    
    if not info_stats:
        ax.text(0.5, 0.5, 'No info score data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Create box plot representation using statistics
    data = {
        'Min': info_stats.get('min', 0),
        'Q25': info_stats.get('q25', 0),
        'Median': info_stats.get('median', 0),
        'Mean': info_stats.get('mean', 0),
        'Q75': info_stats.get('q75', 0),
        'Max': info_stats.get('max', 0)
    }
    
    # Create horizontal bar chart of statistics
    stats_names = list(data.keys())
    stats_values = list(data.values())
    
    y_pos = np.arange(len(stats_names))
    colors = ['red' if name == 'Mean' else 'steelblue' for name in stats_names]
    
    bars = ax.barh(y_pos, stats_values, color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(stats_names)
    ax.set_xlabel('Info Score (R²)')
    ax.set_title(f'Info Score Statistics - Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_xlim(0, 1)
    
    # Add value labels
    for bar, val in zip(bars, stats_values):
        width = bar.get_width()
        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2.,
               f'{val:.3f}', ha='left', va='center')
    
    # Add threshold lines
    ax.axvline(x=0.3, color='orange', linestyle='--', alpha=0.5, label='Poor threshold')
    ax.axvline(x=0.8, color='green', linestyle='--', alpha=0.5, label='Good threshold')
    
    # Add counts above thresholds
    if 'above_0.3' in info_stats and 'above_0.8' in info_stats:
        text = f"Above 0.3: {info_stats['above_0.3']:,}\nAbove 0.8: {info_stats['above_0.8']:,}"
        ax.text(0.95, 0.95, text, transform=ax.transAxes,
               ha='right', va='top', fontsize=9,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(loc='lower right')


def create_chromosome_summary_plots(chr_summary, output_prefix, ref_name, dataset, chromosome):
    """Create comprehensive chromosome-level plots."""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Add main title
    fig.suptitle(f'Chromosome {chromosome} Imputation Summary - Dataset: {dataset}', 
                 fontsize=16, fontweight='bold')
    
    # Create subplots
    ax1 = plt.subplot(2, 2, 1)
    plot_chunk_performance(ax1, chr_summary)
    
    ax2 = plt.subplot(2, 2, 2)
    plot_maf_distribution(ax2, chr_summary)
    
    ax3 = plt.subplot(2, 2, 3)
    plot_r2_by_maf(ax3, chr_summary)
    
    ax4 = plt.subplot(2, 2, 4)
    plot_info_score_distribution(ax4, chr_summary)
    
    # Add summary text box with safe formatting
    chunks = chr_summary.get('chunks_processed', 0)
    variants = chr_summary.get('total_variants', 0)
    well_imp = chr_summary.get('well_imputed_variants', 0)
    mean_r2 = chr_summary.get('mean_r2', 0) or 0
    mean_conc = chr_summary.get('mean_concordance', 0) or 0
    
    summary_text = f"""
    Chunks Processed: {chunks}
    Total Variants: {variants:,}
    Well-Imputed: {well_imp:,}
    Mean R²: {mean_r2:.4f}
    Mean Concordance: {mean_conc:.4f}
    Reference: {ref_name}
    """
    
    fig.text(0.02, 0.02, summary_text.strip(), fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])
    
    return fig


def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level plots')
    parser.add_argument('--chr-summary', required=True,
                       help='Chromosome summary JSON file')
    parser.add_argument('--chunk-plots', nargs='+',
                       help='Chunk plot files (optional, for reference)')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    parser.add_argument('--chromosome', required=True,
                       help='Chromosome identifier')
    
    args = parser.parse_args()
    
    # Load chromosome summary
    logger.info(f"Loading chromosome summary from {args.chr_summary}")
    chr_summary = load_chromosome_summary(args.chr_summary)
    
    # Create plots
    logger.info(f"Creating chromosome-level plots for chromosome {args.chromosome}")
    fig = create_chromosome_summary_plots(
        chr_summary, 
        args.output_prefix,
        args.ref_name,
        args.dataset,
        args.chromosome
    )
    
    # Save to PDF
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_plots.pdf"
    logger.info(f"Saving plots to {output_file}")
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
        
        # Add metadata
        d = pdf.infodict()
        d['Title'] = f'Chromosome {args.chromosome} Imputation Summary'
        d['Author'] = 'Genotype Imputation Pipeline'
        d['Subject'] = f'Imputation QC for {args.dataset}'
        d['Keywords'] = f'Imputation QC Chromosome{args.chromosome} {args.ref_name}'
    
    plt.close(fig)
    
    logger.info(f"Chromosome plots saved successfully to {output_file}")
    
    # Print summary
    print(f"\nChromosome {args.chromosome} Plot Summary:")
    print(f"  Dataset: {args.dataset}")
    print(f"  Reference: {args.ref_name}")
    print(f"  Chunks: {chr_summary.get('chunks_processed', 0)}")
    mean_r2_val = chr_summary.get('mean_r2', 0) or 0
    print(f"  Mean R²: {mean_r2_val:.4f}")
    print(f"  Output: {output_file}")


if __name__ == '__main__':
    main()