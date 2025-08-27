#!/usr/bin/env python3
"""
Generate genome-wide plots for imputation quality metrics.
Creates comprehensive visualizations of imputation performance across all chromosomes.
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


def load_genome_summary(summary_path):
    """Load genome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)


def chromosome_sort_key(chromosome):
    """Generate sort key for chromosome names."""
    chr_str = str(chromosome).replace('chr', '').upper()
    if chr_str.isdigit():
        return (0, int(chr_str))
    elif chr_str == 'X':
        return (1, 23)
    elif chr_str == 'Y':
        return (1, 24)
    elif chr_str in ['MT', 'M']:
        return (1, 25)
    else:
        return (2, chr_str)


def plot_chromosome_comparison(ax, genome_summary):
    """Plot comparison of metrics across chromosomes."""
    chr_details = genome_summary.get('chromosome_details', [])
    
    if not chr_details:
        ax.text(0.5, 0.5, 'No chromosome data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chromosomes
    chr_details_sorted = sorted(chr_details, key=lambda x: chromosome_sort_key(x['chromosome']))
    
    chromosomes = [str(c['chromosome']) for c in chr_details_sorted]
    mean_r2 = [(c.get('mean_r2', 0) or 0) for c in chr_details_sorted]
    concordance = [(c.get('mean_concordance', 0) or 0) for c in chr_details_sorted]
    
    x = np.arange(len(chromosomes))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, mean_r2, width, label='Mean R²', alpha=0.8, color='steelblue')
    bars2 = ax.bar(x + width/2, concordance, width, label='Concordance', alpha=0.8, color='coral')
    
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Score')
    ax.set_title('Performance Metrics by Chromosome')
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45 if len(chromosomes) > 15 else 0, ha='right' if len(chromosomes) > 15 else 'center')
    ax.legend()
    ax.set_ylim(0, 1.1)
    
    # Add genome-wide average lines
    if genome_summary.get('mean_r2'):
        ax.axhline(y=genome_summary['mean_r2'], color='blue', linestyle='--', 
                  alpha=0.5, linewidth=1)
    if genome_summary.get('mean_concordance'):
        ax.axhline(y=genome_summary['mean_concordance'], color='red', linestyle='--', 
                  alpha=0.5, linewidth=1)


def plot_variant_distribution(ax, genome_summary):
    """Plot variant distribution across chromosomes."""
    chr_details = genome_summary.get('chromosome_details', [])
    
    if not chr_details:
        ax.text(0.5, 0.5, 'No chromosome data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chromosomes
    chr_details_sorted = sorted(chr_details, key=lambda x: chromosome_sort_key(x['chromosome']))
    
    chromosomes = [str(c['chromosome']) for c in chr_details_sorted]
    total_variants = [c.get('variants', 0) for c in chr_details_sorted]
    well_imputed = [c.get('well_imputed', 0) for c in chr_details_sorted]
    
    x = np.arange(len(chromosomes))
    
    # Create stacked bar chart
    bars1 = ax.bar(x, well_imputed, label='Well-imputed', alpha=0.8, color='green')
    bars2 = ax.bar(x, [t-w for t, w in zip(total_variants, well_imputed)], 
                   bottom=well_imputed, label='Other variants', alpha=0.8, color='lightgray')
    
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Number of Variants')
    ax.set_title('Variant Distribution Across Chromosomes')
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45 if len(chromosomes) > 15 else 0, ha='right' if len(chromosomes) > 15 else 'center')
    ax.legend()
    
    # Add percentage labels for well-imputed
    for i, (bar, total, well) in enumerate(zip(bars1, total_variants, well_imputed)):
        if total > 0:
            percentage = (well / total) * 100
            if i % max(1, len(chromosomes) // 10) == 0:  # Show labels for subset if many chromosomes
                ax.text(bar.get_x() + bar.get_width()/2., total + max(total_variants)*0.01,
                       f'{percentage:.0f}%', ha='center', va='bottom', fontsize=8)


def plot_genome_wide_maf(ax, genome_summary):
    """Plot genome-wide MAF distribution."""
    maf_bins = genome_summary.get('maf_bins', {})
    
    if not maf_bins:
        ax.text(0.5, 0.5, 'No MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort bins for consistent ordering
    bin_names = sorted(maf_bins.keys())
    counts = [maf_bins[bin_name] for bin_name in bin_names]
    total = sum(counts)
    
    # Check for empty data
    if total == 0:
        ax.text(0.5, 0.5, 'No MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('MAF Distribution (Genome-wide)')
        return
    
    percentages = [(c/total)*100 for c in counts]
    
    # Create bar chart
    x = np.arange(len(bin_names))
    bars = ax.bar(x, percentages, alpha=0.8, color='teal')
    
    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Percentage of Variants (%)')
    ax.set_title('Genome-wide MAF Distribution')
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha='right')
    
    # Add count labels
    for bar, count, pct in zip(bars, counts, percentages):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
               f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=8)
    
    # Add total variants annotation
    ax.text(0.98, 0.98, f'Total variants: {total:,}', transform=ax.transAxes,
           ha='right', va='top', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


def plot_genome_r2_by_maf(ax, genome_summary):
    """Plot genome-wide R² by MAF."""
    r2_by_maf = genome_summary.get('mean_r2_by_maf', {})
    
    if not r2_by_maf:
        ax.text(0.5, 0.5, 'No R² by MAF data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort bins for consistent ordering
    bin_names = sorted(r2_by_maf.keys())
    mean_values = [r2_by_maf[bin_name]['mean'] for bin_name in bin_names]
    n_variants = [r2_by_maf[bin_name]['n_variants'] for bin_name in bin_names]
    
    # Create line plot with markers
    x = np.arange(len(bin_names))
    line = ax.plot(x, mean_values, 'o-', linewidth=2, markersize=8, 
                  color='darkblue', label='Mean R²')
    
    ax.set_xlabel('MAF Bin')
    ax.set_ylabel('Mean R²')
    ax.set_title('Genome-wide R² by MAF Category')
    ax.set_xticks(x)
    ax.set_xticklabels(bin_names, rotation=45, ha='right')
    ax.set_ylim(0, max(mean_values) * 1.1 if mean_values else 1)
    ax.grid(True, alpha=0.3)
    
    # Add value labels
    for i, (val, n) in enumerate(zip(mean_values, n_variants)):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', va='bottom', fontsize=9)
    
    # Add overall mean line
    if genome_summary.get('mean_r2'):
        ax.axhline(y=genome_summary['mean_r2'], color='red', linestyle='--', 
                  alpha=0.7, label=f'Overall Mean: {genome_summary["mean_r2"]:.3f}')
    
    ax.legend()
    
    # Add secondary information
    info_text = f"Total chromosomes: {genome_summary.get('chromosomes_processed', 0)}\n"
    info_text += f"Total chunks: {genome_summary.get('total_chunks', 0)}"
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
           ha='left', va='top', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))


def plot_imputation_summary_stats(ax, genome_summary):
    """Plot summary statistics as a dashboard."""
    ax.axis('off')
    
    # Create dashboard-style display
    # Safely extract values with defaults
    imp_rate = genome_summary.get('imputation_rate', 0) or 0
    well_imp_rate = genome_summary.get('well_imputed_rate', 0) or 0
    mean_r2 = genome_summary.get('mean_r2', 0) or 0
    mean_conc = genome_summary.get('mean_concordance', 0) or 0
    
    stats = [
        ('Chromosomes', genome_summary.get('chromosomes_processed', 0)),
        ('Total Chunks', f"{genome_summary.get('total_chunks', 0):,}"),
        ('Total Variants', f"{genome_summary.get('total_variants', 0):,}"),
        ('Well-Imputed', f"{genome_summary.get('well_imputed_variants', 0):,}"),
        ('Imputation Rate', f"{imp_rate*100:.1f}%"),
        ('Well-Imputed Rate', f"{well_imp_rate*100:.1f}%"),
        ('Mean R²', f"{mean_r2:.4f}"),
        ('Mean Concordance', f"{mean_conc:.4f}"),
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
    
    ax.set_title('Genome-wide Summary Statistics', fontsize=14, weight='bold', pad=20)
    
    # Add dataset info at bottom
    dataset_info = f"Dataset: {genome_summary.get('dataset', 'Unknown')} | Reference: {genome_summary.get('ref_name', 'Unknown')}"
    ax.text(0.5, 0.05, dataset_info, transform=ax.transAxes,
           ha='center', fontsize=11, style='italic')


def plot_chromosome_heatmap(ax, genome_summary):
    """Create heatmap of metrics across chromosomes."""
    chr_details = genome_summary.get('chromosome_details', [])
    
    if not chr_details:
        ax.text(0.5, 0.5, 'No chromosome data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chromosomes
    chr_details_sorted = sorted(chr_details, key=lambda x: chromosome_sort_key(x['chromosome']))
    
    # Prepare data for heatmap
    chromosomes = [str(c['chromosome']) for c in chr_details_sorted]
    metrics = ['Chunks', 'Variants (k)', 'Well-Imp (%)', 'Mean R²', 'Concordance']
    
    # Normalize data for heatmap
    data = []
    for chr_detail in chr_details_sorted:
        # Ensure all values are numeric with proper defaults
        chunks = chr_detail.get('chunks', 0) or 0
        variants = chr_detail.get('variants', 0) or 0
        well_imputed = chr_detail.get('well_imputed', 0) or 0
        mean_r2 = chr_detail.get('mean_r2', 0) or 0
        mean_concordance = chr_detail.get('mean_concordance', 0) or 0
        
        row = [
            float(chunks),
            float(variants) / 1000.0,  # Convert to thousands
            (float(well_imputed) / float(variants)) * 100.0 if variants > 0 else 0.0,
            float(mean_r2),
            float(mean_concordance)
        ]
        data.append(row)
    
    data = np.array(data, dtype=float).T
    
    # Normalize each metric to 0-1 range for visualization
    data_norm = np.zeros_like(data)
    for i in range(data.shape[0]):
        if data[i].size > 0:
            min_val = np.min(data[i])
            max_val = np.max(data[i])
            if min_val is not None and max_val is not None and max_val > min_val:
                data_norm[i] = (data[i] - min_val) / (max_val - min_val)
            else:
                data_norm[i] = 0.5
        else:
            data_norm[i] = 0.5
    
    # Create heatmap
    im = ax.imshow(data_norm, cmap='YlOrRd', aspect='auto')
    
    # Set ticks
    ax.set_xticks(np.arange(len(chromosomes)))
    ax.set_yticks(np.arange(len(metrics)))
    ax.set_xticklabels(chromosomes, rotation=45 if len(chromosomes) > 15 else 0, 
                       ha='right' if len(chromosomes) > 15 else 'center')
    ax.set_yticklabels(metrics)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Normalized Value', rotation=270, labelpad=15)
    
    ax.set_title('Chromosome Metrics Heatmap')
    
    # Add text annotations for small datasets
    if len(chromosomes) <= 10:
        for i in range(len(metrics)):
            for j in range(len(chromosomes)):
                text = ax.text(j, i, f'{data[i, j]:.1f}',
                             ha='center', va='center', color='black' if data_norm[i, j] < 0.5 else 'white',
                             fontsize=8)


def create_genome_summary_plots(genome_summary, output_prefix, ref_name, dataset):
    """Create comprehensive genome-wide plots."""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 16))
    
    # Add main title
    fig.suptitle(f'Genome-wide Imputation Summary - Dataset: {dataset}', 
                 fontsize=18, fontweight='bold')
    
    # Create 2x3 grid of subplots
    ax1 = plt.subplot(2, 3, 1)
    plot_chromosome_comparison(ax1, genome_summary)
    
    ax2 = plt.subplot(2, 3, 2)
    plot_variant_distribution(ax2, genome_summary)
    
    ax3 = plt.subplot(2, 3, 3)
    plot_genome_wide_maf(ax3, genome_summary)
    
    ax4 = plt.subplot(2, 3, 4)
    plot_genome_r2_by_maf(ax4, genome_summary)
    
    ax5 = plt.subplot(2, 3, 5)
    plot_imputation_summary_stats(ax5, genome_summary)
    
    ax6 = plt.subplot(2, 3, 6)
    plot_chromosome_heatmap(ax6, genome_summary)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    
    return fig


def main():
    parser = argparse.ArgumentParser(description='Generate genome-wide plots')
    parser.add_argument('--genome-summary', required=True,
                       help='Genome summary JSON file')
    parser.add_argument('--chr-plots', nargs='+',
                       help='Chromosome plot files (optional, for reference)')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    
    args = parser.parse_args()
    
    # Load genome summary
    logger.info(f"Loading genome summary from {args.genome_summary}")
    genome_summary = load_genome_summary(args.genome_summary)
    
    # Create plots
    logger.info(f"Creating genome-wide plots for dataset {args.dataset}")
    fig = create_genome_summary_plots(
        genome_summary,
        args.output_prefix,
        args.ref_name,
        args.dataset
    )
    
    # Save to PDF
    output_file = f"{args.output_prefix}_{args.ref_name}.genome_plots.pdf"
    logger.info(f"Saving plots to {output_file}")
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
        
        # Add metadata
        d = pdf.infodict()
        d['Title'] = f'Genome-wide Imputation Summary for {args.dataset}'
        d['Author'] = 'Genotype Imputation Pipeline'
        d['Subject'] = f'Genome-wide Imputation QC for {args.dataset}'
        d['Keywords'] = f'Imputation QC Genome {args.ref_name} {args.dataset}'
    
    plt.close(fig)
    
    logger.info(f"Genome-wide plots saved successfully to {output_file}")
    
    # Print summary
    print(f"\nGenome-wide Plot Summary:")
    print(f"  Dataset: {args.dataset}")
    print(f"  Reference: {args.ref_name}")
    print(f"  Chromosomes: {genome_summary.get('chromosomes_processed', 0)}")
    print(f"  Total variants: {genome_summary.get('total_variants', 0):,}")
    mean_r2_val = genome_summary.get('mean_r2', 0) or 0
    print(f"  Mean R²: {mean_r2_val:.4f}")
    print(f"  Output: {output_file}")


if __name__ == '__main__':
    main()