#!/usr/bin/env python3
"""
Generate chromosome-level R² by position plots aggregating chunk data.
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


def plot_r2_along_chromosome(ax, chr_summary):
    """Plot R² values along the chromosome."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chunks by position
    sorted_chunks = sorted(chunk_details, key=lambda x: x.get('start_pos', 0))
    
    positions = []
    r2_values = []
    chunk_boundaries = []
    
    for chunk in sorted_chunks:
        # Use midpoint of chunk for position
        start_pos = chunk.get('start_pos', 0)
        end_pos = chunk.get('end_pos', start_pos + 1000000)
        midpoint = (start_pos + end_pos) / 2
        
        positions.append(midpoint / 1e6)  # Convert to Mb
        r2_values.append(chunk.get('mean_r2', 0) or 0)
        chunk_boundaries.append((start_pos / 1e6, end_pos / 1e6))
    
    # Plot line with markers
    ax.plot(positions, r2_values, 'o-', linewidth=2, markersize=8, 
            color='darkblue', label='Mean R² per chunk')
    
    # Add shaded regions for chunks
    for i, (start, end) in enumerate(chunk_boundaries):
        ax.axvspan(start, end, alpha=0.1, color='gray')
    
    # Add chromosome mean line
    if chr_summary.get('mean_r2'):
        ax.axhline(y=chr_summary['mean_r2'], color='red', linestyle='--', 
                  alpha=0.7, label=f'Chromosome Mean: {chr_summary["mean_r2"]:.3f}')
    
    ax.set_xlabel('Position (Mb)')
    ax.set_ylabel('Mean R²')
    ax.set_title(f'R² Distribution Along Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)
    ax.legend()


def plot_r2_density_by_position(ax, chr_summary):
    """Plot R² density heatmap along chromosome."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chunks by position
    sorted_chunks = sorted(chunk_details, key=lambda x: x.get('start_pos', 0))
    
    # Create 2D density plot data
    positions = []
    r2_distributions = []
    
    for chunk in sorted_chunks:
        start_pos = chunk.get('start_pos', 0)
        end_pos = chunk.get('end_pos', start_pos + 1000000)
        midpoint = (start_pos + end_pos) / 2 / 1e6  # Mb
        
        # Get R² distribution for this chunk
        r2_dist = chunk.get('r2_distribution', {})
        if r2_dist:
            # Create histogram bins
            bins = np.linspace(0, 1, 21)  # 20 bins from 0 to 1
            hist = np.zeros(20)
            
            for bin_name, count in r2_dist.items():
                # Parse bin name like "0.0-0.1" to get index
                try:
                    low = float(bin_name.split('-')[0])
                    idx = int(low * 20)
                    if 0 <= idx < 20:
                        hist[idx] = count
                except:
                    continue
            
            positions.append(midpoint)
            r2_distributions.append(hist)
    
    if positions and r2_distributions:
        # Create heatmap
        data = np.array(r2_distributions).T
        im = ax.imshow(data, aspect='auto', cmap='YlOrRd', 
                      extent=[min(positions), max(positions), 0, 1],
                      origin='lower', interpolation='nearest')
        
        ax.set_xlabel('Position (Mb)')
        ax.set_ylabel('R² Value')
        ax.set_title(f'R² Density Along Chromosome {chr_summary.get("chromosome", "?")}')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Variant Count', rotation=270, labelpad=15)
    else:
        ax.text(0.5, 0.5, 'No R² distribution data available', 
                ha='center', va='center', transform=ax.transAxes)


def plot_variant_density(ax, chr_summary):
    """Plot variant density along chromosome."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chunks by position
    sorted_chunks = sorted(chunk_details, key=lambda x: x.get('start_pos', 0))
    
    positions = []
    total_variants = []
    well_imputed = []
    
    for chunk in sorted_chunks:
        start_pos = chunk.get('start_pos', 0)
        end_pos = chunk.get('end_pos', start_pos + 1000000)
        midpoint = (start_pos + end_pos) / 2 / 1e6  # Mb
        
        positions.append(midpoint)
        total_variants.append(chunk.get('total_variants', 0) or 0)
        well_imputed.append(chunk.get('well_imputed', 0) or 0)
    
    # Create bar chart
    width = (max(positions) - min(positions)) / len(positions) * 0.8 if positions else 1
    
    bars1 = ax.bar(positions, total_variants, width=width, 
                   label='Total Variants', alpha=0.6, color='blue')
    bars2 = ax.bar(positions, well_imputed, width=width, 
                   label='Well-Imputed', alpha=0.8, color='green')
    
    ax.set_xlabel('Position (Mb)')
    ax.set_ylabel('Number of Variants')
    ax.set_title(f'Variant Density Along Chromosome {chr_summary.get("chromosome", "?")}')
    ax.legend()
    
    # Add annotation for sparse regions
    if total_variants:
        mean_density = np.mean(total_variants)
        sparse_threshold = mean_density * 0.5
        
        for pos, count in zip(positions, total_variants):
            if count < sparse_threshold:
                ax.annotate('Low coverage', xy=(pos, count), 
                           xytext=(pos, count + max(total_variants) * 0.1),
                           fontsize=8, ha='center',
                           arrowprops=dict(arrowstyle='->', alpha=0.5))


def plot_quality_metrics_by_position(ax, chr_summary):
    """Plot multiple quality metrics along chromosome."""
    chunk_details = chr_summary.get('chunk_details', [])
    
    if not chunk_details:
        ax.text(0.5, 0.5, 'No chunk data available', 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Sort chunks by position
    sorted_chunks = sorted(chunk_details, key=lambda x: x.get('start_pos', 0))
    
    positions = []
    r2_values = []
    concordance_values = []
    well_imputed_rates = []
    
    for chunk in sorted_chunks:
        start_pos = chunk.get('start_pos', 0)
        end_pos = chunk.get('end_pos', start_pos + 1000000)
        midpoint = (start_pos + end_pos) / 2 / 1e6  # Mb
        
        positions.append(midpoint)
        r2_values.append(chunk.get('mean_r2', 0) or 0)
        concordance_values.append(chunk.get('mean_concordance', 0) or 0)
        
        total_vars = chunk.get('total_variants', 1)
        well_imp = chunk.get('well_imputed', 0)
        well_imputed_rates.append((well_imp / total_vars) if total_vars > 0 else 0)
    
    # Plot multiple metrics
    ax2 = ax.twinx()
    
    line1 = ax.plot(positions, r2_values, 'o-', linewidth=2, markersize=6,
                   color='darkblue', label='Mean R²')
    line2 = ax.plot(positions, concordance_values, 's-', linewidth=2, markersize=6,
                   color='coral', label='Concordance')
    line3 = ax2.plot(positions, well_imputed_rates, '^-', linewidth=2, markersize=6,
                    color='green', label='Well-Imputed Rate')
    
    ax.set_xlabel('Position (Mb)')
    ax.set_ylabel('R² / Concordance', color='black')
    ax2.set_ylabel('Well-Imputed Rate', color='green')
    ax.set_title(f'Quality Metrics Along Chromosome {chr_summary.get("chromosome", "?")}')
    ax.set_ylim(0, 1.1)
    ax2.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)
    
    # Combine legends
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc='lower left')
    
    # Color y-axis labels
    ax2.tick_params(axis='y', labelcolor='green')


def create_chromosome_r2_position_plots(chr_summary, output_prefix):
    """Create comprehensive R² by position plots for chromosome."""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 12))
    
    # Add main title
    fig.suptitle(f'R² Distribution Along Chromosome {chr_summary.get("chromosome", "?")} - Dataset: {chr_summary.get("dataset", "Unknown")}', 
                 fontsize=16, fontweight='bold')
    
    # Create 2x2 grid
    ax1 = plt.subplot(2, 2, 1)
    plot_r2_along_chromosome(ax1, chr_summary)
    
    ax2 = plt.subplot(2, 2, 2)
    plot_variant_density(ax2, chr_summary)
    
    ax3 = plt.subplot(2, 2, 3)
    plot_r2_density_by_position(ax3, chr_summary)
    
    ax4 = plt.subplot(2, 2, 4)
    plot_quality_metrics_by_position(ax4, chr_summary)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    
    return fig


def main():
    parser = argparse.ArgumentParser(description='Generate chromosome R² by position plots')
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
    logger.info(f"Creating R² position plots for chromosome {args.chromosome}")
    fig = create_chromosome_r2_position_plots(chr_summary, args.output_prefix)
    
    # Save to PDF
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_r2_position.pdf"
    logger.info(f"Saving plots to {output_file}")
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
        
        # Add metadata
        d = pdf.infodict()
        d['Title'] = f'R² Distribution Along Chromosome {args.chromosome}'
        d['Author'] = 'Genotype Imputation Pipeline'
        d['Subject'] = f'R² Position Analysis for {args.dataset}'
        d['Keywords'] = f'Imputation R² Position Chromosome {args.chromosome}'
    
    plt.close(fig)
    
    logger.info(f"R² position plots saved successfully to {output_file}")
    
    # Print summary
    print(f"\nR² Position Analysis Summary:")
    print(f"  Chromosome: {args.chromosome}")
    print(f"  Dataset: {args.dataset}")
    mean_r2 = chr_summary.get('mean_r2', 0) or 0
    print(f"  Mean R²: {mean_r2:.4f}")
    print(f"  Output: {output_file}")


if __name__ == '__main__':
    main()