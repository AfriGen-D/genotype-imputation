#!/usr/bin/env python3
"""
Generate genome-wide performance plots for imputation quality assessment.
"""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO)
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

def create_genome_performance_plot(genome_summary):
    """Create comprehensive genome-wide performance visualization."""
    fig = plt.figure(figsize=(16, 12))
    
    # Extract chromosome data
    chr_details = genome_summary.get('chromosome_details', [])
    chr_details_sorted = sorted(chr_details, key=lambda x: chromosome_sort_key(x.get('chromosome', '')))
    
    chromosomes = [chr_data.get('chromosome', '').replace('chr', '') for chr_data in chr_details_sorted]
    mean_r2 = [chr_data.get('mean_r2', 0) or 0 for chr_data in chr_details_sorted]
    total_variants = [chr_data.get('variants', 0) or 0 for chr_data in chr_details_sorted]
    well_imputed = [chr_data.get('well_imputed', 0) or 0 for chr_data in chr_details_sorted]
    
    # Calculate well-imputed rates
    well_imputed_rates = [(w/t)*100 if t > 0 else 0 for w, t in zip(well_imputed, total_variants)]
    
    # Create 2x2 subplot layout
    ax1 = plt.subplot(2, 2, 1)
    bars1 = ax1.bar(chromosomes, mean_r2, alpha=0.8, color='steelblue')
    ax1.set_title('Mean R² by Chromosome', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Mean R²')
    ax1.tick_params(axis='x', rotation=45)
    ax1.axhline(y=genome_summary.get('mean_r2', 0), color='red', linestyle='--', alpha=0.7)
    
    ax2 = plt.subplot(2, 2, 2)
    bars2 = ax2.bar(chromosomes, [v/1e6 for v in total_variants], alpha=0.8, color='forestgreen')
    ax2.set_title('Total Variants by Chromosome (Millions)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Variants (M)')
    ax2.tick_params(axis='x', rotation=45)
    
    ax3 = plt.subplot(2, 2, 3)
    bars3 = ax3.bar(chromosomes, well_imputed_rates, alpha=0.8, color='coral')
    ax3.set_title('Well-Imputed Rate by Chromosome (%)', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Chromosome')
    ax3.set_ylabel('Well-Imputed Rate (%)')
    ax3.tick_params(axis='x', rotation=45)
    overall_rate = genome_summary.get('well_imputed_rate', 0) * 100 if genome_summary.get('well_imputed_rate') else 0
    ax3.axhline(y=overall_rate, color='red', linestyle='--', alpha=0.7)
    
    ax4 = plt.subplot(2, 2, 4)
    ax4.axis('off')
    
    # Summary statistics
    stats_text = f"""
    Genome-Wide Performance Summary
    
    Dataset: {genome_summary.get('dataset', 'Unknown')}
    Reference: {genome_summary.get('ref_name', 'Unknown')}
    
    Total Chromosomes: {genome_summary.get('chromosomes_processed', 0)}
    Total Chunks: {genome_summary.get('total_chunks', 0)}
    Total Variants: {genome_summary.get('total_variants', 0):,}
    Well-Imputed: {genome_summary.get('well_imputed_variants', 0):,}
    
    Overall Mean R²: {genome_summary.get('mean_r2', 0):.4f}
    Well-Imputed Rate: {overall_rate:.1f}%
    """
    
    ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.suptitle(f'Genome-Wide Imputation Performance - {genome_summary.get("dataset", "Unknown")}', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    
    return fig

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-summary', help='Genome summary JSON')
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--ref-name', required=True)
    parser.add_argument('--dataset', required=True)
    # Add all possible arguments to avoid errors
    parser.add_argument('--performance-plot', help='Performance plot')
    parser.add_argument('--r2-plot', help='R2 plot')
    parser.add_argument('--accuracy-plot', help='Accuracy plot')
    parser.add_argument('--maf-plot', help='MAF plot')
    parser.add_argument('--freq-plot', help='Freq plot')
    parser.add_argument('--chr-plot', help='Chr plot')
    
    args = parser.parse_args()
    
    # Determine output file based on script name
    script_name = Path(__file__).stem
    if 'performance' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_performance.pdf"
    elif 'r2_distribution' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_r2_distribution.pdf"
    elif 'accuracy' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_accuracy.pdf"
    elif 'maf_analysis' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_maf_analysis.pdf"
    elif 'freq_comparison' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_freq_comparison.pdf"
    elif 'chr_comparison' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_chr_comparison.pdf"
    elif 'combine_genome' in script_name:
        output_file = f"{args.output_prefix}_{args.ref_name}.genome_all_metrics.pdf"
    else:
        output_file = f"{args.output_prefix}_{args.ref_name}.plot.pdf"
    
    # Load data and create plot
    if args.genome_summary:
        genome_summary = load_genome_summary(args.genome_summary)
        fig = create_genome_performance_plot(genome_summary)
    else:
        # Fallback to placeholder
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f'{script_name}\nNo data provided\n{args.dataset}', 
                ha='center', va='center', fontsize=16, transform=ax.transAxes)
        ax.set_title(f'{script_name.replace("_", " ").title()}')
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    print(f"Generated {output_file}")

if __name__ == '__main__':
    main()
