#!/usr/bin/env python3
"""
Generate genome-wide MAF (Minor Allele Frequency) analysis plots.
"""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def load_genome_summary(summary_path):
    """Load genome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)

def create_genome_maf_plot(genome_summary):
    """Create genome-wide MAF analysis visualization."""
    fig = plt.figure(figsize=(14, 10))
    
    # MAF distribution across genome
    maf_bins = genome_summary.get('maf_bins', {})
    mean_r2_by_maf = genome_summary.get('mean_r2_by_maf', {})
    
    if maf_bins and mean_r2_by_maf:
        # Subplot 1: MAF distribution
        ax1 = plt.subplot(2, 2, 1)
        bin_names = sorted(maf_bins.keys())
        counts = [maf_bins[bin_name] for bin_name in bin_names]
        total = sum(counts)
        percentages = [(c/total)*100 if total > 0 else 0 for c in counts]
        
        bars = ax1.bar(bin_names, percentages, alpha=0.8, color='teal')
        ax1.set_title('Genome-Wide MAF Distribution', fontsize=14, fontweight='bold')
        ax1.set_xlabel('MAF Bin')
        ax1.set_ylabel('Percentage of Variants (%)')
        ax1.tick_params(axis='x', rotation=45)
        
        # Add count labels
        for bar, count, pct in zip(bars, counts, percentages):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=9)
        
        # Subplot 2: R² by MAF
        ax2 = plt.subplot(2, 2, 2)
        maf_bin_names = sorted(mean_r2_by_maf.keys())
        r2_means = [mean_r2_by_maf[bin_name].get('mean', 0) for bin_name in maf_bin_names]
        n_variants = [mean_r2_by_maf[bin_name].get('n_variants', 0) for bin_name in maf_bin_names]
        
        x = np.arange(len(maf_bin_names))
        bars2 = ax2.bar(x, r2_means, alpha=0.8, color='coral')
        ax2.set_title('Mean R² by MAF Category', fontsize=14, fontweight='bold')
        ax2.set_xlabel('MAF Bin')
        ax2.set_ylabel('Mean R²')
        ax2.set_xticks(x)
        ax2.set_xticklabels(maf_bin_names, rotation=45)
        
        # Add sample size annotations
        for i, (r2, n) in enumerate(zip(r2_means, n_variants)):
            ax2.text(i, r2 + 0.02, f'n={n}', ha='center', va='bottom', fontsize=8)
        
        # Subplot 3: Log-scale MAF distribution
        ax3 = plt.subplot(2, 2, 3)
        ax3.bar(bin_names, counts, alpha=0.8, color='forestgreen')
        # Only set log scale if we have positive values
        if any(c > 0 for c in counts):
            ax3.set_yscale('log')
            ax3.set_ylabel('Number of Variants (log)')
            ax3.set_title('MAF Distribution (Log Scale)', fontsize=14, fontweight='bold')
        else:
            ax3.set_ylabel('Number of Variants')
            ax3.set_title('MAF Distribution', fontsize=14, fontweight='bold')
        ax3.set_xlabel('MAF Bin')
        ax3.tick_params(axis='x', rotation=45)
        
        # Subplot 4: Summary statistics
        ax4 = plt.subplot(2, 2, 4)
        ax4.axis('off')
        
        total_vars = genome_summary.get('total_variants', 0)
        rare_variants = maf_bins.get('0.00-0.01', 0) + maf_bins.get('0.01-0.05', 0)
        common_variants = total_vars - rare_variants if total_vars > rare_variants else 0
        
        # Calculate percentages safely
        rare_pct = (rare_variants/total_vars)*100 if total_vars > 0 else 0
        common_pct = (common_variants/total_vars)*100 if total_vars > 0 else 0
        best_r2_cat = max(maf_bin_names, key=lambda x: mean_r2_by_maf[x].get('mean', 0)) if maf_bin_names else 'N/A'
        max_r2 = max(r2_means) if r2_means else 0
        
        stats_text = f"""
        MAF Analysis Summary
        
        Total Variants: {total_vars:,}
        
        Rare Variants (MAF < 0.05): {rare_variants:,}
        ({rare_pct:.1f}%)
        
        Common Variants (MAF ≥ 0.05): {common_variants:,}
        ({common_pct:.1f}%)
        
        Best R² Category: {best_r2_cat}
        (R² = {max_r2:.4f})
        """
        
        ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    else:
        ax1 = plt.subplot(1, 1, 1)
        ax1.text(0.5, 0.5, 'No MAF data available', ha='center', va='center', 
                transform=ax1.transAxes, fontsize=16)
        ax1.set_title('MAF Analysis - No Data')
    
    plt.suptitle(f'Genome-Wide MAF Analysis - {genome_summary.get("dataset", "Unknown")}', 
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
        fig = create_genome_maf_plot(genome_summary)
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
