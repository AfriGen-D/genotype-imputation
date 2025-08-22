#!/usr/bin/env python3
"""
Generate genome-wide chromosome comparison plots.
Similar to r2_snppos analysis but comparing across chromosomes.
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

def create_genome_chr_comparison_plot(genome_summary):
    """Create genome-wide chromosome comparison visualization."""
    fig = plt.figure(figsize=(16, 12))
    
    # Extract chromosome data
    chr_details = genome_summary.get('chromosome_details', [])
    chr_details_sorted = sorted(chr_details, key=lambda x: chromosome_sort_key(x.get('chromosome', '')))
    
    chromosomes = [chr_data.get('chromosome', '').replace('chr', '') for chr_data in chr_details_sorted]
    mean_r2 = [chr_data.get('mean_r2', 0) or 0 for chr_data in chr_details_sorted]
    total_variants = [chr_data.get('variants', 0) or 0 for chr_data in chr_details_sorted]
    well_imputed = [chr_data.get('well_imputed', 0) or 0 for chr_data in chr_details_sorted]
    chunks = [chr_data.get('chunks', 0) or 0 for chr_data in chr_details_sorted]
    
    # Calculate metrics
    well_imputed_rates = [(w/t)*100 if t > 0 else 0 for w, t in zip(well_imputed, total_variants)]
    variants_per_chunk = [(t/c) if c > 0 else 0 for t, c in zip(total_variants, chunks)]
    
    # Create 2x2 subplot layout
    
    # Subplot 1: R² vs chromosome size scatter plot
    ax1 = plt.subplot(2, 2, 1)
    
    # Color by chromosome for easier identification
    colors = plt.cm.viridis(np.linspace(0, 1, len(chromosomes)))
    scatter = ax1.scatter([v/1e6 for v in total_variants], mean_r2, 
                         c=range(len(chromosomes)), cmap='viridis', 
                         s=100, alpha=0.7, edgecolors='black')
    
    # Add chromosome labels
    for i, (chr_name, x, y) in enumerate(zip(chromosomes, [v/1e6 for v in total_variants], mean_r2)):
        ax1.annotate(chr_name, (x, y), xytext=(5, 5), textcoords='offset points', 
                    fontsize=8, alpha=0.8)
    
    ax1.set_title('R² vs Chromosome Size', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Total Variants (Millions)')
    ax1.set_ylabel('Mean R²')
    ax1.grid(True, alpha=0.3)
    
    # Add trend line
    from scipy import stats
    if len(total_variants) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            [v/1e6 for v in total_variants], mean_r2)
        line_x = np.array([min(v/1e6 for v in total_variants), 
                          max(v/1e6 for v in total_variants)])
        line_y = slope * line_x + intercept
        ax1.plot(line_x, line_y, 'r--', alpha=0.8, 
                label=f'Trend (R²={r_value:.3f})')
        ax1.legend()
    
    # Subplot 2: Chromosome metrics comparison
    ax2 = plt.subplot(2, 2, 2)
    
    x = np.arange(len(chromosomes))
    width = 0.35
    
    # Normalize metrics for comparison (0-1 scale)
    norm_r2 = np.array(mean_r2)
    norm_well_imp = np.array(well_imputed_rates) / 100  # Convert percentage to 0-1
    
    bars1 = ax2.bar(x - width/2, norm_r2, width, label='Mean R²', alpha=0.8)
    bars2 = ax2.bar(x + width/2, norm_well_imp, width, label='Well-Imputed Rate', alpha=0.8)
    
    ax2.set_title('Normalized Metrics by Chromosome', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Normalized Value (0-1)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(chromosomes, rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Chromosome heterogeneity analysis
    ax3 = plt.subplot(2, 2, 3)
    
    # Calculate coefficient of variation for each chromosome (using info stats)
    cv_values = []
    for chr_data in chr_details_sorted:
        info_stats = chr_data.get('info_stats', {})
        mean_val = info_stats.get('mean', 0) or 0
        # Estimate std from quartiles (rough approximation)
        q25 = info_stats.get('q25', 0) or 0
        q75 = info_stats.get('q75', 0) or 0
        std_approx = (q75 - q25) / 1.35  # Approximate std from IQR
        cv = (std_approx / mean_val) if mean_val > 0 else 0
        cv_values.append(cv)
    
    bars3 = ax3.bar(chromosomes, cv_values, alpha=0.8, color='coral')
    ax3.set_title('R² Heterogeneity by Chromosome (CV)', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Chromosome')
    ax3.set_ylabel('Coefficient of Variation')
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(True, alpha=0.3)
    
    # Highlight chromosomes with highest variability
    max_cv_idx = np.argmax(cv_values)
    bars3[max_cv_idx].set_color('red')
    ax3.annotate(f'Most variable\n({chromosomes[max_cv_idx]})', 
                xy=(max_cv_idx, cv_values[max_cv_idx]),
                xytext=(10, 10), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    # Subplot 4: Quality summary by chromosome group
    ax4 = plt.subplot(2, 2, 4)
    
    # Group chromosomes by quality
    high_quality_chrs = [chr_name for chr_name, r2 in zip(chromosomes, mean_r2) if r2 >= 0.45]
    medium_quality_chrs = [chr_name for chr_name, r2 in zip(chromosomes, mean_r2) if 0.40 <= r2 < 0.45]
    low_quality_chrs = [chr_name for chr_name, r2 in zip(chromosomes, mean_r2) if r2 < 0.40]
    
    # Create pie chart
    sizes = [len(high_quality_chrs), len(medium_quality_chrs), len(low_quality_chrs)]
    labels = [f'High (R²≥0.45)\n{len(high_quality_chrs)} chr', 
             f'Medium (0.40≤R²<0.45)\n{len(medium_quality_chrs)} chr',
             f'Lower (R²<0.40)\n{len(low_quality_chrs)} chr']
    colors = ['lightgreen', 'orange', 'lightcoral']
    
    # Only show non-zero categories
    non_zero_sizes = []
    non_zero_labels = []
    non_zero_colors = []
    
    for size, label, color in zip(sizes, labels, colors):
        if size > 0:
            non_zero_sizes.append(size)
            non_zero_labels.append(label)
            non_zero_colors.append(color)
    
    if non_zero_sizes:
        wedges, texts, autotexts = ax4.pie(non_zero_sizes, labels=non_zero_labels, 
                                          colors=non_zero_colors, autopct='%1.1f%%',
                                          startangle=90)
        ax4.set_title('Chromosome Quality Distribution', fontsize=14, fontweight='bold')
        
        # Add chromosome names to legend
        legend_text = []
        if high_quality_chrs:
            legend_text.append(f"High: {', '.join(high_quality_chrs)}")
        if medium_quality_chrs:
            legend_text.append(f"Medium: {', '.join(medium_quality_chrs)}")
        if low_quality_chrs:
            legend_text.append(f"Lower: {', '.join(low_quality_chrs)}")
        
        # Add text box with details
        details_text = '\n'.join(legend_text)
        ax4.text(1.3, 0.5, details_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='center', 
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    else:
        ax4.text(0.5, 0.5, 'No chromosome data available', ha='center', va='center',
                transform=ax4.transAxes, fontsize=14)
    
    plt.suptitle(f'Genome-Wide Chromosome Comparison - {genome_summary.get("dataset", "Unknown")}', 
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
        fig = create_genome_chr_comparison_plot(genome_summary)
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
