#!/usr/bin/env python3
"""
Generate chromosome-level accuracy plots for imputation quality metrics.
"""

import json
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def load_chr_summary(summary_path):
    """Load chromosome summary JSON."""
    with open(summary_path, 'r') as f:
        return json.load(f)

def plot_accuracy_metrics(chr_summary, output_file, ref_name, dataset, chromosome):
    """Create accuracy plots for chromosome."""
    with PdfPages(output_file) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Imputation Accuracy - {dataset} {chromosome} ({ref_name})', fontsize=16)
        
        # Plot 1: Chunk accuracy distribution
        ax = axes[0, 0]
        if chr_summary.get('chunk_details'):
            chunks = chr_summary['chunk_details']
            accuracies = [c.get('accuracy', 0) for c in chunks if c.get('accuracy') is not None]
            if accuracies:
                ax.hist(accuracies, bins=20, edgecolor='black', alpha=0.7)
                ax.axvline(sum(accuracies)/len(accuracies), color='red', linestyle='--', label=f'Mean: {sum(accuracies)/len(accuracies):.3f}')
                ax.legend()
        ax.set_xlabel('Accuracy')
        ax.set_ylabel('Number of Chunks')
        ax.set_title('Accuracy Distribution Across Chunks')
        
        # Plot 2: Concordance rates
        ax = axes[0, 1]
        if chr_summary.get('concordance_rate'):
            conc_rate = chr_summary['concordance_rate']
            ax.bar(['Concordance Rate'], [conc_rate], color='green', alpha=0.7)
            ax.set_ylim([0, 1])
            ax.set_ylabel('Rate')
            ax.set_title(f'Overall Concordance Rate: {conc_rate:.3f}')
        else:
            ax.text(0.5, 0.5, 'No concordance data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Concordance Rate')
        
        # Plot 3: Genotype accuracy by MAF bins
        ax = axes[1, 0]
        ax.text(0.5, 0.5, 'MAF-stratified accuracy\n(to be implemented)', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Accuracy by MAF Bins')
        
        # Plot 4: Summary statistics
        ax = axes[1, 1]
        ax.axis('off')
        summary_text = f"""
        Chromosome: {chromosome}
        Reference Panel: {ref_name}
        Total Variants: {chr_summary.get('total_variants', 'N/A')}
        Well Imputed: {chr_summary.get('well_imputed_variants', 'N/A')}
        Mean Accuracy: {chr_summary.get('mean_accuracy', 'N/A')}
        """
        ax.text(0.1, 0.5, summary_text, fontsize=12, verticalalignment='center')
        ax.set_title('Summary Statistics')
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate chromosome-level accuracy plots')
    parser.add_argument('--chr-summary', required=True, help='Chromosome summary JSON file')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--chromosome', required=True, help='Chromosome')
    
    args = parser.parse_args()
    
    # Load chromosome summary
    chr_summary = load_chr_summary(args.chr_summary)
    
    # Generate output filename
    output_file = f"{args.output_prefix}_{args.ref_name}.chr_accuracy.pdf"
    
    # Create plots
    plot_accuracy_metrics(chr_summary, output_file, args.ref_name, args.dataset, args.chromosome)
    
    print(f"Accuracy plots saved to {output_file}")

if __name__ == '__main__':
    main()