#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import sys
import argparse
from scipy import stats

def main():
    parser = argparse.ArgumentParser(description='Plot allele frequency comparison between imputed and reference panels')
    parser.add_argument('freq_file', help='Frequency comparison TSV file')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    parser.add_argument('--level', required=True, help='Analysis level (chunk/chromosome/genome)')
    args = parser.parse_args()
    
    # Read frequency data
    df = pd.read_csv(args.freq_file, sep='\t')
    
    # Create comprehensive comparison plots
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Scatter plot: Imputed vs Reference AF
    ax1 = fig.add_subplot(gs[0, 0])
    if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
        valid_data = df[(df['AF_imputed'] > 0) | (df['AF_reference'] > 0)]
        ax1.scatter(valid_data['AF_reference'], valid_data['AF_imputed'], 
                   alpha=0.3, s=1, color='blue')
        ax1.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='y=x')
        
        # Calculate correlation
        mask = (valid_data['AF_reference'] > 0) & (valid_data['AF_imputed'] > 0)
        if mask.sum() > 1:
            corr = np.corrcoef(valid_data.loc[mask, 'AF_reference'], 
                               valid_data.loc[mask, 'AF_imputed'])[0, 1]
            ax1.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax1.transAxes)
    
    ax1.set_xlabel('Reference Panel AF')
    ax1.set_ylabel('Imputed AF')
    ax1.set_title('Allele Frequency Correlation')
    ax1.grid(True, alpha=0.3)
    
    # 2. Histogram of AF differences
    ax2 = fig.add_subplot(gs[0, 1])
    if 'AF_diff' in df.columns:
        af_diff = df['AF_diff'].dropna()
        ax2.hist(af_diff, bins=50, edgecolor='black', alpha=0.7, color='green')
        ax2.axvline(x=0, color='red', linestyle='--', alpha=0.5)
        ax2.set_xlabel('AF Difference (Imputed - Reference)')
        ax2.set_ylabel('Count')
        ax2.set_title(f'AF Difference Distribution\nMean: {af_diff.mean():.4f}, SD: {af_diff.std():.4f}')
    ax2.grid(True, alpha=0.3)
    
    # 3. MAF comparison
    ax3 = fig.add_subplot(gs[0, 2])
    if 'MAF_imputed' in df.columns and 'MAF_reference' in df.columns:
        maf_bins = np.linspace(0, 0.5, 21)
        
        # Calculate MAF distributions
        imp_maf = df['MAF_imputed'].dropna()
        ref_maf = df['MAF_reference'].dropna()
        
        ax3.hist([ref_maf, imp_maf], bins=maf_bins, label=['Reference', 'Imputed'], 
                alpha=0.6, edgecolor='black')
        ax3.set_xlabel('Minor Allele Frequency')
        ax3.set_ylabel('Count')
        ax3.set_title('MAF Distribution Comparison')
        ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Frequency bins comparison
    ax4 = fig.add_subplot(gs[1, 0])
    if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
        # Bin frequencies for comparison
        bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
        bin_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%', '>50%']
        
        df['ref_bin'] = pd.cut(df['AF_reference'], bins=bins, labels=bin_labels)
        df['imp_bin'] = pd.cut(df['AF_imputed'], bins=bins, labels=bin_labels)
        
        # Count variants in each bin
        ref_counts = df['ref_bin'].value_counts().sort_index()
        imp_counts = df['imp_bin'].value_counts().sort_index()
        
        x = np.arange(len(bin_labels))
        width = 0.35
        
        ax4.bar(x - width/2, ref_counts, width, label='Reference', alpha=0.7)
        ax4.bar(x + width/2, imp_counts, width, label='Imputed', alpha=0.7)
        ax4.set_xlabel('Frequency Bins')
        ax4.set_ylabel('Number of Variants')
        ax4.set_title('Variant Count by Frequency Bin')
        ax4.set_xticks(x)
        ax4.set_xticklabels(bin_labels, rotation=45)
        ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. Bland-Altman plot
    ax5 = fig.add_subplot(gs[1, 1])
    if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
        valid_data = df[(df['AF_imputed'] > 0) & (df['AF_reference'] > 0)]
        if len(valid_data) > 0:
            mean_af = (valid_data['AF_imputed'] + valid_data['AF_reference']) / 2
            diff_af = valid_data['AF_imputed'] - valid_data['AF_reference']
            
            ax5.scatter(mean_af, diff_af, alpha=0.3, s=1)
            ax5.axhline(y=diff_af.mean(), color='red', linestyle='-', label=f'Mean: {diff_af.mean():.4f}')
            ax5.axhline(y=diff_af.mean() + 1.96*diff_af.std(), color='red', linestyle='--', 
                       label=f'±1.96 SD')
            ax5.axhline(y=diff_af.mean() - 1.96*diff_af.std(), color='red', linestyle='--')
            ax5.set_xlabel('Mean AF')
            ax5.set_ylabel('AF Difference (Imputed - Reference)')
            ax5.set_title('Bland-Altman Plot')
            ax5.legend(loc='upper right', fontsize=8)
    ax5.grid(True, alpha=0.3)
    
    # 6. Cumulative distribution
    ax6 = fig.add_subplot(gs[1, 2])
    if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
        # Sort and calculate cumulative distribution
        ref_sorted = np.sort(df['AF_reference'].dropna())
        imp_sorted = np.sort(df['AF_imputed'].dropna())
        
        if len(ref_sorted) > 0 and len(imp_sorted) > 0:
            ref_cdf = np.arange(1, len(ref_sorted) + 1) / len(ref_sorted)
            imp_cdf = np.arange(1, len(imp_sorted) + 1) / len(imp_sorted)
            
            ax6.plot(ref_sorted, ref_cdf, label='Reference', alpha=0.7)
            ax6.plot(imp_sorted, imp_cdf, label='Imputed', alpha=0.7)
            ax6.set_xlabel('Allele Frequency')
            ax6.set_ylabel('Cumulative Proportion')
            ax6.set_title('Cumulative AF Distribution')
            ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 7. QQ plot
    ax7 = fig.add_subplot(gs[2, 0])
    if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
        # Create QQ plot
        valid_data = df[(df['AF_imputed'] > 0) & (df['AF_reference'] > 0)]
        if len(valid_data) > 0:
            ref_quantiles = np.percentile(valid_data['AF_reference'], np.arange(0, 101, 1))
            imp_quantiles = np.percentile(valid_data['AF_imputed'], np.arange(0, 101, 1))
            
            ax7.scatter(ref_quantiles, imp_quantiles, alpha=0.6, s=10)
            ax7.plot([0, 1], [0, 1], 'r--', alpha=0.5)
            ax7.set_xlabel('Reference Quantiles')
            ax7.set_ylabel('Imputed Quantiles')
            ax7.set_title('Q-Q Plot')
    ax7.grid(True, alpha=0.3)
    
    # 8. Frequency by chromosome position (if available)
    ax8 = fig.add_subplot(gs[2, 1:])
    if 'POS' in df.columns and 'AF_diff' in df.columns:
        # Sample data if too many points
        plot_df = df.copy()
        if len(plot_df) > 10000:
            plot_df = plot_df.sample(n=10000)
        
        scatter = ax8.scatter(plot_df['POS'], plot_df['AF_diff'], 
                            c=abs(plot_df['AF_diff']), cmap='RdYlBu_r',
                            alpha=0.5, s=1)
        ax8.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax8.set_xlabel('Genomic Position')
        ax8.set_ylabel('AF Difference (Imputed - Reference)')
        ax8.set_title('AF Difference Along Chromosome')
        plt.colorbar(scatter, ax=ax8, label='|AF Difference|')
    ax8.grid(True, alpha=0.3)
    
    plt.suptitle(f'Allele Frequency Comparison - {args.prefix} ({args.level})', fontsize=16, y=1.02)
    plt.tight_layout()
    
    output_plot = f"{args.prefix}_{args.level}.freq_comparison.png"
    plt.savefig(output_plot, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Calculate and save statistics
    output_stats = f"{args.prefix}_{args.level}.freq_stats.txt"
    with open(output_stats, 'w') as f:
        f.write(f"Frequency Comparison Statistics - {args.prefix} ({args.level})\n")
        f.write("=" * 60 + "\n\n")
        
        if 'AF_imputed' in df.columns and 'AF_reference' in df.columns:
            # Overall statistics
            f.write("Overall Statistics:\n")
            f.write(f"Total variants: {len(df)}\n")
            f.write(f"Imputed-only variants: {(df['AF_reference'] == 0).sum()}\n")
            f.write(f"Reference-only variants: {(df['AF_imputed'] == 0).sum()}\n")
            shared = ((df['AF_imputed'] > 0) & (df['AF_reference'] > 0)).sum()
            f.write(f"Shared variants: {shared}\n\n")
            
            # Correlation statistics
            valid_data = df[(df['AF_imputed'] > 0) & (df['AF_reference'] > 0)]
            if len(valid_data) > 1:
                corr = np.corrcoef(valid_data['AF_reference'], valid_data['AF_imputed'])[0, 1]
                f.write(f"Correlation (shared variants): {corr:.4f}\n")
                
                # Perform statistical tests
                if len(valid_data) > 2:
                    ks_stat, ks_pval = stats.ks_2samp(valid_data['AF_reference'], 
                                                      valid_data['AF_imputed'])
                    f.write(f"Kolmogorov-Smirnov test: statistic={ks_stat:.4f}, p-value={ks_pval:.4e}\n")
            
            # AF difference statistics
            if 'AF_diff' in df.columns:
                af_diff = df['AF_diff'].dropna()
                f.write(f"\nAF Difference Statistics:\n")
                f.write(f"Mean difference: {af_diff.mean():.6f}\n")
                f.write(f"Median difference: {af_diff.median():.6f}\n")
                f.write(f"Std deviation: {af_diff.std():.6f}\n")
                f.write(f"Min difference: {af_diff.min():.6f}\n")
                f.write(f"Max difference: {af_diff.max():.6f}\n")
                f.write(f"Variants with |diff| > 0.1: {(abs(af_diff) > 0.1).sum()}\n")
                f.write(f"Variants with |diff| > 0.2: {(abs(af_diff) > 0.2).sum()}\n")
    
    print(f"Frequency comparison plot saved: {output_plot}")
    print(f"Statistics saved: {output_stats}")
    print(f"Analyzed {len(df)} variants")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"H3ABIONET_CHIPIMPUTATION:CHIPIMPUTATION:REPORT:PLOT_FREQ_COMPARISON":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\n')
        f.write(f'    pandas: {pd.__version__}\n')
        f.write(f'    numpy: {np.__version__}\n')
        
        # Handle scipy version extraction properly
        scipy_version = "1.9.0"  # default
        if hasattr(stats, "__version__"):
            scipy_version = stats.__version__
        elif hasattr(stats, "_version"):
            if hasattr(stats._version, "__version__"):
                scipy_version = stats._version.__version__
        else:
            # Try importing scipy directly
            try:
                import scipy
                scipy_version = scipy.__version__
            except:
                pass
        
        f.write(f'    scipy: {scipy_version}\n')

if __name__ == "__main__":
    main()