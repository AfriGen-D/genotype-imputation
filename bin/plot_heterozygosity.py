#!/usr/bin/env python3
"""
Heterozygosity Plot
Analyzes sample-level heterozygosity for quality control
Detects contamination, inbreeding, or technical issues
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
from pathlib import Path
from scipy import stats

def calculate_sample_heterozygosity(vcf_file):
    """Calculate heterozygosity for each sample in VCF"""
    sample_stats = {}
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    sample_names = header[9:] if len(header) > 9 else []
                    # Initialize stats for each sample
                    for sample in sample_names:
                        sample_stats[sample] = {
                            'het_count': 0,
                            'hom_ref_count': 0,
                            'hom_alt_count': 0,
                            'missing_count': 0,
                            'total_sites': 0
                        }
                continue
            
            fields = line.strip().split('\t')
            format_field = fields[8]
            
            # Check for GT field
            if 'GT' not in format_field:
                continue
            
            format_keys = format_field.split(':')
            gt_idx = format_keys.index('GT')
            
            # Process each sample
            for i, sample_field in enumerate(fields[9:]):
                if i >= len(sample_names):
                    break
                
                sample_name = sample_names[i]
                values = sample_field.split(':')
                
                if len(values) > gt_idx:
                    gt = values[gt_idx]
                    sample_stats[sample_name]['total_sites'] += 1
                    
                    # Parse genotype
                    if gt in ['.', './.', '.|.']:
                        sample_stats[sample_name]['missing_count'] += 1
                    elif '/' in gt or '|' in gt:
                        alleles = gt.replace('|', '/').split('/')
                        if len(alleles) == 2:
                            try:
                                a1, a2 = int(alleles[0]), int(alleles[1])
                                if a1 != a2:
                                    sample_stats[sample_name]['het_count'] += 1
                                elif a1 == 0:
                                    sample_stats[sample_name]['hom_ref_count'] += 1
                                else:
                                    sample_stats[sample_name]['hom_alt_count'] += 1
                            except (ValueError, TypeError):
                                sample_stats[sample_name]['missing_count'] += 1
    
    # Calculate heterozygosity rates
    results = []
    for sample, stats in sample_stats.items():
        n_called = stats['total_sites'] - stats['missing_count']
        
        if n_called > 0:
            het_rate = stats['het_count'] / n_called
            call_rate = n_called / stats['total_sites'] if stats['total_sites'] > 0 else 0
            
            results.append({
                'sample': sample,
                'het_rate': het_rate,
                'het_count': stats['het_count'],
                'hom_ref_count': stats['hom_ref_count'],
                'hom_alt_count': stats['hom_alt_count'],
                'missing_count': stats['missing_count'],
                'total_sites': stats['total_sites'],
                'call_rate': call_rate,
                'n_called': n_called
            })
    
    return pd.DataFrame(results)

def identify_outliers(df, column='het_rate', n_std=3):
    """Identify outlier samples based on heterozygosity"""
    mean = df[column].mean()
    std = df[column].std()
    
    df['z_score'] = (df[column] - mean) / std
    df['is_outlier'] = np.abs(df['z_score']) > n_std
    
    return df, mean, std

def calculate_f_statistic(df):
    """Calculate inbreeding coefficient (F) for each sample"""
    # F = 1 - (observed het / expected het)
    # For simplicity, using population average as expected
    expected_het = df['het_rate'].mean()
    
    f_stats = []
    for _, row in df.iterrows():
        if expected_het > 0:
            f = 1 - (row['het_rate'] / expected_het)
        else:
            f = 0
        f_stats.append(f)
    
    df['f_statistic'] = f_stats
    return df

def main():
    parser = argparse.ArgumentParser(description='Plot sample heterozygosity for QC')
    parser.add_argument('vcf_file', help='Input VCF file')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chr', help='Chromosome', default='')
    parser.add_argument('--n-std', type=float, default=3, help='Number of std deviations for outlier detection')
    
    args = parser.parse_args()
    
    print(f"Calculating heterozygosity from {args.vcf_file}...")
    het_df = calculate_sample_heterozygosity(args.vcf_file)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    if len(het_df) > 0:
        # Identify outliers
        het_df, mean_het, std_het = identify_outliers(het_df, 'het_rate', n_std=args.n_std)
        het_df = calculate_f_statistic(het_df)
        
        # 1. Heterozygosity distribution
        ax = axes[0, 0]
        
        # Histogram
        n, bins, patches = ax.hist(het_df['het_rate'], bins=30, edgecolor='black', 
                                   alpha=0.7, color='steelblue')
        
        # Color outliers
        for i, patch in enumerate(patches):
            bin_center = (bins[i] + bins[i+1]) / 2
            if abs(bin_center - mean_het) > args.n_std * std_het:
                patch.set_facecolor('red')
                patch.set_alpha(0.7)
        
        # Add mean and threshold lines
        ax.axvline(x=mean_het, color='green', linestyle='-', 
                  label=f'Mean = {mean_het:.4f}', linewidth=2)
        ax.axvline(x=mean_het - args.n_std * std_het, color='red', linestyle='--',
                  label=f'±{args.n_std}σ threshold', alpha=0.5)
        ax.axvline(x=mean_het + args.n_std * std_het, color='red', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('Heterozygosity Rate')
        ax.set_ylabel('Number of Samples')
        ax.set_title('Sample Heterozygosity Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Heterozygosity vs Call Rate
        ax = axes[0, 1]
        
        colors = ['red' if outlier else 'blue' for outlier in het_df['is_outlier']]
        ax.scatter(het_df['call_rate'], het_df['het_rate'], 
                  c=colors, alpha=0.6, edgecolors='black')
        
        # Add reference lines
        ax.axhline(y=mean_het, color='green', linestyle='--', alpha=0.3)
        ax.axvline(x=0.95, color='orange', linestyle='--', alpha=0.3,
                  label='95% call rate threshold')
        
        ax.set_xlabel('Call Rate')
        ax.set_ylabel('Heterozygosity Rate')
        ax.set_title('Heterozygosity vs Call Rate')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Highlight problematic samples
        problem_samples = het_df[(het_df['is_outlier']) | (het_df['call_rate'] < 0.95)]
        if len(problem_samples) > 0:
            for _, sample in problem_samples.iterrows():
                ax.annotate(sample['sample'][:10], 
                           xy=(sample['call_rate'], sample['het_rate']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.7)
        
        # 3. Inbreeding coefficient (F) distribution
        ax = axes[1, 0]
        
        ax.hist(het_df['f_statistic'], bins=30, edgecolor='black',
               alpha=0.7, color='coral')
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5,
                  label='F=0 (Random mating)')
        
        # Add interpretation zones
        ax.axvspan(-1, -0.05, alpha=0.2, color='blue', label='Excess heterozygosity')
        ax.axvspan(0.05, 1, alpha=0.2, color='red', label='Inbreeding')
        
        ax.set_xlabel('Inbreeding Coefficient (F)')
        ax.set_ylabel('Number of Samples')
        ax.set_title('Inbreeding Coefficient Distribution')
        ax.set_xlim(-0.5, 0.5)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # 4. Summary statistics and QC recommendations
        ax = axes[1, 1]
        ax.axis('off')
        
        n_outliers = het_df['is_outlier'].sum()
        n_low_call = (het_df['call_rate'] < 0.95).sum()
        
        stats_text = "Heterozygosity QC Summary\n" + "="*35 + "\n\n"
        stats_text += f"Total samples: {len(het_df)}\n"
        stats_text += f"Mean heterozygosity: {mean_het:.4f}\n"
        stats_text += f"Std deviation: {std_het:.4f}\n"
        stats_text += f"Mean call rate: {het_df['call_rate'].mean():.4f}\n\n"
        
        stats_text += "QC Flags:\n"
        stats_text += f"  Heterozygosity outliers: {n_outliers}\n"
        stats_text += f"  Low call rate (<95%): {n_low_call}\n"
        stats_text += f"  Total flagged: {len(problem_samples)}\n\n"
        
        stats_text += "Outlier Samples:\n"
        if n_outliers > 0:
            outlier_samples = het_df[het_df['is_outlier']]
            for _, sample in outlier_samples.head(5).iterrows():
                stats_text += f"  {sample['sample'][:20]}: Het={sample['het_rate']:.3f}"
                stats_text += f" (z={sample['z_score']:.1f})\n"
            if n_outliers > 5:
                stats_text += f"  ... and {n_outliers - 5} more\n"
        else:
            stats_text += "  None detected\n"
        
        stats_text += "\nInterpretation:\n"
        if n_outliers == 0:
            stats_text += "✓ No heterozygosity outliers detected\n"
        else:
            stats_text += "⚠ Review flagged samples for:\n"
            stats_text += "  • Contamination (high het)\n"
            stats_text += "  • Inbreeding (low het)\n"
            stats_text += "  • Technical issues\n"
    else:
        # No data
        for ax in axes.flat:
            ax.text(0.5, 0.5, 'No heterozygosity data available',
                   ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
        
        stats_text = "No samples found in VCF file"
    
    if len(het_df) > 0 and 'stats_text' in locals():
        axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                       fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    # Add main title
    chr_info = f" - Chr {args.chr}" if args.chr else " - Genome-wide"
    fig.suptitle(f'Sample Heterozygosity Analysis\nDataset: {args.sample_id} | Reference: {args.ref_name}{chr_info}',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(args.output_pdf, dpi=150, bbox_inches='tight')
    print(f"Heterozygosity plot saved to {args.output_pdf}")
    
    # Also save outlier list if any
    if len(het_df) > 0 and n_outliers > 0:
        outlier_file = args.output_pdf.replace('.pdf', '_outliers.txt')
        outlier_samples[['sample', 'het_rate', 'z_score']].to_csv(
            outlier_file, sep='\t', index=False
        )
        print(f"Outlier list saved to {outlier_file}")

if __name__ == "__main__":
    main()