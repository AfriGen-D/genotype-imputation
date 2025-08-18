#!/usr/bin/env python3
"""
Concordance vs MAF Plot
Compares imputed genotypes with true genotypes across MAF bins
Used when validation data or masked variants are available
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
from pathlib import Path

def read_vcf_genotypes(vcf_file, sample_col=9):
    """Read genotypes from VCF file"""
    genotypes = []
    positions = []
    refs = []
    alts = []
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt = fields[:5]
            format_field = fields[8]
            sample_field = fields[sample_col] if len(fields) > sample_col else None
            
            if not sample_field:
                continue
            
            # Extract GT field
            format_keys = format_field.split(':')
            if 'GT' not in format_keys:
                continue
            
            gt_idx = format_keys.index('GT')
            sample_values = sample_field.split(':')
            
            if len(sample_values) > gt_idx:
                gt = sample_values[gt_idx]
                genotypes.append(gt)
                positions.append(f"{chrom}:{pos}")
                refs.append(ref)
                alts.append(alt)
    
    return pd.DataFrame({
        'position': positions,
        'genotype': genotypes,
        'ref': refs,
        'alt': alts
    })

def calculate_maf(genotypes_df):
    """Calculate minor allele frequency for each variant"""
    mafs = []
    
    for _, row in genotypes_df.iterrows():
        gt = row['genotype']
        
        # Parse genotype (handle phased and unphased)
        if '|' in gt:
            alleles = gt.split('|')
        elif '/' in gt:
            alleles = gt.split('/')
        else:
            mafs.append(np.nan)
            continue
        
        # Count alleles (assuming biallelic)
        try:
            allele_counts = [int(a) if a != '.' else np.nan for a in alleles]
            if any(pd.isna(allele_counts)):
                mafs.append(np.nan)
            else:
                alt_freq = sum(allele_counts) / (2.0 * len(allele_counts))
                maf = min(alt_freq, 1 - alt_freq)
                mafs.append(maf)
        except:
            mafs.append(np.nan)
    
    genotypes_df['maf'] = mafs
    return genotypes_df

def compare_genotypes(true_df, imputed_df):
    """Compare true and imputed genotypes"""
    # Merge on position
    merged = pd.merge(true_df, imputed_df, on='position', suffixes=('_true', '_imputed'))
    
    # Calculate concordance
    concordance = []
    for _, row in merged.iterrows():
        gt_true = row['genotype_true']
        gt_imputed = row['genotype_imputed']
        
        # Simple concordance check
        if gt_true == gt_imputed:
            concordance.append(1)
        else:
            concordance.append(0)
    
    merged['concordant'] = concordance
    return merged

def bin_by_maf(df, n_bins=10):
    """Bin variants by MAF and calculate concordance per bin"""
    # Filter valid MAF values
    valid_df = df[df['maf_true'].notna()].copy()
    
    if len(valid_df) == 0:
        return None
    
    # Create MAF bins
    maf_bins = np.linspace(0, 0.5, n_bins + 1)
    valid_df['maf_bin'] = pd.cut(valid_df['maf_true'], bins=maf_bins, include_lowest=True)
    
    # Calculate concordance per bin
    concordance_data = []
    for bin_name in valid_df['maf_bin'].unique():
        if pd.isna(bin_name):
            continue
        
        bin_data = valid_df[valid_df['maf_bin'] == bin_name]
        if len(bin_data) == 0:
            continue
        
        concordance_rate = bin_data['concordant'].mean()
        mean_maf = bin_data['maf_true'].mean()
        
        concordance_data.append({
            'maf_bin': str(bin_name),
            'mean_maf': mean_maf,
            'concordance_rate': concordance_rate,
            'n_variants': len(bin_data),
            'n_concordant': bin_data['concordant'].sum()
        })
    
    return pd.DataFrame(concordance_data)

def main():
    parser = argparse.ArgumentParser(description='Plot concordance vs MAF for validation')
    parser.add_argument('imputed_vcf', help='Imputed VCF file')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--true-vcf', help='True/validation VCF file (if available)')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chr', help='Chromosome', default='')
    parser.add_argument('--n-bins', type=int, default=10, help='Number of MAF bins')
    
    args = parser.parse_args()
    
    print(f"Reading imputed VCF: {args.imputed_vcf}")
    imputed_df = read_vcf_genotypes(args.imputed_vcf)
    imputed_df = calculate_maf(imputed_df)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    if args.true_vcf and Path(args.true_vcf).exists():
        print(f"Reading validation VCF: {args.true_vcf}")
        true_df = read_vcf_genotypes(args.true_vcf)
        true_df = calculate_maf(true_df)
        
        # Compare genotypes
        comparison_df = compare_genotypes(true_df, imputed_df)
        
        # Bin by MAF
        concordance_by_maf = bin_by_maf(comparison_df, n_bins=args.n_bins)
        
        if concordance_by_maf is not None and len(concordance_by_maf) > 0:
            # 1. Main concordance vs MAF plot
            ax = axes[0, 0]
            ax.plot(concordance_by_maf['mean_maf'], concordance_by_maf['concordance_rate'],
                   'o-', linewidth=2, markersize=8, color='steelblue')
            ax.set_xlabel('Minor Allele Frequency (MAF)')
            ax.set_ylabel('Concordance Rate')
            ax.set_title('Genotype Concordance by MAF')
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, 0.5)
            ax.set_ylim(0, 1.05)
            
            # 2. Sample size per MAF bin
            ax = axes[0, 1]
            ax.bar(range(len(concordance_by_maf)), concordance_by_maf['n_variants'],
                  color='coral', edgecolor='black', alpha=0.7)
            ax.set_xlabel('MAF Bin')
            ax.set_ylabel('Number of Variants')
            ax.set_title('Variant Count per MAF Bin')
            ax.set_xticks(range(len(concordance_by_maf)))
            ax.set_xticklabels([f"{row['mean_maf']:.3f}" for _, row in concordance_by_maf.iterrows()],
                              rotation=45, ha='right')
            ax.grid(True, alpha=0.3, axis='y')
            
            # 3. Concordance rate distribution
            ax = axes[1, 0]
            overall_concordance = comparison_df['concordant'].mean()
            
            # Stratified concordance rates
            categories = ['Overall', 'MAF < 0.01', 'MAF 0.01-0.05', 'MAF 0.05-0.1', 'MAF > 0.1']
            rates = [overall_concordance]
            
            # Calculate stratified rates
            for maf_range in [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.5)]:
                mask = (comparison_df['maf_true'] >= maf_range[0]) & (comparison_df['maf_true'] < maf_range[1])
                if mask.any():
                    rates.append(comparison_df[mask]['concordant'].mean())
                else:
                    rates.append(0)
            
            colors = ['green' if r > 0.95 else 'orange' if r > 0.9 else 'red' for r in rates]
            bars = ax.bar(range(len(categories)), rates, color=colors, edgecolor='black', alpha=0.7)
            ax.set_ylabel('Concordance Rate')
            ax.set_title('Concordance Rate by MAF Category')
            ax.set_xticks(range(len(categories)))
            ax.set_xticklabels(categories, rotation=45, ha='right')
            ax.axhline(y=0.95, color='green', linestyle='--', alpha=0.5, label='95% threshold')
            ax.axhline(y=0.90, color='orange', linestyle='--', alpha=0.5, label='90% threshold')
            ax.set_ylim(0, 1.05)
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add value labels on bars
            for bar, rate in zip(bars, rates):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                       f'{rate:.3f}', ha='center', va='bottom', fontsize=9)
        
    else:
        # No validation data - show simulated/expected patterns
        ax = axes[0, 0]
        maf_range = np.linspace(0.01, 0.5, 50)
        # Typical concordance pattern (decreases with lower MAF)
        expected_concordance = 0.98 - 0.3 * np.exp(-10 * maf_range)
        ax.plot(maf_range, expected_concordance, '--', linewidth=2, color='gray',
               label='Expected pattern (no validation data)')
        ax.set_xlabel('Minor Allele Frequency (MAF)')
        ax.set_ylabel('Concordance Rate')
        ax.set_title('Expected Concordance Pattern by MAF')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0.6, 1.05)
        
        # Show message in other panels
        for i, ax in enumerate(axes.flat[1:]):
            ax.text(0.5, 0.5, 'Validation data not available\nProvide --true-vcf for concordance analysis',
                   ha='center', va='center', fontsize=11, style='italic')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
    
    # 4. Summary statistics
    ax = axes[1, 1]
    ax.axis('off')
    
    stats_text = "Concordance Analysis Summary\n" + "="*35 + "\n\n"
    
    if args.true_vcf and Path(args.true_vcf).exists() and 'comparison_df' in locals():
        stats_text += f"Total variants compared: {len(comparison_df):,}\n"
        stats_text += f"Overall concordance: {comparison_df['concordant'].mean():.4f}\n\n"
        
        # MAF distribution
        stats_text += "MAF Distribution:\n"
        stats_text += f"  Mean MAF: {comparison_df['maf_true'].mean():.4f}\n"
        stats_text += f"  Median MAF: {comparison_df['maf_true'].median():.4f}\n"
        stats_text += f"  Rare variants (MAF<0.01): {(comparison_df['maf_true'] < 0.01).sum():,}\n"
        stats_text += f"  Low freq (0.01≤MAF<0.05): {((comparison_df['maf_true'] >= 0.01) & (comparison_df['maf_true'] < 0.05)).sum():,}\n"
        stats_text += f"  Common (MAF≥0.05): {(comparison_df['maf_true'] >= 0.05).sum():,}\n"
    else:
        stats_text += "No validation data provided\n\n"
        stats_text += "To perform concordance analysis:\n"
        stats_text += "1. Provide validation VCF with --true-vcf\n"
        stats_text += "2. Or use masked variant approach\n"
        stats_text += "3. Or use cross-validation results\n\n"
        stats_text += "Concordance analysis helps:\n"
        stats_text += "• Validate imputation accuracy\n"
        stats_text += "• Identify MAF-dependent biases\n"
        stats_text += "• Optimize quality thresholds"
    
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    # Add main title
    chr_info = f" - Chr {args.chr}" if args.chr else " - Genome-wide"
    fig.suptitle(f'Concordance Analysis\nSample: {args.sample_id} | Reference: {args.ref_name}{chr_info}',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(args.output_pdf, dpi=150, bbox_inches='tight')
    print(f"Concordance plot saved to {args.output_pdf}")

if __name__ == "__main__":
    main()