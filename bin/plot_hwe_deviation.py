#!/usr/bin/env python3
"""
Hardy-Weinberg Equilibrium Deviation Plot
Detects population structure, technical issues, or selection
Analyzes departure from HWE expectations
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
from pathlib import Path
from scipy import stats

def calculate_allele_frequencies(genotypes):
    """Calculate allele frequencies from genotype counts"""
    n_samples = sum(genotypes.values())
    if n_samples == 0:
        return 0, 0
    
    # Count alleles
    n_ref = 2 * genotypes.get('0/0', 0) + genotypes.get('0/1', 0)
    n_alt = 2 * genotypes.get('1/1', 0) + genotypes.get('0/1', 0)
    total_alleles = 2 * n_samples
    
    if total_alleles == 0:
        return 0, 0
    
    p = n_ref / total_alleles  # ref allele frequency
    q = n_alt / total_alleles  # alt allele frequency
    
    return p, q

def hwe_chi_square(observed, expected):
    """Calculate chi-square test for HWE"""
    chi_sq = 0
    for geno in ['0/0', '0/1', '1/1']:
        obs = observed.get(geno, 0)
        exp = expected.get(geno, 0)
        if exp > 0:
            chi_sq += ((obs - exp) ** 2) / exp
    
    # df = 1 for HWE test
    p_value = 1 - stats.chi2.cdf(chi_sq, df=1)
    return chi_sq, p_value

def calculate_hwe_statistics(vcf_file):
    """Calculate HWE statistics for each variant"""
    hwe_stats = []
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    n_samples = len(header) - 9 if len(header) > 9 else 0
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, variant_id, ref, alt = fields[:5]
            
            # Skip multiallelic variants
            if ',' in alt:
                continue
            
            format_field = fields[8]
            if 'GT' not in format_field:
                continue
            
            format_keys = format_field.split(':')
            gt_idx = format_keys.index('GT')
            
            # Count genotypes
            genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, 'missing': 0}
            
            for sample_field in fields[9:]:
                values = sample_field.split(':')
                if len(values) > gt_idx:
                    gt = values[gt_idx]
                    
                    # Normalize genotype format
                    if gt in ['.', './.', '.|.']:
                        genotype_counts['missing'] += 1
                    else:
                        # Convert phased to unphased for counting
                        gt_norm = gt.replace('|', '/')
                        if gt_norm == '0/0':
                            genotype_counts['0/0'] += 1
                        elif gt_norm in ['0/1', '1/0']:
                            genotype_counts['0/1'] += 1
                        elif gt_norm == '1/1':
                            genotype_counts['1/1'] += 1
            
            # Calculate allele frequencies
            valid_genotypes = {k: v for k, v in genotype_counts.items() if k != 'missing'}
            n_called = sum(valid_genotypes.values())
            
            if n_called < 10:  # Skip variants with too few calls
                continue
            
            p, q = calculate_allele_frequencies(valid_genotypes)
            
            # Calculate expected genotype frequencies under HWE
            expected = {
                '0/0': p * p * n_called,
                '0/1': 2 * p * q * n_called,
                '1/1': q * q * n_called
            }
            
            # Calculate chi-square test
            chi_sq, p_value = hwe_chi_square(valid_genotypes, expected)
            
            # Calculate inbreeding coefficient
            obs_het = genotype_counts['0/1'] / n_called if n_called > 0 else 0
            exp_het = 2 * p * q
            f_stat = 1 - (obs_het / exp_het) if exp_het > 0 else 0
            
            hwe_stats.append({
                'position': f"{chrom}:{pos}",
                'variant_id': variant_id if variant_id != '.' else f"{chrom}_{pos}",
                'p_allele_freq': p,
                'q_allele_freq': q,
                'maf': min(p, q),
                'obs_hom_ref': genotype_counts['0/0'],
                'obs_het': genotype_counts['0/1'],
                'obs_hom_alt': genotype_counts['1/1'],
                'exp_hom_ref': expected['0/0'],
                'exp_het': expected['0/1'],
                'exp_hom_alt': expected['1/1'],
                'chi_square': chi_sq,
                'p_value': p_value,
                'f_statistic': f_stat,
                'n_samples': n_called,
                'call_rate': n_called / (n_called + genotype_counts['missing'])
            })
    
    return pd.DataFrame(hwe_stats)

def main():
    parser = argparse.ArgumentParser(description='Plot Hardy-Weinberg equilibrium deviations')
    parser.add_argument('vcf_file', help='Input VCF file')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chr', help='Chromosome', default='')
    parser.add_argument('--p-threshold', type=float, default=1e-6, 
                       help='P-value threshold for significant HWE deviation')
    
    args = parser.parse_args()
    
    print(f"Calculating HWE statistics from {args.vcf_file}...")
    hwe_df = calculate_hwe_statistics(args.vcf_file)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    if len(hwe_df) > 0:
        # Apply Bonferroni correction
        bonferroni_threshold = 0.05 / len(hwe_df)
        hwe_df['significant'] = hwe_df['p_value'] < args.p_threshold
        hwe_df['bonferroni_significant'] = hwe_df['p_value'] < bonferroni_threshold
        
        # 1. Q-Q plot for HWE p-values
        ax = axes[0, 0]
        
        # Calculate expected and observed -log10(p) values
        observed_p = np.sort(hwe_df['p_value'])
        n = len(observed_p)
        expected_p = np.arange(1, n + 1) / (n + 1)
        
        # Plot Q-Q plot
        ax.scatter(-np.log10(expected_p), -np.log10(observed_p),
                  alpha=0.5, s=10, c='steelblue')
        
        # Add diagonal line for null expectation
        max_val = max(-np.log10(expected_p).max(), -np.log10(observed_p).max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Expected under null')
        
        # Add significance threshold
        ax.axhline(y=-np.log10(args.p_threshold), color='orange', linestyle='--',
                  alpha=0.5, label=f'p = {args.p_threshold}')
        
        ax.set_xlabel('Expected -log₁₀(p)')
        ax.set_ylabel('Observed -log₁₀(p)')
        ax.set_title('HWE P-value Q-Q Plot')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. P-value histogram
        ax = axes[0, 1]
        
        ax.hist(hwe_df['p_value'], bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.axvline(x=args.p_threshold, color='red', linestyle='--',
                  label=f'Significance threshold\n(p = {args.p_threshold})', alpha=0.5)
        ax.axvline(x=bonferroni_threshold, color='orange', linestyle='--',
                  label=f'Bonferroni corrected\n(p = {bonferroni_threshold:.2e})', alpha=0.5)
        
        ax.set_xlabel('HWE P-value')
        ax.set_ylabel('Number of Variants')
        ax.set_title('Distribution of HWE P-values')
        ax.set_yscale('log')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # 3. Inbreeding coefficient by MAF
        ax = axes[1, 0]
        
        # Bin by MAF
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
        hwe_df['maf_bin'] = pd.cut(hwe_df['maf'], bins=maf_bins, include_lowest=True)
        
        # Calculate mean F per MAF bin
        f_by_maf = hwe_df.groupby('maf_bin')['f_statistic'].agg(['mean', 'std', 'count'])
        
        x_pos = np.arange(len(f_by_maf))
        bars = ax.bar(x_pos, f_by_maf['mean'], yerr=f_by_maf['std'],
                     capsize=5, color='green', edgecolor='black', alpha=0.7)
        
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        ax.set_xlabel('MAF Bin')
        ax.set_ylabel('Mean Inbreeding Coefficient (F)')
        ax.set_title('Inbreeding Coefficient by Minor Allele Frequency')
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(idx) for idx in f_by_maf.index], rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add sample counts
        for i, (bar, count) in enumerate(zip(bars, f_by_maf['count'])):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                   f'n={count}', ha='center', va='bottom', fontsize=8)
        
        # 4. Summary statistics
        ax = axes[1, 1]
        ax.axis('off')
        
        n_sig = hwe_df['significant'].sum()
        n_bonf_sig = hwe_df['bonferroni_significant'].sum()
        pct_sig = (n_sig / len(hwe_df)) * 100
        
        stats_text = "HWE Analysis Summary\n" + "="*35 + "\n\n"
        stats_text += f"Total variants tested: {len(hwe_df):,}\n"
        stats_text += f"Mean call rate: {hwe_df['call_rate'].mean():.3f}\n\n"
        
        stats_text += "HWE Deviations:\n"
        stats_text += f"  Significant (p < {args.p_threshold}): {n_sig:,} ({pct_sig:.1f}%)\n"
        stats_text += f"  Bonferroni significant: {n_bonf_sig:,}\n"
        stats_text += f"  Mean chi-square: {hwe_df['chi_square'].mean():.2f}\n\n"
        
        stats_text += "Inbreeding Coefficient:\n"
        stats_text += f"  Mean F: {hwe_df['f_statistic'].mean():.4f}\n"
        stats_text += f"  Std F: {hwe_df['f_statistic'].std():.4f}\n"
        stats_text += f"  Variants with F > 0.1: {(hwe_df['f_statistic'] > 0.1).sum():,}\n"
        stats_text += f"  Variants with F < -0.1: {(hwe_df['f_statistic'] < -0.1).sum():,}\n\n"
        
        stats_text += "Interpretation:\n"
        if pct_sig < 5:
            stats_text += "✓ Low HWE deviation rate\n"
            stats_text += "  Consistent with single population\n"
        elif pct_sig < 10:
            stats_text += "⚠ Moderate HWE deviation rate\n"
            stats_text += "  May indicate:\n"
            stats_text += "  • Mild population structure\n"
            stats_text += "  • Technical artifacts\n"
        else:
            stats_text += "⚠ High HWE deviation rate\n"
            stats_text += "  Likely indicates:\n"
            stats_text += "  • Population stratification\n"
            stats_text += "  • Sample mixing\n"
            stats_text += "  • Genotyping errors\n"
        
        mean_f = hwe_df['f_statistic'].mean()
        if abs(mean_f) < 0.05:
            stats_text += "\n✓ F-statistic near zero\n"
        elif mean_f > 0.05:
            stats_text += "\n⚠ Positive F indicates inbreeding\n"
        else:
            stats_text += "\n⚠ Negative F indicates excess heterozygosity\n"
    else:
        # No data
        for ax in axes.flat:
            ax.text(0.5, 0.5, 'No HWE data available',
                   ha='center', va='center', fontsize=12)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
        
        stats_text = "No variants found for HWE analysis"
    
    if len(hwe_df) > 0 and 'stats_text' in locals():
        axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                       fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    # Add main title
    chr_info = f" - Chr {args.chr}" if args.chr else " - Genome-wide"
    fig.suptitle(f'Hardy-Weinberg Equilibrium Analysis\nDataset: {args.sample_id} | Reference: {args.ref_name}{chr_info}',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(args.output_pdf, dpi=150, bbox_inches='tight')
    print(f"HWE plot saved to {args.output_pdf}")
    
    # Save significant variants list
    if len(hwe_df) > 0 and n_sig > 0:
        sig_file = args.output_pdf.replace('.pdf', '_significant.txt')
        sig_variants = hwe_df[hwe_df['significant']][
            ['variant_id', 'position', 'p_value', 'chi_square', 'f_statistic']
        ]
        sig_variants.to_csv(sig_file, sep='\t', index=False)
        print(f"Significant variants saved to {sig_file}")

if __name__ == "__main__":
    main()