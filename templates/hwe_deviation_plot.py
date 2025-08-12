#!/usr/bin/env python3
"""
Hardy-Weinberg Equilibrium Deviation Plot
Checks for deviations in HWE for imputed vs genotyped variants
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import gzip
from scipy import stats

def read_file(filename):
    """Read potentially gzipped files"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

def calculate_hwe_pvalue(obs_hom_ref, obs_het, obs_hom_alt):
    """
    Calculate HWE p-value using chi-square test
    """
    n = obs_hom_ref + obs_het + obs_hom_alt
    if n == 0:
        return np.nan
    
    # Calculate allele frequencies
    p = (2 * obs_hom_ref + obs_het) / (2 * n)
    q = 1 - p
    
    # Expected counts under HWE
    exp_hom_ref = n * p * p
    exp_het = n * 2 * p * q
    exp_hom_alt = n * q * q
    
    # Chi-square test
    observed = [obs_hom_ref, obs_het, obs_hom_alt]
    expected = [exp_hom_ref, exp_het, exp_hom_alt]
    
    # Avoid division by zero
    if any(e == 0 for e in expected):
        return np.nan
    
    chi2 = sum((o - e) ** 2 / e for o, e in zip(observed, expected))
    p_value = 1 - stats.chi2.cdf(chi2, df=1)
    
    return p_value

# Parse arguments from Nextflow template
info_file = "${info_file}"
vcf_file = "${vcf_file}"  # Optional: VCF with actual genotype counts
output_file = "${output_file}"

# Read data
data = read_file(info_file)

# Create figure
fig, axes = plt.subplots(2, 3, figsize=(14, 10))

# Check if we have HWE data or can calculate it
has_hwe = False
if 'HWE_pval' in data.columns or 'HWE' in data.columns:
    hwe_col = 'HWE_pval' if 'HWE_pval' in data.columns else 'HWE'
    data[hwe_col] = pd.to_numeric(data[hwe_col].replace('-', np.nan), errors='coerce')
    has_hwe = True

# If no HWE but we have MAF, we can show expected vs observed heterozygosity
if not has_hwe and 'MAF' in data.columns:
    data['MAF'] = pd.to_numeric(data['MAF'], errors='coerce')
    data['Expected_Het'] = 2 * data['MAF'] * (1 - data['MAF'])
    
    # This is a proxy for HWE deviation
    if 'Het_Rate' in data.columns:
        data['Het_Rate'] = pd.to_numeric(data['Het_Rate'], errors='coerce')
        data['HWE_deviation'] = data['Het_Rate'] - data['Expected_Het']
        has_hwe = True

if has_hwe and 'Genotyped' in data.columns:
    # Separate by variant type
    genotyped = data[data['Genotyped'] == 'Genotyped'].copy()
    imputed = data[data['Genotyped'] == 'Imputed'].copy()
    
    # Plot 1: Q-Q plot of HWE p-values
    ax = axes[0, 0]
    
    if 'HWE_pval' in data.columns or 'HWE' in data.columns:
        hwe_col = 'HWE_pval' if 'HWE_pval' in data.columns else 'HWE'
        
        # Get -log10 p-values
        gen_pvals = genotyped[hwe_col].dropna()
        imp_pvals = imputed[hwe_col].dropna()
        
        if len(gen_pvals) > 0:
            gen_log = -np.log10(gen_pvals + 1e-300)  # Add small value to avoid log(0)
            theoretical = -np.log10(np.linspace(1/len(gen_log), 1, len(gen_log)))
            ax.scatter(np.sort(theoretical), np.sort(gen_log), alpha=0.5, s=5, 
                      label=f'Genotyped (n={len(gen_log):,})', color='blue')
        
        if len(imp_pvals) > 0:
            imp_log = -np.log10(imp_pvals + 1e-300)
            theoretical = -np.log10(np.linspace(1/len(imp_log), 1, len(imp_log)))
            ax.scatter(np.sort(theoretical), np.sort(imp_log), alpha=0.5, s=5,
                      label=f'Imputed (n={len(imp_log):,})', color='red')
        
        # Add diagonal line
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=1)
        
        # Add significance threshold
        bonferroni = -np.log10(0.05 / len(data))
        ax.axhline(y=bonferroni, color='gray', linestyle=':', alpha=0.5,
                  label=f'Bonferroni (p={0.05/len(data):.2e})')
        
        ax.set_xlabel('Expected -log10(p-value)', fontsize=11)
        ax.set_ylabel('Observed -log10(p-value)', fontsize=11)
        ax.set_title('HWE Q-Q Plot', fontsize=12)
        ax.legend(loc='upper left', fontsize=9)
    
    # Plot 2: Distribution of HWE p-values
    ax = axes[0, 1]
    
    if 'HWE_pval' in data.columns or 'HWE' in data.columns:
        hwe_col = 'HWE_pval' if 'HWE_pval' in data.columns else 'HWE'
        
        bins = np.linspace(0, 1, 21)
        
        if len(genotyped) > 0:
            gen_pvals = genotyped[hwe_col].dropna()
            if len(gen_pvals) > 0:
                ax.hist(gen_pvals, bins=bins, alpha=0.5, label='Genotyped', 
                       color='blue', edgecolor='black', density=True)
        
        if len(imputed) > 0:
            imp_pvals = imputed[hwe_col].dropna()
            if len(imp_pvals) > 0:
                ax.hist(imp_pvals, bins=bins, alpha=0.5, label='Imputed',
                       color='red', edgecolor='black', density=True)
        
        # Add uniform distribution line (expected under null)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Uniform (expected)')
        
        ax.set_xlabel('HWE p-value', fontsize=11)
        ax.set_ylabel('Density', fontsize=11)
        ax.set_title('HWE p-value Distribution', fontsize=12)
        ax.legend()
    
    # Plot 3: Proportion of variants failing HWE by MAF
    ax = axes[0, 2]
    
    if 'MAF' in data.columns and ('HWE_pval' in data.columns or 'HWE' in data.columns):
        hwe_col = 'HWE_pval' if 'HWE_pval' in data.columns else 'HWE'
        
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%']
        
        genotyped['MAF_bin'] = pd.cut(genotyped['MAF'], bins=maf_bins, labels=maf_labels)
        imputed['MAF_bin'] = pd.cut(imputed['MAF'], bins=maf_bins, labels=maf_labels)
        
        # Calculate proportion failing HWE (p < 0.05) by MAF bin
        gen_fail = genotyped.groupby('MAF_bin').apply(
            lambda x: (x[hwe_col] < 0.05).mean() if len(x) > 0 else 0
        )
        imp_fail = imputed.groupby('MAF_bin').apply(
            lambda x: (x[hwe_col] < 0.05).mean() if len(x) > 0 else 0
        )
        
        x = np.arange(len(maf_labels))
        width = 0.35
        
        ax.bar(x - width/2, gen_fail, width, label='Genotyped',
              color='lightblue', edgecolor='black', alpha=0.7)
        ax.bar(x + width/2, imp_fail, width, label='Imputed',
              color='lightcoral', edgecolor='black', alpha=0.7)
        
        ax.axhline(y=0.05, color='gray', linestyle='--', alpha=0.5,
                  label='Expected (5%)')
        
        ax.set_xlabel('MAF Bin', fontsize=11)
        ax.set_ylabel('Proportion Failing HWE (p<0.05)', fontsize=11)
        ax.set_title('HWE Violations by MAF', fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(maf_labels, rotation=45, ha='right')
        ax.legend()
    
    # Plot 4: Inbreeding coefficient (F) distribution
    ax = axes[1, 0]
    
    if 'F' in data.columns or 'Fis' in data.columns:
        f_col = 'F' if 'F' in data.columns else 'Fis'
        
        gen_f = genotyped[f_col].dropna()
        imp_f = imputed[f_col].dropna()
        
        data_to_plot = []
        labels = []
        
        if len(gen_f) > 0:
            data_to_plot.append(gen_f)
            labels.append('Genotyped')
        if len(imp_f) > 0:
            data_to_plot.append(imp_f)
            labels.append('Imputed')
        
        if data_to_plot:
            bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
            colors = ['lightblue', 'lightcoral']
            for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                patch.set_facecolor(color)
            
            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax.set_ylabel('Inbreeding Coefficient (F)', fontsize=11)
            ax.set_title('Inbreeding Coefficient Distribution', fontsize=12)
    else:
        # Calculate pseudo-F from heterozygosity deficit
        if 'MAF' in data.columns:
            ax.text(0.5, 0.5, 'Inbreeding coefficient not available\n' +
                   'Using heterozygosity-based estimation',
                   ha='center', va='center', fontsize=10)
    
    # Plot 5: HWE test statistics summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = 'HWE Test Summary\n' + '='*40 + '\n\n'
    
    if 'HWE_pval' in data.columns or 'HWE' in data.columns:
        hwe_col = 'HWE_pval' if 'HWE_pval' in data.columns else 'HWE'
        
        # Genotyped variants
        gen_pvals = genotyped[hwe_col].dropna()
        if len(gen_pvals) > 0:
            summary_text += 'Genotyped Variants:\n'
            summary_text += f'  Total tested: {len(gen_pvals):,}\n'
            summary_text += f'  Failed HWE (p<0.05): {(gen_pvals < 0.05).sum():,} ({100*(gen_pvals < 0.05).mean():.2f}%)\n'
            summary_text += f'  Failed HWE (p<0.001): {(gen_pvals < 0.001).sum():,} ({100*(gen_pvals < 0.001).mean():.2f}%)\n'
            summary_text += f'  Median p-value: {gen_pvals.median():.3f}\n\n'
        
        # Imputed variants
        imp_pvals = imputed[hwe_col].dropna()
        if len(imp_pvals) > 0:
            summary_text += 'Imputed Variants:\n'
            summary_text += f'  Total tested: {len(imp_pvals):,}\n'
            summary_text += f'  Failed HWE (p<0.05): {(imp_pvals < 0.05).sum():,} ({100*(imp_pvals < 0.05).mean():.2f}%)\n'
            summary_text += f'  Failed HWE (p<0.001): {(imp_pvals < 0.001).sum():,} ({100*(imp_pvals < 0.001).mean():.2f}%)\n'
            summary_text += f'  Median p-value: {imp_pvals.median():.3f}\n\n'
        
        # Statistical comparison
        if len(gen_pvals) > 0 and len(imp_pvals) > 0:
            # KS test to compare distributions
            ks_stat, ks_pval = stats.ks_2samp(gen_pvals, imp_pvals)
            summary_text += 'Distribution Comparison:\n'
            summary_text += f'  KS statistic: {ks_stat:.4f}\n'
            summary_text += f'  KS p-value: {ks_pval:.2e}\n'
    
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', family='monospace')
    
    # Plot 6: MAF vs heterozygosity (HWE expectation)
    ax = axes[1, 2]
    
    if 'MAF' in data.columns:
        # Theoretical HWE curve
        maf_range = np.linspace(0, 0.5, 100)
        expected_het = 2 * maf_range * (1 - maf_range)
        ax.plot(maf_range, expected_het, 'k--', alpha=0.5, linewidth=2, label='HWE expectation')
        
        # Sample points for visualization
        n_sample = min(5000, len(data))
        
        if len(genotyped) > 0:
            gen_sample = genotyped.sample(n=min(n_sample, len(genotyped)))
            gen_sample['Expected_Het'] = 2 * gen_sample['MAF'] * (1 - gen_sample['MAF'])
            ax.scatter(gen_sample['MAF'], gen_sample['Expected_Het'],
                      alpha=0.3, s=2, c='blue', label='Genotyped')
        
        if len(imputed) > 0:
            imp_sample = imputed.sample(n=min(n_sample, len(imputed)))
            imp_sample['Expected_Het'] = 2 * imp_sample['MAF'] * (1 - imp_sample['MAF'])
            ax.scatter(imp_sample['MAF'], imp_sample['Expected_Het'],
                      alpha=0.3, s=2, c='red', label='Imputed')
        
        ax.set_xlabel('Minor Allele Frequency', fontsize=11)
        ax.set_ylabel('Expected Heterozygosity (2pq)', fontsize=11)
        ax.set_title('HWE Expectation Check', fontsize=12)
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 0.5)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)

else:
    # No HWE data available
    for ax in axes.flat:
        ax.text(0.5, 0.5, 'HWE test data not available\n\n' +
                'To calculate HWE statistics:\n' +
                '1. Count genotypes (AA, Aa, aa) per variant\n' +
                '2. Calculate expected counts under HWE\n' +
                '3. Perform chi-square test\n' +
                '4. Compare genotyped vs imputed variants',
                ha='center', va='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_title('HWE Analysis', fontsize=12)

plt.suptitle('Hardy-Weinberg Equilibrium Assessment', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"HWE deviation plot saved to {output_file}")