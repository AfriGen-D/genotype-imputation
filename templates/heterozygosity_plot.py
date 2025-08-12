#!/usr/bin/env python3
"""
Heterozygosity Rate Plot
Compares heterozygosity rates between imputed and genotyped variants
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

# Parse arguments from Nextflow template
info_file = "${info_file}"
vcf_file = "${vcf_file}"  # Optional: actual VCF with genotypes for true het calculation
output_file = "${output_file}"

# Read data
data = read_file(info_file)

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Check if we have genotype status
if 'Genotyped' in data.columns:
    # Separate by variant type
    genotyped = data[data['Genotyped'] == 'Genotyped'].copy()
    imputed = data[data['Genotyped'] == 'Imputed'].copy()
    
    # Calculate expected heterozygosity from MAF (2pq under HWE)
    if 'MAF' in data.columns:
        for df in [genotyped, imputed]:
            df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce')
            df['Expected_Het'] = 2 * df['MAF'] * (1 - df['MAF'])
        
        # Plot 1: Expected heterozygosity distribution
        ax = axes[0, 0]
        
        data_to_plot = []
        labels = []
        colors = []
        
        if len(genotyped) > 0 and 'Expected_Het' in genotyped.columns:
            exp_het_gen = genotyped['Expected_Het'].dropna()
            if len(exp_het_gen) > 0:
                data_to_plot.append(exp_het_gen)
                labels.append(f'Genotyped\n(n={len(exp_het_gen):,})')
                colors.append('lightblue')
        
        if len(imputed) > 0 and 'Expected_Het' in imputed.columns:
            exp_het_imp = imputed['Expected_Het'].dropna()
            if len(exp_het_imp) > 0:
                data_to_plot.append(exp_het_imp)
                labels.append(f'Imputed\n(n={len(exp_het_imp):,})')
                colors.append('lightcoral')
        
        if data_to_plot:
            bp = ax.violinplot(data_to_plot, positions=range(len(data_to_plot)), 
                              widths=0.7, showmeans=True, showmedians=True)
            
            # Color the violin plots
            for pc, color in zip(bp['bodies'], colors):
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
            
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels)
            ax.set_ylabel('Expected Heterozygosity (2pq)', fontsize=11)
            ax.set_title('Expected Heterozygosity Distribution', fontsize=12)
            
            # Add statistics
            for i, (data_vals, label) in enumerate(zip(data_to_plot, labels)):
                mean_val = data_vals.mean()
                median_val = data_vals.median()
                ax.text(i, ax.get_ylim()[1] * 0.95, 
                       f'μ={mean_val:.3f}\nM={median_val:.3f}',
                       ha='center', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Plot 2: Expected het vs MAF
        ax = axes[0, 1]
        
        # Theoretical curve
        maf_range = np.linspace(0, 0.5, 100)
        theoretical_het = 2 * maf_range * (1 - maf_range)
        ax.plot(maf_range, theoretical_het, 'k--', alpha=0.5, linewidth=2, 
               label='Theoretical (HWE)')
        
        # Scatter plot for actual data
        if len(genotyped) > 0:
            # Subsample if too many points
            n_points = min(5000, len(genotyped))
            sample_idx = np.random.choice(len(genotyped), n_points, replace=False)
            genotyped_sample = genotyped.iloc[sample_idx]
            ax.scatter(genotyped_sample['MAF'], genotyped_sample['Expected_Het'],
                      alpha=0.3, s=5, c='blue', label='Genotyped')
        
        if len(imputed) > 0:
            n_points = min(5000, len(imputed))
            sample_idx = np.random.choice(len(imputed), n_points, replace=False)
            imputed_sample = imputed.iloc[sample_idx]
            ax.scatter(imputed_sample['MAF'], imputed_sample['Expected_Het'],
                      alpha=0.3, s=5, c='red', label='Imputed')
        
        ax.set_xlabel('Minor Allele Frequency', fontsize=11)
        ax.set_ylabel('Expected Heterozygosity', fontsize=11)
        ax.set_title('Heterozygosity vs MAF', fontsize=12)
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 0.5)
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Heterozygosity by MAF bins
        ax = axes[1, 0]
        
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%']
        
        genotyped['MAF_bin'] = pd.cut(genotyped['MAF'], bins=maf_bins, labels=maf_labels)
        imputed['MAF_bin'] = pd.cut(imputed['MAF'], bins=maf_bins, labels=maf_labels)
        
        # Calculate mean expected het by MAF bin
        gen_by_maf = genotyped.groupby('MAF_bin')['Expected_Het'].agg(['mean', 'std', 'count'])
        imp_by_maf = imputed.groupby('MAF_bin')['Expected_Het'].agg(['mean', 'std', 'count'])
        
        x = np.arange(len(maf_labels))
        width = 0.35
        
        # Plot bars
        if not gen_by_maf.empty:
            ax.bar(x - width/2, gen_by_maf['mean'], width, yerr=gen_by_maf['std'],
                  label='Genotyped', color='lightblue', edgecolor='black', alpha=0.7, capsize=3)
        
        if not imp_by_maf.empty:
            ax.bar(x + width/2, imp_by_maf['mean'], width, yerr=imp_by_maf['std'],
                  label='Imputed', color='lightcoral', edgecolor='black', alpha=0.7, capsize=3)
        
        ax.set_xlabel('MAF Bin', fontsize=11)
        ax.set_ylabel('Mean Expected Heterozygosity', fontsize=11)
        ax.set_title('Heterozygosity by MAF Category', fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(maf_labels, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Plot 4: Statistical comparison
        ax = axes[1, 1]
        
        # Perform statistical tests
        results_text = 'Statistical Comparison\n' + '='*30 + '\n\n'
        
        if len(genotyped) > 0 and len(imputed) > 0:
            # Overall comparison
            gen_het = genotyped['Expected_Het'].dropna()
            imp_het = imputed['Expected_Het'].dropna()
            
            # T-test
            t_stat, p_value = stats.ttest_ind(gen_het, imp_het)
            results_text += f'Overall Comparison:\n'
            results_text += f'Genotyped mean: {gen_het.mean():.4f}\n'
            results_text += f'Imputed mean: {imp_het.mean():.4f}\n'
            results_text += f'Difference: {gen_het.mean() - imp_het.mean():.4f}\n'
            results_text += f't-statistic: {t_stat:.3f}\n'
            results_text += f'p-value: {p_value:.2e}\n\n'
            
            # By MAF bin comparison
            results_text += 'By MAF Bin (mean difference):\n'
            for maf_bin in maf_labels:
                gen_bin = genotyped[genotyped['MAF_bin'] == maf_bin]['Expected_Het'].dropna()
                imp_bin = imputed[imputed['MAF_bin'] == maf_bin]['Expected_Het'].dropna()
                if len(gen_bin) > 0 and len(imp_bin) > 0:
                    diff = gen_bin.mean() - imp_bin.mean()
                    results_text += f'{maf_bin}: {diff:+.4f}\n'
        
        ax.text(0.1, 0.9, results_text, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', family='monospace')
        ax.set_title('Statistical Analysis', fontsize=12)
        ax.axis('off')
    
    else:
        # No MAF data available
        for ax in axes.flat:
            ax.text(0.5, 0.5, 'MAF data not available\nfor heterozygosity calculation',
                   ha='center', va='center', fontsize=11)
            ax.set_title('Heterozygosity Analysis', fontsize=12)

else:
    # No genotype status available
    for ax in axes.flat:
        ax.text(0.5, 0.5, 'Genotype status not available\n\n' +
                'Cannot separate genotyped vs imputed variants',
                ha='center', va='center', fontsize=11,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_title('Heterozygosity Analysis', fontsize=12)

plt.suptitle('Heterozygosity Rate Comparison', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Heterozygosity plot saved to {output_file}")