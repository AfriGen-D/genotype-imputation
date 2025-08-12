#!/usr/bin/env python3
"""
Calibration Plot
Plots observed vs expected imputation quality metrics to identify systematic bias
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
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
output_file = "${output_file}"

# Read data
data = read_file(info_file)

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Check for calibration data (expected vs observed R-squared)
has_calibration = False

if 'Rsq' in data.columns and 'EmpRsq' in data.columns:
    # Filter and clean data
    calib_data = data[['Rsq', 'EmpRsq']].copy()
    calib_data = calib_data[(calib_data['Rsq'] != '-') & (calib_data['EmpRsq'] != '-')]
    calib_data['Rsq'] = pd.to_numeric(calib_data['Rsq'], errors='coerce')
    calib_data['EmpRsq'] = pd.to_numeric(calib_data['EmpRsq'], errors='coerce')
    calib_data = calib_data.dropna()
    
    if len(calib_data) > 0:
        has_calibration = True

if has_calibration:
    # Plot 1: Main calibration plot
    ax = axes[0, 0]
    
    # Bin the expected R² and calculate mean observed R²
    n_bins = 20
    bins = np.linspace(0, 1, n_bins + 1)
    calib_data['Expected_bin'] = pd.cut(calib_data['Rsq'], bins)
    
    calibration = calib_data.groupby('Expected_bin').agg({
        'Rsq': ['mean', 'std', 'count'],
        'EmpRsq': ['mean', 'std']
    }).dropna()
    
    # Extract values
    exp_mean = calibration['Rsq']['mean'].values
    obs_mean = calibration['EmpRsq']['mean'].values
    obs_std = calibration['EmpRsq']['std'].values
    counts = calibration['Rsq']['count'].values
    
    # Color by sample size
    sizes = np.clip(counts / counts.max() * 100, 20, 200)
    scatter = ax.scatter(exp_mean, obs_mean, s=sizes, c=counts, 
                        cmap='viridis', alpha=0.6, edgecolors='black')
    
    # Add error bars
    ax.errorbar(exp_mean, obs_mean, yerr=obs_std, fmt='none', 
               alpha=0.3, color='gray', capsize=3)
    
    # Perfect calibration line
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=2, label='Perfect calibration')
    
    # Fit regression line
    if len(exp_mean) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(exp_mean, obs_mean)
        x_fit = np.linspace(0, 1, 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, 'b-', alpha=0.7, linewidth=2, 
               label=f'Fitted: y = {slope:.3f}x + {intercept:.3f}')
        
        # Add statistics
        stats_text = f'Slope: {slope:.3f} ± {std_err:.3f}\n'
        stats_text += f'Intercept: {intercept:.3f}\n'
        stats_text += f'R²: {r_value**2:.3f}\n'
        stats_text += f'p-value: {p_value:.2e}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Expected R² (Estimated)', fontsize=11)
    ax.set_ylabel('Observed R² (Empirical)', fontsize=11)
    ax.set_title('Calibration Plot: Expected vs Observed', fontsize=12)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Number of variants', rotation=270, labelpad=15)
    
    # Plot 2: Calibration error distribution
    ax = axes[0, 1]
    
    calibration_error = calib_data['EmpRsq'] - calib_data['Rsq']
    
    ax.hist(calibration_error, bins=50, edgecolor='black', alpha=0.7, color='salmon')
    ax.axvline(x=0, color='r', linestyle='--', linewidth=2, alpha=0.7, label='No bias')
    ax.set_xlabel('Calibration Error (Observed - Expected)', fontsize=11)
    ax.set_ylabel('Number of Variants', fontsize=11)
    ax.set_title('Calibration Error Distribution', fontsize=12)
    ax.legend()
    
    # Add statistics
    mean_error = calibration_error.mean()
    median_error = calibration_error.median()
    std_error = calibration_error.std()
    
    stats_text = f'Mean error: {mean_error:.4f}\n'
    stats_text += f'Median error: {median_error:.4f}\n'
    stats_text += f'Std dev: {std_error:.4f}\n'
    stats_text += f'N = {len(calibration_error):,}'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
           ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 3: Calibration by MAF
    ax = axes[1, 0]
    
    if 'MAF' in data.columns:
        calib_maf = data[['Rsq', 'EmpRsq', 'MAF']].copy()
        calib_maf = calib_maf[(calib_maf['Rsq'] != '-') & (calib_maf['EmpRsq'] != '-')]
        for col in ['Rsq', 'EmpRsq', 'MAF']:
            calib_maf[col] = pd.to_numeric(calib_maf[col], errors='coerce')
        calib_maf = calib_maf.dropna()
        
        # Define MAF bins
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%']
        calib_maf['MAF_bin'] = pd.cut(calib_maf['MAF'], bins=maf_bins, labels=maf_labels)
        
        # Calculate calibration slope for each MAF bin
        slopes = []
        maf_groups = []
        for maf_group, group_data in calib_maf.groupby('MAF_bin'):
            if len(group_data) > 10:  # Need enough data points
                slope, _, _, _, _ = stats.linregress(group_data['Rsq'], group_data['EmpRsq'])
                slopes.append(slope)
                maf_groups.append(maf_group)
        
        if slopes:
            colors = plt.cm.coolwarm(np.linspace(0.2, 0.8, len(slopes)))
            bars = ax.bar(range(len(slopes)), slopes, color=colors, edgecolor='black', alpha=0.7)
            ax.set_xticks(range(len(slopes)))
            ax.set_xticklabels(maf_groups, rotation=45, ha='right')
            ax.axhline(y=1, color='g', linestyle='--', alpha=0.5, label='Perfect calibration')
            ax.set_xlabel('MAF Bin', fontsize=11)
            ax.set_ylabel('Calibration Slope', fontsize=11)
            ax.set_title('Calibration Slope by MAF', fontsize=12)
            ax.legend()
            
            # Add value labels
            for bar, val in zip(bars, slopes):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:.3f}', ha='center', va='bottom' if val > 0 else 'top', fontsize=9)
    else:
        ax.text(0.5, 0.5, 'MAF data not available', ha='center', va='center')
    
    # Plot 4: Q-Q plot
    ax = axes[1, 1]
    
    # Q-Q plot of observed vs expected
    sorted_exp = np.sort(calib_data['Rsq'])
    sorted_obs = np.sort(calib_data['EmpRsq'])
    
    ax.scatter(sorted_exp, sorted_obs, alpha=0.5, s=1)
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, linewidth=2)
    ax.set_xlabel('Expected R² Quantiles', fontsize=11)
    ax.set_ylabel('Observed R² Quantiles', fontsize=11)
    ax.set_title('Q-Q Plot: Expected vs Observed', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    
else:
    # No calibration data available
    for ax in axes.flat:
        ax.text(0.5, 0.5, 'Calibration data not available\n\n' +
                'Requires both estimated (Rsq) and\nempirical (EmpRsq) R² values\n\n' +
                'EmpRsq is typically available only\nfor genotyped variants',
                ha='center', va='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_title('Calibration Analysis', fontsize=12)

plt.suptitle('Imputation Quality Calibration Assessment', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Calibration plot saved to {output_file}")