#!/usr/bin/env python3
"""
Ancestry-aware Imputation Quality Assessment
Evaluates imputation performance stratified by ancestry components
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import gzip
from scipy import stats

# Set style for professional-looking plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def read_file(filename):
    """Helper function to read potentially gzipped files"""
    try:
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as f:
                return pd.read_csv(f, sep='\t')
        else:
            return pd.read_csv(filename, sep='\t')
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return pd.DataFrame()

def calculate_ancestry_specific_metrics(info_data, ancestry_data=None):
    """
    Calculate imputation quality metrics stratified by ancestry
    
    Args:
        info_data: DataFrame with imputation info
        ancestry_data: Optional DataFrame with ancestry proportions per sample
    
    Returns:
        DataFrame with ancestry-specific metrics
    """
    metrics = {}
    
    # Basic metrics without ancestry stratification
    metrics['overall'] = {
        'mean_r2': info_data['Rsq'].mean(),
        'median_r2': info_data['Rsq'].median(),
        'prop_well_imputed': (info_data['Rsq'] >= 0.8).mean(),
        'prop_moderate': ((info_data['Rsq'] >= 0.3) & (info_data['Rsq'] < 0.8)).mean(),
        'prop_poor': (info_data['Rsq'] < 0.3).mean()
    }
    
    if ancestry_data is not None:
        # Calculate metrics for different ancestry groups
        # This would require sample-level imputation data
        pass
    
    return pd.DataFrame(metrics).T

def plot_imputation_by_maf_ancestry(info_data, output_file):
    """
    Create plots showing imputation quality by MAF, potentially stratified by ancestry
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. R-squared distribution by MAF bins
    ax = axes[0, 0]
    maf_bins = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
    info_data['MAF_bin'] = pd.cut(info_data['MAF'], bins=maf_bins, 
                                   labels=['Rare\n(0-0.1%)', 'Very Low\n(0.1-1%)', 
                                          'Low\n(1-5%)', 'Medium\n(5-10%)', 
                                          'Common\n(10-20%)', 'Very Common\n(20-50%)'])
    
    # Box plot of R-squared by MAF bin
    info_data.boxplot(column='Rsq', by='MAF_bin', ax=ax)
    ax.set_xlabel('MAF Bin', fontsize=11)
    ax.set_ylabel('Imputation R²', fontsize=11)
    ax.set_title('Imputation Quality by MAF', fontsize=12)
    ax.axhline(y=0.3, color='r', linestyle='--', alpha=0.5, label='Quality threshold')
    ax.axhline(y=0.8, color='g', linestyle='--', alpha=0.5, label='High quality')
    ax.legend()
    plt.sca(ax)
    plt.xticks(rotation=45)
    
    # 2. Cumulative distribution of R-squared
    ax = axes[0, 1]
    sorted_rsq = np.sort(info_data['Rsq'].dropna())
    cumulative = np.arange(1, len(sorted_rsq) + 1) / len(sorted_rsq)
    ax.plot(sorted_rsq, cumulative, linewidth=2)
    ax.set_xlabel('Imputation R²', fontsize=11)
    ax.set_ylabel('Cumulative Proportion', fontsize=11)
    ax.set_title('Cumulative Distribution of R²', fontsize=12)
    ax.axvline(x=0.3, color='r', linestyle='--', alpha=0.5)
    ax.axvline(x=0.8, color='g', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3)
    
    # 3. Density plot of imputation certainty (dosage distribution)
    ax = axes[1, 0]
    # If we have dosage data, we would plot it here
    # For now, show R-squared density by genotyped status
    for status in info_data['Genotyped'].unique():
        subset = info_data[info_data['Genotyped'] == status]['Rsq'].dropna()
        if len(subset) > 0:
            ax.hist(subset, bins=30, alpha=0.5, label=status, density=True)
    ax.set_xlabel('Imputation R²', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('R² Distribution by Variant Type', fontsize=12)
    ax.legend()
    
    # 4. MAF comparison: imputed vs reference
    ax = axes[1, 1]
    # Scatter plot of MAF if we have both reference and imputed MAF
    if 'ALT_Frq' in info_data.columns:
        ax.scatter(info_data['MAF'], info_data['ALT_Frq'], 
                  c=info_data['Rsq'], cmap='coolwarm', 
                  s=1, alpha=0.5)
        ax.plot([0, 0.5], [0, 0.5], 'k--', alpha=0.3)
        ax.set_xlabel('Reference MAF', fontsize=11)
        ax.set_ylabel('Imputed MAF', fontsize=11)
        ax.set_title('MAF Correlation', fontsize=12)
        cbar = plt.colorbar(ax.collections[0], ax=ax)
        cbar.set_label('R²', rotation=270, labelpad=15)
    else:
        # Alternative plot: R-squared vs number of variants
        maf_grouped = info_data.groupby('MAF_bin').agg({
            'Rsq': ['mean', 'count']
        }).reset_index()
        maf_grouped.columns = ['MAF_bin', 'mean_rsq', 'count']
        
        ax2 = ax.twinx()
        x = range(len(maf_grouped))
        ax.bar(x, maf_grouped['count'], alpha=0.5, color='blue', label='Count')
        ax2.plot(x, maf_grouped['mean_rsq'], 'ro-', label='Mean R²')
        ax.set_xticks(x)
        ax.set_xticklabels(maf_grouped['MAF_bin'], rotation=45, ha='right')
        ax.set_xlabel('MAF Bin', fontsize=11)
        ax.set_ylabel('Number of Variants', fontsize=11, color='blue')
        ax2.set_ylabel('Mean R²', fontsize=11, color='red')
        ax.set_title('Variants and Quality by MAF', fontsize=12)
    
    plt.suptitle('Ancestry-Aware Imputation Quality Assessment', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

def plot_ancestry_stratified_performance(info_data, ancestry_file, output_file):
    """
    Create plots showing imputation performance stratified by ancestry proportions
    This would require sample-level data with ancestry information
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Placeholder for ancestry-specific analysis
    # In practice, this would:
    # 1. Load ancestry proportions for each sample
    # 2. Stratify samples by major ancestry component
    # 3. Calculate imputation metrics per ancestry group
    # 4. Create comparative visualizations
    
    ax = axes[0]
    ax.text(0.5, 0.5, 'Ancestry-stratified analysis\nrequires sample-level data', 
            ha='center', va='center', fontsize=12)
    ax.set_title('Imputation Quality by Ancestry', fontsize=12)
    ax.axis('off')
    
    ax = axes[1]
    # Create a summary statistics table
    summary_stats = calculate_ancestry_specific_metrics(info_data)
    ax.axis('tight')
    ax.axis('off')
    table_data = []
    for metric in ['mean_r2', 'median_r2', 'prop_well_imputed']:
        table_data.append([metric, f"{summary_stats.loc['overall', metric]:.3f}"])
    
    table = ax.table(cellText=table_data, 
                    colLabels=['Metric', 'Value'],
                    cellLoc='center',
                    loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    ax.set_title('Overall Imputation Statistics', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

# Main execution when used as Nextflow template
if __name__ == "__main__":
    # Parse Nextflow template variables
    info_file = "${info_file}"
    ancestry_file = "${ancestry_file}" if "${ancestry_file}" != "" else None
    output_prefix = "${output_prefix}"
    
    # Read imputation info data
    info_data = read_file(info_file)
    
    # Convert numeric columns
    numeric_cols = ['Rsq', 'MAF', 'ALT_Frq', 'AvgCall']
    for col in numeric_cols:
        if col in info_data.columns:
            info_data[col] = pd.to_numeric(info_data[col], errors='coerce')
    
    # Generate plots
    plot_imputation_by_maf_ancestry(info_data, f"{output_prefix}_quality_assessment.pdf")
    
    if ancestry_file:
        plot_ancestry_stratified_performance(info_data, ancestry_file, 
                                            f"{output_prefix}_ancestry_stratified.pdf")
    
    # Generate summary report
    metrics = calculate_ancestry_specific_metrics(info_data)
    metrics.to_csv(f"{output_prefix}_metrics.tsv", sep='\t')
    
    print(f"Ancestry-aware imputation assessment complete. Files generated:")
    print(f"  - {output_prefix}_quality_assessment.pdf")
    if ancestry_file:
        print(f"  - {output_prefix}_ancestry_stratified.pdf")
    print(f"  - {output_prefix}_metrics.tsv")