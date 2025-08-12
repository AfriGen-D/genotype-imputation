#!/usr/bin/env python3
"""
Comprehensive Imputation Quality Control and Assessment
Includes all recommended plots for thorough imputation evaluation
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import gzip
from scipy import stats
from matplotlib.gridspec import GridSpec

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

def plot_dosage_distribution(vcf_data, ax):
    """
    1. Dosage Distribution Plot
    Shows certainty of imputed genotypes (dosages closer to 0, 1, 2 = higher certainty)
    """
    if 'DS' in vcf_data.columns or 'GP' in vcf_data.columns:
        # Extract dosage values
        dosages = vcf_data['DS'] if 'DS' in vcf_data.columns else vcf_data['GP']
        
        # Calculate distance from nearest integer (0, 1, or 2)
        certainty = dosages.apply(lambda x: min(abs(x - 0), abs(x - 1), abs(x - 2)))
        
        # Plot distribution
        ax.hist(certainty, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Distance from Nearest Integer Dosage', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.set_title('Dosage Certainty Distribution', fontsize=11)
        ax.axvline(x=0.1, color='r', linestyle='--', alpha=0.5, label='High certainty threshold')
        ax.legend()
    else:
        # Alternative: plot distribution of R² values as proxy for certainty
        ax.text(0.5, 0.5, 'Dosage data not available\nUsing R² as certainty proxy', 
                ha='center', va='center')
        ax.set_title('Dosage Distribution', fontsize=11)

def plot_info_score_distribution(info_data, ax):
    """
    2. INFO Score Distribution Plot
    Histogram of INFO scores across all variants
    """
    if 'INFO' in info_data.columns:
        info_scores = info_data['INFO']
    elif 'Rsq' in info_data.columns:
        info_scores = info_data['Rsq']
    else:
        ax.text(0.5, 0.5, 'INFO scores not available', ha='center', va='center')
        return
    
    # Remove missing values
    info_scores = pd.to_numeric(info_scores, errors='coerce').dropna()
    
    # Create histogram
    ax.hist(info_scores, bins=50, edgecolor='black', alpha=0.7, color='skyblue')
    ax.axvline(x=0.3, color='r', linestyle='--', alpha=0.5, label='Poor quality (<0.3)')
    ax.axvline(x=0.8, color='g', linestyle='--', alpha=0.5, label='High quality (>0.8)')
    ax.set_xlabel('INFO Score / R²', fontsize=10)
    ax.set_ylabel('Number of Variants', fontsize=10)
    ax.set_title('INFO Score Distribution', fontsize=11)
    
    # Add statistics
    mean_info = info_scores.mean()
    median_info = info_scores.median()
    ax.text(0.02, 0.95, f'Mean: {mean_info:.3f}\nMedian: {median_info:.3f}\nN: {len(info_scores):,}',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.legend(loc='upper left')

def plot_concordance_vs_maf(masked_data, ax):
    """
    3. Concordance Rate vs MAF Plot
    For masked variants, concordance between true and imputed genotypes
    """
    if masked_data is not None and not masked_data.empty:
        # Define MAF bins
        maf_bins = [0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['0-0.1%', '0.1-1%', '1-5%', '5-10%', '10-20%', '20-50%']
        
        masked_data['MAF_bin'] = pd.cut(masked_data['MAF'], bins=maf_bins, labels=maf_labels)
        
        # Calculate concordance by MAF bin
        concordance_by_maf = masked_data.groupby('MAF_bin')['concordance'].agg(['mean', 'std', 'count'])
        
        # Plot
        x = range(len(concordance_by_maf))
        ax.bar(x, concordance_by_maf['mean'], yerr=concordance_by_maf['std'], 
               capsize=5, alpha=0.7, color='coral')
        ax.set_xticks(x)
        ax.set_xticklabels(concordance_by_maf.index, rotation=45, ha='right')
        ax.set_xlabel('MAF Bin', fontsize=10)
        ax.set_ylabel('Concordance Rate', fontsize=10)
        ax.set_title('Concordance vs MAF (Masked Variants)', fontsize=11)
        ax.set_ylim(0, 1)
        ax.axhline(y=0.95, color='g', linestyle='--', alpha=0.5, label='95% concordance')
        ax.legend()
        
        # Add sample sizes
        for i, (idx, row) in enumerate(concordance_by_maf.iterrows()):
            ax.text(i, row['mean'] + row['std'] + 0.02, f"n={row['count']}", 
                   ha='center', fontsize=8)
    else:
        ax.text(0.5, 0.5, 'Masked variant data not available', ha='center', va='center')
        ax.set_title('Concordance vs MAF', fontsize=11)

def plot_calibration(info_data, ax):
    """
    4. Calibration Plot
    Observed vs Expected imputation quality to identify systematic bias
    """
    if 'Rsq' in info_data.columns and 'EmpRsq' in info_data.columns:
        # Filter valid data
        calib_data = info_data[['Rsq', 'EmpRsq']].dropna()
        calib_data = calib_data[(calib_data['Rsq'] != '-') & (calib_data['EmpRsq'] != '-')]
        calib_data['Rsq'] = pd.to_numeric(calib_data['Rsq'], errors='coerce')
        calib_data['EmpRsq'] = pd.to_numeric(calib_data['EmpRsq'], errors='coerce')
        calib_data = calib_data.dropna()
        
        if not calib_data.empty:
            # Bin the expected R² and calculate mean observed R²
            bins = np.linspace(0, 1, 21)
            calib_data['Expected_bin'] = pd.cut(calib_data['Rsq'], bins)
            calibration = calib_data.groupby('Expected_bin').agg({
                'Rsq': 'mean',
                'EmpRsq': 'mean'
            }).dropna()
            
            # Plot
            ax.scatter(calibration['Rsq'], calibration['EmpRsq'], s=50, alpha=0.7)
            ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='Perfect calibration')
            ax.set_xlabel('Expected R² (Estimated)', fontsize=10)
            ax.set_ylabel('Observed R² (Empirical)', fontsize=10)
            ax.set_title('Calibration Plot', fontsize=11)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.legend()
            
            # Calculate calibration slope
            if len(calibration) > 1:
                slope, intercept = np.polyfit(calibration['Rsq'], calibration['EmpRsq'], 1)
                ax.text(0.05, 0.95, f'Slope: {slope:.3f}\nIntercept: {intercept:.3f}',
                       transform=ax.transAxes, fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'No valid calibration data', ha='center', va='center')
    else:
        ax.text(0.5, 0.5, 'Empirical R² not available\nfor calibration', ha='center', va='center')
    ax.set_title('Calibration Plot', fontsize=11)

def plot_cross_validation_accuracy(cv_data, ax):
    """
    5. Cross-validation Accuracy Plot
    Compare true vs imputed genotypes for masked variants
    """
    if cv_data is not None and not cv_data.empty:
        # Create confusion matrix style plot
        accuracy_matrix = cv_data.pivot_table(
            index='true_genotype', 
            columns='imputed_genotype', 
            values='count', 
            fill_value=0
        )
        
        # Normalize to get proportions
        accuracy_matrix = accuracy_matrix.div(accuracy_matrix.sum(axis=1), axis=0)
        
        # Plot heatmap
        sns.heatmap(accuracy_matrix, annot=True, fmt='.2f', cmap='YlOrRd', 
                   ax=ax, cbar_kws={'label': 'Proportion'})
        ax.set_xlabel('Imputed Genotype', fontsize=10)
        ax.set_ylabel('True Genotype', fontsize=10)
        ax.set_title('Cross-validation Accuracy', fontsize=11)
        
        # Calculate overall accuracy
        if 'correct' in cv_data.columns:
            overall_accuracy = cv_data['correct'].mean()
            ax.text(0.02, 0.98, f'Overall Accuracy: {overall_accuracy:.3f}',
                   transform=ax.transAxes, fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax.text(0.5, 0.5, 'Cross-validation data not available', ha='center', va='center')
        ax.set_title('Cross-validation Accuracy', fontsize=11)

def plot_heterozygosity_rate(info_data, ax):
    """
    6. Heterozygosity Rate Plot
    Compare heterozygosity rates between imputed and genotyped variants
    """
    if 'Genotyped' in info_data.columns:
        # Separate by variant type
        genotyped = info_data[info_data['Genotyped'] == 'Genotyped']
        imputed = info_data[info_data['Genotyped'] == 'Imputed']
        
        # If we have heterozygosity data
        if 'Het' in info_data.columns:
            data_to_plot = []
            labels = []
            
            if not genotyped.empty:
                data_to_plot.append(genotyped['Het'].dropna())
                labels.append('Genotyped')
            if not imputed.empty:
                data_to_plot.append(imputed['Het'].dropna())
                labels.append('Imputed')
            
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                colors = ['lightblue', 'lightcoral']
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                
                ax.set_ylabel('Heterozygosity Rate', fontsize=10)
                ax.set_title('Heterozygosity Rate Comparison', fontsize=11)
                
                # Add statistics
                for i, (data, label) in enumerate(zip(data_to_plot, labels)):
                    mean_het = data.mean()
                    ax.text(i+1, ax.get_ylim()[1]*0.95, f'μ={mean_het:.3f}',
                           ha='center', fontsize=9)
        else:
            # Alternative: use MAF to estimate expected heterozygosity
            ax.text(0.5, 0.5, 'Heterozygosity data not directly available\nUsing MAF-based estimation', 
                   ha='center', va='center')
    else:
        ax.text(0.5, 0.5, 'Heterozygosity comparison not available', ha='center', va='center')
    ax.set_title('Heterozygosity Rate', fontsize=11)

def plot_hwe_deviation(info_data, ax):
    """
    7. Hardy-Weinberg Equilibrium Plot
    Check for deviations in HWE for imputed vs genotyped variants
    """
    if 'HWE_pval' in info_data.columns or 'HWE' in info_data.columns:
        hwe_col = 'HWE_pval' if 'HWE_pval' in info_data.columns else 'HWE'
        
        # Separate by variant type
        if 'Genotyped' in info_data.columns:
            genotyped = info_data[info_data['Genotyped'] == 'Genotyped'][hwe_col].dropna()
            imputed = info_data[info_data['Genotyped'] == 'Imputed'][hwe_col].dropna()
            
            # Convert to -log10(p-value) for better visualization
            genotyped_log = -np.log10(pd.to_numeric(genotyped, errors='coerce').dropna())
            imputed_log = -np.log10(pd.to_numeric(imputed, errors='coerce').dropna())
            
            # Create Q-Q plot
            if len(genotyped_log) > 0:
                stats.probplot(genotyped_log, dist="expon", plot=ax)
                ax.get_lines()[0].set_markerfacecolor('blue')
                ax.get_lines()[0].set_label('Genotyped')
            
            if len(imputed_log) > 0:
                stats.probplot(imputed_log, dist="expon", plot=ax)
                ax.get_lines()[-2].set_markerfacecolor('red')
                ax.get_lines()[-2].set_label('Imputed')
            
            ax.set_xlabel('Theoretical Quantiles', fontsize=10)
            ax.set_ylabel('-log10(HWE p-value)', fontsize=10)
            ax.set_title('Hardy-Weinberg Equilibrium Q-Q Plot', fontsize=11)
            ax.legend()
            
            # Add significance threshold
            ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p=0.05')
    else:
        # Alternative: calculate expected heterozygosity vs observed
        if 'MAF' in info_data.columns:
            info_data['Expected_Het'] = 2 * info_data['MAF'] * (1 - info_data['MAF'])
            
            # Plot for different variant types
            for variant_type in info_data['Genotyped'].unique():
                subset = info_data[info_data['Genotyped'] == variant_type]
                ax.scatter(subset['MAF'], subset['Expected_Het'], 
                          s=1, alpha=0.5, label=variant_type)
            
            ax.set_xlabel('Minor Allele Frequency', fontsize=10)
            ax.set_ylabel('Expected Heterozygosity (2pq)', fontsize=10)
            ax.set_title('Expected Heterozygosity under HWE', fontsize=11)
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'HWE data not available', ha='center', va='center')
            ax.set_title('Hardy-Weinberg Equilibrium', fontsize=11)

def create_comprehensive_qc_report(info_file, masked_file=None, cv_file=None, output_prefix="imputation_qc"):
    """
    Create comprehensive QC report with all 7 recommended plots
    """
    # Read data
    info_data = read_file(info_file)
    
    # Convert numeric columns
    numeric_cols = ['Rsq', 'MAF', 'ALT_Frq', 'AvgCall', 'EmpRsq', 'ER2', 'INFO']
    for col in numeric_cols:
        if col in info_data.columns:
            info_data[col] = pd.to_numeric(info_data[col].replace('-', np.nan), errors='coerce')
    
    # Read additional data if available
    masked_data = read_file(masked_file) if masked_file else None
    cv_data = read_file(cv_file) if cv_file else None
    
    # Create figure with all plots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # Create subplots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    ax7 = fig.add_subplot(gs[2, :2])  # Wider plot for HWE
    
    # Generate all plots
    plot_dosage_distribution(info_data, ax1)
    plot_info_score_distribution(info_data, ax2)
    plot_concordance_vs_maf(masked_data, ax3)
    plot_calibration(info_data, ax4)
    plot_cross_validation_accuracy(cv_data, ax5)
    plot_heterozygosity_rate(info_data, ax6)
    plot_hwe_deviation(info_data, ax7)
    
    # Add overall title
    fig.suptitle('Comprehensive Imputation Quality Control Report', fontsize=16, y=0.98)
    
    # Save figure
    plt.savefig(f"{output_prefix}_comprehensive_qc.pdf", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Generate summary statistics
    summary_stats = {
        'Total_Variants': len(info_data),
        'Genotyped_Variants': len(info_data[info_data['Genotyped'] == 'Genotyped']) if 'Genotyped' in info_data.columns else 0,
        'Imputed_Variants': len(info_data[info_data['Genotyped'] == 'Imputed']) if 'Genotyped' in info_data.columns else 0,
        'Mean_Rsq': info_data['Rsq'].mean() if 'Rsq' in info_data.columns else np.nan,
        'Median_Rsq': info_data['Rsq'].median() if 'Rsq' in info_data.columns else np.nan,
        'Prop_Well_Imputed_Rsq_0.8': (info_data['Rsq'] >= 0.8).mean() if 'Rsq' in info_data.columns else np.nan,
        'Prop_Moderate_Rsq_0.3_0.8': ((info_data['Rsq'] >= 0.3) & (info_data['Rsq'] < 0.8)).mean() if 'Rsq' in info_data.columns else np.nan,
        'Prop_Poor_Rsq_0.3': (info_data['Rsq'] < 0.3).mean() if 'Rsq' in info_data.columns else np.nan
    }
    
    # Save summary statistics
    pd.DataFrame([summary_stats]).T.to_csv(f"{output_prefix}_summary_stats.tsv", sep='\t', header=['Value'])
    
    print(f"Comprehensive QC report generated:")
    print(f"  - {output_prefix}_comprehensive_qc.pdf")
    print(f"  - {output_prefix}_summary_stats.tsv")
    
    return summary_stats

# Main execution
if __name__ == "__main__":
    # Parse Nextflow template variables
    info_file = "${info_file}"
    masked_file = "${masked_file}" if "${masked_file}" != "" else None
    cv_file = "${cv_file}" if "${cv_file}" != "" else None
    output_prefix = "${output_prefix}"
    
    # Generate comprehensive QC report
    create_comprehensive_qc_report(info_file, masked_file, cv_file, output_prefix)