#!/usr/bin/env python3
"""
Cross-validation Accuracy Plot
Using masked variants, compares true vs imputed genotypes with confusion matrix visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import gzip

def read_file(filename):
    """Read potentially gzipped files"""
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            return pd.read_csv(f, sep='\t')
    else:
        return pd.read_csv(filename, sep='\t')

# Parse arguments from Nextflow template
cv_file = "${cv_file}"  # Cross-validation results file
info_file = "${info_file}"  # Backup: use info file if CV not available
output_file = "${output_file}"

# Try to read cross-validation data
has_cv_data = False
cv_data = None

if cv_file and cv_file != "":
    try:
        cv_data = read_file(cv_file)
        if 'true_genotype' in cv_data.columns and 'imputed_genotype' in cv_data.columns:
            has_cv_data = True
    except:
        pass

# If no CV data, try to use info file for genotyped variants
if not has_cv_data and info_file and info_file != "":
    try:
        info_data = read_file(info_file)
        # For genotyped variants with empirical R², we can infer accuracy
        if 'Genotyped' in info_data.columns and 'EmpRsq' in info_data.columns:
            genotyped = info_data[info_data['Genotyped'] == 'Genotyped'].copy()
            genotyped['EmpRsq'] = pd.to_numeric(genotyped['EmpRsq'].replace('-', np.nan), errors='coerce')
            genotyped = genotyped.dropna(subset=['EmpRsq'])
            if len(genotyped) > 0:
                # Use EmpRsq as proxy for accuracy
                cv_data = genotyped
                has_cv_data = True
    except:
        pass

# Create figure
fig = plt.figure(figsize=(14, 10))

if has_cv_data and 'true_genotype' in cv_data.columns and 'imputed_genotype' in cv_data.columns:
    # Full cross-validation data available
    
    # Create confusion matrix
    genotype_labels = ['0/0', '0/1', '1/1']
    confusion_matrix = pd.crosstab(cv_data['true_genotype'], cv_data['imputed_genotype'], 
                                   normalize='index')
    
    # Ensure all genotype combinations are present
    for gt in genotype_labels:
        if gt not in confusion_matrix.index:
            confusion_matrix.loc[gt] = 0
        if gt not in confusion_matrix.columns:
            confusion_matrix[gt] = 0
    
    confusion_matrix = confusion_matrix.reindex(index=genotype_labels, columns=genotype_labels, fill_value=0)
    
    # Plot 1: Confusion matrix heatmap
    ax1 = plt.subplot(2, 3, 1)
    sns.heatmap(confusion_matrix, annot=True, fmt='.3f', cmap='YlOrRd', 
               square=True, cbar_kws={'label': 'Proportion'}, ax=ax1,
               vmin=0, vmax=1)
    ax1.set_xlabel('Imputed Genotype', fontsize=11)
    ax1.set_ylabel('True Genotype', fontsize=11)
    ax1.set_title('Genotype Confusion Matrix', fontsize=12)
    
    # Calculate accuracy metrics
    total_correct = np.diag(confusion_matrix.values).sum()
    overall_accuracy = total_correct / confusion_matrix.values.sum() if confusion_matrix.values.sum() > 0 else 0
    
    # Per-genotype accuracy
    per_genotype_acc = np.diag(confusion_matrix.values)
    
    # Plot 2: Per-genotype accuracy
    ax2 = plt.subplot(2, 3, 2)
    bars = ax2.bar(range(len(genotype_labels)), per_genotype_acc, 
                   color=['#2ca02c', '#ff7f0e', '#d62728'], edgecolor='black', alpha=0.7)
    ax2.set_xticks(range(len(genotype_labels)))
    ax2.set_xticklabels(genotype_labels)
    ax2.set_xlabel('True Genotype', fontsize=11)
    ax2.set_ylabel('Accuracy', fontsize=11)
    ax2.set_title('Per-Genotype Accuracy', fontsize=12)
    ax2.set_ylim(0, 1.05)
    ax2.axhline(y=overall_accuracy, color='blue', linestyle='--', alpha=0.5, 
               label=f'Overall: {overall_accuracy:.3f}')
    ax2.legend()
    
    # Add value labels
    for bar, val in zip(bars, per_genotype_acc):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10)
    
    # Plot 3: Accuracy by MAF (if available)
    ax3 = plt.subplot(2, 3, 3)
    if 'MAF' in cv_data.columns:
        cv_data['MAF'] = pd.to_numeric(cv_data['MAF'], errors='coerce')
        cv_data['correct'] = cv_data['true_genotype'] == cv_data['imputed_genotype']
        
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
        maf_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%']
        cv_data['MAF_bin'] = pd.cut(cv_data['MAF'], bins=maf_bins, labels=maf_labels)
        
        accuracy_by_maf = cv_data.groupby('MAF_bin')['correct'].mean()
        counts_by_maf = cv_data.groupby('MAF_bin')['correct'].count()
        
        bars = ax3.bar(range(len(accuracy_by_maf)), accuracy_by_maf.values,
                       color=plt.cm.viridis(np.linspace(0.3, 0.9, len(accuracy_by_maf))),
                       edgecolor='black', alpha=0.7)
        ax3.set_xticks(range(len(accuracy_by_maf)))
        ax3.set_xticklabels(accuracy_by_maf.index, rotation=45, ha='right')
        ax3.set_xlabel('MAF Bin', fontsize=11)
        ax3.set_ylabel('Accuracy', fontsize=11)
        ax3.set_title('Accuracy by MAF', fontsize=12)
        ax3.set_ylim(0, 1.05)
        
        # Add sample sizes
        for bar, val, count in zip(bars, accuracy_by_maf.values, counts_by_maf.values):
            ax3.text(bar.get_x() + bar.get_width()/2., val + 0.01,
                    f'{val:.3f}\n(n={count})', ha='center', va='bottom', fontsize=8)
    else:
        ax3.text(0.5, 0.5, 'MAF data not available', ha='center', va='center')
    
    # Plot 4: Error types distribution
    ax4 = plt.subplot(2, 3, 4)
    error_types = []
    error_counts = []
    
    # Calculate different error types
    for true_gt in genotype_labels:
        for imp_gt in genotype_labels:
            if true_gt != imp_gt:
                count = cv_data[(cv_data['true_genotype'] == true_gt) & 
                              (cv_data['imputed_genotype'] == imp_gt)].shape[0]
                if count > 0:
                    error_types.append(f'{true_gt}→{imp_gt}')
                    error_counts.append(count)
    
    if error_types:
        # Sort by frequency
        sorted_idx = np.argsort(error_counts)[::-1]
        error_types = [error_types[i] for i in sorted_idx]
        error_counts = [error_counts[i] for i in sorted_idx]
        
        # Show top errors
        top_n = min(10, len(error_types))
        ax4.barh(range(top_n), error_counts[:top_n], color='coral', edgecolor='black', alpha=0.7)
        ax4.set_yticks(range(top_n))
        ax4.set_yticklabels(error_types[:top_n])
        ax4.set_xlabel('Number of Errors', fontsize=11)
        ax4.set_title('Most Common Error Types', fontsize=12)
        ax4.invert_yaxis()
    else:
        ax4.text(0.5, 0.5, 'No errors found', ha='center', va='center')
    
elif has_cv_data and 'EmpRsq' in cv_data.columns:
    # Using empirical R² as accuracy proxy
    
    ax1 = plt.subplot(2, 2, 1)
    emp_rsq = pd.to_numeric(cv_data['EmpRsq'], errors='coerce').dropna()
    
    # Distribution of empirical R²
    ax1.hist(emp_rsq, bins=30, edgecolor='black', alpha=0.7, color='skyblue')
    ax1.axvline(x=emp_rsq.mean(), color='r', linestyle='--', linewidth=2, 
               label=f'Mean: {emp_rsq.mean():.3f}')
    ax1.set_xlabel('Empirical R²', fontsize=11)
    ax1.set_ylabel('Number of Variants', fontsize=11)
    ax1.set_title('Cross-validation Accuracy Distribution', fontsize=12)
    ax1.legend()
    
    # Accuracy categories
    ax2 = plt.subplot(2, 2, 2)
    categories = ['Poor\n(<0.3)', 'Moderate\n(0.3-0.8)', 'High\n(≥0.8)']
    counts = [
        (emp_rsq < 0.3).sum(),
        ((emp_rsq >= 0.3) & (emp_rsq < 0.8)).sum(),
        (emp_rsq >= 0.8).sum()
    ]
    colors = ['#d62728', '#ff7f0e', '#2ca02c']
    
    bars = ax2.bar(range(len(categories)), counts, color=colors, edgecolor='black', alpha=0.7)
    ax2.set_xticks(range(len(categories)))
    ax2.set_xticklabels(categories)
    ax2.set_ylabel('Number of Variants', fontsize=11)
    ax2.set_title('Accuracy Categories', fontsize=12)
    
    # Add percentages
    total = sum(counts)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}\n({100*count/total:.1f}%)',
                ha='center', va='bottom', fontsize=10)
    
    # Summary statistics
    ax3 = plt.subplot(2, 2, 3)
    ax3.axis('off')
    
    stats_text = f'Cross-validation Statistics\n\n'
    stats_text += f'Total variants: {len(emp_rsq):,}\n'
    stats_text += f'Mean accuracy (R²): {emp_rsq.mean():.3f}\n'
    stats_text += f'Median accuracy: {emp_rsq.median():.3f}\n'
    stats_text += f'Std deviation: {emp_rsq.std():.3f}\n'
    stats_text += f'Min accuracy: {emp_rsq.min():.3f}\n'
    stats_text += f'Max accuracy: {emp_rsq.max():.3f}\n'
    stats_text += f'\nHigh quality (≥0.8): {(emp_rsq >= 0.8).sum():,} ({100*(emp_rsq >= 0.8).mean():.1f}%)\n'
    stats_text += f'Poor quality (<0.3): {(emp_rsq < 0.3).sum():,} ({100*(emp_rsq < 0.3).mean():.1f}%)'
    
    ax3.text(0.1, 0.9, stats_text, transform=ax3.transAxes, fontsize=11,
            verticalalignment='top', family='monospace')

else:
    # No cross-validation data available
    ax = plt.subplot(1, 1, 1)
    ax.text(0.5, 0.5, 
           'Cross-validation data not available\n\n' +
           'To generate cross-validation metrics:\n\n' +
           '1. Select a subset of known genotypes\n' +
           '2. Mask these genotypes before imputation\n' +
           '3. Impute the masked genotypes\n' +
           '4. Compare imputed vs true genotypes\n' +
           '5. Calculate accuracy metrics\n\n' +
           'This provides unbiased accuracy estimates',
           ha='center', va='center', fontsize=11,
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_title('Cross-validation Analysis', fontsize=13)
    ax.axis('off')

plt.suptitle('Cross-validation Accuracy Assessment', fontsize=14, y=0.98)
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"Cross-validation plot saved to {output_file}")