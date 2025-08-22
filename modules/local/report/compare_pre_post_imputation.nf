process COMPARE_PRE_POST_IMPUTATION {
    tag "$meta.id"
    label 'process_medium'
    label 'python_plotting'
    
    publishDir "/scratch3/users/mamana/results/reports/pre_post_comparison/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), path(pre_vcf), path(pre_index), path(post_vcf), path(post_index)
    
    output:
    tuple val(meta), path("*.comparison_stats.txt")     , emit: stats
    tuple val(meta), path("*.variant_counts.png")       , emit: count_plot
    tuple val(meta), path("*.maf_distribution.png")     , emit: maf_plot
    tuple val(meta), path("*.coverage_improvement.png") , emit: coverage_plot
    tuple val(meta), path("*.variant_gain.png")         , emit: gain_plot
    tuple val(meta), path("*.comparison_report.html")   , emit: report
    path "versions.yml"                                  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create comparison analysis script
    cat <<'PYEOF' > compare_imputation.py
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
import json
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def extract_vcf_stats(vcf_path, label):
    """Extract statistics from VCF file"""
    stats = {
        'label': label,
        'total_variants': 0,
        'snps': 0,
        'indels': 0,
        'multiallelic': 0,
        'maf_bins': defaultdict(int),
        'chr_counts': defaultdict(int),
        'samples': 0,
        'mafs': []
    }
    
    vcf = VCF(vcf_path)
    stats['samples'] = len(vcf.samples)
    
    for variant in vcf:
        stats['total_variants'] += 1
        stats['chr_counts'][variant.CHROM] += 1
        
        # Variant type
        if variant.is_snp:
            stats['snps'] += 1
        elif variant.is_indel:
            stats['indels'] += 1
            
        # Check multiallelic
        if len(variant.ALT) > 1:
            stats['multiallelic'] += 1
        
        # Calculate MAF
        if variant.num_called > 0:
            af = variant.aaf  # alternate allele frequency
            if af is not None and not np.isnan(af):
                maf = min(af, 1 - af)
                stats['mafs'].append(maf)
                
                # Bin MAF
                if maf == 0:
                    stats['maf_bins']['Monomorphic'] += 1
                elif maf < 0.001:
                    stats['maf_bins']['0 < MAF < 0.001'] += 1
                elif maf < 0.005:
                    stats['maf_bins']['0.001 ≤ MAF < 0.005'] += 1
                elif maf < 0.01:
                    stats['maf_bins']['0.005 ≤ MAF < 0.01'] += 1
                elif maf < 0.05:
                    stats['maf_bins']['0.01 ≤ MAF < 0.05'] += 1
                elif maf < 0.1:
                    stats['maf_bins']['0.05 ≤ MAF < 0.1'] += 1
                elif maf < 0.2:
                    stats['maf_bins']['0.1 ≤ MAF < 0.2'] += 1
                else:
                    stats['maf_bins']['MAF ≥ 0.2'] += 1
    
    vcf.close()
    return stats

def plot_variant_counts(pre_stats, post_stats, output_prefix):
    """Create variant count comparison plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Total variant counts
    ax = axes[0, 0]
    categories = ['Total', 'SNPs', 'INDELs', 'Multiallelic']
    pre_counts = [pre_stats['total_variants'], pre_stats['snps'], 
                  pre_stats['indels'], pre_stats['multiallelic']]
    post_counts = [post_stats['total_variants'], post_stats['snps'], 
                   post_stats['indels'], post_stats['multiallelic']]
    
    x = np.arange(len(categories))
    width = 0.35
    ax.bar(x - width/2, pre_counts, width, label='Pre-imputation', color='lightcoral')
    ax.bar(x + width/2, post_counts, width, label='Post-imputation', color='skyblue')
    ax.set_xlabel('Variant Type')
    ax.set_ylabel('Count')
    ax.set_title('Variant Counts: Pre vs Post Imputation')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add percentage increase labels
    for i, (pre, post) in enumerate(zip(pre_counts, post_counts)):
        if pre > 0:
            pct_increase = ((post - pre) / pre) * 100
            ax.text(i, max(pre, post) * 1.05, f'+{pct_increase:.1f}%', 
                   ha='center', fontsize=9, color='green')
    
    # 2. Chromosome distribution
    ax = axes[0, 1]
    chroms = sorted(set(list(pre_stats['chr_counts'].keys()) + 
                       list(post_stats['chr_counts'].keys())),
                   key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').zfill(2)))
    
    pre_chr = [pre_stats['chr_counts'].get(chr, 0) for chr in chroms]
    post_chr = [post_stats['chr_counts'].get(chr, 0) for chr in chroms]
    
    x = np.arange(len(chroms))
    ax.plot(x, pre_chr, 'o-', label='Pre-imputation', color='lightcoral', alpha=0.7)
    ax.plot(x, post_chr, 's-', label='Post-imputation', color='skyblue', alpha=0.7)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Variant Count')
    ax.set_title('Variant Distribution by Chromosome')
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Variant gain by chromosome
    ax = axes[1, 0]
    gain = [post - pre for pre, post in zip(pre_chr, post_chr)]
    colors = ['green' if g >= 0 else 'red' for g in gain]
    bars = ax.bar(x, gain, color=colors, alpha=0.7)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Variant Gain')
    ax.set_title('Imputation Gain by Chromosome')
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha='right')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # 4. Summary statistics table
    ax = axes[1, 1]
    ax.axis('tight')
    ax.axis('off')
    
    summary_data = [
        ['Metric', 'Pre-imputation', 'Post-imputation', 'Change'],
        ['Total Variants', f'{pre_stats["total_variants"]:,}', 
         f'{post_stats["total_variants"]:,}',
         f'+{post_stats["total_variants"] - pre_stats["total_variants"]:,}'],
        ['SNPs', f'{pre_stats["snps"]:,}', 
         f'{post_stats["snps"]:,}',
         f'+{post_stats["snps"] - pre_stats["snps"]:,}'],
        ['INDELs', f'{pre_stats["indels"]:,}', 
         f'{post_stats["indels"]:,}',
         f'+{post_stats["indels"] - pre_stats["indels"]:,}'],
        ['Samples', f'{pre_stats["samples"]}', 
         f'{post_stats["samples"]}', '-']
    ]
    
    table = ax.table(cellText=summary_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    # Style header row
    for i in range(4):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.suptitle(f'Pre/Post-Imputation Comparison: {output_prefix}', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.variant_counts.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_maf_distribution(pre_stats, post_stats, output_prefix):
    """Plot MAF distribution comparison"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. MAF histogram
    ax = axes[0, 0]
    if pre_stats['mafs'] and post_stats['mafs']:
        ax.hist(pre_stats['mafs'], bins=50, alpha=0.5, label='Pre-imputation', 
                color='lightcoral', density=True)
        ax.hist(post_stats['mafs'], bins=50, alpha=0.5, label='Post-imputation', 
                color='skyblue', density=True)
        ax.set_xlabel('Minor Allele Frequency')
        ax.set_ylabel('Density')
        ax.set_title('MAF Distribution')
        ax.legend()
        ax.set_xlim(0, 0.5)
        ax.grid(True, alpha=0.3)
    
    # 2. MAF bins comparison
    ax = axes[0, 1]
    maf_bins = ['Monomorphic', '0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', 
                '0.005 ≤ MAF < 0.01', '0.01 ≤ MAF < 0.05', '0.05 ≤ MAF < 0.1',
                '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']
    
    pre_bin_counts = [pre_stats['maf_bins'].get(bin, 0) for bin in maf_bins]
    post_bin_counts = [post_stats['maf_bins'].get(bin, 0) for bin in maf_bins]
    
    x = np.arange(len(maf_bins))
    width = 0.35
    ax.bar(x - width/2, pre_bin_counts, width, label='Pre-imputation', color='lightcoral')
    ax.bar(x + width/2, post_bin_counts, width, label='Post-imputation', color='skyblue')
    ax.set_xlabel('MAF Bins')
    ax.set_ylabel('Variant Count')
    ax.set_title('Variant Counts by MAF Bins')
    ax.set_xticks(x)
    ax.set_xticklabels(maf_bins, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Cumulative MAF distribution
    ax = axes[1, 0]
    if pre_stats['mafs'] and post_stats['mafs']:
        pre_sorted = np.sort(pre_stats['mafs'])
        post_sorted = np.sort(post_stats['mafs'])
        pre_cumulative = np.arange(1, len(pre_sorted) + 1) / len(pre_sorted)
        post_cumulative = np.arange(1, len(post_sorted) + 1) / len(post_sorted)
        
        ax.plot(pre_sorted, pre_cumulative, label='Pre-imputation', 
                color='lightcoral', alpha=0.7)
        ax.plot(post_sorted, post_cumulative, label='Post-imputation', 
                color='skyblue', alpha=0.7)
        ax.set_xlabel('Minor Allele Frequency')
        ax.set_ylabel('Cumulative Proportion')
        ax.set_title('Cumulative MAF Distribution')
        ax.legend()
        ax.set_xlim(0, 0.5)
        ax.grid(True, alpha=0.3)
    
    # 4. Rare variant enrichment
    ax = axes[1, 1]
    rare_categories = ['Rare\\\\n(MAF < 0.01)', 'Low Frequency\\\\n(0.01 ≤ MAF < 0.05)', 
                       'Common\\\\n(MAF ≥ 0.05)']
    
    pre_rare = sum(pre_bin_counts[1:4])  # MAF < 0.01
    pre_low = pre_bin_counts[4]  # 0.01 ≤ MAF < 0.05
    pre_common = sum(pre_bin_counts[5:])  # MAF ≥ 0.05
    
    post_rare = sum(post_bin_counts[1:4])
    post_low = post_bin_counts[4]
    post_common = sum(post_bin_counts[5:])
    
    pre_counts = [pre_rare, pre_low, pre_common]
    post_counts = [post_rare, post_low, post_common]
    
    x = np.arange(len(rare_categories))
    width = 0.35
    bars1 = ax.bar(x - width/2, pre_counts, width, label='Pre-imputation', color='lightcoral')
    bars2 = ax.bar(x + width/2, post_counts, width, label='Post-imputation', color='skyblue')
    
    ax.set_ylabel('Variant Count')
    ax.set_title('Variant Frequency Categories')
    ax.set_xticks(x)
    ax.set_xticklabels(rare_categories)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add fold change annotations
    for i, (pre, post) in enumerate(zip(pre_counts, post_counts)):
        if pre > 0:
            fold_change = post / pre
            ax.text(i, max(pre, post) * 1.05, f'{fold_change:.1f}x', 
                   ha='center', fontsize=9, color='darkgreen', weight='bold')
    
    plt.suptitle(f'MAF Distribution Analysis: {output_prefix}', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.maf_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_coverage_improvement(pre_stats, post_stats, output_prefix):
    """Plot coverage improvement metrics"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Variant density (variants per Mb)
    ax = axes[0, 0]
    # Assuming average chromosome size for density calculation
    chr_sizes_mb = {'chr1': 249, 'chr2': 242, 'chr3': 198, 'chr4': 190, 'chr5': 182,
                    'chr6': 171, 'chr7': 159, 'chr8': 145, 'chr9': 138, 'chr10': 134,
                    'chr11': 135, 'chr12': 133, 'chr13': 114, 'chr14': 107, 'chr15': 102,
                    'chr16': 90, 'chr17': 83, 'chr18': 80, 'chr19': 59, 'chr20': 64,
                    'chr21': 47, 'chr22': 51, 'chrX': 155, 'chrY': 59}
    
    chroms = sorted(set(list(pre_stats['chr_counts'].keys()) + 
                       list(post_stats['chr_counts'].keys())),
                   key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').zfill(2)))
    
    pre_density = []
    post_density = []
    valid_chroms = []
    
    for chr in chroms:
        if chr in chr_sizes_mb:
            pre_count = pre_stats['chr_counts'].get(chr, 0)
            post_count = post_stats['chr_counts'].get(chr, 0)
            size_mb = chr_sizes_mb[chr]
            pre_density.append(pre_count / size_mb)
            post_density.append(post_count / size_mb)
            valid_chroms.append(chr)
    
    if valid_chroms:
        x = np.arange(len(valid_chroms))
        width = 0.35
        ax.bar(x - width/2, pre_density, width, label='Pre-imputation', 
               color='lightcoral', alpha=0.7)
        ax.bar(x + width/2, post_density, width, label='Post-imputation', 
               color='skyblue', alpha=0.7)
        ax.set_xlabel('Chromosome')
        ax.set_ylabel('Variants per Mb')
        ax.set_title('Variant Density by Chromosome')
        ax.set_xticks(x)
        ax.set_xticklabels(valid_chroms, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # 2. Imputation gain ratio
    ax = axes[0, 1]
    categories = ['Total Variants', 'Rare\\\\n(MAF < 0.01)', 'Low Freq\\\\n(0.01-0.05)', 
                  'Common\\\\n(MAF > 0.05)']
    
    # Calculate gains
    total_gain = (post_stats['total_variants'] / max(pre_stats['total_variants'], 1) - 1) * 100
    
    pre_rare = sum([pre_stats['maf_bins'].get(k, 0) for k in 
                    ['0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', '0.005 ≤ MAF < 0.01']])
    post_rare = sum([post_stats['maf_bins'].get(k, 0) for k in 
                     ['0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', '0.005 ≤ MAF < 0.01']])
    rare_gain = (post_rare / max(pre_rare, 1) - 1) * 100 if pre_rare > 0 else 0
    
    pre_low = pre_stats['maf_bins'].get('0.01 ≤ MAF < 0.05', 0)
    post_low = post_stats['maf_bins'].get('0.01 ≤ MAF < 0.05', 0)
    low_gain = (post_low / max(pre_low, 1) - 1) * 100 if pre_low > 0 else 0
    
    pre_common = sum([pre_stats['maf_bins'].get(k, 0) for k in 
                      ['0.05 ≤ MAF < 0.1', '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']])
    post_common = sum([post_stats['maf_bins'].get(k, 0) for k in 
                       ['0.05 ≤ MAF < 0.1', '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']])
    common_gain = (post_common / max(pre_common, 1) - 1) * 100 if pre_common > 0 else 0
    
    gains = [total_gain, rare_gain, low_gain, common_gain]
    colors = ['green' if g > 0 else 'red' for g in gains]
    
    bars = ax.bar(categories, gains, color=colors, alpha=0.7)
    ax.set_ylabel('Percentage Gain (%)')
    ax.set_title('Imputation Gain by Variant Category')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, gain in zip(bars, gains):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + (5 if height > 0 else -15),
                f'{gain:.1f}%', ha='center', va='bottom' if height > 0 else 'top',
                fontsize=10, weight='bold')
    
    # 3. Variant type proportions
    ax = axes[1, 0]
    labels = ['SNPs', 'INDELs', 'Multiallelic']
    pre_props = [pre_stats['snps'], pre_stats['indels'], pre_stats['multiallelic']]
    post_props = [post_stats['snps'], post_stats['indels'], post_stats['multiallelic']]
    
    x = np.arange(len(labels))
    width = 0.35
    
    # Calculate percentages
    pre_total = sum(pre_props)
    post_total = sum(post_props)
    pre_pct = [p/pre_total*100 if pre_total > 0 else 0 for p in pre_props]
    post_pct = [p/post_total*100 if post_total > 0 else 0 for p in post_props]
    
    bars1 = ax.bar(x - width/2, pre_pct, width, label='Pre-imputation', color='lightcoral')
    bars2 = ax.bar(x + width/2, post_pct, width, label='Post-imputation', color='skyblue')
    
    ax.set_ylabel('Percentage (%)')
    ax.set_title('Variant Type Proportions')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add percentage labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # 4. Summary metrics
    ax = axes[1, 1]
    ax.axis('tight')
    ax.axis('off')
    
    # Calculate key metrics
    total_increase = post_stats['total_variants'] - pre_stats['total_variants']
    pct_increase = (total_increase / pre_stats['total_variants'] * 100) if pre_stats['total_variants'] > 0 else 0
    fold_increase = post_stats['total_variants'] / max(pre_stats['total_variants'], 1)
    
    metrics_data = [
        ['Imputation Summary', 'Value'],
        ['Variants Added', f'{total_increase:,}'],
        ['Percentage Increase', f'{pct_increase:.1f}%'],
        ['Fold Increase', f'{fold_increase:.2f}x'],
        ['Final Variant Count', f'{post_stats["total_variants"]:,}'],
        ['Final SNP Count', f'{post_stats["snps"]:,}'],
        ['Final INDEL Count', f'{post_stats["indels"]:,}']
    ]
    
    table = ax.table(cellText=metrics_data, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.3, 2)
    
    # Style header row
    for i in range(2):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Highlight key metrics
    table[(2, 1)].set_text_props(weight='bold', color='green')
    table[(3, 1)].set_text_props(weight='bold', color='darkblue')
    
    plt.suptitle(f'Coverage Improvement Analysis: {output_prefix}', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.coverage_improvement.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_variant_gain(pre_stats, post_stats, output_prefix):
    """Create detailed variant gain visualization"""
    fig = plt.figure(figsize=(16, 10))
    
    # Create grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Waterfall chart of gains
    ax1 = fig.add_subplot(gs[0, :])
    categories = ['Starting', 'Rare', 'Low Freq', 'Common', 'Final']
    values = [pre_stats['total_variants']]
    
    # Calculate incremental gains
    pre_rare = sum([pre_stats['maf_bins'].get(k, 0) for k in 
                    ['0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', '0.005 ≤ MAF < 0.01']])
    post_rare = sum([post_stats['maf_bins'].get(k, 0) for k in 
                     ['0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', '0.005 ≤ MAF < 0.01']])
    rare_gain = post_rare - pre_rare
    
    pre_low = pre_stats['maf_bins'].get('0.01 ≤ MAF < 0.05', 0)
    post_low = post_stats['maf_bins'].get('0.01 ≤ MAF < 0.05', 0)
    low_gain = post_low - pre_low
    
    pre_common = sum([pre_stats['maf_bins'].get(k, 0) for k in 
                      ['0.05 ≤ MAF < 0.1', '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']])
    post_common = sum([post_stats['maf_bins'].get(k, 0) for k in 
                       ['0.05 ≤ MAF < 0.1', '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']])
    common_gain = post_common - pre_common
    
    # Cumulative values for waterfall
    cumulative = pre_stats['total_variants']
    x_pos = [0]
    
    # Starting bar
    ax1.bar(0, pre_stats['total_variants'], color='lightcoral', alpha=0.7, label='Pre-imputation')
    
    # Gain bars
    gains = [rare_gain, low_gain, common_gain]
    labels = ['Rare\\\\nGain', 'Low Freq\\\\nGain', 'Common\\\\nGain']
    colors = ['#2ecc71', '#3498db', '#9b59b6']
    
    for i, (gain, label, color) in enumerate(zip(gains, labels, colors), 1):
        ax1.bar(i, gain, bottom=cumulative, color=color, alpha=0.7, label=label)
        if gain != 0:
            ax1.text(i, cumulative + gain/2, f'+{gain:,}', ha='center', va='center', 
                    fontsize=9, color='white', weight='bold')
        cumulative += gain
    
    # Final bar
    ax1.bar(4, post_stats['total_variants'], color='skyblue', alpha=0.7, label='Post-imputation')
    
    ax1.set_xticks(range(5))
    ax1.set_xticklabels(['Pre-\\\\nimputation', 'Rare\\\\nGain', 'Low Freq\\\\nGain', 
                         'Common\\\\nGain', 'Post-\\\\nimputation'])
    ax1.set_ylabel('Variant Count')
    ax1.set_title('Waterfall Chart: Variant Gains by Frequency Category', fontsize=12, weight='bold')
    ax1.legend(loc='upper left', ncol=5)
    ax1.grid(True, alpha=0.3)
    
    # 2. Pie chart of gain composition
    ax2 = fig.add_subplot(gs[1, 0])
    if sum(gains) > 0:
        sizes = [g for g in gains if g > 0]
        labels_pie = [l for l, g in zip(['Rare', 'Low Freq', 'Common'], gains) if g > 0]
        colors_pie = [c for c, g in zip(colors, gains) if g > 0]
        
        wedges, texts, autotexts = ax2.pie(sizes, labels=labels_pie, colors=colors_pie,
                                            autopct='%1.1f%%', startangle=90)
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_weight('bold')
        ax2.set_title('Gain Composition\\nby Frequency')
    
    # 3. Chromosome heatmap
    ax3 = fig.add_subplot(gs[1, 1:])
    chroms = sorted(set(list(pre_stats['chr_counts'].keys()) + 
                       list(post_stats['chr_counts'].keys())),
                   key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').zfill(2)))
    
    gain_matrix = []
    for chr in chroms:
        pre_count = pre_stats['chr_counts'].get(chr, 0)
        post_count = post_stats['chr_counts'].get(chr, 0)
        pct_gain = ((post_count - pre_count) / max(pre_count, 1)) * 100 if pre_count > 0 else 0
        gain_matrix.append([pct_gain])
    
    im = ax3.imshow(np.array(gain_matrix).T, cmap='RdYlGn', aspect='auto', vmin=-50, vmax=200)
    ax3.set_xticks(range(len(chroms)))
    ax3.set_xticklabels(chroms, rotation=45, ha='right')
    ax3.set_yticks([0])
    ax3.set_yticklabels(['% Gain'])
    ax3.set_title('Percentage Gain by Chromosome')
    
    # Add text annotations
    for i, chr in enumerate(chroms):
        ax3.text(i, 0, f'{gain_matrix[i][0]:.0f}%', ha='center', va='center',
                color='white' if abs(gain_matrix[i][0]) > 100 else 'black', fontsize=8)
    
    plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
    
    # 4. Before/After comparison
    ax4 = fig.add_subplot(gs[2, :2])
    comparison_data = {
        'Metric': ['Total Variants', 'SNPs', 'INDELs', 'Rare (MAF<0.01)', 
                   'Low Freq (0.01-0.05)', 'Common (MAF>0.05)'],
        'Pre-imputation': [
            pre_stats['total_variants'],
            pre_stats['snps'],
            pre_stats['indels'],
            pre_rare,
            pre_low,
            pre_common
        ],
        'Post-imputation': [
            post_stats['total_variants'],
            post_stats['snps'],
            post_stats['indels'],
            post_rare,
            post_low,
            post_common
        ]
    }
    
    df = pd.DataFrame(comparison_data)
    df['Gain'] = df['Post-imputation'] - df['Pre-imputation']
    df['% Increase'] = (df['Gain'] / df['Pre-imputation'].replace(0, 1)) * 100
    
    # Format for display
    df_display = df.copy()
    for col in ['Pre-imputation', 'Post-imputation', 'Gain']:
        df_display[col] = df_display[col].apply(lambda x: f'{int(x):,}')
    df_display['% Increase'] = df_display['% Increase'].apply(lambda x: f'{x:.1f}%')
    
    # Create table
    ax4.axis('tight')
    ax4.axis('off')
    table = ax4.table(cellText=df_display.values, colLabels=df_display.columns,
                     cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Style header
    for i in range(len(df_display.columns)):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Color code gain column
    for i in range(1, len(df_display) + 1):
        gain_val = df.iloc[i-1]['Gain']
        if gain_val > 0:
            table[(i, 3)].set_text_props(color='green', weight='bold')
            table[(i, 4)].set_text_props(color='green', weight='bold')
    
    ax4.set_title('Detailed Variant Gain Summary', fontsize=11, weight='bold', pad=20)
    
    # 5. Key metrics box
    ax5 = fig.add_subplot(gs[2, 2])
    ax5.axis('off')
    
    total_gain = post_stats['total_variants'] - pre_stats['total_variants']
    fold_change = post_stats['total_variants'] / max(pre_stats['total_variants'], 1)
    
    metrics_text = f"""
    KEY METRICS
    
    Total Gain:
    {total_gain:,} variants
    
    Fold Change:
    {fold_change:.2f}x
    
    Largest Gain:
    {'Rare variants' if rare_gain >= max(low_gain, common_gain) else 'Low freq variants' if low_gain >= common_gain else 'Common variants'}
    
    Average per Chr:
    {total_gain // len(chroms):,} variants
    """
    
    ax5.text(0.5, 0.5, metrics_text, ha='center', va='center', fontsize=11,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.3))
    
    plt.suptitle(f'Comprehensive Variant Gain Analysis: {output_prefix}', 
                fontsize=14, weight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}.variant_gain.png', dpi=150, bbox_inches='tight')
    plt.close()

def generate_html_report(pre_stats, post_stats, output_prefix):
    """Generate comprehensive HTML report"""
    
    total_gain = post_stats['total_variants'] - pre_stats['total_variants']
    pct_increase = (total_gain / pre_stats['total_variants'] * 100) if pre_stats['total_variants'] > 0 else 0
    fold_change = post_stats['total_variants'] / max(pre_stats['total_variants'], 1)
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Pre/Post-Imputation Comparison Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .header {{
                background-color: #2c3e50;
                color: white;
                padding: 20px;
                border-radius: 5px;
                margin-bottom: 20px;
            }}
            .summary-box {{
                background-color: white;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin-bottom: 20px;
            }}
            .metric-card {{
                display: inline-block;
                background-color: #ecf0f1;
                padding: 15px;
                border-radius: 5px;
                margin: 10px;
                min-width: 200px;
                text-align: center;
            }}
            .metric-value {{
                font-size: 24px;
                font-weight: bold;
                color: #2c3e50;
            }}
            .metric-label {{
                font-size: 14px;
                color: #7f8c8d;
                margin-top: 5px;
            }}
            .increase {{
                color: #27ae60;
            }}
            .decrease {{
                color: #e74c3c;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            th {{
                background-color: #34495e;
                color: white;
                padding: 10px;
                text-align: left;
            }}
            td {{
                padding: 10px;
                border-bottom: 1px solid #ecf0f1;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .plot-container {{
                background-color: white;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin-bottom: 20px;
                text-align: center;
            }}
            .plot-title {{
                font-size: 18px;
                font-weight: bold;
                margin-bottom: 10px;
                color: #2c3e50;
            }}
            img {{
                max-width: 100%;
                height: auto;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Pre/Post-Imputation Comparison Report</h1>
            <p>Sample: {output_prefix}</p>
            <p>Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="summary-box">
            <h2>Executive Summary</h2>
            <div>
                <div class="metric-card">
                    <div class="metric-value">{pre_stats['total_variants']:,}</div>
                    <div class="metric-label">Pre-imputation Variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{post_stats['total_variants']:,}</div>
                    <div class="metric-label">Post-imputation Variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value increase">+{total_gain:,}</div>
                    <div class="metric-label">Variants Added</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value increase">{fold_change:.2f}x</div>
                    <div class="metric-label">Fold Increase</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value increase">{pct_increase:.1f}%</div>
                    <div class="metric-label">Percentage Increase</div>
                </div>
            </div>
        </div>
        
        <div class="summary-box">
            <h2>Detailed Comparison</h2>
            <table>
                <tr>
                    <th>Metric</th>
                    <th>Pre-imputation</th>
                    <th>Post-imputation</th>
                    <th>Change</th>
                    <th>% Change</th>
                </tr>
                <tr>
                    <td>Total Variants</td>
                    <td>{pre_stats['total_variants']:,}</td>
                    <td>{post_stats['total_variants']:,}</td>
                    <td class="increase">+{total_gain:,}</td>
                    <td class="increase">{pct_increase:.1f}%</td>
                </tr>
                <tr>
                    <td>SNPs</td>
                    <td>{pre_stats['snps']:,}</td>
                    <td>{post_stats['snps']:,}</td>
                    <td class="increase">+{post_stats['snps'] - pre_stats['snps']:,}</td>
                    <td class="increase">{((post_stats['snps'] - pre_stats['snps'])/max(pre_stats['snps'],1)*100):.1f}%</td>
                </tr>
                <tr>
                    <td>INDELs</td>
                    <td>{pre_stats['indels']:,}</td>
                    <td>{post_stats['indels']:,}</td>
                    <td class="increase">+{post_stats['indels'] - pre_stats['indels']:,}</td>
                    <td class="increase">{((post_stats['indels'] - pre_stats['indels'])/max(pre_stats['indels'],1)*100):.1f}%</td>
                </tr>
                <tr>
                    <td>Multiallelic</td>
                    <td>{pre_stats['multiallelic']:,}</td>
                    <td>{post_stats['multiallelic']:,}</td>
                    <td class="increase">+{post_stats['multiallelic'] - pre_stats['multiallelic']:,}</td>
                    <td class="increase">{((post_stats['multiallelic'] - pre_stats['multiallelic'])/max(pre_stats['multiallelic'],1)*100):.1f}%</td>
                </tr>
            </table>
        </div>
        
        <div class="plot-container">
            <div class="plot-title">Variant Count Comparison</div>
            <img src="{output_prefix}.variant_counts.png" alt="Variant Counts">
        </div>
        
        <div class="plot-container">
            <div class="plot-title">MAF Distribution Analysis</div>
            <img src="{output_prefix}.maf_distribution.png" alt="MAF Distribution">
        </div>
        
        <div class="plot-container">
            <div class="plot-title">Coverage Improvement</div>
            <img src="{output_prefix}.coverage_improvement.png" alt="Coverage Improvement">
        </div>
        
        <div class="plot-container">
            <div class="plot-title">Comprehensive Variant Gain Analysis</div>
            <img src="{output_prefix}.variant_gain.png" alt="Variant Gain">
        </div>
        
    </body>
    </html>
    """
    
    with open(f'{output_prefix}.comparison_report.html', 'w') as f:
        f.write(html_content)

def write_stats_file(pre_stats, post_stats, output_prefix):
    """Write detailed statistics to text file"""
    with open(f'{output_prefix}.comparison_stats.txt', 'w') as f:
        f.write("PRE/POST-IMPUTATION COMPARISON STATISTICS\\n")
        f.write("=" * 60 + "\\n\\n")
        
        f.write("SUMMARY\\n")
        f.write("-" * 30 + "\\n")
        f.write(f"Pre-imputation variants:  {pre_stats['total_variants']:,}\\n")
        f.write(f"Post-imputation variants: {post_stats['total_variants']:,}\\n")
        f.write(f"Variants added:           {post_stats['total_variants'] - pre_stats['total_variants']:,}\\n")
        f.write(f"Fold increase:            {post_stats['total_variants']/max(pre_stats['total_variants'],1):.2f}x\\n")
        f.write(f"Percentage increase:      {((post_stats['total_variants'] - pre_stats['total_variants'])/max(pre_stats['total_variants'],1)*100):.1f}%\\n")
        f.write("\\n")
        
        f.write("VARIANT TYPES\\n")
        f.write("-" * 30 + "\\n")
        f.write(f"{'Type':<15} {'Pre':>12} {'Post':>12} {'Gain':>12}\\n")
        f.write(f"{'SNPs':<15} {pre_stats['snps']:>12,} {post_stats['snps']:>12,} {post_stats['snps']-pre_stats['snps']:>+12,}\\n")
        f.write(f"{'INDELs':<15} {pre_stats['indels']:>12,} {post_stats['indels']:>12,} {post_stats['indels']-pre_stats['indels']:>+12,}\\n")
        f.write(f"{'Multiallelic':<15} {pre_stats['multiallelic']:>12,} {post_stats['multiallelic']:>12,} {post_stats['multiallelic']-pre_stats['multiallelic']:>+12,}\\n")
        f.write("\\n")
        
        f.write("MAF DISTRIBUTION\\n")
        f.write("-" * 30 + "\\n")
        maf_bins = ['Monomorphic', '0 < MAF < 0.001', '0.001 ≤ MAF < 0.005', 
                    '0.005 ≤ MAF < 0.01', '0.01 ≤ MAF < 0.05', '0.05 ≤ MAF < 0.1',
                    '0.1 ≤ MAF < 0.2', 'MAF ≥ 0.2']
        
        f.write(f"{'MAF Bin':<20} {'Pre':>12} {'Post':>12} {'Gain':>12}\\n")
        for bin in maf_bins:
            pre_count = pre_stats['maf_bins'].get(bin, 0)
            post_count = post_stats['maf_bins'].get(bin, 0)
            f.write(f"{bin:<20} {pre_count:>12,} {post_count:>12,} {post_count-pre_count:>+12,}\\n")
        f.write("\\n")
        
        f.write("CHROMOSOME DISTRIBUTION\\n")
        f.write("-" * 30 + "\\n")
        chroms = sorted(set(list(pre_stats['chr_counts'].keys()) + 
                           list(post_stats['chr_counts'].keys())),
                       key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').zfill(2)))
        
        f.write(f"{'Chromosome':<12} {'Pre':>12} {'Post':>12} {'Gain':>12} {'% Increase':>12}\\n")
        for chr in chroms:
            pre_count = pre_stats['chr_counts'].get(chr, 0)
            post_count = post_stats['chr_counts'].get(chr, 0)
            pct_inc = ((post_count - pre_count) / max(pre_count, 1) * 100) if pre_count > 0 else 0
            f.write(f"{chr:<12} {pre_count:>12,} {post_count:>12,} {post_count-pre_count:>+12,} {pct_inc:>11.1f}%\\n")

def main():
    """Main analysis function"""
    prefix = "${prefix}"
    pre_vcf = "${pre_vcf}"
    post_vcf = "${post_vcf}"
    
    print(f"Analyzing pre-imputation VCF: {pre_vcf}")
    pre_stats = extract_vcf_stats(pre_vcf, 'Pre-imputation')
    
    print(f"Analyzing post-imputation VCF: {post_vcf}")
    post_stats = extract_vcf_stats(post_vcf, 'Post-imputation')
    
    print("Generating comparison plots...")
    plot_variant_counts(pre_stats, post_stats, prefix)
    plot_maf_distribution(pre_stats, post_stats, prefix)
    plot_coverage_improvement(pre_stats, post_stats, prefix)
    plot_variant_gain(pre_stats, post_stats, prefix)
    
    print("Writing statistics file...")
    write_stats_file(pre_stats, post_stats, prefix)
    
    print("Generating HTML report...")
    generate_html_report(pre_stats, post_stats, prefix)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()
PYEOF

    # Run the analysis
    python3 compare_imputation.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python3 -c "import seaborn; print(seaborn.__version__)")
        cyvcf2: \$(python3 -c "import cyvcf2; print(cyvcf2.__version__)")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.comparison_stats.txt
    touch ${prefix}.variant_counts.png
    touch ${prefix}.maf_distribution.png
    touch ${prefix}.coverage_improvement.png
    touch ${prefix}.variant_gain.png
    touch ${prefix}.comparison_report.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        pandas: 1.5.0
        matplotlib: 3.6.0
        seaborn: 0.12.0
    END_VERSIONS
    """
}