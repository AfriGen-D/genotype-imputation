#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create aggregated R² score visualizations similar to Terra All of Us
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
import argparse
import json
import glob

def aggregate_info_files(info_pattern):
    """Aggregate multiple info files from different chromosomes"""
    all_data = []
    files = glob.glob(info_pattern)
    
    for info_file in files:
        # Extract chromosome from filename
        chrom = Path(info_file).stem.split('_')[0] if '_' in Path(info_file).stem else 'unknown'
        
        with open(info_file, 'r') as f:
            header = f.readline().strip().split('\t')
            if 'SNP' not in header:
                continue
                
            snp_idx = header.index('SNP')
            maf_idx = header.index('MAF') if 'MAF' in header else -1
            r2_idx = header.index('Rsq') if 'Rsq' in header else header.index('R2') if 'R2' in header else -1
            
            if maf_idx == -1 or r2_idx == -1:
                continue
            
            for line in f:
                parts = line.strip().split('\t')
                try:
                    variant = parts[snp_idx]
                    # Extract position from variant ID (format: chr:pos:ref:alt)
                    if ':' in variant:
                        pos = int(variant.split(':')[1])
                    else:
                        continue
                    
                    maf = float(parts[maf_idx])
                    r2 = float(parts[r2_idx])
                    all_data.append({
                        'chromosome': chrom,
                        'position': pos,
                        'MAF': maf,
                        'R2': r2
                    })
                except (ValueError, IndexError):
                    continue
    
    return pd.DataFrame(all_data)

def plot_aggregated_r2_dashboard(df, output_prefix, title="Aggregated Imputation Quality"):
    """Create comprehensive R² dashboard similar to Terra"""
    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)
    
    # Color scheme
    colors = {
        'excellent': '#2E7D32',  # Dark green
        'good': '#66BB6A',       # Light green
        'moderate': '#FFA726',   # Orange
        'poor': '#EF5350'        # Red
    }
    
    # 1. Overall R² distribution (top-left, spans 2 rows)
    ax1 = fig.add_subplot(gs[0:2, 0])
    r2_bins = np.linspace(0, 1, 51)
    counts, bins = np.histogram(df['R2'], bins=r2_bins)
    
    # Color bars by quality
    bar_colors = []
    for i, (left, right) in enumerate(zip(bins[:-1], bins[1:])):
        mid = (left + right) / 2
        if mid >= 0.8:
            bar_colors.append(colors['excellent'])
        elif mid >= 0.5:
            bar_colors.append(colors['good'])
        elif mid >= 0.3:
            bar_colors.append(colors['moderate'])
        else:
            bar_colors.append(colors['poor'])
    
    ax1.bar(bins[:-1], counts, width=np.diff(bins), align='edge', 
           color=bar_colors, edgecolor='black', linewidth=0.5)
    
    # Add threshold lines
    for threshold, label in [(0.3, 'Poor'), (0.5, 'Moderate'), (0.8, 'Good')]:
        ax1.axvline(threshold, color='black', linestyle='--', alpha=0.5)
        ax1.text(threshold, ax1.get_ylim()[1]*0.95, label, 
                rotation=0, ha='center', fontsize=10)
    
    ax1.set_xlabel('Imputation R²', fontsize=12)
    ax1.set_ylabel('Number of Variants', fontsize=12)
    ax1.set_title('Overall R² Distribution', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. R² by MAF categories (top-middle)
    ax2 = fig.add_subplot(gs[0, 1])
    
    maf_categories = [
        ('Ultra-rare\\n(MAF<0.1%)', df[df['MAF'] < 0.001]),
        ('Rare\\n(0.1-1%)', df[(df['MAF'] >= 0.001) & (df['MAF'] < 0.01)]),
        ('Low\\n(1-5%)', df[(df['MAF'] >= 0.01) & (df['MAF'] < 0.05)]),
        ('Common\\n(≥5%)', df[df['MAF'] >= 0.05])
    ]
    
    positions = []
    medians = []
    boxes = []
    
    for i, (label, data) in enumerate(maf_categories):
        if len(data) > 0:
            positions.append(i)
            box_data = data['R2'].dropna()
            boxes.append(box_data)
            medians.append(np.median(box_data))
    
    bp = ax2.boxplot(boxes, positions=positions, patch_artist=True, 
                     notch=True, showmeans=True)
    
    # Color boxes by median R²
    for patch, median in zip(bp['boxes'], medians):
        if median >= 0.8:
            color = colors['excellent']
        elif median >= 0.5:
            color = colors['good']
        elif median >= 0.3:
            color = colors['moderate']
        else:
            color = colors['poor']
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax2.set_xticks(positions)
    ax2.set_xticklabels([cat[0] for cat in maf_categories if len(cat[1]) > 0])
    ax2.set_ylabel('Imputation R²', fontsize=12)
    ax2.set_title('R² by Allele Frequency', fontsize=14, fontweight='bold')
    ax2.set_ylim([0, 1])
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Quality metrics pie chart (top-right)
    ax3 = fig.add_subplot(gs[0, 2])
    
    quality_thresholds = [
        ('Excellent (R²≥0.8)', (df['R2'] >= 0.8).sum(), colors['excellent']),
        ('Good (0.5≤R²<0.8)', ((df['R2'] >= 0.5) & (df['R2'] < 0.8)).sum(), colors['good']),
        ('Moderate (0.3≤R²<0.5)', ((df['R2'] >= 0.3) & (df['R2'] < 0.5)).sum(), colors['moderate']),
        ('Poor (R²<0.3)', (df['R2'] < 0.3).sum(), colors['poor'])
    ]
    
    labels = []
    sizes = []
    pie_colors = []
    
    for label, count, color in quality_thresholds:
        if count > 0:
            labels.append(f'{label}\\n({count:,} variants)')
            sizes.append(count)
            pie_colors.append(color)
    
    wedges, texts, autotexts = ax3.pie(sizes, labels=labels, colors=pie_colors, 
                                        autopct='%1.1f%%', startangle=90)
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax3.set_title('Imputation Quality Distribution', fontsize=14, fontweight='bold')
    
    # 4. Chromosome-wise mean R² (middle-left, spans 2 columns)
    ax4 = fig.add_subplot(gs[1, 1:3])
    
    if 'chromosome' in df.columns:
        chr_stats = df.groupby('chromosome').agg({
            'R2': ['mean', 'median', 'std', 'count']
        }).reset_index()
        chr_stats.columns = ['chromosome', 'mean_r2', 'median_r2', 'std_r2', 'count']
        
        # Sort chromosomes
        chr_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
        chr_stats['chr_num'] = chr_stats['chromosome'].map(
            {chr: i for i, chr in enumerate(chr_order)}
        )
        chr_stats = chr_stats.sort_values('chr_num')
        
        x = np.arange(len(chr_stats))
        width = 0.35
        
        bars1 = ax4.bar(x - width/2, chr_stats['mean_r2'], width, 
                       label='Mean R²', color='steelblue', alpha=0.8)
        bars2 = ax4.bar(x + width/2, chr_stats['median_r2'], width, 
                       label='Median R²', color='orange', alpha=0.8)
        
        # Add error bars for standard deviation
        ax4.errorbar(x - width/2, chr_stats['mean_r2'], 
                    yerr=chr_stats['std_r2'], fmt='none', 
                    color='black', alpha=0.5, capsize=3)
        
        ax4.set_xlabel('Chromosome', fontsize=12)
        ax4.set_ylabel('R² Score', fontsize=12)
        ax4.set_title('Chromosome-wise Imputation Quality', fontsize=14, fontweight='bold')
        ax4.set_xticks(x)
        ax4.set_xticklabels(chr_stats['chromosome'], rotation=45)
        ax4.legend()
        ax4.set_ylim([0, 1])
        ax4.grid(True, alpha=0.3, axis='y')
        
        # Add count as text above bars
        for i, (mean_bar, count) in enumerate(zip(bars1, chr_stats['count'])):
            ax4.text(mean_bar.get_x() + mean_bar.get_width()/2, 0.95,
                    f'{count:,}', ha='center', va='bottom', fontsize=8)
    
    # 5. MAF vs R² density plot (middle-right)
    ax5 = fig.add_subplot(gs[2, 0])
    
    # Sample for performance
    sample_size = min(50000, len(df))
    df_sample = df.sample(n=sample_size) if len(df) > sample_size else df
    
    hexbin = ax5.hexbin(df_sample['MAF'], df_sample['R2'], 
                        gridsize=50, cmap='YlOrRd', mincnt=1)
    
    # Add mean trend line
    maf_bins = np.linspace(0, 0.5, 21)
    bin_means = []
    bin_centers = []
    
    for i in range(len(maf_bins)-1):
        mask = (df['MAF'] >= maf_bins[i]) & (df['MAF'] < maf_bins[i+1])
        if mask.sum() > 0:
            bin_means.append(df[mask]['R2'].mean())
            bin_centers.append((maf_bins[i] + maf_bins[i+1]) / 2)
    
    ax5.plot(bin_centers, bin_means, 'b-', linewidth=2, label='Mean R²')
    
    ax5.set_xlabel('Minor Allele Frequency', fontsize=12)
    ax5.set_ylabel('Imputation R²', fontsize=12)
    ax5.set_title('MAF vs R² Density', fontsize=14, fontweight='bold')
    ax5.set_xlim([0, 0.5])
    ax5.set_ylim([0, 1])
    ax5.legend()
    
    cb = plt.colorbar(hexbin, ax=ax5)
    cb.set_label('Count', fontsize=10)
    
    # 6. Cumulative distribution by quality (bottom-left)
    ax6 = fig.add_subplot(gs[2, 1])
    
    r2_sorted = np.sort(df['R2'])
    cumulative = np.arange(1, len(r2_sorted) + 1) / len(r2_sorted)
    
    ax6.plot(r2_sorted, cumulative, linewidth=2, color='navy')
    
    # Add shaded regions
    ax6.axvspan(0, 0.3, alpha=0.2, color=colors['poor'], label='Poor')
    ax6.axvspan(0.3, 0.5, alpha=0.2, color=colors['moderate'], label='Moderate')
    ax6.axvspan(0.5, 0.8, alpha=0.2, color=colors['good'], label='Good')
    ax6.axvspan(0.8, 1.0, alpha=0.2, color=colors['excellent'], label='Excellent')
    
    # Add percentage lines
    percentiles = [0.25, 0.5, 0.75, 0.9]
    for p in percentiles:
        r2_val = np.percentile(df['R2'], p * 100)
        ax6.axhline(p, color='gray', linestyle=':', alpha=0.5)
        ax6.axvline(r2_val, color='gray', linestyle=':', alpha=0.5)
        ax6.text(0.02, p, f'{int(p*100)}%', fontsize=8)
        ax6.text(r2_val, 0.02, f'{r2_val:.2f}', rotation=90, fontsize=8)
    
    ax6.set_xlabel('Imputation R²', fontsize=12)
    ax6.set_ylabel('Cumulative Proportion', fontsize=12)
    ax6.set_title('Cumulative R² Distribution', fontsize=14, fontweight='bold')
    ax6.set_xlim([0, 1])
    ax6.set_ylim([0, 1])
    ax6.legend(loc='lower right')
    ax6.grid(True, alpha=0.3)
    
    # 7. Summary statistics table (bottom-middle and right)
    ax7 = fig.add_subplot(gs[2:4, 2])
    ax7.axis('tight')
    ax7.axis('off')
    
    # Calculate statistics
    stats_data = [
        ['Total Variants', f'{len(df):,}'],
        ['Mean R²', f'{df["R2"].mean():.4f}'],
        ['Median R²', f'{df["R2"].median():.4f}'],
        ['Std Dev R²', f'{df["R2"].std():.4f}'],
        ['', ''],
        ['Well Imputed (R²≥0.3)', f'{(df["R2"] >= 0.3).sum():,} ({(df["R2"] >= 0.3).mean():.1%})'],
        ['High Quality (R²≥0.5)', f'{(df["R2"] >= 0.5).sum():,} ({(df["R2"] >= 0.5).mean():.1%})'],
        ['Excellent (R²≥0.8)', f'{(df["R2"] >= 0.8).sum():,} ({(df["R2"] >= 0.8).mean():.1%})'],
        ['', ''],
        ['MAF Categories', ''],
        ['  Rare (<1%)', f'{(df["MAF"] < 0.01).sum():,} variants'],
        ['  Low (1-5%)', f'{((df["MAF"] >= 0.01) & (df["MAF"] < 0.05)).sum():,} variants'],
        ['  Common (≥5%)', f'{(df["MAF"] >= 0.05).sum():,} variants']
    ]
    
    table = ax7.table(cellText=stats_data,
                     colLabels=['Metric', 'Value'],
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.6, 0.4])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    
    # Style the table
    for i in range(len(stats_data) + 1):
        if i == 0:  # Header
            table[(i, 0)].set_facecolor('#E0E0E0')
            table[(i, 1)].set_facecolor('#E0E0E0')
            table[(i, 0)].set_text_props(weight='bold')
            table[(i, 1)].set_text_props(weight='bold')
        elif i in [5, 10]:  # Section separators
            table[(i, 0)].set_facecolor('#F5F5F5')
            table[(i, 1)].set_facecolor('#F5F5F5')
    
    ax7.set_title('Summary Statistics', fontsize=14, fontweight='bold', pad=20)
    
    # 8. Sample size by MAF (bottom-left)
    ax8 = fig.add_subplot(gs[3, 0:2])
    
    maf_bins = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    df['MAF_bin'] = pd.cut(df['MAF'], bins=maf_bins, include_lowest=True)
    maf_counts = df.groupby('MAF_bin').size()
    
    x_pos = np.arange(len(maf_counts))
    bars = ax8.bar(x_pos, maf_counts.values, color='steelblue', edgecolor='black')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax8.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}', ha='center', va='bottom', fontsize=8)
    
    ax8.set_xlabel('MAF Bins', fontsize=12)
    ax8.set_ylabel('Number of Variants (log scale)', fontsize=12)
    ax8.set_title('Variant Count by MAF', fontsize=14, fontweight='bold')
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels([f'{b:.4f}' if b < 0.01 else f'{b:.2f}' 
                         for b in maf_bins[1:]], rotation=45)
    ax8.set_yscale('log')
    ax8.grid(True, alpha=0.3, axis='y')
    
    # Main title
    plt.suptitle(f'{title}\\nComprehensive Imputation Quality Dashboard', 
                fontsize=18, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig(f'{output_prefix}_aggregated_dashboard.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_aggregated_dashboard.png', dpi=150, bbox_inches='tight')
    
    # Save statistics to JSON
    save_aggregated_stats(df, output_prefix)

def save_aggregated_stats(df, output_prefix):
    """Save aggregated statistics to JSON"""
    stats = {
        'summary': {
            'total_variants': len(df),
            'mean_r2': float(df['R2'].mean()),
            'median_r2': float(df['R2'].median()),
            'std_r2': float(df['R2'].std()),
            'min_r2': float(df['R2'].min()),
            'max_r2': float(df['R2'].max())
        },
        'quality_thresholds': {
            'r2_ge_0.3': int((df['R2'] >= 0.3).sum()),
            'r2_ge_0.5': int((df['R2'] >= 0.5).sum()),
            'r2_ge_0.8': int((df['R2'] >= 0.8).sum()),
            'prop_r2_ge_0.3': float((df['R2'] >= 0.3).mean()),
            'prop_r2_ge_0.5': float((df['R2'] >= 0.5).mean()),
            'prop_r2_ge_0.8': float((df['R2'] >= 0.8).mean())
        },
        'maf_categories': {}
    }
    
    # MAF category statistics
    maf_cats = [
        ('ultra_rare', df[df['MAF'] < 0.001]),
        ('rare', df[(df['MAF'] >= 0.001) & (df['MAF'] < 0.01)]),
        ('low_freq', df[(df['MAF'] >= 0.01) & (df['MAF'] < 0.05)]),
        ('common', df[df['MAF'] >= 0.05])
    ]
    
    for name, cat_df in maf_cats:
        if len(cat_df) > 0:
            stats['maf_categories'][name] = {
                'count': int(len(cat_df)),
                'mean_r2': float(cat_df['R2'].mean()),
                'median_r2': float(cat_df['R2'].median()),
                'prop_well_imputed': float((cat_df['R2'] >= 0.3).mean())
            }
    
    with open(f'{output_prefix}_aggregated_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description='Create aggregated R² score visualizations')
    parser.add_argument('info_files', help='Pattern for info files or single file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--title', default='Imputation Quality', help='Dashboard title')
    
    args = parser.parse_args()
    
    # Load data
    if '*' in args.info_files or '?' in args.info_files:
        df = aggregate_info_files(args.info_files)
    else:
        # Single file
        df = aggregate_info_files(args.info_files)
    
    if df.empty:
        print(f"Error: No data found in {args.info_files}")
        sys.exit(1)
    
    # Generate dashboard
    plot_aggregated_r2_dashboard(df, args.output_prefix, args.title)
    
    print(f"Dashboard saved to {args.output_prefix}_aggregated_dashboard.pdf")
    print(f"Statistics saved to {args.output_prefix}_aggregated_stats.json")

if __name__ == '__main__':
    main()