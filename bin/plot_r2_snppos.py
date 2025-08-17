#!/usr/bin/env python3
"""
Plot R² imputation quality vs SNP position from Minimac4 sites VCF file
"""

import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import gzip


def parse_info_field(info_str):
    """Parse VCF INFO field to extract key-value pairs"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def read_sites_vcf(vcf_file):
    """Read sites VCF file and extract relevant information"""
    data = []
    
    # Open file (handle gzipped and plain text)
    if vcf_file.endswith('.gz'):
        opener = gzip.open(vcf_file, 'rt')
    else:
        opener = open(vcf_file, 'r')
    
    with opener as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            info = parse_info_field(parts[7])
            
            # Extract R2 value
            r2 = float(info.get('R2', 0))
            maf = float(info.get('MAF', 0))
            
            data.append({
                'chr': chrom,
                'pos': pos,
                'r2': r2,
                'maf': maf
            })
    
    return pd.DataFrame(data)


def plot_r2_vs_position(df, output_file, sample_id, ref_name):
    """Create R² vs position plot"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot R² vs position
    scatter = ax.scatter(df['pos'], df['r2'], alpha=0.5, s=1, c=df['r2'], 
                        cmap='RdYlGn', vmin=0, vmax=1)
    
    # Add rolling mean if we have enough points
    if len(df) > 10:
        window_size = max(10, len(df) // 100)
        df_sorted = df.sort_values('pos')
        rolling_mean = df_sorted['r2'].rolling(window=window_size, center=True).mean()
        ax.plot(df_sorted['pos'], rolling_mean, 'r-', linewidth=2, 
                label=f'Rolling mean (window={window_size})')
    
    # Add horizontal lines for thresholds
    ax.axhline(y=0.3, color='orange', linestyle='--', alpha=0.5, label='R²=0.3')
    ax.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='R²=0.8')
    
    ax.set_xlabel('Position')
    ax.set_ylabel('R²')
    ax.set_title(f'Imputation Quality (R²) vs SNP Position\n{sample_id} - {ref_name}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    plt.colorbar(scatter, ax=ax, label='R²')
    
    # Add statistics text
    stats_text = f"Mean R²: {df['r2'].mean():.3f}\n"
    stats_text += f"Median R²: {df['r2'].median():.3f}\n"
    stats_text += f"SNPs with R²≥0.3: {(df['r2'] >= 0.3).sum():,} ({(df['r2'] >= 0.3).mean()*100:.1f}%)\n"
    stats_text += f"SNPs with R²≥0.8: {(df['r2'] >= 0.8).sum():,} ({(df['r2'] >= 0.8).mean()*100:.1f}%)"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    
    print(f"Plot saved: {output_file}")
    print(f"Total variants analyzed: {len(df):,}")
    print(f"Mean R²: {df['r2'].mean():.4f}")


def main():
    parser = argparse.ArgumentParser(description='Plot R² vs SNP position from sites VCF')
    parser.add_argument('input_vcf', help='Input sites VCF file from Minimac4')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample ID for plot title')
    parser.add_argument('--ref-name', required=True, help='Reference panel name for plot title')
    
    args = parser.parse_args()
    
    # Read VCF data
    print(f"Reading sites VCF: {args.input_vcf}")
    df = read_sites_vcf(args.input_vcf)
    
    if df.empty:
        print("Warning: No data found in VCF file")
        # Create empty plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_title(f'Imputation Quality (R²) vs SNP Position\n{args.sample_id} - {args.ref_name}')
        plt.savefig(args.output_pdf)
    else:
        # Create plot
        plot_r2_vs_position(df, args.output_pdf, args.sample_id, args.ref_name)


if __name__ == '__main__':
    main()