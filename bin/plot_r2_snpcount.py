#!/usr/bin/env python3
"""
Plot SNP count by R² bins from Minimac4 sites VCF files
"""

import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
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
    """Read sites VCF file and extract R2 information"""
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
            
            info = parse_info_field(parts[7])
            
            # Extract R2 value
            if 'R2' in info:
                r2 = float(info.get('R2', 0))
                data.append(r2)
    
    return data


def main():
    parser = argparse.ArgumentParser(description='Plot SNP count by R² bins')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample ID for plot title')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--info-cutoff', type=float, default=0.3, help='R² cutoff threshold')
    
    args = parser.parse_args()
    
    # Combine all info files
    all_r2 = []
    for vcf_file in glob.glob("*.sites.vcf.gz") + glob.glob("*.info*"):
        try:
            if vcf_file.endswith('.vcf.gz'):
                r2_values = read_sites_vcf(vcf_file)
                all_r2.extend(r2_values)
            else:
                # Try to read as tab-delimited info file
                df = pd.read_csv(vcf_file, sep='\t')
                if 'Rsq' in df.columns:
                    all_r2.extend(df['Rsq'].values)
                elif 'R2' in df.columns:
                    all_r2.extend(df['R2'].values)
        except:
            continue
    
    if not all_r2:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_title('R² vs SNP Count')
        plt.savefig(args.output_pdf)
    else:
        all_r2 = np.array(all_r2)
        
        # Create R² bins
        r2_bins = np.arange(0, 1.05, 0.05)
        counts, bin_edges = np.histogram(all_r2, bins=r2_bins)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Bar plot
        x_pos = np.arange(len(counts))
        bars = ax.bar(x_pos, counts, alpha=0.7)
        
        # Color bars by R² threshold
        for i, bar in enumerate(bars):
            if bin_edges[i] >= args.info_cutoff:
                bar.set_color('green')
            else:
                bar.set_color('red')
        
        # Labels
        ax.set_xlabel('R² Bin')
        ax.set_ylabel('Number of SNPs')
        ax.set_title(f'SNP Count by R² Bin\n{args.sample_id} - {args.ref_name}')
        
        # X-axis labels
        labels = [f"{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}" for i in range(len(counts))]
        ax.set_xticks(x_pos[::2])  # Show every other label to avoid crowding
        ax.set_xticklabels(labels[::2], rotation=45, ha='right')
        
        # Add threshold line
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.text(0.02, 0.98, f"R² threshold: {args.info_cutoff}", 
                transform=ax.transAxes, va='top', fontsize=10)
        
        # Add statistics
        mean_r2 = np.mean(all_r2)
        median_r2 = np.median(all_r2)
        n_pass = np.sum(all_r2 >= args.info_cutoff)
        pct_pass = n_pass / len(all_r2) * 100
        
        stats_text = f"Mean R²: {mean_r2:.3f}\n"
        stats_text += f"Median R²: {median_r2:.3f}\n"
        stats_text += f"Total SNPs: {len(all_r2):,}\n"
        stats_text += f"SNPs with R²≥{args.info_cutoff}: {n_pass:,} ({pct_pass:.1f}%)"
        
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, 
                va='top', ha='right', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(args.output_pdf)
        
        print(f"Plot saved: {args.output_pdf}")
        print(f"Total SNPs: {len(all_r2):,}")
        print(f"Mean R²: {mean_r2:.4f}")


if __name__ == '__main__':
    main()