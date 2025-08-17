#!/usr/bin/env python3
"""
Plot histogram and cumulative distribution of R² from Minimac4 sites VCF files
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
    parser = argparse.ArgumentParser(description='Plot histogram and CDF of R²')
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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax1.set_title('R² Histogram')
        ax2.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax2.set_title('Cumulative Distribution')
        plt.savefig(args.output_pdf)
        plt.close()
        print(f"Empty plot saved: {args.output_pdf}")
    else:
        all_r2 = np.array(all_r2)
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Left plot: Histogram
        n, bins, patches = ax1.hist(all_r2, bins=50, edgecolor='black', alpha=0.7)
        
        # Color bars by threshold
        for i, patch in enumerate(patches):
            if bins[i] >= args.info_cutoff:
                patch.set_facecolor('green')
            else:
                patch.set_facecolor('red')
        
        ax1.axvline(x=args.info_cutoff, color='blue', linestyle='--', 
                    linewidth=2, label=f'R² threshold: {args.info_cutoff}')
        ax1.set_xlabel('R²')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Distribution of R² Values')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Right plot: Cumulative distribution
        sorted_r2 = np.sort(all_r2)
        cumulative = np.arange(1, len(sorted_r2) + 1) / len(sorted_r2)
        
        ax2.plot(sorted_r2, cumulative, 'b-', linewidth=2)
        ax2.axvline(x=args.info_cutoff, color='red', linestyle='--', 
                    linewidth=2, label=f'R² threshold: {args.info_cutoff}')
        ax2.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5, label='Median')
        ax2.set_xlabel('R²')
        ax2.set_ylabel('Cumulative Proportion')
        ax2.set_title('Cumulative Distribution of R² Values')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim([0, 1])
        ax2.set_ylim([0, 1])
        
        # Overall title
        plt.suptitle(f'R² Distribution Analysis\n{args.sample_id} - {args.ref_name}', fontsize=14)
        
        # Add statistics text
        mean_r2 = np.mean(all_r2)
        median_r2 = np.median(all_r2)
        q25 = np.percentile(all_r2, 25)
        q75 = np.percentile(all_r2, 75)
        n_pass = np.sum(all_r2 >= args.info_cutoff)
        pct_pass = n_pass / len(all_r2) * 100
        
        stats_text = f"Statistics:\n"
        stats_text += f"Mean: {mean_r2:.3f}\n"
        stats_text += f"Median: {median_r2:.3f}\n"
        stats_text += f"Q25-Q75: {q25:.3f}-{q75:.3f}\n"
        stats_text += f"Total SNPs: {len(all_r2):,}\n"
        stats_text += f"R²≥{args.info_cutoff}: {n_pass:,} ({pct_pass:.1f}%)"
        
        fig.text(0.02, 0.98, stats_text, transform=fig.transFigure, 
                 va='top', fontsize=9,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        plt.savefig(args.output_pdf)
        plt.close()
        
        print(f"Plot saved: {args.output_pdf}")
        print(f"Total SNPs: {len(all_r2):,}")
        print(f"Mean R²: {mean_r2:.4f}")
        print(f"Median R²: {median_r2:.4f}")
        
        # Ensure file is written
        import os
        if os.path.exists(args.output_pdf):
            print(f"File size: {os.path.getsize(args.output_pdf)} bytes")
        else:
            print(f"ERROR: File {args.output_pdf} was not created!")
            sys.exit(1)


if __name__ == '__main__':
    main()