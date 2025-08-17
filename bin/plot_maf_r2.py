#!/usr/bin/env python3
"""
Plot MAF vs R² from Minimac4 sites VCF files
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
    """Read sites VCF file and extract R2 and MAF information"""
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
            
            # Extract R2 and MAF values
            if 'R2' in info:
                r2 = float(info.get('R2', 0))
                maf = float(info.get('MAF', 0))
                af = float(info.get('AF', 0))
                
                # If MAF is 0 but AF exists, calculate MAF from AF
                if maf == 0 and af > 0:
                    maf = min(af, 1 - af)
                
                data.append({'r2': r2, 'maf': maf})
    
    return pd.DataFrame(data)


def main():
    parser = argparse.ArgumentParser(description='Plot MAF vs R²')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample ID for plot title')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--info-cutoff', type=float, default=0.3, help='R² cutoff threshold')
    
    args = parser.parse_args()
    
    # Combine all info files
    all_data = []
    for vcf_file in glob.glob("*.sites.vcf.gz") + glob.glob("*.info*"):
        try:
            if vcf_file.endswith('.vcf.gz'):
                df = read_sites_vcf(vcf_file)
                if not df.empty:
                    all_data.append(df)
            else:
                # Try to read as tab-delimited info file
                df = pd.read_csv(vcf_file, sep='\t')
                if ('Rsq' in df.columns or 'R2' in df.columns) and ('MAF' in df.columns or 'AF' in df.columns):
                    if 'Rsq' in df.columns:
                        df['r2'] = df['Rsq']
                    elif 'R2' in df.columns:
                        df['r2'] = df['R2']
                    
                    if 'MAF' in df.columns:
                        df['maf'] = df['MAF']
                    elif 'AF' in df.columns:
                        df['maf'] = df['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else x)
                    
                    all_data.append(df[['r2', 'maf']])
        except:
            continue
    
    if not all_data:
        # Create empty plot if no data
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center')
        ax.set_title('MAF vs R²')
        plt.savefig(args.output_pdf)
    else:
        combined = pd.concat(all_data, ignore_index=True)
        combined = combined.dropna()
        
        # Create MAF bins
        maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        combined['maf_bin'] = pd.cut(combined['maf'], bins=maf_bins)
        
        # Calculate statistics for each MAF bin
        stats = combined.groupby('maf_bin').agg({
            'r2': ['mean', 'median', 'std', 'count']
        }).round(3)
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top plot: Scatter plot of MAF vs R²
        scatter = ax1.scatter(combined['maf'], combined['r2'], 
                             alpha=0.3, s=1, c=combined['r2'], 
                             cmap='RdYlGn', vmin=0, vmax=1)
        
        # Add mean R² line for each MAF bin
        for maf_bin in stats.index:
            if pd.notna(maf_bin):
                bin_data = combined[combined['maf_bin'] == maf_bin]
                if len(bin_data) > 0:
                    maf_center = (maf_bin.left + maf_bin.right) / 2
                    mean_r2 = stats.loc[maf_bin, ('r2', 'mean')]
                    ax1.plot(maf_center, mean_r2, 'ko', markersize=8)
        
        ax1.axhline(y=args.info_cutoff, color='blue', linestyle='--', 
                    linewidth=1, label=f'R² threshold: {args.info_cutoff}')
        ax1.set_xlabel('Minor Allele Frequency (MAF)')
        ax1.set_ylabel('R²')
        ax1.set_title(f'MAF vs Imputation Quality (R²)\n{args.sample_id} - {args.ref_name}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax1, label='R²')
        
        # Bottom plot: Box plot by MAF bins
        maf_labels = [f"{b.left:.2f}-{b.right:.2f}" for b in stats.index if pd.notna(b)]
        box_data = [combined[combined['maf_bin'] == b]['r2'].values 
                    for b in stats.index if pd.notna(b)]
        
        bp = ax2.boxplot([d for d in box_data if len(d) > 0], 
                         labels=[l for l, d in zip(maf_labels, box_data) if len(d) > 0], 
                         patch_artist=True)
        
        # Color boxes by mean R²
        for patch, maf_bin in zip(bp['boxes'], [b for b in stats.index if pd.notna(b) and len(combined[combined['maf_bin'] == b]) > 0]):
            mean_r2 = stats.loc[maf_bin, ('r2', 'mean')]
            if mean_r2 >= args.info_cutoff:
                patch.set_facecolor('lightgreen')
            else:
                patch.set_facecolor('lightcoral')
        
        ax2.axhline(y=args.info_cutoff, color='blue', linestyle='--', 
                    linewidth=1, label=f'R² threshold: {args.info_cutoff}')
        ax2.set_xlabel('MAF Bins')
        ax2.set_ylabel('R²')
        ax2.set_title('R² Distribution by MAF Bins')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Add statistics table
        stats_text = "MAF Bin Statistics:\n"
        for i, (maf_bin, row) in enumerate(stats.iterrows()):
            if pd.notna(maf_bin) and i < 5:  # Show first 5 bins
                stats_text += f"{maf_bin.left:.2f}-{maf_bin.right:.2f}: "
                stats_text += f"mean={row[('r2', 'mean')]:.3f}, n={int(row[('r2', 'count')])}\n"
        
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                 va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.savefig(args.output_pdf)
        
        print(f"Plot saved: {args.output_pdf}")
        print(f"Total variants analyzed: {len(combined):,}")
        print(f"Mean R²: {combined['r2'].mean():.4f}")


if __name__ == '__main__':
    main()