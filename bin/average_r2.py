#!/usr/bin/env python3
"""
Calculate average R² statistics from Minimac4 sites VCF files
"""

import sys
import argparse
import pandas as pd
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
                
                data.append({'r2': r2, 'maf': maf, 'file': vcf_file})
    
    return pd.DataFrame(data)


def main():
    parser = argparse.ArgumentParser(description='Calculate average R² statistics')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--output-txt', required=True, help='Output text file for average R²')
    parser.add_argument('--output-csv', required=True, help='Output CSV file for detailed statistics')
    
    args = parser.parse_args()
    
    # Combine all info files
    all_data = []
    file_stats = []
    
    for vcf_file in glob.glob("*.sites.vcf.gz") + glob.glob("*.info*"):
        try:
            if vcf_file.endswith('.vcf.gz'):
                df = read_sites_vcf(vcf_file)
                if not df.empty:
                    # Calculate statistics for this file
                    stats = {
                        'file': vcf_file,
                        'n_variants': len(df),
                        'mean_r2': df['r2'].mean(),
                        'median_r2': df['r2'].median(),
                        'std_r2': df['r2'].std(),
                        'min_r2': df['r2'].min(),
                        'max_r2': df['r2'].max(),
                        'q25_r2': df['r2'].quantile(0.25),
                        'q75_r2': df['r2'].quantile(0.75)
                    }
                    
                    # Count well-imputed variants
                    for threshold in [0.3, 0.5, 0.8]:
                        stats[f'n_r2_ge_{threshold}'] = (df['r2'] >= threshold).sum()
                        stats[f'pct_r2_ge_{threshold}'] = (df['r2'] >= threshold).mean() * 100
                    
                    file_stats.append(stats)
                    all_data.append(df)
            else:
                # Try to read as tab-delimited info file
                df = pd.read_csv(vcf_file, sep='\t')
                if 'Rsq' in df.columns or 'R2' in df.columns:
                    # Get R² column
                    if 'Rsq' in df.columns:
                        r2_col = 'Rsq'
                    else:
                        r2_col = 'R2'
                    
                    # Calculate statistics for this file
                    stats = {
                        'file': vcf_file,
                        'n_variants': len(df),
                        'mean_r2': df[r2_col].mean(),
                        'median_r2': df[r2_col].median(),
                        'std_r2': df[r2_col].std(),
                        'min_r2': df[r2_col].min(),
                        'max_r2': df[r2_col].max(),
                        'q25_r2': df[r2_col].quantile(0.25),
                        'q75_r2': df[r2_col].quantile(0.75)
                    }
                    
                    # Count well-imputed variants
                    for threshold in [0.3, 0.5, 0.8]:
                        stats[f'n_r2_ge_{threshold}'] = (df[r2_col] >= threshold).sum()
                        stats[f'pct_r2_ge_{threshold}'] = (df[r2_col] >= threshold).mean() * 100
                    
                    file_stats.append(stats)
                    
                    # Add to all_data
                    temp_df = pd.DataFrame()
                    temp_df['r2'] = df[r2_col]
                    temp_df['file'] = vcf_file
                    
                    # Try to get MAF
                    if 'MAF' in df.columns:
                        temp_df['maf'] = df['MAF']
                    elif 'AF' in df.columns:
                        temp_df['maf'] = df['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else x)
                    else:
                        temp_df['maf'] = 0
                    
                    all_data.append(temp_df)
        except Exception as e:
            print(f"Warning: Could not process {vcf_file}: {e}", file=sys.stderr)
            continue
    
    if all_data:
        # Combine all data
        combined = pd.concat(all_data, ignore_index=True)
        all_r2 = combined['r2'].values
        
        # Calculate overall statistics
        overall_stats = {
            'Total variants': len(all_r2),
            'Mean R²': np.mean(all_r2),
            'Median R²': np.median(all_r2),
            'Std R²': np.std(all_r2),
            'Min R²': np.min(all_r2),
            'Max R²': np.max(all_r2),
            'Q25 R²': np.quantile(all_r2, 0.25),
            'Q75 R²': np.quantile(all_r2, 0.75)
        }
        
        # Write average R² file
        with open(args.output_txt, 'w') as f:
            f.write(f"Average R² for {args.sample_id} - {args.ref_name}\n")
            f.write("=" * 50 + "\n\n")
            
            for key, value in overall_stats.items():
                if 'Total' in key:
                    f.write(f"{key}: {value:,}\n")
                else:
                    f.write(f"{key}: {value:.4f}\n")
            
            f.write("\nR² Thresholds:\n")
            f.write("-" * 30 + "\n")
            for threshold in [0.3, 0.5, 0.8]:
                n_pass = (all_r2 >= threshold).sum()
                pct_pass = (all_r2 >= threshold).mean() * 100
                f.write(f"R² ≥ {threshold}: {n_pass:,} variants ({pct_pass:.1f}%)\n")
            
            # MAF-specific stats if available
            if 'maf' in combined.columns and combined['maf'].sum() > 0:
                f.write("\nMAF Bins:\n")
                f.write("-" * 30 + "\n")
                
                maf_bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5)]
                for low, high in maf_bins:
                    mask = (combined['maf'] >= low) & (combined['maf'] < high)
                    if mask.any():
                        mean_r2 = combined.loc[mask, 'r2'].mean()
                        n_vars = mask.sum()
                        f.write(f"MAF [{low:.2f}-{high:.2f}): mean R²={mean_r2:.3f}, n={n_vars:,}\n")
        
        # Write summary CSV
        if file_stats:
            summary_df = pd.DataFrame(file_stats)
            summary_df.to_csv(args.output_csv, index=False)
            print(f"Processed {len(file_stats)} info files")
            print(f"Overall mean R²: {overall_stats['Mean R²']:.4f}")
            print(f"Total variants: {overall_stats['Total variants']:,}")
        
    else:
        # No data found
        with open(args.output_txt, 'w') as f:
            f.write("No imputation info data found\n")
        
        pd.DataFrame().to_csv(args.output_csv, index=False)
        print("No data found to process")


if __name__ == '__main__':
    main()