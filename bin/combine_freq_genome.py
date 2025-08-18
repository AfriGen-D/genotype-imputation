#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='Combine chromosome frequency files into genome-wide frequencies')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--input-pattern', default='*.chr_freq.tsv', help='Input file pattern for chromosome frequency files')
    
    args = parser.parse_args()
    
    # Read all chromosome frequency files
    chr_files = glob.glob(args.input_pattern)
    
    print(f"Found {len(chr_files)} chromosome frequency files")
    
    # Combine all chromosome data
    all_dfs = []
    chr_stats = []
    
    for f in sorted(chr_files):
        df = pd.read_csv(f, sep='\t')
        if not df.empty:
            all_dfs.append(df)
            
            # Calculate per-chromosome statistics
            chr_name = f.replace('.chr_freq.tsv', '').split('_')[-1]
            stats = {
                'chromosome': chr_name,
                'n_variants': len(df),
                'n_shared': ((df['AF_imputed'] > 0) & (df['AF_reference'] > 0)).sum() if 'AF_imputed' in df.columns and 'AF_reference' in df.columns else 0,
                'n_imputed_only': (df['AF_reference'] == 0).sum() if 'AF_reference' in df.columns else 0,
                'n_reference_only': (df['AF_imputed'] == 0).sum() if 'AF_imputed' in df.columns else 0
            }
            
            if 'AF_diff' in df.columns:
                stats['mean_af_diff'] = df['AF_diff'].mean()
                stats['median_af_diff'] = df['AF_diff'].median()
                stats['std_af_diff'] = df['AF_diff'].std()
                stats['max_abs_diff'] = df['AF_diff'].abs().max()
            
            chr_stats.append(stats)
    
    if all_dfs:
        # Combine all chromosome data
        genome_df = pd.concat(all_dfs, ignore_index=True)
        
        # Sort by chromosome and position
        if 'CHROM' in genome_df.columns and 'POS' in genome_df.columns:
            # Try to sort chromosomes numerically where possible
            def chr_sort_key(x):
                try:
                    # Extract numeric part if present
                    chr_str = str(x).replace('chr', '')
                    if chr_str.isdigit():
                        return (0, int(chr_str))
                    elif chr_str == 'X':
                        return (1, 23)
                    elif chr_str == 'Y':
                        return (1, 24)
                    elif chr_str == 'MT' or chr_str == 'M':
                        return (1, 25)
                    else:
                        return (2, chr_str)
                except:
                    return (3, str(x))
            
            genome_df['chr_sort'] = genome_df['CHROM'].apply(chr_sort_key)
            genome_df = genome_df.sort_values(['chr_sort', 'POS'])
            genome_df = genome_df.drop('chr_sort', axis=1)
        
        # Save genome-wide frequencies
        genome_df.to_csv(f"{args.output_prefix}.genome_freq.tsv", sep='\t', index=False)
        
        print(f"Combined {len(genome_df)} variants genome-wide")
        
        # Create summary report
        with open(f"{args.output_prefix}.genome_freq_summary.txt", 'w') as f:
            f.write("Genome-wide Frequency Comparison Summary\n")
            f.write("=" * 70 + "\n\n")
            
            # Overall statistics
            f.write("Overall Statistics:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total variants: {len(genome_df):,}\n")
            
            if 'AF_imputed' in genome_df.columns and 'AF_reference' in genome_df.columns:
                shared = ((genome_df['AF_imputed'] > 0) & (genome_df['AF_reference'] > 0)).sum()
                imp_only = (genome_df['AF_reference'] == 0).sum()
                ref_only = (genome_df['AF_imputed'] == 0).sum()
                
                f.write(f"Shared variants: {shared:,} ({100*shared/len(genome_df):.1f}%)\n")
                f.write(f"Imputed-only variants: {imp_only:,} ({100*imp_only/len(genome_df):.1f}%)\n")
                f.write(f"Reference-only variants: {ref_only:,} ({100*ref_only/len(genome_df):.1f}%)\n\n")
                
                # Correlation
                valid_data = genome_df[(genome_df['AF_imputed'] > 0) & (genome_df['AF_reference'] > 0)]
                if len(valid_data) > 1:
                    corr = np.corrcoef(valid_data['AF_reference'], valid_data['AF_imputed'])[0, 1]
                    f.write(f"Overall correlation (shared variants): {corr:.4f}\n\n")
            
            # AF difference statistics
            if 'AF_diff' in genome_df.columns:
                af_diff = genome_df['AF_diff'].dropna()
                f.write("Allele Frequency Difference Statistics:\n")
                f.write("-" * 30 + "\n")
                f.write(f"Mean difference: {af_diff.mean():.6f}\n")
                f.write(f"Median difference: {af_diff.median():.6f}\n")
                f.write(f"Standard deviation: {af_diff.std():.6f}\n")
                f.write(f"Min difference: {af_diff.min():.6f}\n")
                f.write(f"Max difference: {af_diff.max():.6f}\n")
                f.write(f"95% CI: [{np.percentile(af_diff, 2.5):.6f}, {np.percentile(af_diff, 97.5):.6f}]\n")
                
                # Categorize differences
                f.write("\nDifference categories:\n")
                f.write(f"  |diff| < 0.01: {(af_diff.abs() < 0.01).sum():,} ({100*(af_diff.abs() < 0.01).sum()/len(af_diff):.1f}%)\n")
                f.write(f"  0.01 ≤ |diff| < 0.05: {((af_diff.abs() >= 0.01) & (af_diff.abs() < 0.05)).sum():,} ({100*((af_diff.abs() >= 0.01) & (af_diff.abs() < 0.05)).sum()/len(af_diff):.1f}%)\n")
                f.write(f"  0.05 ≤ |diff| < 0.10: {((af_diff.abs() >= 0.05) & (af_diff.abs() < 0.10)).sum():,} ({100*((af_diff.abs() >= 0.05) & (af_diff.abs() < 0.10)).sum()/len(af_diff):.1f}%)\n")
                f.write(f"  0.10 ≤ |diff| < 0.20: {((af_diff.abs() >= 0.10) & (af_diff.abs() < 0.20)).sum():,} ({100*((af_diff.abs() >= 0.10) & (af_diff.abs() < 0.20)).sum()/len(af_diff):.1f}%)\n")
                f.write(f"  |diff| ≥ 0.20: {(af_diff.abs() >= 0.20).sum():,} ({100*(af_diff.abs() >= 0.20).sum()/len(af_diff):.1f}%)\n\n")
            
            # Per-chromosome statistics
            if chr_stats:
                f.write("Per-Chromosome Statistics:\n")
                f.write("-" * 70 + "\n")
                f.write(f"{'Chr':<6} {'Total':<10} {'Shared':<10} {'Imp-only':<10} {'Ref-only':<10}")
                if any('mean_af_diff' in s for s in chr_stats):
                    f.write(f" {'Mean Diff':<12} {'Max |Diff|':<10}")
                f.write("\n")
                
                for stats in chr_stats:
                    f.write(f"{stats['chromosome']:<6} ")
                    f.write(f"{stats['n_variants']:<10,} ")
                    f.write(f"{stats['n_shared']:<10,} ")
                    f.write(f"{stats['n_imputed_only']:<10,} ")
                    f.write(f"{stats['n_reference_only']:<10,} ")
                    if 'mean_af_diff' in stats:
                        f.write(f"{stats['mean_af_diff']:<12.6f} ")
                        f.write(f"{stats['max_abs_diff']:<10.4f}")
                    f.write("\n")
            
            # MAF bin analysis
            if 'MAF_imputed' in genome_df.columns and 'MAF_reference' in genome_df.columns:
                f.write("\nMAF Distribution Analysis:\n")
                f.write("-" * 30 + "\n")
                
                maf_bins = [0, 0.01, 0.05, 0.1, 0.2, 0.5]
                bin_labels = ['0-1%', '1-5%', '5-10%', '10-20%', '20-50%']
                
                genome_df['maf_bin_imp'] = pd.cut(genome_df['MAF_imputed'], bins=maf_bins, labels=bin_labels, include_lowest=True)
                genome_df['maf_bin_ref'] = pd.cut(genome_df['MAF_reference'], bins=maf_bins, labels=bin_labels, include_lowest=True)
                
                f.write("\nImputed MAF distribution:\n")
                for label in bin_labels:
                    count = (genome_df['maf_bin_imp'] == label).sum()
                    pct = 100 * count / len(genome_df) if len(genome_df) > 0 else 0
                    f.write(f"  {label}: {count:,} ({pct:.1f}%)\n")
                
                f.write("\nReference MAF distribution:\n")
                for label in bin_labels:
                    count = (genome_df['maf_bin_ref'] == label).sum()
                    pct = 100 * count / len(genome_df) if len(genome_df) > 0 else 0
                    f.write(f"  {label}: {count:,} ({pct:.1f}%)\n")
        
        print(f"Summary saved to {args.output_prefix}.genome_freq_summary.txt")
        
    else:
        # Create empty output if no data
        pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 
                              'AF_imputed', 'MAF_imputed', 
                              'AF_reference', 'MAF_reference',
                              'AF_diff', 'MAF_diff']).to_csv(
            f"{args.output_prefix}.genome_freq.tsv", sep='\t', index=False
        )
        
        with open(f"{args.output_prefix}.genome_freq_summary.txt", 'w') as f:
            f.write("No frequency data to combine\n")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"COMBINE_FREQ_GENOME":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    pandas: {pd.__version__}\n')
        f.write(f'    numpy: {np.__version__}\n')

if __name__ == "__main__":
    main()