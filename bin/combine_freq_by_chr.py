#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='Combine frequency files by chromosome')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    parser.add_argument('--chromosome', required=True, help='Chromosome identifier')
    parser.add_argument('--pattern', default='*.freq.tsv', help='File pattern for frequency files (default: *.freq.tsv)')
    args = parser.parse_args()
    
    # Read all frequency files for this chromosome
    freq_files = glob.glob(args.pattern)
    
    # Separate imputed and reference files
    imputed_files = [f for f in freq_files if 'imputed' in f]
    ref_files = [f for f in freq_files if 'reference' in f]
    
    # Combine imputed frequencies
    imputed_dfs = []
    for f in imputed_files:
        df = pd.read_csv(f, sep='\t')
        imputed_dfs.append(df)
    
    if imputed_dfs:
        imputed_combined = pd.concat(imputed_dfs, ignore_index=True)
        # Remove duplicates, keeping first occurrence
        imputed_combined = imputed_combined.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])
    else:
        imputed_combined = pd.DataFrame()
    
    # Combine reference frequencies
    ref_dfs = []
    for f in ref_files:
        df = pd.read_csv(f, sep='\t')
        ref_dfs.append(df)
    
    if ref_dfs:
        ref_combined = pd.concat(ref_dfs, ignore_index=True)
        # Remove duplicates, keeping first occurrence
        ref_combined = ref_combined.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])
    else:
        ref_combined = pd.DataFrame()
    
    # Merge imputed and reference on position
    if not imputed_combined.empty and not ref_combined.empty:
        merged = pd.merge(
            imputed_combined[['CHROM', 'POS', 'REF', 'ALT', 'AF', 'MAF']],
            ref_combined[['CHROM', 'POS', 'REF', 'ALT', 'AF', 'MAF']],
            on=['CHROM', 'POS', 'REF', 'ALT'],
            suffixes=('_imputed', '_reference'),
            how='outer'
        )
        
        # Fill missing values with 0
        merged = merged.fillna(0)
        
        # Calculate frequency difference
        merged['AF_diff'] = merged['AF_imputed'] - merged['AF_reference']
        merged['MAF_diff'] = merged['MAF_imputed'] - merged['MAF_reference']
        
        # Save combined frequencies
        output_file = f"{args.prefix}_{args.chromosome}.chr_freq.tsv"
        merged.to_csv(output_file, sep='\t', index=False)
        
        print(f"Combined {len(merged)} variants for chromosome {args.chromosome}")
        print(f"Imputed-only variants: {merged['AF_reference'].eq(0).sum()}")
        print(f"Reference-only variants: {merged['AF_imputed'].eq(0).sum()}")
        print(f"Shared variants: {(merged['AF_imputed'] > 0).sum() - merged['AF_reference'].eq(0).sum()}")
    else:
        # Create empty output if no data
        output_file = f"{args.prefix}_{args.chromosome}.chr_freq.tsv"
        pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 
                              'AF_imputed', 'MAF_imputed', 
                              'AF_reference', 'MAF_reference',
                              'AF_diff', 'MAF_diff']).to_csv(
            output_file, sep='\t', index=False
        )
        print("No frequency data to combine")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"COMBINE_FREQ_BY_CHR":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    pandas: {pd.__version__}\n')
        f.write(f'    numpy: {np.__version__}\n')

if __name__ == "__main__":
    main()