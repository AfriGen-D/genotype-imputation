#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge and compare imputed and reference frequencies')
    parser.add_argument('imputed_freq', help='Imputed frequency file')
    parser.add_argument('ref_freq', help='Reference frequency file')
    parser.add_argument('--prefix', required=True, help='Output prefix')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chunk-id', default='', help='Chunk identifier')
    args = parser.parse_args()
    
    # Read imputed frequencies
    imp_df = pd.read_csv(args.imputed_freq, sep='\t')
    
    # Read reference frequencies
    ref_df = pd.read_csv(args.ref_freq, sep='\t')
    
    # Create a key for merging (CHROM:POS:REF:ALT)
    imp_df['key'] = imp_df['CHROM'].astype(str) + ':' + imp_df['POS'].astype(str) + ':' + \
                     imp_df['REF'] + ':' + imp_df['ALT']
    
    if not ref_df.empty:
        ref_df['key'] = ref_df['CHROM'].astype(str) + ':' + ref_df['POS'].astype(str) + ':' + \
                        ref_df['REF'] + ':' + ref_df['ALT']
        
        # Merge on the key
        merged_df = pd.merge(imp_df, ref_df[['key', 'REF_AF', 'REF_MAF']], 
                            on='key', how='left')
        
        # Fill missing reference values with NA
        merged_df['REF_AF'] = merged_df['REF_AF'].fillna('NA')
        merged_df['REF_MAF'] = merged_df['REF_MAF'].fillna('NA')
        
        # Calculate differences where both values are available
        mask = (merged_df['REF_AF'] != 'NA')
        merged_df['AF_DIFF'] = np.where(mask, 
                                        merged_df['IMP_AF'].astype(float) - merged_df['REF_AF'].astype(float),
                                        'NA')
        merged_df['MAF_DIFF'] = np.where(mask,
                                         merged_df['IMP_MAF'].astype(float) - merged_df['REF_MAF'].astype(float),
                                         'NA')
    else:
        # No reference data available
        merged_df = imp_df.copy()
        merged_df['REF_AF'] = 'NA'
        merged_df['REF_MAF'] = 'NA'
        merged_df['AF_DIFF'] = 'NA'
        merged_df['MAF_DIFF'] = 'NA'
    
    # Select and reorder columns for output
    output_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 
                   'IMP_AF', 'IMP_MAF', 'IMP_R2',
                   'REF_AF', 'REF_MAF', 
                   'AF_DIFF', 'MAF_DIFF']
    
    final_df = merged_df[output_cols]
    
    # Write comparison file
    chunk_suffix = f"_{args.chunk_id}" if args.chunk_id else ""
    output_file = f"{args.prefix}_{args.ref_name}{chunk_suffix}.freq_comparison.tsv"
    final_df.to_csv(output_file, sep='\t', index=False)
    
    # Generate summary statistics
    summary_file = f"{args.prefix}_{args.ref_name}{chunk_suffix}.freq_comparison.summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Frequency Comparison Summary\n")
        f.write(f"============================\n")
        f.write(f"Total variants in imputed data: {len(imp_df)}\n")
        
        if not ref_df.empty:
            f.write(f"Total variants in reference panel: {len(ref_df)}\n")
            
            # Count matched variants
            matched = merged_df[merged_df['REF_AF'] != 'NA']
            f.write(f"Variants found in both: {len(matched)}\n")
            f.write(f"Variants only in imputed: {len(merged_df) - len(matched)}\n")
            
            if len(matched) > 0:
                # Calculate correlation if there are matched variants
                imp_af = matched['IMP_AF'].astype(float)
                ref_af = matched['REF_AF'].astype(float)
                correlation = np.corrcoef(imp_af, ref_af)[0, 1]
                
                # Calculate mean differences
                af_diff = matched['AF_DIFF'].astype(float)
                maf_diff = matched['MAF_DIFF'].astype(float)
                
                f.write(f"\nFrequency Comparison Metrics:\n")
                f.write(f"AF Correlation: {correlation:.4f}\n")
                f.write(f"Mean AF difference: {af_diff.mean():.4f} (SD: {af_diff.std():.4f})\n")
                f.write(f"Mean MAF difference: {maf_diff.mean():.4f} (SD: {maf_diff.std():.4f})\n")
                f.write(f"Max AF difference: {af_diff.abs().max():.4f}\n")
                f.write(f"Max MAF difference: {maf_diff.abs().max():.4f}\n")
        else:
            f.write(f"No reference panel data available for comparison\n")
    
    print(f"Merged frequency comparison for {len(final_df)} variants")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"MERGE_FREQ_COMPARISON":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    pandas: {pd.__version__}\n')
        f.write(f'    numpy: {np.__version__}\n')

if __name__ == "__main__":
    main()