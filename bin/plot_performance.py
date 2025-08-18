#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot performance metrics')
    parser.add_argument('--well-imputed-file', required=True, help='Well imputed variants file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    # Read the well imputed report
    data = pd.read_csv(args.well_imputed_file, sep='\t')
    
    # Create performance plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if 'MAF_BIN' in data.columns and 'COUNT' in data.columns:
        ax.bar(range(len(data)), data['COUNT'])
        ax.set_xticks(range(len(data)))
        ax.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax.set_xlabel('MAF Bin')
        ax.set_ylabel('Number of Well Imputed Variants')
        ax.set_title(f'Imputation Performance by MAF - {args.sample_id}')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_{args.ref_name}.performance.png", dpi=150)
    plt.close()
    
    print(f"Performance plot saved: {args.output_prefix}_{args.ref_name}.performance.png")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"PLOT_PERFORMANCE":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\n')
        f.write(f'    pandas: {pd.__version__}\n')

if __name__ == "__main__":
    main()