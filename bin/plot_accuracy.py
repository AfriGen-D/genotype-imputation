#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot accuracy metrics')
    parser.add_argument('--accuracy-tsv', required=True, help='Accuracy TSV file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    # Read accuracy TSV
    data = pd.read_csv(args.accuracy_tsv, sep='\t')
    
    # Create accuracy plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Mean Rsq by MAF
    if 'MAF_BIN' in data.columns and 'MEAN_RSQ' in data.columns:
        ax1.plot(range(len(data)), data['MEAN_RSQ'], 'o-', linewidth=2, markersize=8)
        ax1.set_xticks(range(len(data)))
        ax1.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax1.set_xlabel('MAF Bin')
        ax1.set_ylabel('Mean Rsq')
        ax1.set_title(f'Imputation Accuracy by MAF')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([0, 1])
    
    # Plot 2: Count distribution
    if 'MAF_BIN' in data.columns and 'COUNT' in data.columns:
        ax2.bar(range(len(data)), data['COUNT'], alpha=0.7)
        ax2.set_xticks(range(len(data)))
        ax2.set_xticklabels(data['MAF_BIN'], rotation=45)
        ax2.set_xlabel('MAF Bin')
        ax2.set_ylabel('Number of Variants')
        ax2.set_title(f'Variant Distribution by MAF')
        ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f'Imputation Accuracy - {args.sample_id} ({args.ref_name})', fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_{args.ref_name}.accuracy.png", dpi=150)
    plt.close()
    
    print(f"Accuracy plot saved: {args.output_prefix}_{args.ref_name}.accuracy.png")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"PLOT_ACCURACY":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\n')
        f.write(f'    pandas: {pd.__version__}\n')

if __name__ == "__main__":
    main()