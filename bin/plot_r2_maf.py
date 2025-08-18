#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Plot R² vs MAF scatter plot')
    parser.add_argument('--acc-info-file', required=True, help='Accuracy INFO file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    # Read accuracy info file
    mafs = []
    rsqs = []
    
    with open(args.acc_info_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    try:
                        maf = float(parts[4])
                        rsq = float(parts[6])
                        mafs.append(maf)
                        rsqs.append(rsq)
                    except (ValueError, IndexError):
                        continue
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if mafs and rsqs:
        # Create hexbin plot for large datasets
        if len(mafs) > 1000:
            hb = ax.hexbin(mafs, rsqs, gridsize=50, cmap='YlOrRd', mincnt=1)
            cb = plt.colorbar(hb)
            cb.set_label('Count')
        else:
            ax.scatter(mafs, rsqs, alpha=0.5, s=10)
        
        # Add reference lines
        ax.axhline(y=0.3, color='r', linestyle='--', alpha=0.5, label='R² = 0.3')
        ax.axhline(y=0.8, color='g', linestyle='--', alpha=0.5, label='R² = 0.8')
        
        ax.set_xlabel('Minor Allele Frequency (MAF)')
        ax.set_ylabel('Imputation Quality (R²)')
        ax.set_title(f'R² vs MAF - {args.sample_id} ({args.ref_name})')
        ax.set_xlim([0, 0.5])
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_{args.ref_name}.r2_maf.png", dpi=150)
    plt.close()
    
    print(f"R² vs MAF plot saved: {args.output_prefix}_{args.ref_name}.r2_maf.png")
    print(f"Total variants plotted: {len(mafs)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"PLOT_R2_MAF":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\n')

if __name__ == "__main__":
    main()