#!/usr/bin/env python3

import sys
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate accuracy report')
    parser.add_argument('--acc-info-file', required=True, help='Accuracy INFO file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    report_file = f"{args.output_prefix}_{args.ref_name}.accuracy.txt"
    tsv_file = f"{args.output_prefix}_{args.ref_name}.accuracy.tsv"
    
    # Calculate accuracy metrics by MAF bins
    maf_bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]
    bin_metrics = {bin_range: {'count': 0, 'rsq_sum': 0} for bin_range in maf_bins}
    
    with open(args.acc_info_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    try:
                        maf = float(parts[3])  # MAF is in column 4 (index 3)
                        rsq = float(parts[5])  # Rsq is in column 6 (index 5)
                        
                        for bin_range in maf_bins:
                            if bin_range[0] <= maf < bin_range[1]:
                                bin_metrics[bin_range]['count'] += 1
                                bin_metrics[bin_range]['rsq_sum'] += rsq
                                break
                    except (ValueError, IndexError):
                        continue
    
    # Write detailed report
    with open(report_file, 'w') as f:
        f.write("Accuracy Report\n")
        f.write("="*50 + "\n")
        f.write(f"Sample: {args.sample_id}\n")
        f.write(f"Reference: {args.ref_name}\n")
        f.write("="*50 + "\n\n")
        
        f.write("MAF Bin\tCount\tMean Rsq\n")
        for bin_range, metrics in sorted(bin_metrics.items()):
            mean_rsq = metrics['rsq_sum'] / metrics['count'] if metrics['count'] > 0 else 0
            f.write(f"{bin_range[0]:.2f}-{bin_range[1]:.2f}\t{metrics['count']}\t{mean_rsq:.4f}\n")
    
    # Write TSV for plotting
    with open(tsv_file, 'w') as f:
        f.write("MAF_BIN\tCOUNT\tMEAN_RSQ\n")
        for bin_range, metrics in sorted(bin_metrics.items()):
            mean_rsq = metrics['rsq_sum'] / metrics['count'] if metrics['count'] > 0 else 0
            f.write(f"{bin_range[0]:.2f}-{bin_range[1]:.2f}\t{metrics['count']}\t{mean_rsq:.4f}\n")
    
    print(f"Accuracy report: {report_file}")
    print(f"Accuracy TSV: {tsv_file}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"REPORT_ACCURACY":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')

if __name__ == "__main__":
    main()