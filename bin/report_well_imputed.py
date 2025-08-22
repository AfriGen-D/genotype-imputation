#!/usr/bin/env python3

import sys
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate well imputed variants report')
    parser.add_argument('--well-info-file', required=True, help='Well imputed INFO file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    
    args = parser.parse_args()
    
    report_file = f"{args.output_prefix}_{args.ref_name}.well_imputed.txt"
    summary_file = f"{args.output_prefix}_{args.ref_name}.well_imputed_summary.txt"
    
    # Process well imputed variants
    maf_bins = [(0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]
    bin_counts = {bin_range: 0 for bin_range in maf_bins}
    total_variants = 0
    
    with open(args.well_info_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    try:
                        maf = float(parts[3])  # MAF is in column 4 (index 3)
                        total_variants += 1
                        
                        for bin_range in maf_bins:
                            if bin_range[0] <= maf < bin_range[1]:
                                bin_counts[bin_range] += 1
                                break
                    except (ValueError, IndexError):
                        continue
    
    # Write detailed report
    with open(report_file, 'w') as f:
        f.write("MAF_BIN\tCOUNT\tPERCENTAGE\n")
        for bin_range, count in sorted(bin_counts.items()):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            f.write(f"{bin_range[0]:.2f}-{bin_range[1]:.2f}\t{count}\t{percentage:.2f}\n")
    
    # Write summary
    with open(summary_file, 'w') as f:
        f.write(f"Total well imputed variants: {total_variants}\n")
        f.write(f"Reference panel: {args.ref_name}\n")
        f.write(f"Sample: {args.sample_id}\n")
    
    print(f"Report generated: {report_file}")
    print(f"Summary generated: {summary_file}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"REPORT_WELL_IMPUTED":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')

if __name__ == "__main__":
    main()