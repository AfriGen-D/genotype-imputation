#!/usr/bin/env python3

import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Filter INFO file by R² threshold')
    parser.add_argument('--input-file', required=True, help='Input INFO file (can be gzipped)')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--r2-threshold', type=float, default=0.3, help='R² threshold for well imputed variants')
    
    args = parser.parse_args()
    
    well_imputed_file = f"{args.output_prefix}_{args.ref_name}.filtered.info"
    accuracy_file = f"{args.output_prefix}_{args.ref_name}.acc.info"
    
    # Read and filter info file
    well_imputed = []
    accuracy = []
    
    with gzip.open(args.input_file, 'rt') if args.input_file.endswith('.gz') else open(args.input_file, 'r') as f:
        header = f.readline()
        
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    try:
                        rsq = float(parts[6])  # Assuming Rsq is in column 7
                        if rsq >= args.r2_threshold:
                            well_imputed.append(line)
                        accuracy.append(line)
                    except (ValueError, IndexError):
                        continue
    
    # Write filtered files
    with open(well_imputed_file, 'w') as f:
        f.write(header)
        f.writelines(well_imputed)
    
    with open(accuracy_file, 'w') as f:
        f.write(header)
        f.writelines(accuracy)
    
    print(f"Well imputed variants (Rsq >= {args.r2_threshold}): {len(well_imputed)}")
    print(f"Total variants for accuracy: {len(accuracy)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"FILTER_INFO_BY_TARGET":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')

if __name__ == "__main__":
    main()