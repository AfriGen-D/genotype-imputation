#!/usr/bin/env python3

import gzip
import sys
import argparse

def parse_info_field(info_str):
    """Parse VCF INFO field and extract key-value pairs"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            # Flag fields without values
            info_dict[item] = True
    return info_dict

def main():
    parser = argparse.ArgumentParser(description='Filter VCF INFO file by R² threshold')
    parser.add_argument('--input-file', required=True, help='Input VCF file (can be gzipped)')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--r2-threshold', type=float, default=0.3, help='R² threshold for well imputed variants')
    
    args = parser.parse_args()
    
    well_imputed_file = f"{args.output_prefix}_{args.ref_name}.filtered.info"
    accuracy_file = f"{args.output_prefix}_{args.ref_name}.acc.info"
    
    # Read and filter VCF file
    well_imputed = []
    accuracy = []
    header_lines = []
    
    open_func = gzip.open if args.input_file.endswith('.gz') else open
    mode = 'rt' if args.input_file.endswith('.gz') else 'r'
    
    with open_func(args.input_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                continue
            
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 8:  # VCF has at least 8 columns
                    info_field = parts[7]
                    info_dict = parse_info_field(info_field)
                    
                    # Extract R2 value
                    if 'R2' in info_dict:
                        try:
                            r2_value = float(info_dict['R2'])
                            
                            # Create tab-delimited info format line
                            # Format: SNP REF ALT MAF AvgCall Rsq Genotyped/Imputed
                            chrom = parts[0]
                            pos = parts[1]
                            snp_id = f"{chrom}:{pos}"
                            ref = parts[3]
                            alt = parts[4]
                            maf = info_dict.get('MAF', '0')
                            avg_cs = info_dict.get('AVG_CS', '1')
                            is_typed = 'Genotyped' if 'TYPED' in info_dict else 'Imputed'
                            
                            info_line = f"{snp_id}\t{ref}\t{alt}\t{maf}\t{avg_cs}\t{r2_value:.6f}\t{is_typed}\n"
                            
                            if r2_value >= args.r2_threshold:
                                well_imputed.append(info_line)
                            accuracy.append(info_line)
                        except (ValueError, KeyError) as e:
                            print(f"Warning: Could not parse R2 value for line: {line.strip()}", file=sys.stderr)
                            continue
    
    # Write header for info format
    header = "SNP\tREF\tALT\tMAF\tAvgCall\tRsq\tGenotyped\n"
    
    # Write filtered files
    with open(well_imputed_file, 'w') as f:
        f.write(header)
        f.writelines(well_imputed)
    
    with open(accuracy_file, 'w') as f:
        f.write(header)
        f.writelines(accuracy)
    
    print(f"Well imputed variants (R² >= {args.r2_threshold}): {len(well_imputed)}")
    print(f"Total variants for accuracy: {len(accuracy)}")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"FILTER_INFO_BY_TARGET":\n')
        f.write(f'    python: {sys.version.split()[0]}\n')

if __name__ == "__main__":
    main()