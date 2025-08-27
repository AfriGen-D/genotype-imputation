#!/usr/bin/env python3
"""
Convert Minimac4 sites VCF file to tab-delimited info format
"""

import sys
import gzip
import argparse

def convert_vcf_to_info(vcf_file, output_file):
    """Convert VCF sites file to Minimac4 info format"""
    
    # Determine if input is gzipped
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    with open_func(vcf_file, mode) as infile, open(output_file, 'w') as outfile:
        # Write header for info format
        outfile.write("SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n")
        
        # Skip VCF headers
        for line in infile:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                continue
                
            # Parse VCF data lines
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
                
            chrom = parts[0]
            pos = parts[1]
            variant_id = parts[2] if parts[2] != '.' else f"{chrom}:{pos}"
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = 'TRUE'
            
            # Extract values from INFO field
            snp = f"{chrom}:{pos}:{ref}:{alt}"
            alt_freq = info_dict.get('AF', '0')
            maf = info_dict.get('MAF', '0')
            avg_call = info_dict.get('AVG_CS', '0')
            rsq = info_dict.get('R2', '0')
            genotyped = '1' if 'TYPED' in info_dict else '0'
            loo_rsq = info_dict.get('ER2', '-')
            emp_r = '-'
            emp_rsq = info_dict.get('ER2', '-')
            dose0 = '-'
            dose1 = '-'
            
            # Write info format line
            outfile.write(f"{snp}\t{ref}\t{alt}\t{alt_freq}\t{maf}\t{avg_call}\t{rsq}\t{genotyped}\t{loo_rsq}\t{emp_r}\t{emp_rsq}\t{dose0}\t{dose1}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert VCF sites file to info format')
    parser.add_argument('vcf_file', help='Input VCF sites file')
    parser.add_argument('output_file', help='Output info file')
    
    args = parser.parse_args()
    
    convert_vcf_to_info(args.vcf_file, args.output_file)
    print(f"Converted {args.vcf_file} to {args.output_file}")

if __name__ == "__main__":
    main()