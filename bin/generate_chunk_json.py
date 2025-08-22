#!/usr/bin/env python3
"""
Generate JSON summary from chunk-level reports for aggregation.
"""

import json
import argparse
import sys
from pathlib import Path
import re

def parse_accuracy_file(accuracy_file):
    """Parse accuracy report file and extract metrics."""
    metrics = {
        'maf_bins': {},
        'r2_by_maf': {}
    }
    
    try:
        with open(accuracy_file, 'r') as f:
            lines = f.readlines()
            
        # Find the data section
        in_data = False
        for line in lines:
            if 'MAF Bin' in line or 'MAF_BIN' in line:
                in_data = True
                continue
            
            if in_data and line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0][0].isdigit():
                    bin_range = parts[0]
                    count = int(parts[1])
                    mean_rsq = float(parts[2])
                    
                    metrics['maf_bins'][bin_range] = count
                    metrics['r2_by_maf'][bin_range] = mean_rsq
                    
    except Exception as e:
        print(f"Warning: Could not parse accuracy file: {e}", file=sys.stderr)
    
    return metrics

def parse_well_imputed_file(well_file):
    """Parse well imputed report file and extract metrics."""
    metrics = {
        'well_imputed_by_maf': {},
        'total_well_imputed': 0
    }
    
    try:
        with open(well_file, 'r') as f:
            lines = f.readlines()
            
        # Find the data section
        in_data = False
        for line in lines:
            if 'MAF_BIN' in line:
                in_data = True
                continue
            
            if in_data and line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[0][0].isdigit():
                    bin_range = parts[0]
                    count = int(parts[1])
                    metrics['well_imputed_by_maf'][bin_range] = count
                    metrics['total_well_imputed'] += count
                    
    except Exception as e:
        print(f"Warning: Could not parse well imputed file: {e}", file=sys.stderr)
    
    return metrics

def parse_summary_file(summary_file):
    """Parse summary file for total counts."""
    metrics = {}
    
    try:
        with open(summary_file, 'r') as f:
            for line in f:
                if 'Total well imputed variants:' in line:
                    metrics['total_well_imputed'] = int(line.split(':')[1].strip())
                elif 'Reference panel:' in line:
                    metrics['ref_name'] = line.split(':')[1].strip()
                elif 'Sample:' in line:
                    metrics['sample_id'] = line.split(':')[1].strip()
    except Exception as e:
        print(f"Warning: Could not parse summary file: {e}", file=sys.stderr)
    
    return metrics

def extract_chunk_info(chunk_id):
    """Extract chromosome and position info from chunk ID."""
    info = {
        'chromosome': None,
        'start': None,
        'end': None
    }
    
    # Pattern: dataset_chr##_start_end
    match = re.search(r'chr(\d+)_(\d+)_(\d+)', chunk_id)
    if match:
        info['chromosome'] = f"chr{match.group(1)}"
        info['start'] = int(match.group(2))
        info['end'] = int(match.group(3))
    
    return info

def main():
    parser = argparse.ArgumentParser(description='Generate JSON summary from chunk reports')
    parser.add_argument('--accuracy-file', required=True, help='Accuracy report file')
    parser.add_argument('--well-imputed-file', required=True, help='Well imputed report file')
    parser.add_argument('--summary-file', required=True, help='Summary file')
    parser.add_argument('--chunk-id', required=True, help='Chunk identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--output', required=True, help='Output JSON file')
    
    args = parser.parse_args()
    
    # Parse all input files
    accuracy_metrics = parse_accuracy_file(args.accuracy_file)
    well_imputed_metrics = parse_well_imputed_file(args.well_imputed_file)
    summary_metrics = parse_summary_file(args.summary_file)
    chunk_info = extract_chunk_info(args.chunk_id)
    
    # Calculate totals
    total_variants = sum(accuracy_metrics['maf_bins'].values())
    total_genotyped = 0  # Would need to extract from VCF
    total_imputed = total_variants - total_genotyped
    
    # Calculate mean R2
    r2_sum = 0
    count_sum = 0
    for bin_range, count in accuracy_metrics['maf_bins'].items():
        if bin_range in accuracy_metrics['r2_by_maf'] and count > 0:
            r2_sum += accuracy_metrics['r2_by_maf'][bin_range] * count
            count_sum += count
    
    mean_r2 = r2_sum / count_sum if count_sum > 0 else 0
    
    # Create comprehensive JSON output
    json_output = {
        'chunk_id': args.chunk_id,
        'ref_name': args.ref_name,
        'chromosome': chunk_info['chromosome'],
        'start': chunk_info['start'],
        'end': chunk_info['end'],
        'total_variants': total_variants,
        'well_imputed_variants': well_imputed_metrics.get('total_well_imputed', 0),
        'total_genotyped': total_genotyped,
        'total_imputed': total_imputed,
        'mean_r2': mean_r2,
        'maf_distribution': accuracy_metrics['maf_bins'],
        'r2_by_maf': accuracy_metrics['r2_by_maf'],
        'well_imputed_by_maf': well_imputed_metrics['well_imputed_by_maf'],
        'concordance': None  # Would need additional data
    }
    
    # Write JSON output
    with open(args.output, 'w') as f:
        json.dump(json_output, f, indent=2)
    
    print(f"JSON summary written to {args.output}")
    print(f"Total variants: {total_variants}")
    print(f"Well imputed: {well_imputed_metrics.get('total_well_imputed', 0)}")
    print(f"Mean R²: {mean_r2:.4f}")

if __name__ == "__main__":
    main()