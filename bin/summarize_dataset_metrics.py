#!/usr/bin/env python3
"""
Summarize metrics across all chunks for a dataset
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import json
import re


def extract_metrics_from_report(report_file):
    """Extract key metrics from a combined report file"""
    metrics = {
        'total_variants': 0,
        'well_imputed_variants': 0,
        'avg_r2': [],
        'avg_maf': [],
        'chunks_processed': 0
    }
    
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            
            # Extract number of chunks
            chunks_match = re.search(r'Number of chunks:\s*(\d+)', content)
            if chunks_match:
                metrics['chunks_processed'] = int(chunks_match.group(1))
            
            # Extract R² values
            r2_matches = re.findall(r'(?:Average R2|Mean R²):\s*([\d.]+)', content, re.IGNORECASE)
            metrics['avg_r2'] = [float(x) for x in r2_matches if x]
            
            # Extract MAF values
            maf_matches = re.findall(r'(?:Average MAF|Mean MAF):\s*([\d.]+)', content, re.IGNORECASE)
            metrics['avg_maf'] = [float(x) for x in maf_matches if x]
            
            # Extract variant counts
            variant_matches = re.findall(r'Total variants:\s*(\d+)', content, re.IGNORECASE)
            if variant_matches:
                metrics['total_variants'] = sum(int(x) for x in variant_matches)
            
            well_imputed_matches = re.findall(r'Well imputed variants:\s*(\d+)', content, re.IGNORECASE)
            if well_imputed_matches:
                metrics['well_imputed_variants'] = sum(int(x) for x in well_imputed_matches)
                
    except Exception as e:
        print(f"Warning: Could not fully parse {report_file}: {e}")
    
    return metrics


def main():
    parser = argparse.ArgumentParser(description='Summarize dataset-level metrics from chunk reports')
    parser.add_argument('combined_report', help='Combined report file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for summary files')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--dataset-id', required=True, help='Dataset identifier')
    
    args = parser.parse_args()
    
    # Extract metrics from combined report
    metrics = extract_metrics_from_report(args.combined_report)
    
    # Calculate summary statistics
    summary = {
        'dataset_id': args.dataset_id,
        'reference_panel': args.ref_name,
        'chunks_processed': metrics['chunks_processed'],
        'total_variants': metrics['total_variants'],
        'well_imputed_variants': metrics['well_imputed_variants'],
        'imputation_rate': (metrics['well_imputed_variants'] / metrics['total_variants'] * 100) 
                          if metrics['total_variants'] > 0 else 0,
        'mean_r2': np.mean(metrics['avg_r2']) if metrics['avg_r2'] else 0,
        'median_r2': np.median(metrics['avg_r2']) if metrics['avg_r2'] else 0,
        'min_r2': np.min(metrics['avg_r2']) if metrics['avg_r2'] else 0,
        'max_r2': np.max(metrics['avg_r2']) if metrics['avg_r2'] else 0,
        'mean_maf': np.mean(metrics['avg_maf']) if metrics['avg_maf'] else 0
    }
    
    # Write summary to JSON
    json_file = f"{args.output_prefix}_{args.ref_name}.summary.json"
    with open(json_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Write summary to text file
    txt_file = f"{args.output_prefix}_{args.ref_name}.summary.txt"
    with open(txt_file, 'w') as f:
        f.write("Dataset-Level Summary Metrics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Dataset ID: {summary['dataset_id']}\n")
        f.write(f"Reference Panel: {summary['reference_panel']}\n")
        f.write(f"Chunks Processed: {summary['chunks_processed']}\n\n")
        
        f.write("Variant Statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total Variants: {summary['total_variants']:,}\n")
        f.write(f"Well Imputed Variants: {summary['well_imputed_variants']:,}\n")
        f.write(f"Imputation Rate: {summary['imputation_rate']:.2f}%\n\n")
        
        f.write("R² Statistics:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Mean R²: {summary['mean_r2']:.4f}\n")
        f.write(f"Median R²: {summary['median_r2']:.4f}\n")
        f.write(f"Min R²: {summary['min_r2']:.4f}\n")
        f.write(f"Max R²: {summary['max_r2']:.4f}\n\n")
        
        f.write(f"Mean MAF: {summary['mean_maf']:.4f}\n")
    
    print(f"Summary metrics written to {json_file} and {txt_file}")


if __name__ == '__main__':
    main()