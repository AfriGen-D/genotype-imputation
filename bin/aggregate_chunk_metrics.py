#!/usr/bin/env python3
"""
Efficiently aggregate metrics from many chunk reports using streaming approach
Designed to handle 1000+ chunks without memory issues
"""

import argparse
import json
import numpy as np
from pathlib import Path
import re
from collections import defaultdict
import sys


class StreamingMetricsAggregator:
    """Aggregate metrics from multiple chunks using streaming approach"""
    
    def __init__(self):
        self.metrics = {
            'chunks_processed': 0,
            'total_variants': 0,
            'well_imputed_variants': 0,
            'r2_sum': 0.0,
            'r2_count': 0,
            'r2_values': [],  # Only store if < 100 chunks for percentiles
            'maf_sum': 0.0,
            'maf_count': 0,
            'variants_by_r2_bin': defaultdict(int),
            'variants_by_maf_bin': defaultdict(int),
            'failed_chunks': []
        }
        self.store_raw_values = True  # Will disable for > 100 chunks
    
    def process_chunk_file(self, filepath):
        """Process a single chunk file and update running statistics"""
        try:
            # Determine file type and extract metrics accordingly
            if filepath.suffix == '.json':
                self._process_json_file(filepath)
            elif filepath.suffix in ['.txt', '.tsv']:
                self._process_text_file(filepath)
            
            self.metrics['chunks_processed'] += 1
            
            # Disable raw value storage if too many chunks
            if self.metrics['chunks_processed'] > 100:
                self.store_raw_values = False
                self.metrics['r2_values'] = []  # Clear to save memory
                
        except Exception as e:
            self.metrics['failed_chunks'].append({
                'file': str(filepath),
                'error': str(e)
            })
    
    def _process_json_file(self, filepath):
        """Extract metrics from JSON format chunk report"""
        with open(filepath, 'r') as f:
            data = json.load(f)
            
            # Update counts
            if 'total_variants' in data:
                self.metrics['total_variants'] += data['total_variants']
            if 'well_imputed_variants' in data:
                self.metrics['well_imputed_variants'] += data['well_imputed_variants']
            
            # Update R² statistics
            if 'mean_r2' in data and data['mean_r2'] is not None:
                self.metrics['r2_sum'] += data['mean_r2']
                self.metrics['r2_count'] += 1
                if self.store_raw_values:
                    self.metrics['r2_values'].append(data['mean_r2'])
            
            # Update MAF statistics
            if 'mean_maf' in data and data['mean_maf'] is not None:
                self.metrics['maf_sum'] += data['mean_maf']
                self.metrics['maf_count'] += 1
    
    def _process_text_file(self, filepath):
        """Extract metrics from text format chunk report"""
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Extract total variants
            variant_match = re.search(r'Total variants:\s*(\d+)', content, re.IGNORECASE)
            if variant_match:
                self.metrics['total_variants'] += int(variant_match.group(1))
            
            # Extract well imputed variants
            well_match = re.search(r'Well imputed variants:\s*(\d+)', content, re.IGNORECASE)
            if well_match:
                self.metrics['well_imputed_variants'] += int(well_match.group(1))
            
            # Extract R² value
            r2_match = re.search(r'(?:Average R2|Mean R²):\s*([\d.]+)', content, re.IGNORECASE)
            if r2_match:
                r2_val = float(r2_match.group(1))
                self.metrics['r2_sum'] += r2_val
                self.metrics['r2_count'] += 1
                if self.store_raw_values:
                    self.metrics['r2_values'].append(r2_val)
            
            # Extract MAF value
            maf_match = re.search(r'(?:Average MAF|Mean MAF):\s*([\d.]+)', content, re.IGNORECASE)
            if maf_match:
                maf_val = float(maf_match.group(1))
                self.metrics['maf_sum'] += maf_val
                self.metrics['maf_count'] += 1
    
    def get_summary(self):
        """Calculate final summary statistics"""
        summary = {
            'chunks_processed': self.metrics['chunks_processed'],
            'chunks_failed': len(self.metrics['failed_chunks']),
            'total_variants': self.metrics['total_variants'],
            'well_imputed_variants': self.metrics['well_imputed_variants'],
            'imputation_rate': 0.0,
            'mean_r2': 0.0,
            'median_r2': None,
            'std_r2': None,
            'percentile_25_r2': None,
            'percentile_75_r2': None,
            'mean_maf': 0.0
        }
        
        # Calculate imputation rate
        if summary['total_variants'] > 0:
            summary['imputation_rate'] = (summary['well_imputed_variants'] / 
                                         summary['total_variants'] * 100)
        
        # Calculate R² statistics
        if self.metrics['r2_count'] > 0:
            summary['mean_r2'] = self.metrics['r2_sum'] / self.metrics['r2_count']
            
            # Calculate percentiles if we have raw values
            if self.metrics['r2_values']:
                summary['median_r2'] = np.median(self.metrics['r2_values'])
                summary['std_r2'] = np.std(self.metrics['r2_values'])
                summary['percentile_25_r2'] = np.percentile(self.metrics['r2_values'], 25)
                summary['percentile_75_r2'] = np.percentile(self.metrics['r2_values'], 75)
        
        # Calculate MAF statistics
        if self.metrics['maf_count'] > 0:
            summary['mean_maf'] = self.metrics['maf_sum'] / self.metrics['maf_count']
        
        return summary


def main():
    parser = argparse.ArgumentParser(
        description='Efficiently aggregate metrics from many chunk reports'
    )
    parser.add_argument('--chunk-reports', nargs='+', required=True,
                       help='List of chunk report files')
    parser.add_argument('--output-prefix', required=True,
                       help='Output prefix for summary files')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset-id', required=True,
                       help='Dataset identifier')
    parser.add_argument('--max-chunks-detail', type=int, default=100,
                       help='Maximum chunks to store detailed metrics for (default: 100)')
    
    args = parser.parse_args()
    
    # Initialize aggregator
    aggregator = StreamingMetricsAggregator()
    
    # Process each chunk file
    print(f"Processing {len(args.chunk_reports)} chunk reports...")
    for i, report_file in enumerate(args.chunk_reports):
        if i % 100 == 0:
            print(f"  Processed {i}/{len(args.chunk_reports)} chunks...")
        aggregator.process_chunk_file(Path(report_file))
    
    # Get summary statistics
    summary = aggregator.get_summary()
    summary['dataset_id'] = args.dataset_id
    summary['reference_panel'] = args.ref_name
    
    # Write summary to JSON
    json_file = f"{args.output_prefix}_{args.ref_name}.summary.json"
    with open(json_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Write human-readable summary
    txt_file = f"{args.output_prefix}_{args.ref_name}.summary.txt"
    with open(txt_file, 'w') as f:
        f.write("Dataset-Level Summary (Aggregated from Chunks)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Dataset ID: {summary['dataset_id']}\n")
        f.write(f"Reference Panel: {summary['reference_panel']}\n")
        f.write(f"Chunks Processed: {summary['chunks_processed']:,}\n")
        if summary['chunks_failed'] > 0:
            f.write(f"Chunks Failed: {summary['chunks_failed']}\n")
        f.write("\n")
        
        f.write("Variant Statistics:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total Variants: {summary['total_variants']:,}\n")
        f.write(f"Well Imputed Variants: {summary['well_imputed_variants']:,}\n")
        f.write(f"Imputation Rate: {summary['imputation_rate']:.2f}%\n\n")
        
        f.write("R² Statistics (across all chunks):\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mean R²: {summary['mean_r2']:.4f}\n")
        if summary['median_r2'] is not None:
            f.write(f"Median R²: {summary['median_r2']:.4f}\n")
            f.write(f"Std Dev R²: {summary['std_r2']:.4f}\n")
            f.write(f"25th Percentile R²: {summary['percentile_25_r2']:.4f}\n")
            f.write(f"75th Percentile R²: {summary['percentile_75_r2']:.4f}\n")
        else:
            f.write("  (Detailed percentiles not available for >100 chunks)\n")
        f.write("\n")
        
        f.write(f"Mean MAF: {summary['mean_maf']:.4f}\n")
        
        if aggregator.metrics['failed_chunks']:
            f.write("\nWarning: Some chunks failed to process:\n")
            for failed in aggregator.metrics['failed_chunks'][:5]:  # Show first 5
                f.write(f"  - {failed['file']}: {failed['error']}\n")
            if len(aggregator.metrics['failed_chunks']) > 5:
                f.write(f"  ... and {len(aggregator.metrics['failed_chunks']) - 5} more\n")
    
    print(f"Summary written to {json_file} and {txt_file}")
    
    # Exit with warning if some chunks failed
    if aggregator.metrics['failed_chunks']:
        sys.exit(1)


if __name__ == '__main__':
    main()