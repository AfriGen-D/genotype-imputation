#!/usr/bin/env python3
"""
Plot dosage distribution to assess imputation certainty
Dosages closer to 0, 1, or 2 indicate higher confidence in genotype calls
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
from pathlib import Path

def read_vcf_dosages(vcf_file):
    """Extract dosage values from VCF file (DS field)"""
    dosages = []
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    sample_idx = 9  # First sample column
                continue
            
            fields = line.strip().split('\t')
            format_field = fields[8]
            
            # Check if DS (dosage) field exists
            if 'DS' not in format_field:
                continue
                
            format_keys = format_field.split(':')
            ds_idx = format_keys.index('DS')
            
            # Process each sample
            for sample_field in fields[9:]:
                values = sample_field.split(':')
                if len(values) > ds_idx and values[ds_idx] != '.':
                    try:
                        dosage = float(values[ds_idx])
                        dosages.append(dosage)
                    except ValueError:
                        continue
    
    return np.array(dosages)

def categorize_dosages(dosages, threshold=0.1):
    """Categorize dosages by certainty level"""
    categories = {
        'High certainty (near 0)': np.sum(dosages <= threshold),
        'High certainty (near 1)': np.sum((dosages >= 1-threshold) & (dosages <= 1+threshold)),
        'High certainty (near 2)': np.sum(dosages >= 2-threshold),
        'Low certainty (0.2-0.8)': np.sum((dosages > 0.2) & (dosages < 0.8)),
        'Low certainty (1.2-1.8)': np.sum((dosages > 1.2) & (dosages < 1.8))
    }
    return categories

def main():
    parser = argparse.ArgumentParser(description='Plot dosage distribution for imputed genotypes')
    parser.add_argument('vcf_file', help='Input VCF file with dosage information')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chunk-id', help='Chunk identifier', default='')
    
    args = parser.parse_args()
    
    print(f"Reading dosages from {args.vcf_file}...")
    dosages = read_vcf_dosages(args.vcf_file)
    
    if len(dosages) == 0:
        print("WARNING: No dosage values found in VCF file")
        # Create empty plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No dosage data available', 
                ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    else:
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Overall dosage distribution histogram
        ax = axes[0, 0]
        ax.hist(dosages, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.5, label='Expected values')
        ax.axvline(x=1, color='red', linestyle='--', alpha=0.5)
        ax.axvline(x=2, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Dosage Value')
        ax.set_ylabel('Frequency')
        ax.set_title('Overall Dosage Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Dosage certainty distribution
        ax = axes[0, 1]
        certainty = np.minimum(
            np.minimum(np.abs(dosages), np.abs(dosages - 1)),
            np.abs(dosages - 2)
        )
        ax.hist(certainty, bins=30, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel('Distance to Nearest Integer (0, 1, or 2)')
        ax.set_ylabel('Frequency')
        ax.set_title('Genotype Certainty Distribution')
        ax.set_xlim(0, 0.5)
        ax.grid(True, alpha=0.3)
        
        # 3. Cumulative distribution
        ax = axes[1, 0]
        sorted_dosages = np.sort(dosages)
        cumulative = np.arange(1, len(sorted_dosages) + 1) / len(sorted_dosages)
        ax.plot(sorted_dosages, cumulative, linewidth=2, color='green')
        ax.set_xlabel('Dosage Value')
        ax.set_ylabel('Cumulative Proportion')
        ax.set_title('Cumulative Dosage Distribution')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 1)
        
        # 4. Summary statistics
        ax = axes[1, 1]
        ax.axis('off')
        
        # Calculate statistics
        categories = categorize_dosages(dosages)
        stats_text = f"Summary Statistics\n" + "="*30 + "\n\n"
        stats_text += f"Total genotypes: {len(dosages):,}\n"
        stats_text += f"Mean dosage: {np.mean(dosages):.3f}\n"
        stats_text += f"Median dosage: {np.median(dosages):.3f}\n"
        stats_text += f"Std deviation: {np.std(dosages):.3f}\n\n"
        
        stats_text += "Certainty Categories:\n"
        total = len(dosages)
        for category, count in categories.items():
            pct = (count / total) * 100 if total > 0 else 0
            stats_text += f"  {category}: {count:,} ({pct:.1f}%)\n"
        
        # High certainty percentage (within 0.1 of 0, 1, or 2)
        high_cert = (categories['High certainty (near 0)'] + 
                    categories['High certainty (near 1)'] + 
                    categories['High certainty (near 2)'])
        high_cert_pct = (high_cert / total * 100) if total > 0 else 0
        stats_text += f"\nOverall high certainty: {high_cert_pct:.1f}%"
        
        ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    # Add main title
    chunk_info = f" - {args.chunk_id}" if args.chunk_id else ""
    fig.suptitle(f'Dosage Distribution Analysis\nSample: {args.sample_id} | Reference: {args.ref_name}{chunk_info}',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(args.output_pdf, dpi=150, bbox_inches='tight')
    print(f"Dosage distribution plot saved to {args.output_pdf}")

if __name__ == "__main__":
    main()