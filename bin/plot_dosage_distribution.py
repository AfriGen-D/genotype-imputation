#!/usr/bin/env python3
"""
Optimized plot dosage distribution with streaming processing for large VCF files
Dosages closer to 0, 1, or 2 indicate higher confidence in genotype calls
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import gzip
from pathlib import Path
from collections import defaultdict

class DosageStats:
    """Incremental statistics calculator for dosage values"""
    
    def __init__(self, n_bins=50, max_dosage=2.5):
        self.n_bins = n_bins
        self.max_dosage = max_dosage
        self.bin_edges = np.linspace(0, max_dosage, n_bins + 1)
        self.bin_counts = np.zeros(n_bins)
        self.certainty_bins = np.zeros(30)  # For certainty distribution
        self.certainty_edges = np.linspace(0, 0.5, 31)
        
        # For incremental statistics
        self.n = 0
        self.sum = 0.0
        self.sum_sq = 0.0
        self.min_val = float('inf')
        self.max_val = float('-inf')
        
        # For categories
        self.categories = defaultdict(int)
        self.threshold = 0.1
        
        # For sampling (keep a subset for cumulative plot)
        self.sample_size = 10000
        self.sample_step = 1
        self.samples = []
        
    def add_dosage(self, dosage):
        """Add a single dosage value to the statistics"""
        self.n += 1
        self.sum += dosage
        self.sum_sq += dosage * dosage
        self.min_val = min(self.min_val, dosage)
        self.max_val = max(self.max_val, dosage)
        
        # Update histogram bins
        bin_idx = np.searchsorted(self.bin_edges, dosage, side='right') - 1
        if 0 <= bin_idx < self.n_bins:
            self.bin_counts[bin_idx] += 1
        
        # Update certainty histogram
        certainty = min(abs(dosage), abs(dosage - 1), abs(dosage - 2))
        cert_idx = np.searchsorted(self.certainty_edges, certainty, side='right') - 1
        if 0 <= cert_idx < len(self.certainty_bins):
            self.certainty_bins[cert_idx] += 1
        
        # Update categories
        if dosage <= self.threshold:
            self.categories['High certainty (near 0)'] += 1
        elif (1 - self.threshold) <= dosage <= (1 + self.threshold):
            self.categories['High certainty (near 1)'] += 1
        elif dosage >= (2 - self.threshold):
            self.categories['High certainty (near 2)'] += 1
        elif 0.2 < dosage < 0.8:
            self.categories['Low certainty (0.2-0.8)'] += 1
        elif 1.2 < dosage < 1.8:
            self.categories['Low certainty (1.2-1.8)'] += 1
        
        # Sample for cumulative plot (reservoir sampling after initial samples)
        if self.n <= self.sample_size:
            self.samples.append(dosage)
        elif self.n % self.sample_step == 0:
            # Randomly replace an element with decreasing probability
            idx = np.random.randint(0, len(self.samples))
            if np.random.random() < self.sample_size / self.n:
                self.samples[idx] = dosage
    
    @property
    def mean(self):
        return self.sum / self.n if self.n > 0 else 0
    
    @property
    def std(self):
        if self.n > 1:
            variance = (self.sum_sq - (self.sum * self.sum) / self.n) / (self.n - 1)
            return np.sqrt(max(0, variance))  # Ensure non-negative
        return 0
    
    @property
    def median_approx(self):
        """Approximate median from histogram"""
        cumsum = np.cumsum(self.bin_counts)
        median_count = self.n / 2
        median_bin = np.searchsorted(cumsum, median_count)
        if median_bin < len(self.bin_edges) - 1:
            return (self.bin_edges[median_bin] + self.bin_edges[median_bin + 1]) / 2
        return self.bin_edges[-1]

def process_vcf_streaming(vcf_file, max_variants=None):
    """Process VCF file in streaming fashion, collecting statistics"""
    stats = DosageStats()
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    variants_processed = 0
    
    print(f"Processing {vcf_file} in streaming mode...")
    
    with opener(vcf_file, 'rt') as f:
        for line_num, line in enumerate(f):
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    n_samples = len(header) - 9
                    print(f"Found {n_samples} samples in VCF")
                continue
            
            # Process variant line
            fields = line.strip().split('\t')
            format_field = fields[8]
            
            # Check if DS (dosage) field exists
            if 'DS' not in format_field:
                continue
                
            format_keys = format_field.split(':')
            try:
                ds_idx = format_keys.index('DS')
            except ValueError:
                continue
            
            # Process each sample's dosage
            for sample_field in fields[9:]:
                values = sample_field.split(':')
                if len(values) > ds_idx and values[ds_idx] != '.':
                    try:
                        dosage = float(values[ds_idx])
                        stats.add_dosage(dosage)
                    except (ValueError, IndexError):
                        continue
            
            variants_processed += 1
            
            # Progress reporting
            if variants_processed % 10000 == 0:
                print(f"  Processed {variants_processed} variants, {stats.n:,} dosages collected...")
            
            # Early stopping for testing
            if max_variants and variants_processed >= max_variants:
                print(f"  Stopped at {max_variants} variants (for testing)")
                break
    
    print(f"Finished processing: {variants_processed} variants, {stats.n:,} total dosages")
    return stats

def create_plots(stats, sample_id, ref_name, chunk_id, output_pdf):
    """Create visualization plots from collected statistics"""
    
    if stats.n == 0:
        print("WARNING: No dosage values found in VCF file")
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
        # Use the pre-computed histogram bins
        bin_centers = (stats.bin_edges[:-1] + stats.bin_edges[1:]) / 2
        ax.bar(bin_centers, stats.bin_counts, width=stats.bin_edges[1] - stats.bin_edges[0],
               edgecolor='black', alpha=0.7, color='steelblue')
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
        cert_centers = (stats.certainty_edges[:-1] + stats.certainty_edges[1:]) / 2
        ax.bar(cert_centers, stats.certainty_bins, 
               width=stats.certainty_edges[1] - stats.certainty_edges[0],
               edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel('Distance to Nearest Integer (0, 1, or 2)')
        ax.set_ylabel('Frequency')
        ax.set_title('Genotype Certainty Distribution')
        ax.set_xlim(0, 0.5)
        ax.grid(True, alpha=0.3)
        
        # 3. Cumulative distribution (from sample)
        ax = axes[1, 0]
        if len(stats.samples) > 0:
            sorted_samples = np.sort(stats.samples)
            cumulative = np.arange(1, len(sorted_samples) + 1) / len(sorted_samples)
            ax.plot(sorted_samples, cumulative, linewidth=2, color='green')
            ax.set_xlabel('Dosage Value')
            ax.set_ylabel('Cumulative Proportion')
            ax.set_title(f'Cumulative Dosage Distribution (sample of {len(stats.samples):,})')
        else:
            ax.text(0.5, 0.5, 'Sample data not available', ha='center', va='center')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 2.2)
        ax.set_ylim(0, 1)
        
        # 4. Summary statistics
        ax = axes[1, 1]
        ax.axis('off')
        
        # Calculate statistics
        stats_text = f"Summary Statistics\n" + "="*30 + "\n\n"
        stats_text += f"Total genotypes: {stats.n:,}\n"
        stats_text += f"Mean dosage: {stats.mean:.3f}\n"
        stats_text += f"Median dosage (approx): {stats.median_approx:.3f}\n"
        stats_text += f"Std deviation: {stats.std:.3f}\n"
        stats_text += f"Min: {stats.min_val:.3f}, Max: {stats.max_val:.3f}\n\n"
        
        stats_text += "Certainty Categories:\n"
        total = stats.n
        for category in ['High certainty (near 0)', 'High certainty (near 1)', 
                        'High certainty (near 2)', 'Low certainty (0.2-0.8)', 
                        'Low certainty (1.2-1.8)']:
            count = stats.categories.get(category, 0)
            pct = (count / total) * 100 if total > 0 else 0
            stats_text += f"  {category}: {count:,} ({pct:.1f}%)\n"
        
        # High certainty percentage
        high_cert = (stats.categories.get('High certainty (near 0)', 0) + 
                    stats.categories.get('High certainty (near 1)', 0) + 
                    stats.categories.get('High certainty (near 2)', 0))
        high_cert_pct = (high_cert / total * 100) if total > 0 else 0
        stats_text += f"\nOverall high certainty: {high_cert_pct:.1f}%"
        
        ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    # Add main title
    chunk_info = f" - {chunk_id}" if chunk_id else ""
    fig.suptitle(f'Dosage Distribution Analysis\nSample: {sample_id} | Reference: {ref_name}{chunk_info}',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=150, bbox_inches='tight')
    print(f"Dosage distribution plot saved to {output_pdf}")

def main():
    parser = argparse.ArgumentParser(description='Plot dosage distribution for imputed genotypes (optimized)')
    parser.add_argument('vcf_file', help='Input VCF file with dosage information')
    parser.add_argument('output_pdf', help='Output PDF file')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    parser.add_argument('--chunk-id', help='Chunk identifier', default='')
    parser.add_argument('--max-variants', type=int, help='Process only first N variants (for testing)')
    
    args = parser.parse_args()
    
    # Process VCF file and collect statistics
    stats = process_vcf_streaming(args.vcf_file, args.max_variants)
    
    # Create plots from statistics
    create_plots(stats, args.sample_id, args.ref_name, args.chunk_id, args.output_pdf)

if __name__ == "__main__":
    main()