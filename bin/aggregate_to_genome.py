#!/usr/bin/env python3
"""
Aggregate chromosome-level summaries to genome-wide level for imputation pipeline.
Combines metrics from all chromosomes into a comprehensive genome summary.
"""

import json
import argparse
from pathlib import Path
import numpy as np
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GenomeAggregator:
    """Aggregates chromosome-level summaries to genome-wide level."""
    
    def __init__(self):
        self.metrics = {
            'dataset': None,
            'ref_name': None,
            'chromosomes_processed': 0,
            'total_chunks': 0,
            'total_variants': 0,
            'well_imputed_variants': 0,
            'total_genotyped': 0,
            'total_imputed': 0,
            'concordance_values': [],
            'r2_weighted_sum': 0.0,
            'r2_weight_total': 0,
            'maf_bins': defaultdict(int),
            'r2_by_maf': defaultdict(list),
            'chromosome_details': [],
            'chromosomes': []
        }
        
    def add_chromosome_summary(self, summary_path):
        """Add metrics from a chromosome summary."""
        try:
            with open(summary_path, 'r') as f:
                chr_data = json.load(f)
                
            # Extract chromosome identifier
            chromosome = chr_data.get('chromosome', Path(summary_path).stem.split('_')[0])
            self.metrics['chromosomes'].append(chromosome)
            
            # Update basic counts
            self.metrics['chromosomes_processed'] += 1
            self.metrics['total_chunks'] += chr_data.get('chunks_processed', 0)
            self.metrics['total_variants'] += chr_data.get('total_variants', 0)
            self.metrics['well_imputed_variants'] += chr_data.get('well_imputed_variants', 0)
            self.metrics['total_genotyped'] += chr_data.get('total_genotyped', 0)
            self.metrics['total_imputed'] += chr_data.get('total_imputed', 0)
            
            # Collect concordance values
            if chr_data.get('mean_concordance') is not None:
                self.metrics['concordance_values'].append(chr_data['mean_concordance'])
            
            # Weighted R2 aggregation
            if chr_data.get('mean_r2') is not None:
                weight = chr_data.get('total_variants', 1)
                self.metrics['r2_weighted_sum'] += chr_data['mean_r2'] * weight
                self.metrics['r2_weight_total'] += weight
            
            # Aggregate MAF bins
            if 'maf_bins' in chr_data:
                for bin_name, count in chr_data['maf_bins'].items():
                    self.metrics['maf_bins'][bin_name] += count
            
            # Aggregate R2 by MAF (using mean values from chromosomes)
            if 'mean_r2_by_maf' in chr_data:
                for bin_name, stats in chr_data['mean_r2_by_maf'].items():
                    if isinstance(stats, dict) and 'mean' in stats:
                        # Weight by number of variants
                        weight = stats.get('n', 1)
                        self.metrics['r2_by_maf'][bin_name].append({
                            'value': stats['mean'],
                            'weight': weight
                        })
            
            # Store chromosome summary
            chr_summary = {
                'chromosome': chromosome,
                'chunks': chr_data.get('chunks_processed', 0),
                'variants': chr_data.get('total_variants', 0),
                'well_imputed': chr_data.get('well_imputed_variants', 0),
                'mean_r2': chr_data.get('mean_r2', 0),
                'mean_concordance': chr_data.get('mean_concordance', 0)
            }
            
            # Add info score stats if available
            if 'info_score_stats' in chr_data:
                chr_summary['info_stats'] = chr_data['info_score_stats']
            
            self.metrics['chromosome_details'].append(chr_summary)
            
            logger.info(f"Added chromosome summary: {chromosome}")
            
        except Exception as e:
            logger.error(f"Error processing chromosome summary {summary_path}: {e}")
    
    def calculate_summary(self):
        """Calculate genome-wide summary statistics."""
        summary = {
            'dataset': self.metrics['dataset'],
            'ref_name': self.metrics['ref_name'],
            'chromosomes_processed': self.metrics['chromosomes_processed'],
            'chromosomes': sorted(self.metrics['chromosomes']),
            'total_chunks': self.metrics['total_chunks'],
            'total_variants': self.metrics['total_variants'],
            'well_imputed_variants': self.metrics['well_imputed_variants'],
            'total_genotyped': self.metrics['total_genotyped'],
            'total_imputed': self.metrics['total_imputed']
        }
        
        # Calculate genome-wide averages
        if self.metrics['concordance_values']:
            summary['mean_concordance'] = np.mean(self.metrics['concordance_values'])
            summary['concordance_std'] = np.std(self.metrics['concordance_values'])
        else:
            summary['mean_concordance'] = None
            summary['concordance_std'] = None
            
        if self.metrics['r2_weight_total'] > 0:
            summary['mean_r2'] = self.metrics['r2_weighted_sum'] / self.metrics['r2_weight_total']
        else:
            summary['mean_r2'] = None
        
        # Calculate weighted R2 by MAF
        summary['mean_r2_by_maf'] = {}
        for bin_name, weighted_values in self.metrics['r2_by_maf'].items():
            if weighted_values:
                total_weight = sum(item['weight'] for item in weighted_values)
                if total_weight > 0:
                    weighted_mean = sum(item['value'] * item['weight'] for item in weighted_values) / total_weight
                    summary['mean_r2_by_maf'][bin_name] = {
                        'mean': weighted_mean,
                        'n_variants': total_weight
                    }
        
        # Convert defaultdicts to regular dicts
        summary['maf_bins'] = dict(self.metrics['maf_bins'])
        
        # Add chromosome details sorted by chromosome number
        summary['chromosome_details'] = sorted(
            self.metrics['chromosome_details'],
            key=lambda x: self._chromosome_sort_key(x['chromosome'])
        )
        
        # Calculate genome-wide info score statistics
        all_info_stats = [chr_detail['info_stats'] 
                         for chr_detail in self.metrics['chromosome_details'] 
                         if 'info_stats' in chr_detail]
        
        if all_info_stats:
            summary['genome_info_stats'] = {
                'mean': np.mean([s['mean'] for s in all_info_stats]),
                'median': np.median([s['median'] for s in all_info_stats]),
                'min': min(s['min'] for s in all_info_stats),
                'max': max(s['max'] for s in all_info_stats),
                'total_above_0.3': sum(s['above_0.3'] for s in all_info_stats),
                'total_above_0.8': sum(s['above_0.8'] for s in all_info_stats)
            }
        
        # Calculate imputation rate
        if summary['total_genotyped'] > 0:
            summary['imputation_rate'] = summary['total_imputed'] / summary['total_genotyped']
        else:
            summary['imputation_rate'] = None
        
        # Calculate well-imputed rate
        if summary['total_variants'] > 0:
            summary['well_imputed_rate'] = summary['well_imputed_variants'] / summary['total_variants']
        else:
            summary['well_imputed_rate'] = None
        
        return summary
    
    def _chromosome_sort_key(self, chromosome):
        """Generate sort key for chromosome names."""
        # Handle numeric chromosomes and X, Y, MT
        chr_str = str(chromosome).replace('chr', '').upper()
        if chr_str.isdigit():
            return (0, int(chr_str))
        elif chr_str == 'X':
            return (1, 23)
        elif chr_str == 'Y':
            return (1, 24)
        elif chr_str in ['MT', 'M']:
            return (1, 25)
        else:
            return (2, chr_str)
    
    def generate_stats_text(self, summary):
        """Generate human-readable statistics text."""
        lines = [
            f"Genome-Wide Imputation Statistics",
            f"=" * 50,
            f"Dataset: {summary.get('dataset', 'Unknown')}",
            f"Reference: {summary.get('ref_name', 'Unknown')}",
            f"",
            f"Summary:",
            f"  Chromosomes processed: {summary.get('chromosomes_processed', 0)}",
            f"  Total chunks: {summary.get('total_chunks', 0)}",
            f"  Total variants: {summary.get('total_variants', 0):,}",
            f"  Well-imputed variants: {summary.get('well_imputed_variants', 0):,}",
            f"  Total genotyped: {summary.get('total_genotyped', 0):,}",
            f"  Total imputed: {summary.get('total_imputed', 0):,}",
            f""
        ]
        
        if summary.get('mean_r2') is not None:
            lines.append(f"  Mean R² (weighted): {summary['mean_r2']:.4f}")
        
        if summary.get('mean_concordance') is not None:
            lines.append(f"  Mean concordance: {summary['mean_concordance']:.4f} "
                        f"(± {summary.get('concordance_std', 0):.4f})")
        
        if summary.get('imputation_rate') is not None:
            lines.append(f"  Imputation rate: {summary['imputation_rate']:.2%}")
        
        if summary.get('well_imputed_rate') is not None:
            lines.append(f"  Well-imputed rate: {summary['well_imputed_rate']:.2%}")
        
        lines.append("")
        
        # MAF distribution
        if summary.get('maf_bins'):
            lines.append("MAF Distribution (genome-wide):")
            total_maf_variants = sum(summary['maf_bins'].values())
            for bin_name in sorted(summary['maf_bins'].keys()):
                count = summary['maf_bins'][bin_name]
                percentage = (count / total_maf_variants * 100) if total_maf_variants > 0 else 0
                lines.append(f"  {bin_name}: {count:,} ({percentage:.1f}%)")
            lines.append("")
        
        # R2 by MAF
        if summary.get('mean_r2_by_maf'):
            lines.append("Mean R² by MAF (genome-wide):")
            for bin_name in sorted(summary['mean_r2_by_maf'].keys()):
                stats = summary['mean_r2_by_maf'][bin_name]
                lines.append(f"  {bin_name}: {stats['mean']:.4f}")
            lines.append("")
        
        # Genome-wide info score statistics
        if summary.get('genome_info_stats'):
            stats = summary['genome_info_stats']
            lines.extend([
                "Genome-wide Info Score Statistics:",
                f"  Mean: {stats['mean']:.4f}",
                f"  Median: {stats['median']:.4f}",
                f"  Min: {stats['min']:.4f}",
                f"  Max: {stats['max']:.4f}",
                f"  Total above 0.3: {stats['total_above_0.3']:,} variants",
                f"  Total above 0.8: {stats['total_above_0.8']:,} variants",
                ""
            ])
        
        # Per-chromosome summary
        lines.append("Per-Chromosome Summary:")
        lines.append(f"{'Chr':<5} {'Chunks':<8} {'Variants':<12} {'Well-Imp':<12} {'Mean R²':<10} {'Concordance':<12}")
        lines.append("-" * 70)
        
        for chr_detail in summary.get('chromosome_details', []):
            chr_name = str(chr_detail['chromosome'])[:4]
            chunks = chr_detail.get('chunks', 0)
            variants = chr_detail.get('variants', 0)
            well_imp = chr_detail.get('well_imputed', 0)
            mean_r2 = chr_detail.get('mean_r2', 0) or 0
            concordance = chr_detail.get('mean_concordance', 0) or 0
            
            lines.append(f"{chr_name:<5} {chunks:<8} {variants:<12,} {well_imp:<12,} "
                        f"{mean_r2:<10.4f} {concordance:<12.4f}")
        
        return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Aggregate chromosome summaries to genome-wide level')
    parser.add_argument('--chr-summaries', nargs='+', required=True,
                       help='Chromosome summary JSON files')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    
    args = parser.parse_args()
    
    # Initialize aggregator
    aggregator = GenomeAggregator()
    aggregator.metrics['dataset'] = args.dataset
    aggregator.metrics['ref_name'] = args.ref_name
    
    # Process chromosome summaries
    logger.info(f"Processing {len(args.chr_summaries)} chromosome summaries...")
    for summary_path in args.chr_summaries:
        aggregator.add_chromosome_summary(summary_path)
    
    # Calculate genome-wide summary
    logger.info("Calculating genome-wide summary...")
    summary = aggregator.calculate_summary()
    
    # Save JSON summary
    json_output = f"{args.output_prefix}_{args.ref_name}.genome_summary.json"
    with open(json_output, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved JSON summary to {json_output}")
    
    # Save text statistics
    text_output = f"{args.output_prefix}_{args.ref_name}.genome_stats.txt"
    stats_text = aggregator.generate_stats_text(summary)
    with open(text_output, 'w') as f:
        f.write(stats_text)
    logger.info(f"Saved text statistics to {text_output}")
    
    # Print summary to stdout
    print(f"\nGenome-Wide Summary for {args.dataset}:")
    print(f"  Chromosomes: {summary['chromosomes_processed']}")
    print(f"  Total chunks: {summary['total_chunks']}")
    print(f"  Total variants: {summary['total_variants']:,}")
    if summary.get('mean_r2') is not None:
        print(f"  Mean R² (weighted): {summary['mean_r2']:.4f}")
    if summary.get('well_imputed_rate') is not None:
        print(f"  Well-imputed rate: {summary['well_imputed_rate']:.2%}")
    
    logger.info("Genome-wide aggregation completed successfully")


if __name__ == '__main__':
    main()