#!/usr/bin/env python3
"""
Aggregate chunk-level metrics to chromosome level for imputation pipeline.
Combines metrics from multiple chunks of the same chromosome.
"""

import json
import argparse
import sys
from pathlib import Path
import numpy as np
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ChromosomeAggregator:
    """Aggregates chunk-level metrics to chromosome level."""
    
    def __init__(self):
        self.metrics = {
            'chromosome': None,
            'dataset': None,
            'ref_name': None,
            'chunks_processed': 0,
            'total_variants': 0,
            'well_imputed_variants': 0,
            'total_genotyped': 0,
            'total_imputed': 0,
            'concordance_sum': 0.0,
            'concordance_count': 0,
            'r2_sum': 0.0,
            'r2_count': 0,
            'maf_bins': defaultdict(int),
            'r2_by_maf': defaultdict(list),
            'info_score_dist': [],
            'chunk_details': []
        }
        
    def add_chunk_report(self, report_path):
        """Add metrics from a chunk report."""
        try:
            with open(report_path, 'r') as f:
                chunk_data = json.load(f)
                
            # Update basic counts
            self.metrics['chunks_processed'] += 1
            self.metrics['total_variants'] += chunk_data.get('total_variants', 0)
            self.metrics['well_imputed_variants'] += chunk_data.get('well_imputed_variants', 0)
            self.metrics['total_genotyped'] += chunk_data.get('total_genotyped', 0)
            self.metrics['total_imputed'] += chunk_data.get('total_imputed', 0)
            
            # Update concordance
            if 'concordance' in chunk_data and chunk_data['concordance'] is not None:
                self.metrics['concordance_sum'] += chunk_data['concordance']
                self.metrics['concordance_count'] += 1
            
            # Update R2
            if 'mean_r2' in chunk_data and chunk_data['mean_r2'] is not None:
                self.metrics['r2_sum'] += chunk_data['mean_r2'] * chunk_data.get('total_variants', 1)
                self.metrics['r2_count'] += chunk_data.get('total_variants', 1)
            
            # Aggregate MAF bins
            if 'maf_distribution' in chunk_data:
                for bin_name, count in chunk_data['maf_distribution'].items():
                    self.metrics['maf_bins'][bin_name] += count
            
            # Aggregate R2 by MAF
            if 'r2_by_maf' in chunk_data:
                for bin_name, r2_values in chunk_data['r2_by_maf'].items():
                    if isinstance(r2_values, list):
                        self.metrics['r2_by_maf'][bin_name].extend(r2_values)
                    elif isinstance(r2_values, (int, float)):
                        self.metrics['r2_by_maf'][bin_name].append(r2_values)
            
            # Store chunk summary
            chunk_summary = {
                'chunk_id': chunk_data.get('chunk_id', Path(report_path).stem),
                'variants': chunk_data.get('total_variants', 0),
                'well_imputed': chunk_data.get('well_imputed_variants', 0),
                'mean_r2': chunk_data.get('mean_r2', 0)
            }
            self.metrics['chunk_details'].append(chunk_summary)
            
            logger.info(f"Added chunk report: {report_path}")
            
        except Exception as e:
            logger.error(f"Error processing chunk report {report_path}: {e}")
    
    def add_info_file(self, info_path):
        """Add imputation info scores from minimac4 info file."""
        try:
            with open(info_path, 'r') as f:
                # Skip header
                header = f.readline().strip().split('\t')
                r2_idx = header.index('Rsq') if 'Rsq' in header else None
                
                if r2_idx is None:
                    logger.warning(f"No Rsq column found in {info_path}")
                    return
                
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) > r2_idx:
                        try:
                            r2_score = float(fields[r2_idx])
                            self.metrics['info_score_dist'].append(r2_score)
                        except ValueError:
                            continue
                            
            logger.info(f"Added info file: {info_path}")
            
        except Exception as e:
            logger.error(f"Error processing info file {info_path}: {e}")
    
    def calculate_summary(self):
        """Calculate summary statistics."""
        summary = dict(self.metrics)
        
        # Calculate averages
        if self.metrics['concordance_count'] > 0:
            summary['mean_concordance'] = self.metrics['concordance_sum'] / self.metrics['concordance_count']
        else:
            summary['mean_concordance'] = None
            
        if self.metrics['r2_count'] > 0:
            summary['mean_r2'] = self.metrics['r2_sum'] / self.metrics['r2_count']
        else:
            summary['mean_r2'] = None
        
        # Calculate R2 statistics by MAF
        summary['mean_r2_by_maf'] = {}
        for bin_name, r2_values in self.metrics['r2_by_maf'].items():
            if r2_values:
                summary['mean_r2_by_maf'][bin_name] = {
                    'mean': np.mean(r2_values),
                    'median': np.median(r2_values),
                    'std': np.std(r2_values),
                    'n': len(r2_values)
                }
        
        # Calculate info score statistics
        if self.metrics['info_score_dist']:
            summary['info_score_stats'] = {
                'mean': np.mean(self.metrics['info_score_dist']),
                'median': np.median(self.metrics['info_score_dist']),
                'min': np.min(self.metrics['info_score_dist']),
                'max': np.max(self.metrics['info_score_dist']),
                'q25': np.percentile(self.metrics['info_score_dist'], 25),
                'q75': np.percentile(self.metrics['info_score_dist'], 75),
                'above_0.3': sum(1 for x in self.metrics['info_score_dist'] if x > 0.3),
                'above_0.8': sum(1 for x in self.metrics['info_score_dist'] if x > 0.8)
            }
        
        # Convert defaultdicts to regular dicts for JSON serialization
        summary['maf_bins'] = dict(summary['maf_bins'])
        summary['r2_by_maf'] = dict(summary['r2_by_maf'])
        
        # Remove raw distributions to keep file size manageable
        del summary['info_score_dist']
        del summary['concordance_sum']
        del summary['concordance_count']
        del summary['r2_sum']
        del summary['r2_count']
        
        return summary
    
    def generate_stats_text(self, summary):
        """Generate human-readable statistics text."""
        lines = [
            f"Chromosome-Level Imputation Statistics",
            f"=" * 50,
            f"Dataset: {summary.get('dataset', 'Unknown')}",
            f"Reference: {summary.get('ref_name', 'Unknown')}",
            f"Chromosome: {summary.get('chromosome', 'Unknown')}",
            f"",
            f"Summary:",
            f"  Chunks processed: {summary.get('chunks_processed', 0)}",
            f"  Total variants: {summary.get('total_variants', 0):,}",
            f"  Well-imputed variants: {summary.get('well_imputed_variants', 0):,}",
            f"  Total genotyped: {summary.get('total_genotyped', 0):,}",
            f"  Total imputed: {summary.get('total_imputed', 0):,}",
            f""
        ]
        
        if summary.get('mean_r2') is not None:
            lines.append(f"  Mean R²: {summary['mean_r2']:.4f}")
        
        if summary.get('mean_concordance') is not None:
            lines.append(f"  Mean concordance: {summary['mean_concordance']:.4f}")
        
        lines.append("")
        
        # MAF distribution
        if summary.get('maf_bins'):
            lines.append("MAF Distribution:")
            for bin_name in sorted(summary['maf_bins'].keys()):
                count = summary['maf_bins'][bin_name]
                lines.append(f"  {bin_name}: {count:,}")
            lines.append("")
        
        # R2 by MAF
        if summary.get('mean_r2_by_maf'):
            lines.append("Mean R² by MAF:")
            for bin_name in sorted(summary['mean_r2_by_maf'].keys()):
                stats = summary['mean_r2_by_maf'][bin_name]
                lines.append(f"  {bin_name}: {stats['mean']:.4f} (n={stats['n']:,})")
            lines.append("")
        
        # Info score statistics
        if summary.get('info_score_stats'):
            stats = summary['info_score_stats']
            lines.extend([
                "Info Score Statistics:",
                f"  Mean: {stats['mean']:.4f}",
                f"  Median: {stats['median']:.4f}",
                f"  Min: {stats['min']:.4f}",
                f"  Max: {stats['max']:.4f}",
                f"  Q25: {stats['q25']:.4f}",
                f"  Q75: {stats['q75']:.4f}",
                f"  Above 0.3: {stats['above_0.3']:,} variants",
                f"  Above 0.8: {stats['above_0.8']:,} variants",
                ""
            ])
        
        # Top performing chunks
        if summary.get('chunk_details'):
            chunks_sorted = sorted(summary['chunk_details'], 
                                 key=lambda x: x.get('mean_r2', 0), 
                                 reverse=True)
            lines.append("Top 5 Chunks by R²:")
            for chunk in chunks_sorted[:5]:
                lines.append(f"  {chunk['chunk_id']}: R²={chunk.get('mean_r2', 0):.4f}, "
                           f"Variants={chunk['variants']:,}")
        
        return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Aggregate chunk-level metrics to chromosome level')
    parser.add_argument('--chunk-reports', nargs='+', required=True,
                       help='Chunk report JSON files')
    parser.add_argument('--chunk-info', nargs='+', 
                       help='Chunk info files from minimac4')
    parser.add_argument('--output-prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    parser.add_argument('--chromosome', required=True,
                       help='Chromosome identifier')
    
    args = parser.parse_args()
    
    # Initialize aggregator
    aggregator = ChromosomeAggregator()
    aggregator.metrics['chromosome'] = args.chromosome
    aggregator.metrics['dataset'] = args.dataset
    aggregator.metrics['ref_name'] = args.ref_name
    
    # Process chunk reports
    logger.info(f"Processing {len(args.chunk_reports)} chunk reports...")
    for report_path in args.chunk_reports:
        aggregator.add_chunk_report(report_path)
    
    # Process info files if provided
    if args.chunk_info:
        logger.info(f"Processing {len(args.chunk_info)} info files...")
        for info_path in args.chunk_info:
            aggregator.add_info_file(info_path)
    
    # Calculate summary
    logger.info("Calculating chromosome-level summary...")
    summary = aggregator.calculate_summary()
    
    # Save JSON summary
    json_output = f"{args.output_prefix}_{args.ref_name}.chr_summary.json"
    with open(json_output, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Saved JSON summary to {json_output}")
    
    # Save text statistics
    text_output = f"{args.output_prefix}_{args.ref_name}.chr_stats.txt"
    stats_text = aggregator.generate_stats_text(summary)
    with open(text_output, 'w') as f:
        f.write(stats_text)
    logger.info(f"Saved text statistics to {text_output}")
    
    # Print summary to stdout
    print(f"\nChromosome {args.chromosome} Summary:")
    print(f"  Chunks: {summary['chunks_processed']}")
    print(f"  Total variants: {summary['total_variants']:,}")
    if summary.get('mean_r2') is not None:
        print(f"  Mean R²: {summary['mean_r2']:.4f}")
    
    logger.info("Chromosome aggregation completed successfully")


if __name__ == '__main__':
    main()