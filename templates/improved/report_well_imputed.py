#!/usr/bin/env python3
"""
Improved version of report_well_imputed.py
Generates reports on well-imputed variants from imputation results.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import logging
from collections import defaultdict
from dataclasses import dataclass, field
import json

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class VariantStats:
    """Statistics for a variant."""
    rsq: float
    maf: float
    info_score: float = 0.0
    empirical_rsq: Optional[float] = None
    genotyped: bool = False
    
    @property
    def is_well_imputed(self) -> bool:
        """Check if variant meets well-imputed criteria."""
        return self.rsq >= 0.8  # Default threshold


@dataclass
class DatasetStats:
    """Statistics for a dataset."""
    total_variants: int = 0
    well_imputed: int = 0
    genotyped: int = 0
    maf_bins: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    rsq_bins: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    
    @property
    def well_imputed_rate(self) -> float:
        """Calculate well-imputed rate."""
        if self.total_variants == 0:
            return 0.0
        return self.well_imputed / self.total_variants
    
    @property
    def genotyped_rate(self) -> float:
        """Calculate genotyped rate."""
        if self.total_variants == 0:
            return 0.0
        return self.genotyped / self.total_variants


class WellImputedReporter:
    """Generate reports on well-imputed variants."""
    
    def __init__(self, rsq_threshold: float = 0.8, maf_threshold: float = 0.01):
        self.rsq_threshold = rsq_threshold
        self.maf_threshold = maf_threshold
        self.dataset_stats: Dict[str, DatasetStats] = defaultdict(DatasetStats)
        self.variants_by_dataset: Dict[str, Dict[str, VariantStats]] = defaultdict(dict)
        
        # Define MAF bins
        self.maf_bins = [
            (0.0, 0.001, "0-0.1%"),
            (0.001, 0.005, "0.1-0.5%"),
            (0.005, 0.01, "0.5-1%"),
            (0.01, 0.05, "1-5%"),
            (0.05, 0.1, "5-10%"),
            (0.1, 0.5, "10-50%")
        ]
        
        # Define R² bins
        self.rsq_bins = [
            (0.0, 0.3, "0-0.3"),
            (0.3, 0.5, "0.3-0.5"),
            (0.5, 0.7, "0.5-0.7"),
            (0.7, 0.8, "0.7-0.8"),
            (0.8, 0.9, "0.8-0.9"),
            (0.9, 1.0, "0.9-1.0")
        ]
    
    def process_info_file(self, info_file: Path, dataset_name: str) -> None:
        """Process a single info file."""
        if not info_file.exists():
            raise FileNotFoundError(f"Info file not found: {info_file}")
        
        logger.info(f"Processing {dataset_name}: {info_file}")
        
        header_indices = {}
        stats = self.dataset_stats[dataset_name]
        
        with open(info_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                
                # Process header
                if "SNP" in line and "Rsq" in line:
                    header_indices = self._parse_header(parts)
                    continue
                
                if not header_indices:
                    logger.warning(f"No header found in {info_file}")
                    break
                
                # Process variant
                variant_stats = self._process_variant(parts, header_indices)
                if variant_stats:
                    variant_id = parts[header_indices.get('snp', 1)]
                    self.variants_by_dataset[dataset_name][variant_id] = variant_stats
                    
                    # Update statistics
                    stats.total_variants += 1
                    
                    if variant_stats.rsq >= self.rsq_threshold:
                        stats.well_imputed += 1
                    
                    if variant_stats.genotyped:
                        stats.genotyped += 1
                    
                    # Bin by MAF
                    for min_maf, max_maf, bin_name in self.maf_bins:
                        if min_maf <= variant_stats.maf < max_maf:
                            stats.maf_bins[bin_name] += 1
                            break
                    
                    # Bin by R²
                    for min_rsq, max_rsq, bin_name in self.rsq_bins:
                        if min_rsq <= variant_stats.rsq <= max_rsq:
                            stats.rsq_bins[bin_name] += 1
                            break
        
        logger.info(f"  Processed {stats.total_variants:,} variants")
        logger.info(f"  Well-imputed: {stats.well_imputed:,} ({stats.well_imputed_rate:.1%})")
    
    def _parse_header(self, header: List[str]) -> Dict[str, int]:
        """Parse header and return column indices."""
        indices = {}
        
        # Map common column names (case-insensitive)
        column_map = {
            'snp': ['snp', 'id', 'variant', 'rsid'],
            'rsq': ['rsq', 'r2', 'info'],
            'maf': ['maf', 'af', 'freq'],
            'emprsq': ['emprsq', 'empirical_rsq', 'emp_r2'],
            'genotyped': ['genotyped', 'typed']
        }
        
        header_lower = [h.lower() for h in header]
        
        for key, possible_names in column_map.items():
            for i, col in enumerate(header_lower):
                if col in possible_names:
                    indices[key] = i
                    break
        
        # Verify required columns
        if 'rsq' not in indices:
            raise ValueError("Required column 'Rsq' not found in header")
        
        return indices
    
    def _process_variant(
        self, 
        parts: List[str], 
        indices: Dict[str, int]
    ) -> Optional[VariantStats]:
        """Process a single variant line."""
        try:
            # Get R²
            rsq_str = parts[indices['rsq']]
            if rsq_str in ['-', 'NA', '']:
                return None
            rsq = float(rsq_str)
            
            # Get MAF
            maf = 0.0
            if 'maf' in indices:
                maf_str = parts[indices['maf']]
                if maf_str not in ['-', 'NA', '']:
                    maf = float(maf_str)
            
            # Skip low MAF variants if threshold is set
            if maf < self.maf_threshold:
                return None
            
            # Get empirical R² if available
            emp_rsq = None
            if 'emprsq' in indices:
                emp_str = parts[indices['emprsq']]
                if emp_str not in ['-', 'NA', '']:
                    emp_rsq = float(emp_str)
            
            # Check if genotyped
            genotyped = False
            if 'genotyped' in indices:
                genotyped_str = parts[indices['genotyped']]
                genotyped = genotyped_str.lower() in ['true', '1', 'yes']
            
            return VariantStats(
                rsq=rsq,
                maf=maf,
                info_score=rsq,
                empirical_rsq=emp_rsq,
                genotyped=genotyped
            )
            
        except (ValueError, IndexError) as e:
            logger.debug(f"Error processing variant: {e}")
            return None
    
    def generate_summary_report(self, output_file: Path) -> None:
        """Generate summary report across all datasets."""
        logger.info(f"Generating summary report: {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("Well-Imputed Variants Summary Report\\n")
            f.write("=" * 70 + "\\n\\n")
            f.write(f"R² Threshold: {self.rsq_threshold}\n")
            f.write(f"MAF Threshold: {self.maf_threshold}\n\n")
            
            # Overall summary
            total_all = sum(s.total_variants for s in self.dataset_stats.values())
            well_imputed_all = sum(s.well_imputed for s in self.dataset_stats.values())
            
            f.write("Overall Summary\\n")
            f.write("-" * 40 + "\\n")
            f.write(f"Total datasets: {len(self.dataset_stats)}\n")
            f.write(f"Total variants: {total_all:,}\n")
            f.write(f"Well-imputed variants: {well_imputed_all:,}\n")
            if total_all > 0:
                f.write(f"Overall well-imputed rate: {100 * well_imputed_all / total_all:.2f}%\n")
            f.write("\\n")
            
            # Per-dataset summary
            f.write("Per-Dataset Summary\\n")
            f.write("-" * 40 + "\\n")
            f.write(f"{'Dataset':<30} {'Total':<12} {'Well-Imputed':<15} {'Rate':<8}\n")
            f.write("-" * 70 + "\\n")
            
            for dataset, stats in sorted(self.dataset_stats.items()):
                f.write(f"{dataset:<30} {stats.total_variants:<12,} "
                       f"{stats.well_imputed:<15,} {stats.well_imputed_rate:>7.1%}\n")
            
            f.write("\\n")
            
            # MAF distribution
            f.write("MAF Distribution (All Datasets Combined)\\n")
            f.write("-" * 40 + "\\n")
            f.write(f"{'MAF Range':<15} {'Count':<12} {'Percentage':<10}\n")
            
            combined_maf = defaultdict(int)
            for stats in self.dataset_stats.values():
                for bin_name, count in stats.maf_bins.items():
                    combined_maf[bin_name] += count
            
            for _, _, bin_name in self.maf_bins:
                count = combined_maf.get(bin_name, 0)
                pct = 100 * count / total_all if total_all > 0 else 0
                f.write(f"{bin_name:<15} {count:<12,} {pct:>9.1f}%\n")
            
            f.write("\\n")
            
            # R² distribution
            f.write("R² Distribution (All Datasets Combined)\\n")
            f.write("-" * 40 + "\\n")
            f.write(f"{'R² Range':<15} {'Count':<12} {'Percentage':<10}\n")
            
            combined_rsq = defaultdict(int)
            for stats in self.dataset_stats.values():
                for bin_name, count in stats.rsq_bins.items():
                    combined_rsq[bin_name] += count
            
            for _, _, bin_name in self.rsq_bins:
                count = combined_rsq.get(bin_name, 0)
                pct = 100 * count / total_all if total_all > 0 else 0
                f.write(f"{bin_name:<15} {count:<12,} {pct:>9.1f}%\n")
    
    def generate_variant_list(self, output_file: Path, dataset: Optional[str] = None) -> None:
        """Generate list of well-imputed variant IDs."""
        logger.info(f"Generating well-imputed variant list: {output_file}")
        
        with open(output_file, 'w') as f:
            if dataset:
                # Single dataset
                if dataset in self.variants_by_dataset:
                    for variant_id, stats in self.variants_by_dataset[dataset].items():
                        if stats.rsq >= self.rsq_threshold:
                            f.write(f"{variant_id}\n")
            else:
                # All datasets - union of well-imputed variants
                all_well_imputed = set()
                for dataset_variants in self.variants_by_dataset.values():
                    for variant_id, stats in dataset_variants.items():
                        if stats.rsq >= self.rsq_threshold:
                            all_well_imputed.add(variant_id)
                
                for variant_id in sorted(all_well_imputed):
                    f.write(f"{variant_id}\n")
    
    def export_json(self, output_file: Path) -> None:
        """Export statistics as JSON."""
        logger.info(f"Exporting statistics to JSON: {output_file}")
        
        data = {
            'parameters': {
                'rsq_threshold': self.rsq_threshold,
                'maf_threshold': self.maf_threshold
            },
            'datasets': {}
        }
        
        for dataset, stats in self.dataset_stats.items():
            data['datasets'][dataset] = {
                'total_variants': stats.total_variants,
                'well_imputed': stats.well_imputed,
                'genotyped': stats.genotyped,
                'well_imputed_rate': stats.well_imputed_rate,
                'genotyped_rate': stats.genotyped_rate,
                'maf_distribution': dict(stats.maf_bins),
                'rsq_distribution': dict(stats.rsq_bins)
            }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reports on well-imputed variants from imputation results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--info-files",
        required=True,
        help="Comma-separated list of info files to process"
    )
    
    parser.add_argument(
        "--datasets",
        required=True,
        help="Comma-separated list of dataset names"
    )
    
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output file prefix for reports"
    )
    
    parser.add_argument(
        "--rsq-threshold",
        type=float,
        default=0.8,
        help="R² threshold for well-imputed variants"
    )
    
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help="Minimum MAF threshold"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Export statistics as JSON"
    )
    
    parser.add_argument(
        "--variant-list",
        action="store_true",
        help="Generate list of well-imputed variant IDs"
    )
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Parse input lists
    info_files = [Path(f.strip()) for f in args.info_files.split(',')]
    datasets = [d.strip() for d in args.datasets.split(',')]
    
    if len(info_files) != len(datasets):
        logger.error(f"Number of info files ({len(info_files)}) must match datasets ({len(datasets)})")
        return 1
    
    # Initialize reporter
    reporter = WellImputedReporter(args.rsq_threshold, args.maf_threshold)
    
    try:
        # Process all files
        for info_file, dataset in zip(info_files, datasets):
            reporter.process_info_file(info_file, dataset)
        
        # Generate reports
        output_prefix = Path(args.output_prefix)
        reporter.generate_summary_report(output_prefix.with_suffix('.summary.txt'))
        
        if args.variant_list:
            reporter.generate_variant_list(output_prefix.with_suffix('.variants.txt'))
        
        if args.json:
            reporter.export_json(output_prefix.with_suffix('.stats.json'))
        
        logger.info("Done!")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())