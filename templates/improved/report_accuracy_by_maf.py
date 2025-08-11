#!/usr/bin/env python3
"""
Improved version of report_accuracy_by_maf.py
Analyzes imputation accuracy stratified by minor allele frequency (MAF).
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, NamedTuple
import logging
from collections import defaultdict
import numpy as np
import json
from dataclasses import dataclass, field

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MAFBin(NamedTuple):
    """Definition of a MAF bin."""
    min_maf: float
    max_maf: float
    label: str


@dataclass
class AccuracyMetrics:
    """Accuracy metrics for variants."""
    count: int = 0
    rsq_sum: float = 0.0
    rsq_squared_sum: float = 0.0
    empirical_rsq_sum: float = 0.0
    empirical_rsq_squared_sum: float = 0.0
    concordance_sum: float = 0.0
    well_imputed_count: int = 0
    genotyped_count: int = 0
    
    @property
    def mean_rsq(self) -> float:
        """Calculate mean R²."""
        return self.rsq_sum / self.count if self.count > 0 else 0.0
    
    @property
    def std_rsq(self) -> float:
        """Calculate standard deviation of R²."""
        if self.count <= 1:
            return 0.0
        mean = self.mean_rsq
        variance = (self.rsq_squared_sum / self.count) - (mean ** 2)
        return np.sqrt(max(0, variance))  # Avoid negative variance due to rounding
    
    @property
    def mean_empirical_rsq(self) -> float:
        """Calculate mean empirical R²."""
        return self.empirical_rsq_sum / self.count if self.count > 0 else 0.0
    
    @property
    def mean_concordance(self) -> float:
        """Calculate mean concordance."""
        return self.concordance_sum / self.count if self.count > 0 else 0.0
    
    @property
    def well_imputed_rate(self) -> float:
        """Calculate well-imputed rate."""
        return self.well_imputed_count / self.count if self.count > 0 else 0.0
    
    def add_variant(
        self, 
        rsq: float, 
        empirical_rsq: Optional[float] = None,
        concordance: Optional[float] = None,
        is_well_imputed: bool = False,
        is_genotyped: bool = False
    ) -> None:
        """Add a variant to the metrics."""
        self.count += 1
        self.rsq_sum += rsq
        self.rsq_squared_sum += rsq ** 2
        
        if empirical_rsq is not None:
            self.empirical_rsq_sum += empirical_rsq
            self.empirical_rsq_squared_sum += empirical_rsq ** 2
        
        if concordance is not None:
            self.concordance_sum += concordance
        
        if is_well_imputed:
            self.well_imputed_count += 1
        
        if is_genotyped:
            self.genotyped_count += 1


class AccuracyByMAFReporter:
    """Analyze and report imputation accuracy by MAF bins."""
    
    # Standard MAF bins used in imputation studies
    STANDARD_MAF_BINS = [
        MAFBin(0.0, 0.001, "0-0.1%"),
        MAFBin(0.001, 0.005, "0.1-0.5%"),
        MAFBin(0.005, 0.01, "0.5-1%"),
        MAFBin(0.01, 0.02, "1-2%"),
        MAFBin(0.02, 0.05, "2-5%"),
        MAFBin(0.05, 0.1, "5-10%"),
        MAFBin(0.1, 0.2, "10-20%"),
        MAFBin(0.2, 0.5, "20-50%")
    ]
    
    def __init__(
        self, 
        rsq_threshold: float = 0.8,
        maf_bins: Optional[List[MAFBin]] = None
    ):
        self.rsq_threshold = rsq_threshold
        self.maf_bins = maf_bins or self.STANDARD_MAF_BINS
        
        # Store metrics by dataset and MAF bin
        self.metrics: Dict[str, Dict[str, AccuracyMetrics]] = defaultdict(
            lambda: defaultdict(AccuracyMetrics)
        )
        
        # Overall metrics by dataset
        self.overall_metrics: Dict[str, AccuracyMetrics] = defaultdict(AccuracyMetrics)
        
        # Track column indices for current file
        self.column_indices: Dict[str, int] = {}
    
    def process_info_file(self, info_file: Path, dataset_name: str) -> None:
        """Process a single info file and compute accuracy metrics."""
        if not info_file.exists():
            raise FileNotFoundError(f"Info file not found: {info_file}")
        
        logger.info(f"Processing {dataset_name}: {info_file}")
        
        line_count = 0
        variant_count = 0
        
        with open(info_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                
                # Process header
                if "SNP" in line and "Rsq" in line:
                    self.column_indices = self._parse_header(parts)
                    continue
                
                if not self.column_indices:
                    logger.warning(f"No header found in {info_file}")
                    break
                
                line_count += 1
                
                # Process variant
                if self._process_variant(parts, dataset_name):
                    variant_count += 1
        
        logger.info(f"  Processed {variant_count:,} variants from {line_count:,} lines")
    
    def _parse_header(self, header: List[str]) -> Dict[str, int]:
        """Parse header and return column indices."""
        indices = {}
        
        # Map column names (case-insensitive)
        column_map = {
            'snp': ['snp', 'id', 'variant', 'rsid'],
            'rsq': ['rsq', 'r2', 'info'],
            'maf': ['maf', 'af', 'freq'],
            'emprsq': ['emprsq', 'empirical_rsq', 'emp_r2', 'EmpRsq'],
            'concordance': ['concordance', 'conc', 'accuracy'],
            'genotyped': ['genotyped', 'typed', 'is_genotyped']
        }
        
        header_lower = [h.lower() for h in header]
        
        for key, possible_names in column_map.items():
            for i, col in enumerate(header_lower):
                if col in possible_names:
                    indices[key] = i
                    logger.debug(f"Found column '{key}' at index {i} ('{header[i]}')")
                    break
        
        # Verify required columns
        if 'rsq' not in indices:
            raise ValueError("Required column 'Rsq' not found in header")
        if 'maf' not in indices:
            logger.warning("MAF column not found - will treat all variants as MAF=0")
        
        return indices
    
    def _process_variant(self, parts: List[str], dataset: str) -> bool:
        """
        Process a single variant line.
        
        Returns:
            True if variant was successfully processed, False otherwise
        """
        try:
            # Get R²
            rsq_str = parts[self.column_indices['rsq']]
            if rsq_str in ['-', 'NA', '']:
                return False
            rsq = float(rsq_str)
            
            # Get MAF
            maf = 0.0
            if 'maf' in self.column_indices:
                maf_str = parts[self.column_indices['maf']]
                if maf_str not in ['-', 'NA', '']:
                    maf = float(maf_str)
            
            # Skip negative MAF (indicates error)
            if maf < 0:
                return False
            
            # Get empirical R² if available
            emp_rsq = None
            if 'emprsq' in self.column_indices:
                emp_str = parts[self.column_indices['emprsq']]
                if emp_str not in ['-', 'NA', '']:
                    emp_rsq = float(emp_str)
            
            # Get concordance if available
            concordance = None
            if 'concordance' in self.column_indices:
                conc_str = parts[self.column_indices['concordance']]
                if conc_str not in ['-', 'NA', '']:
                    concordance = float(conc_str)
            
            # Check if genotyped
            is_genotyped = False
            if 'genotyped' in self.column_indices:
                genotyped_str = parts[self.column_indices['genotyped']]
                is_genotyped = genotyped_str.lower() in ['true', '1', 'yes']
            
            # Determine if well-imputed
            is_well_imputed = rsq >= self.rsq_threshold
            
            # Add to appropriate MAF bin
            bin_label = self._get_maf_bin(maf)
            self.metrics[dataset][bin_label].add_variant(
                rsq, emp_rsq, concordance, is_well_imputed, is_genotyped
            )
            
            # Add to overall metrics
            self.overall_metrics[dataset].add_variant(
                rsq, emp_rsq, concordance, is_well_imputed, is_genotyped
            )
            
            return True
            
        except (ValueError, IndexError) as e:
            logger.debug(f"Error processing variant: {e}")
            return False
    
    def _get_maf_bin(self, maf: float) -> str:
        """Get the MAF bin label for a given MAF value."""
        for maf_bin in self.maf_bins:
            if maf_bin.min_maf <= maf < maf_bin.max_maf:
                return maf_bin.label
        # If MAF is >= 0.5, put in last bin
        if maf >= 0.5:
            return self.maf_bins[-1].label
        return "Unknown"
    
    def generate_report(self, output_file: Path) -> None:
        """Generate detailed accuracy report stratified by MAF."""
        logger.info(f"Generating accuracy report: {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("Imputation Accuracy Report by MAF\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"R² Threshold for well-imputed: {self.rsq_threshold}\n")
            f.write(f"Number of datasets: {len(self.metrics)}\n\n")
            
            # For each dataset
            for dataset in sorted(self.metrics.keys()):
                f.write(f"\nDataset: {dataset}\n")
                f.write("-" * 60 + "\n\n")
                
                # Overall metrics for this dataset
                overall = self.overall_metrics[dataset]
                f.write(f"Overall Statistics:\n")
                f.write(f"  Total variants: {overall.count:,}\n")
                f.write(f"  Mean R²: {overall.mean_rsq:.4f} ± {overall.std_rsq:.4f}\n")
                if overall.empirical_rsq_sum > 0:
                    f.write(f"  Mean Empirical R²: {overall.mean_empirical_rsq:.4f}\n")
                f.write(f"  Well-imputed rate: {overall.well_imputed_rate:.1%}\n")
                if overall.genotyped_count > 0:
                    f.write(f"  Genotyped variants: {overall.genotyped_count:,} "
                           f"({100 * overall.genotyped_count / overall.count:.1f}%)\n")
                f.write("\n")
                
                # Table header
                f.write("Accuracy by MAF:\n")
                f.write(f"{'MAF Range':<12} {'N':<10} {'Mean R²':<12} {'Std R²':<10} "
                       f"{'Well-Imp%':<10}")
                if any(m.empirical_rsq_sum > 0 for m in self.metrics[dataset].values()):
                    f.write(f" {'Emp R²':<10}")
                f.write("\n")
                f.write("-" * 80 + "\n")
                
                # Metrics by MAF bin
                for maf_bin in self.maf_bins:
                    bin_metrics = self.metrics[dataset][maf_bin.label]
                    if bin_metrics.count == 0:
                        continue
                    
                    f.write(f"{maf_bin.label:<12} {bin_metrics.count:<10,} "
                           f"{bin_metrics.mean_rsq:<12.4f} {bin_metrics.std_rsq:<10.4f} "
                           f"{bin_metrics.well_imputed_rate * 100:<10.1f}")
                    
                    if bin_metrics.empirical_rsq_sum > 0:
                        f.write(f" {bin_metrics.mean_empirical_rsq:<10.4f}")
                    
                    f.write("\n")
                
                f.write("\n")
    
    def generate_summary_table(self, output_file: Path) -> None:
        """Generate summary table suitable for publication."""
        logger.info(f"Generating summary table: {output_file}")
        
        with open(output_file, 'w') as f:
            # TSV header
            header = ["Dataset", "MAF_Range", "N_Variants", "Mean_Rsq", "SD_Rsq", 
                     "Well_Imputed_Pct", "Mean_Emp_Rsq"]
            f.write('\t'.join(header) + '\n')
            
            # Data rows
            for dataset in sorted(self.metrics.keys()):
                for maf_bin in self.maf_bins:
                    bin_metrics = self.metrics[dataset][maf_bin.label]
                    if bin_metrics.count == 0:
                        continue
                    
                    row = [
                        dataset,
                        maf_bin.label,
                        str(bin_metrics.count),
                        f"{bin_metrics.mean_rsq:.4f}",
                        f"{bin_metrics.std_rsq:.4f}",
                        f"{bin_metrics.well_imputed_rate * 100:.2f}",
                        f"{bin_metrics.mean_empirical_rsq:.4f}" if bin_metrics.empirical_rsq_sum > 0 else "NA"
                    ]
                    f.write('\t'.join(row) + '\n')
    
    def export_json(self, output_file: Path) -> None:
        """Export metrics as JSON for downstream analysis."""
        logger.info(f"Exporting metrics to JSON: {output_file}")
        
        data = {
            'parameters': {
                'rsq_threshold': self.rsq_threshold,
                'maf_bins': [
                    {'min': b.min_maf, 'max': b.max_maf, 'label': b.label}
                    for b in self.maf_bins
                ]
            },
            'datasets': {}
        }
        
        for dataset in self.metrics:
            dataset_data = {
                'overall': {
                    'n_variants': self.overall_metrics[dataset].count,
                    'mean_rsq': self.overall_metrics[dataset].mean_rsq,
                    'std_rsq': self.overall_metrics[dataset].std_rsq,
                    'well_imputed_rate': self.overall_metrics[dataset].well_imputed_rate,
                    'mean_empirical_rsq': self.overall_metrics[dataset].mean_empirical_rsq
                },
                'by_maf': {}
            }
            
            for maf_bin in self.maf_bins:
                bin_metrics = self.metrics[dataset][maf_bin.label]
                if bin_metrics.count > 0:
                    dataset_data['by_maf'][maf_bin.label] = {
                        'n_variants': bin_metrics.count,
                        'mean_rsq': bin_metrics.mean_rsq,
                        'std_rsq': bin_metrics.std_rsq,
                        'well_imputed_rate': bin_metrics.well_imputed_rate,
                        'mean_empirical_rsq': bin_metrics.mean_empirical_rsq
                    }
            
            data['datasets'][dataset] = dataset_data
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze imputation accuracy stratified by MAF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--info-files",
        required=True,
        help="Comma-separated list of info files"
    )
    
    parser.add_argument(
        "--datasets",
        required=True,
        help="Comma-separated list of dataset names"
    )
    
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output file prefix"
    )
    
    parser.add_argument(
        "--rsq-threshold",
        type=float,
        default=0.8,
        help="R² threshold for well-imputed variants"
    )
    
    parser.add_argument(
        "--custom-maf-bins",
        help="Custom MAF bin boundaries (e.g., '0,0.01,0.05,0.1,0.5')"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Export metrics as JSON"
    )
    
    parser.add_argument(
        "--summary-table",
        action="store_true",
        help="Generate summary table (TSV format)"
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
    
    # Parse custom MAF bins if provided
    maf_bins = None
    if args.custom_maf_bins:
        boundaries = [float(x) for x in args.custom_maf_bins.split(',')]
        maf_bins = []
        for i in range(len(boundaries) - 1):
            label = f"{boundaries[i]*100:.1f}-{boundaries[i+1]*100:.1f}%"
            maf_bins.append(MAFBin(boundaries[i], boundaries[i+1], label))
    
    # Initialize reporter
    reporter = AccuracyByMAFReporter(args.rsq_threshold, maf_bins)
    
    try:
        # Process all files
        for info_file, dataset in zip(info_files, datasets):
            reporter.process_info_file(info_file, dataset)
        
        # Generate reports
        output_prefix = Path(args.output_prefix)
        reporter.generate_report(output_prefix.with_suffix('.report.txt'))
        
        if args.summary_table:
            reporter.generate_summary_table(output_prefix.with_suffix('.summary.tsv'))
        
        if args.json:
            reporter.export_json(output_prefix.with_suffix('.metrics.json'))
        
        logger.info("Done!")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())