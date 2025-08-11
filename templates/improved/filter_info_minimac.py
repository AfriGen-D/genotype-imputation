#!/usr/bin/env python3
"""
Improved version of filter_info_minimac.py
Filters imputation info files based on quality thresholds.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments with improved validation."""
    parser = argparse.ArgumentParser(
        description="Filter imputation info files based on quality thresholds",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--infoFiles",
        required=True,
        help="Comma-separated list of info files to process"
    )
    parser.add_argument(
        "--datasets",
        required=True,
        help="Comma-separated list of dataset names corresponding to info files"
    )
    parser.add_argument(
        "--out_prefix",
        required=True,
        help="Output file prefix for filtered results"
    )
    parser.add_argument(
        "--infoCutoff",
        type=float,
        default=0.3,
        help="Minimum R-squared threshold for well-imputed variants"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    info_files = args.infoFiles.split(',')
    datasets = args.datasets.split(',')
    
    if len(info_files) != len(datasets):
        parser.error(f"Number of info files ({len(info_files)}) must match number of datasets ({len(datasets)})")
    
    # Check if files exist
    for file_path in info_files:
        if not Path(file_path).exists():
            parser.error(f"Info file not found: {file_path}")
    
    return args


class InfoFileProcessor:
    """Process imputation info files with quality filtering."""
    
    def __init__(self, info_cutoff: float, out_prefix: str):
        self.info_cutoff = info_cutoff
        self.out_prefix = out_prefix
        self.header: List[str] = []
        self.info_idx: Optional[int] = None
        self.maf_idx: Optional[int] = None
        self.conc_idx: Optional[int] = None
        self.has_concordance = False
        
    def process_files(self, info_files: List[str], datasets: List[str]) -> Tuple[int, int, int]:
        """
        Process multiple info files and write filtered results.
        
        Returns:
            Tuple of (total_variants, well_imputed_count, concordance_count)
        """
        total_variants = 0
        well_imputed_count = 0
        concordance_count = 0
        
        # Open output files
        with open(f"{self.out_prefix}_well_imputed.tsv", 'w') as well_out, \
             open(f"{self.out_prefix}_well_imputed_snp.tsv", 'w') as snp_out, \
             open(f"{self.out_prefix}_accuracy.tsv", 'w') as acc_out:
            
            for info_file, dataset in zip(info_files, datasets):
                logger.info(f"Processing {dataset}: {info_file}")
                
                file_stats = self._process_single_file(
                    info_file, dataset, well_out, snp_out, acc_out
                )
                
                total_variants += file_stats[0]
                well_imputed_count += file_stats[1]
                concordance_count += file_stats[2]
        
        return total_variants, well_imputed_count, concordance_count
    
    def _process_single_file(
        self, 
        info_file: str, 
        dataset: str,
        well_out, 
        snp_out, 
        acc_out
    ) -> Tuple[int, int, int]:
        """Process a single info file."""
        file_variants = 0
        file_well_imputed = 0
        file_concordance = 0
        
        try:
            with open(info_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    data = line.split('\t')
                    
                    # Process header
                    if "SNP" in line and "Rsq" in line:
                        if not self.header:
                            self._process_header(data, well_out, snp_out, acc_out)
                        continue
                    
                    # Process data line
                    if self.info_idx is None:
                        logger.warning(f"Header not found in {info_file}, skipping...")
                        break
                    
                    stats = self._process_data_line(
                        data, dataset, well_out, snp_out, acc_out
                    )
                    
                    if stats:
                        file_variants += 1
                        if stats[0]:
                            file_well_imputed += 1
                        if stats[1]:
                            file_concordance += 1
                            
        except Exception as e:
            logger.error(f"Error processing {info_file}: {e}")
            raise
        
        logger.info(f"  Processed {file_variants} variants, {file_well_imputed} well-imputed")
        
        return file_variants, file_well_imputed, file_concordance
    
    def _process_header(self, data: List[str], well_out, snp_out, acc_out) -> None:
        """Process header line and set up column indices."""
        self.header = data
        
        try:
            self.info_idx = self.header.index("Rsq")
            self.maf_idx = self.header.index("MAF")
            
            # Write header to output files
            well_out.write('\t'.join(["GROUPS"] + data) + '\n')
            snp_out.write(data[1] + '\n')  # SNP column
            
            # Check for concordance data
            if "EmpRsq" in self.header:
                self.has_concordance = True
                self.conc_idx = self.header.index("EmpRsq")
                acc_out.write('\t'.join(["GROUPS"] + data) + '\n')
                
        except ValueError as e:
            logger.error(f"Required column not found in header: {e}")
            raise
    
    def _process_data_line(
        self, 
        data: List[str], 
        dataset: str,
        well_out, 
        snp_out, 
        acc_out
    ) -> Optional[Tuple[bool, bool]]:
        """
        Process a single data line.
        
        Returns:
            Tuple of (is_well_imputed, has_concordance) or None if skipped
        """
        try:
            # Check MAF filter
            maf = float(data[self.maf_idx])
            if maf < 0.0:
                return None
            
            is_well_imputed = False
            has_concordance = False
            
            # Check R-squared threshold
            rsq_value = data[self.info_idx]
            if rsq_value not in ['-', 'NA', '']:
                rsq = float(rsq_value)
                if rsq >= self.info_cutoff:
                    well_out.write('\t'.join([dataset] + data) + '\n')
                    snp_out.write(data[1] + '\n')  # SNP column
                    is_well_imputed = True
            
            # Check concordance if available
            if self.has_concordance and self.conc_idx is not None:
                conc_value = data[self.conc_idx]
                if conc_value != '-':
                    acc_out.write('\t'.join([dataset] + data) + '\n')
                    has_concordance = True
            
            return is_well_imputed, has_concordance
            
        except (ValueError, IndexError) as e:
            logger.debug(f"Error processing line: {e}")
            return None


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Parse file lists
    info_files = [f.strip() for f in args.infoFiles.split(',')]
    datasets = [d.strip() for d in args.datasets.split(',')]
    
    logger.info(f"Processing {len(info_files)} info files with cutoff {args.infoCutoff}")
    
    # Process files
    processor = InfoFileProcessor(args.infoCutoff, args.out_prefix)
    total, well_imputed, concordance = processor.process_files(info_files, datasets)
    
    # Print summary
    logger.info("=" * 50)
    logger.info(f"Total variants processed: {total:,}")
    logger.info(f"Well-imputed variants (Rsq >= {args.infoCutoff}): {well_imputed:,}")
    if well_imputed > 0:
        logger.info(f"Percentage well-imputed: {100 * well_imputed / total:.2f}%")
    if concordance > 0:
        logger.info(f"Variants with concordance data: {concordance:,}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())