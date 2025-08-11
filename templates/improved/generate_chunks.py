#!/usr/bin/env python
"""
Improved version of generate_chunks.py
Generates genomic chunks for parallel processing of imputation.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import logging
from dataclasses import dataclass

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class GenomicRegion:
    """Represents a genomic region with chromosome and position range."""
    chromosome: str
    start: int
    end: int
    
    def __str__(self) -> str:
        return f"{self.chromosome},{self.start},{self.end}"
    
    def to_bed_format(self) -> str:
        """Return BED format (0-based, half-open)."""
        return f"{self.chromosome}\t{self.start-1}\t{self.end}"


class ChunkGenerator:
    """Generate genomic chunks for parallel processing."""
    
    def __init__(self, chunk_size: int):
        self.chunk_size = chunk_size
        self.chromosome_data: Dict[str, List[int]] = {}
        
    def read_map_file(self, map_file: Path) -> None:
        """
        Read chromosome positions from map file.
        Expected format: chromosome<tab>position
        """
        if not map_file.exists():
            raise FileNotFoundError(f"Map file not found: {map_file}")
        
        logger.info(f"Reading map file: {map_file}")
        
        with open(map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    parts = line.split('\t')
                    if len(parts) < 2:
                        logger.warning(f"Line {line_num}: Invalid format, skipping")
                        continue
                    
                    chrom = parts[0]
                    pos = int(parts[1])
                    
                    if chrom not in self.chromosome_data:
                        self.chromosome_data[chrom] = []
                    self.chromosome_data[chrom].append(pos)
                    
                except ValueError as e:
                    logger.warning(f"Line {line_num}: Error parsing position - {e}")
                    continue
        
        # Sort positions for each chromosome
        for chrom in self.chromosome_data:
            self.chromosome_data[chrom].sort()
            logger.info(f"  Chromosome {chrom}: {len(self.chromosome_data[chrom]):,} positions")
    
    def generate_chunks(
        self, 
        chromosomes: Optional[List[str]] = None,
        specific_region: Optional[str] = None
    ) -> List[GenomicRegion]:
        """
        Generate chunks for specified chromosomes or regions.
        
        Args:
            chromosomes: List of chromosomes to process (None = all)
            specific_region: Specific region in format "chr:start-end"
        
        Returns:
            List of GenomicRegion objects
        """
        chunks = []
        
        if specific_region:
            # Process specific region
            chunks = self._process_specific_region(specific_region)
        else:
            # Process chromosomes
            if chromosomes is None:
                chromosomes = sorted(self.chromosome_data.keys(), key=self._sort_chromosome_key)
            else:
                # Validate requested chromosomes
                available = set(self.chromosome_data.keys())
                requested = set(chromosomes)
                missing = requested - available
                if missing:
                    logger.warning(f"Chromosomes not found in data: {missing}")
                chromosomes = sorted(requested & available, key=self._sort_chromosome_key)
            
            for chrom in chromosomes:
                chunks.extend(self._generate_chunks_for_chromosome(chrom))
        
        return chunks
    
    def _generate_chunks_for_chromosome(self, chrom: str) -> List[GenomicRegion]:
        """Generate chunks for a single chromosome."""
        if chrom not in self.chromosome_data:
            logger.warning(f"Chromosome {chrom} not found in data")
            return []
        
        positions = self.chromosome_data[chrom]
        if not positions:
            return []
        
        min_pos = min(positions)
        max_pos = max(positions)
        
        # Align start to chunk boundary
        chunk_start = min_pos - (min_pos % 10) + 1
        
        chunks = []
        while chunk_start <= max_pos:
            chunk_end = min(chunk_start + self.chunk_size - 1, max_pos)
            chunks.append(GenomicRegion(chrom, chunk_start, chunk_end))
            chunk_start += self.chunk_size
        
        return chunks
    
    def _process_specific_region(self, region: str) -> List[GenomicRegion]:
        """
        Process a specific region string (e.g., "chr1:1000000-2000000").
        Can be a single region or comma-separated list.
        """
        chunks = []
        regions = region.split(',')
        
        for reg in regions:
            try:
                if ':' not in reg:
                    logger.warning(f"Invalid region format: {reg}")
                    continue
                
                chrom, pos_range = reg.split(':')
                start_str, end_str = pos_range.split('-')
                region_start = int(start_str)
                region_end = int(end_str)
                
                # Generate chunks within this region
                chunk_start = region_start
                while chunk_start <= region_end:
                    chunk_end = min(chunk_start + self.chunk_size - 1, region_end)
                    
                    # Only add if chromosome exists in data
                    if chrom in self.chromosome_data:
                        # Verify chunk overlaps with actual data
                        positions = self.chromosome_data[chrom]
                        min_pos = min(positions)
                        max_pos = max(positions)
                        
                        if chunk_start <= max_pos and chunk_end >= min_pos:
                            # Adjust to actual data boundaries
                            adjusted_start = max(chunk_start, min_pos)
                            adjusted_end = min(chunk_end, max_pos)
                            chunks.append(GenomicRegion(chrom, adjusted_start, adjusted_end))
                    
                    chunk_start += self.chunk_size
                    
            except Exception as e:
                logger.warning(f"Error processing region {reg}: {e}")
                continue
        
        return chunks
    
    @staticmethod
    def _sort_chromosome_key(chrom: str) -> Tuple:
        """
        Sort chromosomes in natural order (1, 2, ..., 10, ..., 20, 21, 22, X, Y, MT).
        Handles both 'chr1' and '1' formats.
        """
        # Remove 'chr' prefix if present
        chrom_clean = chrom.replace('chr', '')
        
        # Try to convert to integer for numeric chromosomes
        try:
            return (0, int(chrom_clean))
        except ValueError:
            # Non-numeric chromosomes (X, Y, MT, etc.)
            order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
            return (1, order.get(chrom_clean.upper(), 100))
    
    def write_chunks(self, chunks: List[GenomicRegion], output_file: Path) -> None:
        """Write chunks to output file."""
        logger.info(f"Writing {len(chunks)} chunks to {output_file}")
        
        with open(output_file, 'w') as f:
            for chunk in chunks:
                f.write(str(chunk) + '\n')
    
    def get_statistics(self) -> Dict:
        """Get statistics about the loaded data."""
        stats = {
            'total_chromosomes': len(self.chromosome_data),
            'total_positions': sum(len(positions) for positions in self.chromosome_data.values()),
            'chromosomes': {}
        }
        
        for chrom, positions in self.chromosome_data.items():
            if positions:
                stats['chromosomes'][chrom] = {
                    'count': len(positions),
                    'min': min(positions),
                    'max': max(positions),
                    'range': max(positions) - min(positions)
                }
        
        return stats


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate genomic chunks for parallel imputation processing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "map_file",
        type=Path,
        help="Input map file with chromosome and position columns (tab-separated)"
    )
    
    parser.add_argument(
        "output_file",
        type=Path,
        help="Output file for chunk definitions"
    )
    
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000000,
        help="Size of each chunk in base pairs"
    )
    
    parser.add_argument(
        "--chromosomes",
        type=str,
        help="Comma-separated list of chromosomes to process (default: all)"
    )
    
    parser.add_argument(
        "--region",
        type=str,
        help="Specific region(s) to process (e.g., 'chr1:1000000-2000000')"
    )
    
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print statistics about the input data"
    )
    
    parser.add_argument(
        "--bed-format",
        action="store_true",
        help="Output in BED format (0-based, half-open intervals)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_arguments()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Initialize chunk generator
    generator = ChunkGenerator(args.chunk_size)
    
    try:
        # Read input map file
        generator.read_map_file(args.map_file)
        
        # Print statistics if requested
        if args.stats:
            stats = generator.get_statistics()
            logger.info("=" * 50)
            logger.info("Data Statistics:")
            logger.info(f"Total chromosomes: {stats['total_chromosomes']}")
            logger.info(f"Total positions: {stats['total_positions']:,}")
            for chrom, chrom_stats in sorted(stats['chromosomes'].items()):
                logger.info(f"  {chrom}: {chrom_stats['count']:,} positions, "
                          f"range {chrom_stats['min']:,}-{chrom_stats['max']:,}")
            logger.info("=" * 50)
        
        # Parse chromosomes if specified
        chromosomes = None
        if args.chromosomes:
            chromosomes = [c.strip() for c in args.chromosomes.split(',')]
        
        # Generate chunks
        chunks = generator.generate_chunks(chromosomes, args.region)
        
        if not chunks:
            logger.warning("No chunks generated!")
            return 1
        
        logger.info(f"Generated {len(chunks)} chunks")
        
        # Write output
        if args.bed_format:
            logger.info(f"Writing chunks in BED format to {args.output_file}")
            with open(args.output_file, 'w') as f:
                for chunk in chunks:
                    f.write(chunk.to_bed_format() + '\n')
        else:
            generator.write_chunks(chunks, args.output_file)
        
        logger.info("Done!")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())