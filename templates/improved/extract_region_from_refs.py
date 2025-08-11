#!/usr/bin/env python3
"""
Improved version of extract_region_from_refs.py
Extracts specific genomic regions from reference panel files.
"""

import argparse
import sys
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Set
import logging
import tempfile
import shutil
from dataclasses import dataclass

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class GenomicRegion:
    """Represents a genomic region."""
    chromosome: str
    start: int
    end: int
    
    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"
    
    @classmethod
    def from_string(cls, region_str: str) -> 'GenomicRegion':
        """Parse region from string format (chr:start-end)."""
        if ':' not in region_str or '-' not in region_str:
            raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")
        
        chrom, positions = region_str.split(':')
        start_str, end_str = positions.split('-')
        
        return cls(
            chromosome=chrom,
            start=int(start_str),
            end=int(end_str)
        )


class ReferenceExtractor:
    """Extract regions from reference panel files."""
    
    def __init__(self, bcftools_path: str = "bcftools", tabix_path: str = "tabix"):
        self.bcftools_path = bcftools_path
        self.tabix_path = tabix_path
        self._check_tools()
    
    def _check_tools(self) -> None:
        """Check if required tools are available."""
        tools = [
            (self.bcftools_path, ["--version"]),
            (self.tabix_path, ["--version"])
        ]
        
        for tool, args in tools:
            try:
                result = subprocess.run(
                    [tool] + args,
                    capture_output=True,
                    text=True,
                    check=False
                )
                if result.returncode != 0 and tool != self.tabix_path:
                    raise RuntimeError(f"{tool} not working properly")
                logger.debug(f"{tool} is available")
            except FileNotFoundError:
                raise RuntimeError(f"{tool} not found. Please install it or provide the full path.")
    
    def check_index(self, vcf_file: Path) -> bool:
        """Check if VCF file has an index."""
        # Check for different index formats
        index_extensions = ['.tbi', '.csi', '.idx']
        
        for ext in index_extensions:
            if Path(str(vcf_file) + ext).exists():
                return True
        
        return False
    
    def create_index(self, vcf_file: Path, force: bool = False) -> None:
        """Create index for VCF file if it doesn't exist."""
        if self.check_index(vcf_file) and not force:
            logger.debug(f"Index already exists for {vcf_file}")
            return
        
        logger.info(f"Creating index for {vcf_file}")
        
        # Determine if file is compressed
        is_compressed = str(vcf_file).endswith(('.gz', '.bgz'))
        
        if is_compressed:
            # Use tabix for compressed files
            cmd = [self.tabix_path, "-p", "vcf", str(vcf_file)]
        else:
            # Use bcftools for uncompressed files
            cmd = [self.bcftools_path, "index", str(vcf_file)]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Index created successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create index: {e.stderr}")
            raise
    
    def extract_region(
        self,
        input_file: Path,
        region: GenomicRegion,
        output_file: Path,
        samples: Optional[List[str]] = None,
        exclude_samples: bool = False,
        min_ac: Optional[int] = None,
        max_ac: Optional[int] = None,
        min_af: Optional[float] = None,
        max_af: Optional[float] = None,
        keep_only_pass: bool = False,
        output_type: str = "z"
    ) -> None:
        """
        Extract a genomic region from a VCF/BCF file.
        
        Args:
            input_file: Input VCF/BCF file
            region: Genomic region to extract
            output_file: Output file path
            samples: List of samples to include/exclude
            exclude_samples: If True, exclude listed samples instead of including
            min_ac: Minimum allele count
            max_ac: Maximum allele count
            min_af: Minimum allele frequency
            max_af: Maximum allele frequency
            keep_only_pass: Keep only PASS variants
            output_type: Output type (v=VCF, z=compressed VCF, b=BCF, u=uncompressed BCF)
        """
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        # Create index if needed
        self.create_index(input_file)
        
        # Build bcftools command
        cmd = [
            self.bcftools_path, "view",
            "-r", str(region),
            "-O", output_type,
            "-o", str(output_file)
        ]
        
        # Add sample filtering
        if samples:
            samples_file = None
            if len(samples) > 10:  # Use file for many samples
                samples_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
                samples_file.write('\n'.join(samples))
                samples_file.close()
                
                if exclude_samples:
                    cmd.extend(["-S", f"^{samples_file.name}"])
                else:
                    cmd.extend(["-S", samples_file.name])
            else:
                # Use command line for few samples
                sample_list = ','.join(samples)
                if exclude_samples:
                    cmd.extend(["-s", f"^{sample_list}"])
                else:
                    cmd.extend(["-s", sample_list])
        
        # Add filtering options
        filters = []
        
        if min_ac is not None:
            filters.append(f"AC>={min_ac}")
        if max_ac is not None:
            filters.append(f"AC<={max_ac}")
        if min_af is not None:
            filters.append(f"AF>={min_af}")
        if max_af is not None:
            filters.append(f"AF<={max_af}")
        
        if filters:
            cmd.extend(["-i", " && ".join(filters)])
        
        if keep_only_pass:
            cmd.extend(["-f", "PASS"])
        
        # Add input file
        cmd.append(str(input_file))
        
        # Execute command
        logger.info(f"Extracting region {region} from {input_file}")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Create index for output if compressed
            if output_type in ["z", "b"]:
                self.create_index(output_file)
            
            # Count variants in output
            count_cmd = [self.bcftools_path, "view", "-H", str(output_file)]
            count_result = subprocess.run(count_cmd, capture_output=True, text=True)
            variant_count = len(count_result.stdout.strip().split('\n')) if count_result.stdout else 0
            
            logger.info(f"Extracted {variant_count:,} variants to {output_file}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Extraction failed: {e.stderr}")
            raise
        finally:
            # Clean up temp files
            if samples and len(samples) > 10 and 'samples_file' in locals():
                Path(samples_file.name).unlink(missing_ok=True)
    
    def extract_multiple_regions(
        self,
        input_file: Path,
        regions: List[GenomicRegion],
        output_dir: Path,
        merge_output: bool = False,
        **kwargs
    ) -> List[Path]:
        """
        Extract multiple regions from a reference file.
        
        Args:
            input_file: Input VCF/BCF file
            regions: List of regions to extract
            output_dir: Output directory
            merge_output: If True, merge all regions into a single file
            **kwargs: Additional arguments passed to extract_region
        
        Returns:
            List of output file paths
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files = []
        
        if merge_output:
            # Extract all regions to temporary files then merge
            temp_files = []
            
            for region in regions:
                temp_file = output_dir / f"temp_{region.chromosome}_{region.start}_{region.end}.vcf.gz"
                self.extract_region(input_file, region, temp_file, **kwargs)
                temp_files.append(temp_file)
            
            # Merge files
            merged_file = output_dir / "merged_regions.vcf.gz"
            self.merge_vcf_files(temp_files, merged_file)
            output_files.append(merged_file)
            
            # Clean up temp files
            for temp_file in temp_files:
                temp_file.unlink(missing_ok=True)
                Path(str(temp_file) + ".tbi").unlink(missing_ok=True)
        else:
            # Extract each region to separate file
            for region in regions:
                output_file = output_dir / f"region_{region.chromosome}_{region.start}_{region.end}.vcf.gz"
                self.extract_region(input_file, region, output_file, **kwargs)
                output_files.append(output_file)
        
        return output_files
    
    def merge_vcf_files(self, input_files: List[Path], output_file: Path) -> None:
        """Merge multiple VCF files."""
        logger.info(f"Merging {len(input_files)} files to {output_file}")
        
        cmd = [
            self.bcftools_path, "concat",
            "-a",  # Allow overlaps
            "-O", "z",
            "-o", str(output_file)
        ] + [str(f) for f in input_files]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.create_index(output_file)
            logger.info(f"Merged successfully to {output_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Merge failed: {e.stderr}")
            raise
    
    def get_samples(self, vcf_file: Path) -> List[str]:
        """Get list of samples from VCF file."""
        cmd = [self.bcftools_path, "query", "-l", str(vcf_file)]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            samples = result.stdout.strip().split('\n')
            return [s for s in samples if s]
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to get samples: {e.stderr}")
            raise
    
    def get_variant_count(self, vcf_file: Path, region: Optional[GenomicRegion] = None) -> int:
        """Get variant count in VCF file or region."""
        cmd = [self.bcftools_path, "view", "-H"]
        
        if region:
            cmd.extend(["-r", str(region)])
        
        cmd.append(str(vcf_file))
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            return len(result.stdout.strip().split('\n')) if result.stdout else 0
        except subprocess.CalledProcessError:
            return 0


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract genomic regions from reference panel files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "input_file",
        type=Path,
        help="Input VCF/BCF file"
    )
    
    parser.add_argument(
        "region",
        help="Region to extract (format: chr:start-end) or file with regions"
    )
    
    parser.add_argument(
        "output",
        type=Path,
        help="Output file or directory"
    )
    
    parser.add_argument(
        "--samples",
        help="Comma-separated list of samples to include, or file with sample list"
    )
    
    parser.add_argument(
        "--exclude-samples",
        action="store_true",
        help="Exclude listed samples instead of including them"
    )
    
    parser.add_argument(
        "--min-ac",
        type=int,
        help="Minimum allele count"
    )
    
    parser.add_argument(
        "--max-ac",
        type=int,
        help="Maximum allele count"
    )
    
    parser.add_argument(
        "--min-af",
        type=float,
        help="Minimum allele frequency"
    )
    
    parser.add_argument(
        "--max-af",
        type=float,
        help="Maximum allele frequency"
    )
    
    parser.add_argument(
        "--keep-only-pass",
        action="store_true",
        help="Keep only variants with FILTER=PASS"
    )
    
    parser.add_argument(
        "--output-type",
        choices=["v", "z", "b", "u"],
        default="z",
        help="Output type: v=VCF, z=compressed VCF, b=BCF, u=uncompressed BCF"
    )
    
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Merge multiple regions into single output file"
    )
    
    parser.add_argument(
        "--bcftools",
        default="bcftools",
        help="Path to bcftools executable"
    )
    
    parser.add_argument(
        "--tabix",
        default="tabix",
        help="Path to tabix executable"
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
    
    # Initialize extractor
    extractor = ReferenceExtractor(args.bcftools, args.tabix)
    
    try:
        # Parse regions
        regions = []
        
        # Check if region is a file
        region_path = Path(args.region)
        if region_path.exists():
            # Read regions from file
            with open(region_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        regions.append(GenomicRegion.from_string(line))
            logger.info(f"Read {len(regions)} regions from {region_path}")
        else:
            # Single region from command line
            regions.append(GenomicRegion.from_string(args.region))
        
        # Parse samples if provided
        samples = None
        if args.samples:
            samples_path = Path(args.samples)
            if samples_path.exists():
                # Read samples from file
                with open(samples_path, 'r') as f:
                    samples = [line.strip() for line in f if line.strip()]
                logger.info(f"Read {len(samples)} samples from {samples_path}")
            else:
                # Samples from command line
                samples = [s.strip() for s in args.samples.split(',')]
        
        # Extract regions
        if len(regions) == 1:
            # Single region extraction
            extractor.extract_region(
                args.input_file,
                regions[0],
                args.output,
                samples=samples,
                exclude_samples=args.exclude_samples,
                min_ac=args.min_ac,
                max_ac=args.max_ac,
                min_af=args.min_af,
                max_af=args.max_af,
                keep_only_pass=args.keep_only_pass,
                output_type=args.output_type
            )
        else:
            # Multiple regions
            output_dir = args.output if args.output.is_dir() else args.output.parent
            
            output_files = extractor.extract_multiple_regions(
                args.input_file,
                regions,
                output_dir,
                merge_output=args.merge,
                samples=samples,
                exclude_samples=args.exclude_samples,
                min_ac=args.min_ac,
                max_ac=args.max_ac,
                min_af=args.min_af,
                max_af=args.max_af,
                keep_only_pass=args.keep_only_pass,
                output_type=args.output_type
            )
            
            logger.info(f"Extracted {len(output_files)} files to {output_dir}")
        
        logger.info("Done!")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())