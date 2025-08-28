#!/usr/bin/env python3
"""
Validate that a chunk has sufficient variants for imputation.
Checks both the main region and any buffer regions.
"""

import argparse
import subprocess
import sys
import gzip
import re

def count_variants_in_region(vcf_path, region):
    """Count variants in a specific region of a VCF file"""
    try:
        # Use bcftools to count variants in the region
        cmd = ['bcftools', 'view', '-H', '-r', region, vcf_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Count non-empty lines (each line is a variant)
        variant_count = len([line for line in result.stdout.split('\n') if line.strip()])
        return variant_count
    except subprocess.CalledProcessError as e:
        print(f"Error counting variants: {e}", file=sys.stderr)
        return 0

def parse_vcf_header_for_samples(vcf_path):
    """Get the number of samples from VCF header"""
    try:
        cmd = ['bcftools', 'query', '-l', vcf_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        sample_count = len(result.stdout.strip().split('\n'))
        return sample_count
    except subprocess.CalledProcessError:
        return 0

def validate_chunk_density(vcf_path, ref_vcf_path, chrom, start, end, buffer_size=500000, min_ratio=0.01):
    """
    Validate that a chunk and its buffer regions have sufficient variant density.
    
    Args:
        vcf_path: Path to target VCF file
        ref_vcf_path: Path to reference panel VCF file
        chrom: Chromosome
        start: Start position
        end: End position
        buffer_size: Buffer size on each side (default 500kb)
        min_ratio: Minimum ratio of typed to reference variants
    
    Returns:
        Dictionary with validation results
    """
    
    results = {
        'main_region': {},
        'buffer_before': {},
        'buffer_after': {},
        'combined': {},
        'pass': False,
        'warnings': [],
        'errors': []
    }
    
    # Define regions
    main_region = f"{chrom}:{start}-{end}"
    buffer_before_region = f"{chrom}:{max(1, start - buffer_size)}-{start - 1}" if start > 1 else None
    buffer_after_region = f"{chrom}:{end + 1}-{end + buffer_size}"
    combined_region = f"{chrom}:{max(1, start - buffer_size)}-{end + buffer_size}"
    
    # Count variants in each region
    print(f"Validating chunk {main_region} with {buffer_size}bp buffers...", file=sys.stderr)
    
    # Main region
    target_main = count_variants_in_region(vcf_path, main_region)
    ref_main = count_variants_in_region(ref_vcf_path, main_region)
    ratio_main = target_main / ref_main if ref_main > 0 else 0
    
    results['main_region'] = {
        'region': main_region,
        'target_variants': target_main,
        'ref_variants': ref_main,
        'ratio': ratio_main
    }
    
    # Buffer before (if applicable)
    if buffer_before_region:
        target_before = count_variants_in_region(vcf_path, buffer_before_region)
        ref_before = count_variants_in_region(ref_vcf_path, buffer_before_region)
        ratio_before = target_before / ref_before if ref_before > 0 else 0
        
        results['buffer_before'] = {
            'region': buffer_before_region,
            'target_variants': target_before,
            'ref_variants': ref_before,
            'ratio': ratio_before
        }
    
    # Buffer after
    target_after = count_variants_in_region(vcf_path, buffer_after_region)
    ref_after = count_variants_in_region(ref_vcf_path, buffer_after_region)
    ratio_after = target_after / ref_after if ref_after > 0 else 0
    
    results['buffer_after'] = {
        'region': buffer_after_region,
        'target_variants': target_after,
        'ref_variants': ref_after,
        'ratio': ratio_after
    }
    
    # Combined region (main + buffers)
    target_combined = count_variants_in_region(vcf_path, combined_region)
    ref_combined = count_variants_in_region(ref_vcf_path, combined_region)
    ratio_combined = target_combined / ref_combined if ref_combined > 0 else 0
    
    results['combined'] = {
        'region': combined_region,
        'target_variants': target_combined,
        'ref_variants': ref_combined,
        'ratio': ratio_combined
    }
    
    # Get sample count
    sample_count = parse_vcf_header_for_samples(vcf_path)
    results['sample_count'] = sample_count
    
    # Validation checks
    min_variants_for_phasing = 50  # Eagle typically needs at least 50 variants
    min_variants_for_imputation = 100  # Good imputation needs more variants
    
    # Check main region
    if target_main < min_variants_for_phasing:
        results['errors'].append(f"Main region has only {target_main} variants (minimum: {min_variants_for_phasing})")
    elif target_main < min_variants_for_imputation:
        results['warnings'].append(f"Main region has only {target_main} variants (recommended: {min_variants_for_imputation}+)")
    
    if ratio_main < min_ratio:
        results['errors'].append(f"Main region ratio {ratio_main:.4f} below threshold {min_ratio}")
    
    # Check for empty buffer regions that might cause issues
    if buffer_before_region and target_before == 0 and ref_before > 100:
        results['warnings'].append(f"Buffer before region is empty but has {ref_before} reference variants")
    
    if target_after == 0 and ref_after > 100:
        results['warnings'].append(f"Buffer after region is empty but has {ref_after} reference variants")
    
    # Check combined region
    if ratio_combined < min_ratio:
        results['warnings'].append(f"Combined region ratio {ratio_combined:.4f} below threshold {min_ratio}")
    
    # Special check for the exact issue you encountered
    # If buffer regions have very different densities, warn about potential issues
    if buffer_before_region and ratio_before > 0 and ratio_after == 0:
        results['warnings'].append("Asymmetric variant density: buffer after is empty while buffer before has variants")
    elif buffer_before_region and ratio_before == 0 and ratio_after > 0:
        results['warnings'].append("Asymmetric variant density: buffer before is empty while buffer after has variants")
    
    # Determine if chunk should pass
    # Must have sufficient variants in main region and reasonable overall density
    results['pass'] = (
        len(results['errors']) == 0 and
        target_main >= min_variants_for_phasing and
        ratio_main >= min_ratio
    )
    
    # Additional recommendation for buffer handling
    if not results['pass'] and target_main > 0:
        # If main region has some variants but not enough, suggest alternatives
        if ratio_main < min_ratio:
            results['warnings'].append(f"Consider reducing chunk size or adjusting minRatio parameter")
        if target_main < min_variants_for_phasing:
            results['warnings'].append(f"Consider merging with adjacent chunks to increase variant count")
    
    return results

def format_validation_report(results):
    """Format validation results for output"""
    lines = []
    lines.append("=" * 60)
    lines.append("IMPUTATION CHUNK VALIDATION REPORT")
    lines.append("=" * 60)
    
    # Main region
    main = results['main_region']
    lines.append(f"\nMain Region: {main['region']}")
    lines.append(f"  Target variants: {main['target_variants']}")
    lines.append(f"  Reference variants: {main['ref_variants']}")
    lines.append(f"  Ratio: {main['ratio']:.4f}")
    
    # Buffer regions
    if results['buffer_before']:
        before = results['buffer_before']
        lines.append(f"\nBuffer Before: {before['region']}")
        lines.append(f"  Target variants: {before['target_variants']}")
        lines.append(f"  Reference variants: {before['ref_variants']}")
        lines.append(f"  Ratio: {before['ratio']:.4f}")
    
    after = results['buffer_after']
    lines.append(f"\nBuffer After: {after['region']}")
    lines.append(f"  Target variants: {after['target_variants']}")
    lines.append(f"  Reference variants: {after['ref_variants']}")
    lines.append(f"  Ratio: {after['ratio']:.4f}")
    
    # Combined
    combined = results['combined']
    lines.append(f"\nCombined Region: {combined['region']}")
    lines.append(f"  Target variants: {combined['target_variants']}")
    lines.append(f"  Reference variants: {combined['ref_variants']}")
    lines.append(f"  Ratio: {combined['ratio']:.4f}")
    
    # Sample info
    lines.append(f"\nSample Count: {results['sample_count']}")
    
    # Errors and warnings
    if results['errors']:
        lines.append("\nERRORS:")
        for error in results['errors']:
            lines.append(f"  ✗ {error}")
    
    if results['warnings']:
        lines.append("\nWARNINGS:")
        for warning in results['warnings']:
            lines.append(f"  ⚠ {warning}")
    
    # Final verdict
    lines.append("\n" + "=" * 60)
    if results['pass']:
        lines.append("VERDICT: PASS ✓ - Chunk suitable for imputation")
    else:
        lines.append("VERDICT: FAIL ✗ - Chunk should be skipped or merged")
    lines.append("=" * 60)
    
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(
        description='Validate chunk suitability for imputation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate a chunk with default settings
  %(prog)s --vcf target.vcf.gz --ref-vcf ref.vcf.gz --chrom chr9 --start 45046581 --end 50046580
  
  # Validate with custom buffer and ratio
  %(prog)s --vcf target.vcf.gz --ref-vcf ref.vcf.gz --chrom chr9 --start 45046581 --end 50046580 \\
           --buffer 250000 --min-ratio 0.001
        """
    )
    
    parser.add_argument('--vcf', required=True, help='Target VCF file (phased)')
    parser.add_argument('--ref-vcf', required=True, help='Reference panel VCF file')
    parser.add_argument('--chrom', required=True, help='Chromosome')
    parser.add_argument('--start', type=int, required=True, help='Start position')
    parser.add_argument('--end', type=int, required=True, help='End position')
    parser.add_argument('--buffer', type=int, default=500000, help='Buffer size (default: 500000)')
    parser.add_argument('--min-ratio', type=float, default=0.01, help='Minimum variant ratio (default: 0.01)')
    parser.add_argument('--output', help='Output file for report (default: stdout)')
    parser.add_argument('--status-file', help='Write PASS/FAIL status to this file')
    
    args = parser.parse_args()
    
    # Validate the chunk
    results = validate_chunk_density(
        args.vcf,
        args.ref_vcf,
        args.chrom,
        args.start,
        args.end,
        args.buffer,
        args.min_ratio
    )
    
    # Format report
    report = format_validation_report(results)
    
    # Output report
    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
    else:
        print(report)
    
    # Write status file if requested
    if args.status_file:
        with open(args.status_file, 'w') as f:
            f.write("PASS\n" if results['pass'] else "FAIL\n")
    
    # Exit with appropriate code
    sys.exit(0 if results['pass'] else 1)

if __name__ == '__main__':
    main()