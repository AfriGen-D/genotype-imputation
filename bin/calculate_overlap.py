#!/usr/bin/env python3

import argparse
import sys

def read_map_file(filepath):
    """Read a map file and return a set of positions per chromosome"""
    chrom_positions = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    pos = int(parts[1])
                    if chrom not in chrom_positions:
                        chrom_positions[chrom] = set()
                    chrom_positions[chrom].add(pos)
    return chrom_positions

def calculate_overlap(chunk_map, ref_map, min_ratio, output_file, chunk_id, region):
    """Calculate overlap between chunk and reference panel"""
    
    # Read map files
    chunk_positions = read_map_file(chunk_map)
    ref_positions = read_map_file(ref_map)
    
    # Calculate overlap for each chromosome
    results = []
    total_chunk_vars = 0
    total_ref_vars = 0
    total_overlap = 0
    
    for chrom in set(list(chunk_positions.keys()) + list(ref_positions.keys())):
        chunk_vars = chunk_positions.get(chrom, set())
        ref_vars = ref_positions.get(chrom, set())
        overlap = chunk_vars.intersection(ref_vars)
        
        total_chunk_vars += len(chunk_vars)
        total_ref_vars += len(ref_vars)
        total_overlap += len(overlap)
        
        if len(ref_vars) > 0:
            ratio = len(overlap) / len(ref_vars)
        else:
            ratio = 0.0
            
        results.append({
            'chrom': chrom,
            'chunk_vars': len(chunk_vars),
            'ref_vars': len(ref_vars),
            'overlap': len(overlap),
            'ratio': ratio
        })
    
    # Calculate overall ratio
    if total_ref_vars > 0:
        overall_ratio = total_overlap / total_ref_vars
    else:
        overall_ratio = 0.0
    
    # Write output file
    with open(output_file, 'w') as f:
        f.write(f"Chunk: {chunk_id}\n")
        f.write(f"Region: {region}\n")
        f.write(f"Chunk variants: {total_chunk_vars}\n")
        f.write(f"Reference variants in region: {total_ref_vars}\n")
        f.write(f"Overlapping variants: {total_overlap}\n")
        f.write(f"Overlap ratio: {overall_ratio:.4f}\n")
        f.write(f"Min ratio threshold: {min_ratio}\n")
        f.write("\n")
        
        # Per-chromosome details
        f.write("Per-chromosome breakdown:\n")
        for result in results:
            f.write(f"  {result['chrom']}: chunk={result['chunk_vars']}, ref={result['ref_vars']}, "
                   f"overlap={result['overlap']}, ratio={result['ratio']:.4f}\n")
        
        f.write("\n")
        if overall_ratio < min_ratio:
            f.write(f"WARNING: Overlap ratio {overall_ratio:.4f} is below minimum threshold {min_ratio}\n")
            status = "FAIL"
        else:
            f.write(f"PASS: Overlap ratio meets threshold\n")
            status = "PASS"
    
    # Print summary to stdout
    print(f"Chunk {chunk_id} overlap check: {status}")
    print(f"  Chunk variants: {total_chunk_vars}")
    print(f"  Reference variants: {total_ref_vars}")
    print(f"  Overlapping: {total_overlap}")
    print(f"  Ratio: {overall_ratio:.4f} (threshold: {min_ratio})")
    
    # Write status file for workflow filtering
    # Need both good ratio AND minimum absolute overlap count
    min_overlap_count = 50  # Minimum number of overlapping variants for Eagle
    with open(output_file.replace('.overlap.txt', '.overlap.status'), 'w') as f:
        if overall_ratio >= min_ratio and total_overlap >= min_overlap_count:
            f.write("PASS\n")
        else:
            f.write("FAIL\n")
    
    # Handle different overlap scenarios
    if overall_ratio == 0:
        if total_chunk_vars == 0:
            print(f"WARNING: No variants found in chunk {chunk_id}")
        elif total_ref_vars == 0:
            print(f"WARNING: No reference variants found in region {region}")
        else:
            print(f"WARNING: No overlapping variants between chunk and reference panel")
        print(f"  This chunk will be skipped from phasing/imputation")
        # Don't exit with error - let workflow filter handle it
    elif overall_ratio < min_ratio:
        print(f"WARNING: Low overlap ratio - phasing/imputation may fail")
        print(f"  This chunk will be skipped from phasing/imputation")
    elif total_overlap < 100:  # Additional check for minimum absolute overlap
        print(f"WARNING: Too few overlapping variants ({total_overlap}) for reliable phasing")
        print(f"  Consider this chunk may fail in Eagle phasing")
    
    # Check if we have enough overlap for phasing
    # Eagle typically needs at least 100 overlapping variants for good phasing
    if total_overlap < 50:
        print(f"WARNING: Very low overlap count ({total_overlap} variants)")
        print(f"  Eagle phasing is likely to fail with genetic map issues")
    
    return overall_ratio

def main():
    parser = argparse.ArgumentParser(description='Calculate overlap between chunk and reference panel')
    parser.add_argument('--chunk-map', required=True, help='Chunk map file')
    parser.add_argument('--ref-map', required=True, help='Reference panel map file')
    parser.add_argument('--min-ratio', type=float, default=0.001, help='Minimum overlap ratio')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--chunk-id', required=True, help='Chunk identifier')
    parser.add_argument('--region', default='', help='Genomic region')
    
    args = parser.parse_args()
    
    calculate_overlap(
        args.chunk_map,
        args.ref_map,
        args.min_ratio,
        args.output,
        args.chunk_id,
        args.region
    )

if __name__ == '__main__':
    main()