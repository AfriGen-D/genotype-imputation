#!/usr/bin/env python3

import re
import sys
from collections import defaultdict

print("=" * 80)
print("CHROMOSOME MISMATCH ANALYSIS")
print("=" * 80)
print()

# Parse the Nextflow log to find processing patterns
log_file = ".nextflow.log"

# Track statistics per chromosome
chr_stats = defaultdict(lambda: {
    'total_chunks': 0,
    'passed_mismatch': 0,
    'failed_mismatch': 0,
    'reached_phasing': 0,
    'reached_imputation': 0,
    'completed_imputation': 0
})

# Expected chunks per chromosome based on size
expected_chunks = {
    'chr1': 10, 'chr2': 10, 'chr3': 8, 'chr4': 8, 'chr5': 7,
    'chr6': 7, 'chr7': 7, 'chr8': 6, 'chr9': 6, 'chr10': 6,
    'chr11': 6, 'chr12': 6, 'chr13': 4, 'chr14': 4, 'chr15': 4,
    'chr16': 4, 'chr17': 4, 'chr18': 3, 'chr19': 3, 'chr20': 3,
    'chr21': 2, 'chr22': 2
}

print("Reading Nextflow log file...")
with open(log_file, 'r') as f:
    for line in f:
        # Extract chromosome from chunk IDs
        chunk_match = re.search(r'(chr[0-9XY]+)_\d+_\d+', line)
        if chunk_match:
            chrom = chunk_match.group(1)
            
            # Track different stages
            if 'GENERATE_CHUNKS_VCF' in line and 'Task completed' in line:
                chr_stats[chrom]['total_chunks'] += 1
            elif 'CHECK_MISMATCH' in line:
                if 'Task completed' in line:
                    chr_stats[chrom]['passed_mismatch'] += 1
                elif 'FAILED' in line or 'ERROR' in line:
                    chr_stats[chrom]['failed_mismatch'] += 1
            elif 'EAGLE_PHASING' in line and 'Submitted' in line:
                chr_stats[chrom]['reached_phasing'] += 1
            elif 'IMPUTE_MINIMAC4' in line:
                if 'Submitted' in line:
                    chr_stats[chrom]['reached_imputation'] += 1
                elif 'Task completed' in line:
                    chr_stats[chrom]['completed_imputation'] += 1

# Calculate mismatch issues
print("\nChromosome Analysis:")
print("-" * 80)
print(f"{'Chromosome':<12} {'Expected':<10} {'Imputed':<10} {'Missing':<10} {'Loss %':<10}")
print("-" * 80)

problematic_chromosomes = []
total_expected = 0
total_completed = 0

for chrom in sorted(expected_chunks.keys(), key=lambda x: int(x[3:]) if x[3:].isdigit() else 99):
    expected = expected_chunks.get(chrom, 0)
    completed = chr_stats[chrom]['completed_imputation']
    missing = expected - completed
    loss_pct = (missing / expected * 100) if expected > 0 else 0
    
    total_expected += expected
    total_completed += completed
    
    # Flag chromosomes with >10% chunk loss
    flag = " ***" if loss_pct > 10 else ""
    if loss_pct > 10:
        problematic_chromosomes.append((chrom, loss_pct, missing))
    
    print(f"{chrom:<12} {expected:<10} {completed:<10} {missing:<10} {loss_pct:>6.1f}%{flag}")

print("-" * 80)
print(f"{'TOTAL':<12} {total_expected:<10} {total_completed:<10} {total_expected - total_completed:<10} "
      f"{(total_expected - total_completed) / total_expected * 100:>6.1f}%")
print("=" * 80)

if problematic_chromosomes:
    print("\n⚠️  CHROMOSOMES WITH HIGH MISMATCH RATES:")
    print("-" * 80)
    for chrom, loss_pct, missing in sorted(problematic_chromosomes, key=lambda x: x[1], reverse=True):
        print(f"  • {chrom}: {loss_pct:.1f}% chunk loss ({missing} chunks failed)")
        
        # Provide specific insights
        if chrom in ['chr9', 'chr19']:
            print(f"    → Known for high GC content and repetitive regions")
        elif chrom in ['chr6']:
            print(f"    → Contains HLA region with high population diversity")
        elif chrom in ['chr17']:
            print(f"    → Contains inversion polymorphisms common in African populations")
    print()

print("\n📊 SUMMARY:")
print("-" * 80)
print(f"Total chunks expected: {total_expected}")
print(f"Total chunks completed: {total_completed}")
print(f"Total chunks lost: {total_expected - total_completed}")
print(f"Overall success rate: {total_completed / total_expected * 100:.1f}%")
print()

# Identify the specific failed chunks
print("🔍 INVESTIGATING FAILED CHUNKS:")
print("-" * 80)

# Find chunks that started but didn't complete
failed_chunks = []
with open(log_file, 'r') as f:
    content = f.read()
    
    # Find all CHECK_MISMATCH submissions
    mismatch_chunks = re.findall(r'CHECK_MISMATCH.*?(chr\d+_\d+_\d+)', content)
    
    # Find all completed imputations
    completed_chunks = re.findall(r'IMPUTE_MINIMAC4.*?Task completed.*?(chr\d+_\d+_\d+)', content)
    completed_set = set(completed_chunks)
    
    # Find the difference
    for chunk in set(mismatch_chunks):
        if chunk not in completed_set:
            failed_chunks.append(chunk)

if failed_chunks:
    print("Chunks that failed QC/mismatch checks:")
    for chunk in sorted(failed_chunks)[:10]:  # Show first 10
        chrom = chunk.split('_')[0]
        region = f"{chunk.split('_')[1]}-{chunk.split('_')[2]}"
        print(f"  • {chunk} (Region: {chrom}:{region})")
else:
    print("No specific failed chunks identified in log")

print("=" * 80)