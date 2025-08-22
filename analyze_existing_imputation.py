#!/usr/bin/env python3
"""
Standalone script to analyze existing imputation results
Generates pre/post comparison and R2 window analysis from existing files
"""

import os
import sys
import argparse
import glob
from pathlib import Path

def find_imputation_files(work_dir, pattern="*.info"):
    """Find imputation info files in work directory"""
    files = []
    for root, dirs, filenames in os.walk(work_dir):
        for filename in filenames:
            if filename.endswith('.info') or filename.endswith('.vcf.gz'):
                files.append(os.path.join(root, filename))
    return files

def find_vcf_pairs(work_dir):
    """Find pre and post-imputation VCF pairs"""
    pairs = []
    
    # Look for patterns that indicate pre/post VCFs
    pre_patterns = ['*target*.vcf.gz', '*input*.vcf.gz', '*qc*.vcf.gz']
    post_patterns = ['*imputed*.vcf.gz', '*output*.vcf.gz', '*minimac*.vcf.gz']
    
    pre_vcfs = []
    post_vcfs = []
    
    for pattern in pre_patterns:
        pre_vcfs.extend(glob.glob(os.path.join(work_dir, '**', pattern), recursive=True))
    
    for pattern in post_patterns:
        post_vcfs.extend(glob.glob(os.path.join(work_dir, '**', pattern), recursive=True))
    
    # Try to match based on chromosome
    for pre in pre_vcfs:
        pre_name = os.path.basename(pre)
        # Extract chromosome from filename
        import re
        chr_match = re.search(r'chr(\d+|X|Y)', pre_name)
        if chr_match:
            chr_id = chr_match.group(0)
            # Find matching post-imputation file
            for post in post_vcfs:
                if chr_id in os.path.basename(post):
                    pairs.append((pre, post))
                    break
    
    return pairs

def generate_r2_analysis(info_file, output_dir):
    """Generate R2 genomic window analysis from info file"""
    
    script = f"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read info file
try:
    df = pd.read_csv('{info_file}', sep='\\t', comment='#')
    
    # Look for R2 column
    r2_cols = ['Rsq', 'R2', 'RSQ', 'r2', 'INFO']
    r2_col = None
    for col in r2_cols:
        if col in df.columns:
            r2_col = col
            break
    
    if r2_col:
        # Basic statistics
        print(f"R2 Statistics for {os.path.basename(info_file)}:")
        print(f"  Mean R2: {{df[r2_col].mean():.3f}}")
        print(f"  Median R2: {{df[r2_col].median():.3f}}")
        print(f"  Min R2: {{df[r2_col].min():.3f}}")
        print(f"  Max R2: {{df[r2_col].max():.3f}}")
        print(f"  Variants with R2 > 0.8: {{(df[r2_col] > 0.8).sum()}} ({{(df[r2_col] > 0.8).sum()/len(df)*100:.1f}}%)")
        print(f"  Variants with R2 < 0.3: {{(df[r2_col] < 0.3).sum()}} ({{(df[r2_col] < 0.3).sum()/len(df)*100:.1f}}%)")
        
        # Create plot
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 2, 1)
        plt.hist(df[r2_col], bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('R2')
        plt.ylabel('Frequency')
        plt.title('R2 Distribution')
        plt.axvline(0.3, color='red', linestyle='--', label='Poor (0.3)')
        plt.axvline(0.8, color='green', linestyle='--', label='Good (0.8)')
        plt.legend()
        
        plt.subplot(2, 2, 2)
        # Cumulative distribution
        sorted_r2 = np.sort(df[r2_col])
        cumulative = np.arange(1, len(sorted_r2) + 1) / len(sorted_r2)
        plt.plot(sorted_r2, cumulative)
        plt.xlabel('R2')
        plt.ylabel('Cumulative Proportion')
        plt.title('Cumulative R2 Distribution')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        output_file = os.path.join('{output_dir}', f'{{os.path.basename(info_file).replace(".info", "")}}_r2_analysis.png')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"\\nPlot saved to: {{output_file}}")
        
except Exception as e:
    print(f"Error analyzing {{info_file}}: {{e}}")
"""
    
    # Execute the analysis
    exec(script)

def generate_comparison_report(pre_vcf, post_vcf, output_dir):
    """Generate pre/post imputation comparison"""
    
    print(f"\nComparing:")
    print(f"  Pre:  {pre_vcf}")
    print(f"  Post: {post_vcf}")
    
    # Create a simple comparison summary
    report_file = os.path.join(output_dir, 'comparison_summary.txt')
    
    with open(report_file, 'w') as f:
        f.write("PRE/POST IMPUTATION COMPARISON\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Pre-imputation:  {pre_vcf}\n")
        f.write(f"Post-imputation: {post_vcf}\n\n")
        
        # Use bcftools to get variant counts if available
        import subprocess
        try:
            pre_count = subprocess.check_output(
                f"bcftools view -H {pre_vcf} 2>/dev/null | wc -l", 
                shell=True, text=True
            ).strip()
            post_count = subprocess.check_output(
                f"bcftools view -H {post_vcf} 2>/dev/null | wc -l", 
                shell=True, text=True
            ).strip()
            
            f.write(f"Pre-imputation variants:  {pre_count}\n")
            f.write(f"Post-imputation variants: {post_count}\n")
            
            if pre_count.isdigit() and post_count.isdigit():
                gain = int(post_count) - int(pre_count)
                fold = float(post_count) / float(pre_count) if int(pre_count) > 0 else 0
                f.write(f"Variants gained: {gain}\n")
                f.write(f"Fold increase: {fold:.2f}x\n")
        except:
            f.write("Could not count variants (bcftools not available)\n")
    
    print(f"Report saved to: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze existing imputation results')
    parser.add_argument('--work-dir', default='/scratch3/users/mamana/nextflow-work',
                       help='Nextflow work directory')
    parser.add_argument('--output-dir', default='/scratch3/users/mamana/results/reports/standalone',
                       help='Output directory for reports')
    parser.add_argument('--analysis-type', choices=['r2', 'comparison', 'both'], 
                       default='both', help='Type of analysis to perform')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("IMPUTATION ANALYSIS TOOL")
    print("=" * 60)
    
    if args.analysis_type in ['r2', 'both']:
        print("\n1. R2 Genomic Window Analysis")
        print("-" * 30)
        
        # Find info files
        info_files = find_imputation_files(args.work_dir, '*.info')
        
        if info_files:
            print(f"Found {len(info_files)} info files")
            # Analyze up to 5 most recent files
            for info_file in sorted(info_files, key=os.path.getmtime, reverse=True)[:5]:
                print(f"\nAnalyzing: {os.path.basename(info_file)}")
                generate_r2_analysis(info_file, args.output_dir)
        else:
            print("No info files found")
    
    if args.analysis_type in ['comparison', 'both']:
        print("\n2. Pre/Post Imputation Comparison")
        print("-" * 30)
        
        # Find VCF pairs
        vcf_pairs = find_vcf_pairs(args.work_dir)
        
        if vcf_pairs:
            print(f"Found {len(vcf_pairs)} VCF pairs")
            for pre_vcf, post_vcf in vcf_pairs[:3]:  # Analyze up to 3 pairs
                generate_comparison_report(pre_vcf, post_vcf, args.output_dir)
        else:
            print("No matching VCF pairs found")
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print(f"Results saved to: {args.output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()