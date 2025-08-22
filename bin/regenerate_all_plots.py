#!/usr/bin/env python3
"""
Regenerate all plots using proper data-driven plotting scripts.
This script identifies and runs the correct plotting scripts for chromosome and genome-level data.
"""

import os
import subprocess
import sys
import json
import argparse
from pathlib import Path

def run_command(cmd, description):
    """Run a command and report success/failure."""
    print(f"  Running: {description}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"    Success: {description}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"    Error in {description}: {e.stderr}")
        return False

def regenerate_chromosome_plots(results_dir, dataset, ref_name):
    """Regenerate chromosome-level plots."""
    print(f"\n=== Regenerating Chromosome-Level Plots ===")
    
    chr_dir = Path(results_dir) / dataset / ref_name / "chromosome"
    plots_dir = chr_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Find all chromosome JSON files
    json_files = list(chr_dir.glob("*_chr*_*.chr_summary.json"))
    print(f"Found {len(json_files)} chromosome summary files")
    
    container_path = "/users/mamana/genotype-imputation/singularity_cache/mamana-python-plotting-1.1.0.img"
    bin_dir = "/users/mamana/genotype-imputation/bin"
    
    success_count = 0
    
    for json_file in sorted(json_files):
        # Extract chromosome from filename
        filename_parts = json_file.stem.split('_')
        chr_part = None
        for part in filename_parts:
            if part.startswith('chr'):
                chr_part = part
                break
        
        if not chr_part:
            print(f"    Could not extract chromosome from {json_file.name}")
            continue
            
        chr_num = chr_part
        
        # Generate prefix (remove .chr_summary.json)
        prefix = json_file.stem.replace('.chr_summary', '')
        
        print(f"  Processing chromosome {chr_num}")
        
        # Change to chromosome directory
        os.chdir(str(chr_dir))
        
        # Run real chromosome performance plotting
        cmd = f"""singularity exec {container_path} python3 {bin_dir}/plot_chr_performance.py \\
            --chr-summary {json_file.name} \\
            --output-prefix {prefix} \\
            --ref-name {ref_name} \\
            --dataset {dataset} \\
            --chromosome {chr_num}"""
        
        if run_command(cmd, f"chromosome {chr_num} performance plot"):
            success_count += 1
            # Move plot to plots directory
            pdf_file = f"{prefix}_{ref_name}.chr_performance.pdf"
            if os.path.exists(pdf_file):
                os.rename(pdf_file, plots_dir / pdf_file)
                print(f"    Moved {pdf_file} to plots/")
    
    print(f"Successfully generated {success_count} chromosome performance plots")
    return success_count

def regenerate_genome_plots(results_dir, dataset, ref_name):
    """Regenerate genome-level plots."""
    print(f"\n=== Regenerating Genome-Level Plots ===")
    
    genome_dir = Path(results_dir) / dataset / ref_name / "genome"
    plots_dir = genome_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Find genome summary JSON
    json_files = list(genome_dir.glob("*.genome_summary.json"))
    
    if not json_files:
        print("    No genome summary JSON found")
        return 0
        
    json_file = json_files[0]
    print(f"Found genome summary: {json_file.name}")
    
    container_path = "/users/mamana/genotype-imputation/singularity_cache/mamana-python-plotting-1.1.0.img"
    bin_dir = "/users/mamana/genotype-imputation/bin"
    
    # Generate prefix
    prefix = json_file.stem.replace('.genome_summary', '')
    
    # Change to genome directory
    os.chdir(str(genome_dir))
    
    success_count = 0
    
    # List of genome plotting scripts to run
    genome_scripts = [
        "plot_genome_summary.py",
        "plot_genome_performance.py", 
        "plot_genome_maf_analysis.py",
        "plot_genome_r2_distribution.py",
        "plot_genome_chr_comparison.py",
        "plot_genome_r2_position.py",
        "plot_genome_maf_r2.py",
        "plot_genome_r2_snpcount.py"
    ]
    
    for script in genome_scripts:
        script_name = script.replace('.py', '').replace('plot_', '')
        cmd = f"""singularity exec {container_path} python3 {bin_dir}/{script} \\
            --genome-summary {json_file.name} \\
            --output-prefix {prefix} \\
            --ref-name {ref_name} \\
            --dataset {dataset}"""
        
        if run_command(cmd, f"genome {script_name} plots"):
            success_count += 1
            
    # Move all generated plots to plots directory
    for pdf_pattern in ["*.genome_*.pdf", "*_genome.pdf"]:
        for pdf_file in Path(".").glob(pdf_pattern):
            target = plots_dir / pdf_file.name
            pdf_file.rename(target)
            print(f"    Moved {pdf_file.name} to plots/")
    
    print(f"Successfully generated {success_count} genome-level plot sets")
    return success_count

def main():
    parser = argparse.ArgumentParser(description='Regenerate all imputation quality plots')
    parser.add_argument('--results-dir', default='/scratch3/users/mamana/results',
                       help='Results directory path')
    parser.add_argument('--dataset', default='awigen_500_b38',
                       help='Dataset name')
    parser.add_argument('--ref-name', default='H3AR6x',
                       help='Reference panel name')
    
    args = parser.parse_args()
    
    print("=== Plot Regeneration Script ===")
    print(f"Dataset: {args.dataset}")
    print(f"Reference: {args.ref_name}")
    print(f"Results directory: {args.results_dir}")
    
    # Check if directories exist
    base_dir = Path(args.results_dir) / args.dataset / args.ref_name
    if not base_dir.exists():
        print(f"Error: Base directory does not exist: {base_dir}")
        sys.exit(1)
    
    original_cwd = os.getcwd()
    total_success = 0
    
    try:
        # Regenerate chromosome plots
        chr_success = regenerate_chromosome_plots(args.results_dir, args.dataset, args.ref_name)
        total_success += chr_success
        
        # Regenerate genome plots  
        genome_success = regenerate_genome_plots(args.results_dir, args.dataset, args.ref_name)
        total_success += genome_success
        
        print(f"\n=== Summary ===")
        print(f"Total plots generated: {total_success}")
        print(f"Chromosome plots: {chr_success}")
        print(f"Genome plots: {genome_success}")
        
        if total_success > 0:
            print("\nPlots regenerated successfully! Check the plots/ directories for results.")
        else:
            print("\nNo plots were generated. Check the error messages above.")
            
    finally:
        os.chdir(original_cwd)

if __name__ == '__main__':
    main()