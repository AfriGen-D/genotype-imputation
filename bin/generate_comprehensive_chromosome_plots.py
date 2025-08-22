#!/usr/bin/env python3

import json
import subprocess
import argparse
import sys
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def check_file_exists(file_path):
    """Check if a file exists and is readable."""
    if not os.path.exists(file_path):
        logger.error(f"File does not exist: {file_path}")
        return False
    if not os.access(file_path, os.R_OK):
        logger.error(f"File is not readable: {file_path}")
        return False
    return True

def run_plot_script(script_name, chr_summary_file, output_file, use_container=True):
    """Run a plotting script with error handling."""
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        logger.error(f"Plot script not found: {script_path}")
        return False
    
    try:
        if use_container:
            cmd = [
                'singularity', 'exec', 'mamana/python-plotting:1.1.0',
                'python3', str(script_path),
                '--chr-summary', chr_summary_file,
                '--output', output_file
            ]
        else:
            cmd = [
                'python3', str(script_path),
                '--chr-summary', chr_summary_file,
                '--output', output_file
            ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout:
            logger.info(f"Script output: {result.stdout.strip()}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run {script_name}: {e}")
        if e.stdout:
            logger.error(f"Stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"Stderr: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running {script_name}: {e}")
        return False

def create_comprehensive_chromosome_plots(chr_summary_file, output_dir, dataset, ref_name, chromosome, use_container=True):
    """Create all chromosome-level plots using the new comprehensive scripts."""
    
    # Validate input file
    if not check_file_exists(chr_summary_file):
        return False
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define plot types and their corresponding scripts
    plot_configs = [
        {
            'name': 'Original Summary',
            'script': 'plot_chromosome_summary.py',
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_performance.pdf',
            'description': 'Multi-panel performance dashboard (original 4-panel layout)'
        },
        {
            'name': 'MAF Analysis',
            'script': 'plot_chr_maf_analysis.py', 
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_maf_analysis.pdf',
            'description': 'Comprehensive MAF distribution and R² by frequency analysis'
        },
        {
            'name': 'R² Distribution',
            'script': 'plot_chr_r2_distribution.py',
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_r2_distribution.pdf', 
            'description': 'R² distribution histograms, cumulative plots, and quality thresholds'
        },
        {
            'name': 'R² by Position',
            'script': 'plot_chr_r2_position.py',
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_r2_position.pdf',
            'description': 'R² values by genomic position with variant density analysis'
        },
        {
            'name': 'MAF vs R²', 
            'script': 'plot_chr_maf_r2.py',
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_maf_r2.pdf',
            'description': 'MAF vs R² relationship analysis with trend lines and efficiency metrics'
        },
        {
            'name': 'SNP Count Analysis',
            'script': 'plot_chr_r2_snpcount.py',
            'output': f'{dataset}_{chromosome}_{ref_name}.chr_r2_snpcount.pdf',
            'description': 'R² vs SNP count correlation and chunk performance ranking'
        }
    ]
    
    successful_plots = []
    failed_plots = []
    
    logger.info(f"Creating comprehensive chromosome plots for {chromosome}")
    logger.info(f"Dataset: {dataset}, Reference: {ref_name}")
    logger.info(f"Output directory: {output_dir}")
    
    # Create each plot type
    for config in plot_configs:
        output_file = os.path.join(output_dir, config['output'])
        
        logger.info(f"Creating {config['name']} plot...")
        logger.info(f"Description: {config['description']}")
        
        # Handle the original script differently (it has different arguments)
        if config['script'] == 'plot_chromosome_summary.py':
            success = run_original_chromosome_script(
                chr_summary_file, output_file, dataset, ref_name, chromosome, use_container
            )
        else:
            success = run_plot_script(config['script'], chr_summary_file, output_file, use_container)
        
        if success:
            successful_plots.append({
                'name': config['name'],
                'file': output_file,
                'description': config['description']
            })
            logger.info(f"✓ Successfully created {config['name']} plot")
        else:
            failed_plots.append({
                'name': config['name'], 
                'script': config['script'],
                'description': config['description']
            })
            logger.error(f"✗ Failed to create {config['name']} plot")
    
    # Generate summary report
    generate_plot_summary(successful_plots, failed_plots, output_dir, dataset, ref_name, chromosome)
    
    return len(successful_plots), len(failed_plots)

def run_original_chromosome_script(chr_summary_file, output_file, dataset, ref_name, chromosome, use_container=True):
    """Run the original chromosome plot script with its specific arguments."""
    script_path = Path(__file__).parent / 'plot_chromosome_summary.py'
    
    if not script_path.exists():
        logger.error(f"Original chromosome script not found: {script_path}")
        return False
    
    try:
        # Extract output prefix from output file
        output_prefix = output_file.replace('.chr_plots.pdf', '').replace('.chr_performance.pdf', '')
        
        if use_container:
            cmd = [
                'singularity', 'exec', 'mamana/python-plotting:1.1.0',
                'python3', str(script_path),
                '--chr-summary', chr_summary_file,
                '--output-prefix', output_prefix,
                '--ref-name', ref_name,
                '--dataset', dataset,
                '--chromosome', chromosome
            ]
        else:
            cmd = [
                'python3', str(script_path),
                '--chr-summary', chr_summary_file,
                '--output-prefix', output_prefix,
                '--ref-name', ref_name,
                '--dataset', dataset,
                '--chromosome', chromosome
            ]
        
        logger.info(f"Running original script: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout:
            logger.info(f"Original script output: {result.stdout.strip()}")
        
        # The original script creates a file with a specific naming pattern
        original_output = f"{output_prefix}_{ref_name}.chr_plots.pdf"
        if os.path.exists(original_output) and original_output != output_file:
            # Rename to match our expected output
            os.rename(original_output, output_file)
            logger.info(f"Renamed {original_output} to {output_file}")
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to run original chromosome script: {e}")
        if e.stdout:
            logger.error(f"Stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"Stderr: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running original chromosome script: {e}")
        return False

def generate_plot_summary(successful_plots, failed_plots, output_dir, dataset, ref_name, chromosome):
    """Generate a summary report of created plots."""
    
    summary_file = os.path.join(output_dir, f'{dataset}_{chromosome}_{ref_name}.plot_summary.txt')
    
    try:
        with open(summary_file, 'w') as f:
            f.write(f"Chromosome {chromosome} Plot Generation Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Dataset: {dataset}\n")
            f.write(f"Reference Panel: {ref_name}\n")
            f.write(f"Chromosome: {chromosome}\n")
            f.write(f"Output Directory: {output_dir}\n\n")
            
            f.write(f"Successfully Created Plots ({len(successful_plots)}):\n")
            f.write("-" * 30 + "\n")
            for plot in successful_plots:
                f.write(f"• {plot['name']}\n")
                f.write(f"  File: {os.path.basename(plot['file'])}\n")
                f.write(f"  Description: {plot['description']}\n\n")
            
            if failed_plots:
                f.write(f"Failed Plots ({len(failed_plots)}):\n")
                f.write("-" * 20 + "\n")
                for plot in failed_plots:
                    f.write(f"• {plot['name']} (script: {plot['script']})\n")
                    f.write(f"  Description: {plot['description']}\n\n")
            
            f.write("Plot Descriptions:\n")
            f.write("-" * 18 + "\n")
            f.write("1. Original Summary: Traditional 4-panel layout with chunk performance, MAF distribution, R² by MAF, and info score statistics\n\n")
            f.write("2. MAF Analysis: 6-panel comprehensive analysis including MAF distribution (linear and log scale), R² by MAF with confidence intervals, violin plots, chunk variability, and summary statistics\n\n")
            f.write("3. R² Distribution: 6-panel R² quality analysis with histograms, cumulative distributions, quality threshold breakdowns, box plots by MAF, chunk variability, and detailed statistics\n\n")
            f.write("4. R² by Position: 4-panel genomic position analysis with scatter plots of R² vs position, variant density maps, imputation success rates, and correlation statistics\n\n")
            f.write("5. MAF vs R²: 6-panel relationship analysis with scatter plots, trend lines, confidence intervals, efficiency metrics, category comparisons, and correlation statistics\n\n")
            f.write("6. SNP Count Analysis: 6-panel variant count analysis with R² vs count correlations, count distributions, quartile analysis, efficiency plots, performance rankings, and summary statistics\n\n")
        
        logger.info(f"Plot summary saved to {summary_file}")
        
    except Exception as e:
        logger.error(f"Failed to generate plot summary: {e}")

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive chromosome-level plots')
    parser.add_argument('--chr-summary', required=True,
                       help='Path to chromosome summary JSON file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for plots')
    parser.add_argument('--dataset', required=True,
                       help='Dataset name')
    parser.add_argument('--ref-name', required=True,
                       help='Reference panel name')
    parser.add_argument('--chromosome', required=True,
                       help='Chromosome identifier (e.g., chr1, chr2)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not os.path.exists(args.chr_summary):
        logger.error(f"Chromosome summary file not found: {args.chr_summary}")
        sys.exit(1)
    
    # Create comprehensive plots
    successful, failed = create_comprehensive_chromosome_plots(
        args.chr_summary,
        args.output_dir,
        args.dataset,
        args.ref_name,
        args.chromosome,
        use_container=True  # Default to using container
    )
    
    # Print final summary
    print(f"\nChromosome {args.chromosome} Plot Generation Complete")
    print(f"{'='*50}")
    print(f"Successfully created: {successful} plots")
    print(f"Failed to create: {failed} plots")
    
    if successful > 0:
        print(f"\nPlots saved to: {args.output_dir}")
        print(f"Summary report: {args.dataset}_{args.chromosome}_{args.ref_name}.plot_summary.txt")
    
    if failed > 0:
        print(f"\nSome plots failed to generate. Check the logs above for details.")
        sys.exit(1)
    else:
        print(f"\n✓ All chromosome plots generated successfully!")

if __name__ == '__main__':
    main()