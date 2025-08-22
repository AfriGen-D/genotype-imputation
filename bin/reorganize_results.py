#!/usr/bin/env python3
"""
Reorganize imputation results into a hierarchical structure:
dataset -> reference_panel -> {chunks, chromosome, genome, reports}
"""

import os
import shutil
import argparse
import re
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_metadata(filename):
    """Extract dataset, chromosome, and reference panel from filename."""
    metadata = {
        'dataset': None,
        'chromosome': None,
        'ref_panel': None,
        'chunk_start': None,
        'chunk_end': None
    }
    
    # Pattern: dataset_chr##_start_end_refpanel.extension
    # Example: awigen_500_b38_chr1_788531_25788530_H3AR6x.accuracy.txt
    
    # Extract reference panel (usually at the end before extension)
    ref_match = re.search(r'_(H3AR6x|1000G|HRC|CAAPA|other_ref)', filename)
    if ref_match:
        metadata['ref_panel'] = ref_match.group(1)
    
    # Extract chromosome
    chr_match = re.search(r'_chr(\d+|X|Y|MT)_', filename)
    if chr_match:
        metadata['chromosome'] = f"chr{chr_match.group(1)}"
    
    # Extract dataset (everything before _chr)
    if chr_match:
        dataset_part = filename[:chr_match.start()]
        metadata['dataset'] = dataset_part
    
    # Extract chunk coordinates if present
    coord_match = re.search(r'_chr\d+_(\d+)_(\d+)', filename)
    if coord_match:
        metadata['chunk_start'] = coord_match.group(1)
        metadata['chunk_end'] = coord_match.group(2)
    
    return metadata

def get_analysis_level(filepath, metadata):
    """Determine if file is chunk, chromosome, or genome level."""
    filename = os.path.basename(filepath)
    parent_dir = os.path.basename(os.path.dirname(filepath))
    
    # Check for specific patterns
    if 'genome' in parent_dir or 'final' in parent_dir:
        return 'genome'
    elif 'chromosome' in parent_dir or '.chr_' in filename:
        return 'chromosome'
    elif metadata['chunk_start'] and metadata['chunk_end']:
        return 'chunks'
    else:
        # Default based on file naming
        if '.chr_' in filename or '_chr_' in filename:
            return 'chromosome'
        elif 'genome' in filename or 'final' in filename:
            return 'genome'
        else:
            return 'chunks'

def reorganize_results(input_dir, output_dir, dry_run=False):
    """Reorganize results into hierarchical structure."""
    
    # Find all result files
    result_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(('.txt', '.tsv', '.pdf', '.html', '.json', '.png', '.jpg')):
                result_files.append(os.path.join(root, file))
    
    logger.info(f"Found {len(result_files)} result files to reorganize")
    
    # Track unique combinations
    datasets = set()
    ref_panels = set()
    
    # Process each file
    for filepath in result_files:
        filename = os.path.basename(filepath)
        
        # Skip version files
        if filename == 'versions.yml':
            continue
        
        # Extract metadata
        metadata = extract_metadata(filename)
        
        # Skip if we couldn't extract required metadata
        if not metadata['dataset'] or not metadata['ref_panel']:
            logger.warning(f"Could not extract metadata from: {filename}")
            continue
        
        datasets.add(metadata['dataset'])
        ref_panels.add(metadata['ref_panel'])
        
        # Determine analysis level
        level = get_analysis_level(filepath, metadata)
        
        # Build new path
        new_dir = os.path.join(
            output_dir,
            metadata['dataset'],
            metadata['ref_panel'],
            level
        )
        
        # Add chromosome subdirectory for chunk and chromosome levels
        if level in ['chunks', 'chromosome'] and metadata['chromosome']:
            new_dir = os.path.join(new_dir, metadata['chromosome'])
        
        # Create directory
        if not dry_run:
            os.makedirs(new_dir, exist_ok=True)
        
        # Copy or move file
        new_path = os.path.join(new_dir, filename)
        
        if dry_run:
            logger.info(f"Would move: {filepath}")
            logger.info(f"        to: {new_path}")
        else:
            # Check if source exists and is not the same as destination
            if os.path.exists(filepath) and filepath != new_path:
                # Create hard link instead of copying to save space
                try:
                    if os.path.exists(new_path):
                        os.remove(new_path)
                    os.link(filepath, new_path)
                    logger.debug(f"Linked {filename} to {new_dir}")
                except OSError:
                    # Fall back to copy if hard link fails (e.g., across filesystems)
                    shutil.copy2(filepath, new_path)
                    logger.debug(f"Copied {filename} to {new_dir}")
    
    # Create report summary structure
    if not dry_run:
        for dataset in datasets:
            for ref_panel in ref_panels:
                report_dir = os.path.join(output_dir, dataset, ref_panel, 'reports')
                os.makedirs(report_dir, exist_ok=True)
                
                # Look for final reports
                final_pattern = f"{dataset}_{ref_panel}.final_report"
                for ext in ['.html', '.pdf']:
                    final_files = [f for f in result_files if final_pattern + ext in f]
                    for final_file in final_files:
                        if os.path.exists(final_file):
                            dest = os.path.join(report_dir, os.path.basename(final_file))
                            if not os.path.exists(dest):
                                try:
                                    os.link(final_file, dest)
                                except OSError:
                                    shutil.copy2(final_file, dest)
    
    logger.info(f"Reorganization complete!")
    logger.info(f"Datasets found: {', '.join(sorted(datasets))}")
    logger.info(f"Reference panels found: {', '.join(sorted(ref_panels))}")
    
    # Print tree structure
    if not dry_run:
        print("\nNew directory structure:")
        for dataset in sorted(datasets):
            print(f"└── {dataset}/")
            for ref_panel in sorted(ref_panels):
                ref_path = os.path.join(output_dir, dataset, ref_panel)
                if os.path.exists(ref_path):
                    print(f"    └── {ref_panel}/")
                    for subdir in ['chunks', 'chromosome', 'genome', 'reports']:
                        subpath = os.path.join(ref_path, subdir)
                        if os.path.exists(subpath):
                            count = sum(1 for _ in Path(subpath).rglob('*') if _.is_file())
                            print(f"        ├── {subdir}/ ({count} files)")

def main():
    parser = argparse.ArgumentParser(description='Reorganize imputation results')
    parser.add_argument('--input-dir', default='/scratch3/users/mamana/results',
                       help='Input results directory')
    parser.add_argument('--output-dir', default='/scratch3/users/mamana/results_organized',
                       help='Output directory for reorganized results')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually doing it')
    
    args = parser.parse_args()
    
    reorganize_results(args.input_dir, args.output_dir, args.dry_run)

if __name__ == "__main__":
    main()