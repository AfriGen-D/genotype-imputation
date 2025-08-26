#!/usr/bin/env python3
"""
Migration script to convert old pipeline configurations to new v2.0 format
"""

import argparse
import json
import sys
from pathlib import Path
import re

def parse_old_config(config_file):
    """Parse old Nextflow config format"""
    params = {}
    
    with open(config_file, 'r') as f:
        content = f.read()
        
    # Extract params block
    params_match = re.search(r'params\s*\{([^}]+)\}', content, re.DOTALL)
    if not params_match:
        print(f"Warning: No params block found in {config_file}")
        return params
    
    params_content = params_match.group(1)
    
    # Parse individual parameters
    for line in params_content.split('\n'):
        line = line.strip()
        if '=' in line and not line.startswith('//'):
            # Remove comments
            line = line.split('//')[0].strip()
            
            # Parse key-value
            parts = line.split('=', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip().rstrip(';').strip()
                
                # Remove quotes
                if value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]
                elif value.startswith("'") and value.endswith("'"):
                    value = value[1:-1]
                
                params[key] = value
    
    return params

def map_parameters(old_params):
    """Map old parameter names to new ones"""
    
    # Parameter mapping dictionary
    param_map = {
        # Input/output
        'target_datasets': 'input',
        'ref_panels': 'reference_panels',
        'outDir': 'outdir',
        
        # QC parameters
        'site_miss': 'qc_max_missing',
        'hwe': 'qc_hwe_pvalue',
        'mac': 'qc_min_ac',
        'min_ac': 'qc_min_ac',
        
        # Imputation parameters  
        'NE': 'impute_ne',
        'impute_iter': None,  # Removed in v2
        'impute_burnin': None,  # Removed in v2
        'impute_info_cutoff': 'impute_info_cutoff',
        'chunk_size': 'chunk_size',
        'buffer_size': 'impute_buffer',
        
        # Tools
        'eagle_genetic_map': 'eagle_genetic_map',
        'reference_genome': 'reference_genome',
        
        # Chromosomes
        'chromosomes': 'chromosomes',
    }
    
    new_params = {}
    unmapped = []
    
    for key, value in old_params.items():
        if key in param_map:
            new_key = param_map[key]
            if new_key:  # Skip None mappings
                new_params[new_key] = value
        else:
            unmapped.append(key)
            new_params[key] = value  # Keep unmapped params
    
    # Set defaults for new required parameters
    if 'phasing_tool' not in new_params:
        new_params['phasing_tool'] = 'eagle'
    if 'imputation_tool' not in new_params:
        new_params['imputation_tool'] = 'minimac4'
    if 'report_level' not in new_params:
        new_params['report_level'] = 'detailed'
    if 'genome_build' not in new_params:
        # Try to infer from chromosome format
        if 'chromosomes' in new_params:
            if new_params['chromosomes'].startswith('chr'):
                new_params['genome_build'] = 'b38'
            else:
                new_params['genome_build'] = 'b37'
        else:
            new_params['genome_build'] = 'b38'
    
    return new_params, unmapped

def convert_reference_panels(ref_panels_str):
    """Convert reference panel format"""
    # Parse old format: ref_name:m3vcf:vcf or similar
    panels = []
    
    if isinstance(ref_panels_str, str):
        # Try to parse as a path to a config file
        if Path(ref_panels_str).exists():
            with open(ref_panels_str, 'r') as f:
                ref_config = json.load(f)
                for name, paths in ref_config.items():
                    panels.append([name, paths.get('m3vcf', ''), paths.get('vcf', '')])
        else:
            # Try to parse as inline format
            for panel in ref_panels_str.split(','):
                parts = panel.strip().split(':')
                if len(parts) >= 3:
                    panels.append([parts[0], parts[1], parts[2]])
    
    return panels

def create_samplesheet(target_datasets_str, output_path):
    """Create CSV samplesheet from old target_datasets format"""
    samples = []
    
    # Parse old format
    if isinstance(target_datasets_str, str):
        for dataset in target_datasets_str.split(','):
            parts = dataset.strip().split(':')
            if len(parts) >= 2:
                sample_name = parts[0]
                vcf_path = parts[1]
                samples.append({
                    'sample': sample_name,
                    'vcf': vcf_path,
                    'population': 'Unknown',
                    'sex': 'Unknown'
                })
    
    # Write CSV
    with open(output_path, 'w') as f:
        f.write('sample,vcf,population,sex\\n')
        for sample in samples:
            f.write(f"{sample['sample']},{sample['vcf']},{sample['population']},{sample['sex']}\\n")
    
    return output_path

def generate_new_config(params, output_file):
    """Generate new format configuration file"""
    
    config_content = f"""/*
========================================================================================
    Migrated Configuration for Imputation Pipeline v2.0
========================================================================================
*/

params {{
    // Input/output options
    input                      = '{params.get('input', 'samplesheet.csv')}'
    outdir                     = '{params.get('outdir', './results')}'
    
    // Reference genome options
    genome_build               = '{params.get('genome_build', 'b38')}'
    reference_panels           = {json.dumps(params.get('reference_panels', []))}
    eagle_genetic_map          = {params.get('eagle_genetic_map', 'null')}
    reference_genome           = {params.get('reference_genome', 'null')}
    chromosomes                = '{params.get('chromosomes', 'ALL')}'
    
    // QC options
    qc_min_ac                  = {params.get('qc_min_ac', 2)}
    qc_min_maf                 = {params.get('qc_min_maf', 0.01)}
    qc_max_missing             = {params.get('qc_max_missing', 0.05)}
    qc_hwe_pvalue              = {params.get('qc_hwe_pvalue', '1e-6')}
    
    // Phasing options
    phasing_tool               = '{params.get('phasing_tool', 'eagle')}'
    chunk_size                 = {params.get('chunk_size', 5000000)}
    
    // Imputation options
    imputation_tool            = '{params.get('imputation_tool', 'minimac4')}'
    impute_info_cutoff         = {params.get('impute_info_cutoff', 0.3)}
    impute_window              = {params.get('impute_window', 500000)}
    impute_ne                  = {params.get('impute_ne', 20000)}
    impute_buffer              = {params.get('impute_buffer', 250000)}
    
    // Reporting options
    report_level               = '{params.get('report_level', 'detailed')}'
    generate_plots             = {params.get('generate_plots', 'true')}
    concordance_analysis       = {params.get('concordance_analysis', 'false')}
    
    // Computational resources
    max_cpus                   = {params.get('max_cpus', 16)}
    max_memory                 = '{params.get('max_memory', '128.GB')}'
    max_time                   = '{params.get('max_time', '240.h')}'
}}
"""
    
    with open(output_file, 'w') as f:
        f.write(config_content)
    
    return output_file

def main():
    parser = argparse.ArgumentParser(
        description='Migrate old pipeline configuration to v2.0 format'
    )
    parser.add_argument(
        'config_file',
        help='Path to old configuration file'
    )
    parser.add_argument(
        '-o', '--output',
        default='migrated.config',
        help='Output configuration file (default: migrated.config)'
    )
    parser.add_argument(
        '--create-samplesheet',
        action='store_true',
        help='Create samplesheet from target_datasets'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Parse old configuration
    print(f"Parsing {args.config_file}...")
    old_params = parse_old_config(args.config_file)
    
    if args.verbose:
        print(f"Found {len(old_params)} parameters")
        print("Old parameters:", json.dumps(old_params, indent=2))
    
    # Map parameters
    print("Mapping parameters to v2.0 format...")
    new_params, unmapped = map_parameters(old_params)
    
    # Convert reference panels if present
    if 'reference_panels' in new_params:
        new_params['reference_panels'] = convert_reference_panels(
            new_params['reference_panels']
        )
    
    # Create samplesheet if requested
    if args.create_samplesheet and 'target_datasets' in old_params:
        samplesheet_path = Path(args.output).stem + '_samplesheet.csv'
        print(f"Creating samplesheet: {samplesheet_path}")
        create_samplesheet(old_params['target_datasets'], samplesheet_path)
        new_params['input'] = samplesheet_path
    
    # Generate new configuration
    print(f"Generating new configuration: {args.output}")
    generate_new_config(new_params, args.output)
    
    # Report results
    print("\\nMigration complete!")
    print(f"  - Mapped {len(new_params)} parameters")
    if unmapped:
        print(f"  - Unmapped parameters kept as-is: {', '.join(unmapped)}")
    print(f"  - New configuration saved to: {args.output}")
    
    if args.create_samplesheet:
        print(f"  - Samplesheet created: {new_params['input']}")
    
    print("\\nNext steps:")
    print("1. Review the migrated configuration")
    print("2. Update any file paths as needed")
    print("3. Test with: nextflow run workflows/main.nf -c", args.output)

if __name__ == '__main__':
    main()