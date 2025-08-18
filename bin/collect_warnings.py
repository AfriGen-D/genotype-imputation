#!/usr/bin/env python3
"""
Collect and summarize warnings from pipeline execution
Parses log files for warnings, errors, and important messages
"""

import sys
import argparse
import re
from pathlib import Path
from datetime import datetime
import json

def parse_log_file(log_file):
    """Parse Nextflow log file for warnings and errors"""
    warnings = []
    errors = []
    info_messages = []
    
    warning_patterns = [
        r'WARN.*',
        r'WARNING:.*',
        r'Warning:.*',
        r'.*\[warn\].*',
        r'.*failed due to zero genetic distance.*',
        r'.*Eagle failed.*',
        r'.*Missing output file.*',
        r'.*low confidence.*',
        r'.*poor quality.*'
    ]
    
    error_patterns = [
        r'ERROR.*',
        r'Error:.*',
        r'.*\[error\].*',
        r'.*failed.*',
        r'.*exception.*',
        r'.*terminated.*'
    ]
    
    info_patterns = [
        r'.*Successfully completed.*',
        r'.*chunks? processed.*',
        r'.*variants? imputed.*',
        r'.*Process.*completed.*'
    ]
    
    with open(log_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Check for warnings
            for pattern in warning_patterns:
                if re.search(pattern, line, re.IGNORECASE):
                    warnings.append({
                        'line': line_num,
                        'message': line,
                        'type': 'warning'
                    })
                    break
            
            # Check for errors (exclude false positives)
            for pattern in error_patterns:
                if re.search(pattern, line, re.IGNORECASE):
                    # Filter out non-critical errors
                    if not any(skip in line.lower() for skip in ['error rate', 'error threshold', 'standard error']):
                        errors.append({
                            'line': line_num,
                            'message': line,
                            'type': 'error'
                        })
                        break
            
            # Check for important info
            for pattern in info_patterns:
                if re.search(pattern, line, re.IGNORECASE):
                    info_messages.append({
                        'line': line_num,
                        'message': line,
                        'type': 'info'
                    })
                    break
    
    return warnings, errors, info_messages

def categorize_warnings(warnings):
    """Categorize warnings by type"""
    categories = {
        'phasing': [],
        'imputation': [],
        'file_io': [],
        'quality': [],
        'memory': [],
        'other': []
    }
    
    for warning in warnings:
        msg = warning['message'].lower()
        
        if any(word in msg for word in ['eagle', 'phasing', 'phase', 'genetic distance']):
            categories['phasing'].append(warning)
        elif any(word in msg for word in ['impute', 'minimac', 'r2', 'maf']):
            categories['imputation'].append(warning)
        elif any(word in msg for word in ['file', 'missing', 'output', 'input', 'path']):
            categories['file_io'].append(warning)
        elif any(word in msg for word in ['quality', 'confidence', 'threshold', 'filter']):
            categories['quality'].append(warning)
        elif any(word in msg for word in ['memory', 'ram', 'heap', 'resource']):
            categories['memory'].append(warning)
        else:
            categories['other'].append(warning)
    
    return categories

def generate_summary(warnings, errors, info_messages):
    """Generate a summary of warnings and errors"""
    summary = {
        'timestamp': datetime.now().isoformat(),
        'total_warnings': len(warnings),
        'total_errors': len(errors),
        'total_info': len(info_messages),
        'categories': {},
        'recommendations': []
    }
    
    # Categorize warnings
    categorized = categorize_warnings(warnings)
    for category, items in categorized.items():
        if items:
            summary['categories'][category] = len(items)
    
    # Generate recommendations based on warnings
    if categorized['phasing']:
        summary['recommendations'].append(
            "Consider adjusting chunk sizes or using alternative phasing strategies for problematic regions"
        )
    
    if categorized['imputation']:
        summary['recommendations'].append(
            "Review imputation quality thresholds and consider filtering low-confidence variants"
        )
    
    if categorized['file_io']:
        summary['recommendations'].append(
            "Check file paths and ensure all required input files are available"
        )
    
    if categorized['memory']:
        summary['recommendations'].append(
            "Consider increasing memory allocation or reducing chunk sizes"
        )
    
    if len(errors) > 0:
        summary['recommendations'].append(
            f"Critical: {len(errors)} errors detected - review error messages for details"
        )
    
    return summary

def write_outputs(warnings, errors, info_messages, summary, args):
    """Write warning reports to files"""
    
    # Write detailed warnings file
    with open(args.output_warnings, 'w') as f:
        f.write(f"Pipeline Warnings Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Sample: {args.sample_id} | Reference: {args.ref_name}\n")
        f.write("=" * 80 + "\n\n")
        
        if warnings:
            f.write(f"WARNINGS ({len(warnings)} total)\n")
            f.write("-" * 40 + "\n")
            
            categorized = categorize_warnings(warnings)
            for category, items in categorized.items():
                if items:
                    f.write(f"\n{category.upper()} ({len(items)} warnings):\n")
                    for item in items[:10]:  # Show first 10 of each category
                        f.write(f"  Line {item['line']}: {item['message'][:150]}\n")
                    if len(items) > 10:
                        f.write(f"  ... and {len(items) - 10} more\n")
        else:
            f.write("No warnings detected\n")
        
        if errors:
            f.write(f"\n\nERRORS ({len(errors)} total)\n")
            f.write("-" * 40 + "\n")
            for error in errors[:20]:  # Show first 20 errors
                f.write(f"  Line {error['line']}: {error['message'][:150]}\n")
            if len(errors) > 20:
                f.write(f"  ... and {len(errors) - 20} more\n")
        
        f.write(f"\n\nRECOMMENDATIONS\n")
        f.write("-" * 40 + "\n")
        for rec in summary['recommendations']:
            f.write(f"• {rec}\n")
    
    # Write JSON summary
    with open(args.output_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Write TSV for easy parsing
    tsv_file = args.output_warnings.replace('.txt', '.tsv')
    with open(tsv_file, 'w') as f:
        f.write("type\tcategory\tline\tmessage\n")
        
        categorized = categorize_warnings(warnings)
        for category, items in categorized.items():
            for item in items:
                f.write(f"warning\t{category}\t{item['line']}\t{item['message']}\n")
        
        for error in errors:
            f.write(f"error\tgeneral\t{error['line']}\t{error['message']}\n")

def main():
    parser = argparse.ArgumentParser(description='Collect and summarize pipeline warnings')
    parser.add_argument('log_file', help='Nextflow log file or execution trace')
    parser.add_argument('--output-warnings', required=True, help='Output warnings file')
    parser.add_argument('--output-summary', required=True, help='Output summary JSON file')
    parser.add_argument('--sample-id', required=True, help='Sample identifier')
    parser.add_argument('--ref-name', required=True, help='Reference panel name')
    
    args = parser.parse_args()
    
    if not Path(args.log_file).exists():
        # Create empty outputs if log file doesn't exist
        print(f"Warning: Log file {args.log_file} not found. Creating empty reports.")
        
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_warnings': 0,
            'total_errors': 0,
            'total_info': 0,
            'categories': {},
            'recommendations': ['Log file not available - ensure pipeline logging is enabled']
        }
        
        with open(args.output_warnings, 'w') as f:
            f.write(f"Pipeline Warnings Report\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Sample: {args.sample_id} | Reference: {args.ref_name}\n")
            f.write("=" * 80 + "\n\n")
            f.write("Log file not available\n")
        
        with open(args.output_summary, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return
    
    print(f"Parsing log file: {args.log_file}")
    warnings, errors, info_messages = parse_log_file(args.log_file)
    
    print(f"Found {len(warnings)} warnings, {len(errors)} errors, {len(info_messages)} info messages")
    
    summary = generate_summary(warnings, errors, info_messages)
    write_outputs(warnings, errors, info_messages, summary, args)
    
    print(f"Warning report saved to {args.output_warnings}")
    print(f"Summary saved to {args.output_summary}")

if __name__ == "__main__":
    main()