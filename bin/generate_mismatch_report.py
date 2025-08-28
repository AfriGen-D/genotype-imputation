#!/usr/bin/env python3
"""
Generate comprehensive mismatch report and determine if pipeline should continue
"""

import sys
import json
import argparse
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description='Check global mismatch statistics')
    parser.add_argument('--statuses', nargs='+', required=True, 
                       help='List of mismatch status values')
    parser.add_argument('--max-failed', type=int, default=100,
                       help='Maximum number of failed chunks allowed')
    parser.add_argument('--max-percent', type=float, default=0.15,
                       help='Maximum failure percentage allowed')
    parser.add_argument('--report', default='global_mismatch_report.txt',
                       help='Output report file')
    parser.add_argument('--status', default='pipeline_status.txt',
                       help='Output status file')
    
    args = parser.parse_args()
    
    # Parse statuses
    statuses = args.statuses
    
    # Count status types
    total_chunks = len(statuses)
    passed_chunks = sum(1 for s in statuses if s == "PASS")
    warned_chunks = sum(1 for s in statuses if s == "WARN")
    failed_chunks = sum(1 for s in statuses if s == "FAIL")
    
    # Calculate rates
    failure_rate = failed_chunks / total_chunks if total_chunks > 0 else 0
    warning_rate = warned_chunks / total_chunks if total_chunks > 0 else 0
    pass_rate = passed_chunks / total_chunks if total_chunks > 0 else 0
    
    # Determine pipeline status
    should_stop = False
    stop_reasons = []
    
    if failed_chunks > args.max_failed:
        should_stop = True
        stop_reasons.append(f"Too many failed chunks: {failed_chunks} > {args.max_failed}")
    
    if failure_rate > args.max_percent:
        should_stop = True
        stop_reasons.append(f"Failure rate too high: {failure_rate*100:.1f}% > {args.max_percent*100:.0f}%")
    
    # Generate detailed report
    with open(args.report, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GLOBAL MISMATCH CHECK REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        # Summary statistics
        f.write("📊 CHUNK STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total chunks analyzed:     {total_chunks:>6}\n")
        f.write(f"Chunks passed:            {passed_chunks:>6} ({pass_rate*100:>5.1f}%)\n")
        f.write(f"Chunks with warnings:     {warned_chunks:>6} ({warning_rate*100:>5.1f}%)\n")
        f.write(f"Chunks failed:            {failed_chunks:>6} ({failure_rate*100:>5.1f}%)\n")
        f.write("\n")
        
        # Threshold information
        f.write("⚙️  CONFIGURED THRESHOLDS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Max failed chunks allowed:    {args.max_failed}\n")
        f.write(f"Max failure percentage:       {args.max_percent*100:.0f}%\n")
        f.write("\n")
        
        # Pipeline decision
        f.write("🚦 PIPELINE DECISION\n")
        f.write("-" * 40 + "\n")
        
        if should_stop:
            f.write("❌ STATUS: STOP PIPELINE\n\n")
            f.write("Reasons:\n")
            for reason in stop_reasons:
                f.write(f"  • {reason}\n")
            f.write("\n")
            
            # Diagnostic information
            f.write("🔍 DIAGNOSTIC INFORMATION\n")
            f.write("-" * 40 + "\n")
            f.write("This level of mismatch indicates:\n")
            
            if failure_rate > 0.5:
                f.write("  ⚠️  SEVERE INCOMPATIBILITY (>50% failure)\n")
                f.write("  Likely causes:\n")
                f.write("    • Wrong genome build (b37 vs b38)\n")
                f.write("    • Wrong species or organism\n")
                f.write("    • Corrupted input files\n")
            elif failure_rate > 0.3:
                f.write("  ⚠️  MAJOR INCOMPATIBILITY (>30% failure)\n")
                f.write("  Likely causes:\n")
                f.write("    • Population mismatch with reference panel\n")
                f.write("    • Wrong reference panel selected\n")
                f.write("    • Chromosome naming convention mismatch\n")
            else:
                f.write("  ⚠️  SIGNIFICANT INCOMPATIBILITY (>15% failure)\n")
                f.write("  Likely causes:\n")
                f.write("    • Poor quality genotyping data\n")
                f.write("    • Reference panel missing populations\n")
                f.write("    • Structural variations in population\n")
            
            f.write("\n")
            f.write("📋 RECOMMENDED ACTIONS\n")
            f.write("-" * 40 + "\n")
            f.write("1. Verify genome build consistency\n")
            f.write("2. Check reference panel compatibility\n")
            f.write("3. Review input data quality metrics\n")
            f.write("4. Consider using a different reference panel\n")
            f.write("5. Check for sample/population stratification\n")
            
        else:
            f.write("✅ STATUS: CONTINUE PIPELINE\n\n")
            f.write(f"Failed chunks ({failed_chunks}) within acceptable limits\n\n")
            
            if failed_chunks > 0:
                f.write("📝 NOTES\n")
                f.write("-" * 40 + "\n")
                f.write(f"{failed_chunks} chunks will be excluded from imputation.\n")
                f.write("This is acceptable for:\n")
                f.write("  • Regions with population-specific variants\n")
                f.write("  • Structural variation hotspots\n")
                f.write("  • Low-complexity or repetitive regions\n")
                f.write("  • Regions with poor genotyping quality\n")
            
            if warned_chunks > 0:
                f.write(f"\n{warned_chunks} chunks proceeding with warnings.\n")
                f.write("These have marginal match quality but are still usable.\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    # Write status file
    with open(args.status, 'w') as f:
        f.write("STOP" if should_stop else "CONTINUE")
    
    # Exit with appropriate code
    if should_stop:
        print(f"\n{'='*80}", file=sys.stderr)
        print("PIPELINE TERMINATION REQUIRED", file=sys.stderr)
        print(f"{'='*80}", file=sys.stderr)
        print(f"Failed chunks: {failed_chunks}/{total_chunks} ({failure_rate*100:.1f}%)", file=sys.stderr)
        print(f"Exceeds threshold - stopping to prevent wasted computation", file=sys.stderr)
        print(f"See {args.report} for detailed analysis", file=sys.stderr)
        print(f"{'='*80}\n", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"✅ Mismatch check passed: {failed_chunks}/{total_chunks} chunks failed ({failure_rate*100:.1f}%)")
        print(f"Pipeline will continue with {total_chunks - failed_chunks} chunks")

if __name__ == "__main__":
    main()