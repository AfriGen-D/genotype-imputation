process CHECK_GLOBAL_MISMATCH {
    tag "global_mismatch_check"
    label 'process_single'
    
    container 'mamana/python-plotting:1.1.0'
    
    input:
    val mismatch_statuses
    
    output:
    path "global_mismatch_report.txt", emit: report
    path "pipeline_status.txt"       , emit: status
    path "versions.yml"              , emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def max_failed_chunks = params.max_failed_chunks ?: 100
    def max_failed_percent = params.max_failed_percent ?: 0.15
    """
    #!/usr/bin/env python3
    
    import sys
    import json
    
    # Parse input statuses from Nextflow list
    status_str = '''${mismatch_statuses}'''
    # Convert string representation of list to actual list
    status_str = status_str.strip().strip('[').strip(']')
    statuses = [s.strip().strip("'").strip('"') for s in status_str.split(',') if s.strip()]
    
    # Count status types
    total_chunks = len(statuses)
    passed_chunks = sum(1 for s in statuses if s == "PASS")
    warned_chunks = sum(1 for s in statuses if s == "WARN")
    failed_chunks = sum(1 for s in statuses if s == "FAIL")
    
    # Calculate failure rate
    failure_rate = failed_chunks / total_chunks if total_chunks > 0 else 0
    
    # Generate report
    with open("global_mismatch_report.txt", "w") as f:
        f.write("=" * 80 + "\\n")
        f.write("GLOBAL MISMATCH CHECK REPORT\\n")
        f.write("=" * 80 + "\\n\\n")
        
        f.write(f"Total chunks analyzed: {total_chunks}\\n")
        f.write(f"Chunks passed: {passed_chunks} ({passed_chunks/total_chunks*100:.1f}%)\\n")
        f.write(f"Chunks with warnings: {warned_chunks} ({warned_chunks/total_chunks*100:.1f}%)\\n")
        f.write(f"Chunks failed: {failed_chunks} ({failure_rate*100:.1f}%)\\n")
        f.write("\\n")
        
        f.write("Thresholds:\\n")
        f.write(f"  Max failed chunks allowed: ${max_failed_chunks}\\n")
        f.write(f"  Max failure percentage allowed: {${max_failed_percent}*100:.0f}%\\n")
        f.write("\\n")
        
        # Determine pipeline status
        if failed_chunks > ${max_failed_chunks}:
            status = "STOP"
            reason = f"Too many failed chunks: {failed_chunks} > ${max_failed_chunks}"
            f.write(f"❌ PIPELINE STATUS: STOP\\n")
            f.write(f"   Reason: {reason}\\n\\n")
            f.write("CRITICAL ERROR: Dataset appears incompatible with reference panel\\n")
            f.write("Possible causes:\\n")
            f.write("  • Wrong genome build (b37 vs b38)\\n")
            f.write("  • Population mismatch with reference panel\\n")
            f.write("  • Poor quality genotyping data\\n")
            f.write("  • Incorrect reference panel selection\\n")
        elif failure_rate > ${max_failed_percent}:
            status = "STOP"
            reason = f"Failure rate too high: {failure_rate*100:.1f}% > {${max_failed_percent}*100:.0f}%"
            f.write(f"❌ PIPELINE STATUS: STOP\\n")
            f.write(f"   Reason: {reason}\\n\\n")
            f.write("CRITICAL ERROR: High proportion of chunks failing\\n")
            f.write("Recommendations:\\n")
            f.write("  • Check genome build consistency\\n")
            f.write("  • Verify reference panel compatibility\\n")
            f.write("  • Review input data quality\\n")
        else:
            status = "CONTINUE"
            f.write(f"✅ PIPELINE STATUS: CONTINUE\\n")
            f.write(f"   Failed chunks within acceptable limits\\n\\n")
            if failed_chunks > 0:
                f.write(f"Note: {failed_chunks} chunks will be excluded from imputation\\n")
                f.write("This is acceptable for regions with:\\n")
                f.write("  • Population-specific variants\\n")
                f.write("  • Structural variations\\n")
                f.write("  • Low-quality genotyping\\n")
        
        f.write("\\n" + "=" * 80 + "\\n")
    
    # Write status for pipeline control
    with open("pipeline_status.txt", "w") as f:
        f.write(status)
    
    # Exit with error if pipeline should stop
    if status == "STOP":
        print(f"\\n{'='*80}", file=sys.stderr)
        print(f"PIPELINE TERMINATION REQUIRED", file=sys.stderr)
        print(f"{'='*80}", file=sys.stderr)
        print(f"Failed chunks: {failed_chunks}/{total_chunks} ({failure_rate*100:.1f}%)", file=sys.stderr)
        print(f"Threshold exceeded - stopping pipeline to prevent wasted computation", file=sys.stderr)
        print(f"{'='*80}\\n", file=sys.stderr)
        sys.exit(1)
    
    # Create versions file
    with open("versions.yml", "w") as f:
        f.write('\\"${task.process}\\":\\n')
        f.write('    python: ' + sys.version.split()[0] + '\\n')
    """
    
    stub:
    """
    echo "CONTINUE" > pipeline_status.txt
    touch global_mismatch_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}