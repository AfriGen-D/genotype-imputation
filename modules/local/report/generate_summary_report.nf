process GENERATE_SUMMARY_REPORT {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    tuple val(meta), path(overlap_results), path(mismatch_results), path(imputation_logs)
    
    output:
    tuple val(meta), path("${prefix}_imputation_summary.txt"), emit: summary
    tuple val(meta), path("${prefix}_skipped_chunks.tsv"), emit: skipped
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import os
    import glob
    import re
    from datetime import datetime
    
    # Output files
    summary_file = "${prefix}_imputation_summary.txt"
    skipped_file = "${prefix}_skipped_chunks.tsv"
    
    # Initialize counters
    total_chunks = 0
    passed_overlap = 0
    failed_overlap = []
    passed_mismatch = 0
    failed_mismatch = []
    successfully_imputed = 0
    
    # Parse overlap results
    for overlap_file in glob.glob("*overlap.txt"):
        total_chunks += 1
        chunk_id = overlap_file.replace('.overlap.txt', '')
        
        # Check status file
        status_file = overlap_file.replace('.overlap.txt', '.status')
        if os.path.exists(status_file):
            with open(status_file, 'r') as f:
                status = f.read().strip()
                if status == "PASS":
                    passed_overlap += 1
                else:
                    with open(overlap_file, 'r') as of:
                        lines = of.readlines()
                        reason = "Insufficient overlap with reference panel"
                        if len(lines) > 1:
                            # Extract overlap ratio from file
                            for line in lines:
                                if 'Overall ratio' in line:
                                    ratio = line.split(':')[1].strip()
                                    reason = f"Insufficient overlap (ratio: {ratio})"
                                    break
                        failed_overlap.append((chunk_id, reason))
    
    # Parse mismatch results
    for mismatch_file in glob.glob("*mismatch.txt"):
        chunk_id = mismatch_file.replace('.mismatch.txt', '')
        status_file = mismatch_file.replace('.mismatch.txt', '.mismatch.status')
        
        if os.path.exists(status_file):
            with open(status_file, 'r') as f:
                status = f.read().strip()
                if status == "PASS":
                    passed_mismatch += 1
                else:
                    with open(mismatch_file, 'r') as mf:
                        lines = mf.readlines()
                        reason = "High allele mismatch with reference panel"
                        for line in lines:
                            if 'Mismatch rate' in line:
                                rate = line.split(':')[1].strip()
                                reason = f"High allele mismatch (rate: {rate})"
                                break
                        failed_mismatch.append((chunk_id, reason))
    
    # Count successfully imputed chunks
    for log_file in glob.glob("*.log"):
        with open(log_file, 'r') as f:
            content = f.read()
            if "Imputation successful" in content or "minimac4" in content.lower():
                successfully_imputed += 1
    
    # Generate summary report
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\\n")
        f.write("IMPUTATION PIPELINE SUMMARY REPORT\\n")
        f.write("=" * 80 + "\\n\\n")
        
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        f.write(f"Sample: ${meta.id}\\n\\n")
        
        f.write("CHUNK PROCESSING SUMMARY\\n")
        f.write("-" * 40 + "\\n")
        f.write(f"Total chunks generated: {total_chunks}\\n")
        f.write(f"Chunks passing overlap check: {passed_overlap} ({passed_overlap/total_chunks*100:.1f}%)\\n")
        f.write(f"Chunks failing overlap check: {len(failed_overlap)} ({len(failed_overlap)/total_chunks*100:.1f}%)\\n")
        f.write(f"Chunks passing mismatch check: {passed_mismatch}\\n")
        f.write(f"Chunks failing mismatch check: {len(failed_mismatch)}\\n")
        f.write(f"Successfully imputed chunks: {successfully_imputed}\\n\\n")
        
        if failed_overlap:
            f.write("CHUNKS SKIPPED DUE TO INSUFFICIENT OVERLAP\\n")
            f.write("-" * 40 + "\\n")
            for chunk, reason in failed_overlap:
                f.write(f"  - {chunk}: {reason}\\n")
            f.write("\\n")
        
        if failed_mismatch:
            f.write("CHUNKS SKIPPED DUE TO HIGH ALLELE MISMATCH\\n")
            f.write("-" * 40 + "\\n")
            for chunk, reason in failed_mismatch:
                f.write(f"  - {chunk}: {reason}\\n")
            f.write("\\n")
        
        f.write("QUALITY CONTROL PARAMETERS\\n")
        f.write("-" * 40 + "\\n")
        f.write(f"Minimum overlap ratio: ${params.minRatio}\\n")
        f.write(f"Minimum overlap count: 50 variants\\n")
        f.write(f"Maximum mismatch rate: ${params.max_mismatch_rate}\\n")
        f.write(f"R² threshold: ${params.r2_threshold}\\n\\n")
        
        f.write("RECOMMENDATIONS\\n")
        f.write("-" * 40 + "\\n")
        if len(failed_overlap) > total_chunks * 0.3:
            f.write("⚠️  High proportion of chunks failed overlap check.\\n")
            f.write("   Consider using a reference panel with better coverage for this region.\\n")
        if len(failed_mismatch) > 0:
            f.write("⚠️  Some chunks had high allele mismatch rates.\\n")
            f.write("   Check reference genome version and strand orientation.\\n")
        if successfully_imputed == total_chunks:
            f.write("✓  All chunks successfully imputed!\\n")
        
        f.write("\\n" + "=" * 80 + "\\n")
    
    # Generate skipped chunks TSV for downstream analysis
    with open(skipped_file, 'w') as f:
        f.write("chunk_id\\treason\\tstatus\\n")
        for chunk, reason in failed_overlap:
            f.write(f"{chunk}\\t{reason}\\tskipped_overlap\\n")
        for chunk, reason in failed_mismatch:
            f.write(f"{chunk}\\t{reason}\\tskipped_mismatch\\n")
    
    # Write versions
    import sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Sample: ${meta.id}" > ${prefix}_imputation_summary.txt
    echo "chunk_id\treason\tstatus" > ${prefix}_skipped_chunks.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}