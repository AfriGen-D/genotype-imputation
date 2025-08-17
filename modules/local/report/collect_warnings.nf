process COLLECT_WARNINGS {
    tag "$meta.id"
    label 'process_single'
    
    container 'quay.io/biocontainers/python:3.11'
    
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    tuple val(meta), path(log_file)
    
    output:
    tuple val(meta), path("*_warnings.txt"), emit: warnings
    tuple val(meta), path("*_skipped_chunks_summary.txt"), emit: summary
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import re
    from collections import defaultdict
    
    warnings_file = "${prefix}_warnings.txt"
    summary_file = "${prefix}_skipped_chunks_summary.txt"
    
    # Parse log file for warnings
    warnings = []
    skipped_overlap = []
    skipped_mismatch = []
    
    with open("${log_file}", 'r') as f:
        for line in f:
            if 'WARN' in line:
                warnings.append(line.strip())
                
                # Parse specific warning types
                if 'Skipping chunk' in line and 'insufficient overlap' in line:
                    match = re.search(r'Skipping chunk (\\S+) due to insufficient overlap', line)
                    if match:
                        skipped_overlap.append(match.group(1))
                        
                elif 'Skipping chunk' in line and 'high allele mismatch' in line:
                    match = re.search(r'Skipping chunk (\\S+) due to high allele mismatch', line)
                    if match:
                        skipped_mismatch.append(match.group(1))
    
    # Write all warnings
    with open(warnings_file, 'w') as f:
        f.write("PIPELINE WARNINGS\\n")
        f.write("=" * 50 + "\\n\\n")
        
        if warnings:
            for warning in warnings:
                f.write(warning + "\\n")
        else:
            f.write("No warnings generated during pipeline execution.\\n")
    
    # Write summary of skipped chunks
    with open(summary_file, 'w') as f:
        f.write("SKIPPED CHUNKS SUMMARY\\n")
        f.write("=" * 50 + "\\n\\n")
        
        if skipped_overlap:
            f.write(f"Chunks skipped due to insufficient overlap ({len(skipped_overlap)}):\\n")
            for chunk in sorted(skipped_overlap):
                f.write(f"  - {chunk}\\n")
            f.write("\\n")
        
        if skipped_mismatch:
            f.write(f"Chunks skipped due to high allele mismatch ({len(skipped_mismatch)}):\\n")
            for chunk in sorted(skipped_mismatch):
                f.write(f"  - {chunk}\\n")
            f.write("\\n")
        
        if not skipped_overlap and not skipped_mismatch:
            f.write("No chunks were skipped during processing.\\n")
        else:
            total_skipped = len(skipped_overlap) + len(skipped_mismatch)
            f.write(f"\\nTotal chunks skipped: {total_skipped}\\n")
            f.write(f"  - Due to overlap issues: {len(skipped_overlap)}\\n")
            f.write(f"  - Due to mismatch issues: {len(skipped_mismatch)}\\n")
    
    # Write versions
    import sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "No warnings" > ${prefix}_warnings.txt
    echo "No skipped chunks" > ${prefix}_skipped_chunks_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}