process COLLECT_WARNINGS {
    tag "$meta.id"
    label 'process_single'
    
    container 'mamana/python-plotting:1.0.0'
    
    publishDir "${params.outdir}/reports/${meta.id}", mode: 'copy'
    
    input:
    tuple val(meta), val(ref_name), path(log_file)
    
    output:
    tuple val(meta), val(ref_name), path("*_warnings.txt"), emit: warnings
    tuple val(meta), val(ref_name), path("*_warnings_summary.json"), emit: summary
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def warnings_out = "${prefix}_${ref_name}_warnings.txt"
    def summary_out = "${prefix}_${ref_name}_warnings_summary.json"
    """
    # Copy the Python script from bin directory
    cp ${projectDir}/bin/collect_warnings.py .
    
    # Run the warning collection script
    python3 collect_warnings.py \\
        ${log_file} \\
        --output-warnings ${warnings_out} \\
        --output-summary ${summary_out} \\
        --sample-id "${meta.id}" \\
        --ref-name "${ref_name}"
    
    # Ensure filesystem sync
    sync
    
    # Verify files exist
    if [ ! -f "${warnings_out}" ]; then
        echo "ERROR: Expected output file ${warnings_out} was not created!"
        exit 1
    fi
    
    if [ ! -f "${summary_out}" ]; then
        echo "ERROR: Expected output file ${summary_out} was not created!"
        exit 1
    fi
    
    echo "Output files verified:"
    ls -la ${warnings_out} ${summary_out}
    
    # Ensure files have proper permissions
    chmod 644 ${warnings_out} ${summary_out}
    
    # Small delay to ensure filesystem operations complete
    sleep 1
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def warnings_out = "${prefix}_${ref_name}_warnings.txt"
    def summary_out = "${prefix}_${ref_name}_warnings_summary.json"
    """
    echo "No warnings" > ${warnings_out}
    echo '{"total_warnings": 0, "total_errors": 0}' > ${summary_out}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}